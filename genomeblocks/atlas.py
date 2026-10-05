"""GIGGLE-like enrichment over a binned bin × track sparse index.

Tile the genome at fixed resolution and store one bit per (bin, track) cell
in a single CSR matrix. A query becomes one SpMV — ``M[query_bins].sum(0)``
returns per-track shared-bin counts in one C-level operation, vectorized
across thousands of tracks.

Three search modes (all share the same overlap-count primitive):

  * ``Atlas.search(query)``        Fisher 2x2 vs. genome null (bin units)
  * ``Atlas.search(query, ref=r)`` direct 2x2 query vs. ref per track
  * ``Atlas.bootstrap(query, n)``  shuffled-position null (chrom-aware)

Sized for ChIP-Atlas-scale collections (~25k mouse ChIP-seq tracks):
~2.5 GB resident at 1 kb resolution, sub-second per query.

References:
  Layer et al. *Nat. Methods* 15, 123-126 (2018).
"""
from __future__ import annotations

import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np

from ._table import TableMixin
from .genome import read_sizes

if TYPE_CHECKING:                                   # names used in annotations only
    import pandas as pd
    from scipy.sparse import csr_matrix
# scipy.sparse is imported where the matrix is built, so `import genomeblocks`
# does not pay for scipy. Loci.enrich / Loci.enrich_mc live in loci.py and call
# Atlas.search / Atlas.bootstrap.


_BED_EXT = (".bed", ".bed.gz", ".narrowPeak", ".narrowPeak.gz",
            ".broadPeak", ".broadPeak.gz")


class Atlas(TableMixin):
    """Sparse bin x track index for fast multi-file overlap enrichment.

    As a table it is its track table (name, n_peaks, n_bins + metadata):
    ``len``, ``shape``, ``columns``, ``head()``, ``describe()``,
    ``to_pandas()`` and the Arrow / dataframe protocols. Every query set is
    anything :func:`~genomeblocks.as_loci` takes."""

    def __init__(
        self,
        *,
        bin_size: int,
        chrom_names: Sequence[str],
        chrom_sizes: Dict[str, int],
        chrom_offsets: Dict[str, int],
        n_bins: int,
        track_names: Sequence[str],
        track_n_peaks: np.ndarray,
        track_n_bins: np.ndarray,
        M: "csr_matrix",
    ):
        self.bin_size = int(bin_size)
        self.chrom_names = list(chrom_names)
        self.chrom_sizes = dict(chrom_sizes)
        self.chrom_offsets = dict(chrom_offsets)
        self.n_bins = int(n_bins)
        self.track_names = list(track_names)
        self.track_n_peaks = np.asarray(track_n_peaks, dtype=np.int64)
        self.track_n_bins = np.asarray(track_n_bins, dtype=np.int64)
        self.M = M
        self.meta = None  # Optional[pandas.DataFrame] aligned to track_names
        self._chrom_max_bins = {
            c: (s + self.bin_size - 1) // self.bin_size
            for c, s in self.chrom_sizes.items()
        }

    # ---- construction --------------------------------------------------

    @classmethod
    def make(
        cls,
        paths: Union[str, Iterable[str]],
        *,
        chromsizes: Union[str, Dict[str, int], "object"],
        names: Optional[Sequence[str]] = None,
        bin_size: int = 1000,
        workers: Optional[int] = None,
        verbose: bool = True,
        meta: Union[str, "pd.DataFrame", None] = None,
        meta_columns: Optional[Sequence[str]] = None,
        meta_id_col: Optional[str] = None,
        meta_sep: str = "\t",
        name_pattern: Optional[str] = None,
    ) -> "Atlas":
        files = _resolve_paths(paths)
        if not files:
            raise ValueError(f"No BED files matched: {paths!r}")
        missing = [p for p in files if not os.path.isfile(p)]
        if missing:
            hint = (" Did you pass `os.listdir(dir)`? That returns bare "
                    "filenames -- pass the directory string directly, or "
                    "join it to each name.") if not os.path.isabs(missing[0]) else ""
            raise FileNotFoundError(
                f"{len(missing)}/{len(files)} input files not found. "
                f"First missing: {missing[0]!r}.{hint}"
            )
        if names is None:
            names = [_basename(p, pattern=name_pattern) for p in files]
        elif len(names) != len(files):
            raise ValueError("names length must match number of input files")

        chrom_sizes = {str(k): int(v) for k, v in read_sizes(chromsizes).items()}
        if not chrom_sizes:
            raise ValueError("chromsizes is empty: pass a {chrom: length} dict, a .chrom.sizes path, "
                             "a Genome with sizes or a cooler")
        chrom_names = list(chrom_sizes.keys())
        chrom_offsets: Dict[str, int] = {}
        cum = 0
        for c in chrom_names:
            chrom_offsets[c] = cum
            cum += (chrom_sizes[c] + bin_size - 1) // bin_size
        n_bins = cum

        if workers is None:
            workers = max(1, (os.cpu_count() or 2) // 2)

        n_tracks = len(files)
        per_track_bins: List[Optional[np.ndarray]] = [None] * n_tracks
        per_track_n_peaks = np.zeros(n_tracks, dtype=np.int64)

        args = [
            (i, p, chrom_offsets, chrom_sizes, bin_size)
            for i, p in enumerate(files)
        ]

        if workers > 1 and n_tracks > 1:
            from tqdm import tqdm
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futures = {ex.submit(_build_worker, a): a[0] for a in args}
                it = as_completed(futures)
                if verbose:
                    it = tqdm(it, total=n_tracks, desc="[atlas index]")
                for fut in it:
                    tid, bins, n_peaks = fut.result()
                    per_track_bins[tid] = bins
                    per_track_n_peaks[tid] = n_peaks
        else:
            it = args
            if verbose:
                from tqdm import tqdm
                it = tqdm(it, total=n_tracks, desc="[atlas index]")
            for a in it:
                tid, bins, n_peaks = _build_worker(a)
                per_track_bins[tid] = bins
                per_track_n_peaks[tid] = n_peaks

        # Assemble CSC (bins, tracks): tracks become columns, sorted bin ids
        # become row indices. Single-pass concat keeps peak memory at ~nnz*4.
        lens = np.array(
            [0 if b is None else len(b) for b in per_track_bins],
            dtype=np.int64,
        )
        nnz = int(lens.sum())
        if nnz == 0:
            raise ValueError(
                "Atlas index is empty -- no BED intervals fell on supplied chroms.",
            )

        idx_dtype = np.int64 if n_bins > 2_147_483_647 else np.int32
        indptr = np.zeros(n_tracks + 1, dtype=np.int64)
        np.cumsum(lens, out=indptr[1:])
        indices = np.empty(nnz, dtype=idx_dtype)
        for i, b in enumerate(per_track_bins):
            if b is None or len(b) == 0:
                continue
            indices[indptr[i]:indptr[i + 1]] = b
        data = np.ones(nnz, dtype=np.uint8)
        track_n_bins = lens.astype(np.int64)
        del per_track_bins

        from scipy.sparse import csc_matrix
        M = csc_matrix(
            (data, indices, indptr),
            shape=(n_bins, n_tracks),
        ).tocsr()

        atlas = cls(
            bin_size=bin_size,
            chrom_names=chrom_names,
            chrom_sizes=chrom_sizes,
            chrom_offsets=chrom_offsets,
            n_bins=n_bins,
            track_names=list(names),
            track_n_peaks=per_track_n_peaks,
            track_n_bins=track_n_bins,
            M=M,
        )
        if meta is not None:
            atlas.attach_meta(
                meta,
                id_col=meta_id_col,
                columns=meta_columns,
                sep=meta_sep,
            )
        return atlas

    # ---- metadata ------------------------------------------------------

    def attach_meta(
        self,
        meta: Union[str, "pd.DataFrame"],
        *,
        id_col: Optional[str] = None,
        columns: Optional[Sequence[str]] = None,
        sep: str = "\t",
    ) -> "Atlas":
        """Attach per-track metadata aligned to ``self.track_names``.

        ``meta`` can be a path to a TSV/CSV or an in-memory ``DataFrame``.
        ``columns`` supplies headers for a header-less file (also implies
        ``header=None`` on read). ``id_col`` selects the column that joins
        against ``track_names``; defaults to the first column.

        Tracks with no row in the metadata get NaN; metadata rows for
        unknown tracks are dropped.
        """
        import pandas as pd
        if isinstance(meta, str):
            df = pd.read_csv(
                meta, sep=sep,
                header=None if columns is not None else "infer",
                names=list(columns) if columns is not None else None,
                dtype=str,
            )
        else:
            df = meta.copy()
        if id_col is None:
            id_col = df.columns[0]
        df = df.set_index(id_col)
        df.index = df.index.astype(str)
        self.meta = df.reindex(self.track_names)
        return self

    # ---- bin conversion ------------------------------------------------

    def _intervals_to_bin_ranges(self, loci) -> np.ndarray:
        """Map each interval to ``[abs_start_bin, abs_end_bin)`` (columns only).

        ``loci`` is anything :func:`~genomeblocks.as_loci` takes. Intervals on
        chromosomes missing from the atlas are dropped. Returns an ``(n, 2)``
        int64 array; empty if no interval qualifies.
        """
        from .interop import as_loci
        L = as_loci(loci)
        if not len(L):
            return np.zeros((0, 2), dtype=np.int64)
        bs = self.bin_size
        names = L.genome.names
        off = np.array([self.chrom_offsets.get(n, -1) for n in names], np.int64)[L.codes]
        cmax = np.array([self._chrom_max_bins.get(n, 0) for n in names], np.int64)[L.codes]
        s = np.maximum(L.starts // bs, 0)
        e = np.maximum((L.ends - 1) // bs + 1, s + 1)
        s, e = np.minimum(s, cmax), np.minimum(e, cmax)
        ok = (off >= 0) & (s < e)
        return np.column_stack([off[ok] + s[ok], off[ok] + e[ok]])

    def _bins_union(self, ranges: np.ndarray) -> np.ndarray:
        """Sorted unique bin ids covered by any of the input ranges."""
        if len(ranges) == 0:
            return np.zeros(0, dtype=np.int64)
        starts = ranges[:, 0]
        lengths = ranges[:, 1] - ranges[:, 0]
        return np.unique(_explode_ranges(starts, lengths))

    def _overlap_counts(self, q_bins: np.ndarray) -> np.ndarray:
        if len(q_bins) == 0:
            return np.zeros(len(self.track_names), dtype=np.int64)
        sub = self.M[q_bins]
        return np.asarray(sub.sum(axis=0)).ravel().astype(np.int64)

    # ---- GIGGLE search -------------------------------------------------

    def search(
        self,
        query,
        *,
        ref=None,
        alternative: str = "two-sided",
    ):
        """GIGGLE-style enrichment of each track over the query.

        ``query`` / ``ref`` are anything :func:`~genomeblocks.as_loci` takes.
        With ``ref`` supplied, the c/d cells come from the reference set
        instead of the genome null -- "factors more enriched in query than
        in ref". Returns a ``pandas.DataFrame`` sorted by ``giggle_score``.
        """
        import pandas as pd

        q_bins = self._bins_union(self._intervals_to_bin_ranges(query))
        if len(q_bins) == 0:
            raise ValueError("Query has no bins on the atlas's chromosomes.")
        a = self._overlap_counts(q_bins)
        n_q = int(len(q_bins))

        if ref is not None:
            r_bins = self._bins_union(self._intervals_to_bin_ranges(ref))
            if len(r_bins) == 0:
                raise ValueError("Reference has no bins on the atlas's chromosomes.")
            c = self._overlap_counts(r_bins)
            n_r = int(len(r_bins))
            b = n_q - a
            d = n_r - c
        else:
            n_r = None
            c = self.track_n_bins - a
            b = n_q - a
            d = self.n_bins - n_q - c

        a, b, c, d = (np.maximum(x, 0).astype(np.int64) for x in (a, b, c, d))
        log_or, p, score = _fisher_vec(a, b, c, d, alternative=alternative)

        n_tracks = len(self.track_names)
        out = {"name": self.track_names, "n_query_bins": np.full(n_tracks, n_q, dtype=np.int64)}
        if n_r is not None:
            out["n_ref_bins"] = np.full(n_tracks, n_r, dtype=np.int64)
        out.update({
            "track_n_bins": self.track_n_bins,
            "track_n_peaks": self.track_n_peaks,
            "overlaps": a,
            "log2_odds": log_or,
            "p": p,
            "giggle_score": score,
        })
        df = pd.DataFrame(out)
        df = self._with_meta(df)
        return df.sort_values("giggle_score", ascending=False).reset_index(drop=True)

    # ---- Monte Carlo bootstrap ----------------------------------------

    def bootstrap(
        self,
        query,
        *,
        n: int = 10,
        pool=None,
        sample: Optional[int] = None,
        replace: bool = False,
        keep_chrom: bool = True,
        seed: Optional[int] = None,
        verbose: bool = True,
    ):
        """Compare observed overlaps to a resampled null.

        Three knobs that interact:

        * ``pool=None`` (default) -- null is per-interval position shuffle
          (``keep_chrom`` preserves the original chromosome). Each group
          shuffles independently.
        * ``pool=<Loci>`` -- null is sampled from a curated CRE universe.
          Use this when query is a subset of pool and you want to control
          for the universe's bias (LOLA / regioneR style).
        * ``sample=<int>`` -- per iteration, subsample BOTH every query
          group and the pool down to ``sample`` regions. The pool draw is
          **shared across all groups** within each iteration, so cross-
          group comparisons are paired against the same null. Use this
          for differently-sized query groups.

        ``query`` may be a single interval set (anything
        :func:`~genomeblocks.as_loci` takes) or a ``{group: intervals}``
        dict. With a dict of groups, the result is a long-format DataFrame
        with a ``group`` column.

        With ``sample`` set, ``observed`` and ``expected`` are means over
        the bootstrap distribution of subsamples; ``z`` is the paired z
        of (obs - null) per iteration.
        """
        import pandas as pd

        if n <= 0:
            raise ValueError("n must be positive.")
        rng = np.random.default_rng(seed)

        from .interop import as_loci
        single = not (isinstance(query, dict) and not _is_columns(query))
        groups = {"query": as_loci(query)} if single else {k: as_loci(v) for k, v in query.items()}
        pool = as_loci(pool) if pool is not None else None
        if not groups:
            raise ValueError("query dict is empty.")

        # per-group bin ranges (one row per interval)
        group_ranges: Dict[str, np.ndarray] = {}
        for g, q in groups.items():
            r = self._intervals_to_bin_ranges(q)
            if len(r) == 0:
                raise ValueError(
                    f"query group {g!r} has no bins on atlas chroms.",
                )
            group_ranges[g] = r

        if pool is not None:
            pool_ranges = self._intervals_to_bin_ranges(pool)
            if len(pool_ranges) == 0:
                raise ValueError("pool has no bins on atlas chromosomes.")
            n_pool = int(len(pool_ranges))
        else:
            pool_ranges = None
            n_pool = 0
            chrom_idx = {c: i for i, c in enumerate(self.chrom_names)}
            chrom_offsets_arr = np.array(
                [self.chrom_offsets[c] for c in self.chrom_names],
                dtype=np.int64,
            )
            chrom_bins = np.array(
                [self._chrom_max_bins[c] for c in self.chrom_names],
                dtype=np.int64,
            )
            chrom_p = chrom_bins / chrom_bins.sum()

        # If subsampling, validate sizes once
        if sample is not None:
            if sample <= 0:
                raise ValueError("sample must be positive.")
            for g, r in group_ranges.items():
                if len(r) < sample and not replace:
                    raise ValueError(
                        f"query group {g!r} has {len(r)} regions but "
                        f"sample={sample}. Pass replace=True or lower sample.",
                    )
            if pool is not None and n_pool < sample and not replace:
                raise ValueError(
                    f"pool has {n_pool} regions but sample={sample}. "
                    f"Pass replace=True or lower sample.",
                )

        n_tracks = len(self.track_names)

        # Per-group accumulators
        obs_acc = {g: np.zeros((n, n_tracks), dtype=np.int64) for g in groups}
        null_acc = {g: np.zeros((n, n_tracks), dtype=np.int64) for g in groups}

        # If query isn't subsampled, observed is fixed -- compute once
        fixed_observed: Dict[str, np.ndarray] = {}
        if sample is None:
            for g, r in group_ranges.items():
                bins_u = np.unique(_explode_ranges(r[:, 0], r[:, 1] - r[:, 0]))
                fixed_observed[g] = self._overlap_counts(bins_u)

        # Pre-stage the per-group genome-shuffle inputs (only if pool is None)
        gs: Dict[str, dict] = {}
        if pool is None:
            bs = self.bin_size
            for g, q in groups.items():
                ci = np.array([chrom_idx.get(n, -1) for n in q.genome.names], np.int64)[q.codes]
                ok = ci >= 0
                s_b = q.starts[ok] // bs
                e_b = np.maximum((q.ends[ok] - 1) // bs + 1, s_b + 1)
                gs[g] = {
                    "c_ids": ci[ok],
                    "lens": np.minimum(e_b - s_b, chrom_bins[ci[ok]]).astype(np.int64),
                }

        it = range(n)
        if verbose:
            from tqdm import tqdm
            it = tqdm(it, desc="[atlas bootstrap]")

        for i in it:
            # ---- null draw(s) for this iteration ----
            if pool is not None and sample is not None:
                # ONE shared pool draw, broadcast to every group
                idx = rng.choice(n_pool, size=sample, replace=replace)
                sampled = pool_ranges[idx]
                bins_u = np.unique(
                    _explode_ranges(sampled[:, 0],
                                    sampled[:, 1] - sampled[:, 0]),
                )
                shared_null = self._overlap_counts(bins_u)
                for g in groups:
                    null_acc[g][i] = shared_null
            elif pool is not None:
                # No subsample: each group draws len(group) from pool
                # (independent draws -- not paired across groups)
                for g, r in group_ranges.items():
                    k = len(r)
                    if k > n_pool and not replace:
                        raise ValueError(
                            f"query group {g!r} has {k} regions but pool "
                            f"has {n_pool}. Pass replace=True or set sample.",
                        )
                    idx = rng.choice(n_pool, size=k, replace=replace)
                    sampled = pool_ranges[idx]
                    bins_u = np.unique(
                        _explode_ranges(sampled[:, 0],
                                        sampled[:, 1] - sampled[:, 0]),
                    )
                    null_acc[g][i] = self._overlap_counts(bins_u)
            else:
                # Genome-shuffle null per group (independent)
                for g in groups:
                    c_ids_arr = gs[g]["c_ids"]
                    lens_arr = gs[g]["lens"]
                    cs = c_ids_arr if keep_chrom else rng.choice(
                        len(self.chrom_names), size=len(c_ids_arr), p=chrom_p,
                    )
                    cb = chrom_bins[cs]
                    free = np.maximum(cb - lens_arr, 0)
                    local = rng.integers(0, free + 1)
                    starts = chrom_offsets_arr[cs] + local
                    q_bins = np.unique(_explode_ranges(starts, lens_arr))
                    null_acc[g][i] = self._overlap_counts(q_bins)

            # ---- observed for this iteration ----
            if sample is None:
                for g in groups:
                    obs_acc[g][i] = fixed_observed[g]
            else:
                for g, r in group_ranges.items():
                    idx = rng.choice(len(r), size=sample, replace=replace)
                    sampled = r[idx]
                    bins_u = np.unique(
                        _explode_ranges(sampled[:, 0],
                                        sampled[:, 1] - sampled[:, 0]),
                    )
                    obs_acc[g][i] = self._overlap_counts(bins_u)

        # ---- per-group summary ----
        out_frames: List["pd.DataFrame"] = []
        for g in groups:
            obs = obs_acc[g]
            null = null_acc[g]
            mean_obs = obs.mean(axis=0)
            mean_null = null.mean(axis=0)
            obs_std = obs.std(axis=0, ddof=1) if n > 1 else np.zeros(n_tracks)
            null_std = null.std(axis=0, ddof=1) if n > 1 else np.zeros(n_tracks)
            diff = obs - null
            mean_diff = diff.mean(axis=0)
            std_diff = (diff.std(axis=0, ddof=1)
                        if n > 1 else np.zeros(n_tracks))
            z = np.divide(mean_diff, std_diff,
                          out=np.zeros_like(mean_diff, dtype=float),
                          where=std_diff > 0)
            log2fc = np.log2((mean_obs + 1.0) / (mean_null + 1.0))
            n_le_zero = (diff <= 0).sum(axis=0)
            p_emp = (1 + n_le_zero) / (n + 1)

            df = pd.DataFrame({
                "name": self.track_names,
                "observed": mean_obs,
                "expected": mean_null,
                "obs_std": obs_std,
                "null_std": null_std,
                "log2fc": log2fc,
                "z": z,
                "p_emp": p_emp,
                "track_n_bins": self.track_n_bins,
            })
            df = self._with_meta(df)
            df = df.sort_values("z", ascending=False).reset_index(drop=True)
            if not single:
                df.insert(0, "group", g)
            out_frames.append(df)

        return out_frames[0] if single else pd.concat(out_frames, ignore_index=True)

    # ---- helpers -----------------------------------------------------

    def _with_meta(self, df):
        if self.meta is None or self.meta.empty:
            return df
        return df.merge(
            self.meta, left_on="name", right_index=True, how="left",
        )

    # ---- persistence ---------------------------------------------------

    def save(self, path: str) -> None:
        arrays = dict(
            bin_size=np.int64(self.bin_size),
            n_bins=np.int64(self.n_bins),
            chrom_names=np.array(self.chrom_names),
            chrom_lens=np.array(
                [self.chrom_sizes[c] for c in self.chrom_names],
                dtype=np.int64,
            ),
            chrom_offsets=np.array(
                [self.chrom_offsets[c] for c in self.chrom_names],
                dtype=np.int64,
            ),
            track_names=np.array(self.track_names),
            track_n_peaks=self.track_n_peaks,
            track_n_bins=self.track_n_bins,
            M_data=self.M.data,
            M_indices=self.M.indices,
            M_indptr=self.M.indptr,
            M_shape=np.array(self.M.shape, dtype=np.int64),
        )
        if self.meta is not None and not self.meta.empty:
            cols = list(self.meta.columns)
            arrays["meta_columns"] = np.array([str(c) for c in cols], dtype=str)
            for col in cols:
                arrays[f"meta__{col}"] = np.array(
                    self.meta[col].fillna("").tolist(), dtype=str,
                )
        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path: str) -> "Atlas":
        d = np.load(path, allow_pickle=False)
        chrom_names = [str(x) for x in d["chrom_names"].tolist()]
        chrom_lens = d["chrom_lens"].astype(np.int64)
        chrom_offsets = d["chrom_offsets"].astype(np.int64)
        chrom_sizes = {c: int(l) for c, l in zip(chrom_names, chrom_lens)}
        chrom_off = {c: int(o) for c, o in zip(chrom_names, chrom_offsets)}
        from scipy.sparse import csr_matrix
        M = csr_matrix(
            (d["M_data"], d["M_indices"], d["M_indptr"]),
            shape=tuple(int(x) for x in d["M_shape"].tolist()),
        )
        track_names = [str(x) for x in d["track_names"].tolist()]
        atlas = cls(
            bin_size=int(d["bin_size"]),
            chrom_names=chrom_names,
            chrom_sizes=chrom_sizes,
            chrom_offsets=chrom_off,
            n_bins=int(d["n_bins"]),
            track_names=track_names,
            track_n_peaks=d["track_n_peaks"].astype(np.int64),
            track_n_bins=d["track_n_bins"].astype(np.int64),
            M=M,
        )
        if "meta_columns" in d.files:
            import pandas as pd
            cols = [str(x) for x in d["meta_columns"].tolist()]
            data = {c: [str(x) for x in d[f"meta__{c}"].tolist()] for c in cols}
            df = pd.DataFrame(data, index=track_names)
            atlas.meta = df.replace("", np.nan)
        return atlas

    # ---- the track table -------------------------------------------------

    @property
    def columns(self) -> list:
        meta = [] if self.meta is None else [str(c) for c in self.meta.columns]
        return ["name", "n_peaks", "n_bins"] + meta

    def to_pandas(self):
        """One row per track: name, n_peaks, n_bins and the metadata columns."""
        import pandas as pd
        df = pd.DataFrame({"name": self.track_names, "n_peaks": self.track_n_peaks, "n_bins": self.track_n_bins})
        if self.meta is not None and not self.meta.empty:
            m = self.meta.reset_index(drop=True)
            for c in m.columns:
                df[str(c)] = m[c].to_numpy()
        return df

    def to_arrow(self):
        import pyarrow as pa
        return pa.Table.from_pandas(self.to_pandas(), preserve_index=False)

    def head(self, n: int = 5):
        return self.to_pandas().head(n)

    def tail(self, n: int = 5):
        return self.to_pandas().tail(n)

    def describe(self):
        import pandas as pd
        nb = self.track_n_bins
        rows = [("tracks", len(self)), ("chromosomes", len(self.chrom_names)), ("bin size", self.bin_size),
                ("bins", self.n_bins), ("nonzero cells", int(self.M.nnz)),
                ("density", float(self.M.nnz) / max(self.n_bins * max(len(self), 1), 1)),
                ("bins per track min", int(nb.min()) if len(nb) else 0),
                ("bins per track median", float(np.median(nb)) if len(nb) else 0.0),
                ("bins per track max", int(nb.max()) if len(nb) else 0),
                ("metadata columns", ", ".join(self.columns[3:]) or "—")]
        return pd.DataFrame(rows, columns=["", "value"]).set_index("")

    summary = describe

    def __iter__(self):
        return iter(self.track_names)

    def __getitem__(self, key):
        """``atlas['SRX...']`` or ``atlas[i]``: that track's row as a dict."""
        if isinstance(key, str):
            if key not in self.track_names:
                raise KeyError(f"no track named {key!r}")
            key = self.track_names.index(key)
        row = self.to_pandas().iloc[int(key)]
        return {k: (v.item() if hasattr(v, "item") else v) for k, v in row.items()}

    def _repr_html_(self):
        from ._display import table_html
        import pandas as pd
        df = self.to_pandas()
        n = len(df)
        show, gap = (pd.concat([df.head(6), df.tail(3)]), 5) if n > 10 else (df, None)
        return table_html(f"Atlas · {n:,} tracks · {self.n_bins:,} bins of {self.bin_size:,} bp · "
                          f"{self.M.nnz:,} cells", list(df.columns),
                          [list(r) for r in show.itertuples(index=False)], gap_after=gap)

    # ---- dunders -------------------------------------------------------

    def __len__(self) -> int:
        return len(self.track_names)

    def __str__(self) -> str:
        return (f"Atlas(tracks={len(self.track_names)}, "
                f"bins={self.n_bins:,}, bin_size={self.bin_size}, "
                f"nnz={self.M.nnz:,})")
    __repr__ = __str__


def _is_columns(d: dict) -> bool:
    """A dict of interval columns (chrom / start / end) rather than of groups."""
    keys = {str(k).lower() for k in d}
    return bool(keys & {"chrom", "chromosome", "chr", "seqnames", "seqname", "start", "end"})


# ---- vectorized stat ---------------------------------------------------


def _fisher_vec(a, b, c, d, *, alternative: str = "two-sided"):
    """Vectorized Fisher 2x2 over many tables.

    Returns ``(log2_odds_ratio, p, giggle_score)``. ``p`` uses the
    hypergeometric distribution (exact); the GIGGLE score is
    ``-log10(p) * log2(OR)`` -- sign flips naturally for depletion.
    """
    from scipy.stats import hypergeom
    a = np.asarray(a); b = np.asarray(b); c = np.asarray(c); d = np.asarray(d)
    N = a + b + c + d
    K = a + c
    n = a + b
    p_enrich = hypergeom.sf(a - 1, N, K, n)
    p_deplete = hypergeom.cdf(a, N, K, n)
    if alternative == "greater":
        p = p_enrich
    elif alternative == "less":
        p = p_deplete
    else:
        p = np.minimum(np.minimum(p_enrich, p_deplete) * 2.0, 1.0)
    p = np.clip(np.nan_to_num(p, nan=1.0), 1e-300, 1.0)
    log_or = np.log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
    score = -np.log10(p) * log_or
    return log_or, p, score


# ---- range helpers -----------------------------------------------------


def _explode_ranges(starts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Concatenate ``arange(s, s+l)`` for each (s, l) pair, vectorized."""
    if len(starts) == 0:
        return np.zeros(0, dtype=np.int64)
    total = int(lengths.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64)
    out_starts = np.empty(len(starts), dtype=np.int64)
    out_starts[0] = 0
    if len(starts) > 1:
        np.cumsum(lengths[:-1], out=out_starts[1:])
    within = np.arange(total, dtype=np.int64) - np.repeat(out_starts, lengths)
    return np.repeat(starts.astype(np.int64), lengths) + within


# ---- workers ----------------------------------------------------------


def _build_worker(args):
    """Read one BED file, return (track_id, sorted_unique_bin_ids, n_peaks)."""
    import pandas as pd

    tid, path, chrom_offsets, chrom_sizes, bin_size = args
    try:
        df = pd.read_csv(
            path, sep="\t", header=None, comment="#",
            usecols=[0, 1, 2], names=["chrom", "start", "end"],
            on_bad_lines="skip", engine="c",
            dtype={0: str},
        )
    except pd.errors.EmptyDataError:
        return tid, np.zeros(0, dtype=np.int64), 0

    df = df.dropna(subset=["start", "end"])
    if len(df) == 0:
        return tid, np.zeros(0, dtype=np.int64), 0

    df = df[df["chrom"].isin(chrom_offsets)]
    if len(df) == 0:
        return tid, np.zeros(0, dtype=np.int64), 0

    starts = df["start"].to_numpy(dtype=np.int64)
    ends = df["end"].to_numpy(dtype=np.int64)
    n_peaks = int(len(df))

    chrom_max = {
        c: (s + bin_size - 1) // bin_size for c, s in chrom_sizes.items()
    }
    offs = df["chrom"].map(chrom_offsets).to_numpy(dtype=np.int64)
    cb_max = df["chrom"].map(chrom_max).to_numpy(dtype=np.int64)

    s_bins = np.maximum(starts // bin_size, 0)
    e_bins = np.maximum((ends - 1) // bin_size + 1, s_bins + 1)
    s_bins = np.minimum(s_bins, cb_max)
    e_bins = np.minimum(e_bins, cb_max)
    valid = s_bins < e_bins
    if not valid.any():
        return tid, np.zeros(0, dtype=np.int64), n_peaks

    s_abs = s_bins[valid] + offs[valid]
    spans = e_bins[valid] - s_bins[valid]
    return tid, np.unique(_explode_ranges(s_abs, spans)), n_peaks


# ---- file / chromsize helpers ----------------------------------------


def _basename(p: str, pattern: Optional[str] = None) -> str:
    """Extract a track id from a file path.

    Default: strip the BED-family extension. If ``pattern`` is given, it's
    used as ``re.search`` against the basename and the first group (or the
    full match if there are no groups) is returned. Useful for ChIP-Atlas
    files like ``SRX10895570.05.bed.gz`` -- pass ``r'^[^.]+'`` to keep
    only the SRX accession.
    """
    b = os.path.basename(p)
    if pattern is not None:
        import re
        m = re.search(pattern, b)
        if m:
            return m.group(1) if m.groups() else m.group(0)
    for ext in sorted(_BED_EXT, key=len, reverse=True):
        if b.endswith(ext):
            return b[: -len(ext)]
    return b


def _resolve_paths(paths) -> List[str]:
    if isinstance(paths, str):
        if any(c in paths for c in "*?["):
            return sorted(glob.glob(paths))
        if os.path.isdir(paths):
            return sorted(
                os.path.join(paths, f)
                for f in os.listdir(paths)
                if f.endswith(_BED_EXT)
            )
        return [paths]
    return list(paths)
