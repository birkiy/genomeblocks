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
from typing import Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix

from .loci import Loci
from .locus import Locus


_BED_EXT = (".bed", ".bed.gz", ".narrowPeak", ".narrowPeak.gz",
            ".broadPeak", ".broadPeak.gz")


class Atlas:
    """Sparse bin x track index for fast multi-file overlap enrichment."""

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
        M: csr_matrix,
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

        chrom_sizes = _resolve_chromsizes(chromsizes)
        chrom_names = list(chrom_sizes.keys())
        chrom_offsets: Dict[str, int] = {}
        cum = 0
        for c in chrom_names:
            chrom_offsets[c] = cum
            cum += (chrom_sizes[c] + bin_size - 1) // bin_size
        n_bins = cum

        if workers is None:
            workers = max(1, (os.cpu_count() or 2) - 1)

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

    def _intervals_to_bin_ranges(self, loci: Loci) -> np.ndarray:
        """Map each interval to ``[abs_start_bin, abs_end_bin)``.

        Intervals on chroms missing from the atlas are silently dropped.
        Returns an ``(n, 2)`` int64 array; empty if no interval qualifies.
        """
        bs = self.bin_size
        co = self.chrom_offsets
        cmax = self._chrom_max_bins
        out_s: List[int] = []
        out_e: List[int] = []
        for l in loci:
            off = co.get(l.chrom)
            if off is None:
                continue
            cb = cmax[l.chrom]
            s = int(l.start) // bs
            e = (int(l.end) - 1) // bs + 1
            if s < 0:
                s = 0
            if e <= s:
                e = s + 1
            if s > cb:
                s = cb
            if e > cb:
                e = cb
            if s >= e:
                continue
            out_s.append(off + s)
            out_e.append(off + e)
        if not out_s:
            return np.zeros((0, 2), dtype=np.int64)
        return np.column_stack(
            [np.asarray(out_s, dtype=np.int64),
             np.asarray(out_e, dtype=np.int64)],
        )

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
        query: Loci,
        *,
        ref: Optional[Loci] = None,
        alternative: str = "two-sided",
    ):
        """GIGGLE-style enrichment of each track over the query.

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
            other_label, other_val = "n_ref_bins", n_r
        else:
            c = self.track_n_bins - a
            b = n_q - a
            d = self.n_bins - n_q - c
            other_label, other_val = "n_track_bins", None  # filled below

        a, b, c, d = (np.maximum(x, 0).astype(np.int64) for x in (a, b, c, d))
        log_or, p, score = _fisher_vec(a, b, c, d, alternative=alternative)

        n_tracks = len(self.track_names)
        out = {
            "name": self.track_names,
            "n_query_bins": np.full(n_tracks, n_q, dtype=np.int64),
            other_label: (np.full(n_tracks, other_val, dtype=np.int64)
                          if other_val is not None else self.track_n_bins),
            "track_n_bins": self.track_n_bins,
            "track_n_peaks": self.track_n_peaks,
            "overlaps": a,
            "log2_odds": log_or,
            "p": p,
            "giggle_score": score,
        }
        df = pd.DataFrame(out)
        df = self._with_meta(df)
        return df.sort_values("giggle_score", ascending=False).reset_index(drop=True)

    # ---- Monte Carlo bootstrap ----------------------------------------

    def bootstrap(
        self,
        query: Loci,
        *,
        n: int = 10,
        keep_chrom: bool = True,
        seed: Optional[int] = None,
        verbose: bool = True,
    ):
        """Compare observed shared-bin counts to a shuffled-position null.

        ``keep_chrom=True`` (default) preserves each interval's chromosome
        and only randomizes the start bin -- the standard ChIP-seq
        permutation. ``keep_chrom=False`` resamples chroms by length.
        """
        import pandas as pd

        if n <= 0:
            raise ValueError("n must be positive.")
        rng = np.random.default_rng(seed)

        chrom_idx = {c: i for i, c in enumerate(self.chrom_names)}
        chrom_offsets_arr = np.array(
            [self.chrom_offsets[c] for c in self.chrom_names],
            dtype=np.int64,
        )
        chrom_bins = np.array(
            [self._chrom_max_bins[c] for c in self.chrom_names],
            dtype=np.int64,
        )

        bs = self.bin_size
        c_ids: List[int] = []
        lens: List[int] = []
        for l in query:
            ci = chrom_idx.get(l.chrom)
            if ci is None:
                continue
            s_b = int(l.start) // bs
            e_b = max((int(l.end) - 1) // bs + 1, s_b + 1)
            length = min(e_b - s_b, int(chrom_bins[ci]))
            c_ids.append(ci)
            lens.append(length)
        if not c_ids:
            raise ValueError("Empty query (after chrom filter).")
        c_ids_arr = np.asarray(c_ids, dtype=np.int64)
        lens_arr = np.asarray(lens, dtype=np.int64)

        observed = self._overlap_counts(
            self._bins_union(self._intervals_to_bin_ranges(query)),
        )

        chrom_p = chrom_bins / chrom_bins.sum()
        n_intervals = len(c_ids_arr)
        n_tracks = len(self.track_names)
        boot = np.zeros((n, n_tracks), dtype=np.int64)

        it = range(n)
        if verbose:
            from tqdm import tqdm
            it = tqdm(it, desc="[atlas bootstrap]")
        for i in it:
            cs = c_ids_arr if keep_chrom else rng.choice(
                len(self.chrom_names), size=n_intervals, p=chrom_p,
            )
            cb = chrom_bins[cs]
            free = np.maximum(cb - lens_arr, 0)
            local_starts = rng.integers(0, free + 1)
            starts = chrom_offsets_arr[cs] + local_starts
            q_bins = np.unique(_explode_ranges(starts, lens_arr))
            boot[i] = self._overlap_counts(q_bins)

        mean = boot.mean(axis=0)
        std = boot.std(axis=0, ddof=1) if n > 1 else np.zeros(n_tracks)
        z = np.divide(observed - mean, std,
                      out=np.zeros_like(mean, dtype=float), where=std > 0)
        log2fc = np.log2((observed + 1.0) / (mean + 1.0))
        more_extreme = (np.abs(boot - mean[None, :])
                        >= np.abs(observed - mean)[None, :]).sum(axis=0)
        p_emp = (1 + more_extreme) / (n + 1)

        df = pd.DataFrame({
            "name": self.track_names,
            "observed": observed,
            "expected": mean,
            "std": std,
            "log2fc": log2fc,
            "z": z,
            "p_emp": p_emp,
            "track_n_bins": self.track_n_bins,
        })
        df = self._with_meta(df)
        return df.sort_values("z", ascending=False).reset_index(drop=True)

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
            atlas.meta = df.replace("", pd.NA)
        return atlas

    # ---- dunders -------------------------------------------------------

    def __len__(self) -> int:
        return len(self.track_names)

    def __str__(self) -> str:
        return (f"Atlas(tracks={len(self.track_names)}, "
                f"bins={self.n_bins:,}, bin_size={self.bin_size}, "
                f"nnz={self.M.nnz:,})")
    __repr__ = __str__


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


def _resolve_chromsizes(chromsizes) -> Dict[str, int]:
    if isinstance(chromsizes, dict):
        return {k: int(v) for k, v in chromsizes.items()}
    if isinstance(chromsizes, str):
        sizes = {}
        with open(chromsizes) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                sizes[parts[0]] = int(parts[1])
        return sizes
    if hasattr(chromsizes, "chromsizes"):
        return {k: int(v) for k, v in dict(chromsizes.chromsizes).items()}
    raise TypeError(f"Unsupported chromsizes type: {type(chromsizes)}")


# ---- fluent attachments to Loci ---------------------------------------


def _loci_enrich(self: Loci, atlas: Atlas, *,
                 ref: Optional[Loci] = None, **kw):
    return atlas.search(self, ref=ref, **kw)


def _loci_enrich_mc(self: Loci, atlas: Atlas, *, n: int = 10, **kw):
    return atlas.bootstrap(self, n=n, **kw)


Loci.enrich = _loci_enrich
Loci.enrich_mc = _loci_enrich_mc
