"""Pairs: BEDPE as a table, and contact counting from pairs files.

A :class:`Pairs` holds the two anchors of every loop as two row-aligned
:class:`~genomeblocks.Loci` (``P.a``, ``P.b``) plus extra columns (name,
score, ...). This module is the one BEDPE reader in genomeblocks:
``Architecture.make`` and the browser call it under the hood.

``count_pairs`` / ``count_pairs_2d`` stream HiC-Pro allValidPairs, 4DN
``.pairs`` or Juicer files once and count contacts per window.
"""
from __future__ import annotations

from typing import Iterator, Optional, Tuple

import numpy as np

from .genome import Genome
from .loci import SCODE, STRANDS, Loci

_BEDPE = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]


class Pairs:
    """BEDPE rows: anchor ``a`` and anchor ``b`` (two aligned Loci) + columns."""

    def __init__(self, a: Loci, b: Loci, cols=None, *, filename: Optional[str] = None):
        if len(a) != len(b):
            raise ValueError("both anchors need one row per pair")
        if b.genome is not a.genome:
            b = a._check(b)
        self.a, self.b = a, b
        self.cols = {k: np.asarray(v) for k, v in (cols or {}).items()}
        self.filename = filename

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, filename: str, *, genome: Optional[Genome] = None, min_score: Optional[float] = None,
             max_distance: Optional[float] = None, backend: Optional[str] = None) -> "Pairs":
        """Read a BEDPE (``.gz`` too): the six anchor columns, then name, score,
        strand1, strand2 when present; further columns are kept as ``col11`` ..."""
        from .backends.tables import read_columns, sniff
        _, ncol = sniff(filename)
        if ncol and ncol < 6:
            raise ValueError(f"{filename}: a BEDPE needs at least 6 columns, found {ncol}")
        use = list(range(max(ncol, 6)))
        d = read_columns(filename, use, ints=[1, 2, 4, 5], floats=[7] if ncol > 7 else [], backend=backend)
        g = genome if genome is not None else Genome()
        s1 = np.array([SCODE.get(x, 0) for x in d[8]], np.int8) if 8 in d else None
        s2 = np.array([SCODE.get(x, 0) for x in d[9]], np.int8) if 9 in d else None
        a = Loci(g.encode(d[0]), d[1], d[2], s1, genome=g)
        b = Loci(g.encode(d[3]), d[4], d[5], s2, genome=g)
        cols = {}
        if 6 in d:
            cols["name"] = d[6]
        if 7 in d:
            cols["score"] = d[7]
        for i in range(10, ncol):
            cols[f"col{i + 1}"] = d[i]
        P = cls(a, b, cols, filename=str(filename))
        return P.filter(min_score=min_score, max_distance=max_distance)

    @classmethod
    def from_frame(cls, df, *, genome: Optional[Genome] = None) -> "Pairs":
        """From a pandas / polars / pyarrow frame with BEDPE columns (chrom1, start1,
        end1, chrom2, start2, end2 by name, else the first six)."""
        from .interop import _as_pandas
        pdf = _as_pandas(df)
        cols = list(pdf.columns)
        low = {str(c).lower(): c for c in cols}
        need = _BEDPE[:6]
        if all(n in low for n in need):
            names = [low[n] for n in need]
        elif len(cols) >= 6:
            names = cols[:6]
        else:
            raise ValueError(f"a BEDPE frame needs the columns {list(need)} (by name, or as its first six "
                             f"columns); got {cols}")
        g = genome if genome is not None else Genome()

        def strand(k):
            c = low.get(k)
            return None if c is None else np.array([SCODE.get(str(x), 0) for x in pdf[c]], np.int8)
        a = Loci(g.encode(pdf[names[0]].astype(str).to_numpy(object)), pdf[names[1]].to_numpy(np.int64),
                 pdf[names[2]].to_numpy(np.int64), strand("strand1"), genome=g)
        b = Loci(g.encode(pdf[names[3]].astype(str).to_numpy(object)), pdf[names[4]].to_numpy(np.int64),
                 pdf[names[5]].to_numpy(np.int64), strand("strand2"), genome=g)
        used = set(names) | {low.get("strand1"), low.get("strand2")}
        return cls(a, b, {str(k): pdf[k].to_numpy() for k in cols if k not in used})

    def take(self, idx) -> "Pairs":
        return Pairs(self.a.take(idx), self.b.take(idx), {k: v[idx] for k, v in self.cols.items()},
                     filename=self.filename)

    # ── shape ─────────────────────────────────────────────────────────────
    def __len__(self):
        return len(self.a)

    def __getitem__(self, key):
        """``P[i]`` the two anchors of pair i; ``P['score']`` a column; else a selection."""
        if isinstance(key, (int, np.integer)):
            return self.a[key], self.b[key]
        if isinstance(key, str):
            return self.cols[key]
        return self.take(key)

    def __iter__(self):
        for i in range(len(self)):
            yield self.a[i], self.b[i]

    @property
    def genome(self) -> Genome:
        return self.a.genome

    @property
    def mids1(self) -> np.ndarray:
        return self.a.centers

    @property
    def mids2(self) -> np.ndarray:
        return self.b.centers

    @property
    def is_cis(self) -> np.ndarray:
        return self.a.codes == self.b.codes

    @property
    def distance(self) -> np.ndarray:
        """|mid2 - mid1| for cis pairs, inf for trans."""
        d = np.abs(self.mids2 - self.mids1).astype(float)
        d[~self.is_cis] = np.inf
        return d

    # ── selections ────────────────────────────────────────────────────────
    def filter(self, *, min_score: Optional[float] = None, max_distance: Optional[float] = None) -> "Pairs":
        keep = np.ones(len(self), bool)
        if min_score is not None and "score" in self.cols:
            keep &= np.asarray(self.cols["score"], float) >= min_score
        if max_distance is not None:
            keep &= self.distance <= max_distance
        return self if keep.all() else self.take(keep)

    def anchors_overlap(self, loci, *, r: int = 0, backend: Optional[str] = None):
        """(mask anchor a overlaps ``loci``, mask anchor b overlaps ``loci``), anchors widened by ``r``."""
        L = self.a._check(loci)
        a, b = (self.a.slop(r), self.b.slop(r)) if r else (self.a, self.b)
        return a.overlap_any(L, backend=backend), b.overlap_any(L, backend=backend)

    def overlapping(self, loci, *, r: int = 0, both: bool = False, backend: Optional[str] = None) -> "Pairs":
        """Pairs with one anchor (or ``both``) overlapping ``loci`` — bedtools pairtobed."""
        m1, m2 = self.anchors_overlap(loci, r=r, backend=backend)
        return self.take((m1 & m2) if both else (m1 | m2))

    # ── export ────────────────────────────────────────────────────────────
    def to_pandas(self):
        import pandas as pd
        d = {"chrom1": self.a.chroms.astype(str), "start1": self.a.starts, "end1": self.a.ends,
             "chrom2": self.b.chroms.astype(str), "start2": self.b.starts, "end2": self.b.ends}
        d.update({k: v for k, v in self.cols.items() if k in ("name", "score")})
        d["strand1"], d["strand2"] = STRANDS[self.a.strands], STRANDS[self.b.strands]
        d.update({k: v for k, v in self.cols.items() if k not in ("name", "score")})
        return pd.DataFrame(d)

    def to_polars(self):
        import polars as pl
        return pl.from_pandas(self.to_pandas())

    def to_arrow(self):
        import pyarrow as pa
        return pa.Table.from_pandas(self.to_pandas(), preserve_index=False)

    def to_bedpe(self, path: str) -> str:
        df = self.to_pandas()
        if "name" not in df:
            df.insert(6, "name", ".")
        if "score" not in df:
            df.insert(7, "score", 0)
        df.to_csv(path, sep="\t", header=False, index=False)
        return path

    def save(self, path: str):
        import pyarrow.parquet as pq
        pq.write_table(self.to_arrow(), path)

    @classmethod
    def load(cls, path: str, *, genome: Optional[Genome] = None) -> "Pairs":
        import pyarrow.parquet as pq
        return cls.from_frame(pq.read_table(path), genome=genome)

    def __repr__(self):
        n = len(self)
        cis = int(self.is_cis.sum()) if n else 0
        return f"Pairs(n={n:,}, cis={cis:,}, trans={n - cis:,}{', cols=[' + ', '.join(self.cols) + ']' if self.cols else ''})"


def read_bedpe(filename: str, **kw) -> Pairs:
    """A BEDPE file as :class:`Pairs` (same as ``Pairs.make``)."""
    return Pairs.make(filename, **kw)


def as_pairs(x, *, genome: Optional[Genome] = None) -> Pairs:
    """Pairs from a BEDPE path, a Pairs, or a frame with BEDPE columns."""
    if isinstance(x, Pairs):
        return x
    if isinstance(x, str) or hasattr(x, "__fspath__"):
        p = str(x)
        return Pairs.load(p, genome=genome) if p.endswith(".parquet") else Pairs.make(p, genome=genome)
    return Pairs.from_frame(x, genome=genome)


# ── pairs files: contacts per window ───────────────────────────────────────

# (chrom1_col, pos1_col, chrom2_col, pos2_col) — 0-based columns
PAIRS_FORMAT_COLUMNS = {
    "allvalidpairs": (1, 2, 4, 5),  # HiC-Pro: readID chr1 pos1 strand1 chr2 pos2 strand2 ...
    "pairs":         (1, 2, 3, 4),  # pairtools / 4DN: readID chrom1 pos1 chrom2 pos2 strand1 strand2 ...
    "juicer":        (2, 3, 6, 7),  # Juicer medium: readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2
}


def _detect_pairs_format(filename: str) -> str:
    import gzip
    op = gzip.open if str(filename).endswith(".gz") else open
    with op(filename, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            num = lambda s: s.lstrip("-").isdigit()                      # noqa: E731
            # Juicer medium: readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 ... with the
            # strands as +/- or as SAM flags (0 forward, 16 reverse)
            if (len(fields) >= 8 and fields[1] in ("+", "-", "0", "16") and fields[5] in ("+", "-", "0", "16")
                    and num(fields[3]) and num(fields[7]) and not num(fields[2]) and not num(fields[6])):
                return "juicer"
            if len(fields) >= 7 and fields[3] in ("+", "-"):
                return "allvalidpairs"
            return "pairs"
    raise ValueError(f"No data lines in {filename}")


def read_pairs_chunks(filename: str, *, format: str = "auto", columns: Optional[Tuple[int, int, int, int]] = None,
                      chunksize: int = 2_000_000) -> Iterator:
    """Stream a pairs file as DataFrames of chrom1, pos1, chrom2, pos2 (categorical chroms)."""
    import pandas as pd
    if columns is not None:
        usecols = list(columns)
    else:
        if format == "auto":
            format = _detect_pairs_format(filename)
        if format not in PAIRS_FORMAT_COLUMNS:
            raise ValueError(f"Unknown pairs format: {format!r}. Known: {sorted(PAIRS_FORMAT_COLUMNS)}, "
                             f"or pass columns=(c1, p1, c2, p2).")
        usecols = list(PAIRS_FORMAT_COLUMNS[format])
    order = sorted(range(4), key=lambda k: usecols[k])
    names = ["chrom1", "pos1", "chrom2", "pos2"]
    reader = pd.read_csv(filename, sep="\t", header=None, comment="#", usecols=usecols,
                         dtype={usecols[0]: "category", usecols[1]: "int64",
                                usecols[2]: "category", usecols[3]: "int64"},
                         chunksize=chunksize, engine="c")
    for ch in reader:
        ch.columns = [names[k] for k in order]
        yield ch[names]


class _Windows:
    """Windows sorted on the genome axis for vectorised position lookup."""

    def __init__(self, loci: Loci):
        self.loci = loci
        g = loci.genome
        self.key = loci.codes.astype(np.int64) * (1 << 40) + loci.starts
        self.order = np.argsort(self.key, kind="stable")
        self.key, self.ends = self.key[self.order], loci.ends[self.order]
        self.codes = loci.codes[self.order]
        self.genome = g

    def locate(self, chroms, pos) -> np.ndarray:
        """Row of the window holding each (chrom, pos); -1 if none."""
        chroms = chroms.astype("category")
        lut = np.array([self.genome.code.get(str(c), -1) for c in chroms.cat.categories] + [-1], np.int64)
        c = lut[chroms.cat.codes.to_numpy()]
        i = np.searchsorted(self.key, c * (1 << 40) + pos, side="right") - 1
        j = np.maximum(i, 0)
        ok = (c >= 0) & (i >= 0) & (self.codes[j] == c) & (pos < self.ends[j])
        return np.where(ok, self.order[j], -1)


def count_pairs(loci, pairs_file: str, *, target_chrom: Optional[str] = None, format: str = "auto",
                columns: Optional[Tuple[int, int, int, int]] = None, chunksize: int = 2_000_000,
                verbose: bool = True):
    """Contacts landing in each window of ``loci``, split by partner chromosome.

    One streaming pass. For each pair, each end inside a window counts once
    for that window under the other end's chromosome (a cis pair with both
    ends in windows counts once per end). Windows on a chromosome must not
    overlap. Returns a DataFrame aligned to the rows of ``loci``: chrom, start,
    end, uid, then 'count' (with ``target_chrom``) or one column per partner
    chromosome with at least one contact.
    """
    import pandas as pd
    from .interop import as_loci
    loci = as_loci(loci)
    W = _Windows(loci)
    n = len(loci)
    names: dict = {}
    mat = np.zeros((n, 0), np.int64)
    target = np.zeros(n, np.int64)
    n_seen = 0
    for chunk in read_pairs_chunks(pairs_file, format=format, columns=columns, chunksize=chunksize):
        n_seen += len(chunk)
        for self_c, self_p, other_c in (("chrom1", "pos1", "chrom2"), ("chrom2", "pos2", "chrom1")):
            w = W.locate(chunk[self_c], chunk[self_p].to_numpy(np.int64))
            other = chunk[other_c].astype("category")
            other_names = [str(c) for c in other.cat.categories] + ["nan"]
            if target_chrom is not None:
                hit = (w >= 0) & np.array([c == target_chrom for c in other_names])[other.cat.codes.to_numpy()]
                target += np.bincount(w[hit], minlength=n)
                continue
            p_lut = np.array([names.setdefault(c, len(names)) for c in other_names], np.int64)
            if len(names) > mat.shape[1]:
                mat = np.pad(mat, ((0, 0), (0, len(names) - mat.shape[1])))
            hit = w >= 0
            p = p_lut[other.cat.codes.to_numpy()][hit]
            P = mat.shape[1]
            mat += np.bincount(w[hit] * P + p, minlength=n * P).reshape(n, P)
        if verbose:
            print(f"[INFO] processed {n_seen:,} pairs", end="\r")
    if verbose:
        print(f"\n[INFO] processed {n_seen:,} pairs total")
    out = pd.DataFrame({"chrom": loci.chroms.astype(str), "start": loci.starts, "end": loci.ends,
                        "uid": loci.uid})
    if target_chrom is not None:
        out["count"] = target
    else:
        for chrom in sorted(names):
            col = mat[:, names[chrom]]
            if col.any():
                out[chrom] = col
    return out


def count_pairs_2d(loci_a, pairs_file: str, *, loci_b=None, format: str = "auto",
                   columns: Optional[Tuple[int, int, int, int]] = None, chunksize: int = 2_000_000,
                   verbose: bool = True):
    """Window-to-window contact counts as a ``scipy.sparse.csr_matrix``
    (rows of ``loci_a`` x rows of ``loci_b``; ``loci_b=None`` uses ``loci_a``,
    which makes the matrix symmetric). One streaming pass; both orientations
    of every pair are counted, so a pair inside one window adds 2 to its
    diagonal cell in the symmetric case."""
    import scipy.sparse as sp
    from .interop import as_loci
    loci_a = as_loci(loci_a)
    same = loci_b is None
    loci_b = loci_a if same else loci_a._check(loci_b)
    Wa = _Windows(loci_a)
    Wb = Wa if same else _Windows(loci_b)
    n_a, n_b = len(loci_a), len(loci_b)
    rows_acc, cols_acc, data_acc = [], [], []
    n_seen = 0
    for chunk in read_pairs_chunks(pairs_file, format=format, columns=columns, chunksize=chunksize):
        n_seen += len(chunk)
        p1 = chunk["pos1"].to_numpy(np.int64)
        p2 = chunk["pos2"].to_numpy(np.int64)
        a1, a2 = Wa.locate(chunk["chrom1"], p1), Wa.locate(chunk["chrom2"], p2)
        b1, b2 = (a1, a2) if same else (Wb.locate(chunk["chrom1"], p1), Wb.locate(chunk["chrom2"], p2))
        keys = []
        for ga, gb in ((a1, b2), (a2, b1)):
            ok = (ga >= 0) & (gb >= 0)
            keys.append(ga[ok] * n_b + gb[ok])
        keys = np.concatenate(keys)
        if keys.size:
            uk, uc = np.unique(keys, return_counts=True)
            rows_acc.append(uk // n_b)
            cols_acc.append(uk % n_b)
            data_acc.append(uc.astype(np.int64))
        if verbose:
            print(f"[INFO] processed {n_seen:,} pairs", end="\r")
    if verbose:
        print(f"\n[INFO] processed {n_seen:,} pairs total")
    if not data_acc:
        return sp.csr_matrix((n_a, n_b), dtype=np.int64)
    coo = sp.coo_matrix((np.concatenate(data_acc), (np.concatenate(rows_acc), np.concatenate(cols_acc))),
                        shape=(n_a, n_b), dtype=np.int64)
    coo.sum_duplicates()
    return coo.tocsr()


def pair_2d_block(mat, loci_a, loci_b, chrom_a: str, chrom_b: str):
    """Dense block of a ``count_pairs_2d`` matrix for one chromosome pair:
    (block, windows of loci_a on chrom_a, windows of loci_b on chrom_b)."""
    ra = np.flatnonzero(loci_a.codes == loci_a.genome.code.get(chrom_a, -1))
    cb = np.flatnonzero(loci_b.codes == loci_b.genome.code.get(chrom_b, -1))
    if not len(ra) or not len(cb):
        return np.zeros((len(ra), len(cb)), np.int64), loci_a.take(ra), loci_b.take(cb)
    return mat[ra, :][:, cb].toarray(), loci_a.take(ra), loci_b.take(cb)


def pair_2d_to_frame(mat, loci_a, loci_b):
    """Non-zero cells of a ``count_pairs_2d`` matrix as a long DataFrame
    (chrom1 start1 end1 chrom2 start2 end2 count)."""
    import pandas as pd
    coo = mat.tocoo()
    r, c = coo.row.astype(np.int64), coo.col.astype(np.int64)
    return pd.DataFrame({"chrom1": loci_a.chroms[r].astype(str), "start1": loci_a.starts[r],
                         "end1": loci_a.ends[r], "chrom2": loci_b.chroms[c].astype(str),
                         "start2": loci_b.starts[c], "end2": loci_b.ends[c], "count": coo.data})
