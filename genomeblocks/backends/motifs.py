"""Motif backends: one motif library, three scanning engines.

A :class:`Library` holds motifs as count matrices (W x 4, A C G T) whatever
they were read from — a JASPAR / jaspar16 / TRANSFAC / uniprobe / MEME file
(parsed here, in pure Python), Biopython ``Bio.motifs`` objects, lightmotif
motifs, or plain arrays. Each motif is scored with one log-odds matrix,

    log2( (count + 0.1) / (column total + 0.4) / 0.25 )

computed bit for bit the way lightmotif computes it (float32, its summation
order, the platform ``log2f``), and every engine scans that same matrix, so
MOODS (C++, the default), lightmotif (SIMD) and Biopython (numpy) report
the same hits. A hit is a window position whose score is at least the
threshold; windows containing N never hit. :func:`threshold_from_pvalue`
turns a p-value into a per-motif score cutoff, the same whichever engine
scans.

Sequences are scanned as one block of windows: each motif is scanned once
over the concatenation and the hits are split back per window (a hit that
would straddle two windows is dropped), so per-window results equal scanning
each window on its own.
"""
from __future__ import annotations

import re
from typing import List, Optional, Sequence

import numpy as np

from . import resolve
from .._table import TableMixin

PSEUDOCOUNT = 0.1
_ACGT_FROM_LM = [0, 1, 3, 2]                    # lightmotif columns are A C T G N
_BAD = re.compile(r"[^ACGTN]")
_LOG2F = None


def _log2f(x: np.ndarray) -> np.ndarray:
    """float32 log2 through the platform C library — what lightmotif's Rust
    ``f32::log2`` calls — so the matrices match it bit for bit; a correctly
    rounded float64 log2 when no C library can be loaded."""
    global _LOG2F
    if _LOG2F is None:
        try:
            import ctypes
            import ctypes.util
            libm = ctypes.CDLL(ctypes.util.find_library("m") or "libm.so.6")
            f = libm.log2f
            f.restype, f.argtypes = ctypes.c_float, [ctypes.c_float]
            _LOG2F = np.frompyfunc(lambda v: f(v), 1, 1)
        except (OSError, AttributeError, TypeError):
            _LOG2F = False
    x = np.asarray(x, np.float32)
    if _LOG2F is False:
        return np.log2(x.astype(np.float64)).astype(np.float32)
    return np.asarray(_LOG2F(x), dtype=np.float32)


def logodds_matrix(counts, pseudocount: float = PSEUDOCOUNT) -> np.ndarray:
    """(W x 4) log2-odds of a (W x 4) A C G T count matrix against a uniform
    background — lightmotif's ``normalize(p).log_odds()``, reproduced exactly
    (float32 arithmetic summed in its A C T G order), returned as float64."""
    c = np.asarray(counts, np.float32)
    p = np.float32(pseudocount)
    a, cc, g, t = c[:, 0] + p, c[:, 1] + p, c[:, 2] + p, c[:, 3] + p
    total = ((a + cc) + t) + g
    freq = np.stack([a / total, cc / total, g / total, t / total], 1).astype(np.float32)
    return _log2f(freq / np.float32(0.25)).astype(np.float64)


def _dp_tail(q: np.ndarray, floor: int = 0):
    """Distribution of the shifted integer score S = sum_j (q[j, b_j] - min_j) under
    a uniform background, kept only from ``floor`` up (lower partial sums that can
    no longer reach ``floor`` are dropped). Returns (P(S = floor + i))_i."""
    mins = q.min(1)
    d = q - mins[:, None]
    rng = d.max(1)
    rem = np.concatenate([np.cumsum(rng[::-1])[::-1][1:], [0]])       # most the later columns can add
    lo, dist = 0, np.ones(1)
    for j in range(len(d)):
        hi = lo + len(dist) - 1 + int(rng[j])
        new_lo = max(lo, floor - int(rem[j]))
        new = np.zeros(hi - new_lo + 1)
        for b in range(4):
            a = lo + int(d[j, b]) - new_lo                             # where dist[0] lands
            s0 = max(0, -a)
            if s0 < len(dist):
                new[a + s0:a + len(dist)] += 0.25 * dist[s0:]
        lo, dist = new_lo, new
    return lo, dist


def threshold_from_pvalue(lo_matrix, pvalue: float, *, resolution: float = 1e-3) -> float:
    """The smallest score ``s`` with P(score >= s) <= ``pvalue`` for one (W x 4)
    log-odds matrix, under a uniform background — exact on a ``resolution`` grid
    (agrees with MOODS' ``threshold_from_p`` to a few thousandths of a bit).

    A coarse pass bounds the answer from below, so the exact pass only builds
    the upper tail of the score distribution. A motif too short to reach
    ``pvalue`` gets its best possible score (only perfect matches count)."""
    m = np.asarray(lo_matrix, np.float64)
    q = np.rint(m / resolution).astype(np.int64)
    base = int(q.min(1).sum())
    coarse = max(resolution * 20, 0.02)
    qc = np.rint(m / coarse).astype(np.int64)
    lo_c, dist_c = _dp_tail(qc)
    tail_c = np.cumsum(dist_c[::-1])[::-1]
    kc = np.flatnonzero(tail_c <= pvalue)
    if len(kc):
        t_c = (lo_c + kc[0] + int(qc.min(1).sum())) * coarse
        floor = int(np.floor((t_c - len(m) * coarse - coarse) / resolution)) - base - 1
        floor = max(floor, 0)
    else:
        floor = 0
    lo, dist = _dp_tail(q, floor)
    tail = np.cumsum(dist[::-1])[::-1]
    k = np.flatnonzero(tail <= pvalue)
    if floor and (not len(k) or k[0] == 0):                            # bound not below the answer: no pruning
        lo, dist = _dp_tail(q, 0)
        tail = np.cumsum(dist[::-1])[::-1]
        k = np.flatnonzero(tail <= pvalue)
    if not len(k):
        return float(m.max(1).sum())
    return float((lo + k[0] + base) * resolution)


# ── the library ────────────────────────────────────────────────────────────

class Library(TableMixin):
    """Motifs as (W x 4) count matrices plus their names.

    As a table (name, description, width, consensus): ``len``, ``shape``,
    ``columns``, ``head()``, ``describe()``, ``to_pandas()`` and the Arrow /
    dataframe protocols; ``lib['CTCF']`` / ``lib[0]`` is a count matrix."""

    def __init__(self, names: Sequence[str], counts: Sequence[np.ndarray],
                 descriptions: Optional[Sequence[str]] = None, *, pseudocount: float = PSEUDOCOUNT,
                 _logodds: Optional[List[np.ndarray]] = None):
        self.names = list(names)
        self.counts = [np.asarray(c, np.float64) for c in counts]
        self.descriptions = list(descriptions) if descriptions is not None else [""] * len(self.names)
        self.pseudocount = pseudocount
        self._lo = _logodds

    def __len__(self):
        return len(self.names)

    def __repr__(self):
        head = ", ".join(self.names[:4]) + (", ..." if len(self) > 4 else "")
        return f"Library({len(self)} motifs: {head})"

    @property
    def widths(self) -> np.ndarray:
        return np.array([len(c) for c in self.counts], np.int64)

    # ── like a DataFrame ─────────────────────────────────────────────────
    @property
    def columns(self) -> list:
        return ["name", "description", "width", "consensus"]

    def consensus(self, i: int) -> str:
        """The most frequent base at each position of motif ``i``."""
        return "".join("ACGT"[k] for k in np.asarray(self.counts[i]).argmax(1))

    def to_pandas(self):
        import pandas as pd
        return pd.DataFrame({"name": self.names, "description": self.descriptions, "width": self.widths,
                             "consensus": [self.consensus(i) for i in range(len(self))]})

    def to_arrow(self):
        import pyarrow as pa
        return pa.Table.from_pandas(self.to_pandas(), preserve_index=False)

    def head(self, n: int = 5):
        return self.to_pandas().head(n)

    def tail(self, n: int = 5):
        return self.to_pandas().tail(n)

    def describe(self):
        import pandas as pd
        w = self.widths
        rows = [("motifs", len(self)), ("width min", int(w.min()) if len(w) else 0),
                ("width median", float(np.median(w)) if len(w) else 0.0),
                ("width max", int(w.max()) if len(w) else 0), ("pseudocount", self.pseudocount),
                ("log-odds", "log2((count + p) / (total + 4p) / 0.25)")]
        return pd.DataFrame(rows, columns=["", "value"]).set_index("").rename_axis(None)

    summary = describe

    def __iter__(self):
        return iter(self.names)

    def __getitem__(self, key):
        """``lib[i]`` or ``lib['name']``: the (W x 4) count matrix of that motif."""
        if isinstance(key, (int, np.integer)):
            return self.counts[int(key)]
        if isinstance(key, str):
            if key in self.names:
                return self.counts[self.names.index(key)]
            hits = self.indices(key, "substring")
            if len(hits) == 1:
                return self.counts[hits[0]]
            raise KeyError(f"no motif named {key!r}" + (f" ({len(hits)} match the substring: use select())"
                                                        if hits else ""))
        raise TypeError("index a Library with a motif name or a position")

    def _repr_html_(self):
        from .._display import table_html
        df = self.to_pandas()
        n = len(df)
        import pandas as pd
        show, gap = (pd.concat([df.head(6), df.tail(3)]), 5) if n > 10 else (df, None)
        return table_html(f"Library · {n:,} motifs", list(df.columns),
                          [list(r) for r in show.itertuples(index=False)], gap_after=gap,
                          note=f"pseudocount {self.pseudocount} · lib[name] gives the count matrix")

    def logodds(self, i: int) -> np.ndarray:
        """(W x 4) log2-odds matrix of motif ``i`` (A C G T columns)."""
        if self._lo is None:
            self._lo = [None] * len(self)
        if self._lo[i] is None:
            self._lo[i] = logodds_matrix(self.counts[i], self.pseudocount)
        return self._lo[i]

    def pfm(self, i: int, pseudo: float = 0.01) -> np.ndarray:
        """(4 x W) probability matrix (rows A C G T)."""
        c = self.counts[i] + pseudo
        return (c / c.sum(1, keepdims=True)).T

    def select(self, motifs=None, match: str = "substring") -> "Library":
        """The motifs named in ``motifs`` (``None`` keeps all). ``match='substring'``
        is case-insensitive and also searches the descriptions (TF symbols in
        JASPAR headers); ``'exact'`` compares names verbatim."""
        idx = self.indices(motifs, match)
        return self.take(idx)

    def indices(self, motifs=None, match: str = "substring") -> List[int]:
        if motifs is None:
            return list(range(len(self)))
        if isinstance(motifs, str):
            motifs = [motifs]
        if match == "exact":
            wanted = set(motifs)
            return [i for i, n in enumerate(self.names) if n in wanted]
        if match == "substring":
            low = [m.lower() for m in motifs]
            return [i for i, (n, d) in enumerate(zip(self.names, self.descriptions))
                    if any(x in n.lower() or (d and x in d.lower()) for x in low)]
        raise ValueError(f"match must be 'exact' or 'substring', got {match!r}")

    def take(self, idx) -> "Library":
        lo = None if self._lo is None else [self._lo[i] for i in idx]
        return Library([self.names[i] for i in idx], [self.counts[i] for i in idx],
                       [self.descriptions[i] for i in idx], pseudocount=self.pseudocount, _logodds=lo)

    # ── export ──────────────────────────────────────────────────────────
    def to_biopython(self):
        """A list of ``Bio.motifs.Motif`` (counts), named like the library."""
        from Bio import motifs as bm
        out = []
        for name, c, d in zip(self.names, self.counts, self.descriptions):
            m = bm.Motif(alphabet="ACGT", counts={b: c[:, k].tolist() for k, b in enumerate("ACGT")})
            m.name, m.matrix_id = d or name, name
            out.append(m)
        return out

    def to_moods(self):
        """MOODS matrices: one 4 x W log-odds list of lists per motif (A C G T rows)."""
        return [self.logodds(i).T.tolist() for i in range(len(self))]

    def to_meme(self, path: str, pseudo: float = 0.01):
        """Write the library as MEME minimal format (probabilities)."""
        write_meme({n: self.pfm(i, pseudo) for i, n in enumerate(self.names)}, path)

    def to_jaspar(self, path: str):
        """Write the library as JASPAR (bracketed, jaspar16) counts."""
        with open(path, "w") as f:
            for n, c, d in zip(self.names, self.counts, self.descriptions):
                f.write(f">{n} {d}".rstrip() + "\n")
                for k, b in enumerate("ACGT"):
                    f.write(f"{b}  [ " + " ".join(f"{v:g}" for v in c[:, k]) + " ]\n")


def write_meme(pfms: dict, path: str, *, alphabet: str = "ACGT", bg=(0.25, 0.25, 0.25, 0.25)):
    """Write ``{name: (4 x W) probability matrix}`` as MEME minimal format."""
    bg_str = " ".join(f"{a} {p:.4f}" for a, p in zip(alphabet, bg))
    with open(path, "w") as f:
        f.write(f"MEME version 4\n\nALPHABET= {alphabet}\n\nstrands: + -\n\n"
                f"Background letter frequencies\n{bg_str}\n\n")
        for name, pfm in pfms.items():
            pfm = np.asarray(pfm, float)
            f.write(f"MOTIF {name}\nletter-probability matrix: alength= 4 w= {pfm.shape[1]}\n")
            for j in range(pfm.shape[1]):
                f.write(" ".join(f"{pfm[i, j]:.6f}" for i in range(4)) + "\n")
            f.write("\n")


def _read_meme(path, nsites: float = 100.0):
    """MEME (minimal or full) -> (names, counts, descriptions). Probabilities are
    scaled by the motif's ``nsites`` (100 when absent) to give counts."""
    names, counts, descs = [], [], []
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MOTIF"):
            parts = line.split()
            name = parts[1] if len(parts) > 1 else f"motif{len(names) + 1}"
            desc = " ".join(parts[2:])
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("letter-probability"):
                i += 1
            if i >= len(lines):
                break
            hdr = lines[i]
            w = int(re.search(r"w=\s*(\d+)", hdr).group(1))
            m = re.search(r"nsites=\s*([\d.eE+-]+)", hdr)
            n = float(m.group(1)) if m else nsites
            rows = []
            i += 1
            while len(rows) < w and i < len(lines):
                vals = lines[i].split()
                if vals:
                    rows.append([float(v) for v in vals[:4]])
                i += 1
            names.append(name)
            counts.append(np.array(rows) * n)
            descs.append(desc)
        else:
            i += 1
    return names, counts, descs


def load_motifs(src, format: str = "jaspar") -> Library:
    """A :class:`Library` from a motif file or motif objects.

    ``src`` may be a path (``format`` = 'jaspar' (raw counts), 'jaspar16'
    (bracketed), 'transfac', 'uniprobe' or 'meme'; ``.gz`` too), a
    ``Library``, one or many Biopython ``Bio.motifs.Motif``, lightmotif
    motifs, or a ``{name: matrix}`` dict of counts / probabilities (4 x W or
    W x 4, rows or columns in A C G T order). Probability matrices (columns
    summing to 1) are scaled to 100 sites.
    """
    if isinstance(src, Library):
        return src
    if isinstance(src, str) or hasattr(src, "__fspath__"):
        import os
        path = str(src)
        if not os.path.exists(path):
            raise FileNotFoundError(f"no such motif file: {path!r}")
        if format not in FORMATS:
            raise ValueError(f"unknown motif format {format!r}; use one of {', '.join(FORMATS)}")
        if format == "meme":
            return Library(*_read_meme(path))
        reader = {"jaspar": _read_jaspar, "jaspar16": _read_jaspar, "transfac": _read_transfac,
                  "uniprobe": _read_uniprobe}[format]
        names, counts, descs = reader(path)
        if not names:
            raise ValueError(f"{path}: no {format} motifs found; pass format= for another layout "
                             f"('jaspar', 'jaspar16', 'transfac', 'uniprobe', 'meme')")
        return Library(names, counts, descs)
    if isinstance(src, np.ndarray):
        src = [src]
    if isinstance(src, (list, tuple)) and src and all(isinstance(a, (np.ndarray, list, tuple)) for a in src):
        src = {f"motif{k + 1}": a for k, a in enumerate(src)}       # bare matrices
    if isinstance(src, dict):
        names, counts = [], []
        for k, v in src.items():
            a = np.asarray(v, np.float64)
            if a.shape[0] == 4 and a.shape[1] != 4:
                a = a.T
            if np.allclose(a.sum(1), 1.0):
                a = a * 100.0
            names.append(str(k))
            counts.append(a)
        return Library(names, counts)
    items = list(src) if not hasattr(src, "counts") else [src]
    names, counts, descs = [], [], []
    for m in items:
        c = m.counts
        if hasattr(c, "keys"):                                       # Biopython
            a = np.array([c[b] for b in "ACGT"], np.float64).T
            names.append(getattr(m, "matrix_id", None) or getattr(m, "name", None) or f"motif{len(names) + 1}")
            descs.append(getattr(m, "name", "") or "")
        else:                                                         # lightmotif
            a = np.asarray(c, np.float64)[:, _ACGT_FROM_LM]
            names.append(_lm_name(m, len(names)))
            descs.append(getattr(m, "description", "") or "")
        counts.append(a)
    return Library(names, counts, descs)


FORMATS = ("jaspar", "jaspar16", "transfac", "uniprobe", "meme")


def _open_text(path):
    import gzip
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _num(tok: str, path, what) -> float:
    try:
        return float(tok)
    except ValueError:
        raise ValueError(f"{path}: {what}: {tok!r} is not a number") from None


def _read_jaspar(path):
    """JASPAR, raw (four rows of counts) or bracketed (``A [ ... ]``), one or many
    motifs: ``>ID description`` then the A C G T rows. Names and descriptions
    follow lightmotif: the ID, then the rest of the header line."""
    names, counts, descs = [], [], []
    name, desc, rows = None, "", []

    def flush():
        if name is None:
            return
        if len(rows) != 4:
            raise ValueError(f"{path}: motif {name!r} has {len(rows)} rows, expected A C G T")
        by = {}
        for k, (lab, vals) in enumerate(rows):
            by[lab or "ACGT"[k]] = vals
        if set(by) != set("ACGT") or len({len(v) for v in by.values()}) != 1:
            raise ValueError(f"{path}: motif {name!r} needs four equal-length rows A C G T")
        names.append(name)
        descs.append(desc)
        counts.append(np.array([by[b] for b in "ACGT"], np.float64).T)

    with _open_text(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                head = line[1:].strip().split(None, 1)
                name = head[0] if head else f"motif{len(names) + 1}"
                desc = head[1].strip() if len(head) > 1 else ""
                rows = []
                continue
            if name is None:
                continue
            lab = None
            if line[0].upper() in "ACGT" and not line[0].isdigit() and (len(line) == 1 or not line[1].isalnum()):
                lab, line = line[0].upper(), line[1:]
            line = line.replace("[", " ").replace("]", " ")
            rows.append((lab, [_num(t, path, f"motif {name!r}") for t in line.split()]))
    flush()
    return names, counts, descs


def _read_transfac(path):
    """TRANSFAC matrices (``P0`` / ``PO`` header, numbered rows, ``//`` between
    motifs). Name: NA, else ID, else AC (as lightmotif); description: DE.
    Rows of frequencies are scaled to 100 sites."""
    names, counts, descs = [], [], []
    cur = {}

    def flush():
        if cur.get("rows"):
            order = cur.get("order") or list("ACGT")
            rows = np.array(cur["rows"], np.float64)
            m = np.zeros((len(rows), 4))
            for k, b in enumerate(order[:rows.shape[1]]):
                if b in "ACGT":
                    m[:, "ACGT".index(b)] = rows[:, k]
            if np.any(m % 1):                                     # frequencies, not counts
                m = m / np.maximum(m.sum(1, keepdims=True), 1e-12) * 100.0
            names.append(cur.get("NA") or cur.get("ID") or cur.get("AC") or f"motif{len(names) + 1}")
            descs.append(cur.get("DE", ""))
            counts.append(m)
        cur.clear()

    with _open_text(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("//"):
                flush()
                continue
            key, rest = line[:2], line[2:].strip()
            if key in ("AC", "ID", "NA", "DE") and key not in cur:
                cur[key] = rest
            elif key in ("P0", "PO"):
                cur["order"] = [x.upper() for x in rest.split()]
                cur["rows"] = []
            elif "rows" in cur and line.split()[0].isdigit():
                vals = line.split()[1:]
                n = len(cur.get("order") or "ACGT")
                cur["rows"].append([_num(v, path, "matrix row") for v in vals[:n]])
    flush()
    return names, counts, descs


def _read_uniprobe(path):
    """UniPROBE: an optional name line, then ``A:`` ``C:`` ``G:`` ``T:`` rows of
    probabilities (one or many motifs, blank lines between). Scaled to 100 sites."""
    names, counts, descs = [], [], []
    name, rows = None, {}

    def flush():
        if len(rows) == 4:
            m = np.array([rows[b] for b in "ACGT"], np.float64).T
            m = m / np.maximum(m.sum(1, keepdims=True), 1e-12) * 100.0
            names.append(name or f"motif{len(names) + 1}")
            descs.append("")
            counts.append(m)
        elif rows:
            raise ValueError(f"{path}: motif {name!r} needs rows A: C: G: T:")

    with _open_text(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            lab = line[:2].upper()
            if len(lab) == 2 and lab[1] == ":" and lab[0] in "ACGT":
                rows[lab[0]] = [_num(v, path, "uniprobe row") for v in line[2:].split()]
            else:
                flush()
                name, rows = line, {}
    flush()
    return names, counts, descs


def _lm_name(m, k: int) -> str:
    """A lightmotif motif's name: NA, else ID, else AC, else motif<k+1>."""
    for attr in ("name", "id", "accession"):
        v = getattr(m, attr, None)
        if v:
            return str(v)
    return f"motif{k + 1}"


# ── scanning ───────────────────────────────────────────────────────────────

def _rc(mat):
    """Reverse complement of a (W x 4) A C G T matrix."""
    return mat[::-1, ::-1]


class Block:
    """Windows (strings) laid end to end once, scanned per motif by one engine."""

    def __init__(self, seqs: Sequence[str], backend: Optional[str] = None):
        self.backend = resolve("motifs", backend)
        self.n = len(seqs)
        seqs = [s.upper() for s in seqs]
        keep = [i for i, s in enumerate(seqs) if s and not _BAD.search(s)]
        self.rows = np.asarray(keep, np.int64)                       # block window -> input row
        lens = np.array([len(seqs[i]) for i in keep], np.int64)
        self.offsets = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        self.text = "".join(seqs[i] for i in keep)
        self._engine = None

    def _prepare(self):
        if self._engine is None and self.text:
            if self.backend == "lightmotif":
                import lightmotif
                self._engine = lightmotif.stripe(self.text)
            elif self.backend == "biopython":
                from Bio.Seq import Seq
                self._engine = Seq(self.text)
            else:
                self._engine = self.text
        return self._engine

    def _positions(self, mat, threshold):
        """Start positions in the block of hits of one (W x 4) log-odds matrix."""
        eng = self._prepare()
        if eng is None:
            return np.zeros(0, np.int64)
        if self.backend == "lightmotif":
            import lightmotif
            sm = lightmotif.ScoringMatrix({b: mat[:, k].tolist() for k, b in enumerate("ACGT")}
                                          | {"N": [float("-inf")] * len(mat)})
            return np.fromiter((h.position for h in lightmotif.scan(sm, eng, threshold=threshold)), np.int64)
        if self.backend == "moods":
            import MOODS.scan
            sc = MOODS.scan.Scanner(7)
            sc.set_motifs([mat.T.tolist()], [0.25] * 4, [float(threshold)])
            return np.fromiter((h.pos for h in sc.scan(eng)[0]), np.int64)
        from Bio.motifs.matrix import PositionSpecificScoringMatrix
        ps = PositionSpecificScoringMatrix("ACGT", {b: mat[:, k].tolist() for k, b in enumerate("ACGT")})
        with np.errstate(invalid="ignore"):
            s = np.asarray(ps.calculate(eng), np.float64)
        return np.flatnonzero(s >= threshold).astype(np.int64)

    def hits(self, mat, threshold, both: bool = False):
        """(input row, position in the window, strand) of every hit; strand 0 '+', 1 '-'."""
        w = len(mat)
        out_r, out_p, out_s = [], [], []
        for strand, m in ((0, mat), (1, _rc(mat))) if both else ((0, mat),):
            pos = self._positions(np.ascontiguousarray(m), threshold)
            if not len(pos):
                continue
            k = np.searchsorted(self.offsets, pos, side="right") - 1
            ok = (k < len(self.rows)) & (pos + w <= self.offsets[np.minimum(k + 1, len(self.rows))])
            out_r.append(self.rows[k[ok]])
            out_p.append(pos[ok] - self.offsets[k[ok]])
            out_s.append(np.full(int(ok.sum()), strand, np.int8))
        if not out_r:
            z = np.zeros(0, np.int64)
            return z, z.copy(), np.zeros(0, np.int8)
        return np.concatenate(out_r), np.concatenate(out_p), np.concatenate(out_s)

    def counts(self, mat, threshold, both: bool = False) -> np.ndarray:
        """Hits per input window."""
        r, _, _ = self.hits(mat, threshold, both)
        return np.bincount(r, minlength=self.n).astype(np.int64)

    def moods_counts(self, mats, thresholds, both: bool = False) -> np.ndarray:
        """MOODS scans many matrices in one pass: (windows x motifs) counts."""
        import MOODS.scan
        out = np.zeros((self.n, len(mats)), np.int64)
        if not self.text:
            return out
        thr = np.broadcast_to(np.asarray(thresholds, np.float64), (len(mats),))
        all_m, owner, widths, all_t = [], [], [], []
        for j, m in enumerate(mats):
            for mm in ((m, _rc(m)) if both else (m,)):
                all_m.append(np.ascontiguousarray(mm).T.tolist())
                owner.append(j)
                widths.append(len(m))
                all_t.append(float(thr[j]))
        sc = MOODS.scan.Scanner(7)
        sc.set_motifs(all_m, [0.25] * 4, all_t)
        res = sc.scan(self.text)
        for k, hits in enumerate(res):
            if not hits:
                continue
            pos = np.fromiter((h.pos for h in hits), np.int64)
            i = np.searchsorted(self.offsets, pos, side="right") - 1
            ok = (i < len(self.rows)) & (pos + widths[k] <= self.offsets[np.minimum(i + 1, len(self.rows))])
            out[:, owner[k]] += np.bincount(self.rows[i[ok]], minlength=self.n)
        return out
