"""Prototype: a columnar Loci that keeps genomeblocks' public Loci API.

Storage is four numpy columns (chromosome code, start, end, strand code)
plus the chromosome name list. ``loci[i]`` and iteration hand out
``LocusView`` objects: real ``Locus`` subclasses (so ``isinstance`` checks
and every Locus method keep working) whose fields read from and write to the
columns. Whole-set operations run on the arrays; everything else in
genomeblocks keeps working through the views, unchanged.
"""
from __future__ import annotations

import numpy as np

from genomeblocks.loci import Loci
from genomeblocks.locus import Locus

from .npintervals import STRIDE, Intervals, merge as _np_merge, overlaps_any

_STRANDS = (".", "+", "-")
_SCODE = {".": 0, "+": 1, "-": 2}


class LocusView(Locus):
    """A Locus whose four fields live in a ColumnarLoci (read/write-through)."""

    def __init__(self, owner, i):                 # no dataclass __init__: nothing is copied
        self.__dict__["_o"] = owner
        self.__dict__["_i"] = i

    def _set(self, name, value):
        o, i = self.__dict__["_o"], self.__dict__["_i"]
        if name == "chrom":
            o.codes[i] = o._code_for(value)
        elif name == "strand":
            o.strands[i] = _SCODE[value]
        else:
            getattr(o, name + "s")[i] = value
        o._dirty()

    chrom = property(lambda s: s._o.names[s._o.codes[s._i]], lambda s, v: s._set("chrom", v))
    start = property(lambda s: int(s._o.starts[s._i]), lambda s, v: s._set("start", v))
    end = property(lambda s: int(s._o.ends[s._i]), lambda s, v: s._set("end", v))
    strand = property(lambda s: _STRANDS[s._o.strands[s._i]], lambda s, v: s._set("strand", v))

    def __setattr__(self, name, value):           # route field writes through the properties
        if name in ("chrom", "start", "end", "strand"):
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value


class ColumnarLoci(Loci):
    """Loci stored as columns. Same API; objects are made on access."""

    def __init__(self, iterable=(), *, filename=None, columns=None):
        list.__init__(self)                         # the list part stays empty
        self.filename = filename
        self._uids = None
        self._cgr = None
        if columns is not None:
            self.codes, self.starts, self.ends, self.strands, self.names = columns
        else:
            items = list(iterable)
            iv = Intervals.from_loci(items) if items else Intervals([], [], [], [])
            self.codes, self.starts, self.ends, self.names = iv.codes, iv.starts, iv.ends, iv.names
            self.strands = np.fromiter((_SCODE.get(l.strand, 0) for l in items), np.int8, len(items))

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, filename, filetype=None):
        import pandas as pd
        df = pd.read_csv(filename, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 5],
                         names=["chrom", "start", "end", "strand"],
                         dtype={"chrom": "category", "strand": "category"})
        names = [str(c) for c in df["chrom"].cat.categories]
        order = np.argsort(names)
        rank = np.empty(len(names), np.int64); rank[order] = np.arange(len(names))
        codes = rank[df["chrom"].cat.codes.to_numpy()]
        strands = df["strand"].map(_SCODE).fillna(0).to_numpy(np.int8).copy()   # pandas 3 arrays are read-only
        return cls(filename=filename, columns=(codes, df["start"].to_numpy(np.int64).copy(),
                                               df["end"].to_numpy(np.int64).copy(), strands,
                                               [names[i] for i in order]))

    @classmethod
    def from_loci(cls, loci):
        return cls(loci)

    def _take(self, idx):
        return ColumnarLoci(columns=(self.codes[idx], self.starts[idx], self.ends[idx],
                                     self.strands[idx], self.names))

    def _code_for(self, name):
        if name not in self.names:
            self.names.append(name)
        return self.names.index(name)

    def _dirty(self):
        self._uids = None
        self._cgr = None

    def _iv(self):
        return Intervals(self.codes, self.starts, self.ends, self.names)

    # ── list behaviour ────────────────────────────────────────────────────
    def __len__(self):
        return len(self.starts)

    def __bool__(self):
        return len(self.starts) > 0

    def __iter__(self):
        for i in range(len(self.starts)):
            yield LocusView(self, i)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            n = len(self.starts)
            if not -n <= key < n:
                raise IndexError(key)
            return LocusView(self, int(key) % n)
        if isinstance(key, str):
            try:
                return LocusView(self, self.uids[key])
            except KeyError:
                raise KeyError(f"No Locus with UID '{key}' found.")
        if isinstance(key, slice):
            return self._take(np.arange(len(self.starts))[key])
        return self._take(np.asarray(key))

    def append(self, l):
        self.extend([l])

    def extend(self, items):
        other = items if isinstance(items, ColumnarLoci) else ColumnarLoci(items)
        if not len(other):
            return
        lut = np.array([self._code_for(n) for n in other.names], np.int64)
        self.codes = np.concatenate([self.codes, lut[other.codes]])
        self.starts = np.concatenate([self.starts, other.starts])
        self.ends = np.concatenate([self.ends, other.ends])
        self.strands = np.concatenate([self.strands, other.strands])
        self._dirty()

    def copy(self):
        return ColumnarLoci(filename=self.filename,
                            columns=(self.codes.copy(), self.starts.copy(), self.ends.copy(),
                                     self.strands.copy(), list(self.names)))

    def __str__(self):
        return f"Loci(n={len(self)}, columnar)"
    __repr__ = __str__

    # ── indexes ───────────────────────────────────────────────────────────
    @property
    def uids(self):
        if self._uids is None:
            ch = np.asarray(self.names, dtype=object)[self.codes]
            st = np.asarray(_STRANDS, dtype=object)[self.strands]
            self._uids = {f"{c}:{s}-{e}({d})": i for i, (c, s, e, d) in
                          enumerate(zip(ch.tolist(), self.starts.tolist(), self.ends.tolist(), st.tolist()))}
        return self._uids

    def _build_cgr(self):
        from genomeblocks.loci import _get_cgranges, _PyIntervalIndex
        cg = _get_cgranges()
        idx = cg.cgranges() if cg is not None else _PyIntervalIndex()
        names = self.names
        for i, (c, s, e) in enumerate(zip(self.codes.tolist(), self.starts.tolist(), self.ends.tolist())):
            idx.add(names[c], s, e, i)
        idx.index()
        self._cgr = idx

    # ── whole-set operations on the arrays ────────────────────────────────
    @staticmethod
    def _cols(o):
        return o._iv() if isinstance(o, ColumnarLoci) else Intervals.from_loci(list(o))

    def intersect(self, o):
        return self._take(overlaps_any(self._iv(), self._cols(o)))

    def difference(self, o):
        return self._take(~overlaps_any(self._iv(), self._cols(o)))

    def __and__(self, o): return self.intersect(o)
    def __sub__(self, o): return self.difference(o)
    __truediv__ = __sub__

    def __add__(self, o):
        out = self.copy()
        out.extend(o)
        return out
    __or__ = __add__

    def __xor__(self, o):
        return (self - o) + (o - self)

    def slop(self, n):
        return ColumnarLoci(columns=(self.codes, np.maximum(self.starts - n, 0), self.ends + n,
                                     self.strands, self.names))

    def sort(self):
        rank = np.argsort(np.argsort(np.asarray(self.names, dtype=object)))
        return self._take(np.lexsort((self.starts, rank[self.codes])))

    def merge(self):
        if not len(self):
            return ColumnarLoci()
        s = self._take(np.lexsort((self.starts, self.codes)))
        m = _np_merge(s._iv())
        # keep the strand of the first interval of each block, as genomeblocks does
        g = s.codes * STRIDE + s.starts
        first = np.searchsorted(g, m.codes * STRIDE + m.starts)
        out = ColumnarLoci(columns=(m.codes, m.starts, m.ends, s.strands[first], m.names))
        return out.sort()                            # chromosome-name order, like today

    def to_frame(self, names=None):
        import pandas as pd
        return pd.DataFrame({"Chr": np.asarray(self.names, dtype=object)[self.codes],
                             "Start": self.starts, "End": self.ends,
                             "Strand": np.asarray(_STRANDS, dtype=object)[self.strands],
                             "Name": names if names is not None else list(self.uids)})


def annotations_columnar(genes, L):
    """Genes.annotations with one vectorised overlap test per region class."""
    import pandas as pd
    a = L._iv()
    annot = genes.annot
    label = np.full(len(L), "Intergenic", dtype=object)
    for key, name in (("body", "Intronic"), ("exon", "Exonic"), ("utr3", "3UTR"),
                      ("utr5", "5UTR"), ("prom", "Promoter-TSS")):     # lowest priority first
        label[overlaps_any(a, Intervals.from_loci(list(annot[key])))] = name
    return pd.DataFrame({"uid": list(L.uids), "annotation": label})


def bin_ranges_columnar(atlas, L):
    """Atlas._intervals_to_bin_ranges on the columns (no per-locus loop)."""
    off = np.array([atlas.chrom_offsets.get(n, -1) for n in L.names], np.int64)[L.codes]
    cmax = np.array([atlas._chrom_max_bins.get(n, 0) for n in L.names], np.int64)[L.codes]
    bs = atlas.bin_size
    s = np.maximum(L.starts // bs, 0)
    e = np.maximum((L.ends - 1) // bs + 1, s + 1)
    s, e = np.minimum(s, cmax), np.minimum(e, cmax)
    ok = (off >= 0) & (s < e)
    return np.column_stack([off[ok] + s[ok], off[ok] + e[ok]])
