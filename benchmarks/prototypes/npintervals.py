"""Prototype: a dependency-free, columnar interval engine in numpy.

Intervals are three parallel arrays (chromosome code, start, end) plus the
list of chromosome names. Every operation works on all chromosomes at once by
placing them on one genome-wide axis (``code * STRIDE + pos``), so there is no
per-chromosome Python loop.

This is what "implement it from scratch" would look like for the whole-set
operations (intersect / difference / merge); point lookups stay with an index
such as cgranges.
"""
from __future__ import annotations

import numpy as np

STRIDE = np.int64(1) << 40          # > any chromosome length (1.1e12)


class Intervals:
    """Columnar interval set: codes[i] indexes into names."""
    __slots__ = ("codes", "starts", "ends", "names")

    def __init__(self, codes, starts, ends, names):
        self.codes = np.asarray(codes, np.int64)
        self.starts = np.asarray(starts, np.int64)
        self.ends = np.asarray(ends, np.int64)
        self.names = list(names)

    def __len__(self):
        return len(self.starts)

    # ── conversion ────────────────────────────────────────────────────────
    @classmethod
    def from_loci(cls, loci):
        """From a list of Locus objects (one pass in Python)."""
        n = len(loci)
        chroms = [l.chrom for l in loci]
        names, codes = np.unique(np.asarray(chroms, dtype=object), return_inverse=True)
        starts = np.fromiter((l.start for l in loci), np.int64, n)
        ends = np.fromiter((l.end for l in loci), np.int64, n)
        return cls(codes, starts, ends, list(names))

    @classmethod
    def from_frame(cls, df, chrom="chrom", start="start", end="end"):
        cat = df[chrom].astype("category")
        return cls(cat.cat.codes.to_numpy(), df[start].to_numpy(), df[end].to_numpy(),
                   [str(c) for c in cat.cat.categories])

    def recode(self, names):
        """Same intervals, coded against another name list (unknown -> -1)."""
        lut = {n: i for i, n in enumerate(names)}
        m = np.array([lut.get(n, -1) for n in self.names], np.int64)
        return Intervals(m[self.codes], self.starts, self.ends, names)

    def to_loci(self, Locus, Loci):
        names = self.names
        return Loci(Locus(names[c], int(s), int(e))
                    for c, s, e in zip(self.codes.tolist(), self.starts.tolist(), self.ends.tolist()))

    # ── helpers ───────────────────────────────────────────────────────────
    def _g(self):
        return self.codes * STRIDE + self.starts, self.codes * STRIDE + self.ends


def overlaps_any(a: Intervals, b: Intervals) -> np.ndarray:
    """Boolean mask: does each interval of ``a`` overlap any interval of ``b``?

    Sort ``b`` by start on the genome axis and keep a running maximum of its
    ends. For an interval [s, e) of ``a``, the candidates are the ``b``
    intervals starting before ``e`` (one searchsorted); it overlaps iff the
    largest end among them exceeds ``s``. O((n + m) log m), fully vectorised.
    """
    b = b.recode(a.names) if b.names != a.names else b
    keep = b.codes >= 0
    bs, be = b.codes[keep] * STRIDE + b.starts[keep], b.codes[keep] * STRIDE + b.ends[keep]
    if bs.size == 0:
        return np.zeros(len(a), bool)
    o = np.argsort(bs, kind="stable")
    bs, run_max = bs[o], np.maximum.accumulate(be[o])
    as_, ae = a._g()
    k = np.searchsorted(bs, ae, side="left")          # b starting before a.end
    hit = k > 0
    hit[hit] = run_max[k[hit] - 1] > as_[hit]
    return hit


def merge(a: Intervals) -> Intervals:
    """bedtools-merge semantics (overlapping or book-ended intervals fuse)."""
    if len(a) == 0:
        return a
    gs, ge = a._g()
    o = np.lexsort((gs,))
    gs, ge, codes = gs[o], ge[o], a.codes[o]
    run_end = np.maximum.accumulate(ge)
    new = np.empty(len(gs), bool)
    new[0] = True
    new[1:] = gs[1:] > run_end[:-1]                    # gap -> new block
    first = np.flatnonzero(new)
    last = np.append(first[1:], len(gs)) - 1
    ms, me = gs[first], run_end[last]
    c = codes[first]
    return Intervals(c, ms - c * STRIDE, me - c * STRIDE, a.names)


def slop(a: Intervals, n: int) -> Intervals:
    return Intervals(a.codes, np.maximum(a.starts - n, 0), a.ends + n, a.names)
