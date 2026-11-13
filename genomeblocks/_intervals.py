"""Vectorised interval kernels on the genome-wide axis.

Every chromosome is laid end to end on one axis (``code * STRIDE + pos``), so
an operation over the whole genome is a handful of numpy calls — no
per-chromosome or per-interval Python loop. All inputs share one Genome, so
codes are directly comparable.
"""
from __future__ import annotations

import numpy as np

STRIDE = np.int64(1) << 40          # > any chromosome length (1.1e12)


def gpos(codes, pos):
    return np.asarray(codes, np.int64) * STRIDE + np.asarray(pos, np.int64)


def overlaps_any(qc, qs, qe, rc, rs, re_):
    """mask[i] = query i overlaps at least one reference interval (half-open).

    Sort the references by start, keep a running max of their ends; a query
    [s, e) overlaps iff the largest end among references starting before e
    exceeds s. Works for references of any length. O((n+m) log m).
    """
    if len(rs) == 0:
        return np.zeros(len(qs), bool)
    bs, be = gpos(rc, rs), gpos(rc, re_)
    o = np.argsort(bs, kind="stable")
    bs, run_max = bs[o], np.maximum.accumulate(be[o])
    k = np.searchsorted(bs, gpos(qc, qe), side="left")
    hit = k > 0
    hit[hit] = run_max[k[hit] - 1] > gpos(qc, qs)[hit]
    return hit


def overlap_pairs(qc, qs, qe, rc, rs, re_):
    """All (query, reference) index pairs that overlap (half-open).

    References are sorted once; a reference can only overlap [s, e) if it
    starts in (s - maxlen, e), so each query's candidates are one
    searchsorted range. Returns (qi, ri) sorted by qi.
    """
    if len(qs) == 0 or len(rs) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    rg = gpos(rc, rs)
    o = np.argsort(rg, kind="stable")
    rg_s, rend_s = rg[o], gpos(rc, re_)[o]
    maxlen = max(int((np.asarray(re_) - np.asarray(rs)).max()), 0)
    qgs, qge = gpos(qc, qs), gpos(qc, qe)
    lo = np.searchsorted(rg_s, qgs - maxlen, side="right")
    hi = np.searchsorted(rg_s, qge, side="left")
    n = np.maximum(hi - lo, 0)
    qi = np.repeat(np.arange(len(qs)), n)
    start = np.repeat(lo - np.concatenate([[0], np.cumsum(n)[:-1]]), n)
    k = start + np.arange(n.sum())
    keep = rend_s[k] > qgs[qi]
    return qi[keep], o[k[keep]]


def merge(codes, starts, ends):
    """bedtools-merge semantics: overlapping or book-ended intervals fuse.

    Returns (codes, starts, ends, first) where ``first`` is the input index of
    the first interval (in sorted order) of each merged block.
    """
    if len(starts) == 0:
        z = np.zeros(0, np.int64)
        return z, z, z, z
    gs, ge = gpos(codes, starts), gpos(codes, ends)
    o = np.argsort(gs, kind="stable")
    gs, ge = gs[o], ge[o]
    run_end = np.maximum.accumulate(ge)
    new = np.empty(len(gs), bool)
    new[0] = True
    new[1:] = gs[1:] > run_end[:-1]
    first = np.flatnonzero(new)
    last = np.append(first[1:], len(gs)) - 1
    c = np.asarray(codes)[o][first]
    return c, gs[first] - c * STRIDE, run_end[last] - c * STRIDE, o[first]


def count_overlaps(qc, qs, qe, rc, rs, re_):
    """n[i] = number of reference intervals overlapping query i (half-open).

    A reference [s, e) overlaps [a, b) iff s < b and e > a; every reference
    ending at or before a also starts before b, so the count is
    #(starts < b) - #(ends <= a) — two searchsorted calls. The rule also
    gives genomeblocks' zero-length behaviour ([p, p) overlaps [s, e) iff
    s < p < e). Chromosomes are apart on the genome axis, so references on
    other chromosomes cancel out.
    """
    if len(rs) == 0 or len(qs) == 0:
        return np.zeros(len(qs), np.int64)
    rs_s = np.sort(gpos(rc, rs))
    re_s = np.sort(gpos(rc, re_))
    return (np.searchsorted(rs_s, gpos(qc, qe), side="left")
            - np.searchsorted(re_s, gpos(qc, qs), side="right")).astype(np.int64)


def _merged_axis(codes, starts, ends):
    """Merged reference blocks on the genome axis: (block starts, block ends,
    cumulative covered bases before each block). Zero-length intervals cover no
    base, so they are left out (they would otherwise split a piece or a gap)."""
    c, s, e, _ = merge(*_nonempty(codes, starts, ends))
    bs, be = gpos(c, s), gpos(c, e)
    cum = np.concatenate([[0], np.cumsum(be - bs)])
    return bs, be, cum


def covered_bases(qc, qs, qe, rc, rs, re_):
    """bases[i] = bases of query i covered by the reference set (overlapping
    references counted once: bedtools coverage's covered-bases column)."""
    if len(rs) == 0 or len(qs) == 0:
        return np.zeros(len(qs), np.int64)
    bs, be, cum = _merged_axis(rc, rs, re_)
    a, b = gpos(qc, qs), gpos(qc, qe)
    j0 = np.searchsorted(be, a, side="right")             # first block ending after a
    j1 = np.searchsorted(bs, b, side="left")              # first block starting at or after b
    n = j1 - j0
    out = cum[np.maximum(j1, j0)] - cum[j0]
    has = n > 0
    out[has] -= np.maximum(a[has] - bs[j0[has]], 0)        # clip the first block at a
    out[has] -= np.maximum(be[j1[has] - 1] - b[has], 0)    # and the last at b
    return np.maximum(out, 0).astype(np.int64)


def subtract(qc, qs, qe, rc, rs, re_):
    """bedtools subtract: the parts of each query interval no reference covers.

    Returns (row, starts, ends): one row per surviving piece, ``row`` the query
    row it came from, pieces in query order and left to right.
    """
    if len(qs) == 0:
        z = np.zeros(0, np.int64)
        return z, z, z
    qs, qe = np.asarray(qs, np.int64), np.asarray(qe, np.int64)
    if len(rs) == 0:
        return np.arange(len(qs)), qs.copy(), qe.copy()
    bs, be, _ = _merged_axis(rc, rs, re_)
    off = np.asarray(qc, np.int64) * STRIDE
    a, b = off + qs, off + qe
    j0 = np.searchsorted(be, a, side="right")
    j1 = np.searchsorted(bs, b, side="left")
    n = np.maximum(j1 - j0, 0)                              # blocks inside each query; n + 1 candidate pieces
    row = np.repeat(np.arange(len(qs)), n + 1)
    k = np.arange(len(row)) - np.repeat(np.cumsum(n + 1) - (n + 1), n + 1)
    j = j0[row] + k
    ps = np.where(k == 0, a[row], be[np.minimum(j - 1, len(be) - 1)])
    pe = np.where(k == n[row], b[row], bs[np.minimum(j, len(bs) - 1)])
    ps, pe = np.maximum(ps, a[row]), np.minimum(pe, b[row])
    keep = pe > ps
    row, ps, pe = row[keep], ps[keep], pe[keep]
    return row, ps - off[row], pe - off[row]


def _nonempty(codes, starts, ends):
    codes, starts, ends = (np.asarray(x, np.int64) for x in (codes, starts, ends))
    keep = ends > starts
    return codes[keep], starts[keep], ends[keep]


def complement(codes, starts, ends, sizes):
    """bedtools complement: the gaps between merged intervals within each
    chromosome ``[0, size)``; ``sizes`` holds one length per genome code (a
    chromosome with no interval is one whole gap, a code with size < 0 is left
    out). Returns (codes, starts, ends) in (code, start) order.

    The merged blocks are already in that order, and the gap ending at each
    block comes right before it, so only each chromosome's closing gap (and the
    chromosomes with no block) are inserted — no sort of the gaps."""
    sizes = np.asarray(sizes, np.int64)
    c, s, e = _nonempty(codes, starts, ends)
    if len(c):
        c, s, e, _ = merge(c, s, e)
    first = np.r_[True, c[1:] != c[:-1]] if len(c) else np.zeros(0, bool)
    gs = np.where(first, 0, np.r_[0, e[:-1]]) if len(c) else np.zeros(0, np.int64)
    gc, ge = c.astype(np.int64), np.minimum(s, sizes[c]) if len(c) else np.zeros(0, np.int64)
    # each chromosome's closing gap goes after its last block's gap
    last = np.flatnonzero(np.r_[c[1:] != c[:-1], True]) if len(c) else np.zeros(0, np.int64)
    at = last + 1
    tc = gc[last]
    gc = np.insert(gc, at, tc)
    gs = np.insert(gs, at, e[last])
    ge = np.insert(ge, at, sizes[tc])
    # chromosomes with a size and no block: one whole gap, at their place in code order
    known = np.flatnonzero(sizes >= 0)
    empty = np.setdiff1d(known, np.unique(c))
    if len(empty):
        at = np.searchsorted(gc, empty, side="left")
        gc, gs, ge = np.insert(gc, at, empty), np.insert(gs, at, 0), np.insert(ge, at, sizes[empty])
    keep = (ge > gs) & (sizes[gc] >= 0)
    return gc[keep], gs[keep], ge[keep]


def running_argmax(a: np.ndarray) -> np.ndarray:
    """idx[i] = argmax(a[:i+1]) (first index of the running maximum)."""
    if len(a) == 0:
        return np.zeros(0, np.int64)
    run = np.maximum.accumulate(a)
    is_new = np.r_[True, a[1:] > run[:-1]]
    idx = np.where(is_new, np.arange(len(a)), 0)
    return np.maximum.accumulate(idx)


def nearest(qc, qs, qe, rc, rs, re_):
    """(row in the reference, distance) of the nearest reference interval per query.

    Overlaps have distance 0 (the overlapping interval with the lowest start
    wins); otherwise the closer of the left and right neighbour, the left one
    on a tie. Distance is the gap in bases: book-ended intervals are 0 apart.
    Queries with nothing on their chromosome get (-1, -1).
    """
    qc, qs, qe = np.asarray(qc), np.asarray(qs, np.int64), np.asarray(qe, np.int64)
    rc, rs, re_ = np.asarray(rc), np.asarray(rs, np.int64), np.asarray(re_, np.int64)
    n = len(qs)
    best = np.full(n, -1, np.int64)
    dist = np.full(n, np.iinfo(np.int64).max, np.int64)
    if n == 0 or len(rs) == 0:
        return best, np.full(n, -1, np.int64)
    og = gpos(rc, rs)
    order = np.lexsort((rs, np.asarray(rc, np.int64)))
    og, oe = og[order], gpos(rc, re_)[order]
    oc = np.asarray(rc)[order]
    gqs, gqe = gpos(qc, qs), gpos(qc, qe)
    # left: the interval with the largest end among those starting before q.end
    k = np.searchsorted(og, gqe, side="left")
    run_max = np.maximum.accumulate(oe)
    arg_max = running_argmax(oe)
    has = k > 0
    j = np.where(has, arg_max[np.maximum(k - 1, 0)], -1)
    same = has & (oc[np.maximum(j, 0)] == qc)
    d_left = np.where(same, np.maximum(gqs - run_max[np.maximum(k - 1, 0)], 0), dist)
    # overlaps: distance 0, the overlapping interval with the lowest start wins
    pq, pr_ = overlap_pairs(qc, qs, qe, rc, rs, re_)
    ov = np.zeros(n, bool)
    if len(pq):
        srt = np.lexsort((pr_, rs[pr_], pq))
        pq, pr_ = pq[srt], pr_[srt]
        firstq = np.r_[True, pq[1:] != pq[:-1]]
        best[pq[firstq]] = pr_[firstq]
        dist[pq[firstq]] = 0
        ov[pq[firstq]] = True
    # right: the first interval starting at or after q.end
    ok2 = k < len(og)
    j2 = np.where(ok2, np.minimum(k, len(og) - 1), 0)
    same2 = ok2 & (oc[j2] == qc)
    d_right = np.where(same2, og[j2] - gqe, dist)
    use_left = ~ov & same & (d_left <= d_right)
    use_right = ~ov & same2 & ~use_left
    best[use_left] = order[j[use_left]]
    dist[use_left] = d_left[use_left]
    best[use_right] = order[j2[use_right]]
    dist[use_right] = d_right[use_right]
    dist[best < 0] = -1
    return best, dist


class PointIndex:
    """Single-window lookups on the genome axis in O(log n + k).

    Rows are sorted by start once; a running maximum of their ends bounds
    the scan on the left (no interval before the first position whose
    running max exceeds the query start can reach it). Built lazily by
    ``Loci`` and cached until the rows change.
    """

    __slots__ = ("order", "gs", "ge", "run")

    def __init__(self, codes, starts, ends):
        gs, ge = gpos(codes, starts), gpos(codes, ends)
        self.order = np.argsort(gs, kind="stable")
        self.gs, self.ge = gs[self.order], ge[self.order]
        self.run = np.maximum.accumulate(self.ge) if len(self.ge) else self.ge

    def query(self, code: int, start: int, end: int) -> np.ndarray:
        """Rows (sorted) overlapping [start, end) on chromosome ``code``."""
        if not len(self.gs):
            return np.zeros(0, np.int64)
        qs = np.int64(code) * STRIDE + np.int64(start)
        qe = np.int64(code) * STRIDE + np.int64(end)
        lo = int(np.searchsorted(self.run, qs, side="right"))
        hi = int(np.searchsorted(self.gs, qe, side="left"))
        if hi <= lo:
            return np.zeros(0, np.int64)
        k = np.flatnonzero(self.ge[lo:hi] > qs) + lo
        return np.sort(self.order[k])
