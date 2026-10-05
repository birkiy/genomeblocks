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
