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
