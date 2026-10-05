"""Interval backends: overlap, nearest, merge and point lookups.

Every function takes columnar ``Loci`` that share one Genome and returns row
numbers, so the caller never sees which engine ran. The answers are
normalised to genomeblocks' rules whatever the engine:

* intervals are half-open; two intervals overlap when ``s1 < e2 and s2 < e1``;
* ``overlap_pairs`` is sorted by (query row, reference row);
* ``nearest`` gives distance 0 to an overlap and otherwise the gap in bases
  (book-ended intervals are 0 apart), -1 when the chromosome has nothing;
  on an exact tie the engines may pick different, equally near, rows;
* ``merge`` fuses overlapping and book-ended intervals (bedtools merge).

Zero-length intervals are points between two bases: ``[p, p)`` overlaps
``[s, e)`` exactly when ``s < p < e``. Every engine's pairs are filtered with
that rule, so the engines agree on them (bedtools, which widens empty
intervals on its side, may still miss pairs that need the widened base).

Support matrix::

                 genomeblocks cgranges ncls bioframe pyranges bedtools
    overlap          yes        yes    yes    yes      yes      yes
    nearest          yes         -      -     yes      yes      yes
    merge            yes         -      -     yes      yes      yes
    point lookup     yes        yes    yes    yes      yes      yes
"""
from __future__ import annotations

import numpy as np

from .. import _intervals as K
from . import resolve, unsupported

_NEAREST = ("genomeblocks", "bioframe", "pyranges", "bedtools")
_MERGE = ("genomeblocks", "bioframe", "pyranges", "bedtools")


# ── helpers: hand columns to the engines with codes as chromosome labels ──
# Both sides share one Genome, so the integer code is a safe chromosome name
# for every engine (no decoding, no clashes, identical on both sides).

def _frame(L, cols=("chrom", "start", "end")):
    import pandas as pd
    c, s, e = cols
    return pd.DataFrame({c: L.codes.astype(str), s: L.starts, e: L.ends})


def _sorted_pairs(qi, ri):
    qi, ri = np.asarray(qi, np.int64), np.asarray(ri, np.int64)
    o = np.lexsort((ri, qi))
    return qi[o], ri[o]


def _gap(q, r, qi, ri):
    """Distance in bases between query rows qi and reference rows ri (0 on overlap)."""
    d = np.maximum(r.starts[ri] - q.ends[qi], q.starts[qi] - r.ends[ri])
    return np.maximum(d, 0)


# ── overlap ────────────────────────────────────────────────────────────────

def overlap_pairs(q, r, *, backend=None):
    """(query rows, reference rows) of every overlapping pair."""
    b = resolve("intervals", backend)
    if len(q) == 0 or len(r) == 0:
        z = np.zeros(0, np.int64)
        return z, z.copy()
    if b == "genomeblocks":
        return _sorted_pairs(*K.overlap_pairs(q.codes, q.starts, q.ends, r.codes, r.starts, r.ends))
    qi, ri = _PAIRS[b](q, r)
    return _sorted_pairs(*_half_open(q, r, qi, ri))


def _half_open(q, r, qi, ri):
    """Keep the pairs that satisfy genomeblocks' rule (s1 < e2 and s2 < e1 on
    the same chromosome): engines differ on zero-length intervals."""
    qi, ri = np.asarray(qi, np.int64), np.asarray(ri, np.int64)
    if not len(qi):
        return qi, ri
    ok = ((q.codes[qi] == r.codes[ri]) & (q.starts[qi] < r.ends[ri]) & (r.starts[ri] < q.ends[qi]))
    return qi[ok], ri[ok]


def overlap_any(q, r, *, backend=None):
    """Boolean mask over the query rows: overlaps at least one reference row."""
    b = resolve("intervals", backend)
    if b == "genomeblocks":
        return K.overlaps_any(q.codes, q.starts, q.ends, r.codes, r.starts, r.ends)
    mask = np.zeros(len(q), bool)
    qi, _ = overlap_pairs(q, r, backend=b)
    mask[qi] = True
    return mask


def _pairs_cgranges(q, r):
    import cgranges
    idx = cgranges.cgranges()
    for i, (c, s, e) in enumerate(zip(r.codes.tolist(), r.starts.tolist(), r.ends.tolist())):
        idx.add(str(c), s, e, i)
    idx.index()
    qi, ri = [], []
    for i, (c, s, e) in enumerate(zip(q.codes.tolist(), q.starts.tolist(), q.ends.tolist())):
        for _, _, j in idx.overlap(str(c), s, e):
            qi.append(i)
            ri.append(j)
    return qi, ri


def _pairs_ncls(q, r):
    from ncls import NCLS
    qi, ri = [], []
    for c in np.intersect1d(np.unique(q.codes), np.unique(r.codes)):
        rm, qm = np.flatnonzero(r.codes == c), np.flatnonzero(q.codes == c)
        if not len(rm):
            continue
        tree = NCLS(r.starts[rm], r.ends[rm], rm.astype(np.int64))
        a, b = tree.all_overlaps_both(q.starts[qm], q.ends[qm], qm.astype(np.int64))
        qi.append(a)
        ri.append(b)
    if not qi:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(qi), np.concatenate(ri)


def _pairs_bioframe(q, r):
    import bioframe as bf
    ov = bf.overlap(_frame(q), _frame(r), how="inner", return_index=True, return_input=False)
    return ov["index"].to_numpy(np.int64), ov["index_"].to_numpy(np.int64)


def _pyranges(L, name):
    import pyranges as pr
    df = _frame(L, ("Chromosome", "Start", "End"))
    df[name] = np.arange(len(L))
    return pr.PyRanges(df)


def _pr_df(gr):
    return gr.df if hasattr(gr, "df") else gr                # pyranges 0.x / 1.x


def _pairs_pyranges(q, r):
    gq, gr = _pyranges(q, "__q"), _pyranges(r, "__r")
    j = gq.join_ranges(gr) if hasattr(gq, "join_ranges") else gq.join(gr)
    df = _pr_df(j)
    if not len(df):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return df["__q"].to_numpy(np.int64), df["__r"].to_numpy(np.int64)


def _bedtool(L):
    import pybedtools
    df = _frame(L)
    df["name"] = np.arange(len(L))
    return pybedtools.BedTool.from_dataframe(df)


def _read_bt(bt):
    """A BedTool result as a header-less DataFrame of strings ('' when empty)."""
    import pandas as pd
    fn = bt.fn
    try:
        return pd.read_csv(fn, sep="\t", header=None, dtype=str)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _pairs_bedtools(q, r):
    res = _read_bt(_bedtool(q).intersect(_bedtool(r), wa=True, wb=True))
    if not len(res):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return res[3].to_numpy(np.int64), res[7].to_numpy(np.int64)   # half-open filter in overlap_pairs


_PAIRS = {"cgranges": _pairs_cgranges, "ncls": _pairs_ncls, "bioframe": _pairs_bioframe,
          "pyranges": _pairs_pyranges, "bedtools": _pairs_bedtools}


# ── nearest ────────────────────────────────────────────────────────────────

def nearest(q, r, *, backend=None):
    """(reference row, distance) of the nearest reference interval per query row."""
    b = resolve("intervals", backend)
    if b == "genomeblocks":
        return K.nearest(q.codes, q.starts, q.ends, r.codes, r.starts, r.ends)
    if b not in _NEAREST:
        raise unsupported("intervals", b, "nearest", _NEAREST)
    best = np.full(len(q), -1, np.int64)
    if len(q) and len(r):
        qi, ri = _NEAREST_IMPL[b](q, r)
        qi, ri = np.asarray(qi, np.int64), np.asarray(ri, np.int64)
        ok = ri >= 0
        best[qi[ok]] = ri[ok]
        # engines break ties among overlapping intervals their own way; use the
        # genomeblocks rule (lowest start, then lowest row) on the engine's overlaps
        pq, pr_ = overlap_pairs(q, r, backend=b)
        if len(pq):
            o = np.lexsort((pr_, r.starts[pr_], pq))
            pq, pr_ = pq[o], pr_[o]
            first = np.r_[True, pq[1:] != pq[:-1]]
            best[pq[first]] = pr_[first]
    dist = np.full(len(q), -1, np.int64)
    hit = best >= 0
    dist[hit] = _gap(q, r, np.flatnonzero(hit), best[hit])
    return best, dist


def _nearest_bioframe(q, r):
    import bioframe as bf
    cl = bf.closest(_frame(q), _frame(r), k=1, return_index=True, return_distance=False,
                    return_input=False)
    ri = cl["index_"].to_numpy(dtype=float, na_value=np.nan)
    return cl["index"].to_numpy(np.int64), np.where(np.isnan(ri), -1, ri).astype(np.int64)


def _nearest_pyranges(q, r):
    gq, gr = _pyranges(q, "__q"), _pyranges(r, "__r")
    nb = gq.nearest_ranges(gr) if hasattr(gq, "nearest_ranges") else gq.nearest(gr)
    df = _pr_df(nb)
    if not len(df):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return df["__q"].to_numpy(np.int64), df["__r"].to_numpy(np.int64)


def _nearest_bedtools(q, r):
    res = _read_bt(_bedtool(q).sort().closest(_bedtool(r).sort(), t="first"))
    if not len(res):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    ri = res[7].replace(".", "-1").to_numpy(np.int64)
    return res[3].to_numpy(np.int64), ri


_NEAREST_IMPL = {"bioframe": _nearest_bioframe, "pyranges": _nearest_pyranges,
                 "bedtools": _nearest_bedtools}


# ── merge ──────────────────────────────────────────────────────────────────

def merge(L, *, backend=None):
    """Merged blocks of ``L``: (codes, starts, ends, first) where ``first`` is the
    row of the first interval (lowest start) of each block."""
    b = resolve("intervals", backend)
    if b == "genomeblocks" or len(L) == 0:
        return K.merge(L.codes, L.starts, L.ends)
    if b not in _MERGE:
        raise unsupported("intervals", b, "merge", _MERGE)
    c, s, e = _MERGE_IMPL[b](L)
    c, s, e = np.asarray(c, np.int64), np.asarray(s, np.int64), np.asarray(e, np.int64)
    o = np.lexsort((s, c))
    c, s, e = c[o], s[o], e[o]
    # first input row of each block: the lowest-start row at the block start
    g = K.gpos(L.codes, L.starts)
    order = np.argsort(g, kind="stable")
    first = order[np.searchsorted(g[order], K.gpos(c, s), side="left")]
    return c, s, e, first


def _merge_bioframe(L):
    import bioframe as bf
    m = bf.merge(_frame(L), min_dist=0)
    return m["chrom"].astype(np.int64), m["start"], m["end"]


def _merge_pyranges(L):
    gr = _pyranges(L, "__i")
    m = gr.merge_overlaps() if hasattr(gr, "merge_overlaps") else gr.merge()
    df = _pr_df(m)
    return df["Chromosome"].astype(str).astype(np.int64), df["Start"], df["End"]


def _merge_bedtools(L):
    res = _read_bt(_bedtool(L).sort().merge())
    return res[0].astype(np.int64), res[1].astype(np.int64), res[2].astype(np.int64)


_MERGE_IMPL = {"bioframe": _merge_bioframe, "pyranges": _merge_pyranges, "bedtools": _merge_bedtools}


# ── point lookups ──────────────────────────────────────────────────────────

def point_rows(L, code: int, start: int, end: int, *, backend=None) -> np.ndarray:
    """Rows of ``L`` overlapping one window (sorted). Indexes are cached on ``L``."""
    b = resolve("intervals", backend)
    if b == "genomeblocks":
        idx = L._index_cache("genomeblocks", lambda: K.PointIndex(L.codes, L.starts, L.ends))
        return idx.query(code, start, end)
    if b == "cgranges":
        def build():
            import cgranges
            ix = cgranges.cgranges()
            for i, (c, s, e) in enumerate(zip(L.codes.tolist(), L.starts.tolist(), L.ends.tolist())):
                ix.add(str(c), s, e, i)
            ix.index()
            return ix
        ix = L._index_cache("cgranges", build)
        rows = np.array(sorted(j for *_, j in ix.overlap(str(code), int(start), int(end))), np.int64)
        return _point_filter(L, rows, start, end)
    if b == "ncls":
        def build():
            from ncls import NCLS
            trees = {}
            for c in np.unique(L.codes):
                m = np.flatnonzero(L.codes == c)
                if len(m):
                    trees[int(c)] = NCLS(L.starts[m], L.ends[m], m.astype(np.int64))
            return trees
        tree = L._index_cache("ncls", build).get(int(code))
        if tree is None:
            return np.zeros(0, np.int64)
        rows = np.array(sorted(j for *_, j in tree.find_overlap(int(start), int(end))), np.int64)
        return _point_filter(L, rows, start, end)
    one = type(L)([code], [start], [end], genome=L.genome)
    _, ri = overlap_pairs(one, L, backend=b)
    return np.sort(ri)


def _point_filter(L, rows, start, end):
    """Rows that overlap ``[start, end)`` under the half-open rule."""
    if not len(rows):
        return rows
    ok = (L.starts[rows] < end) & (start < L.ends[rows])
    return rows[ok]
