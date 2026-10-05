"""Super-enhancers, ROSE-style, on columnar tables.

    peaks ─▶ stitch (peaks within ``stitch`` bp become one region) ─▶ score = mean signal x width
          ─▶ slope-1 knee on the ranked scores ─▶ SEs

``knee`` is the same slope-1 rule as ``Architecture.elbow`` (prime hubs), so
SEs and hubs are cut the same way.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .loci import Loci


def knee(values) -> Tuple[int, np.ndarray]:
    """Slope-1 knee on a ranked curve. Returns (cutoff, order): ``order`` are the
    indices of the positive values sorted descending (stable), and the first
    ``cutoff`` of them are above the knee."""
    from scipy.ndimage import uniform_filter1d
    v = np.asarray(values, float)
    idx = np.flatnonzero(v > 0)
    order = idx[np.argsort(-v[idx], kind="stable")]
    y = v[order]
    if len(y) < 3:
        return len(y), order
    ys = y[::-1]
    m = len(ys)
    xn = np.arange(m) / (m - 1)
    yn = (ys - ys[0]) / (ys[-1] - ys[0])
    size = max(11, m // 200) if m >= 50 else 1                 # smooth long curves only
    yn_s = uniform_filter1d(yn, size=size, mode="nearest") if size > 1 else yn
    crossed = np.flatnonzero(np.gradient(yn_s, xn) >= 1.0)
    i = int(crossed[0]) if len(crossed) else m
    return m - i, order


def stitch_peaks(peaks, stitch: int = 12_500, *, backend: Optional[str] = None) -> Loci:
    """ROSE stitching: peaks whose gap is at most ``stitch`` bp form one region
    running from the first peak's start to the last peak's end (genome order).

    ``peaks`` is anything :func:`~genomeblocks.as_loci` takes; ``backend``
    picks the interval engine."""
    from .backends.intervals import overlap_pairs
    from .interop import as_loci
    peaks = as_loci(peaks)
    if not len(peaks):
        return Loci(genome=peaks.genome, is_sorted=True)
    clusters = peaks.slop(stitch // 2).merge(backend=backend)
    ci, pi = overlap_pairs(clusters, peaks, backend=backend)      # peak -> its cluster
    starts = np.full(len(clusters), np.iinfo(np.int64).max, np.int64)
    ends = np.full(len(clusters), -1, np.int64)
    np.minimum.at(starts, ci, peaks.starts[pi])
    np.maximum.at(ends, ci, peaks.ends[pi])
    n = np.bincount(ci, minlength=len(clusters))
    return Loci(clusters.codes, starts, ends, genome=peaks.genome, cols={"n_peaks": n}, is_sorted=True)


def call_se(peaks, bigwigs, *, stitch: int = 12_500, workers: int = 1, verbose: bool = False,
            return_all: bool = False, backend: Optional[str] = None, **signal_kw):
    """ROSE with genomeblocks: stitch peaks, score mean signal (averaged over
    ``bigwigs``) x width, keep the stitched regions above the slope-1 knee.

    ``peaks``: anything :func:`~genomeblocks.as_loci` takes. ``bigwigs``: one
    or several bigWig paths / open handles / a ``{name: path}`` dict (what
    :meth:`Loci.signal` takes; extra keywords such as ``backend=`` for the
    bigWig engine go to it).

    Returns the SEs in genome order with columns ``n_peaks``, ``score`` and
    ``rank`` (1 = strongest). ``return_all=True`` also returns every stitched
    region with its score (for the hockey-stick plot): ``(se, all_regions)``.
    """
    from .interop import as_loci
    peaks = as_loci(peaks)
    st = stitch_peaks(peaks, stitch)
    if not len(st):
        se = st.take([])
        return (se, st) if return_all else se
    signal_kw.setdefault("progress", False)
    signal_kw.setdefault("verbose", False)
    cube = st.signal(bigwigs, span=True, n_bins=1, workers=workers, backend=backend, **signal_kw)
    sig = np.nan_to_num(cube[:, :, 0]).mean(1)
    score = sig * st.lengths
    cut, order = knee(score)
    keep = order[:cut]
    rank = np.zeros(len(st), np.int64)
    rank[order] = np.arange(1, len(order) + 1)
    st.cols["score"], st.cols["rank"] = score, rank
    se = st.take(np.sort(keep))
    if verbose:
        print(f"[INFO] {len(peaks):,} peaks → {len(st):,} stitched regions → {len(se):,} SEs "
              f"(median {int(np.median(se.lengths)) if len(se) else 0:,} bp)")
    return (se, st) if return_all else se


def nearest_gene_within(regions, tss, maxd: int) -> np.ndarray:
    """For each region: the row in ``tss`` (1-bp TSS Loci) whose TSS is closest to
    the region centre, or -1 when that TSS lies more than ``maxd`` from the region
    (0 when it is inside it)."""
    from .interop import as_loci
    regions = as_loci(regions)
    tss = regions._check(as_loci(tss))
    out = np.full(len(regions), -1, np.int64)
    order = np.lexsort((tss.starts, tss.codes))
    tc, tp = tss.codes[order], tss.starts[order]
    for c in np.unique(regions.codes):
        m = np.flatnonzero(regions.codes == c)
        lo, hi = np.searchsorted(tc, c), np.searchsorted(tc, c, side="right")
        if lo == hi:
            continue
        p = tp[lo:hi]
        cen = regions.centers[m]
        k = np.clip(np.searchsorted(p, cen), 1, len(p) - 1) if len(p) > 1 else np.zeros(len(m), int)
        if len(p) > 1:
            left, right = p[k - 1], p[k]
            k = np.where(np.abs(cen - left) <= np.abs(right - cen), k - 1, k)
        t = p[k]
        s, e = regions.starts[m], regions.ends[m]
        d = np.where((s <= t) & (t < e), 0, np.minimum(np.abs(t - s), np.abs(t - e)))
        hit = d <= maxd
        out[m[hit]] = order[lo + k[hit]]
    return out
