"""Super-enhancers, ROSE-style, on columnar tables.

    peaks ─▶ stitch (merge within ``stitch`` bp) ─▶ score = mean signal x width
          ─▶ slope-1 knee on the ranked scores ─▶ SEs

``knee`` is the same slope-1 rule as ``Architecture.elbow`` (prime hubs), so
SEs and hubs are cut the same way.
"""
from __future__ import annotations

from typing import Sequence, Tuple

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
    yn_s = uniform_filter1d(yn, size=max(11, m // 200), mode="nearest")
    crossed = np.flatnonzero(np.gradient(yn_s, xn) >= 1.0)
    i = int(crossed[0]) if len(crossed) else m
    return m - i, order


def call_se(peaks: Loci, bigwigs: Sequence[str], *, stitch: int = 12_500, workers: int = 1,
            verbose: bool = False) -> Loci:
    """ROSE with genomeblocks: stitch peaks, score mean signal (averaged over
    ``bigwigs``) x width, keep the stitched regions above the slope-1 knee.

    Returns the SEs in genome order with columns ``score`` and ``rank``
    (1 = strongest). ``all_regions`` on the result holds every stitched region
    with its score (handy for the hockey-stick plot).
    """
    bigwigs = [bigwigs] if isinstance(bigwigs, str) else list(bigwigs)
    st = peaks.slop(stitch // 2).merge().slop(-stitch // 2)
    sig = np.nan_to_num(st.signal(bigwigs, span=True, n_bins=1, workers=workers, progress=False,
                                  verbose=False)[:, :, 0]).mean(1)
    score = sig * st.lengths
    cut, order = knee(score)
    keep = order[:cut]
    rank = np.zeros(len(st), np.int64)
    rank[order] = np.arange(1, len(order) + 1)
    st.cols["score"], st.cols["rank"] = score, rank
    se = st.take(np.sort(keep))
    se.all_regions = st
    if verbose:
        print(f"[INFO] {len(peaks):,} peaks → {len(st):,} stitched regions → {len(se):,} SEs "
              f"(median {int(np.median(se.lengths)) if len(se) else 0:,} bp)")
    return se


def nearest_gene_within(regions: Loci, tss: Loci, maxd: int) -> np.ndarray:
    """For each region: the row in ``tss`` (1-bp TSS Loci) whose TSS is closest to
    the region centre, or -1 when that TSS lies more than ``maxd`` from the region
    (0 when it is inside it)."""
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
