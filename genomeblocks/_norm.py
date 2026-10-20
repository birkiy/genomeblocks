"""Count normalisation shared by :class:`~genomeblocks.Expression` and the
signal cubes: library sizes, CPM, edgeR's TMM factors and ``filterByExpr``
rule, DESeq2's median-of-ratios size factors.

Pure numpy, feature rows x sample columns. Every engine would give the same
numbers, so none of this is a backend; the rules follow the R sources
(edgeR ``calcNormFactors`` / ``cpm`` / ``filterByExpr``, DESeq2
``estimateSizeFactorsForMatrix``) closely enough that the test suite pins
them against values computed in R.
"""
from __future__ import annotations

from typing import Optional

import numpy as np


def _matrix(counts) -> np.ndarray:
    x = np.asarray(counts, dtype=float)
    if x.ndim != 2:
        raise ValueError(f"expected a features x samples matrix, got shape {x.shape}")
    return x


def lib_sizes(counts) -> np.ndarray:
    """Column totals, ignoring NaN rows (features that were not quantified)."""
    return np.nansum(_matrix(counts), axis=0)


def tmm_factors(counts, *, trim_lfc: float = 0.3, trim_mag: float = 0.05, ref: Optional[int] = None,
                lib_size=None) -> np.ndarray:
    """edgeR's trimmed mean of M-values (Robinson & Oshlack 2010), one factor
    per sample, scaled to a geometric mean of 1.

    The reference sample is the one whose upper quartile of counts per
    library size is closest to the mean upper quartile (``ref=`` picks one).
    Against it, every other sample's log-ratios M are trimmed by rank
    (``trim_lfc`` on M, ``trim_mag`` on the mean log-abundance A) and
    averaged with the asymptotic-variance weights of ``calcNormFactors``.
    Rows holding a NaN or only zeros take no part.
    """
    x = _matrix(counts)
    x = x[~np.isnan(x).any(axis=1)]
    lib = np.asarray(lib_size, float) if lib_size is not None else x.sum(axis=0)
    m = x.shape[1]
    if m == 0:
        return np.zeros(0)
    if not (lib > 0).all():
        raise ValueError("tmm_factors: every sample needs a positive library size")
    x = x[(x > 0).any(axis=1)]
    if ref is None:
        if x.shape[0] == 0:
            return np.ones(m)
        f75 = np.percentile(x / lib, 75, axis=0)
        if np.median(f75) < 1e-20:
            ref = int(np.argmax(np.sqrt(x).sum(axis=0)))
        else:
            ref = int(np.argmin(np.abs(f75 - f75.mean())))
    from scipy.stats import rankdata
    out = np.ones(m)
    r, nr = x[:, ref], lib[ref]
    for j in range(m):
        if j == ref:
            continue
        o, no = x[:, j], lib[j]
        with np.errstate(divide="ignore", invalid="ignore"):
            logr = np.log2((o / no) / (r / nr))
            abse = (np.log2(o / no) + np.log2(r / nr)) / 2
            v = (no - o) / no / o + (nr - r) / nr / r
        fin = np.isfinite(logr) & np.isfinite(abse) & (abse > -1e10)
        logr, abse, v = logr[fin], abse[fin], v[fin]
        n = len(logr)
        if n == 0 or np.max(np.abs(logr)) < 1e-6:
            continue
        lo_l = np.floor(n * trim_lfc) + 1
        hi_l = n + 1 - lo_l
        lo_s = np.floor(n * trim_mag) + 1
        hi_s = n + 1 - lo_s
        rl, rs = rankdata(logr), rankdata(abse)
        keep = (rl >= lo_l) & (rl <= hi_l) & (rs >= lo_s) & (rs <= hi_s)
        if not keep.any():
            continue
        f = np.sum(logr[keep] / v[keep]) / np.sum(1 / v[keep])
        out[j] = 2 ** f if np.isfinite(f) else 1.0
    return out / np.exp(np.mean(np.log(out)))


def size_factors(counts) -> np.ndarray:
    """DESeq2's median-of-ratios size factors (``estimateSizeFactorsForMatrix``):
    for each sample, the median over features of count / geometric mean of
    that feature across samples, over features with a finite geometric mean
    and a positive count in the sample."""
    x = _matrix(counts)
    x = x[~np.isnan(x).any(axis=1)]
    with np.errstate(divide="ignore"):
        logx = np.log(x)
    loggeo = logx.mean(axis=1)
    ok = np.isfinite(loggeo)
    if not ok.any():
        raise ValueError("size_factors: no feature has a positive count in every sample")
    out = np.empty(x.shape[1])
    for j in range(x.shape[1]):
        sel = ok & (x[:, j] > 0)
        out[j] = np.exp(np.median(logx[sel, j] - loggeo[sel]))
    return out


def cpm(counts, *, lib_size=None, factors=None, log: bool = False, prior: float = 2.0) -> np.ndarray:
    """Counts per million by (normalised) library size, as edgeR's ``cpm``.

    ``factors`` (TMM or size factors) scale the library sizes. With
    ``log=True`` a prior count (scaled by each library's size, edgeR's
    ``prior.count``) is added before the log2.
    """
    x = _matrix(counts)
    lib = np.asarray(lib_size, float) if lib_size is not None else np.nansum(x, axis=0)
    if factors is not None:
        lib = lib * np.asarray(factors, float)
    if not log:
        return x / lib * 1e6
    ps = lib / lib.mean() * prior
    return np.log2((x + ps) / (lib + 2 * ps) * 1e6)


def filter_by_expr(counts, *, group=None, lib_size=None, min_count: float = 10, min_total_count: float = 15,
                   large_n: int = 10, min_prop: float = 0.7) -> np.ndarray:
    """edgeR's ``filterByExpr`` rule as a boolean mask over the features.

    A feature is kept when its CPM is at least ``min_count / median library
    size`` in at least as many samples as the smallest group has (shrunk
    towards ``min_prop`` of that size above ``large_n``), and its total count
    is at least ``min_total_count``. Without ``group`` every sample must pass.
    NaN rows are never kept.
    """
    x = _matrix(counts)
    n = x.shape[1]
    lib = np.asarray(lib_size, float) if lib_size is not None else np.nansum(x, axis=0)
    if group is None:
        min_n = n
    else:
        _, cnt = np.unique(np.asarray(group).astype(str), return_counts=True)
        min_n = int(cnt.min())
    if min_n > large_n:
        min_n = large_n + (min_n - large_n) * min_prop
    cutoff = min_count / np.median(lib) * 1e6
    c = x / lib * 1e6
    tol = 1e-14
    with np.errstate(invalid="ignore"):
        keep = (np.sum(c >= cutoff, axis=1) >= min_n - tol) & (np.sum(x, axis=1) >= min_total_count - tol)
    return keep & ~np.isnan(x).any(axis=1)
