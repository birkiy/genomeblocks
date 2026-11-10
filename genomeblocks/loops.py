"""HiChIP loops: FitHiChIP's peak-to-all model on a cooler, in numpy.

:func:`fithichip` calls the loops FitHiChIP 11.0 calls with ``IntType=3``
(peak-to-all), ``BiasType=1`` (coverage bias), ``UseP2PBackgrnd=0`` and
``MergeInt=1`` — its defaults otherwise — and reads the pixels straight from
a cooler instead of writing, filtering and sorting text files (FitHiChIP took
22 min to 19 h on the samples of ``benchmarks/bench_loops.py``; this takes
3 s to 1.5 min). The model, step by step:

1. **Coverage** of a 5 kb bin = the contacts of every pixel whose first bin it
   is (cis and trans) plus the cis pixels whose second bin it is; the diagonal
   counts twice. Bins shorter than 2 bp get 0, and a chromosome has bins only
   if it is the first chromosome of some pixel.
2. **Peak bins**: a bin [start, end] touching a peak [start, end] (closed
   intervals, as GRanges in FitHiChIP).
3. **Bias** = coverage / mean nonzero coverage of the bins of the same class
   (peak bins, other bins).
4. **Tested pairs**: cis pixels with a count, start distance in
   [``low``, ``high``], at least one peak bin. Like FitHiChIP, the first bin of
   the coverage table sorted by name (``sort -k1,1 -k2,2n``) has no features
   and its pairs are not tested.
5. **Distance model**: the pairs sorted by distance are cut into ``nbins``
   groups holding equal numbers of contacts (whole distances); each group's
   prior = contacts / possible peak-to-all pairs at those distances (zero
   contacts included) / all contacts, against its mean distance. R
   ``smooth.spline(cv=TRUE)`` (:mod:`genomeblocks._smoothspline`) fits it, is
   evaluated on the bin-size grid and fitted once more: p(d).
6. **Bias regression** in each group, on pairs with both biases in
   ``bias_range``: log10 count = a + b1 log10 bias1 + b2 log10 bias2; a, b1, b2
   are smoothed over distance with the same two splines.
7. **Test**: expected count E = 10^(a(d) + b1(d) log10 bias1 + b2(d) log10
   bias2), or N p(d) when a bias is 0; P = P(X >= count) for
   X ~ Binomial(floor(sum E), E / floor(sum E)); Benjamini-Hochberg q; loops:
   q < ``qvalue``.
8. **Merging** (FitHiChIP's ``CombineNearbyInteraction.py --percent 100
   --Neigh 2``): loops touching in 8-connectivity form components; in each,
   loops are taken by q, then count, unless both bins are within ``neigh`` bins
   of an already taken loop.

On a synthetic cooler and five HiChIP samples (0.3–52 M tested pairs) this
gives FitHiChIP's tested pairs, loops and merged loops; P and q agree to 1e-7
(relative) where FitHiChIP's table of tested pairs is at hand — FitHiChIP
itself moves them by as much between two machines — and to 1e-5 on the
synthetic cooler (``benchmarks/bench_loops.py``).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from .bedpe import Pairs, open_cooler
from .genome import Genome, read_sizes
from .loci import Loci


@dataclass
class FitHiChIP:
    """The result of :func:`fithichip`.

    ``tested``, ``loops`` and ``merged`` are :class:`~genomeblocks.Pairs` (bin
    pairs as anchors a / b). ``tested`` carries ``count``, ``bias1``, ``bias2``,
    ``is_peak1``, ``is_peak2``, ``distance``, ``p`` (distance prior),
    ``expected``, ``p_bias``, ``dbinom``, ``P`` and ``Q``; ``loops`` is its
    rows with Q < ``qvalue``; ``merged`` the merged loops with ``count``,
    ``P``, ``Q``, the component's bin span (``span1_start`` ...), ``sum_count``
    and ``strong_conn`` (the fraction of the span's bin pairs that are loops).
    ``groups`` is the equal-occupancy table (FitHiChIP's ``Bin_Info.log``),
    ``coefficients`` the regression per group."""
    tested: Pairs
    loops: Pairs
    merged: Pairs
    groups: "object"
    coefficients: "object"
    binsize: int
    qvalue: float
    timings: dict = field(default_factory=dict)

    def __repr__(self):
        return (f"FitHiChIP(tested={len(self.tested):,}, loops={len(self.loops):,} (q < {self.qvalue}), "
                f"merged={len(self.merged):,}, binsize={self.binsize})")

    def write(self, outdir: str, prefix: str, *, tested: bool = False) -> dict:
        """FitHiChIP's browser files — ``{prefix}.interactions_FitHiC_Q{q}_IGV.bedpe`` and
        ``Merge_Nearby_Interactions/{prefix}.interactions_FitHiC_Q{q}_MergeNearContacts_IGV.bedpe``, byte for
        byte in FitHiChIP's format — plus the loop tables as TSV and, with ``tested=True``, every tested pair as
        parquet. Returns the paths."""
        out = Path(outdir)
        (out / "Merge_Nearby_Interactions").mkdir(parents=True, exist_ok=True)
        stem = f"{prefix}.interactions_FitHiC_Q{self.qvalue}"
        paths = {"igv": out / f"{stem}_IGV.bedpe",
                 "merged_igv": out / "Merge_Nearby_Interactions" / f"{stem}_MergeNearContacts_IGV.bedpe",
                 "loops": out / f"{stem}.loops.tsv",
                 "merged": out / "Merge_Nearby_Interactions" / f"{stem}_MergeNearContacts.tsv"}
        _write_igv(self.loops, paths["igv"])
        _write_igv(self.merged, paths["merged_igv"])
        self.loops.to_pandas().to_csv(paths["loops"], sep="\t", index=False)
        self.merged.to_pandas().to_csv(paths["merged"], sep="\t", index=False)
        if tested:
            paths["tested"] = out / f"{prefix}.interactions_FitHiC.parquet"
            self.tested.save(str(paths["tested"]))
        return {k: str(v) for k, v in paths.items()}


def fithichip(cool: str, peaks, *, resolution: Optional[int] = None, chromsizes=None, low: int = 5_000,
              high: int = 10_000_000, qvalue: float = 0.01, nbins: int = 200,
              bias_range: Tuple[float, float] = (0.2, 5.0), neigh: int = 2, verbose: bool = True) -> FitHiChIP:
    """FitHiChIP peak-to-all loops (coverage bias regression) from a cooler and anchor peaks.

    Args:
        cool: a .cool, a cooler URI, or an .mcool with ``resolution``. Raw counts (``pixels/count``) are used;
            any balancing weights are ignored, as in FitHiChIP's coverage-bias mode.
        peaks: the anchors, anything :func:`~genomeblocks.as_loci` reads (narrowPeak, BED, a Loci, a frame).
        chromsizes: the chromosomes to use (FitHiChIP's ``ChrSizeFile``); default: the cooler's.
        low, high: distance range of the tested pairs (bin starts), inclusive.
        qvalue: FDR threshold of the loops (strict ``<``).
        nbins: number of equal-occupancy groups of the distance model.
        bias_range: biases of the pairs used to fit the bias regression.
        neigh: merging neighbourhood in bins.
    """
    from scipy import stats

    T = {}
    t0 = time.perf_counter()
    clr = open_cooler(cool, resolution, caller="fithichip")
    binsize = int(clr.binsize)
    names = list(clr.chromnames)
    g = Genome.from_sizes(clr.chromsizes)
    with clr.open("r") as h5:
        bchromi = h5["bins/chrom"][:]
        bstart = h5["bins/start"][:].astype(np.int64)
        bend = h5["bins/end"][:].astype(np.int64)
        p1 = h5["pixels/bin1_id"][:]
        p2 = h5["pixels/bin2_id"][:]
        cnt = h5["pixels/count"][:]
    sizes = read_sizes(chromsizes) if chromsizes is not None else dict(clr.chromsizes)
    T["read"] = time.perf_counter() - t0

    # ── coverage, peak bins, bias ──
    t0 = time.perf_counter()
    nb = len(bstart)
    bcode = np.array([g.code[c] for c in names], np.int32)[bchromi]
    cis = bchromi[p1] == bchromi[p2]
    cov = np.bincount(p1, weights=cnt, minlength=nb) + np.bincount(p2[cis], weights=cnt[cis], minlength=nb)
    cov[(bend - bstart) < 2] = 0
    first_end = np.zeros(len(names), bool)
    first_end[np.unique(bchromi[p1])] = True
    used = np.array([c in sizes for c in names]) & first_end
    has_row = used[bchromi]                                   # the bins of FitHiChIP's coverage table
    bins = Loci(bcode, bstart, bend, genome=g)
    from .interop import as_loci
    pk = as_loci(peaks, genome=g)
    pk = pk.take(np.isin(pk.chroms, [c for c in names if c in sizes]))
    is_peak = bins.overlap_any(pk.slop(1)) & has_row         # closed intervals: touching counts
    bias = np.zeros(nb)
    for cls in (True, False):
        sel = has_row & (is_peak == cls)
        nz = sel & (cov > 0)
        if nz.any():
            bias[sel] = cov[sel] / cov[nz].mean()
    feat = has_row.copy()
    first_chrom = min((c for c, u in zip(names, used) if u), key=lambda s: s.encode(), default=None)
    if first_chrom is not None:
        rows = np.flatnonzero(bchromi == names.index(first_chrom))
        feat[rows[np.argmin(bstart[rows])]] = False
    T["coverage"] = time.perf_counter() - t0

    # ── tested pairs, sorted by distance ──
    t0 = time.perf_counter()
    dist = np.abs(bstart[p2] - bstart[p1])
    sel = cis & (cnt > 0) & feat[p1] & feat[p2] & (dist >= low) & (dist <= high) & (is_peak[p1] | is_peak[p2])
    i1, i2, cc, d = p1[sel], p2[sel], cnt[sel].astype(np.int64), dist[sel]
    del p1, p2, cnt, cis, dist, sel
    o = np.argsort(d, kind="stable")
    i1, i2, cc, d = i1[o], i2[o], cc[o], d[o]
    n = len(d)
    if n == 0:
        raise ValueError("no peak-to-all pairs with contacts in the distance range")
    T["pairs"] = time.perf_counter() - t0

    # ── distance model ──
    t0 = time.perf_counter()
    ud, ufirst, ucount = np.unique(d, return_index=True, return_counts=True)
    ucc = np.add.reduceat(cc, ufirst)
    possible = _possible_pairs(bchromi[has_row], bstart[has_row] // binsize, is_peak[has_row], ud // binsize)
    total = int(cc.sum())
    per_group = np.floor(total / (nbins if total >= nbins else n))
    groups = []
    acc = nel = nall = tsum = 0.0
    for k in range(len(ud)):
        nel += ucount[k]
        acc += ucc[k]
        nall += possible[k]
        tsum += possible[k] * ud[k]
        if acc >= per_group:
            groups.append((int(tsum / nall) if nall else np.nan, nel, nall, acc))
            acc = nel = nall = tsum = 0.0
    import pandas as pd
    gr = pd.DataFrame(groups, columns=["AvgIntDist", "NumLociNonZeroContact", "NumLoci_AllPossible", "NumContact"])
    gr["AvgContact"] = gr["NumContact"] / gr["NumLoci_AllPossible"]
    gr["PriorProb"] = gr["AvgContact"] / total
    grid = np.arange(d.min(), d.max() + 1, binsize, dtype=float)
    p_u = _spline_twice(gr["AvgIntDist"].to_numpy(float), gr["PriorProb"].to_numpy(float), grid, ud)
    T["distance model"] = time.perf_counter() - t0

    # ── bias regression ──
    t0 = time.perf_counter()
    b1, b2 = bias[i1], bias[i2]
    lo, hi = bias_range
    ends = np.cumsum(gr["NumLociNonZeroContact"].to_numpy().astype(np.int64))
    coef = []
    for s, e in zip(np.r_[0, ends[:-1]], ends):
        x1, x2, c = b1[s:e], b2[s:e], cc[s:e]
        f = (x1 >= lo) & (x1 <= hi) & (x2 >= lo) & (x2 <= hi)
        coef.append(_lm2(np.log10(c[f].astype(float)), np.log10(x1[f]), np.log10(x2[f])))
    co = pd.DataFrame(coef, columns=["Intercept", "LogBias1", "LogBias2"])
    co.insert(0, "Dist", gr["AvgIntDist"].to_numpy())
    cu = {c: _spline_twice(co["Dist"].to_numpy(float), co[c].to_numpy(float), grid, ud)
          for c in ("Intercept", "LogBias1", "LogBias2")}
    T["bias regression"] = time.perf_counter() - t0

    # ── expected counts, binomial test, BH ──
    t0 = time.perf_counter()
    grp = np.repeat(np.arange(len(ud)), ucount)
    p_d = p_u[grp]
    pos = (b1 > 0) & (b2 > 0)
    E = np.empty(n)
    E[~pos] = total * p_d[~pos]
    gp = grp[pos]
    E[pos] = 10 ** (cu["Intercept"][gp] + cu["LogBias1"][gp] * np.log10(b1[pos]) + cu["LogBias2"][gp] * np.log10(b2[pos]))
    E[~np.isfinite(E)] = 2.0 ** 15
    s_e = E.sum()
    ntot = int(s_e) if s_e < 2 ** 31 else 2 ** 31 - 1
    pr = E / ntot
    db = stats.binom.pmf(cc, ntot, pr)
    P = db + stats.binom.sf(cc, ntot, pr)
    Q = _bh(P)
    T["test"] = time.perf_counter() - t0

    # ── tables, merging ──
    t0 = time.perf_counter()
    tested = Pairs(Loci(bcode[i1], bstart[i1], bend[i1], genome=g), Loci(bcode[i2], bstart[i2], bend[i2], genome=g),
                   {"count": cc, "bias1": b1, "bias2": b2, "is_peak1": is_peak[i1], "is_peak2": is_peak[i2],
                    "distance": d, "p": p_d, "expected": E, "p_bias": pr, "dbinom": db, "P": P, "Q": Q})
    loops = tested.take(np.flatnonzero(Q < qvalue))
    merged = _merge(loops, binsize, neigh)
    T["merge"] = time.perf_counter() - t0
    res = FitHiChIP(tested, loops, merged, gr, co, binsize, qvalue, T)
    if verbose:
        print(f"[fithichip] {n:,} peak-to-all pairs, {len(loops):,} loops (q < {qvalue}), {len(merged):,} merged | "
              + ", ".join(f"{k} {v:.1f}s" for k, v in T.items()))
    return res


# ── pieces ─────────────────────────────────────────────────────────────────

def _spline_twice(x, y, grid, at):
    """FitHiChIP's curve: smooth.spline(cv=TRUE) of the finite (x, y), evaluated on ``grid``, fitted again,
    evaluated at ``at``."""
    from ._smoothspline import smooth_spline
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= 2:
        raise ValueError("fewer than 3 finite points for the spline fit")
    s1 = smooth_spline(x[ok], y[ok])
    return smooth_spline(grid, s1.predict(grid)).predict(at)


def _possible_pairs(chrom, binidx, is_peak, shifts):
    """Peak-to-all locus pairs (zero contacts included) at each bin shift, over all chromosomes: peak-peak
    pairs once, peak-other pairs in both directions (FitHiChIP's ``AllPossibleContacts`` for IntType 3)."""
    shifts = np.asarray(shifts, np.int64)
    pad = int(shifts.max()) + 1 if len(shifts) else 1
    o = np.lexsort((binidx, chrom))
    chrom, binidx, is_peak = chrom[o], binidx[o], is_peak[o]
    starts = np.r_[0, np.flatnonzero(np.diff(chrom)) + 1]
    lens = np.diff(np.r_[starts, len(chrom)])
    span = np.maximum.reduceat(binidx, starts) + 1
    base = np.r_[0, np.cumsum(span + pad)[:-1]]
    pos = np.repeat(base, lens) + binidx                    # chromosomes laid end to end, > max shift apart
    P = np.zeros(int(base[-1] + span[-1] + pad), bool)
    N = np.zeros_like(P)
    P[pos[is_peak]] = True
    N[pos[~is_peak]] = True
    out = np.empty(len(shifts), np.int64)
    for k, s in enumerate(shifts):
        if s == 0:
            out[k] = np.count_nonzero(P)
        else:
            out[k] = (np.count_nonzero(P[:-s] & P[s:]) + np.count_nonzero(P[:-s] & N[s:])
                      + np.count_nonzero(N[:-s] & P[s:]))
    return out


def _lm2(y, x1, x2):
    """Coefficients of R ``lm(y ~ x1 + x2)``; NA for a column that is aliased (as R's pivoting, tol 1e-7)."""
    X = np.column_stack([np.ones(len(y)), x1, x2])
    keep = []
    for j in range(3):
        cols = keep + [j]
        sub = X[:, cols]
        if len(y) >= len(cols) and np.linalg.matrix_rank(sub, tol=1e-7 * max(np.abs(sub).max(), 1e-300)) == len(cols):
            keep.append(j)
    beta = np.full(3, np.nan)
    if keep:
        beta[keep] = np.linalg.lstsq(X[:, keep], y, rcond=None)[0]
    return beta


def _bh(p):
    """R ``p.adjust(p, "BH")``."""
    n = len(p)
    o = np.argsort(-p, kind="stable")
    q = np.minimum.accumulate((n / np.arange(n, 0, -1)) * p[o])
    out = np.empty(n)
    out[o] = np.minimum(1.0, q)
    return out


def _q15(q):
    """The values as FitHiChIP's text tables hold them (15 significant digits)."""
    return np.array([float(f"{v:.15g}") for v in q])


def _merge(loops: Pairs, binsize: int, neigh: int = 2) -> Pairs:
    """FitHiChIP ``CombineNearbyInteraction.py --percent 100 --Neigh neigh`` (8-connectivity)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    g = loops.genome
    if len(loops) == 0:
        return Pairs(Loci(genome=g), Loci(genome=g), {k: np.zeros(0) for k in
                     ("count", "P", "Q", "span1_start", "span1_end", "span2_start", "span2_end", "sum_count",
                      "strong_conn")})
    ch = loops.a.codes.astype(np.int64)
    cc = loops.cols["count"]
    # FitHiChIP keys a bin by end // binsize and keeps the first loop of a key in its file order
    o = np.lexsort((-cc, loops.b.starts, loops.a.starts, ch))
    k1 = loops.a.ends[o] // binsize
    k2 = loops.b.ends[o] // binsize
    key = (ch[o] << 42) | (k1 << 21) | k2
    key, first = np.unique(key, return_index=True)
    rows = o[first]
    ch, k1, k2 = ch[rows], k1[first], k2[first]
    cc, P, Q = cc[rows], loops.cols["P"][rows], loops.cols["Q"][rows]
    n = len(rows)
    src, dst = [], []
    for d1 in (-1, 0, 1):
        for d2 in (-1, 0, 1):
            if d1 or d2:
                nk = (ch << 42) | ((k1 + d1) << 21) | (k2 + d2)
                j = np.minimum(np.searchsorted(key, nk), n - 1)
                hit = (key[j] == nk) & (k1 + d1 >= 0) & (k2 + d2 >= 0)
                src.append(np.flatnonzero(hit))
                dst.append(j[hit])
    src, dst = np.concatenate(src), np.concatenate(dst)
    _, comp = connected_components(coo_matrix((np.ones(len(src)), (src, dst)), shape=(n, n)), directed=False)
    taken, keep = set(), np.zeros(n, bool)
    r = range(-neigh, neigh + 1)
    for i in np.lexsort((k2, k1, -cc, _q15(Q), comp)):
        c, a, b = comp[i], k1[i], k2[i]
        if not any((c, a + x, b + y) in taken for x in r for y in r):
            taken.add((c, a, b))
            keep[i] = True
    # component span, summed counts, fraction of the span's bin pairs that are loops
    lo1 = np.full(comp.max() + 1, np.iinfo(np.int64).max)
    lo2 = lo1.copy()
    hi1 = np.full(comp.max() + 1, -1)
    hi2 = hi1.copy()
    np.minimum.at(lo1, comp, k1)
    np.minimum.at(lo2, comp, k2)
    np.maximum.at(hi1, comp, k1)
    np.maximum.at(hi2, comp, k2)
    sumc = np.bincount(comp, weights=cc)
    nrow = hi1 - lo1 + 1
    rc = np.repeat(np.arange(len(lo1)), nrow)
    rb = np.repeat(lo1 - np.r_[0, np.cumsum(nrow)[:-1]], nrow) + np.arange(nrow.sum())
    cch = np.zeros(len(lo1), np.int64)
    cch[comp] = ch
    lo_key = (cch[rc] << 42) | (rb << 21) | lo2[rc]
    hi_key = (cch[rc] << 42) | (rb << 21) | hi2[rc]
    inbox = np.searchsorted(key, hi_key, "right") - np.searchsorted(key, lo_key, "left")
    strong = np.add.reduceat(inbox, np.r_[0, np.cumsum(nrow)[:-1]]) / (nrow * (hi2 - lo2 + 1))
    s = np.flatnonzero(keep)
    cs = comp[s]
    a = Loci(ch[s].astype(np.int32), (k1[s] - 1) * binsize, k1[s] * binsize, genome=g)
    b = Loci(ch[s].astype(np.int32), (k2[s] - 1) * binsize, k2[s] * binsize, genome=g)
    return Pairs(a, b, {"count": cc[s], "P": P[s], "Q": Q[s], "span1_start": (lo1[cs] - 1) * binsize,
                        "span1_end": hi1[cs] * binsize, "span2_start": (lo2[cs] - 1) * binsize,
                        "span2_end": hi2[cs] * binsize, "sum_count": sumc[cs].astype(np.int64),
                        "strong_conn": strong[cs]})


def _awk_num(v: float) -> str:
    """A number as awk prints it (integers as integers, else %.6g)."""
    return str(int(v)) if float(v).is_integer() and abs(v) < 2 ** 53 else f"{v:.6g}"


def _write_igv(P: Pairs, path) -> None:
    """FitHiChIP's IGV bedpe: bin midpoints +- 1, '-', min(500, -log10 q) (500 for q = 0), '.', '.', sorted
    ``sort -k1,1 -k2,2n -k5,5n`` (chromosome names in byte order)."""
    q = _q15(P.cols["Q"])
    score = [min(500.0, -math.log(v) / math.log(10)) if v > 0 else 500 for v in q]     # awk: C's log
    c1, c2 = P.a.chroms.astype(str), P.b.chroms.astype(str)
    m1 = (P.a.starts + P.a.ends) / 2
    m2 = (P.b.starts + P.b.ends) / 2
    order = sorted(range(len(P)), key=lambda i: (c1[i].encode(), m1[i], m2[i]))
    with open(path, "w") as f:
        for i in order:
            f.write(f"{c1[i]}\t{_awk_num(m1[i] - 1)}\t{_awk_num(m1[i] + 1)}\t{c2[i]}\t{_awk_num(m2[i] - 1)}\t"
                    f"{_awk_num(m2[i] + 1)}\t-\t{_awk_num(score[i])}\t.\t.\n")
