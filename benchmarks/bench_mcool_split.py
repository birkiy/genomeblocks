#!/usr/bin/env python3
"""add_mcool: how a pixel's count is shared among the edges in it.

Two samples, A and B, each with its own graph (its CREs, its loops) and one union graph (the CREs of
both, the loops of both). The pixel counts are a sample's reads either way, so a rule that is fair to
the union should leave each sample's CRE strengths where its own graph puts them, apart from the
contact that the other sample's loops genuinely add. Rules compared, each scored on the super-hubs
(slope-1 knee on CRE strength, ``se.knee``) a sample loses between its own and the union graph, and
on contact conservation (summed edge weight / summed counts of the pixels touched):

- equal split, ``diagonal=True``    the 2.0.1 rule
- no split                          every edge gets the full pixel count
- each pixel once per CRE           a CRE takes the count of every distinct pixel it touches, once
- centre binning                    the equal split, CREs binned by their centre
- ATAC-weighted                     ``weights=`` the sample's accessibility, ``diagonal=True``
- peak mask                         ``weights=`` the sample's 0/1 peak mask, ``diagonal=True``
- equal split                       ``diagonal=False`` (the new default)
- peak mask, no diagonal            ``weights=`` the peak mask, ``diagonal=False``

The data are synthetic and seeded (one 100 Mb chromosome, 5 kb bins, written to a temporary
directory), so the bench runs anywhere in seconds; it also times ``add_mcool`` with and without
``weights=``. The rows with ``part = "tcell"`` are the same comparison measured on H3K27ac HiChIP of
naive CD4, Th17 and Treg cells (5 kb coolers, FitHiChIP loops, ATAC-seq CREs; own graph = the state's
CREs and loops, union = the CREs and loops of the three states). Those data are not part of this
repository; the numbers are recorded here as measured.

    python bench_mcool_split.py        # writes results/mcool_split.json
"""
from __future__ import annotations

import contextlib
import io
import tempfile

import numpy as np
import pandas as pd

from common import Recorder, timeit

from genomeblocks import Architecture, Loci, Pairs, as_loci
from genomeblocks.se import knee

BIN, NBINS, SEED = 5000, 20_000, 7
LOOP, OFF_LOOP, SELF = 2.0, 0.6, 600         # count scale of a loop pixel, any other pixel, a bin's self-contact


# ── synthetic samples ───────────────────────────────────────────────────────
def synthesize(rng):
    """CREs (some sharing a bin, some crossing a bin edge), per-sample peaks and accessibility, loops."""
    cre_bins = np.sort(rng.choice(np.arange(1, NBINS - 1), 12_000, replace=False))
    per_bin = rng.choice([1, 2, 3], size=len(cre_bins), p=[0.86, 0.11, 0.03])
    b = np.repeat(cre_bins, per_bin)
    s = b * BIN + rng.integers(0, BIN, len(b))
    ln = np.clip(rng.lognormal(np.log(900), 0.4, len(b)), 200, 3000).astype(np.int64)
    cre = as_loci(pd.DataFrame({"chrom": "chr1", "start": s, "end": s + ln})).sort().merge()
    n = len(cre)
    peak = {"A": rng.random(n) < 0.72, "B": rng.random(n) < 0.72}
    none = ~(peak["A"] | peak["B"])
    pick_a = rng.random(n) < 0.5
    peak["A"] |= none & pick_a
    peak["B"] |= none & ~pick_a
    atac = {k: rng.lognormal(0, 0.6, n) * np.where(p, 1.0, 0.2) for k, p in peak.items()}
    activity = np.zeros(NBINS)
    activity[cre_bins] = rng.lognormal(0, 1.0, len(cre_bins))

    def loops(n_loops):
        i = rng.choice(cre_bins, n_loops, p=activity[cre_bins] / activity[cre_bins].sum())
        d = np.minimum(np.floor(rng.pareto(0.6, n_loops) + 1), 400).astype(np.int64)
        j = cre_bins[np.minimum(np.searchsorted(cre_bins, i + d), len(cre_bins) - 1)]
        k = np.unique(np.minimum(i, j) * NBINS + np.maximum(i, j))
        k = k[k // NBINS != k % NBINS]
        return np.c_[k // NBINS, k % NBINS]

    shared = loops(9000)
    lp = {"A": np.unique(np.r_[shared, loops(6000)], axis=0), "B": np.unique(np.r_[shared, loops(6000)], axis=0)}
    return cre, peak, atac, activity, lp


def as_pairs(bins):
    mid = bins * BIN + BIN // 2
    return Pairs.from_frame(pd.DataFrame({"chrom1": "chr1", "start1": mid[:, 0] - 1, "end1": mid[:, 0] + 1,
                                          "chrom2": "chr1", "start2": mid[:, 1] - 1, "end2": mid[:, 1] + 1}))


def write_cooler(path, keys, counts):
    import cooler
    bins = pd.DataFrame({"chrom": "chr1", "start": np.arange(NBINS) * BIN, "end": np.arange(1, NBINS + 1) * BIN})
    keep = counts > 0
    cooler.create_cooler(path, bins, pd.DataFrame({"bin1_id": keys[keep] // NBINS, "bin2_id": keys[keep] % NBINS,
                                                   "count": counts[keep]}))


def pixel_counts(rng, union, loops, activity):
    """Counts of every pixel the union skeleton touches, per sample: loop pixels of the sample are enriched,
    a bin's self-contact (the diagonal) is far larger than a loop pixel."""
    b = union.loci.starts // BIN
    k = np.unique(np.minimum(b[union.src], b[union.tgt]) * NBINS + np.maximum(b[union.src], b[union.tgt]))
    i, j = k // NBINS, k % NBINS
    out = {}
    for smp, lp in loops.items():
        is_loop = np.isin(k, lp[:, 0] * NBINS + lp[:, 1])
        lam = 40 * np.sqrt(activity[i] * activity[j]) * np.maximum(j - i, 1) ** -0.3 * np.where(is_loop, LOOP, OFF_LOOP)
        lam[i == j] = SELF * activity[i[i == j]]
        out[smp] = rng.poisson(lam)
    return k, out


# ── rules ───────────────────────────────────────────────────────────────────
def quiet(fn):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn()


def strength(A, w):
    n = len(A.loci)
    return np.bincount(A.src, w, n) + np.bincount(A.tgt, w, n)


def rule_strength(rule, A, cool, keys, counts, peak, atac):
    """(CRE strength, contact conservation) of graph A under one rule; peak / atac are per-row of A.loci."""
    b = A.loci.starts // BIN
    e = np.minimum(b[A.src], b[A.tgt]) * NBINS + np.maximum(b[A.src], b[A.tgt])
    c = np.zeros(len(e))
    hit = np.searchsorted(keys, e)
    ok = keys[np.minimum(hit, len(keys) - 1)] == e
    c[ok] = counts[hit[ok]]
    diag = b[A.src] == b[A.tgt]
    if rule == "no split":
        w = c
    elif rule == "each pixel once per CRE":
        n, npx = len(A.loci), int(e.max()) + 1
        src, tgt = A.src.astype(np.int64), A.tgt.astype(np.int64)
        pairs = np.unique(np.r_[src * npx + e, tgt * npx + e])
        hit = np.searchsorted(keys, pairs % npx)
        ok = keys[np.minimum(hit, len(keys) - 1)] == pairs % npx
        S = np.bincount(pairs[ok] // npx, counts[hit[ok]], n)
        return S, S.sum() / (2 * c_pixels(e, c, np.ones(len(e), bool)))
    elif rule == "centre binning":
        cen = A.loci.centers
        C = Architecture.from_edges(Loci(A.loci.codes, cen, cen + 1, genome=A.loci.genome), A.src, A.tgt)
        w = quiet(lambda: C.add_mcool(cool, diagonal=True)).ep["w"]
        bc = cen // BIN
        ec = np.minimum(bc[A.src], bc[A.tgt]) * NBINS + np.maximum(bc[A.src], bc[A.tgt])
        hit = np.searchsorted(keys, ec)
        okc = keys[np.minimum(hit, len(keys) - 1)] == ec
        cc = np.zeros(len(ec))
        cc[okc] = counts[hit[okc]]
        return strength(A, w), w.sum() / c_pixels(ec, cc, np.ones(len(ec), bool))
    else:
        kw = {"equal split, diagonal=True": dict(diagonal=True), "ATAC-weighted": dict(weights=atac, diagonal=True),
              "peak mask": dict(weights=peak, diagonal=True), "equal split, diagonal=False": {},
              "peak mask, no diagonal": dict(weights=peak)}[rule]
        w = quiet(lambda: A.copy().add_mcool(cool, **kw)).ep["w"]
        return strength(A, w), w.sum() / c_pixels(e, c, np.ones(len(e), bool) if kw.get("diagonal") else ~diag)
    return strength(A, w), w.sum() / c_pixels(e, c, np.ones(len(e), bool))


def c_pixels(e, c, use):
    """Summed count of the distinct pixels the edges in ``use`` touch."""
    _, first = np.unique(e[use], return_index=True)
    return c[use][first].sum()


def hubs(S):
    cut, order = knee(S)
    h = np.zeros(len(S), bool)
    h[order[:cut]] = True
    return h, float(S[order[cut - 1]])


RULES = ["equal split, diagonal=True", "no split", "each pixel once per CRE", "centre binning", "ATAC-weighted",
         "peak mask", "equal split, diagonal=False", "peak mask, no diagonal"]

# H3K27ac HiChIP of three T-cell states (see the docstring): own -> union super-hub cut-off, super-hubs, own
# super-hubs lost in the union, contact conservation (own, union), and the IKZF2 promoter CRE in Treg
# (strength own / union and whether it is a super-hub in the union).
TCELL = [  # rule, state, cut own, cut union, hubs own, hubs union, lost, conservation own, union
    ("equal split, diagonal=True", "Naive", 2961, 3137.5, 2073, 2341, 76, 1.0, 1.0),
    ("equal split, diagonal=True", "Th17", 3469, 3515, 1945, 2520, 43, 1.0, 1.0),
    ("equal split, diagonal=True", "Treg", 1930, 2191, 2096, 2379, 171, 1.0, 1.0),
    ("no split", "Naive", 3949, 3930, 1785, 2370, 0, 1.202, 1.267),
    ("no split", "Th17", 4380, 4541, 1796, 2431, 3, 1.187, 1.261),
    ("no split", "Treg", 2726, 2830, 1614, 2291, 2, 1.185, 1.262),
    ("each pixel once per CRE", "Naive", 3237, 3499, 2045, 2353, 68, 1.097, 1.126),
    ("each pixel once per CRE", "Th17", 3661, 3947, 2049, 2509, 51, 1.089, 1.122),
    ("each pixel once per CRE", "Treg", 2336, 2377, 1781, 2499, 0, 1.088, 1.123),
    ("centre binning", "Naive", 3118, 3106, 2085, 2536, 30, 1.0, 1.0),
    ("centre binning", "Th17", 3340, 3228, 2293, 3101, 31, 1.0, 1.0),
    ("centre binning", "Treg", 2140, 2206, 2006, 2584, 49, 1.0, 1.0),
    ("ATAC-weighted", "Naive", 2970, 3090, 2053, 2375, 31, 1.0, 1.0),
    ("ATAC-weighted", "Th17", 3670, 3762, 1777, 2259, 25, 1.0, 1.0),
    ("ATAC-weighted", "Treg", 2145, 2076, 1777, 2596, 18, 1.0, 1.0),
    ("peak mask", "Naive", 2961, 3123.5, 2073, 2350, 35, 1.0, 1.0),
    ("peak mask", "Th17", 3469, 3267, 1945, 2798, 0, 1.0, 1.0),
    ("peak mask", "Treg", 1930, 2185, 2096, 2405, 121, 1.0, 1.0),
    ("equal split, diagonal=False", "Naive", 2958.5, 3043, 1781, 2035, 61, 1.0, 1.0),
    ("equal split, diagonal=False", "Th17", 3116, 3069, 1892, 2504, 58, 1.0, 1.0),
    ("equal split, diagonal=False", "Treg", 1902, 2022, 1795, 2225, 81, 1.0, 1.0),
    ("peak mask, no diagonal", "Naive", 2958.5, 3046, 1781, 2058, 8, 1.0, 1.0),
    ("peak mask, no diagonal", "Th17", 3116, 3032, 1892, 2576, 0, 1.0, 1.0),
    ("peak mask, no diagonal", "Treg", 1902, 2012, 1795, 2273, 7, 1.0, 1.0),
]
TCELL_IKZF2 = {  # rule: (strength own, strength union, union super-hub)
    "equal split, diagonal=True": (2168, 2179.5, False), "no split": (4461, 4507, True),
    "each pixel once per CRE": (4336, 4359, True), "centre binning": (2168, 2179, False),
    "ATAC-weighted": (3245.6, 3262.8, True), "peak mask": (2168, 2179.5, False),
    "equal split, diagonal=False": (2168, 2179.5, True), "peak mask, no diagonal": (2168, 2179.5, True)}


if __name__ == "__main__":
    rec = Recorder("mcool_split")
    rng = np.random.default_rng(SEED)
    cre, peak, atac, activity, lp = synthesize(rng)
    union = Architecture.make(cre, as_pairs(np.r_[lp["A"], lp["B"]]), r=BIN // 2, verbose=False)
    keys, counts = pixel_counts(rng, union, lp, activity)
    with tempfile.TemporaryDirectory() as td:
        cool = {}
        for smp in ("A", "B"):
            cool[smp] = f"{td}/{smp}.cool"
            write_cooler(cool[smp], keys, counts[smp])
        b = union.loci.starts // BIN
        rec.add(part="synthetic", step="data", n_cres=len(cre), cres_sharing_a_bin=int((np.bincount(b)[b] > 1).sum()),
                cres_crossing_a_bin_edge=int((b != (union.loci.ends - 1) // BIN).sum()),
                union_edges=union.n_links, diagonal_edges=int((b[union.src] == b[union.tgt]).sum()),
                loops_A=len(lp["A"]), loops_B=len(lp["B"]))
        diag = b[union.src] == b[union.tgt]
        for smp in ("A", "B"):
            w = quiet(lambda: union.copy().add_mcool(cool[smp], diagonal=True)).ep["w"]
            rec.add(part="synthetic", step="self-contact", sample=smp, diagonal_edge_share=round(float(diag.mean()), 4),
                    diagonal_contact_share=round(float(w[diag].sum() / w.sum()), 4))
        for smp in ("A", "B"):
            mine = np.flatnonzero(peak[smp])
            own = Architecture.make(cre.take(mine), as_pairs(lp[smp]), r=BIN // 2, verbose=False)
            for rule in RULES:
                So, co = rule_strength(rule, own, cool[smp], keys, counts[smp], np.ones(len(mine)), atac[smp][mine])
                Su, cu = rule_strength(rule, union, cool[smp], keys, counts[smp], peak[smp].astype(float), atac[smp])
                So_u = np.zeros(len(cre))
                So_u[mine] = So                                       # own strengths on the union rows
                ho, cut_o = hubs(So_u)
                hu, cut_u = hubs(Su)
                both = So_u > 0
                closed_in_bin = np.bincount(b[~peak[smp]], minlength=NBINS)[b] > 0
                diluted = both & closed_in_bin                        # own CREs sharing their bin with a CRE without a peak
                rec.add(part="synthetic", step="rule", sample=smp, rule=rule, cutoff_own=cut_o, cutoff_union=cut_u,
                        hubs_own=int(ho.sum()), hubs_union=int(hu.sum()), lost=int((ho & ~hu).sum()),
                        own_cres_weaker_in_union=int((Su[both] < So_u[both] * (1 - 1e-9)).sum()),
                        own_cres_sharing_a_bin_with_a_closed_cre=int(diluted.sum()),
                        their_median_union_over_own=round(float(np.median(Su[diluted] / So_u[diluted])), 3),
                        conservation_own=round(float(co), 4), conservation_union=round(float(cu), 4))
        w_args = {"default (equal, diagonal=False)": {}, "weights= peak mask": {"weights": peak["A"].astype(float)},
                  "diagonal=True": {"diagonal": True}}
        for label, kw in w_args.items():
            t = timeit(lambda: quiet(lambda: union.copy().add_mcool(cool["A"], **kw)), repeat=5)
            rec.add(part="synthetic", step="time add_mcool", call=label, seconds=t["median"], runs=t["runs"],
                    n_edges=union.n_links)
    for rule, state, cut_o, cut_u, h_o, h_u, lost, c_o, c_u in TCELL:
        rec.add(part="tcell", step="rule", sample=state, rule=rule, cutoff_own=cut_o, cutoff_union=cut_u, hubs_own=h_o,
                hubs_union=h_u, lost=lost, conservation_own=c_o, conservation_union=c_u)
    for rule, (s_o, s_u, hub) in TCELL_IKZF2.items():
        rec.add(part="tcell", step="IKZF2 promoter (Treg)", rule=rule, strength_own=s_o, strength_union=s_u,
                union_super_hub=hub)
    rec.save(note="part 'tcell': recorded from H3K27ac HiChIP of naive CD4, Th17 and Treg cells (data not in the repository)")
