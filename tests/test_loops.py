"""HiChIP loops (FitHiChIP's model): R smooth.spline parity, the pieces FitHiChIP defines, planted loops."""
import numpy as np
import pandas as pd
import pytest

import genomeblocks as gb
from genomeblocks import Genome, Loci, Pairs
from genomeblocks import loops as L
from genomeblocks._smoothspline import smooth_spline

# R 4.1.2 and 4.5.3 give these to the last digit (smooth.spline(x, y, cv = TRUE); predict(fit, at)$y)
X1 = np.round(np.exp(np.linspace(np.log(5000), np.log(1e7), 60)))
Y1 = 1e-3 * (X1 / 5000) ** -1.1 * (1 + 0.1 * np.sin(np.arange(60)))
AT1 = [5000, 12345, 1e5, 3.3e6, 1e7, 2e7]
R1 = [0.0010139903366975117, 0.0003871459901856799, 3.3421836702286927e-05, 8.0230179320009181e-07,
      2.4869979577386692e-07, -1.922744662520739e-07]
X3 = np.repeat([1, 2, 3, 5, 8, 13, 21, 34, 55, 89], 3).astype(float)        # ties: R averages them, weight 3
Y3 = np.log(X3) + 0.3 * np.cos(np.arange(30))
AT3 = [0, 1, 4.5, 50, 89, 100]
R3 = [0.5286278523169583, 0.68708725564703976, 1.2349976120639479, 4.05800530910075, 4.3120145844744391,
      4.2904903139550923]


@pytest.mark.parametrize("x, y, at, expected, spar", [(X1, Y1, AT1, R1, 0.27193223216623275),
                                                      (X3, Y3, AT3, R3, 0.75144611643228898)])
def test_smooth_spline_matches_r(x, y, at, expected, spar):
    f = smooth_spline(x, y)
    assert f.spar == pytest.approx(spar, abs=1e-6)
    np.testing.assert_allclose(f.predict(at), expected, rtol=1e-7)          # incl. linear extrapolation


def test_smooth_spline_needs_four_distinct_x():
    with pytest.raises(ValueError, match="four unique"):
        smooth_spline([1, 1, 2, 3], [1, 2, 3, 4])


def test_bh_matches_r_p_adjust():
    p = np.array([0.01, 0.04, 0.03, 0.04, 0.2, 1e-5, 0.5, 0.03])
    r = [0.040000000000000001, 0.05333333333333333, 0.05333333333333333, 0.05333333333333333, 0.22857142857142856,
         8.0000000000000007e-05, 0.5, 0.05333333333333333]
    np.testing.assert_allclose(L._bh(p), r, rtol=1e-15)


def test_possible_pairs_counts_every_peak_to_all_pair():
    rng = np.random.default_rng(1)
    chrom = np.repeat([0, 1], [60, 40])
    idx = np.r_[np.sort(rng.choice(70, 60, replace=False)), np.arange(40)]     # chromosome 0 has holes
    peak = rng.random(100) < 0.2
    shifts = np.array([1, 2, 5, 17, 39, 69])
    got = L._possible_pairs(chrom, idx, peak, shifts)
    want = []
    for s in shifts:
        n = 0
        for c in (0, 1):
            on = dict(zip(idx[chrom == c], peak[chrom == c]))
            n += sum(1 for i, pi in on.items() if i + s in on and (pi or on[i + s]))
        want.append(n)
    assert got.tolist() == want


def _loops(genome, rows, binsize=5000):
    """Pairs of bins (b1, b2) on chr1 with count / P / Q columns."""
    b1, b2, cnt, q = (np.array(v) for v in zip(*rows))
    a = Loci(np.zeros(len(b1), np.int32), b1 * binsize, (b1 + 1) * binsize, genome=genome)
    b = Loci(np.zeros(len(b1), np.int32), b2 * binsize, (b2 + 1) * binsize, genome=genome)
    return Pairs(a, b, {"count": cnt, "P": q, "Q": q})


def test_merge_keeps_the_best_loop_of_each_neighbourhood():
    g = Genome(["chr1"])
    loops = _loops(g, [(10, 50, 9, 1e-9), (10, 51, 30, 1e-3), (10, 52, 30, 1e-3), (10, 53, 4, 1e-8),   # one chain
                       (20, 80, 5, 1e-4), (20, 81, 9, 1e-4)])                                         # a q tie
    m = L._merge(loops, 5000, neigh=2)
    kept = sorted(zip(m.a.starts // 5000, m.b.starts // 5000))
    assert kept == [(10, 50), (10, 53), (20, 81)]                  # 53 is 3 bins from 50; the tie goes to the count
    chain = m.b.starts // 5000 == 50
    assert m.cols["sum_count"][chain].tolist() == [73]
    assert m.cols["span2_start"][chain].tolist() == [50 * 5000] and m.cols["span2_end"][chain].tolist() == [54 * 5000]
    assert m.cols["strong_conn"][chain].tolist() == [1.0]


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """Two chromosomes of 5 kb bins: distance-decay background scaled by a per-bin coverage, peaks, and 3x3
    blobs of extra contacts between peak bins (the loops)."""
    cooler = pytest.importorskip("cooler")
    rng = np.random.default_rng(7)
    sizes = {"chrA": 3_000_000, "chrB": 2_000_002}                     # chrB ends in a 2 bp bin
    bins = pd.concat([pd.DataFrame({"chrom": c, "start": np.arange(0, n, 5000)}) for c, n in sizes.items()],
                     ignore_index=True)
    bins["end"] = np.minimum(bins.start + 5000, bins.chrom.map(sizes))
    nb = len(bins)
    chrom = bins.chrom.to_numpy()
    peak_bins = np.sort(rng.choice(nb, 120, replace=False))
    w = rng.lognormal(0, 0.3, nb)
    w[peak_bins] *= 3
    rows = []
    for c in sizes:
        idx = np.flatnonzero(chrom == c)
        for k in range(1, 400):
            i, j = idx[:-k], idx[k:]
            lam = 8.0 * (k ** -1.0) * w[i] * w[j]
            n = rng.poisson(lam)
            rows.append(np.column_stack([i[n > 0], j[n > 0], n[n > 0]]))
        rows.append(np.column_stack([idx, idx, rng.poisson(20, len(idx))]))
    loops = []
    pk = peak_bins[chrom[peak_bins] == "chrA"]
    for i in range(0, len(pk) - 8, 8):
        a, b = pk[i], pk[i + 5]
        if 40 <= b - a <= 300:
            loops.append((a, b))
            for da in (-1, 0, 1):
                for db in (-1, 0, 1):
                    rows.append(np.array([[a + da, b + db, 60 if da == db == 0 else 25]]))
    px = pd.DataFrame(np.vstack(rows), columns=["bin1_id", "bin2_id", "count"])
    px = px.groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()
    px = px[px["count"] > 0]
    d = tmp_path_factory.mktemp("loops")
    cooler.create_cooler(str(d / "t.cool"), bins, px)
    peaks = bins.iloc[peak_bins][["chrom", "start", "end"]].assign(start=lambda t: t.start + 1000,
                                                                     end=lambda t: t.start + 1500)
    peaks.to_csv(d / "peaks.bed", sep="\t", header=False, index=False)
    return str(d / "t.cool"), str(d / "peaks.bed"), bins, loops


def test_fithichip_finds_the_planted_loops(planted):
    cool, peaks, bins, planted_loops = planted
    r = gb.fithichip(cool, peaks, low=5000, high=2_000_000, verbose=False)
    assert len(planted_loops) >= 3
    centre = {(bins.start[a], bins.start[b]) for a, b in planted_loops}
    found = set(zip(r.loops.a.starts.tolist(), r.loops.b.starts.tolist()))
    merged = set(zip(r.merged.a.starts.tolist(), r.merged.b.starts.tolist()))
    assert centre <= found
    assert centre <= merged                                        # one loop per blob: its strongest pixel
    assert (r.loops.cols["Q"] < 0.01).all()


def test_fithichip_tested_pairs_follow_fithichip(planted):
    cool, peaks, bins, _ = planted
    r = gb.fithichip(cool, peaks, low=5000, high=2_000_000, verbose=False)
    t = r.tested
    chroms1 = t.a.chroms
    assert (t.cols["is_peak1"] | t.cols["is_peak2"]).all()
    assert ((t.cols["distance"] >= 5000) & (t.cols["distance"] <= 2_000_000)).all()
    assert np.all(np.diff(t.cols["distance"]) >= 0)
    # FitHiChIP's feature table lacks the first bin of `sort -k1,1 -k2,2n`: chrA:0-5000 is never tested
    assert not ((chroms1 == "chrA") & (t.a.starts == 0)).any()
    assert ((chroms1 == "chrB") & (t.a.starts == 0)).any() or not (t.b.chroms == "chrB").any()


def test_fithichip_writes_fithichip_igv_files(planted, tmp_path):
    cool, peaks, _, _ = planted
    r = gb.fithichip(cool, peaks, low=5000, high=2_000_000, verbose=False)
    paths = r.write(str(tmp_path), "s")
    lines = open(paths["igv"]).read().splitlines()
    assert len(lines) == len(r.loops)
    f = lines[0].split("\t")
    assert f[0] == "chrA" and int(f[2]) - int(f[1]) == 2 and f[6] == "-" and f[8:] == [".", "."]
    q = np.sort(r.loops.cols["Q"])[0]
    assert max(float(x.split("\t")[7]) for x in lines) == pytest.approx(min(500, -np.log10(q)), rel=1e-5)
    assert len(open(paths["merged_igv"]).read().splitlines()) == len(r.merged)

