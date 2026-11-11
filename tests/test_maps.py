"""MAPS loops: VGAM / R parity of the pieces, MAPS's quirks, and MAPS itself on a small synthetic pairs file."""
import numpy as np
import pytest

from genomeblocks import _maps as M
from genomeblocks._maps import maps, pospoisson_fit

B = 5000


def test_pospoisson_fit_matches_vgam():
    # R 4.5.3, VGAM 1.1-14: vglm(y ~ x1 + x2, family = pospoisson()); coef(fit); fit@iter; fitted(fit)[c(1, 2, 40)]
    i = np.arange(1, 41)
    X = np.column_stack([np.ones(40), np.sin(i), np.log(1 + i)])
    y = 1 + (i * i) % 7
    coef, eta, it = pospoisson_fit(X, y)
    np.testing.assert_allclose(coef, [1.0856051765245596, -0.16763126774845785, -0.0068249937708628914], rtol=1e-12)
    assert it == 4
    lam = np.exp(eta)
    np.testing.assert_allclose((-lam / np.expm1(-lam))[[0, 1, 39]],
                               [2.7740672082925819, 2.743539052671442, 2.7643729691771757], rtol=1e-13)


def test_pospoisson_fit_is_the_ztp_mle():
    rng = np.random.default_rng(3)
    X = np.column_stack([np.ones(3000), rng.normal(size=3000)])
    lam = np.exp(0.3 + 0.8 * X[:, 1])
    y = rng.poisson(lam)
    while (y == 0).any():                                     # zero-truncated draws
        z = y == 0
        y[z] = rng.poisson(lam[z])
    coef, eta, _ = pospoisson_fit(X, y)
    lam = np.exp(eta)
    score = X.T @ (y - lam - lam / np.expm1(lam))             # gradient of the ZTP log-likelihood
    assert np.abs(score).max() < 1e-3 * len(y) ** 0.5
    np.testing.assert_allclose(coef, [0.3, 0.8], atol=0.06)


def test_p_value_is_maps_formula():
    # R: ppois(c, E, lower.tail = FALSE) / ppois(0, E, lower.tail = FALSE) — P(Y > c), truncated mean as rate
    got = M._ztp_p(np.array([1, 3, 12, 80.]), np.array([0.5, 2.5, 10, 37.2]))
    np.testing.assert_allclose(got, [0.2292529587316009, 0.26410273990892624, 0.20845298735611048,
                                     3.5820806958078883e-10], rtol=1e-13)


def test_numbers_print_as_r_write_table():
    # R: options(scipen = 999); write.table(...) of the same values
    vals = [1e-20, 123456.7, 1, 2.5e-7, 1 / 3, 1.5e20, 0.000123456789012345678]
    want = ["0.00000000000000000001", "123456.7", "1", "0.00000025", "0.333333333333333", "150000000000000000000",
            "0.000123456789012346"]
    assert [M._r_num(v) for v in vals] == want


def test_peak_bins_follow_maps_including_the_first_peak_quirk():
    # bins floor(s/B) .. ceil(e/B) - 1; a peak ending on a bin edge does not reach the next bin
    got = M._peak_bins(np.array([1000, 20000]), np.array([6000, 25000]), B, 10)
    assert np.flatnonzero(got).tolist() == [0, 1, 4]
    # first peak spans exactly 4 bins (MAPS.py L110-117): the bin before it is dropped, even if covered
    got = M._peak_bins(np.array([50_100, 45_000]), np.array([69_000, 46_000]), B, 20)
    assert np.flatnonzero(got).tolist() == [10, 11, 12, 13]
    got = M._peak_bins(np.array([50_100, 45_000]), np.array([74_000, 46_000]), B, 20)   # 5 bins: no quirk
    assert np.flatnonzero(got).tolist() == [9, 10, 11, 12, 13, 14]


def test_clusters_are_connected_components_numbered_by_first_row():
    from genomeblocks import Genome, Loci, Pairs
    g = Genome(["chr1"])
    # bins (b1, b2): a chain 0-1-2 (0 and 2 are 6 bins apart, linked through 1), a pair, a singleton
    b = np.array([[100, 300], [500, 900], [103, 303], [106, 306], [502, 901], [800, 1000]])
    n = len(b)
    a = Loci(np.zeros(n, np.int32), b[:, 0] * B, b[:, 0] * B + B, genome=g)
    c = Loci(np.zeros(n, np.int32), b[:, 1] * B, b[:, 1] * B + B, genome=g)
    p = np.array([1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 0.0])
    T = Pairs(a, c, {"set": np.array(["AND"] * n, object), "peak1": np.ones(n), "peak2": np.ones(n),
                     "count": np.full(n, 20), "dist": b[:, 1] - b[:, 0], "expected": np.ones(n),
                     "ratio": np.full(n, 20.), "P": p, "fdr": p})
    S = M._label(T, np.arange(n), ["chr1"], B, 15_000)
    # singletons first, then clusters, each in row order
    assert (S.a.starts // B).tolist() == [800, 100, 500, 103, 106, 502]
    assert S.cols["label"].tolist() == [0, 1, 2, 1, 1, 2]
    assert S.cols["ClusterSize"].tolist() == [1, 3, 2, 3, 3, 2]
    assert S.cols["CountNei"].tolist() == [1, 2, 2, 3, 2, 2]
    # p = 0 counts as the smallest other p; a cluster's score is the sum of -log10 p
    np.testing.assert_allclose(S.cols["NegLog10P"], [9, 20, 15, 20, 20, 15])
    assert S.cols["ClusterType"].tolist() == ["Singleton", "BroadPeak", "SharpPeak", "BroadPeak", "BroadPeak",
                                             "SharpPeak"]


# ── MAPS on a synthetic pairs file ─────────────────────────────────────────

# planted loops: clusters (MAPS's classify_peaks stops with an error when every call is a singleton) and singletons
PLANTED = {"chr1": [(15, 45), (15, 46), (16, 45), (105, 165), (205, 255), (35, 72)],
           "chr2": [(25, 85), (25, 86), (305, 345), (305, 347)]}


def synthetic(d):
    """Two 2 Mb chromosomes: features, peaks (every 10th bin, plus MAPS's first-peak quirk on chr2), short pairs
    per bin (more at peaks), long pairs with a deterministic distance / short-count decay at every bin pair with a
    peak, +60 contacts at PLANTED (AND and XOR). Everything deterministic, so MAPS's own output can be pinned."""
    nb = 400
    bins = np.arange(nb)
    is_pk = bins % 10 == 5
    feats, peaks, lines = [], [], []
    for ci, c in enumerate(["chr1", "chr2"]):
        for b in bins:
            el = 0 if b < 2 else 1000 + (b * 37) % 3000
            feats.append(f"{c}\t{b * B}\t{b * B + B}\t{el}\t{0.35 + 0.003 * ((b * 13) % 100):.4f}\t"
                         f"{0.6 + 0.008 * ((b * 7) % 50):.4f}\n")
        pk = is_pk.copy()
        if c == "chr2":
            peaks.append(f"{c}\t{50 * B + 100}\t{53 * B + 4000}\n")   # 4 bins: bin 49 is then no peak bin
            peaks.append(f"{c}\t{49 * B + 100}\t{49 * B + 900}\n")
            pk[50:54] = True
        for b in np.flatnonzero(is_pk):
            peaks.append(f"{c}\t{b * B + 1000}\t{b * B + 2500}\n")
        s = 3 + (bins * 11 + ci * 5) % 40 + 50 * pk
        for b in bins:
            mid = b * B + 2500
            for k in range(s[b]):
                lines.append(f"s{ci}.{b}.{k}\t{c}\t{mid - 150 - k}\t{c}\t{mid + 150}\t+\t-\n")
            lines.append(f"x{ci}.{b}\t{c}\t{mid - 150}\t{c}\t{mid + 150}\t+\t+\n")      # same strand: not short
        i, j = np.triu_indices(nb, 1)
        keep = (j - i <= 210) & (pk[i] | pk[j])
        i, j = i[keep], j[keep]
        lam = 4.0 * ((s[i] + 1.0) * (s[j] + 1.0)) ** 0.35 / (1.0 + j - i) ** 0.9
        u = ((i * 7919 + j * 104729 + ci * 13) % 1000) / 1000
        cnt = np.floor(lam + u).astype(int)
        for a, b2 in PLANTED[c]:
            cnt[(i == a) & (j == b2)] += 60
        for a, b2, n in zip(i, j, cnt):
            for k in range(n):
                lines.append(f"l{ci}.{a}.{b2}.{k}\t{c}\t{a * B + 2001 + k}\t{c}\t{b2 * B + 3001}\t+\t-\n")
    lines.append("t1\tchr1\t100000\tchr2\t200000\t+\t-\n")
    (d / "features.txt").write_text("".join(feats))
    (d / "peaks.bed").write_text("".join(peaks))
    (d / "syn.pairs").write_text("## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
                                 + "".join(lines))
    return d / "syn.pairs", d / "peaks.bed", d / "features.txt"


# MAPS (ijuric/MAPS 21222de; Python 3.11, pandas 3.1.1; R 4.5.3, VGAM 1.1-14, data.table 1.18.6.1) on the
# feather-style files of synthetic() (callers/maps/derive_inputs.py), N_CHROMS=2, no sex chromosomes:
# bin pairs per regression and the sig3Dinteractions.bedpe it wrote.
MAPS_TESTED = {("chr1", "AND"): 558, ("chr1", "XOR"): 8681, ("chr2", "AND"): 649, ("chr2", "XOR"): 9277}
MAPS_BEDPE = """\
chr1	start1	end1	chr2	start2	end2	count	expected	fdr	ClusterLabel	ClusterSize	ClusterType	ClusterNegLog10P	ClusterSummit
chr1	525000	530000	chr1	825000	830000	62	1.78345483404144	0.000000000000000000000000000000000000000000000000000000000000000000000394948742199232	chr1_01	1	Singleton	72.1500934638795	1
chr1	1025000	1030000	chr1	1275000	1280000	62	1.8149413022015	0.000000000000000000000000000000000000000000000000000000000000000000000573026519148153	chr1_02	1	Singleton	71.687429482086	1
chr1	175000	180000	chr1	360000	365000	63	1.9615456628265	0.0000000000000000000000000000000000000000000000000000000000000000000614738182733505	chr1_03	1	Singleton	71.149879567282	1
chr1	75000	80000	chr1	225000	230000	63	2.68710774433829	0.000000000000000000000000000000000000000000000000000000000000332880694565696	chr1_1	3	BroadPeak	202.620340863266	0
chr1	75000	80000	chr1	230000	235000	62	1.95782467061351	0.000000000000000000000000000000000000000000000000000000000000000000596025705146517	chr1_1	3	BroadPeak	202.620340863266	0
chr1	80000	85000	chr1	225000	230000	62	1.92102579187494	0.000000000000000000000000000000000000000000000000000000000000000000282232318757111	chr1_1	3	BroadPeak	202.620340863266	1
chr2	125000	130000	chr2	425000	430000	62	1.77151546658523	0.000000000000000000000000000000000000000000000000000000000000000000000305159330342015	chr2_1	2	SharpPeak	150.848352540839	0
chr2	1525000	1530000	chr2	1725000	1730000	63	1.99579677830798	0.00000000000000000000000000000000000000000000000000000000000000000000669043711232968	chr2_2	2	SharpPeak	144.765625750457	0
chr2	125000	130000	chr2	430000	435000	62	1.40242513927304	0.00000000000000000000000000000000000000000000000000000000000000000000000000279751909373254	chr2_1	2	SharpPeak	150.848352540839	1
chr2	1525000	1530000	chr2	1735000	1740000	61	1.56050687950954	0.0000000000000000000000000000000000000000000000000000000000000000000000385933091849265	chr2_2	2	SharpPeak	144.765625750457	1
"""


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    d = tmp_path_factory.mktemp("maps")
    pairs, peaks, feats = synthetic(d)
    r = maps(pairs, peaks, features=feats, chroms=["chr1", "chr2"], verbose=False)
    return r, d


def test_maps_finds_the_planted_loops(syn):
    r, _ = syn
    got = {(c, s1 // B, s2 // B) for c, s1, s2 in zip(r.loops.a.chroms, r.loops.a.starts, r.loops.b.starts)}
    for c, pl in PLANTED.items():
        for a, b in pl:
            assert (c, a, b) in got
    t = r.tested
    xor = (t.a.chroms == "chr1") & (t.a.starts == 35 * B) & (t.b.starts == 72 * B)
    assert t.cols["set"][xor].tolist() == ["XOR"]
    # every call passes MAPS's cutoffs
    L = r.significant
    assert (L.cols["count"] >= 12).all() and (L.cols["ratio"] >= 2).all() and (L.cols["fdr"] < 0.01).all()


def test_maps_tests_the_bin_pairs_maps_tests(syn):
    r, _ = syn
    t = r.tested
    got = {}
    for c in ("chr1", "chr2"):
        for s in ("AND", "XOR"):
            got[(c, s)] = int(((t.a.chroms == c) & (t.cols["set"] == s)).sum())
    assert got == MAPS_TESTED
    assert (t.cols["dist"] >= 2).all() and (t.cols["dist"] <= 200).all()
    assert not ((t.a.starts < 2 * B) | (t.b.starts < 2 * B)).any()             # effective length 0
    chr2_49 = (t.a.chroms == "chr2") & ((t.a.starts == 49 * B) | (t.b.starts == 49 * B))
    assert not (t.cols["peak1"][chr2_49 & (t.a.starts == 49 * B)]).any()        # the first-peak quirk


def test_maps_writes_maps_bedpe(syn, tmp_path):
    r, _ = syn
    path = r.write(str(tmp_path), "syn.5k")["bedpe"]
    text = open(path).read()
    head = text.splitlines()[0].split("\t")
    assert head == ["chr1", "start1", "end1", "chr2", "start2", "end2", "count", "expected", "fdr", "ClusterLabel",
                    "ClusterSize", "ClusterType", "ClusterNegLog10P", "ClusterSummit"]
    assert len(text.splitlines()) == len(r.loops) + 1
    _same_bedpe(text, MAPS_BEDPE)


def test_maps_from_a_cooler_with_given_short_counts(syn):
    """Long counts from a cooler's pixels + short counts passed in = the pairs-file result (the synthetic short
    pairs fall in one bin, so the cooler's distance >= 1 pixels hold only long pairs)."""
    import cooler
    import pandas as pd
    r, d = syn
    p = pd.read_csv(d / "syn.pairs", sep="\t", comment="#", header=None, names=["id", "c1", "p1", "c2", "p2", "s1", "s2"])
    p = p[p.c1 == p.c2]
    off = {"chr1": 0, "chr2": 400}
    b1 = (p.p1 - (p.s1 == "+")) // B + p.c1.map(off)
    b2 = (p.p2 - (p.s2 == "+")) // B + p.c1.map(off)
    px = pd.DataFrame({"bin1_id": np.minimum(b1, b2), "bin2_id": np.maximum(b1, b2)}).value_counts().sort_index()
    px = px.rename("count").reset_index()
    bins = pd.DataFrame({"chrom": np.repeat(["chr1", "chr2"], 400), "start": np.tile(np.arange(400) * B, 2)})
    bins["end"] = bins.start + B
    uri = str(d / "syn.cool")
    cooler.create_cooler(uri, bins, px, ordered=True)
    _, short, _ = M._scan_pairs(str(d / "syn.pairs"), ["chr1", "chr2"], B, 1000, 200, "auto", None, None, 10 ** 6,
                                want_long=False, verbose=False)
    c = maps(uri, d / "peaks.bed", features=d / "features.txt", chroms=["chr1", "chr2"], short=short, verbose=False)
    for k in ("count", "expected", "P", "fdr"):
        np.testing.assert_array_equal(c.tested.cols[k], r.tested.cols[k])
    assert c.loops.cols["ClusterLabel"].tolist() == r.loops.cols["ClusterLabel"].tolist()
    with pytest.raises(ValueError, match="short"):
        maps(uri, d / "peaks.bed", features=d / "features.txt", chroms=["chr1", "chr2"], verbose=False)


def _same_bedpe(got, want):
    """Same rows and labels; numbers equal to 1e-12 (R's and scipy's incomplete gamma differ in the last digits)."""
    g = [ln.split("\t") for ln in got.strip().splitlines()]
    w = [ln.split("\t") for ln in want.strip().splitlines()]
    assert len(g) == len(w) and g[0] == w[0]
    for a, b in zip(g[1:], w[1:]):
        for k, (x, y) in enumerate(zip(a, b)):
            if k in (7, 8, 12):
                assert float(x) == pytest.approx(float(y), rel=1e-12)
            else:
                assert x == y
