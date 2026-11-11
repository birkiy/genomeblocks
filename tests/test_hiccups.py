"""HiCCUPS (Juicer's CPU path): the kernel, thresholds, Poisson, Java numbers and HashSet order against values
juicer_tools 2.13.06 itself produced, the float32 shortcut's invariant, and planted loops end to end.

Reference values come from small Java harnesses run against juicer_tools_2.13.06.jar (kept with the parity run in
fastloops/callers/hiccups/javatest/): Kernel.java calls GPUController.runCPUVersion (by reflection) on the matrix
built by ``_kernel_block``; Fdr.java calls HiCCUPSUtils.calculateThresholdAndFDR on the histogram of
``_fdr_hist``; Pois.java prints PoissonDistributionImpl.probability; HashOrder.java prints
``new ArrayList<>(new HashSet<>(list))`` of Feature2D objects; FloatStr.java prints String.valueOf(float).
"""
import numpy as np
import pandas as pd
import pytest

import genomeblocks._hiccups as H

F32 = np.float32


# ── the background kernel ────────────────────────────────────────────────────

def _kernel_block(row0, col0, size=512):
    """Kernel.java's matrix: counts by distance plus two loops, a divisive norm with NaN bins, Juicer-normalised
    (float32(n / (v_i v_j)) above the diagonal, NaN records dropped)."""
    I = np.arange(row0, row0 + size)[:, None]
    J = np.arange(col0, col0 + size)[None, :]
    n = (np.floor(60.0 / (1.0 + (J - I) / 2.5)).astype(np.int64) + np.where((I * 31 + J * 17) % 7 == 0, 2, 0)
         + np.where((np.abs(I - 150) <= 1) & (np.abs(J - 190) <= 1), 40, 0)
         + np.where((np.abs(I - 300) <= 1) & (np.abs(J - 760) <= 1), 9, 0))
    with np.errstate(invalid="ignore"):
        B = (n / (_kr(I).astype(np.float64) * _kr(J).astype(np.float64))).astype(F32).astype(np.float64)
    B[np.isnan(B) | (I > J)] = 0
    st = H._Block.__new__(H._Strip)
    st.B, st.R0, st.C0, st.R1, st.C1 = B, row0, col0, row0 + size, col0 + size
    st.tables()
    return st


def _kr(I):
    return np.where(I % 97 == 5, np.nan, 0.5 + ((I * 13) % 10) / 10.0).astype(F32)


D32 = (60.0 / (1.0 + np.arange(2048) / 2.5) + 0.25).astype(F32)
V32 = _kr(np.arange(2048))
# juicer_tools 2.13.06 runCPUVersion, window 7, peak width 4: row, col, expectedBL, Donut, H, V, binBL, Donut, H, V,
# observed. Tiles: rows/cols from 0 (diagonal tile) and rows from 0 x columns from 472 (off-diagonal tile).
JUICER_KERNEL = {(0, 0): """
100 107 1.8031758 3.3488698 1.3258535 5.4681063 2 5 1 7 15
100 108 4.9798465 5.7711782 6.3939953 8.043875 7 7 8 9 14
150 190 1.0427247 0.9875892 1.2761728 1.2935507 0 0 1 1 43
150 191 1.831797 1.6348648 2.043948 2.0889795 2 2 3 3 43
149 189 5.8272114 5.6628327 6.915672 6.7856016 7 7 8 8 43
200 230 1.3704958 1.297868 1.2034461 1.5530779 1 1 0 1 4
250 400 0.29990864 0.16884483 0.29773328 0.3795584 0 0 0 0 0
300 480 0.09139839 0.082931526 0.017343448 0.017343448 0 0 0 0 0
60 70 3.3384378 3.2464883 4.1531982 4.051689 5 5 6 6 12""", (0, 472): """
300 760 0.09440004 0.083038375 0.01735909 0.017359084 0 0 0 0 9
301 760 0.15103398 0.128022 0.03029948 0.0 0 0 0 0 9
299 761 0.44367877 0.33989298 0.068368204 0.2285451 0 0 0 0 9
100 600 0.06328346 0.07112661 0.14958684 0.14958681 0 0 0 0 2
480 500 2.019533 1.9426945 2.5125928 2.5846443 3 2 4 4 8
450 900 0.066550024 0.07110596 0.15969926 0.15969928 0 0 0 0 0"""}


@pytest.mark.parametrize("tile", list(JUICER_KERNEL))
def test_kernel_is_juicers_to_the_bit(tile):
    ref = np.array([[float(x) for x in l.split()] for l in JUICER_KERNEL[tile].strip().splitlines()])
    R, C = ref[:, 0].astype(np.int64), ref[:, 1].astype(np.int64)
    st = _kernel_block(*tile)
    ex = H._exact_at(st, R, C, D32, V32, 4, 7)
    np.testing.assert_array_equal(ex["e"], ref[:, 2:6].T.astype(F32))       # float32 values, bit for bit
    np.testing.assert_array_equal(ex["bins"], ref[:, 6:10].T)
    o, bins, flag = H._fast(st, R, C, D32, V32, H._dtables(D32, 4, int((C - R).max())), 4, 7)
    np.testing.assert_array_equal(o, ref[:, 10])
    np.testing.assert_array_equal(bins[:, ~flag], ref[:, 6:10].T[:, ~flag])


def test_kernel_masks_nan_norm_rows():
    st = _kernel_block(0, 0)                    # row 102 has a NaN norm: every expected NaN, chunk -1
    ex = H._exact_at(st, np.array([102]), np.array([140]), D32, V32, 4, 7)
    assert np.isnan(ex["e"]).all() and (ex["bins"] == -1).all()


@pytest.mark.parametrize("window, pw", [(7, 4), (5, 2), (3, 1)])
def test_fast_sums_give_the_exact_chunks_unless_flagged(window, pw):
    """The float64 summed-area shortcut agrees with the float32 kernel on every pixel it does not flag, and flags
    few (the invariant that makes the shortcut exact)."""
    st = _kernel_block(0, 0)
    I, J = np.meshgrid(np.arange(20, 492), np.arange(20, 492), indexing="ij")
    keep = (J - I > pw + 2) & ~np.isnan(V32[I]) & ~np.isnan(V32[J])
    R, C = I[keep], J[keep]
    o, bins, flag = H._fast(st, R, C, D32, V32, H._dtables(D32, pw, int((C - R).max())), pw, window)
    ex = H._exact_at(st, R, C, D32, V32, pw, window)
    np.testing.assert_array_equal(bins[:, ~flag], ex["bins"][:, ~flag])
    assert flag.mean() < 1e-3


# ── thresholds, Poisson, numbers ─────────────────────────────────────────────

def _fdr_hist():
    i = np.arange(H.W1)[:, None]
    j = np.arange(H.W2)[None, :]
    h = ((1000000 * (i + 1)) >> np.minimum(j // (i + 2), 62)) + np.where((j % 97 == 0) & (j < 300), 3, 0)
    h[7] = 0
    return np.repeat(h[None].astype(np.int64), 4, axis=0)


# Fdr.java, 10% FDR (calculateThresholdAndFDR gets 1 / 0.1): the threshold per chunk, and a few q-values.
# Chunk 7 has no pixel (stays 10000); chunk 39 hits Juicer's "j == 0" case (Poisson(2^(40/3)) lies mostly above
# 10000 counts) and gets 10000 - 2.
JUICER_THRESHOLDS = [4, 4, 5, 5, 6, 8, 9, 10000, 13, 16, 20, 24, 30, 36, 45, 56, 69, 86, 107, 133, 166, 208, 261, 326,
                     409, 513, 645, 1035, 1246, 1507, 1828, 2225, 2717, 3327, 4085, 5027, 6202, 7666, 9494, 9998]
JUICER_Q = {(0, 0): 1.0, (0, 3): 0.3566931486129761, (0, 10): 2.843869515345432e-05, (5, 8): 0.11013404279947281,
            (8, 20): 0.0010117586934939027, (20, 150): 3.3680598735809326, (39, 9999): 0.0, (7, 0): 0.0}


def test_fdr_thresholds_and_q_values_are_juicers():
    thr, q = H._fdr_tables(_fdr_hist(), 1.0 / 0.1)
    np.testing.assert_array_equal(thr[0], np.array(JUICER_THRESHOLDS, F32))
    for (i, j), v in JUICER_Q.items():
        assert q[0, i, j] == F32(v)


def test_poisson_pmf_is_commons_math_2_2():
    # PoissonDistributionImpl(2^((i+1)/3)).probability(k) as juicer_tools computes it (Pois.java)
    ref = {(0, 0): 0.2836764218990301, (0, 1): 0.35740989530944706, (0, 5): 0.007505137509018255,
           (14, 30): 0.06814215338823745, (20, 100): 0.0014497717284445608, (39, 8000): 4.780896350562382e-126,
           (39, 9999): 2.4702034994823162e-05}
    for (i, k), v in ref.items():
        got = H._poisson_pmf(i, 10000)[k]
        assert got == pytest.approx(v, rel=1e-11)
        assert F32(got) == F32(v)                       # what Juicer uses


def test_java_float_strings():
    vals = [1.8532424292061478e-05, 5.125570723407691e-08, 3012.1572265625, 0.0013010199181735516, 393.0,
            0.0010000000474974513, 0.0009999900357797742, 1e7, 9999999.0, 123456.703125, 1.0, 0.10000000149011612,
            np.nan, np.inf, 0.0]
    java = ["1.8532424E-5", "5.1255707E-8", "3012.1572", "0.0013010199", "393.0", "0.001", "9.9999E-4", "1.0E7",
            "9999999.0", "123456.7", "1.0", "0.1", "NaN", "Infinity", "0.0"]
    assert H._jfloats(np.array(vals, F32)) == java


def test_hashset_order_is_javas():
    """``new ArrayList<>(new HashSet<>(features))`` of Feature2D objects (HashOrder.java, chr7, 13 features with
    duplicates)."""
    rows = []
    for i in range(13):
        res = 10000 if i % 3 == 0 else 5000
        s1 = 1000000 + ((i * 7919) % 200) * 5000
        rows += [(s1, s1 + 50000 + ((i * 104729) % 300) * 5000, res)] * (2 if i % 5 == 0 else 1)
    a = np.array(rows)
    out = H._hashset({"chrom": np.zeros(len(a), np.int64), "start1": a[:, 0], "start2": a[:, 1], "res": a[:, 2]},
                     ["chr7"])
    java = ("1165000,2230000,5000 1000000,1050000,10000 1545000,1690000,5000 1380000,2010000,5000 "
            "1760000,2970000,5000 1595000,1790000,5000 1975000,2750000,5000 1140000,1430000,10000 "
            "1355000,2710000,10000 1570000,2490000,10000 1190000,1530000,5000 1950000,3450000,5000 "
            "1785000,2270000,10000")
    assert " ".join(f"{x},{y},{r}" for x, y, r in zip(out["start1"], out["start2"], out["res"])) == java


# ── post-processing rules ────────────────────────────────────────────────────

def _table(rows, res):
    a = np.array(rows, dtype=np.float64).reshape(-1, 3)
    n = len(a)
    return {"chrom": np.zeros(n, np.int64), "start1": a[:, 0].astype(np.int64), "start2": a[:, 1].astype(np.int64),
            "res": np.full(n, res), "observed": a[:, 2].astype(F32), "e": np.ones((4, n), F32),
            "bins": np.zeros((4, n), np.int64), "fdr": np.zeros((4, n), F32), "peak": a[:, 2].astype(F32)}


def test_coalesce_grows_clusters_from_the_strongest_pixel():
    t = _table([(100000, 300000, 20), (105000, 305000, 30), (110000, 310000, 25), (160000, 300000, 40)], 5000)
    out = H._coalesce(t, 5000, 20000)
    # seeds by observed: (160k, 300k) alone (the others are > 20 kb away); then (105k, 305k) takes (110k, 310k)
    # (7.07 kb), centroid (107500, 307500) -> 107500 (int), radius round(3535.5) = 3536, then (100k, 300k)
    # (10.6 kb <= 20 kb + 3536): centroid (105000, 305000), radius round(7071.07) = 7071
    assert out["start1"].tolist() == [160000, 105000]
    assert out["numCollapsed"].tolist() == [1, 3]
    assert out["centroid1"].tolist() == [162500, 107500] and out["centroid2"].tolist() == [302500, 307500]
    assert out["radius"].tolist() == [0, 7071]


def test_merge_keeps_juicers_resolution_rules():
    five = _table([(1000000, 1500000, 50), (3000000, 3500000, 50), (5000000, 5080000, 20), (7000000, 7600000, 120)],
                  5000)
    ten = _table([(1010000, 1510000, 60), (9000000, 9500000, 60)], 10000)
    tf = _table([(1000000, 1525000, 80), (12000000, 12500000, 80)], 25000)
    m = H._merge_resolutions({5000: five, 10000: ten, 25000: tf}, ["chr1"])
    got = set(zip(m["start1"].tolist(), m["res"].tolist()))
    # 5 kb: kept with a 10 kb loop within 20 kb, within 110 kb of the diagonal, or observed > 100; 10 kb: kept
    # without a 5 kb loop within 20 kb; 25 kb: kept without a merged loop within 50 kb
    assert got == {(1000000, 5000), (5000000, 5000), (7000000, 5000), (9000000, 10000), (12000000, 25000)}


def test_expected_is_normalised_per_chromosome():
    """Juicer's chromosome factor: the expected of a chromosome sums to its normalised contacts."""
    rng = np.random.default_rng(3)
    data = {"names": ["a", "b"], "sizes": {"a": 3_000_000, "b": 2_000_000}, "pix": {}, "v": {}}
    for c, W in (("a", 600), ("b", 400)):
        i = rng.integers(0, W, 200000)
        j = np.minimum(i + rng.geometric(0.02, len(i)), W - 1)
        df = pd.DataFrame({"i": i, "j": j}).value_counts().reset_index()
        data["pix"][c] = (df["i"].to_numpy(), df["j"].to_numpy(), df["count"].to_numpy(np.float64))
        data["v"][c] = np.ones(W)
    e = H._juicer_expected(data, 5000)
    for c, W in (("a", 600), ("b", 400)):
        b1, b2, n = data["pix"][c]
        k = np.arange(W)
        assert np.sum((W - k) * e[c].astype(np.float64)) == pytest.approx(n.sum(), rel=1e-5)


# ── end to end ───────────────────────────────────────────────────────────────

LOOPS = [(2_500_000, 2_800_000), (4_500_000, 5_000_000), (7_500_000, 7_650_000), (10_000_000, 11_500_000)]


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """One 14 Mb chromosome: distance decay, per-bin bias, four 3x3-bin loops; a 5 kb cooler coarsened into an
    .mcool with 10 and 25 kb, balanced."""
    import cooler
    d = tmp_path_factory.mktemp("hiccups")
    rng = np.random.default_rng(11)
    L, res = 14_000_000, 5000
    W = L // res
    bias = rng.lognormal(0, 0.2, W)
    i = rng.integers(0, W, 4_000_000)
    j = i + np.exp(rng.uniform(0, np.log(3000), len(i))).astype(np.int64)
    keep = (j < W) & (rng.random(len(i)) < bias[i] * bias[np.minimum(j, W - 1)] / bias.max() ** 2)
    i, j = i[keep], j[keep]
    for a, b in LOOPS:
        k = 400
        i = np.r_[i, a // res + rng.integers(-1, 2, k)]
        j = np.r_[j, b // res + rng.integers(-1, 2, k)]
    px = (pd.DataFrame({"bin1_id": i, "bin2_id": j}).value_counts().rename("count").reset_index()
          .sort_values(["bin1_id", "bin2_id"]))
    bins = pd.DataFrame({"chrom": "chr1", "start": np.arange(W) * res, "end": np.arange(1, W + 1) * res})
    base = str(d / "t.5000.cool")
    cooler.create_cooler(base, bins, px)
    out = str(d / "t.mcool")
    cooler.zoomify_cooler(base, out, [5000, 10000, 25000], chunksize=10_000_000)
    for r in (5000, 10000, 25000):
        cooler.balance_cooler(cooler.Cooler(f"{out}::resolutions/{r}"), store=True)
    return out


def test_hiccups_finds_the_planted_loops(planted, tmp_path):
    r = H.hiccups(planted, verbose=False)
    m = r.loops
    assert set(r.pixels) == {5000, 10000, 25000}
    for a, b in LOOPS:
        near = (np.abs(m.a.starts - a) <= 20000) & (np.abs(m.b.starts - b) <= 20000)
        assert near.any(), (a, b)
    assert len(m) <= 2 * len(LOOPS)
    hit = m.cols["count"] > 1.75 * m.cols["expectedDonut"]           # the post-processing filters held
    assert hit.all() and (m.cols["count"] > 1.5 * m.cols["expectedH"]).all()
    assert set(np.unique(m.cols["resolution"])) <= {5000, 10000, 25000}

    paths = r.write(tmp_path)
    lines = open(paths["merged"]).read().splitlines()
    assert lines[0].split("\t")[11:] == H._FINAL and lines[1] == "# juicer_tools version 2.13.06"
    assert len(lines) == len(m) + 2
    first = lines[2].split("\t")
    assert first[10] == "0,255,255" and first[11] == f"{float(m.cols['count'][0])}"
    thr = open(paths["thresholds_5000"]).read().splitlines()
    assert len(thr) == H.W1 and thr[0].startswith("0\t")


def test_hiccups_restrict_limits_the_search(planted):
    full = H.hiccups(planted, resolutions=(25000,), restrict=False, verbose=False)
    cpu = H.hiccups(planted, resolutions=(25000,), verbose=False)
    # at 25 kb the CPU path searches only the diagonal tiles (472 bins = 11.8 Mb); the planted loops are closer
    assert len(cpu.enriched[25000]) <= len(full.enriched[25000])
    span = (cpu.enriched[25000].b.starts - cpu.enriched[25000].a.starts) // 25000
    assert span.max() < 472


def test_hiccups_workers_do_not_change_the_result(planted):
    one = H.hiccups(planted, resolutions=(10000,), verbose=False)
    four = H.hiccups(planted, resolutions=(10000,), workers=4, verbose=False)
    a, b = one.enriched[10000], four.enriched[10000]
    np.testing.assert_array_equal(a.a.starts, b.a.starts)
    np.testing.assert_array_equal(a.cols["expectedDonut"], b.cols["expectedDonut"])


def test_hiccups_rejects_resolutions_without_defaults(planted):
    with pytest.raises(ValueError, match="peak_width"):
        H.hiccups(planted, resolutions=(5000, 7000))
