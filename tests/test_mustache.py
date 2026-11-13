"""Mustache's loop caller: its block layout and distance rules, the band scale space against scipy on the full
block, the local normalisation and BH against Mustache / statsmodels, and planted loops on a synthetic cooler
compared with Mustache 1.3.3's own calls."""
import math

import numpy as np
import pytest

from genomeblocks import _mustache as M


@pytest.mark.parametrize("given, res, want", [(None, 5000, 2_000_000), (None, 10000, 2_000_000),
                                              (None, 25000, 5_000_000), (None, 1000, 2_000_000),
                                              (None, 500, 1_000_000), (100_000, 5000, 1_000_000),
                                              (3_000_000, 5000, 3_000_000), (60_000_000, 5000, 50_000_000),
                                              (20_000_000, 5000, 10_000_000)])
def test_distance_rules_follow_mustache_main(given, res, want):
    assert M._distance_filter(given, res) == want


def test_blocks_overlap_by_the_distance_and_end_at_the_last_contact():
    assert M._blocks(1500, 400) == ([0], [1500])
    assert M._blocks(5000, 400) == ([0, 1600, 3000], [2000, 3600, 5000])
    assert M._blocks(9000, 1500) == ([0, 1500, 3000, 4500, 6000], [3000, 4500, 6000, 7500, 9000])


def test_sigmas_are_mustache_scales():
    sig = M._sigmas(1.6)
    assert len(sig) == 12
    assert [s for s, _ in sig][:3] == [1.6, 1.6 * 2 ** 0.1, 1.6 * 2 ** 0.2]
    # radius = ceil(2 sigma), as Mustache sets the kernel width 2 ceil(2 sigma) + 1
    assert [int(t * s + 0.5) for s, t in sig] == [math.ceil(2 * s) for s, _ in sig]


def test_band_scale_space_equals_scipy_on_the_full_block():
    from scipy.ndimage import gaussian_filter, maximum_filter
    rng = np.random.default_rng(3)
    C, D = 300, 60
    x = rng.integers(0, C, 20000)
    y = x + rng.integers(0, D + 1, len(x))
    k = np.unique((x * C + y)[y < C])
    x, y = k // C, k % C
    cm = np.zeros((C, C))
    cm[x, y] = rng.standard_normal(len(x))
    t = y - x >= 4
    r, c = x[t], y[t]
    dd = np.arange(C)[None, :] - np.arange(C)[:, None]
    cm[(dd <= 4) | (dd > D)] = 2
    sigs = [M._sigmas(o) for o in (1.6, 3.2)]
    R = max(int(tt * s + 0.5) for sig in sigs for s, tt in sig)
    got = list(M._band_dog(cm, D, R, r, c, sigs))
    assert len(got) == 22
    for o, sig in enumerate(sigs):
        G = [gaussian_filter(cm, s, truncate=tt) for s, tt in sig]
        for j, (L, Mx) in enumerate(got[11 * o:11 * (o + 1)]):
            Lf = G[j] - G[j + 1]
            assert np.array_equal(L, Lf[r, c])                          # bit for bit
            assert np.array_equal(Mx, maximum_filter(Lf, footprint=np.ones((3, 3)), mode="constant")[r, c])


def _diag_pixels():
    """A 3000-bin chromosome at 5 kb with pixels on diagonals 4..400 (deterministic, no RNG)."""
    x, y = np.meshgrid(np.arange(3000), np.arange(4, 401))
    x, y = x.ravel(), (x + y).ravel()
    keep = (y < 3000) & ((x * 7 + y * 13) % 5 != 0)
    x, y = x[keep], y[keep]
    o = np.lexsort((y, x))
    x, y = x[o], y[o]
    v = 1.0 / (1 + (y - x)) * (1 + 0.5 * np.sin(x / 37.0) + 0.25 * np.cos(y / 11.0))
    return x.astype(np.int64), y.astype(np.int64), v


# mustache.normalize_sparse(x, y, v, 5000, 400) (Mustache 1.3.3 at df71aa1, numpy 1.26.4) on _diag_pixels():
# the values at these rows of the row-major pixel list. 3000 bins: the local (2 Mb window) branch.
NORM_ROWS = [0, 1000, 123456, 400000, 777777, 887005]
NORM_REF = [0.5341371443418756, 0.6452453409836872, -1.4141872330641119, 1.2493015841742687,
            -0.5260468052265151, -1.4477639974174017]
# ... and on the pixels of _diag_pixels() with x, y < 700: the whole-diagonal branch ((700 - 400) bins <= 2 Mb)
NORM700_ROWS = [0, 1000, 50000, 123456, 157905]
NORM700_REF = [0.5873971432452267, 0.6531991679893306, -1.914100509694317, -1.001626095683078,
               0.3912586146381252]


def test_local_normalisation_matches_mustache():
    x, y, v = _diag_pixels()
    z = M._normalize(x, y, v, 5000, 400)
    np.testing.assert_allclose(z[NORM_ROWS], NORM_REF, rtol=1e-12)


def test_whole_diagonal_normalisation_matches_mustache():
    x, y, v = _diag_pixels()
    k = y < 700
    x, y, v = x[k], y[k], v[k]
    z = M._normalize(x, y, v, 5000, 400)
    np.testing.assert_allclose(z[NORM700_ROWS], NORM700_REF, rtol=1e-12)
    z = M._normalize(x, y, v, 5000, 399)                   # Mustache leaves the last diagonal (D) as it is
    last = y - x == 399
    assert last.any() and np.array_equal(z[last], v[last]) and not np.array_equal(z[~last], v[~last])


def test_bh_matches_statsmodels():
    p = np.array([0.01, 0.04, 0.03, 0.04, 0.2, 1e-5, 0.5, 0.03, 1.0])
    # statsmodels 0.15.0 multipletests(p, method="fdr_bh")[1]
    want = [0.045000000000000005, 0.060000000000000005, 0.060000000000000005, 0.060000000000000005,
            0.2571428571428572, 9e-05, 0.5625, 0.060000000000000005, 1.0]
    assert M._bh(p).tolist() == want                                     # the same arithmetic: bit for bit


def synthetic_cool(path):
    """chr1 (10 Mb, 5 kb bins): Poisson distance decay, five planted 5x5 blobs; chr2 (0.8 Mb) ignored by
    Mustache's > 1 Mb rule. Unit weights. Returns the planted bin pairs."""
    import cooler
    import pandas as pd
    rng = np.random.default_rng(7)
    res, n = 5000, 2000
    planted = [(300, 360), (700, 790), (1100, 1250), (1400, 1460), (1700, 1900)]
    i, j = np.meshgrid(np.arange(n), np.arange(0, 401), indexing="ij")
    i, j = i.ravel(), (i + j).ravel()
    keep = j < n
    i, j = i[keep], j[keep]
    lam = 40.0 * (j - i + 1.0) ** -0.75
    for a, b in planted:
        lam += 6 * lam * np.exp(-((i - a) ** 2 + (j - b) ** 2) / (2 * 1.2 ** 2))
    cnt = rng.poisson(lam)
    m = cnt > 0
    bins = pd.DataFrame({"chrom": ["chr1"] * n + ["chr2"] * 160,
                         "start": np.r_[np.arange(n), np.arange(160)] * res,
                         "end": np.r_[np.arange(1, n + 1), np.arange(1, 161)] * res})
    bins["weight"] = 1.0
    pix = pd.DataFrame({"bin1_id": i[m], "bin2_id": j[m], "count": cnt[m]})
    i2 = np.arange(n, n + 150)
    pix = pd.concat([pix, pd.DataFrame({"bin1_id": i2, "bin2_id": i2 + 5, "count": 3})], ignore_index=True)
    cooler.create_cooler(str(path), bins, pix, ordered=True)
    return planted


# Mustache 1.3.3 (df71aa1; numpy 1.26.4, scipy 1.16.3): ``mustache -f syn.cool -r 5000 -o syn.tsv -p 1`` on
# synthetic_cool(), the whole of syn.tsv. Its FDR column depends on the CPU: the local normalisation's
# np.convolve sums run through BLAS ddot, and OpenBLAS picks a ddot kernel per CPU, so the integer-count ties
# round differently. One output per kernel family (OPENBLAS_CORETYPE=SkylakeX on AVX-512 CPUs; Haswell, Zen
# and Sandybridge on AVX2 / AVX CPUs; Prescott, the generic one); the port writes Mustache's file under each.
_HEAD = "BIN1_CHR\tBIN1_START\tBIN1_END\tBIN2_CHROMOSOME\tBIN2_START\tBIN2_END\tFDR\tDETECTION_SCALE\n"
_LOOPS = [("1500000\t1505000", "1800000\t1805000", 0), ("3500000\t3505000", "3950000\t3955000", 1),
          ("7000000\t7005000", "7300000\t7305000", 0)]
_FDR = {"SkylakeX": ("3.897234377525649e-05", "0.0027665667848083184"),
        "Haswell": ("3.896939018144252e-05", "0.0027663060180849186"),
        "Prescott": ("3.898689013148715e-05", "0.0027673450454184976")}
MUSTACHE_SYN = {k: _HEAD + "".join(f"chr1\t{a}\tchr1\t{b}\t{q[i]}\t2.111212657236631\n" for a, b, i in _LOOPS)
                for k, q in _FDR.items()}


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    d = tmp_path_factory.mktemp("mustache")
    planted = synthetic_cool(d / "syn.cool")
    return str(d / "syn.cool"), planted


def test_mustache_finds_the_planted_loops(syn):
    cool, planted = syn
    r = M.mustache(cool, threads=2, verbose=False)
    got = set(zip((r.loops.a.starts // 5000).tolist(), (r.loops.b.starts // 5000).tolist()))
    near = {(a + da, b + db) for a, b in planted for da in (-1, 0, 1) for db in (-1, 0, 1)}
    assert got <= near                                        # nothing but the planted blobs
    assert {p for p in planted if p[1] - p[0] <= 100} <= got  # the blobs near the diagonal, at their centre
    assert (r.loops.cols["FDR"] < 0.2).all()
    assert set(r.loops.a.chroms.astype(str)) == {"chr1"}     # chr2 is under 1 Mb


def test_mustache_writes_mustaches_output(syn, tmp_path):
    cool, _ = syn
    r = M.mustache(cool, threads=1, verbose=False)
    path = r.write(str(tmp_path), "syn")["loops"]
    assert open(path).read() in MUSTACHE_SYN.values()          # Mustache's file, for this CPU's BLAS kernel


def test_weights_given_as_an_array_equal_a_weight_column(syn, tmp_path):
    import shutil
    import h5py
    cool, _ = syn
    copy = tmp_path / "w.cool"
    shutil.copy(cool, copy)
    with h5py.File(copy, "r+") as h5:
        nb = h5["bins/start"].shape[0]
        w = 1 + 0.5 * np.sin(np.arange(nb) / 7.0)
        h5["bins"].create_dataset("w2", data=w)
    a = M.mustache(str(copy), weight="w2", threads=1, verbose=False)
    b = M.mustache(str(copy), weight=w, threads=1, verbose=False)
    assert len(a.loops) and a.loops.to_pandas().equals(b.loops.to_pandas())
    with pytest.raises(ValueError, match="one value per bin"):
        M.mustache(cool, weight=np.ones(3), verbose=False)
    with pytest.raises(ValueError, match="no bins/KR"):
        M.mustache(cool, weight="KR", verbose=False)
