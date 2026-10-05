"""bigWig signal: the three backends agree, windows at chromosome ends, tmm."""
import numpy as np
import pytest

import genomeblocks as gb
from genomeblocks import Loci, as_loci
from genomeblocks.backends.bigwig import _reduce, open_bigwig
from genomeblocks.signal import tmm

from conftest import CHROM_SIZES, installed_backends

BW_BACKENDS = installed_backends("bigwig")


def _brute(values, chrom, a, b, n_bins, stat):
    v = values.get(chrom)
    if v is None:
        return np.zeros(n_bins)
    seg = v[max(a, 0):min(b, len(v))]
    seg = np.concatenate([np.full(max(0, -a), np.nan), seg, np.full(max(0, b - len(v)), np.nan)])
    return _reduce(seg, n_bins, stat, 0.0)


@pytest.fixture
def values(bw_path):
    h = open_bigwig(bw_path)
    try:
        return {c: h.values(c, 0, n) for c, n in CHROM_SIZES.items()}
    finally:
        h.close()


@pytest.mark.parametrize("backend", BW_BACKENDS)
@pytest.mark.parametrize("stat", ["mean", "min", "max", "sum", "std", "coverage"])
def test_handles_bin_alike(bw_path, values, backend, stat):
    h = open_bigwig(bw_path, backend=backend)
    try:
        for a, b, nb in ((1000, 1300, 7), (500, 900, 3), (2000, 2100, 300), (0, 19_999, 13), (-50, 100, 10),
                         (19_900, 20_200, 6)):
            got = h.stats_array("chr1", a, b, n_bins=nb, stat=stat)
            assert np.allclose(got, _brute(values, "chr1", a, b, nb, stat), atol=1e-3), (a, b, nb)
        assert h.stats_array("chrZ", 0, 100, n_bins=3).tolist() == [0.0, 0.0, 0.0]
        assert np.isnan(h.values("chrZ", 0, 5)).all()
        assert np.isnan(h.values("chr1", -3, 2)[:3]).all()
    finally:
        h.close()


@pytest.mark.parametrize("backend", BW_BACKENDS)
def test_signal_cube_matches_brute_force_at_the_edges(bw_path, values, backend):
    L = as_loci([("chr1", 20, 80), ("chr1", 19_950, 20_000), ("chr2", 10, 20), ("chr1", 2000, 2100)])
    for nb, fl, span in ((20, 300, False), (200, 3000, False), (7, 100, False), (5, 0, True)):
        S = L.signal([bw_path], n_bins=nb, flank=fl, span=span, backend=backend, verbose=False, progress=False,
                     dtype=np.float64)
        assert S.shape == (4, 1, nb)
        for i in range(len(L)):
            if span:
                a, b = int(L.starts[i]), int(L.ends[i])
            else:
                c = (L.starts[i] + L.ends[i]) // 2
                a, b = int(c - fl), int(c + fl)
            assert np.allclose(S[i, 0], _brute(values, L.chroms[i], a, b, nb, "mean"), atol=1e-4), (i, nb, fl)


def test_signal_inputs_and_handles(bw_path, bw2_path, cre):
    import pyBigWig
    S1 = cre.signal({"a": bw_path, "b": bw2_path}, n_bins=10, flank=500, verbose=False, progress=False)
    S2 = cre.signal([bw_path, bw2_path], n_bins=10, flank=500, verbose=False, progress=False)
    assert S1.shape == (7, 2, 10) and np.array_equal(S1, S2)
    h = pyBigWig.open(bw_path)
    S3 = cre.signal([h], n_bins=10, flank=500, verbose=False, progress=False)
    assert np.allclose(S3[:, 0], S1[:, 0], atol=1e-4)
    with pytest.raises(ValueError, match="open pybigwig handle"):
        open_bigwig(h, backend="python")
    assert open_bigwig(h).backend == "pybigwig"
    S4 = cre.to_pandas().pipe(lambda df: gb.signal.signal(df, bw_path, n_bins=10, flank=500, verbose=False, progress=False))
    assert np.array_equal(S4[:, 0], S1[:, 0])
    with pytest.raises(ValueError, match="unknown stat"):
        cre.signal(bw_path, n_bins=4, flank=100, agg="median", verbose=False, progress=False)


def test_multiprocess_matches_sequential(bw_path, bw2_path):
    L = Loci.tile_genome(CHROM_SIZES, 500)
    a = L.signal([bw_path, bw2_path], n_bins=8, flank=400, verbose=False, progress=False)
    b = L.signal([bw_path, bw2_path], n_bins=8, flank=400, workers=2, verbose=False, progress=False)
    assert np.array_equal(a, b)


def test_tmm():
    rng = np.random.default_rng(0)
    base = np.abs(rng.normal(size=(50, 1, 4))) + 0.1
    cube = np.concatenate([base, 3.0 * base, 0.5 * base], axis=1)        # three scaled copies of one track
    out = tmm(cube)
    assert out.shape == cube.shape and np.isfinite(out).all()
    assert np.allclose(out[:, 0], out[:, 1]) and np.allclose(out[:, 0], out[:, 2])   # scaling undone
    lib = out.mean(axis=2).sum(axis=0)
    assert np.allclose(lib, lib[0])                                        # equal library sizes after TMM + CPM
    cube = np.abs(rng.normal(size=(50, 3, 4))) + 0.1
    cube[:, 2] = 0
    with pytest.warns(RuntimeWarning, match="no signal"):
        out = tmm(cube)
    assert np.isfinite(out).all() and (out[:, 2] == 0).all()
    with pytest.raises(ValueError, match="every track is empty"):
        tmm(np.zeros((5, 2, 3)))
