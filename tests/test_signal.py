import numpy as np
import pytest

from genomeblocks import tmm
from genomeblocks.signal import _even_ranges, _tmm_norm_factors
from genomeblocks.signal_draw import _bcast, _resolve_groups   # plotting helpers


def test_bcast_scalar_and_sequence():
    assert _bcast(5, 3, "x") == [5, 5, 5]
    assert _bcast("Blues", 3, "cmap") == ["Blues", "Blues", "Blues"]   # str = scalar
    assert _bcast([1, 2, 3], 3, "x") == [1, 2, 3]


def test_bcast_length_mismatch_raises():
    with pytest.raises(ValueError):
        _bcast([1, 2], 3, "x")


def test_even_ranges_partition():
    rngs = _even_ranges(10, 3)
    assert rngs[0][0] == 0 and rngs[-1][1] == 10
    # contiguous and covering
    for (a, b), (c, d) in zip(rngs, rngs[1:]):
        assert b == c
    assert sum(b - a for a, b in rngs) == 10


def test_resolve_groups_from_masks():
    S = np.random.rand(6, 1, 4).astype(np.float32)
    gidx = [np.array([True, True, True, False, False, False]),
            np.array([False, False, False, True, True, True])]
    groups = _resolve_groups(S, gidx, sort=None)
    assert sorted(groups[0].tolist()) == [0, 1, 2]
    assert sorted(groups[1].tolist()) == [3, 4, 5]


def test_tmm_preserves_shape():
    cube = np.abs(np.random.rand(20, 3, 10)).astype(np.float32) + 0.1
    out = tmm(cube)
    assert out.shape == cube.shape
    assert np.isfinite(out).all()


def test_tmm_norm_factors_geomean_one():
    # vendored edgeR TMM (replaces the conorm dependency): one factor per
    # column (sample), scaled to a geometric mean of 1.
    rng = np.random.default_rng(1)
    data = np.abs(rng.normal(20, 6, (300, 4))) + 1.0
    f = _tmm_norm_factors(data)
    assert f.shape == (4,)
    assert np.isclose(np.exp(np.log(f).mean()), 1.0)
