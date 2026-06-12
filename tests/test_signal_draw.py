import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.figure import Figure

from genomeblocks import Loci
from genomeblocks.locus import Locus
from genomeblocks.signal_draw import plot_heatmap, plot_profiles


def _data():
    up = Loci([Locus("chr1", i * 1000, i * 1000 + 200) for i in range(4)])
    down = Loci([Locus("chr1", 50000 + i * 1000, 50000 + i * 1000 + 200) for i in range(4)])
    regions = up + down
    S = np.random.rand(len(regions), 2, 20).astype(np.float32)
    return regions, S, {"up": up, "down": down}


def test_plot_heatmap_returns_figure():
    regions, S, groups = _data()
    fig = plot_heatmap(regions, S, groups=groups, sets=["up", "down"],
                       samples=["t0", "t1"], profile=False)
    assert isinstance(fig, Figure)


def test_plot_heatmap_no_groups_is_one_block():
    regions, S, _ = _data()
    fig = plot_heatmap(regions, S, profile=False)
    assert isinstance(fig, Figure)


def test_plot_profiles_returns_figure():
    regions, S, groups = _data()
    fig = plot_profiles(regions, S, groups=groups)
    assert isinstance(fig, Figure)


def test_plot_heatmap_attached_to_loci():
    regions, S, groups = _data()
    fig = regions.plot_heatmap(S, groups=groups, profile=False)
    assert isinstance(fig, Figure)


def test_unknown_group_name_raises():
    regions, S, groups = _data()
    import pytest
    with pytest.raises(ValueError):
        plot_heatmap(regions, S, groups=groups, sets=["nope"], profile=False)
