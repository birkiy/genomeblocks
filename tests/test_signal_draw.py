"""Heatmaps and profiles (headless)."""
import matplotlib.pyplot as plt
import numpy as np
import pytest

import genomeblocks as gb
from genomeblocks import signal_draw as sd


@pytest.fixture
def cube(cre):
    return np.random.default_rng(0).random((len(cre), 2, 12))


def test_plot_heatmap_groups_and_samples(cre, cube):
    fig = sd.plot_heatmap(cre, cube, groups={"up": np.arange(3), "down": cre.tail(4)}, samples=["a", "b"], vmax=[1, 1])
    assert fig.axes
    fig = sd.plot_heatmap(cre.to_pandas(), cube, groups={"mask": np.array([True] * 3 + [False] * 4)}, cmap="Reds")
    assert fig.axes
    with pytest.raises(ValueError, match="samples must have length 2"):
        sd.plot_heatmap(cre, cube, samples=["only one"])
    with pytest.raises(ValueError, match="vmax must have one entry per track, got 3 for 2 tracks"):
        sd.plot_heatmap(cre, cube, vmax=[1, 2, 3])
    with pytest.raises(ValueError, match="rows"):
        sd.plot_heatmap(cre, cube[:3])


def test_plot_profiles(cre, cube):
    fig = sd.plot_profiles(cre, cube, groups={"a": [0, 1], "b": [2, 3, 4]}, track=1)
    assert fig.axes


def _xticklabels(fig):
    labels = [t.get_text() for t in fig.axes[-1].get_xticklabels()]
    plt.close(fig)
    return labels


def test_profile_and_heatmap_tick_labels_agree(cre, cube):
    assert _xticklabels(sd.plot_profiles(cre, cube, height=500)) == ["-0.5kb", "center", "+0.5kb"]
    for h in (500, 3000, 1500, 250):
        assert _xticklabels(sd.plot_profiles(cre, cube, height=h)) == _xticklabels(sd.plot_heatmap(cre, cube, height=h))


def test_compare_heatmap(cre, bw_path, bw2_path):
    a, b = cre.head(4), cre.tail(5)
    fig, union, S, groups = sd.compare_heatmap(a, b, [bw_path, bw2_path], n_bins=10, flank=500, normalize=False)
    assert fig.axes and len(union) == len(a | b) - len(a & b) and S.shape[0] == len(union)
    fig, *_ = gb.compare_heatmap(a.to_pandas(), b, [bw_path], n_bins=10, flank=500, normalize=True, samples=["x"])
    assert fig.axes
    S = np.random.default_rng(0).random((len(union), 4, 10))
    merged = {"x": [0, 1], "y": [2, 3]}
    fig, _, S2, _ = sd.compare_heatmap(a, b, [], S=S, samples=merged, cmap=["Reds", "Blues"], flank=500)
    assert fig.axes and S2.shape[1] == 2
    for kw in ({"cmap": ["Reds"] * 4}, {"vmax": [1, 2, 3, 4]}, {"ymax": [1, 2, 3, 4]}):
        with pytest.raises(ValueError, match="one entry per plotted column .merged sample., got 4 for 2 columns"):
            sd.compare_heatmap(a, b, [], S=S, samples=merged, flank=500, **kw)


def test_plot_motif_heatmap(cre):
    M = np.random.default_rng(0).poisson(0.3, size=(len(cre), 2, 20)).astype(float)
    fig = gb.plot_motif_heatmap(cre, M, ["CTCF", "AR"], r=500, groups={"a": np.arange(3), "b": np.arange(3, 7)})
    assert fig.axes
    with pytest.raises(ValueError, match="rows, motifs, bins"):
        gb.plot_motif_heatmap(cre, M[:, 0])


def test_group_mask_inputs(cre):
    assert np.flatnonzero(sd.group_mask(cre, [0, 2])).tolist() == [0, 2]
    assert np.flatnonzero(sd.group_mask(cre, cre.head(2))).tolist() == [0, 1]
    assert np.flatnonzero(sd.group_mask(cre, cre.head(2).to_pandas())).tolist() == [0, 1]
    assert np.flatnonzero(sd.group_mask(cre, np.array([True] + [False] * 6))).tolist() == [0]
    with pytest.raises(ValueError, match="group mask"):
        sd.group_mask(cre, np.array([True, False]))
