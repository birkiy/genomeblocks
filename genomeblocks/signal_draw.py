"""Visualization for genomic signal — heatmaps and average profiles.

Split out of ``signal.py`` so signal extraction stays matplotlib-free. Grouping
uses a plain ``groups: dict[str, Loci]`` (no Tags layer)::

    S = loci.signal(bigwigs)
    plot_heatmap(loci, S, groups={"up": up_loci, "down": down_loci})
"""
from __future__ import annotations
from typing import Dict, List, Sequence
from collections.abc import Sequence as SequenceABC

import numpy as np

from .loci import Loci
# matplotlib is imported inside the plotting functions: loci.py imports this
# module to attach Loci.plot_heatmap, so a top-level pyplot import would make
# every `from genomeblocks import Loci` pay for matplotlib.
# NOTE: signal/tmm are imported lazily inside compare_heatmap. loci.py imports
# this module at its tail, so a top-level `from .signal import ...` here would
# deadlock when `genomeblocks.signal` is the first thing imported.


def _bcast(x, n, name):
    """Convert scalar/str to list [x]*n; sequence of length n passes through"""
    if isinstance(x, str) or not isinstance(x, SequenceABC):
        return [x] * n
    if len(x) != n:
        raise ValueError(f"{name} must have length {n}")
    return list(x)


def _resolve_groups(S, gidx, sort):
    """Convert boolean masks to sorted index arrays."""
    groups = []
    for idx in gidx:
        indices = np.where(idx)[0]
        if sort == "group":
            order = np.argsort(S[indices].mean(axis=(1, 2)))[::-1]
            indices = indices[order]
        groups.append(indices)
    if sort == "global":
        vals = S.mean(axis=(1, 2))
        groups = [g[np.argsort(vals[g])[::-1]] for g in groups]
    return groups


def plot_heatmap(
    loci: Loci,
    S: np.ndarray,
    *,
    groups: Dict[str, "Loci"] | None = None,
    sets: List[str] | None = None,
    samples: List[str] | None = None,
    colors: Dict[str, tuple] | None = None,
    ymax: float = 10,
    ymin: float = 0,
    height: int = 3000,
    cmap: str = "Blues",
    vmax: float = 10,
    profile: bool = True,
    sort: str | None = "group",
    dpi: int = 100,
    ):
    """
    Plot heatmap of genomic signals.

    Args:
        loci: Genomic loci
        S: Signal array of shape (regions x tracks x bins)
        groups: Dict mapping group name -> Loci (if None, all loci are one group)
        sets: List of group names controlling row order (if None, uses all keys)
        samples: List of sample names (column labels)
        colors: Dictionary mapping set names to colors
        ymax: Maximum y-axis value for profile plots (scalar or per-track list)
        ymin: Minimum y-axis value for profile plots (scalar or per-track list)
        height: Region height in base pairs (for x-axis labels)
        cmap: Colormap for heatmaps (scalar or per-track list)
        vmax: Maximum value for heatmap color scaling (scalar or per-track list)
        profile: Include average profile above heatmaps
        sort: Sorting mode - "group" (per-group), "global", or None (no sort)
        dpi: Figure DPI

    Returns:
        matplotlib Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    if groups is None:
        sets = ["all"]
        gidx = [np.ones(len(loci), dtype=bool)]
    else:
        if sets is None:
            sets = list(groups.keys())
        gidx = []
        for name in sets:
            if name not in groups:
                raise ValueError(f"Group '{name}' not found in groups dict")
            g_uids = {l.uid for l in groups[name]}
            idx = np.array([loc.uid in g_uids for loc in loci])
            gidx.append(idx)

    if samples is None:
        samples = [f"track_{i}" for i in range(S.shape[1])]

    n = len(samples)
    cmaps = _bcast(cmap, n, "cmap")
    vms = _bcast(vmax, n, "vmax")
    ys = _bcast(ymax, n, "ymax")
    yl = _bcast(ymin, n, "ymin")

    if colors is None:
        colors = {k: plt.get_cmap('tab10')(i) for i, k in enumerate(sets)}

    g_ = _resolve_groups(S, gidx, sort)
    nb = S.shape[-1]

    rows = ([max(sum(len(g) for g in g_) // 4, 1)] if profile else []) + [max(len(g), 1) for g in g_]
    fig = plt.figure(figsize=(3 * n, 10), dpi=dpi)
    gs = gridspec.GridSpec(len(rows), n, height_ratios=rows)
    plt.subplots_adjust(hspace=0.05, wspace=0.3)

    for i, s in enumerate(samples):
        if profile:
            ax = fig.add_subplot(gs[0, i])
            for j, idx in enumerate(g_):
                if len(idx) > 0:
                    ax.plot(S[idx, i, :].mean(0), color=colors[sets[j]], lw=2,
                            label=sets[j])
            ax.set_ylim(yl[i], ys[i])
            ax.set_xticks([])
            if i == 0:
                ax.set_ylabel("CPM signal")
            if i == n - 1:
                ax.legend(fontsize=7, frameon=False)
            ax.set_title(s)

        for j, idx in enumerate(g_):
            ax = fig.add_subplot(gs[j + (1 if profile else 0), i])
            if len(idx) > 0:
                ax.imshow(S[idx, i, :], aspect='auto', cmap=cmaps[i],
                          vmin=0, vmax=vms[i])
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_ylabel(sets[j], rotation=0, ha='right', va='center')
            if j == len(g_) - 1:
                ax.set_xticks([0, nb // 2, nb])
                kb = round(height / 1000, 1)
                ax.set_xticklabels([f"-{kb}kb", "center", f"+{kb}kb"])
    return fig


def plot_profiles(
    loci: Loci,
    S: np.ndarray,
    *,
    groups: Dict[str, "Loci"] | None = None,
    sets: List[str] | None = None,
    colors: Dict[str, tuple] | None = None,
    ylim: float | None = None,
    dpi: int = 100,
    height: int = 3000,
    ):
    """
    Plot average signal profiles by set.

    Args:
        loci: Genomic loci
        S: Signal array of shape (regions × tracks × bins)
        groups: Dict mapping group name -> Loci (if None, all loci are one group)
        sets: Group names / order (if None, uses all keys in groups)
        colors: Dictionary mapping set names to colors
        ylim: Y-axis limit
        dpi: Figure DPI
        height: Region height in base pairs

    Returns:
        matplotlib Figure
    """
    import matplotlib.pyplot as plt

    # Build uid to index mapping for efficient lookup
    uid_to_idx = {loc.uid: i for i, loc in enumerate(loci)}

    if groups is None:
        sets = ["all"]
        gidx = [np.ones(len(loci), dtype=bool)]
    else:
        if sets is None:
            sets = list(groups.keys())
        gidx = []
        for name in sets:
            if name not in groups:
                raise ValueError(f"Group '{name}' not found in groups dict")
            g_uids = {l.uid for l in groups[name]}
            idx = np.array([loc.uid in g_uids for loc in loci])
            gidx.append(idx)

    if colors is None:
        colors = {k: plt.get_cmap('tab10')(i) for i, k in enumerate(sets)}
    fig, axs = plt.subplots(1, len(sets), figsize=(len(sets) * 3, 3), dpi=dpi, squeeze=False)
    axs = axs[0]

    for ax, idx, lab in zip(axs, gidx, sets):
        m = S[idx, 0, :].mean(0)
        ax.plot(m, color=colors[lab], lw=2)
        ax.set_title(lab)
        ax.set_xticks([0, S.shape[-1] // 2, S.shape[-1]])
        kb = height // 1000
        ax.set_xticklabels([f"-{kb}kb", "center", f"+{kb}kb"])
        ax.set_ylim(0, (ylim if ylim is not None else np.percentile(m, 99)))
        ax.set_yticks([])

    axs[0].set_ylabel("signal")
    return fig


def compare_heatmap(
    a: Loci,
    b: Loci,
    bigwigs: Sequence[str],
    *,
    a_name: str = "A",
    b_name: str = "B",
    common_name: str = "common",
    sets: List[str] | None = None,
    samples: List[str] | Dict[str, List[int]] | None = None,
    n_bins: int = 200,
    flank: int = 3_000,
    agg: str = "mean",
    normalize: bool = True,
    cmap: str = "Blues",
    vmax: float = 10,
    ymax: float = 10,
    ymin: float = 0,
    profile: bool = True,
    sort: str | None = "group",
    colors: Dict[str, tuple] | None = None,
    dpi: int = 100,
    S: np.ndarray | None = None,
    signal_kw: dict | None = None,
    ):
    """
    Compare two Loci sets as a heatmap grid.

    Computes a-specific (a - b), b-specific (b - a), and common (a & b),
    extracts signal from bigwigs, and plots a heatmap grid.

    Args:
        a, b: Two Loci objects to compare
        bigwigs: BigWig file paths (columns of the heatmap)
        a_name, b_name: Labels for the specific rows
        common_name: Label for the common row
        sets: Row order. Defaults to [a_name, common_name, b_name].
              Pass a subset or reorder to customise.
        samples: Column labels. Either a list of strings (one per bigwig,
                 no merging) or a dict mapping column names to lists of
                 bigwig indices to average together, e.g.
                 ``{"ATAC": [0, 1], "H3K4me3": [2, 3]}``.
        n_bins, flank, agg: Passed to signal() when S is not provided
        normalize: Apply TMM normalization before merging (default True)
        cmap, vmax, ymax, ymin: Passed to plot_heatmap (scalar or per-track).
                 When samples is a dict, lengths must match the merged count.
        profile: Include average profile above heatmaps
        sort: Sorting mode - "group" (per-group), "global", or None
        colors: Dict mapping group names to colors
        dpi: Figure DPI
        S: Pre-computed signal array for the union loci. Skips extraction
           and normalization.
        signal_kw: Extra kwargs passed to signal()

    Returns:
        (fig, union_loci, S, groups)
    """
    a_specific = a - b
    b_specific = b - a
    common = a & b
    union = a_specific | common | b_specific

    groups = {
        a_name: a_specific,
        common_name: common,
        b_name: b_specific,
    }

    # --- signal extraction + normalization + merging ----
    from .signal import signal, tmm  # lazy: avoids an import cycle (see top)
    if S is None:
        kw = dict(n_bins=n_bins, flank=flank, agg=agg)
        if signal_kw:
            kw.update(signal_kw)
        S = signal(union, bigwigs, **kw)
        S = np.nan_to_num(S)

        if normalize:
            S = tmm(S)
            S = np.nan_to_num(S)

    # merge tracks when samples is a dict
    if isinstance(samples, dict):
        merged = np.stack(
            [S[:, idx, :].mean(axis=1) for idx in samples.values()],
            axis=1,
        )
        S = merged
        sample_labels = list(samples.keys())
    else:
        sample_labels = samples  # list[str] or None

    if sets is None:
        sets = [a_name, common_name, b_name]

    height = flank

    fig = plot_heatmap(
        union, S,
        groups=groups,
        sets=sets,
        samples=sample_labels,
        colors=colors,
        ymax=ymax, ymin=ymin,
        height=height,
        cmap=cmap, vmax=vmax,
        profile=profile,
        sort=sort,
        dpi=dpi,
    )

    return fig, union, S, groups



Loci.plot_heatmap = plot_heatmap
Loci.plot_profiles = plot_profiles
