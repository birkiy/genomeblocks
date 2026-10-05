"""Heatmaps and average profiles of a signal cube.

Kept apart from ``signal.py`` so extraction stays matplotlib-free. Groups of
rows are a plain ``{name: rows}`` dict, where ``rows`` is a Loci (or anything
:func:`~genomeblocks.as_loci` reads — matched to the plotted loci by
coordinates), a boolean mask, or row numbers::

    S = cre.signal(bigwigs)
    plot_heatmap(cre, S, groups={"up": up, "down": cre["lfc"] < 0})
"""
from __future__ import annotations

from collections.abc import Sequence as SequenceABC
from typing import Dict, List, Sequence

import numpy as np


def _bcast(x, n, name, unit="track", units="tracks"):
    """Scalar / str -> [x] * n; a sequence of length n passes through."""
    if isinstance(x, str) or not isinstance(x, SequenceABC):
        return [x] * n
    if len(x) != n:
        raise ValueError(f"{name} must have one entry per {unit}, got {len(x)} for {n} {units}")
    return list(x)


def _flank_ticklabels(height) -> List[str]:
    """x tick labels at the left edge, the centre and the right edge, in kb."""
    kb = round(height / 1000, 1)
    return [f"-{kb}kb", "center", f"+{kb}kb"]


def group_mask(loci, rows) -> np.ndarray:
    """Boolean mask over ``loci`` for one group: a mask, row numbers, or intervals
    (matched by chromosome, start and end)."""
    n = len(loci)
    a = None if hasattr(rows, "codes") else (np.asarray(rows) if isinstance(rows, (np.ndarray, list)) else None)
    if a is not None and a.dtype == bool:
        if len(a) != n:
            raise ValueError(f"a group mask has {len(a)} values for {n} rows")
        return a
    if a is not None and a.dtype.kind in "iu":
        m = np.zeros(n, bool)
        m[a] = True
        return m
    from .interop import as_loci
    G = loci._check(as_loci(rows))
    import pandas as pd
    have = pd.MultiIndex.from_arrays([G.codes.astype(np.int64), G.starts, G.ends])
    return pd.MultiIndex.from_arrays([loci.codes.astype(np.int64), loci.starts, loci.ends]).isin(have)


def _resolve_groups(S, gidx, sort):
    """Boolean masks -> index arrays, sorted by mean signal within (or across) groups."""
    groups = []
    for idx in gidx:
        indices = np.where(idx)[0]
        if sort == "group":
            indices = indices[np.argsort(S[indices].mean(axis=(1, 2)))[::-1]]
        groups.append(indices)
    if sort == "global":
        vals = S.mean(axis=(1, 2))
        groups = [g[np.argsort(vals[g])[::-1]] for g in groups]
    return groups


def _group_masks(loci, groups, sets):
    if groups is None:
        return ["all"], [np.ones(len(loci), dtype=bool)]
    sets = list(groups) if sets is None else list(sets)
    for name in sets:
        if name not in groups:
            raise ValueError(f"Group {name!r} not found in groups (have: {', '.join(map(str, groups))})")
    return sets, [group_mask(loci, groups[name]) for name in sets]


def plot_heatmap(loci, S: np.ndarray, *, groups: Dict[str, object] | None = None, sets: List[str] | None = None,
                 samples: List[str] | None = None, colors: Dict[str, tuple] | None = None, ymax=10, ymin=0,
                 height: int = 3000, cmap="Blues", vmax=10, profile: bool = True, sort: str | None = "group",
                 dpi: int = 100, ylabel: str = "CPM signal"):
    """deeptools-style heatmap of a (rows x tracks x bins) cube, one column per track.

    Args:
        loci: the loci the cube was extracted for (anything ``as_loci`` takes).
        S: the cube.
        groups: ``{name: rows}`` blocks of rows (Loci / mask / row numbers);
            None = one block.
        sets: block order (default: the keys of ``groups``).
        samples: column titles.
        ymax, ymin, cmap, vmax: scalar or one per track.
        height: bp from the centre to the edge (x tick labels).
        profile: draw the mean profile of each block above the heatmaps.
        sort: 'group' (by mean signal within each block), 'global', or None.
    Returns:
        matplotlib Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from .interop import as_loci
    loci = as_loci(loci)
    S = np.asarray(S)
    if S.ndim != 3 or S.shape[0] != len(loci):
        raise ValueError(f"S must be (rows, tracks, bins) with {len(loci)} rows, got shape {S.shape}")
    sets, gidx = _group_masks(loci, groups, sets)
    n = S.shape[1]
    if samples is None:
        samples = [f"track_{i}" for i in range(n)]
    samples = list(samples)
    if len(samples) != n:
        raise ValueError(f"samples must have length {n} (one label per track), got {len(samples)}")
    cmaps, vms = _bcast(cmap, n, "cmap"), _bcast(vmax, n, "vmax")
    ys, yl = _bcast(ymax, n, "ymax"), _bcast(ymin, n, "ymin")
    if colors is None:
        colors = {k: plt.get_cmap("tab10")(i) for i, k in enumerate(sets)}
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
                    ax.plot(S[idx, i, :].mean(0), color=colors[sets[j]], lw=2, label=sets[j])
            ax.set_ylim(yl[i], ys[i])
            ax.set_xticks([])
            if i == 0:
                ax.set_ylabel(ylabel)
            if i == n - 1:
                ax.legend(fontsize=7, frameon=False)
            ax.set_title(s)
        for j, idx in enumerate(g_):
            ax = fig.add_subplot(gs[j + (1 if profile else 0), i])
            if len(idx) > 0:
                ax.imshow(S[idx, i, :], aspect="auto", cmap=cmaps[i], vmin=0, vmax=vms[i])
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_ylabel(sets[j], rotation=0, ha="right", va="center")
            if j == len(g_) - 1:
                ax.set_xticks([0, nb // 2, nb])
                ax.set_xticklabels(_flank_ticklabels(height))
    return fig


def plot_profiles(loci, S: np.ndarray, *, groups: Dict[str, object] | None = None, sets: List[str] | None = None,
                  colors: Dict[str, tuple] | None = None, ylim: float | None = None, dpi: int = 100,
                  height: int = 3000, track: int = 0):
    """Mean profile of one track (``track``) per block of rows, side by side."""
    import matplotlib.pyplot as plt
    from .interop import as_loci
    loci = as_loci(loci)
    S = np.asarray(S)
    sets, gidx = _group_masks(loci, groups, sets)
    if colors is None:
        colors = {k: plt.get_cmap("tab10")(i) for i, k in enumerate(sets)}
    fig, axs = plt.subplots(1, len(sets), figsize=(len(sets) * 3, 3), dpi=dpi, squeeze=False)
    for ax, idx, lab in zip(axs[0], gidx, sets):
        m = S[idx, track, :].mean(0) if idx.any() else np.zeros(S.shape[-1])
        ax.plot(m, color=colors[lab], lw=2)
        ax.set_title(lab)
        ax.set_xticks([0, S.shape[-1] // 2, S.shape[-1]])
        ax.set_xticklabels(_flank_ticklabels(height))
        top = ylim if ylim is not None else (np.percentile(m, 99) if m.any() else 1.0)
        ax.set_ylim(0, top)
        ax.set_yticks([])
    axs[0][0].set_ylabel("signal")
    return fig


def compare_heatmap(a, b, bigwigs: Sequence[str], *, a_name: str = "A", b_name: str = "B",
                    common_name: str = "common", sets: List[str] | None = None, samples=None, n_bins: int = 200,
                    flank: int = 3_000, agg: str = "mean", normalize: bool = True, cmap="Blues", vmax=10,
                    ymax=10, ymin=0, profile: bool = True, sort: str | None = "group", colors=None, dpi: int = 100,
                    S: np.ndarray | None = None, signal_kw: dict | None = None):
    """Two interval sets as one heatmap grid: a-only, shared, b-only.

    ``a`` and ``b`` are anything ``as_loci`` takes. The cube is extracted for
    their union (a-only, then common, then b-only rows), TMM-normalised when
    ``normalize``, and ``samples`` may be a ``{column: [track indices]}`` dict
    to average replicate tracks into one column.

    Returns:
        (fig, union_loci, S, groups) — ``groups`` maps each block name to its rows.
    """
    from .interop import as_loci
    from .signal import signal, tmm
    a = as_loci(a)
    b = a._check(as_loci(b))
    a_only, b_only, common = a - b, b - a, a & b
    union = a_only + common + b_only
    n1, n2 = len(a_only), len(a_only) + len(common)
    rows = np.arange(len(union))
    groups = {a_name: rows < n1, common_name: (rows >= n1) & (rows < n2), b_name: rows >= n2}
    if S is None:
        kw = dict(n_bins=n_bins, flank=flank, agg=agg)
        kw.update(signal_kw or {})
        S = np.nan_to_num(signal(union, bigwigs, **kw))
        if normalize:
            S = np.nan_to_num(tmm(S))
    if isinstance(samples, dict):
        S = np.stack([S[:, idx, :].mean(axis=1) for idx in samples.values()], axis=1)
        labels = list(samples)
        per = dict(unit="plotted column (merged sample)", units="columns")
        cmap, vmax = _bcast(cmap, len(labels), "cmap", **per), _bcast(vmax, len(labels), "vmax", **per)
        ymax, ymin = _bcast(ymax, len(labels), "ymax", **per), _bcast(ymin, len(labels), "ymin", **per)
    else:
        labels = samples
    fig = plot_heatmap(union, S, groups=groups, sets=sets or [a_name, common_name, b_name], samples=labels,
                       colors=colors, ymax=ymax, ymin=ymin, height=flank, cmap=cmap, vmax=vmax, profile=profile,
                       sort=sort, dpi=dpi)
    return fig, union, S, groups
