from __future__ import annotations
import time
import shutil
from multiprocessing import Process, Value, cpu_count
from multiprocessing.shared_memory import SharedMemory
from typing import List, Tuple, Sequence, Dict
from collections.abc import Sequence as SequenceABC

import numpy as np
import pyBigWig as bw
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import gridspec

from .loci import Loci
from .tags import Tags


def plan_workers(n_tracks: int, n_loci: int, *, cores: int | None = None,
                max_bw_parallel: int = 6) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """
    Plan the distribution of work across multiple processes.
    
    Args:
        n_tracks: Number of bigwig tracks
        n_loci: Number of genomic loci
        cores: Number of CPU cores to use (default: all available)
        max_bw_parallel: Maximum number of bigwig files to process in parallel
    
    Returns:
        List of work chunks as [((track_start,track_end), (loci_start,loci_end)), ...]
    """
    c = cores or cpu_count()
    t_chunks = min(n_tracks, max_bw_parallel, c)
    t_step = (n_tracks + t_chunks - 1) // t_chunks
    l_per_t = max(1, c // t_chunks)
    l_step = (n_loci + l_per_t - 1) // l_per_t

    plan = []
    for ti in range(t_chunks):
        t_lo = ti * t_step
        t_hi = min(n_tracks, (ti + 1) * t_step)
        for li in range(l_per_t):
            l_lo = li * l_step
            l_hi = min(n_loci, (li + 1) * l_step)
            plan.append(((t_lo, t_hi), (l_lo, l_hi)))
    return plan


def _worker(
    t_b: Tuple[int, int],
    l_b: Tuple[int, int],
    bwp: Sequence[str],
    loci: Loci,
    n_bins: int,
    flank: int,
    agg: str,
    shm_name: str,
    shape: Tuple[int, int, int],
    dt: np.dtype,
    ctr: Value, # type: ignore
    span: bool = False,
    ):
    """Worker process for parallel signal extraction from bigwig files."""
    t_lo, t_hi = t_b
    l_lo, l_hi = l_b
    shm = SharedMemory(name=shm_name)
    cube = np.ndarray(shape, dtype=dt, buffer=shm.buf)

    hs = [bw.open(p) for p in bwp[t_lo:t_hi]]
    ch = hs[0].chroms()

    try:
        for r in range(l_lo, l_hi):
            loc = loci[r]
            chrom = loc.chrom
            if chrom not in ch:  # skip unknown chrom
                with ctr.get_lock():
                    ctr.value += (t_hi - t_lo)
                continue

            if span:
                L = max(0, loc.start)
                R = min(ch[chrom], loc.end)
                pre, post = 0, 0
                core = max(1, n_bins)
            else:
                c = (loc.start + loc.end) // 2
                L, R = c - flank, c + flank
                pre = max(0, -L)
                post = max(0, R - ch.get(chrom, 0))
                L = max(0, L)
                R = min(ch.get(chrom, 0), R)
                core = n_bins - (pre + post)

            if core > 0:
                for off, h in enumerate(hs):
                    t = t_lo + off
                    xs = h.stats(chrom, L, R, nBins=core, type=agg)
                    xs = np.nan_to_num(xs).astype(dt, copy=False)
                    if span:
                        arr = xs
                    else:
                        arr = np.concatenate((
                            np.zeros(pre, dtype=dt),
                            xs,
                            np.zeros(post, dtype=dt)
                        ))
                    cube[r, t, :len(arr)] = arr

            with ctr.get_lock():
                ctr.value += (t_hi - t_lo)
    finally:
        for h in hs:
            h.close()
        shm.close()


def signal(
    loci: Loci,
    bigwigs: Sequence[str],
    *,
    n_bins: int = 200,
    flank: int = 3_000,
    agg: str = "mean",
    dtype: str | np.dtype = np.float32,
    progress: bool = True,
    max_bw_parallel: int = 6,
    max_workers: int | None = None,
    span: bool = False,
    verbose: bool = True
    ) -> np.ndarray:
    """
    Extract signal from bigwig files for given genomic loci.

    Args:
        loci: Genomic loci to extract signal from
        bigwigs: List of bigwig file paths
        n_bins: Number of bins to divide each region into
        flank: Number of base pairs to include on each side of locus center
        agg: Aggregation method ('mean', 'max', 'min', etc.)
        dtype: Numpy dtype for the output array
        progress: Show progress bar
        max_bw_parallel: Maximum number of bigwig files to process in parallel
        max_workers: Maximum number of worker processes
        span: Use full region span instead of center±flank

    Returns:
        numpy array of shape (n_loci, n_tracks, n_bins)
    """
    if verbose:
        print(f"[INFO] Extracting {len(bigwigs)} bigwigs for {len(loci)} loci "
            f"into {n_bins} bins (span={span}, agg='{agg}'). 🪏")


    n_loci, n_tracks = len(loci), len(bigwigs)
    if n_loci == 0:
        raise ValueError("No loci provided.")

    bytes_need = n_loci * n_tracks * n_bins * np.dtype(dtype).itemsize
    if bytes_need > 0.5 * shutil.disk_usage("/").free:
        raise MemoryError("Cube may exceed safe RAM; try disk-chunk mode.")

    shm = SharedMemory(create=True, size=bytes_need)
    cube = np.ndarray((n_loci, n_tracks, n_bins), dtype=dtype, buffer=shm.buf)
    cube.fill(0)

    plan = plan_workers(
        n_tracks,
        n_loci,
        cores=(max_workers or cpu_count()),
        max_bw_parallel=max_bw_parallel,
    )

    ctr = Value('i', 0)
    procs = []
    for t_b, l_b in plan:
        p = Process(target=_worker, args=(
            t_b, l_b, bigwigs, loci, n_bins, flank, agg,
            shm.name, cube.shape, dtype, ctr, span
        ))
        p.start()
        procs.append(p)

    if progress:
        tot = n_loci * n_tracks
        with tqdm(total=tot, dynamic_ncols=True) as bar:
            last = 0
            while any(p.is_alive() for p in procs):
                v = ctr.value
                bar.update(v - last)
                last = v
                time.sleep(0.2)
            bar.update(ctr.value - last)

    for p in procs:
        p.join()
        if p.exitcode != 0:
            shm.close()
            shm.unlink()
            raise RuntimeError(f"worker {p.pid} exit {p.exitcode}")

    out = cube.copy()
    shm.close()
    shm.unlink()
    return out


def tmm(cube: np.ndarray) -> np.ndarray:
    """TMM-normalize a signal cube (regions × tracks × bins).

    Computes per-track TMM normalization factors and library-size
    scaling so that tracks become comparable.

    Args:
        cube: Signal array of shape (regions, tracks, bins)

    Returns:
        Normalized copy of the cube (same shape).
    """
    import conorm
    means = np.nanmean(cube, axis=2)                        # (regions, tracks)
    factors = conorm.tmm_norm_factors(means)                 # (tracks,)
    lib_size = means.sum(0)                                  # (tracks,)
    scale = 1.0 / (factors * lib_size / 1_000_000)
    return cube * scale[None, :, None]


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
    tags: Tags | None = None,
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
        tags: Tags object for grouping loci (if None, all loci treated as one group)
        groups: Dict mapping group names to Loci objects (alternative to tags).
                Mutually exclusive with tags.
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
    if groups is not None and tags is not None:
        raise ValueError("Cannot specify both 'tags' and 'groups'")

    if groups is not None:
        tags = Tags.make(loci, verbose=False)
        tags.add(groups)
        if sets is None:
            sets = list(groups.keys())

    if tags is None:
        sets = ["all"]
        gidx = [np.ones(len(loci), dtype=bool)]
    else:
        if sets is None:
            sets = sorted(tags.keys())
        gidx = []
        for tag_name in sets:
            if tag_name not in tags:
                raise ValueError(f"Tag '{tag_name}' not found in Tags object")
            tag_view = tags[tag_name]
            tag_uids = tag_view.uids if hasattr(tag_view, 'uids') else set(tag_view)
            idx = np.array([loc.uid in tag_uids for loc in loci])
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
    tags: Tags | None = None,
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
        tags: Tags object for grouping loci (if None, all loci treated as one group)
        sets: List of tag names to use as groups (if None, uses all tags in Tags object)
        colors: Dictionary mapping set names to colors
        ylim: Y-axis limit
        dpi: Figure DPI
        height: Region height in base pairs
    
    Returns:
        matplotlib Figure
    """
    # Build uid to index mapping for efficient lookup
    uid_to_idx = {loc.uid: i for i, loc in enumerate(loci)}
    
    if tags is None:
        # No tags provided - treat all loci as one group
        sets = ["all"]
        gidx = [np.ones(len(loci), dtype=bool)]
    else:
        # Use tags to group loci
        if sets is None:
            sets = sorted(tags.keys())
        
        gidx = []
        for tag_name in sets:
            if tag_name not in tags:
                raise ValueError(f"Tag '{tag_name}' not found in Tags object")
            tag_view = tags[tag_name]
            # Get UIDs from the tag and create boolean index
            tag_uids = tag_view.uids if hasattr(tag_view, 'uids') else set(tag_view)
            idx = np.array([loc.uid in tag_uids for loc in loci])
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
        (fig, union_loci, S, tags)
    """
    a_specific = a - b
    b_specific = b - a
    common = a & b
    union = a_specific | common | b_specific

    tags = Tags.make(union, verbose=False)
    tags.add({
        a_name: a_specific,
        common_name: common,
        b_name: b_specific,
    })

    # --- signal extraction + normalization + merging ----
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
        tags=tags,
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

    return fig, union, S, tags


Loci.signal = signal

Loci.plot_heatmap = plot_heatmap
Loci.plot_profiles = plot_profiles