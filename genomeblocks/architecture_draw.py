"""Region drawing for :class:`~genomeblocks.Architecture` graphs.

Plain matplotlib on top of the tables: the CREs of a region become nodes
(optionally merged when closer than ``merge_distance``), the edges among
them become lines, and the positions come from the graph backend's layout
(:func:`genomeblocks.backends.graph.layout_edges`) — graph-tool when
installed, else scipy's spectral layout, or igraph / networkx on request::

    A.draw("chr8:127.5-128.5 Mb")
    A.draw(("chr8", 127_500_000, 128_500_000), layout="genomic", merge_distance=5_000,
           vertex_size_by="strength", vertex_color="annot", backend="networkx")

Kept separate from ``architecture.py`` so the graph object stays free of
matplotlib; everything here is imported lazily.
"""
from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np


def _scale(values, lo_hi: Tuple[float, float]) -> np.ndarray:
    v = np.asarray(values, float)
    lo, hi = lo_hi
    if not len(v):
        return v
    if np.nanmax(v) > np.nanmin(v):
        return lo + (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v)) * (hi - lo)
    return np.full(len(v), (lo + hi) / 2.0)


def _arc(x0, x1, n=40):
    """A half circle from x0 to x1 above the line y = 0."""
    c, r = (x0 + x1) / 2.0, abs(x1 - x0) / 2.0
    th = np.linspace(0, np.pi, n)
    return np.column_stack([c + r * np.cos(th), r * np.sin(th)])


def draw(A, region, *, layout: str = "spring", backend: Optional[str] = None, merge_distance: Optional[int] = None,
         vertex_size_by=None, edge_width_by: str = "w", vertex_size_range: Tuple[float, float] = (40, 400),
         edge_width_range: Tuple[float, float] = (0.6, 4.0), vertex_color=None, cmap: str = "viridis",
         edge_color: str = "#b8b8b8", figsize: Tuple[float, float] = (8, 6), show_labels: bool = True,
         label_prop: str = "gene", font_size: int = 8, linked_only: bool = True, seed: int = 0,
         title: bool = True, ax=None):
    """Draw the CREs of ``region`` and the edges among them.

    Args:
        A: the Architecture.
        region: ``'chr1:1,000-2,000'``, ``'chr8:127.5-128.5 Mb'``, ``(chrom, start, end)`` or a Locus.
        layout: ``'spring'`` (the graph backend's force-directed or spectral
            layout), ``'circular'``, or ``'genomic'`` (nodes on a line at their
            genomic position, edges as arcs — a browser-like view).
        backend: graph backend for the spring layout (``None`` = default).
        merge_distance: merge CREs whose centres are within this many bp into
            one node (edge weights and size values are summed).
        vertex_size_by: a vertex column name or an array (one value per Loci
            row) that scales node sizes; default: number of CREs in the node.
        edge_width_by: edge column that scales line widths (default ``'w'``).
        vertex_color: a colour, a vertex column name (numbers are mapped
            through ``cmap``, strings get one colour per value) or an array
            aligned to the Loci rows.
        show_labels / label_prop: label nodes with that vertex column (unique
            non-empty values joined by ``/``), else with their coordinates.
        linked_only: skip CREs without edges (default).
        ax: draw into an existing Axes.

    Returns:
        the matplotlib Axes.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from .backends.graph import layout_edges
    from .locus import parse_region

    chrom, start, end = parse_region(region)
    L = A.loci
    rows = L.overlap_rows(chrom, start, end)
    if linked_only:
        rows = rows[A.degree[rows] > 0]
    if not len(rows):
        raise ValueError(f"no CREs{' with links' if linked_only else ''} in {chrom}:{start:,}-{end:,}")
    rows = rows[np.argsort(L.starts[rows], kind="stable")]
    centers = L.centers[rows].astype(float)

    # nodes: one per CRE, or one per cluster of nearby CREs
    if merge_distance is not None and len(rows) > 1:
        cluster = np.concatenate([[0], np.cumsum(np.diff(centers) > merge_distance)])
    else:
        cluster = np.arange(len(rows))
    m = int(cluster.max()) + 1
    count = np.bincount(cluster, minlength=m).astype(float)
    x = np.bincount(cluster, weights=centers, minlength=m) / count
    node_start = np.full(m, np.iinfo(np.int64).max, np.int64)
    node_end = np.zeros(m, np.int64)
    np.minimum.at(node_start, cluster, L.starts[rows])
    np.maximum.at(node_end, cluster, L.ends[rows])

    # edges among the nodes (parallel edges summed, self edges dropped)
    pos_of = np.full(len(L), -1, np.int64)
    pos_of[rows] = cluster
    keep = (pos_of[A.src] >= 0) & (pos_of[A.tgt] >= 0)
    s, t = pos_of[A.src[keep]], pos_of[A.tgt[keep]]
    if edge_width_by is not None and edge_width_by in A.ep:
        w = np.asarray(A.ep[edge_width_by], float)[keep]
    else:
        w = np.ones(int(keep.sum()))
    inter = s != t
    s, t, w = s[inter], t[inter], np.nan_to_num(w[inter])
    key = np.minimum(s, t).astype(np.int64) * m + np.maximum(s, t)
    uk, inv = np.unique(key, return_inverse=True)
    es, et = uk // m, uk % m
    ew = np.bincount(inv.ravel(), weights=w, minlength=len(uk)) if len(uk) else np.zeros(0)

    # sizes
    if vertex_size_by is None:
        size_val = count
    else:
        v = np.asarray(A.vp[vertex_size_by] if isinstance(vertex_size_by, str) else vertex_size_by, float)
        size_val = np.bincount(cluster, weights=np.nan_to_num(v[rows]), minlength=m)
    sizes = _scale(size_val, vertex_size_range)
    widths = _scale(ew, edge_width_range)

    # colours
    colors = "#4A90E2"
    legend = None
    if vertex_color is not None:
        vals = None
        if isinstance(vertex_color, str) and vertex_color in A.vp:
            vals = np.asarray(A.vp[vertex_color])[rows]
        elif not isinstance(vertex_color, str):
            vals = np.asarray(vertex_color)[rows]
        if vals is None:
            colors = vertex_color
        elif vals.dtype.kind in "biuf":
            node_val = np.bincount(cluster, weights=np.nan_to_num(vals.astype(float)), minlength=m) / count
            colors = plt.get_cmap(cmap)(_scale(node_val, (0.0, 1.0)))
        else:
            labels_c = np.array([""] * m, dtype=object)
            for k, c in enumerate(cluster.tolist()):
                if not labels_c[c] and str(vals[k]):
                    labels_c[c] = str(vals[k])
            uniq = sorted(set(labels_c.tolist()))
            palette = plt.get_cmap("tab10")
            cmap_s = {u: palette(i % 10) for i, u in enumerate(uniq)}
            colors = [cmap_s[u] for u in labels_c.tolist()]
            legend = cmap_s

    # labels
    labels: Optional[Sequence[str]] = None
    if show_labels:
        if label_prop in A.vp:
            raw = np.asarray(A.vp[label_prop], dtype=object)[rows]
            labels = []
            for c in range(m):
                parts = []
                for v in raw[cluster == c].tolist():
                    sv = str(v) if v is not None else ""
                    if sv and sv != "nan" and sv not in parts:
                        parts.append(sv)
                labels.append("/".join(parts) if parts else f"{chrom}:{node_start[c]:,}-{node_end[c]:,}")
        else:
            labels = [f"{chrom}:{node_start[c]:,}-{node_end[c]:,}" for c in range(m)]

    pos = layout_edges(m, es, et, kind=layout, backend=backend, seed=seed, x=x)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if layout == "genomic":
        segs = [_arc(pos[a, 0], pos[b, 0]) for a, b in zip(es.tolist(), et.tolist())]
        if segs:
            ax.add_collection(LineCollection(segs, colors=edge_color, linewidths=widths, zorder=1))
        ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=colors, zorder=2, edgecolors="white", linewidths=0.5)
        span = max(end - start, 1)
        ax.set_xlim(start - 0.02 * span, end + 0.02 * span)
        top = max((abs(pos[b, 0] - pos[a, 0]) / 2.0 for a, b in zip(es.tolist(), et.tolist())), default=span / 10)
        ax.set_ylim(-0.15 * top, top * 1.15)
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        unit, div = ("Mb", 1e6) if span >= 1e6 else (("kb", 1e3) if span >= 1e3 else ("bp", 1))
        ax.set_xlabel(f"{chrom} ({unit})")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / div:,.2f}"))
        if labels is not None:
            for c in range(m):
                ax.annotate(labels[c], (pos[c, 0], 0), xytext=(0, -10), textcoords="offset points",
                            ha="center", va="top", fontsize=font_size, rotation=45, clip_on=True)
    else:
        segs = [[pos[a], pos[b]] for a, b in zip(es.tolist(), et.tolist())]
        if segs:
            ax.add_collection(LineCollection(segs, colors=edge_color, linewidths=widths, zorder=1))
        ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=colors, zorder=2, edgecolors="white", linewidths=0.5)
        if labels is not None:
            for c in range(m):
                ax.annotate(labels[c], pos[c], xytext=(0, 6 + np.sqrt(sizes[c]) / 2), textcoords="offset points",
                            ha="center", va="bottom", fontsize=font_size)
        pad = 0.15
        lo, hi = pos.min(0), pos.max(0)
        rng_ = np.maximum(hi - lo, 1e-9)
        ax.set_xlim(lo[0] - pad * rng_[0], hi[0] + pad * rng_[0])
        ax.set_ylim(lo[1] - pad * rng_[1], hi[1] + pad * rng_[1] + 0.1 * rng_[1])
        ax.set_aspect("equal", adjustable="datalim")
        ax.axis("off")
    if legend:
        from matplotlib.lines import Line2D
        ax.legend([Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8) for c in legend.values()],
                  list(legend), fontsize=font_size, loc="upper right", frameon=False)
    if title:
        merged = f" ({len(rows)} CREs merged at {merge_distance:,} bp)" if merge_distance is not None else ""
        ax.set_title(f"{A.name}: {chrom}:{start:,}-{end:,}\n{m} nodes{merged}, {len(uk)} edges", fontsize=11)
    return ax
