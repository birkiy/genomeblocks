"""IGV-like genomic region browser.

Lightweight, fully vectorial region viewer built on matplotlib.  Given a
genomic region and a dict of named tracks, ``browser()`` lays out one
sub-axes per track (sharing x) plus a coordinate ruler, dispatching to the
appropriate drawer by track type:

    - ``.bw`` / ``.bigwig`` (or a list)   → binned coverage (a list of bigwigs
                                            is averaged into one track — replicate grouping)
    - ``.narrowPeak`` / ``.bed`` / Loci   → interval rectangles
    - ``.bedpe`` / list[Pair]             → half-sine arcs between anchors
    - ``Genes``                           → stacked gene models (exon/CDS)

All output is SVG-clean (no rasterized patches), and only binned values
are pulled from bigwig files so large regions stay cheap.
"""
from __future__ import annotations
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.ticker import FuncFormatter

from .locus import Locus
from .loci import Loci
from .genes import Genes
from .bedpe import Pair, read_bedpe


# ── defaults ─────────────────────────────────────────────────────────────────

_DEFAULT_HEIGHTS = {
    'bed':        0.2,
    'narrowPeak': 0.2,
    'bw':         0.5,
    'bedpe':      1.5,
    'genes':      1.5,
}

_DEFAULT_COLORS = {
    'bed':        '#444444',
    'narrowPeak': '#444444',
    'bw':         '#4c78a8',
    'bedpe':      '#888888',
    'genes':      '#000000',
}


# ── region parsing ───────────────────────────────────────────────────────────

def _parse_region(region) -> Tuple[str, int, int]:
    if isinstance(region, Locus):
        return region.chrom, int(region.start), int(region.end)
    if isinstance(region, (tuple, list)) and len(region) == 3:
        return str(region[0]), int(region[1]), int(region[2])
    if isinstance(region, str):
        chrom, rest = region.split(':')
        rest = rest.replace(',', '').replace(' ', '')
        s, e = rest.split('-')
        return chrom, int(s), int(e)
    raise ValueError(f"Cannot parse region: {region!r}")


# ── type detection ───────────────────────────────────────────────────────────

def _detect_track_type(track: Any) -> str:
    if isinstance(track, str):
        p = track.lower()
        if p.endswith(('.bw', '.bigwig')):  return 'bw'
        if p.endswith('.narrowpeak'):       return 'narrowPeak'
        if p.endswith('.bed'):              return 'bed'
        if p.endswith('.bedpe'):            return 'bedpe'
        raise ValueError(f"Unknown track extension: {track}")
    if isinstance(track, Genes): return 'genes'
    if isinstance(track, Loci):  return 'bed'
    if isinstance(track, (list, tuple)) and track:
        # list[Pair] → bedpe; list of bigwig paths → one averaged bw track
        # (replicate grouping), mirroring how the heatmap averages columns.
        if isinstance(track[0], Pair):
            return 'bedpe'
        if all(isinstance(t, str) and t.lower().endswith(('.bw', '.bigwig')) for t in track):
            return 'bw'
    raise ValueError(f"Cannot detect track type for: {type(track).__name__}")


# ── drawers ──────────────────────────────────────────────────────────────────

def _draw_intervals(ax, track, chrom, start, end, color):
    """Draw Loci / bed / narrowPeak as a single row of rectangles."""
    loci = Loci.make(track) if isinstance(track, str) else track
    y0, y1 = 0.15, 0.85
    patches = []
    for l in loci:
        if l.chrom != chrom: continue
        if l.end < start or l.start > end: continue
        x0 = max(l.start, start)
        x1 = min(l.end, end)
        if x1 <= x0: continue
        patches.append(Rectangle((x0, y0), x1 - x0, y1 - y0))
    if patches:
        ax.add_collection(PatchCollection(patches, facecolor=color, edgecolor='none'))
    ax.set_ylim(0, 1)
    ax.set_yticks([])


def _draw_bigwig(ax, track, chrom, start, end, color, n_bins, ymax):
    """Draw a bigwig coverage track via binned stats (fill_between + outline).

    ``track`` may be a single bigwig (path or handle) or a list of bigwigs —
    in which case their per-bin means are averaged into one track (replicate
    grouping), the same averaging the heatmap does across signal columns.
    """
    from .signal import _bw_open  # lazy — avoid import cost if no bw track
    bws = track if isinstance(track, (list, tuple)) else [track]
    ys = []
    for tr in bws:
        opened = isinstance(tr, str)
        h = _bw_open(tr) if opened else tr
        try:
            ys.append(h.stats_array(chrom, start, end, n_bins=n_bins,
                                    stat='mean', missing=0.0).astype(np.float64, copy=False))
        finally:
            if opened:
                h.close()
    y = ys[0] if len(ys) == 1 else np.mean(ys, axis=0)
    x = np.linspace(start, end, n_bins, endpoint=False) + (end - start) / (2 * n_bins)

    ax.fill_between(x, 0.0, y, facecolor=color, linewidth=0, step='mid')
    ax.plot(x, y, color=color, linewidth=0.4, drawstyle='steps-mid')

    peak = float(ymax) if ymax is not None else (float(y.max()) * 1.05 if y.size and y.max() > 0 else 1.0)
    ax.set_ylim(0, peak)
    ax.set_yticks([0, peak])
    ax.tick_params(axis='y', labelsize=6, length=2, pad=1)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def _draw_bedpe(ax, track, chrom, start, end, color, arc_points=40,
                max_arc_height=None):
    """Draw bedpe as half-sine arcs connecting the two anchor midpoints.

    The y-axis is capped at ``max_arc_height`` (defaults to half the view
    span), giving IGV-like proportions: loops fully contained in the view
    fit the panel, while long-range loops are clipped at the panel top so
    their takeoff angle stays informative.
    """
    pairs = read_bedpe(track, verbose=False) if isinstance(track, str) else track
    span = max(1, end - start)
    cap = float(max_arc_height) if max_arc_height is not None else span * 0.5

    scores = [p.score for p in pairs if p.chrom1 == chrom and p.chrom2 == chrom]
    smax = max(scores) if scores else 0.0

    theta = np.linspace(0.0, np.pi, arc_points)
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    for p in pairs:
        if p.chrom1 != chrom or p.chrom2 != chrom: continue
        m1, m2 = sorted((p.mid1, p.mid2))
        if m2 < start or m1 > end: continue
        cx = 0.5 * (m1 + m2)
        r = 0.5 * (m2 - m1)
        if r <= 0: continue
        # Peak proportional to anchor span; panel clips anything > cap.
        xs = cx + r * cos_t
        ys = r * sin_t
        lw = 0.5 + 1.5 * (p.score / smax if smax > 0 else 0.0)
        ax.plot(xs, ys, color=color, linewidth=lw,
                solid_capstyle='round', clip_on=True)

    ax.set_ylim(0, cap * 1.05)
    ax.set_yticks([])


def _stack_genes(visible_with_tx, gap_frac, span):
    """Assign each shown transcript to a stacking row to avoid overlap.

    ``visible_with_tx`` is a list of ``(gene, [transcript, ...])`` — caller
    pre-filters transcripts (e.g. cap at k per gene) before stacking.

    Returns ``(row_of_transcript_id, n_rows)``.
    """
    gap = span * gap_frac
    rows = []  # rightmost end per row
    row_of = {}
    for _, tx_list in visible_with_tx:
        for t in tx_list:
            placed = False
            for i, rend in enumerate(rows):
                if t.start > rend + gap:
                    rows[i] = t.end
                    row_of[id(t)] = i
                    placed = True
                    break
            if not placed:
                row_of[id(t)] = len(rows)
                rows.append(t.end)
    return row_of, max(1, len(rows))


def _select_transcripts(gene, max_per_gene):
    """Pick up to ``max_per_gene`` transcripts to display, longest-first by
    genomic span (end − start). ``max_per_gene=None`` keeps all; ``1``
    collapses each gene to its longest isoform."""
    tx_list = list((gene.transcripts or {}).values())
    if max_per_gene is None or len(tx_list) <= max_per_gene:
        return tx_list
    tx_list.sort(key=lambda t: t.end - t.start, reverse=True)
    return tx_list[:max_per_gene]


def _draw_genes(ax, genes, chrom, start, end, color, show_labels=True,
                max_transcripts_per_gene=None):
    """Draw a Genes object: body line, exon boxes, taller CDS boxes, label.

    Set ``max_transcripts_per_gene=1`` to collapse each gene to its longest
    transcript (cleaner view for dense regions); ``None`` (default) shows
    every isoform.
    """
    visible = [g for g in genes.values()
               if g.chrom == chrom and g.end >= start and g.start <= end]
    visible.sort(key=lambda g: (g.start, -g.end))

    visible_with_tx = [(g, _select_transcripts(g, max_transcripts_per_gene))
                       for g in visible]
    row_of, n_rows = _stack_genes(visible_with_tx, gap_frac=0.1, span=end - start)

    exon_boxes, cds_boxes = [], []

    for g, tx_list in visible_with_tx:
        for t in tx_list:

            r = row_of[id(t)]
            yc = (n_rows - 1 - r) + 0.5  # flipped: row 0 on top

            # Thin body line (intron backbone) across full gene span
            gx0 = max(t.start, start)
            gx1 = min(t.end, end)
            ax.plot([gx0, gx1], [yc, yc], color=color, linewidth=1, zorder=1)


            for e in t.exons:
                if e.end < start or e.start > end: continue
                x0 = max(e.start, start)
                w = min(e.end, end) - x0
                if w <= 0: continue
                exon_boxes.append(Rectangle((x0, yc - 0.12), w, 0.24))

            for c in t.cds:
                if c.end < start or c.start > end: continue
                x0 = max(c.start, start)
                w = min(c.end, end) - x0
                if w <= 0: continue
                cds_boxes.append(Rectangle((x0, yc - 0.24), w, 0.48))

        if show_labels and tx_list:
            lbl = g.gene_name or g.gene_id
            lbl_x = max(g.start, start)
            arrow = '→' if g.strand == '+' else ('←' if g.strand == '-' else '')
            ax.text(lbl_x-100, yc + 0.14, f"{arrow} {lbl}".strip(),
                    fontsize=7, color='#222', ha='left', va='bottom',
                    clip_on=True)

    if exon_boxes:
        ax.add_collection(PatchCollection(exon_boxes, facecolor=color, edgecolor='none', zorder=2))
    if cds_boxes:
        ax.add_collection(PatchCollection(cds_boxes,  facecolor=color, edgecolor='none', zorder=3))

    ax.set_ylim(0, n_rows+1)
    ax.set_yticks([])


# ── coordinate ruler ─────────────────────────────────────────────────────────

def _format_ruler(ax, start, end):
    span = end - start
    if span >= 1_000_000:
        unit, div = 'Mb', 1_000_000
    elif span >= 1_000:
        unit, div = 'kb', 1_000
    else:
        unit, div = 'bp', 1
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v/div:,.2f} {unit}"))


# ── public entry point ──────────────────────────────────────────────────────

def browser(
    region: Union[Locus, Tuple[str, int, int], str],
    tracks: Dict[str, Any],
    *,
    figsize: Tuple[float, Optional[float]] = (10, None),
    dpi: int = 100,
    track_heights: Optional[Dict[str, float]] = None,
    colors: Optional[Dict[str, str]] = None,
    bw_n_bins: int = 1000,
    bw_ymax: Optional[Dict[str, float]] = None,
    label_fontsize: int = 8,
    hspace: float = 0.15,
    genes_max_transcripts: Optional[int] = None,
):
    """Plot an IGV-like browser view of a genomic region.

    Parameters
    ----------
    region : Locus | (chrom, start, end) | 'chr1:1,000-2,000'
        The region to display.
    tracks : dict[str, Any]
        Ordered mapping of track name → track source.  Type is auto-detected:
        path strings dispatch by extension (``.bw``, ``.narrowPeak``, ``.bed``,
        ``.bedpe``); ``Loci`` / ``Genes`` / ``list[Pair]`` objects are also
        accepted directly.
    figsize : (w, h)
        Figure size.  If ``h`` is None it is derived from track heights.
    dpi : int
    track_heights : dict[str, float], optional
        Per-track height overrides.  Defaults depend on track type.
    colors : dict[str, str], optional
        Per-track color overrides.
    bw_n_bins : int
        Number of bins to request per bigwig track.
    bw_ymax : dict[str, float], optional
        Per-track y-axis maximum for bigwig tracks (auto-scaled otherwise).
    label_fontsize : int
        Y-label font size.
    hspace : float
        Vertical spacing between tracks.

    Returns
    -------
    (fig, axes_by_name) : (matplotlib.figure.Figure, dict[str, Axes])
        ``axes_by_name`` includes an ``'_axis'`` entry for the ruler.
    """
    chrom, start, end = _parse_region(region)
    if end <= start:
        raise ValueError(f"Invalid region: end ({end}) must be > start ({start})")

    types = {name: _detect_track_type(tr) for name, tr in tracks.items()}
    heights_cfg = track_heights or {}
    heights = [heights_cfg.get(n, _DEFAULT_HEIGHTS[types[n]]) for n in tracks]
    ruler_h = 0.25

    w, h = figsize
    if h is None:
        h = max(2.0, sum(heights) + ruler_h + 0.4)

    fig = plt.figure(figsize=(w, h), dpi=dpi)
    gs = gridspec.GridSpec(
        len(tracks) + 1, 1,
        height_ratios=heights + [ruler_h],
        hspace=hspace,
    )

    colors_cfg = colors or {}
    ymax_cfg = bw_ymax or {}
    axes_by_name: Dict[str, plt.Axes] = {}
    first_ax = None

    for i, (name, track) in enumerate(tracks.items()):
        ax = fig.add_subplot(gs[i], sharex=first_ax)
        if first_ax is None:
            first_ax = ax
        axes_by_name[name] = ax
        tt = types[name]
        col = colors_cfg.get(name, _DEFAULT_COLORS[tt])

        if tt in ('bed', 'narrowPeak'):
            _draw_intervals(ax, track, chrom, start, end, color=col)
        elif tt == 'bw':
            _draw_bigwig(ax, track, chrom, start, end, color=col,
                         n_bins=bw_n_bins, ymax=ymax_cfg.get(name))
        elif tt == 'bedpe':
            _draw_bedpe(ax, track, chrom, start, end, color=col)
        elif tt == 'genes':
            _draw_genes(ax, track, chrom, start, end, color=col,
                        max_transcripts_per_gene=genes_max_transcripts)

        ax.set_xlim(start, end)
        ax.set_ylabel(name, rotation=0, ha='right', va='center',
                      fontsize=label_fontsize, labelpad=10)
        for side in ('top', 'right', 'bottom'):
            ax.spines[side].set_visible(False)
        if tt != 'bw':
            ax.spines['left'].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    # Coordinate ruler along the bottom
    ax_ruler = fig.add_subplot(gs[-1], sharex=first_ax)
    axes_by_name['_axis'] = ax_ruler
    ax_ruler.set_xlim(start, end)
    ax_ruler.set_ylim(0, 1)
    ax_ruler.set_yticks([])
    for side in ('top', 'right', 'left'):
        ax_ruler.spines[side].set_visible(False)
    ax_ruler.tick_params(axis='x', labelsize=7, length=3)
    _format_ruler(ax_ruler, start, end)
    ax_ruler.set_xlabel(f"{chrom}:{start:,}-{end:,}",
                        fontsize=label_fontsize)

    return fig, axes_by_name
