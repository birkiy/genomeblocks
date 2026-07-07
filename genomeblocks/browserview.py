"""IGV-like genomic region browser.

Lightweight, fully vectorial region viewer built on matplotlib.  Given a
genomic region and a dict of named tracks, ``browser()`` lays out one
sub-axes per track (sharing x) plus a coordinate ruler, dispatching to the
appropriate drawer by track type:

    - ``.bw`` / ``.bigwig`` (or a list)   → binned coverage (a list of bigwigs
                                            is averaged into one track — replicate grouping)
    - ``.bam``                            → per-base coverage with IGV-style
                                            reference-mismatch coloring (needs ``reference``)
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
from collections import defaultdict

from .locus import Locus
from .loci import Loci
from .genes import Genes
from .bedpe import Pair, read_bedpe


# ── defaults ─────────────────────────────────────────────────────────────────

_DEFAULT_HEIGHTS = {
    'bed':        0.2,
    'narrowPeak': 0.2,
    'bw':         0.5,
    'bam':        0.6,
    'bedpe':      1.5,
    'genes':      1.5,
    'sequence':   0.2,
}

_DEFAULT_COLORS = {
    'bed':        '#444444',
    'narrowPeak': '#444444',
    'bw':         '#4c78a8',
    'bam':        '#a6a6a6',
    'bedpe':      '#888888',
    'genes':      '#000000',
}

# IGV-standard nucleotide colors, used both for BAM mismatch bars and the
# reference sequence track (A green, C blue, G orange, T red, else gray).
_NUC_COLORS = {
    'A': '#009600',
    'C': '#0000ff',
    'G': '#d17105',
    'T': '#ff0000',
    'N': '#b0b0b0',
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
        if p.endswith('.bam'):              return 'bam'
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


def _draw_bam(ax, track, chrom, start, end, color, *, reference, min_baseq,
              allele_freq, ymax):
    """Draw a BAM coverage track with IGV-style mismatch coloring.

    Total per-base depth is a gray filled step (like the bigwig drawer); on
    top, positions where the non-reference allele fraction reaches
    ``allele_freq`` get the mismatching reads drawn as stacked, nucleotide-
    colored bars (the matched fraction stays gray from the fill underneath).
    Mismatch detection is vectorized and only the handful of passing positions
    are rendered, so wide views stay cheap.
    """
    from .bam import pileup_counts, reference_seq

    counts = pileup_counts(track, chrom, start, end, min_baseq=min_baseq)  # (4, n)
    total = counts.sum(axis=0)
    n = total.size
    xc = np.arange(start, end) + 0.5
    ax.fill_between(xc, 0.0, total, facecolor=color, linewidth=0, step='mid')

    if reference is not None and n:
        ref = reference_seq(reference, chrom, start, end)
        ref = (ref + 'N' * n)[:n]                       # pad short edges with N
        b2i = np.full(256, -1, dtype=np.int64)
        for i, b in enumerate(b'ACGT'):
            b2i[b] = i
        ref_idx = b2i[np.frombuffer(ref.encode('ascii'), dtype=np.uint8)]
        matched = np.zeros(n, dtype=np.int64)
        valid = ref_idx >= 0
        matched[valid] = counts[ref_idx[valid], np.nonzero(valid)[0]]
        mism = total - matched
        with np.errstate(invalid='ignore', divide='ignore'):
            frac = np.where(total > 0, mism / total, 0.0)
        boxes = defaultdict(list)
        for p in np.nonzero((total > 0) & (frac >= allele_freq))[0]:
            y = 0.0
            for bi, base in enumerate('ACGT'):
                if bi == ref_idx[p]:
                    continue
                c = counts[bi, p]
                if c <= 0:
                    continue
                boxes[base].append(Rectangle((start + p, y), 1.0, c))
                y += c
        for base, patches in boxes.items():
            ax.add_collection(PatchCollection(patches, facecolor=_NUC_COLORS[base],
                                              edgecolor='none', zorder=3))

    peak = float(ymax) if ymax is not None else (
        float(total.max()) * 1.05 if n and total.max() > 0 else 1.0)
    ax.set_ylim(0, peak)
    ax.set_yticks([0, peak])
    ax.tick_params(axis='y', labelsize=6, length=2, pad=1)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def _draw_sequence(ax, chrom, start, end, reference):
    """Draw the reference bases beneath the tracks, IGV-style.

    Base-level zoom (≤ 200 bp) shows colored letters; a bit wider (≤ 5 kb)
    shows a colored strip; wider still it stays blank (letters/strip would be
    illegible and expensive).
    """
    from .bam import reference_seq
    seq = reference_seq(reference, chrom, start, end)
    span = end - start
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    if span <= 0 or not seq:
        return
    if span <= 200:
        fs = min(8.0, max(3.0, 900.0 / span))
        for i, b in enumerate(seq):
            ax.text(start + i + 0.5, 0.5, b, ha='center', va='center',
                    fontsize=fs, family='monospace', clip_on=True,
                    color=_NUC_COLORS.get(b, '#888888'))
    elif span <= 5000:
        boxes = defaultdict(list)
        for i, b in enumerate(seq):
            boxes[b].append(Rectangle((start + i, 0.15), 1.0, 0.7))
        for b, patches in boxes.items():
            ax.add_collection(PatchCollection(patches, facecolor=_NUC_COLORS.get(b, '#888888'),
                                              edgecolor='none'))


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
    bw_ymax: Union[float, Dict[str, float], None] = None,
    bw_share: Optional[Sequence[Sequence[str]]] = None,
    reference: Optional[str] = None,
    show_sequence: bool = True,
    bam_min_baseq: int = 15,
    bam_allele_freq: float = 0.2,
    bam_ymax: Union[float, Dict[str, float], None] = None,
    bam_share: Optional[Sequence[Sequence[str]]] = None,
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
        path strings dispatch by extension (``.bw``, ``.bam``, ``.narrowPeak``,
        ``.bed``, ``.bedpe``); ``Loci`` / ``Genes`` / ``list[Pair]`` objects are
        also accepted directly.  ``.bam`` files need a coordinate-sorted BAM
        with a ``.bai`` index alongside.
    figsize : (w, h)
        Figure size.  If ``h`` is None it is derived from track heights.
    dpi : int
    track_heights : dict[str, float], optional
        Per-track height overrides.  Defaults depend on track type.
    colors : dict[str, str], optional
        Per-track color overrides.
    bw_n_bins : int
        Number of bins to request per bigwig track.
    bw_ymax : float | dict[str, float], optional
        Y-axis maximum for bigwig tracks. A scalar applies to every bigwig
        track; a dict sets it per track. Auto-scaled (per track) otherwise.
    bw_share : list[list[str]], optional
        Groups of bigwig track names that should share one y-scale (the group's
        region max), so tracks are directly comparable — e.g.
        ``[["AR 0h", "AR 4h"]]`` scales both AR tracks together. An explicit
        ``bw_ymax`` for a track still takes precedence.
    reference : str, optional
        Path to an indexed FASTA (``.fa`` + ``.fai``).  Required for ``.bam``
        mismatch coloring; also drives the reference-sequence track drawn just
        above the ruler (see ``show_sequence``).
    show_sequence : bool
        When ``reference`` is given, append a reference-sequence track at the
        bottom.  Shows colored letters at base-level zoom (≤ 200 bp), a color
        strip up to 5 kb, blank beyond.  Access its axes at ``'_sequence'``.
    bam_min_baseq : int
        Minimum base quality for a read to count toward BAM coverage (matches
        IGV's default of 15).
    bam_allele_freq : float
        Fraction of non-reference reads at a position before it is colored as a
        mismatch (IGV's ``0.2`` default); below this the bar stays fully gray.
    bam_ymax : float | dict[str, float], optional
        Coverage y-axis maximum for BAM tracks; scalar for all, or per-track
        dict.  Auto-scaled per track otherwise.
    bam_share : list[list[str]], optional
        Groups of BAM track names sharing one coverage y-scale (their region
        max), so depths are directly comparable — e.g. one group per condition.
        An explicit ``bam_ymax`` for a track still takes precedence.
    label_fontsize : int
        Y-label font size.
    hspace : float
        Vertical spacing between tracks.

    Returns
    -------
    (fig, axes_by_name) : (matplotlib.figure.Figure, dict[str, Axes])
        ``axes_by_name`` includes an ``'_axis'`` entry for the ruler and, when a
        reference sequence track is drawn, a ``'_sequence'`` entry.
    """
    chrom, start, end = _parse_region(region)
    if end <= start:
        raise ValueError(f"Invalid region: end ({end}) must be > start ({start})")

    types = {name: _detect_track_type(tr) for name, tr in tracks.items()}
    if any(t == 'bam' for t in types.values()) and reference is None:
        raise ValueError("BAM tracks need a `reference` FASTA for mismatch "
                         "coloring; pass reference='genome.fa'.")

    heights_cfg = track_heights or {}
    heights = [heights_cfg.get(n, _DEFAULT_HEIGHTS[types[n]]) for n in tracks]
    ruler_h = 0.25
    seq_on = reference is not None and show_sequence
    seq_h = heights_cfg.get('_sequence', _DEFAULT_HEIGHTS['sequence'])
    tail_heights = ([seq_h] if seq_on else []) + [ruler_h]

    w, h = figsize
    if h is None:
        h = max(2.0, sum(heights) + sum(tail_heights) + 0.4)

    fig = plt.figure(figsize=(w, h), dpi=dpi)
    gs = gridspec.GridSpec(
        len(tracks) + len(tail_heights), 1,
        height_ratios=heights + tail_heights,
        hspace=hspace,
    )

    colors_cfg = colors or {}
    # Resolve per-bigwig y-limits. bw_ymax may be a scalar (all bw tracks) or a
    # per-track dict; bw_share lists groups that share one y-scale (their region
    # max) so e.g. 0h/4h tracks are directly comparable. Explicit bw_ymax wins.
    bw_names = [n for n in tracks if types[n] == 'bw']
    if isinstance(bw_ymax, (int, float)):
        ymax_cfg = {n: float(bw_ymax) for n in bw_names}
    else:
        ymax_cfg = dict(bw_ymax or {})
    if bw_share:
        from .signal import _bw_open

        def _region_max(track):
            bws = track if isinstance(track, (list, tuple)) else [track]
            m = 0.0
            for tr in bws:
                h = _bw_open(tr)
                try:
                    v = h.stats_array(chrom, start, end, n_bins=bw_n_bins,
                                      stat='mean', missing=0.0)
                    m = max(m, float(np.nan_to_num(v).max()))
                finally:
                    h.close()
            return m

        for grp in bw_share:
            members = [n for n in grp if n in tracks and types[n] == 'bw']
            if not members:
                continue
            gmax = max(_region_max(tracks[n]) for n in members) * 1.05
            for n in members:
                ymax_cfg.setdefault(n, gmax)

    # Resolve per-BAM coverage y-limits, same scalar/dict/share logic as bigwig.
    bam_names = [n for n in tracks if types[n] == 'bam']
    if isinstance(bam_ymax, (int, float)):
        bam_ymax_cfg = {n: float(bam_ymax) for n in bam_names}
    else:
        bam_ymax_cfg = dict(bam_ymax or {})
    if bam_share:
        from .bam import pileup_counts

        def _bam_region_max(track):
            return float(pileup_counts(track, chrom, start, end,
                                       min_baseq=bam_min_baseq).sum(axis=0).max())

        for grp in bam_share:
            members = [n for n in grp if n in tracks and types[n] == 'bam']
            if not members:
                continue
            gmax = max(_bam_region_max(tracks[n]) for n in members) * 1.05
            for n in members:
                bam_ymax_cfg.setdefault(n, gmax)

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
        elif tt == 'bam':
            _draw_bam(ax, track, chrom, start, end, color=col,
                      reference=reference, min_baseq=bam_min_baseq,
                      allele_freq=bam_allele_freq, ymax=bam_ymax_cfg.get(name))
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
        if tt not in ('bw', 'bam'):
            ax.spines['left'].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    # Reference sequence track (letters or color strip) just above the ruler.
    if seq_on:
        ax_seq = fig.add_subplot(gs[len(tracks)], sharex=first_ax)
        axes_by_name['_sequence'] = ax_seq
        _draw_sequence(ax_seq, chrom, start, end, reference)
        ax_seq.set_xlim(start, end)
        ax_seq.set_ylabel('sequence', rotation=0, ha='right', va='center',
                          fontsize=label_fontsize, labelpad=10)
        for side in ('top', 'right', 'bottom', 'left'):
            ax_seq.spines[side].set_visible(False)
        ax_seq.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

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
