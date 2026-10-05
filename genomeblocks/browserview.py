"""IGV-like genomic region browser.

Lightweight, fully vectorial region viewer built on matplotlib. Given a
genomic region and a dict of named tracks, ``browser()`` lays out one
sub-axes per track (sharing x) plus a coordinate ruler, dispatching to the
appropriate drawer by track type:

    - bigWig path / open handle (or a list)  → binned coverage (a list is averaged
                                               into one track — replicate grouping)
    - ``.bam``                               → per-base coverage with IGV-style
                                               reference-mismatch colouring (needs ``reference``)
    - intervals (Loci, BED / narrowPeak path, a frame, a list of regions —
      anything :func:`~genomeblocks.as_loci` takes)  → interval rectangles
    - ``Pairs`` / ``.bedpe`` path             → half-sine arcs between anchors
    - ``Genes``                              → stacked gene models (exon / CDS)

All output is SVG-clean (no rasterised patches), and only binned values are
pulled from bigWig files so large regions stay cheap.
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from .locus import Locus, parse_region
from .genes import Genes
from .bedpe import Pairs, as_pairs


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

# IGV-standard nucleotide colours, used both for BAM mismatch bars and the
# reference sequence track (A green, C blue, G orange, T red, else gray).
_NUC_COLORS = {
    'A': '#009600',
    'C': '#0000ff',
    'G': '#d17105',
    'T': '#ff0000',
    'N': '#b0b0b0',
}


# ── type detection ───────────────────────────────────────────────────────────

def _is_bigwig(x) -> bool:
    from .backends.bigwig import handle_backend
    if isinstance(x, (str, os.PathLike)):
        return str(x).lower().removesuffix('.gz').endswith(('.bw', '.bigwig'))
    return handle_backend(x) is not None


def _detect_track_type(track: Any) -> str:
    """'genes', 'bedpe', 'bw', 'bam', 'narrowPeak' or 'bed' (anything as_loci takes)."""
    if isinstance(track, Genes):
        return 'genes'
    if isinstance(track, Pairs):
        return 'bedpe'
    if isinstance(track, (str, os.PathLike)):
        p = str(track).lower().removesuffix('.gz')
        if p.endswith(('.bw', '.bigwig')):
            return 'bw'
        if p.endswith('.bam'):
            return 'bam'
        if p.endswith('.bedpe'):
            return 'bedpe'
        if p.endswith('.narrowpeak'):
            return 'narrowPeak'
        return 'bed'                                    # a BED-like path or a region string
    if _is_bigwig(track):
        return 'bw'
    if isinstance(track, (list, tuple)) and track and all(_is_bigwig(t) for t in track):
        return 'bw'
    return 'bed'                                        # Loci, frames, lists of regions: as_loci decides


# ── drawers ──────────────────────────────────────────────────────────────────

def _draw_intervals(ax, track, chrom, start, end, color):
    """Draw an interval set as a single row of rectangles (columns, no Python loop)."""
    from .interop import as_loci
    L = as_loci(track)
    rows = L.overlap_rows(chrom, start, end)
    y0, y1 = 0.15, 0.85
    if len(rows):
        x0 = np.maximum(L.starts[rows], start)
        x1 = np.minimum(L.ends[rows], end)
        ok = x1 > x0
        patches = [Rectangle((a, y0), b - a, y1 - y0) for a, b in zip(x0[ok].tolist(), x1[ok].tolist())]
        if patches:
            ax.add_collection(PatchCollection(patches, facecolor=color, edgecolor='none'))
    ax.set_ylim(0, 1)
    ax.set_yticks([])


def _bigwig_bins(track, chrom, start, end, n_bins, backend=None) -> np.ndarray:
    """Mean per bin of one bigWig, or the average over a list of them."""
    from .backends.bigwig import open_bigwig
    bws = track if isinstance(track, (list, tuple)) else [track]
    ys = []
    for tr in bws:
        h = open_bigwig(tr, backend=backend)              # a path or an open handle (left open)
        try:
            ys.append(h.stats_array(chrom, start, end, n_bins=n_bins, stat='mean', missing=0.0))
        finally:
            h.close()
    return ys[0] if len(ys) == 1 else np.mean(ys, axis=0)


def _draw_bigwig(ax, track, chrom, start, end, color, n_bins, ymax, backend=None):
    """Draw a bigWig coverage track via binned stats (fill_between + outline).

    ``track`` may be a single bigWig (path or open handle) or a list of them —
    in which case their per-bin means are averaged into one track (replicate
    grouping), the same averaging the heatmap does across signal columns.
    """
    y = _bigwig_bins(track, chrom, start, end, n_bins, backend)
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
    """Draw a BAM coverage track with IGV-style mismatch colouring.

    Total per-base depth is a gray filled step (like the bigWig drawer); on
    top, positions where the non-reference allele fraction reaches
    ``allele_freq`` get the mismatching reads drawn as stacked, nucleotide-
    coloured bars (the matched fraction stays gray from the fill underneath).
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

    Base-level zoom (≤ 200 bp) shows coloured letters; a bit wider (≤ 5 kb)
    shows a coloured strip; wider still it stays blank.
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


def _draw_bedpe(ax, track, chrom, start, end, color, arc_points=40, max_arc_height=None):
    """Draw pairs as half-sine arcs connecting the two anchor midpoints.

    The y-axis is capped at ``max_arc_height`` (defaults to half the view
    span), giving IGV-like proportions: loops fully contained in the view fit
    the panel, while long-range loops are clipped at the panel top so their
    takeoff angle stays informative. Line width follows the score column.
    """
    P = as_pairs(track)
    span = max(1, end - start)
    cap = float(max_arc_height) if max_arc_height is not None else span * 0.5
    code = P.genome.code.get(chrom)
    if code is not None and len(P):
        cis = P.is_cis & (P.a.codes == code)
        m1 = np.minimum(P.mids1, P.mids2)[cis].astype(float)
        m2 = np.maximum(P.mids1, P.mids2)[cis].astype(float)
        sc = (np.asarray(P.cols['score'], float) if 'score' in P.cols else np.ones(len(P)))[cis]
        inview = (m2 >= start) & (m1 <= end) & (m2 > m1)
        m1, m2, sc = m1[inview], m2[inview], np.nan_to_num(sc[inview])
        if len(m1):
            theta = np.linspace(0.0, np.pi, arc_points)
            cx, r = 0.5 * (m1 + m2), 0.5 * (m2 - m1)
            xs = cx[:, None] + r[:, None] * np.cos(theta)[None, :]
            ys = r[:, None] * np.sin(theta)[None, :]
            segs = np.stack([xs, ys], axis=2)
            smax = sc.max()
            lw = 0.5 + 1.5 * (sc / smax if smax > 0 else np.zeros_like(sc))
            ax.add_collection(LineCollection(segs, colors=color, linewidths=lw, capstyle='round'))
    ax.set_ylim(0, cap * 1.05)
    ax.set_yticks([])


def _draw_genes(ax, genes, chrom, start, end, color, show_labels=True, max_transcripts_per_gene=None):
    """Draw a Genes object from its three tables: body line, exon boxes, taller
    CDS boxes, a label per gene.

    Set ``max_transcripts_per_gene=1`` to collapse each gene to one isoform
    (the canonical one when selected, else the longest); ``None`` shows every
    isoform.
    """
    G, T, F = genes.genes, genes.transcripts, genes.features
    gi = G.overlap_rows(chrom, start, end)
    ax.set_yticks([])
    if not len(gi):
        ax.set_ylim(0, 1)
        return
    gi = gi[np.lexsort((-G.ends[gi], G.starts[gi]))]
    gene_rank = {g: k for k, g in enumerate(gi.tolist())}
    tx = np.flatnonzero(np.isin(T.cols['gene'], gi))
    if not len(tx):
        ax.set_ylim(0, 1)
        return
    tlen = (T.ends - T.starts)[tx]
    is_canon = np.zeros(len(T), bool)
    canon = G.cols.get('canonical')
    if canon is not None:
        c = np.asarray(canon)
        is_canon[c[c >= 0]] = True
    g_rank = np.array([gene_rank[g] for g in T.cols['gene'][tx].tolist()])
    order = np.lexsort((-tlen, ~is_canon[tx], g_rank))         # gene, canonical first, then longest
    tx = tx[order]
    if max_transcripts_per_gene is not None:
        g_of = T.cols['gene'][tx]
        first = np.r_[True, g_of[1:] != g_of[:-1]]
        starts_of_run = np.where(first, np.arange(len(g_of)), 0)
        rank = np.arange(len(g_of)) - np.maximum.accumulate(starts_of_run)
        tx = tx[rank < max_transcripts_per_gene]

    # greedy stacking: a transcript goes to the first row whose rightmost end leaves a gap
    gap = (end - start) * 0.1
    rows_end: list = []
    row_of = np.empty(len(tx), np.int64)
    for k, (ts, te) in enumerate(zip(T.starts[tx].tolist(), T.ends[tx].tolist())):
        for i, rend in enumerate(rows_end):
            if ts > rend + gap:
                rows_end[i] = te
                row_of[k] = i
                break
        else:
            row_of[k] = len(rows_end)
            rows_end.append(te)
    n_rows = max(1, len(rows_end))
    yc = (n_rows - 1 - row_of) + 0.5                              # row 0 on top
    row_of_tx = np.full(len(T), -1, np.int64)
    row_of_tx[tx] = np.arange(len(tx))

    for k, t in enumerate(tx.tolist()):
        gx0, gx1 = max(int(T.starts[t]), start), min(int(T.ends[t]), end)
        ax.plot([gx0, gx1], [yc[k], yc[k]], color=color, linewidth=1, zorder=1)

    fm = row_of_tx[F.cols['transcript']] >= 0
    for kind, half, z in ((0, 0.12, 2), (1, 0.24, 3)):            # exons thin, CDS tall
        sel = fm & (F.cols['kind'] == kind)
        if not sel.any():
            continue
        x0 = np.maximum(F.starts[sel], start)
        x1 = np.minimum(F.ends[sel], end)
        yk = yc[row_of_tx[F.cols['transcript'][sel]]]
        ok = x1 > x0
        patches = [Rectangle((a, y - half), b - a, 2 * half)
                   for a, b, y in zip(x0[ok].tolist(), x1[ok].tolist(), yk[ok].tolist())]
        if patches:
            ax.add_collection(PatchCollection(patches, facecolor=color, edgecolor='none', zorder=z))

    if show_labels:
        names = G.cols['gene_name']
        ids = G.cols['gene_id']
        g_of = T.cols['gene'][tx]
        first = np.r_[True, g_of[1:] != g_of[:-1]]
        for k in np.flatnonzero(first).tolist():
            g = int(g_of[k])
            lbl = str(names[g] or ids[g])
            arrow = '→' if G.strands[g] == 1 else ('←' if G.strands[g] == 2 else '')
            ax.text(max(int(G.starts[g]), start) - 100, yc[k] + 0.14, f"{arrow} {lbl}".strip(),
                    fontsize=7, color='#222', ha='left', va='bottom', clip_on=True)
    ax.set_ylim(0, n_rows + 1)


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
    backend: Optional[str] = None,
):
    """Plot an IGV-like browser view of a genomic region.

    Parameters
    ----------
    region : Locus | (chrom, start, end) | 'chr1:1,000-2,000' | 'chr8:127.7-128.1 Mb'
        The region to display.
    tracks : dict[str, Any]
        Ordered mapping of track name → track source. Type is auto-detected:
        bigWig paths or open handles (a list is averaged), ``.bam`` paths,
        ``.bedpe`` paths or ``Pairs``, ``Genes``, and any interval set
        :func:`~genomeblocks.as_loci` takes (Loci, BED / narrowPeak paths,
        frames, lists of regions). ``.bam`` files need a coordinate-sorted BAM
        with a ``.bai`` index alongside.
    figsize : (w, h)
        Figure size. If ``h`` is None it is derived from track heights.
    track_heights, colors : dict, optional
        Per-track height / colour overrides. Defaults depend on track type.
    bw_n_bins : int
        Number of bins to request per bigWig track.
    bw_ymax : float | dict[str, float], optional
        Y-axis maximum for bigWig tracks (scalar for all, dict per track).
    bw_share : list[list[str]], optional
        Groups of bigWig track names that share one y-scale (the group's
        region max), e.g. ``[["AR 0h", "AR 4h"]]``. An explicit ``bw_ymax`` wins.
    reference : str, optional
        Path to an indexed FASTA (``.fa`` + ``.fai``). Required for ``.bam``
        mismatch colouring; also drives the reference-sequence track.
    show_sequence : bool
        When ``reference`` is given, append a reference-sequence track at the
        bottom (letters ≤ 200 bp, a colour strip ≤ 5 kb). Axes at ``'_sequence'``.
    bam_min_baseq, bam_allele_freq, bam_ymax, bam_share
        BAM coverage settings (IGV defaults 15 and 0.2); scaling as for bigWigs.
    genes_max_transcripts : int, optional
        Isoforms shown per gene (``1`` collapses to the canonical / longest).
    backend : str, optional
        bigWig engine ('pybigtools', 'pybigwig', 'python'); default automatic.

    Returns
    -------
    (fig, axes_by_name) : (matplotlib.figure.Figure, dict[str, Axes])
        ``axes_by_name`` includes an ``'_axis'`` entry for the ruler and, when a
        reference sequence track is drawn, a ``'_sequence'`` entry.
    """
    chrom, start, end = parse_region(region)
    if end <= start:
        raise ValueError(f"Invalid region: end ({end}) must be > start ({start})")

    types = {name: _detect_track_type(tr) for name, tr in tracks.items()}
    if any(t == 'bam' for t in types.values()) and reference is None:
        raise ValueError("BAM tracks need a `reference` FASTA for mismatch "
                         "colouring; pass reference='genome.fa'.")

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
    gs = gridspec.GridSpec(len(tracks) + len(tail_heights), 1, height_ratios=heights + tail_heights,
                           hspace=hspace)

    colors_cfg = colors or {}
    # per-bigWig y-limits: scalar / dict / shared groups; an explicit bw_ymax wins
    bw_names = [n for n in tracks if types[n] == 'bw']
    if isinstance(bw_ymax, (int, float)):
        ymax_cfg = {n: float(bw_ymax) for n in bw_names}
    else:
        ymax_cfg = dict(bw_ymax or {})
    if bw_share:
        for grp in bw_share:
            members = [n for n in grp if n in tracks and types[n] == 'bw']
            if not members:
                continue
            gmax = max(float(np.nan_to_num(_bigwig_bins(tracks[n], chrom, start, end, bw_n_bins, backend)).max())
                       for n in members) * 1.05
            for n in members:
                ymax_cfg.setdefault(n, gmax)

    bam_names = [n for n in tracks if types[n] == 'bam']
    if isinstance(bam_ymax, (int, float)):
        bam_ymax_cfg = {n: float(bam_ymax) for n in bam_names}
    else:
        bam_ymax_cfg = dict(bam_ymax or {})
    if bam_share:
        from .bam import pileup_counts

        def _bam_region_max(track):
            return float(pileup_counts(track, chrom, start, end, min_baseq=bam_min_baseq).sum(axis=0).max())

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
            _draw_bigwig(ax, track, chrom, start, end, color=col, n_bins=bw_n_bins, ymax=ymax_cfg.get(name),
                         backend=backend)
        elif tt == 'bam':
            _draw_bam(ax, track, chrom, start, end, color=col, reference=reference, min_baseq=bam_min_baseq,
                      allele_freq=bam_allele_freq, ymax=bam_ymax_cfg.get(name))
        elif tt == 'bedpe':
            _draw_bedpe(ax, track, chrom, start, end, color=col)
        elif tt == 'genes':
            _draw_genes(ax, track, chrom, start, end, color=col, max_transcripts_per_gene=genes_max_transcripts)

        ax.set_xlim(start, end)
        ax.set_ylabel(name, rotation=0, ha='right', va='center', fontsize=label_fontsize, labelpad=10)
        for side in ('top', 'right', 'bottom'):
            ax.spines[side].set_visible(False)
        if tt not in ('bw', 'bam'):
            ax.spines['left'].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    if seq_on:
        ax_seq = fig.add_subplot(gs[len(tracks)], sharex=first_ax)
        axes_by_name['_sequence'] = ax_seq
        _draw_sequence(ax_seq, chrom, start, end, reference)
        ax_seq.set_xlim(start, end)
        ax_seq.set_ylabel('sequence', rotation=0, ha='right', va='center', fontsize=label_fontsize, labelpad=10)
        for side in ('top', 'right', 'bottom', 'left'):
            ax_seq.spines[side].set_visible(False)
        ax_seq.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    ax_ruler = fig.add_subplot(gs[-1], sharex=first_ax)
    axes_by_name['_axis'] = ax_ruler
    ax_ruler.set_xlim(start, end)
    ax_ruler.set_ylim(0, 1)
    ax_ruler.set_yticks([])
    for side in ('top', 'right', 'left'):
        ax_ruler.spines[side].set_visible(False)
    ax_ruler.tick_params(axis='x', labelsize=7, length=3)
    _format_ruler(ax_ruler, start, end)
    ax_ruler.set_xlabel(f"{chrom}:{start:,}-{end:,}", fontsize=label_fontsize)

    return fig, axes_by_name
