---
title: browser
parent: API Reference
layout: default
nav_order: 7
---

# `genomeblocks.browserview`
{: .no_toc }

The `browser()` function lives in the `genomeblocks.browserview` module and is
re-exported at the top level as `genomeblocks.browser`. Import it either way:

```python
from genomeblocks import browser            # the function (recommended)
from genomeblocks.browserview import browser  # equivalent
```

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `browser(region, tracks, ...)`

```python
browser(
    region: Locus | tuple[str, int, int] | str,
    tracks: dict[str, Any],
    *,
    figsize: tuple[float, float | None] = (10, None),
    dpi: int = 100,
    track_heights: dict[str, float] | None = None,
    colors: dict[str, str] | None = None,
    bw_n_bins: int = 1000,
    bw_ymax: dict[str, float] | None = None,
    label_fontsize: int = 8,
    hspace: float = 0.15,
) -> (matplotlib.figure.Figure, dict[str, matplotlib.axes.Axes])
```

Returns `(fig, axes_by_name)`. The axes dict includes an entry `'_axis'` for the coordinate ruler.

---

## Region parsing

Accepted shapes:

- `Locus(chrom, start, end)` instance
- `(chrom, start, end)` tuple/list
- `"chr6:122,600,000-122,800,000"` string (commas and spaces stripped)

---

## Track sources and detection

| Source type | Detected as |
|---|---|
| `str` ending `.bw` / `.bigwig` | `bw` |
| `str` ending `.narrowPeak` | `narrowPeak` |
| `str` ending `.bed` | `bed` |
| `str` ending `.bedpe` | `bedpe` |
| `Loci` instance | `bed` |
| `Genes` instance | `genes` |
| `list[Pair]` | `bedpe` |

Per-type defaults:

| Type | Height | Color |
|---|---|---|
| `bed` / `narrowPeak` | `0.2` | `#444444` |
| `bw` | `0.5` | `#4c78a8` |
| `bedpe` | `1.5` | `#888888` |
| `genes` | `1.5` | `#000000` |

---

## Internal drawers

Low-level renderers — you can import them if you want to add a track type:

| Name | Signature |
|---|---|
| `_parse_region(region)` | `→ (chrom, start, end)` |
| `_detect_track_type(track)` | `→ str` |
| `_draw_intervals(ax, track, chrom, start, end, color)` | — |
| `_draw_bigwig(ax, track, chrom, start, end, color, n_bins, ymax)` | — |
| `_draw_bedpe(ax, track, chrom, start, end, color, arc_points=40, max_arc_height=None)` | — |
| `_draw_genes(ax, genes, chrom, start, end, color, show_labels=True)` | — |
| `_stack_genes(visible, gap_frac, span)` | `→ (row_of_id, n_rows)` |
| `_format_ruler(ax, start, end)` | — |
