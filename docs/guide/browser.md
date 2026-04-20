---
title: Browser
parent: User Guide
layout: default
nav_order: 6
---

# Browser
{: .no_toc }

An IGV-like, fully-vectorial, matplotlib-native region viewer. One call dispatches across bigWig / narrowPeak / BED / BEDPE / Genes and renders publication-ready SVG or PDF.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## One-call API

```python
from genomeblocks import browser, Loci, Genes
from genomeblocks.bedpe import read_bedpe

fig, axes = browser(
    region=("chr6", 122_600_000, 122_800_000),
    tracks={
        "HiChIP loops": read_bedpe("loops.bedpe"),
        "ATAC signal":  "ATAC.bw",
        "ATAC peaks":   Loci.make("peaks.narrowPeak"),
        "Genes":        Genes.make("gencode.v38.gtf"),
    },
    figsize=(12, None),      # height auto-sized to track count
    bw_n_bins=120,
    colors={"HiChIP loops": "#DA0000",
            "ATAC signal":  "#4c78a8"},
    bw_ymax={"ATAC signal": 150.0},
)
fig.savefig("region.svg")
```

Track *type* is auto-detected from:

| Source | Rendered as |
|---|---|
| `*.bw`, `*.bigwig` | binned coverage (fill_between + outline) |
| `*.narrowPeak`, `*.bed` / `Loci` | rectangle rows |
| `*.bedpe` / `list[Pair]` | half-sine arcs between anchors |
| `Genes` | stacked gene models (thin body line, exon boxes, taller CDS boxes, optional label & strand arrow) |

---

## Region specification

Any of:

```python
region = "chr6:122,600,000-122,800,000"
region = ("chr6", 122_600_000, 122_800_000)
region = some_locus                   # Locus → .chrom/.start/.end
```

Commas and spaces are stripped from string forms.

---

## Styling knobs

- `track_heights={name: float}` — per-track height in inches (defaults: bw=0.5, peaks=0.2, bedpe=1.5, genes=1.5).
- `colors={name: "#hex"}` — per-track color override.
- `bw_n_bins=N` — bin count per bigWig track (larger → finer detail, slower reads).
- `bw_ymax={name: y}` — fix a per-track y-axis; otherwise auto-scaled to the visible max.
- `label_fontsize`, `hspace` — cosmetic fine-tuning.

Returns `(fig, axes_by_name)` — you can grab any axis by name and overlay annotations (highlight boxes, TSS arrows, etc.) before saving.

---

## Why SVG-clean?

The browser avoids `imshow` / rasterized patches entirely: coverage tracks use `fill_between` + `steps-mid` outlines, intervals use `Rectangle`/`PatchCollection`, and loops use vector `plot()` arcs. The output is small, editable in Inkscape / Illustrator, and reads well in Jupyter at any zoom level.

---

## Real example

See [`examples/browser_example.py`](https://github.com/birkiy/genomeblocks/blob/main/examples/browser_example.py) for a rendering of the Nanog locus (mm10, chr6:122,286,666-122,902,344) with HiChIP loops, ATAC signal, ATAC peaks, and GENCODE protein-coding genes — the output is [`examples/browser_Nanog.svg`](https://github.com/birkiy/genomeblocks/blob/main/examples/browser_Nanog.svg).

---

## Tips

- For very wide regions, drop `bw_n_bins` to 120–240 to keep vector size tiny without losing the profile shape.
- `max_arc_height` (passed through `**kwargs` to `_draw_bedpe`) caps arc height; defaults to half the view span. Loops wider than the region clip at the top, preserving the takeoff angle.
- Gene stacking uses a greedy per-row algorithm with a `span × 0.1` gap; pass `show_labels=False` to hide names in dense regions.
