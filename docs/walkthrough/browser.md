---
title: "5. Genome browser"
parent: "Example: AR & FOXA1"
layout: default
nav_order: 5
---

# Genome browser
{: .no_toc }

An IGV-like view of every condition at one locus — signal tracks (with replicate
averaging and shared y-axes), peak calls, and gene models.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What the browser is for

Heatmaps and enrichment summarise thousands of regions; a **browser view** zooms
into one locus to *see* the data — the per-base coverage, where peaks were
called, and which genes are nearby. It is how you sanity-check a result and how
you build a figure for a specific gene.

## Tracks are a dict, types are auto-detected

`browser(region, tracks)` takes a region and an **ordered dict** of named
tracks. Each value's type is inferred from its extension or Python type:

| track value | rendered as |
|---|---|
| `*.bw` path (or a **list** of them) | binned coverage (a list is **averaged**) |
| `*.bed` / `*.narrowPeak` path, or a `Loci` | interval rectangles |
| a `Genes` object | stacked gene models (exons / CDS) |
| `*.bedpe` path / `list[Pair]` | arc track |

```python
from genomeblocks.browser import browser

region = "chr19:50,792,009-50,923,669"
tracks = {
    "ATAC 0h":   [BW["ATAC_0h_r1"], BW["ATAC_0h_r2"]],   # list -> averaged
    "ATAC 4h":   [BW["ATAC_4h_r1"], BW["ATAC_4h_r2"]],
    "AR 0h":     BW["AR_0h"],
    "AR 4h":     BW["AR_4h"],
    "FOXA1 0h":  BW["FOXA1_0h"],
    "FOXA1 4h":  BW["FOXA1_4h"],
    "AR 4h peaks":    PEAK["AR_4h"],
    "FOXA1 4h peaks": PEAK["FOXA1_4h"],
    "genes":     genes,
}
```

The region accepts a `"chr:start-end"` string (commas allowed), a
`(chrom, start, end)` tuple, or a `Locus`.

## Replicate averaging, the browser way

The two ATAC replicates are passed as a **list of bigWig paths** — the browser
extracts each and averages their per-bin means into a single track. This is the
same operation the [heatmap](signal-heatmaps) does over its track columns, so the
two figures show the ATAC replicates consistently.

## Shared y-axes for honest comparison

By default each bigWig track auto-scales to its own maximum — which makes a
low-signal 0 h track look as tall as a high-signal 4 h track and hides the
induction. `bw_share` groups tracks that should share one y-scale (the group's
region maximum):

```python
fig, axes = browser(region, tracks, bw_n_bins=2000, figsize=(11, None),
                    bw_share=[["ATAC 0h", "ATAC 4h"],
                              ["AR 0h", "AR 4h"],
                              ["FOXA1 0h", "FOXA1 4h"]])
```

Now AR 0 h and AR 4 h sit on the same axis, so the DHT-induced AR gain — and the
ATAC and FOXA1 changes — are read directly off the heights.

| control | effect |
|---|---|
| `bw_n_bins` | bins per bigWig track (horizontal resolution) |
| `bw_share` | list of name-groups that share a y-scale |
| `bw_ymax` | a fixed y-max: a scalar for all bigWig tracks, or a per-track dict |
| `figsize=(w, None)` | width fixed; height derived from the track heights |

`browser` returns `(fig, axes_by_name)`, so you can grab any track's axis by name
to annotate it (highlight a peak, mark a TSS) before saving.

## Reading the result

At this prostate locus (the KLK locus on chr19) you can see the logic of the
whole analysis in one panel: ATAC marks the accessible landscape, FOXA1 occupies
sites at 0 h and 4 h, and AR appears/strengthens at 4 h — strongest where FOXA1
is already bound (the AR+F sites) — with the called peaks and gene models lined
up underneath.

---

That completes the walkthrough. The full, runnable notebook is at
[`examples/ar_foxa1_lncap/`](https://github.com/birkiy/genomeblocks/tree/main/examples/ar_foxa1_lncap);
for per-module reference see the [User Guide](../guide/).
