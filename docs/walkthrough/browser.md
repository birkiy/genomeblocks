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

`gb.browser(region, tracks)` takes a region and an **ordered dict** of named
tracks. Each value's type is inferred from its extension or Python type:

| track value | rendered as |
|---|---|
| `*.bw` path or an open bigWig handle (or a **list** of them) | binned coverage (a list is **averaged**) |
| any interval input — a `Loci`, a `*.bed` / `*.narrowPeak` path, a frame, a list of regions | interval rectangles |
| a `Genes` object | stacked gene models (exons / CDS) |
| a `Pairs` or a `*.bedpe` path | arc track |
| a `*.bam` path (with `reference=`) | per-base coverage with mismatch colouring |

```python
import genomeblocks as gb

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
    "AR+F":      ARpF,                                     # a Loci table is a track too
    "genes":     genes,
}
```

The region accepts a `"chr:start-end"` string (commas and `kb` / `Mb` units
allowed), a `(chrom, start, end)` tuple, or a `Locus`. Like every coordinate in
genomeblocks it is 0-based and half-open: `chr19:50,792,009-50,923,669` is the
131,660 bases from position 50,792,009 up to but not including 50,923,669.

## Replicate averaging, the browser way

The two ATAC replicates are passed as a **list of bigWig paths** — the browser
extracts each and averages their per-bin means into a single track. This is the
same operation the [heatmap]({{ '/walkthrough/signal-heatmaps/' | relative_url }}) does over its track columns, so the
two figures show the ATAC replicates consistently.

## Shared y-axes for honest comparison

By default each bigWig track auto-scales to its own maximum — which makes a
low-signal 0 h track look as tall as a high-signal 4 h track and hides the
induction. `bw_share` groups tracks that should share one y-scale (the group's
region maximum):

```python
fig, axes = gb.browser(region, tracks, bw_n_bins=2000, figsize=(11, None),
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
| `backend` | the bigWig reader (`pybigtools`, `pybigwig`, `python`) |

`browser` returns `(fig, axes_by_name)`, so you can grab any track's axis by name
to annotate it (highlight a peak, mark a TSS) before saving; `axes["_axis"]` is
the coordinate ruler.

## Reading the result

At this prostate locus (the KLK locus on chr19) you can see the logic of the
whole analysis in one panel: ATAC marks the accessible landscape, FOXA1 occupies
sites at 0 h and 4 h, and AR appears/strengthens at 4 h — strongest where FOXA1
is already bound (the AR+F sites) — with the called peaks and gene models lined
up underneath.

{: .tip }
> For an interactive version of the same view, `gb.igv_html(path, regions=..., loci=..., genes=..., signal=...)`
> writes a one-file IGV page and `gb.View(...)` a one-file interactive view;
> see the [Browser guide]({{ '/guide/browser/' | relative_url }}).

---

That completes the walkthrough. The full, runnable notebook is at
[`examples/ar_foxa1_lncap/`](https://github.com/birkiy/genomeblocks/tree/main/examples/ar_foxa1_lncap);
for per-module reference see the [User Guide]({{ '/guide/' | relative_url }}).
