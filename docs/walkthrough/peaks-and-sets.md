---
title: "1. Peaks, set algebra & categories"
parent: "Example: AR & FOXA1"
layout: default
nav_order: 1
---

# Peaks, set algebra & categories
{: .no_toc }

How `genomeblocks` represents genomic intervals, and how set operations turn raw
peak calls into the **accessible chromatin** background and the **AR+F / AR−F**
categories.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Loci`: the interval table

A [`Loci`]({{ '/guide/loci/' | relative_url }}) is a table of intervals: numpy
columns `chrom` (as integer codes against a `Genome`), `start`, `end`, `strand`
plus any extra columns, 0-based and half-open. The row number is the join key
for everything downstream — a signal cube, a motif matrix, an annotation table
all have one row per locus in the same order. Every row also has a stable
`uid` (`chr:start-end(strand)`), which labels results.

Load a BED / narrowPeak file with `Loci.make`:

```python
import genomeblocks as gb
from genomeblocks import Loci

peaks = {k: Loci.make(v) for k, v in PEAK.items()}   # PEAK = {name: path}
for k, v in peaks.items():
    print(f"{k:12} {len(v):>8,} peaks")
```

`Loci.make` reads `chrom / start / end` (and the strand column when present),
so ChIP-Atlas `.05.bed` files load directly; `keep=True` also keeps `name`,
`score` and the narrowPeak columns. Rows come back sorted in genome order.

```python
peaks["AR_4h"]                 # Loci(n=3,645, chroms=..., sorted)
peaks["AR_4h"].head(3)         # the first rows, still a Loci
peaks["AR_4h"].to_pandas()     # chrom, start, end, strand as a DataFrame
```

## Set algebra — and one subtlety that matters

`Loci` overloads Python's set operators. The distinction between them is the
single most important idea in this page:

| op | method | result |
|---|---|---|
| `a & b` | `intersect` | the **`a` rows** that overlap any `b` row (left-hand intervals kept) |
| `a - b` (and `/`) | `difference` | the **`a` rows** that overlap **no** `b` row |
| `a + b` (and `\|`) | concatenation | every row from both — **duplicates kept, nothing merged** |
| `a.merge()` | — | fuse overlapping or book-ended intervals into one, in genome order |

{: .warning }
> `+` (and `|`) **concatenate** — they do *not* deduplicate or merge. To build a
> true *union of regions* follow with `.merge()`. And `a & b` returns the
> original `a` rows that overlap `b` (peak-level membership), **not** the
> geometric intersection rectangle. This is what you want for **"which of my
> peaks fall in this other set."**

The right-hand side of any of these is anything
[`as_loci`]({{ '/interoperability/' | relative_url }}) takes — another `Loci`, a
BED path, a pandas / polars frame — so `peaks["AR_4h"] & "atac.narrowPeak"`
works as written.

## Accessible chromatin = the ATAC union

Open chromatin (ATAC-seq peaks) is our universe of "places a factor *could*
bind." We pool all four ATAC peak sets (both timepoints, both replicates) and
merge them into one non-redundant region set:

```python
accessible = (peaks["ATAC_0h_r1"] + peaks["ATAC_0h_r2"]
              + peaks["ATAC_4h_r1"] + peaks["ATAC_4h_r2"]).merge()
```

`+` stacks the hundreds of thousands of peaks; `.merge()` collapses overlaps
(the same enhancer called in two replicates) into a single interval and returns
the result sorted. `accessible` is the background pool we reuse for motif
enrichment later.

## Inside vs outside accessible chromatin

A quick sanity check on the biology: transcription factors mostly bind open
chromatin. We count, for each factor/timepoint, how many of its peaks overlap
`accessible`:

```python
for name in ["AR_0h", "AR_4h", "FOXA1_0h", "FOXA1_4h"]:
    s = peaks[name]
    inside = len(s & accessible)        # peaks overlapping accessible
    total  = len(s)
    print(f"{name:10} {100*inside/total:5.1f}% inside accessible "
          f"({inside:,}/{total:,})")
```

`len(s & accessible)` is the count of `s` peaks that land in open chromatin
(`s.overlap_any(accessible).sum()` gives the same number as a boolean mask).
Plotting `inside` vs `total - inside` as a stacked bar shows the expected
picture: the large majority of AR and FOXA1 binding sits in ATAC-accessible
regions.

## Keep only accessible peaks

We restrict the three peak sets that define our categories to accessible
chromatin, so every downstream comparison is on equal (open-chromatin) footing:

```python
F0 = peaks["FOXA1_0h"] & accessible
F4 = peaks["FOXA1_4h"] & accessible
A4 = peaks["AR_4h"]    & accessible
```

## A Venn over genomic intervals

Three peak sets do not share identical intervals, so a Venn needs a common
*universe*: merge all peaks into one region set, then ask which input sets
touch each region. `overlap_any` returns one boolean per universe row, and the
row numbers where it is `True` are the region ids `matplotlib_venn` counts:

```python
import numpy as np
from matplotlib_venn import venn3

universe = (F0 + F4 + A4).merge()
ids = [set(np.flatnonzero(universe.overlap_any(s)).tolist()) for s in (F0, F4, A4)]
venn3(ids, set_labels=["FOXA1 0h", "FOXA1 4h", "AR 4h"])
```

## Defining AR+F and AR−F

The Venn motivates the split. "FOXA1-bound" means a FOXA1 peak at **either**
timepoint, so we union `F0` and `F4` first, then partition the 4 h AR peaks:

```python
F_any = (F0 + F4).merge()          # FOXA1-bound at 0 h or 4 h
ARpF  = A4 & F_any                 # AR 4h that IS FOXA1-bound  -> FOXA1-dependent
ARmF  = A4 - F_any                 # AR 4h that is NOT FOXA1-bound -> FOXA1-independent
```

`ARpF` and `ARmF` are disjoint and together equal `A4`. They are the two
`Loci` tables carried through the rest of the walkthrough.

{: .tip }
> Because `&` and `-` keep the **AR** rows, `ARpF` and `ARmF` are genuine AR
> peak tables you can feed straight into signal extraction, annotation, motif
> scanning, and the browser — no coordinate bookkeeping required. Overlap
> queries run on genomeblocks' own numpy kernel by default; `backend=` on
> `intersect` / `difference` / `merge` switches to bioframe, pyranges, ncls,
> cgranges or bedtools (see [Backends]({{ '/backends/' | relative_url }})).

Next: **[Signal heatmaps →]({{ '/walkthrough/signal-heatmaps/' | relative_url }})**
