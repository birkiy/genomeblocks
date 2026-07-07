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

## `Loci`: the interval container

A [`Loci`](../guide/loci) is a list of `Locus` intervals (`chrom, start, end,
strand`) with a fast overlap index ([cgranges](https://github.com/lh3/cgranges))
built lazily behind the `.cgr` property. Every locus has a stable `.uid`
(`chr:start-end(strand)`) used as its identity in set operations and joins.

Load a BED/narrowPeak file with `Loci.make`:

```python
from genomeblocks import Loci

peaks = {k: Loci.make(v) for k, v in PEAK.items()}   # PEAK = {name: path}
for k, v in peaks.items():
    print(f"{k:12} {len(v):>8,} peaks")
```

`Loci.make` reads the first three columns as `chrom/start/end` (and column 6 as
strand if present), so ChIP-Atlas `.05.bed` files load directly.

## Set algebra — and one subtlety that matters

`Loci` overloads Python's set operators. The distinction between them is the
single most important idea in this page:

| op | method | result |
|---|---|---|
| `a & b` | `intersect` | the **`a` peaks** that overlap any `b` peak (LHS intervals kept) |
| `a - b` (and `\`) | `difference` | the **`a` peaks** that overlap **no** `b` peak |
| `a + b` (and `\|`) | concatenation | every locus from both — **duplicates kept, nothing merged** |
| `a.sort().merge()` | — | sort, then fuse overlapping/adjacent intervals into one |

{: .warning }
> `+` (and `|`) **concatenate** — they do *not* deduplicate or merge. To build a
> true *union of regions* you must follow with `.sort().merge()`. And `a & b`
> returns the original `a` intervals that overlap `b` (peak-level membership),
> **not** the geometric intersection rectangle. This is what you want for **"which
> of my peaks fall in this other set."**

## Accessible chromatin = the ATAC union

Open chromatin (ATAC-seq peaks) is our universe of "places a factor *could*
bind." We pool all four ATAC peak sets (both timepoints, both replicates) and
merge them into one non-redundant region set:

```python
accessible = ( \
    peaks["ATAC_0h_r1"] \
    + peaks["ATAC_0h_r2"] \
    + peaks["ATAC_4h_r1"] \
    + peaks["ATAC_4h_r2"]
).sort().merge()
```

`+` stacks the ~hundreds-of-thousands of peaks; `.sort().merge()` collapses
overlaps (e.g. the same enhancer called in two replicates) into a single
interval. The result, `accessible`, is the background pool we reuse for motif
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

`len(s & accessible)` is the count of `s` peaks that land in open chromatin.
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

A Venn of three peak sets is a plain set intersection.

```python
from matplotlib_venn import venn3

venn3((F0, F4, A4), set_labels=["FOXA1 0h", "FOXA1 4h", "AR 4h"])
```

## Defining AR+F and AR−F

The Venn motivates the split. "FOXA1-bound" means a FOXA1 peak at **either**
timepoint, so we union `F0` and `F4` first, then partition the 4 h AR peaks:

```python
F_any = (F0 + F4).sort().merge()   # FOXA1-bound at 0 h or 4 h
ARpF  = A4 & F_any                 # AR 4h that IS FOXA1-bound  -> FOXA1-dependent
ARmF  = A4 - F_any                 # AR 4h that is NOT FOXA1-bound -> FOXA1-independent
```

`ARpF` and `ARmF` are disjoint and together equal `A4`. They are the two
`Loci` sets carried through the rest of the walkthrough.

{: .tip }
> Because `&` and `-` keep the **AR** intervals, `ARpF` and `ARmF` are genuine AR
> peak sets you can feed straight into signal extraction, annotation, motif
> scanning, and the browser — no coordinate bookkeeping required.

Next: **[Signal heatmaps →](signal-heatmaps)**
