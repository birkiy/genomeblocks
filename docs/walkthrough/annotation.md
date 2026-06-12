---
title: "3. Genomic annotation"
parent: "Example: AR & FOXA1"
layout: default
nav_order: 3
---

# Genomic annotation
{: .no_toc }

Labelling each AR set by its gene context — promoter, UTR, exon, intron,
intergenic — with a `Genes` model parsed from a GTF.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Why annotate

Where a regulatory element sits relative to genes is a first clue to how it
acts. Distal enhancers cluster in introns and intergenic space; promoter-proximal
sites suggest direct core-promoter regulation. Comparing the **annotation
breakdown** of AR+F vs AR−F asks whether the two cistromes occupy different parts
of the genome.

## Building a gene model

`Genes.make` parses a GTF (here GENCODE v49, protein-coding) into a `Genes`
object — a dictionary of genes, each with transcripts, exons, CDS and UTRs:

```python
from genomeblocks import Genes

genes = Genes.make(GTF)            # GTF path; promoter_r defaults to 1000 bp
```

From these features `Genes` lazily derives the interval sets it needs for
annotation — promoters (TSS ± `promoter_r`), exons, 5′/3′ UTRs, and gene bodies —
and merges each.

## Annotating a `Loci` set

`genes.annotations(loci)` returns a `pandas.DataFrame` with one row per input
locus (`uid`) and its **region class**:

```python
annot_pf = genes.annotations(ARpF)
annot_mf = genes.annotations(ARmF)
```

Each locus is assigned the **first** matching class in a fixed priority order, so
every locus gets exactly one label:

| class | meaning |
|---|---|
| `Promoter-TSS` | overlaps a TSS ± `promoter_r` window |
| `5UTR` / `3UTR` | overlaps a 5′ / 3′ UTR |
| `Exonic` | overlaps an exon |
| `Intronic` | inside a gene body but not the above |
| `Intergenic` | none of the above |

The priority order means a peak that touches both a promoter and an intron is
called `Promoter-TSS` — the most specific, regulatorily-meaningful class wins.

## Pie charts per set

A `value_counts()` on the `annotation` column gives the composition of each set,
drawn as side-by-side pies:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(9, 4))
for ax, df, title in [(axes[0], annot_pf, "AR+F"), (axes[1], annot_mf, "AR-F")]:
    vc = df["annotation"].value_counts()
    ax.pie(vc.values, labels=vc.index, autopct="%1.0f%%", textprops={"fontsize": 7})
    ax.set_title(f"{title}  (n={len(df):,})")
```

A typical enhancer-dominated cistrome is mostly `Intronic` + `Intergenic`; a
shift in the promoter fraction between AR+F and AR−F is the kind of difference
this view surfaces.

{: .note }
> `annotations()` is the lightweight, region-class labeller. `genes` also exposes
> `nearest_genes(loci)` to attach the closest gene to each peak — useful when you
> want to name targets rather than just classify location.

Next: **[Motif & ChIP-Atlas enrichment →](enrichment)**
