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

`Genes.make` parses a GTF or GFF3 (here GENCODE v49, protein-coding) into a
[`Genes`]({{ '/guide/genes/' | relative_url }}) object — three linked `Loci`
tables: `genes`, `transcripts` (each pointing at its gene row) and `features`
(exons, CDS and UTRs, each pointing at its transcript row):

```python
from genomeblocks import Genes

genes = Genes.make(GTF)            # GTF path; promoter_r defaults to 1000 bp
genes                              # Genes(<n> genes, <n> transcripts, <n> exons)
genes.head()                       # the genes table: chrom, start, end, strand, gene_id, gene_name, gene_type
```

GTF coordinates are 1-based; the tables are **0-based, half-open** like every
other table, so a gene annotated `1001-5000` is stored as `start=1000, end=5000`.
A gene's TSS is the 1-bp interval `[t, t+1)` with `t = start` on `+` and
`t = end - 1` on `-`; `genes.get_tss()` returns them as a `Loci` with
`gene_name`, `gene_id` and `gene` (the row) columns.

From these tables `Genes` lazily derives the interval sets it needs for
annotation — promoters (TSS ± `promoter_r`), exons, 5′/3′ UTRs, and gene
bodies — and merges each.

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
`genes.labels(loci)` gives the same answer as an integer code per row (an
index into `genomeblocks.genes.LABELS`) when you want to add it as a column:
`ARpF["annotation"] = genes.labels(ARpF)`.

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
> `nearest_genes(loci)` — a DataFrame of `uid`, `gene_name` and `distance` to the
> closest promoter window — useful when you want to name targets rather than
> just classify location. Both take any interval input (`as_loci`), so a
> pandas frame of peaks works as well as a `Loci`.

Next: **[Motif & ChIP-Atlas enrichment →]({{ '/walkthrough/enrichment/' | relative_url }})**
