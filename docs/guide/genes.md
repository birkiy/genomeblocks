---
title: Genes
parent: User Guide
layout: default
nav_order: 3
---

# Genes
{: .no_toc }

GTF/GFF and UCSC RefSeq parsing, plus a small toolkit for assigning CREs to genes (annotation classes, nearest gene, enhancer-to-gene bundling).
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Data model

```
Genes                       # dict[gene_id → Gene]
 └── Gene                   # Locus with .gene_id, .gene_name, .gene_type, .transcripts, .tss
      └── Transcript        # Locus with .transcript_id, .exons, .cds, .utr
           ├── Exon         # Locus with .exon_number
           ├── CDS          # Exon subclass
           └── UTR          # Exon subclass with .type ("5'" or "3'")
```

Everything inherits from `Locus`, so every object carries a UID, `.length`, `.center`, and participates in the interval APIs.

---

## Parsing GTF / GFF

```python
from genomeblocks import Genes

genes = Genes.make("gencode.v38.annotation.gtf",
                   gene_name_key="gene_name",   # override for non-GENCODE sources
                   gene_type_key="gene_type",
                   chr_map={"1": "chr1", "2": "chr2", ...},  # optional rename
                   promoter_r=1000)             # ±1 kb around TSS
```

The parser handles feature types `gene`, `transcript`, `exon`, `CDS`, `five_prime_UTR`, `three_prime_UTR`, and the legacy `UTR` (auto-classified 5'/3' from CDS position). Unknown feature types are collected and reported once at the end.

---

## Parsing UCSC RefSeq

If you have a UCSC `refGene.txt` / `ncbiRefSeq.txt` dump:

```python
genes = Genes.make_ucsc("ncbiRefSeq.txt",
                        chr_map=None,
                        promoter_r=1000,
                        keep_alt_contigs=False)
```

- Alt-contig transcripts (e.g. `chr6_GL000251v2_alt`) are dropped by default. Pass `keep_alt_contigs=True` to keep them under a `gene_name__chrom` key.
- Transcript-ID collisions within the same gene (e.g. MHC paralogs) get a `__N` suffix — no row is silently dropped.

---

## TSS helpers

```python
tss_by_gene       = genes.get_tss()                     # {gene_name → Locus}
tss_by_transcript = genes.get_tss_transcripts()         # {tx_id    → Locus}
tss_coding        = genes.get_tss(gene_type="protein_coding")
```

`get_tss_transcripts()` captures alternative promoters — TP53, CDKN2A, TCF7L2, etc. have multiple functional TSSs. Use `nearest_transcripts()` downstream when alt-promoter resolution matters.

---

## Annotating Loci

Four annotations wrap the common use cases:

### 1. Region class per CRE

```python
df = genes.annotations(cre)
# Columns: uid, annotation
# annotation ∈ {Promoter-TSS, 5UTR, 3UTR, Exonic, Intronic, Intergenic}
```

Priority order is Promoter-TSS → 5UTR → 3UTR → Exonic → Intronic → Intergenic.

### 2. Nearest gene

```python
df = genes.nearest_genes(cre)
# Columns: Chr, Start, End, Name (cre uid), Name_b (gene name), Distance
```

TSS is slopped by `promoter_r` before the nearest lookup, so a CRE inside a promoter window is reported as 0-distance to that gene.

### 3. Nearest transcript

```python
df = genes.nearest_transcripts(cre)   # same shape; Name_b = transcript_id
```

### 4. ROSE-style enhancer→gene bundle

```python
df = genes.enhancer_to_genes(cre, prox=50_000, level="gene")
# Columns: uid, overlap (comma-separated), proximal (within ±50 kb of TSS), closest (single name)
```

Unlike the original ROSE, `closest` uses real interval-to-interval distance (`pyranges.nearest`), not enhancer-center → TSS. Pass `level="transcript"` to resolve alt-promoters.

---

## Pre-computed annotation Loci

`genes.annot` (lazy-built dict) exposes the building blocks used by `annotations()`:

```python
genes.annot["body"]   # all gene bodies
genes.annot["prom"]   # TSS slopped by promoter_r, sorted and merged
genes.annot["exon"]   # all exons, sorted and merged
genes.annot["utr5"]   # all 5' UTRs, sorted and merged
genes.annot["utr3"]   # all 3' UTRs, sorted and merged
```

They are ordinary `Loci`, so you can intersect, tag, or use them as sources for `Architecture.make_spread(...)`.

---

## Printing

```python
print(genes.table())
# Name                 Count
# ------------------------------
# Transcripts          227463
# Exons                1398532
# CDS                  736812
# UTR                  358301
```

---

## End-to-end: annotate CREs and drop a CSV

```python
from genomeblocks import Loci, Genes

cre = Loci.make("cre.bed")
genes = Genes.make("gencode.v38.annotation.gtf")

annot = genes.annotations(cre)           # region class
near  = genes.nearest_genes(cre)         # nearest gene

out = (annot
       .merge(near[["Name", "Name_b"]].rename(columns={"Name": "uid", "Name_b": "nearest_gene"}),
              on="uid", how="left"))
out.to_csv("cre_annotated.csv", index=False)
```
