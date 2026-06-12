---
title: genes
parent: API Reference
layout: default
nav_order: 4
---

# `genomeblocks.genes`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Transcript(Locus)`

```python
@dataclass
class Transcript(Locus):
    transcript_id: str = ""
    exons: Loci
    cds:   Loci
    utr:   Loci
```

Methods: `add_exon(e)`, `add_cds(c)`, `add_utr(u)`.

---

## `Gene(Locus)`

```python
@dataclass
class Gene(Locus):
    gene_id: str = ""
    gene_name: str = ""
    gene_type: str | None = ""
    transcripts: dict[str, Transcript]
    # computed in __post_init__:
    tss: Locus   # TSS as a 1-bp Locus (strand-aware)
```

Method: `add_transript(t_id, t)`.

---

## `Genes(dict)`

Subclass of `dict[str, Gene]`.

### Factories

```python
Genes.make(filename, gene_name_key="gene_name",
                     gene_type_key="gene_type",
                     chr_map=None,
                     promoter_r=1000) -> Genes
# GTF/GFF parser. Handles gene/transcript/exon/CDS/five_prime_UTR/three_prime_UTR/UTR.

Genes.make_ucsc(filename, chr_map=None,
                          promoter_r=1000,
                          keep_alt_contigs=False) -> Genes
# UCSC RefSeq table parser (column order: bin, name, chrom, strand, txStart, ...).
```

### Properties

| Name | Description |
|---|---|
| `filename` | Source file path. |
| `_promoter_r` | Promoter half-window (bp). |
| `annot` | Lazy `{'body','prom','exon','utr5','utr3'}` dict of merged/sorted `Loci`. |

### Methods

```python
Genes.get_tss(gene_type=None) -> dict[str, Locus]
# {gene_name → TSS Locus}; optional gene_type filter.

Genes.annotations(loci) -> pandas.DataFrame
# Columns: uid, annotation ∈ {Promoter-TSS, 5UTR, 3UTR, Exonic, Intronic, Intergenic}

Genes.nearest_genes(loci) -> pandas.DataFrame
# pyranges.nearest against slopped-TSS Loci; Name_b = gene name, Distance in bp.

Genes.table() -> str
# Summary of counts.
```
