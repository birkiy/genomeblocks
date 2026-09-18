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
    tss_score:   float | None = None   # set by Genes.select_isoforms()
    tss_support: bool  | None = None   # set by Genes.select_isoforms()
    # property:
    tss: Locus   # TSS as a 1-bp Locus (strand-aware)
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
    canonical: str | None = None   # transcript key set by Genes.select_isoforms()
    # computed in __post_init__:
    tss: Locus   # TSS as a 1-bp Locus (strand-aware)
    # property:
    canonical_transcript: Transcript | None
```

Methods: `add_transript(t_id, t)`, `set_span(start, end)` (moves the body and refreshes `tss`).

---

## `Genes(dict)`

Subclass of `dict[str, Gene]`.

### Factories

```python
Genes.make(filename, gene_name_key="gene_name",
                     gene_type_key="gene_type",
                     chr_map=None,
                     promoter_r=1000,
                     cre=None, bw=None, r=None, kw=None) -> Genes
# GTF/GFF parser. Handles gene/transcript/exon/CDS/five_prime_UTR/three_prime_UTR/UTR.

Genes.make_ucsc(filename, chr_map=None,
                          promoter_r=1000,
                          keep_alt_contigs=False,
                          cre=None, bw=None, r=None, kw=None) -> Genes
# UCSC RefSeq table parser (column order: bin, name, chrom, strand, txStart, ...).
```

`cre` (BED/narrowPeak path, `Loci`, `Locus`, or a list of any of those) and/or `bw`
(bigwig path or list of paths) run `select_isoforms()` right after parsing, so each
gene follows its longest **ATAC-supported** isoform instead of its longest annotated
one. `r` is the TSS half-window and defaults to `promoter_r`; `kw` is a dict
forwarded to `select_isoforms()` (`agg`, `min_signal`, `min_frac`, `rank`,
`collapse`, `verbose`). Peaks/bigwigs must use the same chromosome names as the
parsed genes (i.e. post `chr_map`).

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

Genes.select_isoforms(cre=None, bw=None, *, r=None, agg="max",
                      min_signal=0.0, min_frac=0.5, rank="longest",
                      collapse=True, verbose=True) -> Genes
# Keep the isoforms whose TSS is open, and point each gene at the best of them.
```

`select_isoforms()` marks every transcript (`tss_score`, `tss_support`), sets
`Gene.canonical`, and — with `collapse=True` — moves each gene's body/TSS onto its
canonical isoform (so `annot`, `get_tss()`, `annotations()` and `nearest_genes()`
follow it). No transcript is dropped. A transcript is supported when it passes each
criterion given:

| Argument | Effect |
|---|---|
| `cre` | TSS window (±`r`, default `promoter_r`) must overlap an ATAC/DNase peak. A `Loci` is used as-is, so its interval index is reused. |
| `bw` | TSS window must score `> min_signal` **and** reach `min_frac` of the best TSS score among that gene's candidates. Several bigwigs → a TSS keeps its highest score (open in *any* sample). |
| `rank` | `'longest'` (default) picks the longest supported isoform, `'signal'` the strongest TSS; the other value breaks ties. |

Genes with no supported isoform fall back to the same ranking over all their
isoforms, so every gene still gets a canonical transcript.
