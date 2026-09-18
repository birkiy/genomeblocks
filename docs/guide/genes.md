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
 └── Gene                   # Locus with .gene_id, .gene_name, .gene_type, .transcripts, .tss, .canonical
      └── Transcript        # Locus with .transcript_id, .exons, .cds, .utr, .tss
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

## Picking the isoform your cells actually use

A GTF gene spans the **union** of its isoforms, so `gene.start` / `gene.end` / `gene.tss`
follow the longest *annotated* transcript. For genes with a long, rarely used isoform
(TGFBR3 is the classic case) that TSS can sit tens of kilobases away from the promoter
that is open in your cell type — which then skews promoter annotation, nearest-gene
calls and the browser view.

Hand `Genes.make()` an ATAC-seq peak file and/or a bigwig and it keeps the isoforms
whose TSS is actually accessible, then points the gene at the longest of those:

```python
genes = Genes.make("gencode.v38.annotation.gtf",
                   promoter_r=1000,
                   cre="Th17_atac_peaks.narrowPeak",   # path, Loci, Locus, or a list
                   bw="Th17_atac.bw")                  # bigwig, or a list of them
```

`cre` takes a `Loci` (or a `Locus`, or a mixed list of paths and `Loci`) just as
happily as a file path — handy when the peaks are already in memory, reused from
another step, or filtered first:

```python
peaks = Loci.make("Th17_atac_peaks.narrowPeak")
genes = Genes.make("gencode.v38.annotation.gtf", cre=peaks - blacklist)
```

The TSS half-window is `r`, which defaults to `promoter_r` — one window for both
promoter annotation and isoform support, unless you override it. Either evidence
argument works on its own — peaks alone, signal alone, or both (peaks gate first,
then the signal cut, which also makes the bigwig pass much cheaper). Tuning knobs go
in `kw`, or call the method directly on an already-parsed object:

```python
genes.select_isoforms(peaks, ["naive.bw", "th17.bw", "treg.bw"],
                      r=500,            # default: the object's promoter_r
                      min_signal=0.0,   # absolute floor a TSS score must exceed
                      min_frac=0.5,     # ...and ≥50% of the gene's best TSS score
                      rank="longest",   # or "signal" for the strongest TSS
                      collapse=True)    # move gene body/TSS onto the chosen isoform
```

With several bigwigs a TSS keeps its **highest** score, so an isoform open in any one
of your samples counts as supported.

What you get back:

```python
g = genes["ENSG00000069702"]        # TGFBR3
g.canonical                         # 'ENST00000212355' — the chosen transcript key
g.canonical_transcript.tss_score    # ATAC signal at its TSS window
g.start, g.end, g.tss               # now follow that isoform (collapse=True)
[t.tss_support for t in g.transcripts.values()]   # per-isoform support calls
```

- **Nothing is dropped.** All isoforms stay in `gene.transcripts`; only the gene's own
  span, its TSS and `canonical` change.
- Genes with no supported isoform **fall back** to the longest annotated one, so every
  gene keeps a canonical transcript.
- Because the gene span moved, `annot`, `get_tss()`, `annotations()` and
  `nearest_genes()` all follow the supported isoform. `annot['exon'/'utr5'/'utr3']`
  still pool every isoform.
- Peaks and bigwigs must use the **same chromosome names** as the annotation (i.e.
  after `chr_map`); a `[WARN]` is printed if nothing at all comes out supported.
- In the browser, `browser(..., genes_max_transcripts=1)` draws the canonical isoform
  rather than the longest one.

---

## TSS helpers

```python
tss_by_gene = genes.get_tss()                          # {gene_name → Locus}
tss_coding  = genes.get_tss(gene_type="protein_coding")
```

---

## Annotating Loci

Two annotations wrap the common use cases:

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

{: .note }
> `annotations()` + `nearest_genes()` are exactly what `Architecture.annotate()`
> uses under the hood to tie each CRE to a region class and a gene.

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

They are ordinary `Loci`, so you can intersect them with CRE sets or use them as promoter sources when building an `Architecture`.

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
