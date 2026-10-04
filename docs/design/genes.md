---
title: Genes
parent: Design
layout: default
nav_order: 2
---

# Genes
{: .no_toc }

`Genes` parses a GTF (or a UCSC table) into a dictionary of gene objects, then
builds a cached index of promoter, UTR, exon and gene-body intervals that
labels any set of regions.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## From GTF lines to objects
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/genes-model.svg %}
</div><figcaption>
<strong>One pass over the file builds a tree.</strong> <code>Genes</code> is a <code>dict</code> keyed by <code>gene_id</code>. Each <code>Gene</code> holds its <code>Transcript</code>s, and each transcript holds its <code>Exon</code>, <code>CDS</code> and <code>UTR</code> records. Every node is a <code>Locus</code>, so any of them can go straight into <code>Loci</code> operations.
</figcaption></figure>

- Attributes are read from `gene_name_key` / `gene_type_key` (GENCODE defaults),
  and `chr_map` renames chromosomes while parsing.
- `Gene.tss` and `Transcript.tss` are 1-bp loci on the strand-aware start.
- `Genes.make_ucsc` reads UCSC `refGene`-style tables into the same objects.

## Labelling regions: the annotation index
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/genes-annot.svg %}
</div><figcaption>
<strong>First hit wins.</strong> <code>genes.annot</code> holds five merged <code>Loci</code>: promoters (every TSS ± <code>promoter_r</code>), 5′ UTRs, 3′ UTRs, exons and gene bodies. <code>annotations(cres)</code> tests each CRE against them in that order and labels it with the first set it overlaps. A CRE that overlaps none of them is <em>Intergenic</em>.
</figcaption></figure>

The index is built on first access and cached on the object, so labelling
several CRE sets costs one build. `nearest_genes(cres)` is a separate path: a
`nearest` join (pyranges) against every TSS widened by `promoter_r`, returning
the closest gene name per CRE.

## Picking the isoform your cells use
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/genes-isoforms.svg %}
</div><figcaption>
<strong>Open chromatin decides the TSS.</strong> An isoform is <em>supported</em> when its TSS window (± <code>r</code>) overlaps a peak. When bigWigs are also given, its TSS signal must exceed <code>min_signal</code> and reach <code>min_frac</code> of the gene's best TSS. The longest supported isoform becomes canonical (<code>rank="signal"</code> picks the strongest TSS instead). With <code>collapse=True</code> the gene's span and TSS move to it, which every downstream label and nearest-gene call then reads.
</figcaption></figure>

Genes with no supported isoform fall back to the same ranking over all their
isoforms, so every gene keeps a canonical transcript, and no transcript is
removed. The per-isoform evidence stays on the objects
(`Transcript.tss_score`, `Transcript.tss_support`, `Gene.canonical`).

## Costs
{: .sec-navy }

- **Parsing** is a pure-Python line loop that creates one object per GTF
  record, so it is the slowest step of the classic pipeline. The
  [columnar Genes]({{ '/design/columnar/' | relative_url }}) parses with polars (or pandas) into three
  tables instead.
- **Labelling** is five indexed lookups per CRE against merged interval sets.
- **Isoform selection** is one indexed lookup per transcript TSS, plus one
  bigWig summary per TSS when `bw` is given.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/#genes' | relative_url }}) page.
