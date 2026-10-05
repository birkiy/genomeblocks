---
title: Columnar
parent: Design
layout: default
nav_order: 9
---

# Columnar
{: .no_toc }

`genomeblocks.columnar` (new in 1.1) stores CREs, genes and the contact
graph as numpy tables that share one genome and line up by row. Objects still
appear when you look at a single row, and the classic functions accept the
tables.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

```python
import genomeblocks.columnar as gbc
cre   = gbc.Loci.make("atac.narrowPeak")          # sorted, one block per chromosome
genes = gbc.Genes.make("gencode.gtf")              # GTF or GFF3
A = (gbc.Architecture.make(cre, "loops.bedpe")
       .add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes).strength())
A.chrom("chr8"); A.cis; A.trans                    # views, no copy
```

## Tables that line up by row
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-rows.svg %}
</div><figcaption>
<strong>The row number is the join key.</strong> Row <em>i</em> of the CREs is row <em>i</em> of every label, signal and graph column, and the Architecture's edge table stores row numbers, not uids. Joining a result back to its CREs is array indexing, with no dictionary in between. Genes are three tables linked the same way: transcripts point to their gene's row, features to their transcript's row.
</figcaption></figure>

- **`Genome`** is the one chromosome dictionary of a session. Tables store
  chromosomes as small integer codes into it, so codes from any two tables are
  directly comparable. Codes only grow (a new name gets the next code), and
  sorting uses a separate natural order (chr1, chr2, …, chr10, …, chrX).
- **`Loci`** holds `codes`, `starts`, `ends`, `strands` plus optional extra
  columns. `L[i]` returns a `LocusView`, a real `Locus` whose fields read and
  write the columns. `make` sorts into genome order, so each chromosome is one
  contiguous block of rows.
- **Interval kernels** place every chromosome on one axis
  (`code · 2⁴⁰ + position`), so overlap, merge and nearest are a handful of
  numpy calls over the whole genome. `overlaps_any` keeps a running maximum of
  reference ends, which works for references of any length.
- **`Genes`** parses with polars when installed and pandas otherwise (GTF and
  GFF3), and `select_isoforms` applies the classic rule vectorised.

## One edge table, sorted into blocks
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-edges.svg %}
</div><figcaption>
<strong>Views are offsets, not copies.</strong> Edges are kept sorted so each chromosome's cis edges form one contiguous block, in genome order, and every trans edge sits in a final block. <code>A.chrom("chr2")</code>, <code>A.cis</code> and <code>A.trans</code> are slices of the same arrays. Edits through a view write through to the full graph.
</figcaption></figure>

- **Neighbours** come from one adjacency index over all edges, so a CRE's trans
  partners are always included.
- **`A.graph()`** builds a graph-tool `Graph` from the arrays on demand (vertex
  *i* is row *i*) and caches it; components, PageRank and the rest run there.
- **Trans loops** are first-class: `normalize` gives them the mean trans weight
  as their expectation, and `prune` keeps them.
- **`save(path)` / `load(path)`** write parquet tables plus a small JSON file,
  readable from pandas, polars or R.
- **`to_legacy()`** returns a classic `Architecture` for the drawing helpers.

## HiChIP short-range tracks and super-enhancers
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-se.svg %}
</div><figcaption>
<strong>A ChIP track from HiChIP pairs.</strong> Pairs closer than about 1 kb are mostly undigested ChIP fragments, so their 5′ ends behave like ChIP-seq reads. <code>hichip</code> turns them into a BED for MACS3 and a 147-bp-fragment coverage bigWig. <code>se.call_se</code> stitches peaks, scores each stitched region by mean signal × width, and cuts at the same slope-1 knee as prime hubs.
</figcaption></figure>

The short-range step gives the same ends and a byte-identical bedGraph as the
usual `awk | sort | bedtools genomecov` recipe.

## Sharing results
{: .sec-purple }

- **`View`** writes a one-file HTML browser (about 40 kB of viewer plus the
  data as gzipped typed arrays). It shows CREs, loops, anchor contact profiles,
  bigWig signal, intervals and gene models, with gene search and a CRE
  inspector. It needs no server.
- **`igv_html`** writes the same data as an igv.js page for people who prefer
  IGV.

## Optional dependencies
{: .sec-navy }

`pip install "genomeblocks[columnar]"` adds polars (faster GTF/GFF3 and pairs
parsing) and pyarrow (parquet save/load). Without polars every reader falls
back to pandas with identical results.

Measured against the classic modules on the
[Benchmarks]({{ '/benchmarks/#columnar' | relative_url }}) page.
