---
title: Tables
parent: Design
layout: default
nav_order: 9
permalink: /design/tables/
---

# Tables
{: .no_toc }

genomeblocks stores CREs, genes, loops, the contact graph and the enrichment
index as numpy tables that share one `Genome` and line up by row. A single
row appears as a `Locus` when you look at it; whole-set work runs on the
columns through the backends.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

```python
import genomeblocks as gb
cre   = gb.Loci.make("atac.narrowPeak", keep=True)   # sorted, one block per chromosome
genes = gb.Genes.make("gencode.gtf")                 # GTF or GFF3, 0-based tables
A = (gb.Architecture.make(cre, "loops.bedpe")
       .add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes).strength())
A.chrom("chr8"); A.cis; A.trans                      # views, no copy
A.save("arch/"); cre.save("cre.parquet")             # parquet, readable from polars or R
```

## Tables that line up by row
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-rows.svg %}
</div><figcaption>
<strong>The row number is the join key.</strong> Row <em>i</em> of the CREs is row <em>i</em> of every label, signal and graph column, and the Architecture's edge table stores row numbers, not uids. Joining a result back to its CREs is array indexing, with no dictionary in between. Genes are three tables linked the same way: transcripts point to their gene's row, features to their transcript's row.
</figcaption></figure>

### Columns and codes

A `Loci` is four aligned arrays plus a dict of extra columns:

| Column | dtype | Meaning |
|---|---|---|
| `codes` | int32 | chromosome code into the table's `Genome` |
| `starts` | int64 | 0-based start (BED convention) |
| `ends` | int64 | end, exclusive |
| `strands` | int8 | 0 `.`, 1 `+`, 2 `-` |
| `cols` | dict of arrays | `name`, `score`, `gene_name`, ... one value per row |

```python
cre = gb.Loci.make("peaks.bed", keep=True)
cre                       # -> Loci(n=7, chroms=2, sorted, cols=[name, score])
cre.codes                 # -> array([0, 0, 0, 0, 0, 1, 1], dtype=int32)
cre.strands               # -> array([1, 2, 1, 0, 2, 1, 1], dtype=int8)
cre.genome                # -> Genome(2 chroms: chr1, chr2)
cre["start"][:3]          # -> array([ 900, 1900, 4900])    a column
cre[0]                    # -> Locus[0](chr1:900-1100(+))   a row
cre[0].name, cre[0].score # -> ('p1', 10.0)                 extra columns on the row too
```

`cre[i]` is a `LocusView`: a real `Locus` whose `chrom`, `start`, `end` and
`strand` read and write the columns, so nothing is copied to look at one row
and an edit through it lands in the table. Derived columns (`chroms`,
`centers`, `lengths`, `strand`, `uid`, `names`) are computed from the arrays
on request; `uid` and the `uid → row` dict are built once and dropped when
the rows change.

### Sorted rows are chromosome blocks

`Loci.make` sorts into genome order, so each chromosome is one contiguous
block of rows. `chrom_offsets` gives the blocks and `by_chrom` returns a
slice of the same arrays:

```python
cre.is_sorted             # -> True
cre.chrom_offsets         # -> {'chr1': (0, 5), 'chr2': (5, 7)}
v = cre.by_chrom("chr2")  # a view on sorted Loci
import numpy as np
np.shares_memory(cre.starts, v.starts)      # -> True
```

`take` with a slice (and therefore `head`, `tail`, `by_chrom`) shares memory
with the parent; `take` with a mask or an index array copies, as numpy does.
`slop` reuses `codes`, `strands` and the extra columns and allocates only the
shifted starts and ends. Sorting is by a natural-order rank
(`chr1, chr2, ..., chr10, ..., chrX`) kept on the `Genome`, so the order never
depends on which name was seen first.

### Interval kernels on one axis

Whole-set operations place every chromosome on one axis
(`code · 2⁴⁰ + position`), so overlap, merge and nearest over the whole
genome are a few sorts and `searchsorted` calls. `overlaps_any` keeps a
running maximum of reference ends, which works for references of any length;
`overlap_pairs` bounds each query's candidates by the longest reference. Point
lookups (`overlap_rows`) use a cached index built on first use. These
kernels are the default `intervals` backend; the others receive the same
columns and are normalised to the same rules
([Backends]({{ '/backends/' | relative_url }})).

## One Genome per table, re-coded on contact
{: .sec-green }

A `Genome` maps chromosome names to integer codes (and holds sizes when it
knows them). Codes are append-only: a new name gets the next free code and
existing codes never change. Tables made from one another (`take`, `slop`,
`merge`, an Architecture's vertices, a `Genes`' three tables, a `Pairs`' two
anchors) share their Genome, so their codes are directly comparable.

Tables read separately each get their own Genome. An operation between two of
them brings the second onto the first through `Loci._check`: a lookup table
with one entry per chromosome of the other Genome, applied to its codes in
O(n). Names the first Genome has not seen are added, which is why the codes
of the result can differ from the input's:

```python
other = gb.as_loci([("chr2", 450, 650), ("chrX", 1, 2), ("chr1", 1000, 1050)])
other.genome.names, other.codes       # -> (['chr2', 'chrX', 'chr1'], array([0, 1, 2], dtype=int32))
cre.overlap_pairs(other)              # -> (array([0, 5]), array([2, 0]))   row 0 of cre ↔ row 2 of other
cre.genome.names                      # -> ['chr1', 'chr2', 'chrX']          chrX was added
cre._check(other).codes               # -> array([1, 2, 0], dtype=int32)     other, on cre's Genome
```

Pass one Genome to several constructors to skip even that:

```python
g = gb.Genome.from_sizes("hg38.chrom.sizes")        # or Genome.from_fasta("hg38.fa")
cre  = gb.as_loci("peaks.bed", genome=g)
loops = gb.Pairs.make("loops.bedpe", genome=g)
cre.genome is loops.genome                          # -> True
```

There is no session-wide Genome and nothing to set up before the first call;
every table is self-contained and two notebooks' tables never interfere.

### Genes and Pairs are tables of tables

`Genes` is three `Loci` on one Genome: `genes` (gene_id, gene_name,
gene_type), `transcripts` (transcript_id, `gene` → genes row) and `features`
(kind, `transcript` → transcripts row, exon_number). A GTF's 1-based starts
become 0-based while parsing, so genes and CREs compare directly, and a gene's
TSS is the 1-bp interval `[t, t+1)` with `t = start` on `+` and `end - 1` on
`-`.

```python
genes = gb.Genes.make("genes.gtf")
genes.transcripts["gene"]      # -> array([0, 0, 1, 2])          the gene row of each transcript
genes.features["transcript"]   # -> array([0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 3, 3])
genes.features["kind"]         # -> 0 exon, 1 CDS, 2 5'UTR, 3 3'UTR
genes.get_tss().to_records()[:2]
# -> [('chr1', 1000, 1001, '+'), ('chr1', 10999, 11000, '-')]
```

`Pairs` is two row-aligned `Loci`, `P.a` and `P.b`, plus columns (`name`,
`score`, ...). `P.a.genome is P.b.genome`, and `Architecture.make` re-codes
both anchors onto the CREs' Genome before mapping them.

## One edge table, sorted into blocks
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-edges.svg %}
</div><figcaption>
<strong>Views are offsets, not copies.</strong> Edges are kept sorted so each chromosome's cis edges form one contiguous block, in genome order, and every trans edge sits in a final block. <code>A.chrom("chr2")</code>, <code>A.cis</code> and <code>A.trans</code> are slices of the same arrays. Edits through a view write through to the full graph. <code>A.graph()</code> builds the graph of the graph backend (graph-tool when installed, else scipy; igraph or networkx on request).
</figcaption></figure>

An `Architecture` is the CRE `Loci` (vertex *i* = row *i*, vertex columns in
`A.vp`) plus an edge table: `src`, `tgt` (int32, `src < tgt`) and edge
columns in `A.ep` (`w`, `n`, `d`, ...). Edges are stored in canonical order,
cis blocks by chromosome rank then the trans block, sorted by `(src, tgt)`
within a block. `blocks` is the offsets table:

```python
A = gb.Architecture.make(cre, "loops.bedpe", r=100, verbose=False).normalize(verbose=False)
A               # -> Architecture(name='Skeleton', loci=6, links=4 [3 cis · 1 trans], edge_props=[w, d, n], vertex_props=[])
A.src, A.tgt    # -> (array([0, 1, 5, 0], dtype=int32), array([2, 4, 6, 5], dtype=int32))
A.blocks        # -> {'chr1': (0, 2), 'chr2': (2, 3), 'trans': (3, 4)}
A.is_cis        # -> array([ True,  True,  True, False])
```

### Views share memory and write through

```python
c1 = A.chrom("chr1")
c1              # -> Architecture(name='Skeleton:chr1', loci=4, links=2 [2 cis · 0 trans], ...)
np.shares_memory(A.src, c1.src), np.shares_memory(A.ep["w"], c1.ep["w"])   # -> (True, True)
A.cis.n_links, A.trans.n_links          # -> (3, 1)
c1.ep["w"][0] = 50.0
A.ep["w"]                               # -> array([50.,  3.,  2.,  1.])   written through
[(c, v.n_links) for c, v in A.chroms()] # -> [('chr1', 2), ('chr2', 1)]    one view per cis block
```

`A.vp` is shared by every view (same rows, same columns). `subgraph` and
`region` select edges with an index array, so they are copies of the selected
edges, still over the same `Loci`. `copy()` copies everything.

### Neighbours, graphs, algorithms

- **Neighbours** come from one CSR adjacency built over *all* edges, both
  directions, so a CRE's trans partners are always included:
  `A.neighbor_rows(0)` → `(array([2, 5]), array([0, 3]))` (partner rows, edge
  ids), `A.neighbors(0)` → the partner `Loci`.
- **`A.graph(backend=)`** builds the engine's graph from the arrays
  (vertex *i* = row *i*, edge columns as edge properties). The graph-tool
  graph is built once and cached with its edge columns re-synced on every
  call; the others are built on demand.
- **`components`**, **`pagerank`** store their result as a vertex column
  (`A.vp["component"]`, `A.vp["pagerank"]`) and give the same vector on every
  graph backend.
- **Trans loops** are first-class: `normalize` gives them the mean trans
  weight as their expectation, `prune` keeps them, and `dmax` in `make`
  applies to cis loops only.

### Save and load as parquet

`A.save(path)` writes a directory: `vertices.parquet` (the `Loci` columns plus
`vp:<name>` for each vertex column), `edges.parquet` (`src`, `tgt`,
`ep:<name>` for each edge column) and `meta.json`. `Genes.save` writes
`genes.parquet`, `transcripts.parquet`, `features.parquet` and `meta.json`
(with `"coordinates": "0-based half-open"`). `Loci.save`, `Pairs.save` and
`Atlas.save` write single files. Everything is plain parquet:

```python
A.components(); A.save("arch/")
sorted(os.listdir("arch/"))   # -> ['edges.parquet', 'meta.json', 'vertices.parquet']
import pyarrow.parquet as pq
pq.read_table("arch/edges.parquet").schema.names      # -> ['src', 'tgt', 'ep:w', 'ep:d', 'ep:n']
pq.read_table("arch/vertices.parquet").schema.names
# -> ['chrom', 'start', 'end', 'strand', 'name', 'score', 'vp:component']
B = gb.Architecture.load("arch/")                     # same edges, vp['component'] back
import polars as pl
pl.read_parquet("arch/edges.parquet").shape           # -> (4, 5)    any parquet reader
pl.read_parquet("cre.parquet").schema
# -> chrom: Categorical, start: Int64, end: Int64, strand: Categorical, name: String, score: Float64
```

`load(path, genome=g)` reads a table back onto a Genome of your choice; with
none, it gets its own.

## HiChIP short-range tracks and super-enhancers
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-se.svg %}
</div><figcaption>
<strong>A ChIP track from HiChIP pairs.</strong> Pairs closer than about 1 kb are mostly undigested ChIP fragments, so their 5′ ends behave like ChIP-seq reads. <code>hichip</code> turns them into a BED for MACS3 and a 147-bp-fragment coverage bigWig. <code>se.call_se</code> stitches peaks, scores each stitched region by mean signal × width, and cuts at the same slope-1 knee as prime hubs.
</figcaption></figure>

Both modules work on the columns. `hichip.shortrange_ends` reads a HiC-Pro
allValidPairs file (polars when installed, else pandas in chunks), keeps the
cis pairs within `max_dist` and returns the stranded 5′ ends as a sorted
`Loci`; `fragments` extends them to `extsize`-bp fragments clipped to the
chromosome; `coverage` sums them into bedGraph runs and `to_bigwig` writes
those. `shortrange_track` does the whole recipe minus peak calling:

```python
from genomeblocks import hichip, se
ends = hichip.shortrange_ends("sample.allValidPairs", max_dist=1000)   # Loci(n=..., sorted)
out = hichip.shortrange_track("sample.allValidPairs", "sample", "hg38.chrom.sizes")
sorted(out)         # -> ['bed', 'bigwig', 'ends']   BED6 for MACS3 + coverage bigWig
```

The short-range step gives the same ends and a byte-identical bedGraph as the
usual `awk | sort | bedtools genomecov` recipe.

`se.stitch_peaks` is `slop(stitch / 2).merge()` on the interval backend plus
one `overlap_pairs` to find each cluster's first start and last end, so it
returns ROSE's stitched regions with an `n_peaks` column; `se.call_se` scores
them from one `Loci.signal` call (`span=True, n_bins=1`) and cuts at the
knee:

```python
st = se.stitch_peaks(cre, stitch=1000)
st.to_records()[:2], st["n_peaks"]
# -> ([('chr1', 900, 2100, '.'), ('chr1', 4900, 5100, '.')], array([2, 1, 2, 1, 1]))
ses = se.call_se(cre, "h3k27ac.bw", stitch=1000)       # Loci(n=3, chroms=1, cols=[n_peaks, score, rank])
ses, every = se.call_se(cre, "h3k27ac.bw", stitch=1000, return_all=True)
```

## Sharing results
{: .sec-purple }

- **`gb.View`** writes a one-file HTML browser: a small canvas viewer plus the
  tables as gzipped typed arrays. It shows CREs, loops, anchor contact
  profiles, bigWig signal, intervals and gene models, with gene search and a
  CRE inspector, and needs no server. Tracks are added one call at a time
  (`v.signal(...)`, `v.loops(...)`, `v.genes()`) and `v.save("page.html")`
  writes the page.
- **`gb.igv_html`** writes the same data as an igv.js page for people who
  prefer IGV: `igv_html("share.html", regions=[...], loci={"CREs": cre},
  genes=genes, signal={"ATAC": "atac.bw"}, architecture=A)`.

Both read `Genes` as 0-based tables directly (`to_bed12`, `representative`),
so what they draw is what the tables hold.

## Costs
{: .sec-navy }

A `Loci` of a million rows is about 21 MB (four columns); an edge table is
8 bytes per edge plus 8 per edge column. Views cost two integers. Everything
that touches a whole table is a numpy pass or a `searchsorted`; the only
Python-level loops are the writers of text formats. The
[Benchmarks]({{ '/benchmarks/' | relative_url }}) page measures each piece
against the tools people would otherwise use, and against every backend.
