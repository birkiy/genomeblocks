---
title: Quickstart
layout: default
nav_order: 3
---

# Quickstart
{: .no_toc }

A 10-minute tour of `genomeblocks` 2.0. We load ATAC peaks as a table, do set algebra on them, annotate them with a GTF, build a chromatin-contact graph from loops and a Hi-C matrix, extract bigWig signal into a cube, draw a browser view, and hand everything to pandas and polars.
{: .fs-5 .fw-300 }

The outputs shown (`# ->`) come from the tiny synthetic dataset in the test suite (`tests/conftest.py`): 7 peaks, 3 genes, 4 loops, two bigWigs and a one-resolution `.cool`. Swap in your own files; nothing else changes.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Load peaks as `Loci`

```python
import genomeblocks as gb
from genomeblocks import Loci, as_loci

peaks = Loci.make("atac_peaks.narrowPeak", keep=True)   # keep the narrowPeak columns too
peaks
# -> Loci(n=7, chroms=2, sorted, cols=[name, score, signalValue, pValue, qValue, peak])
peaks.shape, peaks.columns
# -> (7, 10), ['chrom', 'start', 'end', 'strand', 'name', 'score', 'signalValue', 'pValue', 'qValue', 'peak']

peaks[0]                      # -> Locus[0](chr1:900-1100(.))   a view on row 0
peaks[0].uid                  # -> 'chr1:900-1100(.)'
peaks[0].signalValue          # -> 1.0                           extra columns are attributes
peaks["start"][:3]            # -> array([ 900, 1900, 4900])     columns are numpy arrays
peaks.describe()              # rows, chromosomes, bases covered, length quantiles, strands
```

A `Loci` is a table: `codes` / `starts` / `ends` / `strands` columns on a `Genome`, plus `cols` for everything else. `Loci.make` sorts into genome order by default. The same table comes from anything interval-like through `as_loci`, which every function below calls on its inputs:

```python
import pandas as pd

as_loci("chr1:1,000-2,000")                               # -> Loci(n=1, chroms=1)
as_loci([("chr1", 100, 200), ("chr2", 50, 80)])           # -> Loci(n=2, chroms=2)
df = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [100, 300], "end": [500, 700]})
as_loci(df)                                               # -> Loci(n=2, chroms=1)
Loci.from_frame(df)                                       # the same, by name
```

---

## 2. Interval algebra

```python
# Widen each peak by ±200 bp, then fuse overlapping or book-ended rows.
cre = peaks.slop(200).merge()
cre                                   # -> Loci(n=7, chroms=2, sorted)

se = Loci.make("H3K27ac_SE.bed")      # three super-enhancer regions
cre & se                              # -> Loci(n=4, chroms=1, sorted)   rows of cre overlapping se
cre - se                              # -> Loci(n=3, chroms=2, sorted)   rows of cre not overlapping se
cre | se                              # -> Loci(n=10, chroms=2)          concatenation
cre ^ se                              # -> Loci(n=4, chroms=2)           rows of either not overlapping the other

cre.overlap_any(se)                   # -> array([ True,  True, False,  True,  True, False, False])
j, dist = cre.nearest(se)             # row in se and gap in bp, per row of cre (0 = overlap)
cre.overlaps("chr1:900-2,000")        # -> Loci(n=2, chroms=1)           a point lookup (cached index)
```

`&` keeps whole rows of the left side (bedtools `intersect -u`), `-` drops them (`-v`). Results carry the left side's extra columns. Every whole-set operation takes `backend=`; see [Backends]({{ '/backends/' | relative_url }}) and the [Loci guide]({{ '/guide/loci/' | relative_url }}).

---

## 3. Annotate with a GTF

`Genes` holds three linked tables — `genes`, `transcripts`, `features` — each a `Loci` in the **same 0-based, half-open coordinates as every other table**. A GTF's 1-based start becomes `start - 1` while parsing (genePred tables are already 0-based).

```python
from genomeblocks import Genes

genes = Genes.make("genes.gtf", promoter_r=1000)         # GTF or GFF3, .gz too
genes                                                    # -> Genes(3 genes, 4 transcripts, 7 exons)
genes.to_pandas()
#   chrom  start    end strand gene_id gene_name       gene_type
# 0  chr1   1000   5000      +      G1    GENE_A  protein_coding
# 1  chr1  10000  11000      -      G2    GENE_B          lncRNA
# 2  chr2  20000  30000      +      G3    GENE_C  protein_coding
genes.to_pandas("transcripts")                           # transcript_id, gene (row), gene_id, gene_name
genes.to_pandas("features")                              # kind, transcript (row), exon_number, feature, ...

g = genes["GENE_A"]                 # -> Gene(GENE_A, chr1:1,000-5,000(+), protein_coding)
g.tss                               # -> Locus(chrom='chr1', start=1000, end=1001, strand='+')
g.transcripts                       # -> Loci(n=2, chroms=1, cols=[transcript_id, gene])
g.exons                             # -> Loci(n=4, chroms=1, cols=[kind, transcript, exon_number])
```

A gene's TSS is the 1-bp interval `[t, t+1)` with `t = start` on `+` and `t = end - 1` on `-`:

```python
genes.get_tss().to_pandas()
#   chrom  start    end strand gene_name gene_id  gene
# 0  chr1   1000   1001      +    GENE_A      G1     0
# 1  chr1  10999  11000      -    GENE_B      G2     1
# 2  chr2  20000  20001      +    GENE_C      G3     2
```

Region classes and nearest genes for any interval set, one row per input row:

```python
annot = genes.annotations(cre)      # DataFrame(uid, annotation)
annot.groupby("annotation").size()
# annotation
# 3UTR            1
# Intergenic      2
# Promoter-TSS    4

genes.nearest_genes(cre)            # DataFrame(uid, gene_name, distance) to the nearest promoter window
```

Promoter windows are `get_tss().slop(promoter_r)`; use them to split a CRE set, and the split feeds the heatmaps below as a plain dict:

```python
prom = genes.get_tss().slop(1000)
groups = {"promoters": cre & prom, "enhancers": cre - prom}
{k: len(v) for k, v in groups.items()}        # -> {'promoters': 4, 'enhancers': 3}
```

### Isoforms your cells use
{: .no_toc }

`select_isoforms` points each gene at the longest isoform whose TSS is open in your peaks and/or bigWig, moves the gene body onto it, and records the choice:

```python
genes.select_isoforms(cre, "ATAC.bw", r=500)
# [INFO] Isoform support: 2/3 genes with an open TSS (3/4 isoforms), 1 fell back to the longest isoform.
genes["GENE_A"].canonical          # -> 'T1'
genes.to_pandas()[["gene_name", "canonical"]]
#   gene_name canonical
# 0    GENE_A        T1
# 1    GENE_B        T2
# 2    GENE_C        T3
```

The same happens at parse time with `Genes.make(path, cre=..., bw=...)`. See the [Genes guide]({{ '/guide/genes/' | relative_url }}).

---

## 4. Loops as `Pairs`, contacts as an `Architecture`

A BEDPE is a `Pairs` table: anchor `a` and anchor `b` are two row-aligned `Loci`, the other columns ride along.

```python
from genomeblocks import Pairs, Architecture

loops = Pairs.make("loops.bedpe")
loops                                   # -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
loops.to_pandas()
#   chrom1  start1  end1 chrom2  start2   end2 name  score strand1 strand2
# 0   chr1     900  1100   chr1    4900   5100   l1    5.0       +       -
# 1   chr1    1900  2100   chr1   10900  11100   l2    3.0       +       +
# 2   chr1     900  1100   chr2     500    600   l3    1.0       -       +
# 3   chr2     500   600   chr2    5000   5100   l4    2.0       +       +
loops.distance                          # -> array([4000., 9000.,   inf, 4500.])   inf = trans
loops.overlapping(cre, r=100)           # pairs with an anchor on a CRE (bedtools pairtobed)
```

`Architecture.make` draws an edge between every pair of CREs within `r` bp of the two anchor midpoints. The vertices **are** the `Loci` rows; the edges are a table (`src`, `tgt`, then the `ep` columns). `add_mcool` and `normalize` work on the graph's own loci — they take no `loci` argument:

```python
A = (Architecture.make(cre, loops, r=200)             # or a .bedpe path, or a BEDPE frame
       .add_mcool("hic.cool")                          # .mcool: add resolution=5000
       .normalize())                                   # O/E: ep.d (distance), ep.n (observed / expected)
# [INFO] 4 loops | 4 mapped (100.0%) | loci=6, links=4 (1 trans)
# [INFO] Set distributed weights for 4/4 edges from cooler. [w]
# [INFO] Power-law fit: alpha=0.802, C=4.143e+03 on 3 cis edges; 1 trans edges use the mean trans weight → ep.n
A
# -> Architecture(name='Skeleton', loci=6, links=4 [3 cis · 1 trans], edge_props=[w, d, n], vertex_props=[])
A.to_pandas()                                          # the edge table: src, tgt, uid1, uid2, ..., cis, w, d, n
```

Annotate the vertices, sum the edge weights per vertex, and pick the hubs — the `vp` columns line up with `cre` by row:

```python
A.annotate(genes).strength()
A.vp.annot          # -> array(['Promoter-TSS', 'Promoter-TSS', '3UTR', 'Promoter-TSS', 'Promoter-TSS', 'Intergenic', 'Intergenic'])
A.vp.gene           # -> array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_A', ''])
A.vp.strength       # -> array([2.12, 1.08, 1.12, 0.  , 1.08, 1.82, 0.82])

res = A.prime_hubs(key="n")
sorted(res["prime_genes"])             # -> ['GENE_A']
```

Views are zero-copy slices of the edge table, and graph algorithms run on whichever graph backend is installed (graph-tool when present, scipy otherwise; igraph and networkx on request):

```python
A.cis, A.trans, A.chrom("chr1")        # Skeleton:cis · Skeleton:trans · Skeleton:chr1
A.components()                         # -> array([ 0,  1,  0, -1,  1,  0,  0])   -1 = no links
A.pagerank()                           # one value per vertex, stored as vp.pagerank
A.to_networkx()                        # -> Graph with 7 nodes and 4 edges
ax = A.draw("chr1:0-12 kb", layout="genomic")      # the region's subgraph on matplotlib
```

Full API: [Architecture guide]({{ '/guide/architecture/' | relative_url }}).

---

## 5. Signal cubes and heatmaps

`loci.signal(bigwigs)` returns one `(rows, tracks, bins)` array — row `i` is row `i` of the `Loci`:

```python
S = cre.signal(["ATAC.bw", "H3K27ac.bw"], n_bins=50, flank=1000)
S.shape, S.dtype                       # -> (7, 2, 50), dtype('float32')
N = gb.tmm(S)                          # TMM-normalised cube, same shape

fig = cre.plot_heatmap(S, groups=groups, samples=["ATAC", "H3K27ac"], height=1000, vmax=8)
fig.savefig("heatmap.pdf")
fig = cre.plot_profiles(S, groups=groups, height=1000)    # mean profile of track 0 per group
```

`compare_heatmap` extracts and draws two interval sets as one grid (a-only, shared, b-only):

```python
fig, union, S2, blocks = gb.compare_heatmap(
    cre, se, ["ATAC.bw", "H3K27ac.bw"],
    a_name="CRE-only", b_name="SE-only", n_bins=50, flank=1000)
union                                  # -> Loci(n=8, chroms=2)
{k: int(v.sum()) for k, v in blocks.items()}   # -> {'CRE-only': 3, 'common': 4, 'SE-only': 1}
```

The cube also leaves as a labelled `xarray.DataArray` (`genomeblocks.interop.cube_to_xarray(S, cre, tracks)`), a long DataFrame or an AnnData. See the [Signal guide]({{ '/guide/signal/' | relative_url }}).

---

## 6. Browser view, IGV page, interactive view

Three ways to look at a region. All of them take the tables as they are:

```python
# a matplotlib figure: bigWig paths, Loci, Pairs / BEDPE paths, Genes, BAM paths
fig, axes = gb.browser("chr1:0-12 kb",
                       {"ATAC": "ATAC.bw", "CREs": cre, "loops": loops, "genes": genes},
                       figsize=(10, None))
fig.savefig("browser.svg")
list(axes)                             # -> ['ATAC', 'CREs', 'loops', 'genes', '_axis']

# one HTML file that opens igv.js on the regions, with the tracks embedded
gb.igv_html("igv.html", regions=["chr1:0-12 kb"], loci={"CREs": cre}, genes=genes,
            signal={"ATAC": "ATAC.bw"}, architecture=A)

# one interactive HTML file over an Architecture: signal, CREs, loops and genes
v = gb.View(A, genes=genes).signal("ATAC", "ATAC.bw").cres().loops().genes()
v.save("view.html")
```

A full real-data example is in [`examples/ar_foxa1_lncap/`](https://github.com/birkiy/genomeblocks/tree/main/examples/ar_foxa1_lncap) and the [browser walkthrough]({{ '/walkthrough/browser/' | relative_url }}). See the [Browser guide]({{ '/guide/browser/' | relative_url }}).

---

## 7. Hand off to pandas, polars and the rest

Every table converts explicitly, and the protocols let other libraries take it as it is:

```python
cre.to_pandas()                        # chrom (categorical, genome order), start, end, strand, + cols
cre.to_polars()                        # through Arrow; numeric columns are not copied
A.to_polars()                          # the edge table
cre.to_pyranges(); cre.to_bioframe(); cre.to_bed("cre.bed")
cre.save("cre.parquet"); Loci.load("cre.parquet")

import polars as pl, duckdb
pl.DataFrame(cre).shape                # -> (7, 4)                      Arrow C stream
duckdb.sql("select chrom, count(*) n from cre group by chrom order by chrom").df()
#   chrom  n
# 0  chr1  5
# 1  chr2  2

Loci.from_frame(cre.to_polars()).equals(cre)       # -> True   round trips are exact
```

See [Interoperability]({{ '/interoperability/' | relative_url }}) for every input and output.

---

## Next steps

- [Concepts →]({{ '/concepts/' | relative_url }}) — the mental model behind the tables.
- [Backends →]({{ '/backends/' | relative_url }}) — which engine runs what, and how to pick one.
- [Example: AR & FOXA1 →]({{ '/walkthrough/' | relative_url }}) — a complete real-data walkthrough.
- [Architecture guide →]({{ '/guide/architecture/' | relative_url }}) — build and mine chromatin networks.
