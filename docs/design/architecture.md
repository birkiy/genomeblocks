---
title: Architecture
parent: Design
layout: default
nav_order: 7
---

# Architecture
{: .no_toc }

`Architecture` is a chromatin-contact graph stored as two tables: the CRE
`Loci` are the vertices (row *i* = vertex *i*) and a sorted edge table holds
`src`, `tgt` and the edge columns. Loops make edges, Hi-C contacts weigh
them, and graph algorithms run through the `graph` backend.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## Two tables
{: .sec-purple }

| Part | What it holds |
|---|---|
| `A.loci` | the CRE `Loci`, genome-sorted; vertex *i* is row *i* |
| `A.src`, `A.tgt` | `int32` row numbers, `src < tgt`, one entry per edge |
| `A.ep` | edge columns: `w` (Hi-C weight), `d` (distance), `n` (observed / expected), ... — float arrays aligned to the edges |
| `A.vp` | vertex columns: `annot`, `gene`, `strength`, `component`, ... — arrays aligned to the Loci rows |

```python
import genomeblocks as gb

cre = gb.Loci.make("peaks.bed")
A = gb.Architecture.make(cre, "loops.bedpe", r=100)
# -> [INFO] 4 loops | 4 mapped (100.0%) | loci=6, links=4 (1 trans)
A.loci is cre
# -> True
A.src, A.tgt
# -> (array([0, 1, 5, 0], dtype=int32), array([2, 4, 6, 5], dtype=int32))
A.blocks
# -> {'chr1': (0, 2), 'chr2': (2, 3), 'trans': (3, 4)}
A.degree, A.n_loci, A.n_links, A.n_trans
# -> (array([2, 1, 1, 0, 1, 2, 1]), 6, 4, 1)
```

Because vertex ids are row numbers, every vertex column, annotation and
signal cube built from the same `Loci` lines up with the graph by position:
there is no uid → vertex dictionary. `A[uid]` and `A[uid1, uid2]` still work
through `loci.uids` for lookups by name.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-edges.svg %}
</div><figcaption>
<strong>Views are offsets, not copies.</strong> Edges are kept sorted so each chromosome's cis edges form one contiguous block, in genome order, and every trans edge sits in a final block. <code>A.chrom("chr2")</code>, <code>A.cis</code> and <code>A.trans</code> are slices of the same arrays; they share the vertex columns, and edits through a view write through to the full graph.
</figcaption></figure>

- **Canonical order.** The constructor sorts edges by (block, `src`, `tgt`),
  where a cis edge's block is its chromosome's natural rank and every trans
  edge comes after; `blocks` finds the block boundaries once with `diff`.
- **Views.** `chrom(c)`, `cis`, `trans`, `region(...)`, `subgraph(rows=,
  vp=, values=)` and `A & B` return an `Architecture` over the same `Loci`
  with a slice (or a selection) of the edge arrays.
- **Neighbours** come from one CSR adjacency over all edges (built on first
  use), so a CRE's trans partners are always included; `neighbors(x)`
  returns them with the edge columns as a DataFrame.

```python
A.chrom("chr1").src, A.cis.n_links, A.trans.src, A.trans.tgt
# -> (array([0, 1], dtype=int32), 3, array([0], dtype=int32), array([5], dtype=int32))
A.neighbors(0)
# ->    row                uid chrom    cis    w
# -> 0    2  chr1:4900-5100(+)  chr1   True  0.0
# -> 1    5    chr2:500-600(+)  chr2  False  0.0
```

## The pipeline
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-pipeline.svg %}
</div><figcaption>
<strong>Six chained calls, each writing a column.</strong> Edge columns (<code>ep.w</code>, <code>ep.d</code>, <code>ep.n</code>) and vertex columns (<code>vp.annot</code>, <code>vp.gene</code>, <code>vp.strength</code>) are plain arrays in two dicts. Each step reads what the previous one wrote, so the steps run in this order; a step whose input column is missing raises and names the call that writes it.
</figcaption></figure>

## From loops to weighted edges
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-weights.svg %}
</div><figcaption>
<strong>Loops make edges, Hi-C weighs them, distance normalises them.</strong> <code>make</code> takes every CRE within ± <code>r</code> of each loop anchor's midpoint and links every CRE under one anchor to every CRE under the other. <code>add_mcool</code> maps each CRE to its cooler bin and gives each edge the pixel count of its bin pair, divided among the edges that share that pixel. <code>normalize</code> divides by the expected contact at that distance.
</figcaption></figure>

- **`make(loci, bedpe, r=2500, dmax=1e9, trans=True)`** reads the loops with
  [`as_pairs`]({{ '/design/bedpe/' | relative_url }}) (a BEDPE path, a `Pairs`
  or a frame), re-codes the anchors onto the CREs' `Genome`, drops cis loops
  longer than `dmax` (trans loops stay unless `trans=False`), builds one
  `Loci` of widened anchor midpoints per side, and asks the `intervals`
  backend for `overlap_pairs(anchors, loci)`. The two hit lists are crossed
  per loop (`_cross`, a `repeat`-based group product), self-pairs are
  dropped, and `unique` on `lo · n + hi` keys makes each CRE pair one edge.
  Every interval backend gives the same edge list.
- **`add_mcool(mcool, resolution=)`** maps every CRE to its cooler bin from
  the bin offsets, forms one `(bin1, bin2)` key per edge, and reads the
  pixel table one chromosome block at a time (rows of the lower bin), so
  memory stays bounded and trans pixels are found in the same pass. A
  pixel's count is split evenly across the edges in it, so dense clusters of
  CREs do not multiply one contact. Edges without a pixel keep `w = 0`.
- **`normalize(source="w", name="n")`** fits `w ≈ C · d^−α` on cis edges with
  positive distance and weight and writes `ep.d` and `ep.n = w / E(d)`. The
  fit needs at least three such edges at two distinct distances; without
  them `normalize` warns (a `UserWarning` naming the usable cis-edge count),
  `A.fit` holds `alpha = NaN` and the cis `n` is 0. Trans edges have no
  distance: `ep.d = inf` and their expectation is the mean trans weight.
- **`prune()`** removes co-located cis edges (`d = 0`); trans edges stay.

```python
A.ep["w"][:] = [5.0, 3.0, 2.0, 1.0]
A.normalize()
# -> [INFO] Power-law fit: alpha=0.397, C=1.006e+02 on 3 cis edges; 1 trans edges use the mean trans weight → ep.n
A.ep["d"], A.ep["n"].round(3)
# -> (array([4000., 9000., 4500.,   inf]), array([1.341, 1.11 , 0.562, 1.   ]))
```

## Genes and hubs
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-hubs.svg %}
</div><figcaption>
<strong>Contacts assign genes; strength finds hubs.</strong> <code>annotate</code> labels each CRE with its region class (see <a href="{{ '/design/genes/' | relative_url }}">Genes</a>). A promoter CRE takes its nearest gene. Every other CRE takes the gene of the promoter it contacts most strongly, by the edge column <code>key</code> (O/E by default), cis or trans. <code>strength</code> sums each vertex's incident O/E. <code>prime_hubs</code> ranks CREs by that strength, cuts at the point where the scaled curve's slope reaches 1, and returns the hubs' genes.
</figcaption></figure>

- `annotate(genes)` is `genes.labels(A.loci)` plus `genes.nearest_tss` on
  the promoter rows; the gene assignment is vectorised over the edge arrays:
  one `lexsort` by (vertex, −weight) and the first row per vertex wins. Use a
  different `name=` per weight key (`key="n_CM", name="gene_CM"`) to keep
  several assignments side by side.
- `strength(key)` is two `bincount`s over `src` and `tgt`.
- `elbow(key)` / `prime_hubs` use `se.knee`, the same slope-1 rule on the
  smoothed, 0–1-scaled ranked curve that `call_se` uses for
  super-enhancers.
- `support(genes, r=5000)` finds the linked CREs within `r` of every gene's
  TSS in one `overlap_pairs` call, as `{gene_name: [uid, ...]}` (or row
  arrays with `rows=True`).

```python
genes = gb.Genes.make("genes.gtf")
A.annotate(genes, verbose=False).strength(verbose=False)
A.vp["annot"][:3], A.vp["gene"]
# -> (array(['Promoter-TSS', 'Promoter-TSS', '3UTR'], dtype=object),
#     array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_A', ''], dtype=object))
A.support(genes, r=5000)
# -> {'GENE_A': ['chr1:900-1100(+)', 'chr1:1900-2100(-)', 'chr1:4900-5100(+)'], 'GENE_B': ['chr1:10900-11100(-)']}
```

The chr2 CRE in row 5 is Intergenic and gets GENE_A through its trans edge
to the GENE_A promoter: assignments follow every edge, whatever the block.

## The graph backend seam
{: .sec-purple }

Nothing in the tables depends on a graph library. When an algorithm needs a
graph object, `backends.graph.native(A, backend)` builds one from the arrays
(vertex *i* = Loci row *i*, every edge, edge columns as edge attributes,
`vprops` copied on request) for the engine in use:

| Backend | Object | Default? |
|---|---|---|
| `graph-tool` | `graph_tool.Graph` (built once, cached; edge columns re-synced on every `A.graph()`) | first choice when installed (conda) |
| `scipy` | symmetric CSR adjacency, `scipy.sparse.csgraph` | the pip default |
| `igraph` | `igraph.Graph` | on request |
| `networkx` | `networkx.Graph` (parallel edges merged, their numeric columns summed, as in the adjacency) | on request |

graph-tool and scipy are the only two in the automatic list, because they
give identical answers; `components` renumbers every engine's labels by the
first row of each component (rows without links get -1), and `pagerank`
returns one array per row for every engine, so results do not depend on the
engine. The scipy PageRank is a power iteration on the column-stochastic
matrix with dangling mass spread uniformly, matching the others to `1e-6`.

```python
A.components()
# -> array([ 0,  1,  0, -1,  1,  0,  0])
for b in ("scipy", "igraph", "networkx"):
    print(b, A.components(backend=b).tolist(), A.pagerank("w", backend=b).round(4).tolist())
# -> scipy [0, 1, 0, -1, 1, 0, 0] [0.2178, 0.1626, 0.1786, 0.0244, 0.1626, 0.1466, 0.1074]
# -> igraph [0, 1, 0, -1, 1, 0, 0] [0.2178, 0.1626, 0.1786, 0.0244, 0.1626, 0.1466, 0.1074]
# -> networkx [0, 1, 0, -1, 1, 0, 0] [0.2178, 0.1626, 0.1786, 0.0244, 0.1626, 0.1466, 0.1074]
type(A.graph()).__name__                    # scipy here: graph-tool is not installed
# -> 'csr_matrix'
A.components(backend="graph-tool")                  # without graph-tool
# -> ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool
```

`to_graph_tool`, `to_networkx`, `to_igraph`, `to_scipy(weight)` and
`to_anndata()` (vertices as `obs`, one sparse adjacency per edge column in
`obsp`, the layout scanpy's graph tools read) are the exports; `backends.graph.layout`
gives 2-D positions (the engine's force-directed layout, scipy's spectral
layout, a circle, or genomic position) for drawing.

## Drawing a region
{: .sec-purple }

`A.draw(region, layout="spring" | "circular" | "genomic", backend=None,
merge_distance=None, vertex_size_by=, vertex_color=, edge_width_by="w")`
(`architecture_draw.draw`) takes the linked CREs of a region from
`loci.overlap_rows`, optionally merges CREs closer than `merge_distance`
into one node, sums parallel edges, asks `backends.graph.layout_edges` for
positions, and draws nodes and edges with plain matplotlib. The genomic
layout puts nodes on a line at their position and edges as arcs, a
browser-like view. No graph library is needed for the circular and genomic
layouts; the spring layout uses whichever graph backend is in use.

```python
ax = A.draw("chr1:0-12 kb", layout="genomic")
ax.get_title().splitlines()[-1]
# -> '4 nodes, 2 edges'
ax = A.draw(("chr1", 0, 12_000), merge_distance=1500, vertex_size_by="strength", vertex_color="annot", backend="networkx")
ax.get_title().splitlines()[-1]
# -> '3 nodes (4 CREs merged at 1,500 bp), 2 edges'
```

## Tables out, storing and combining
{: .sec-purple }

- **As a table** an `Architecture` is its edge table: `shape`, `columns`,
  `head()`, `describe()`, `to_pandas()` / `edges_frame()` (`src`, `tgt`,
  uids, chromosomes, `cis`, edge columns), `to_polars`, `to_arrow` and the
  `TableMixin` protocols. `vertices_frame()` is the vertex side with every
  `vp` column; `block_counts()` is the chromosome × chromosome edge count.
- **`from_edges(loci, src, tgt, **ep)`, `from_frame(loci, df)`,
  `from_scipy(loci, M)`** build a graph from arrays, an edge table (rows or
  uids) or a sparse matrix.
- **`save(dir)` / `load(dir)`** write `vertices.parquet` (the Loci plus
  `vp:` columns), `edges.parquet` (`src`, `tgt`, `ep:` columns) and
  `meta.json`, readable from pandas, polars or R.
- **`G | H`, `G & H`** combine graphs on the same `Loci` by edge key;
  `copy()` is deep.

```python
B = gb.Architecture.from_frame(cre, A.edges_frame())
list(B) == list(A)
# -> True
A.save("arch"); back = gb.Architecture.load("arch")
list(back) == list(A), list(back.vp)
# -> (True, ['annot', 'gene', 'strength', 'component', 'pagerank'])
import polars as pl
pl.DataFrame(A).shape
# -> (4, 10)
```

## Costs
{: .sec-purple }

- **`make`:** two `overlap_pairs` (anchors against CREs), one group cross
  product, one `unique`; no Python loop over loops.
- **`add_mcool`:** one pass over the pixel table, chromosome block by block,
  with a `searchsorted` per block.
- **`normalize`, `strength`, `annotate`:** a `curve_fit` on the cis edges,
  `bincount`s, one `lexsort` over the edges.
- **Views:** two offsets, no copy; the adjacency is one `argsort` of
  `2 · n_links` entries, built once.
- **Graph algorithms:** the engine's cost plus one array copy in and out;
  the graph-tool object is cached.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
