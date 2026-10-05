---
title: Architecture
parent: User Guide
layout: default
nav_order: 4
---

# Architecture
{: .no_toc }

Chromatin-contact networks as two tables: build from loop files, overlay Hi-C
from `.cool` / `.mcool`, normalize to O/E, annotate vertices with gene context,
mine hub genes, and hand the graph to graph-tool, igraph, networkx or scipy.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is
{: .sec-purple }

`gb.Architecture` is a graph stored as two tables that share one key:

| Table | What it holds | Key |
|---|---|---|
| vertices | a [`Loci`]({{ '/guide/loci/' | relative_url }}) — `A.loci` — plus vertex columns in `A.vp` | row number: vertex *i* is `A.loci` row *i* |
| edges | `A.src`, `A.tgt` (int32 row numbers, `src < tgt`) plus edge columns in `A.ep` | edge number |

Because vertices *are* Loci rows, a signal cube, a motif matrix or any other
array built from the same Loci lines up with the vertex columns by position.
There is no uid-to-vertex dictionary.

The edge columns the pipeline creates:

| Edge column | Set by | Meaning |
|---|---|---|
| `ep.w` | `make` (zeros), `add_mcool` | raw contact weight |
| `ep.d` | `normalize` | genomic distance between the two CRE centres; `inf` for trans edges |
| `ep.n` | `normalize` | observed / expected weight |

The vertex columns the pipeline creates: `vp.annot` and `vp.gene`
(`annotate`), `vp.strength` (`strength`), `vp.component` (`components`),
`vp.pagerank` (`pagerank`). `A.ep.w` and `A.ep["w"]` are the same thing.

Edges are kept sorted so that the cis edges of each chromosome form one
contiguous block in genome order, followed by one block of trans edges:

```
| cis chr1 | cis chr2 | ... | cis chrX | trans |
```

`A.blocks` gives the edge ranges, and the views below are zero-copy slices of
them. As a table the Architecture is its edge table: `shape`, `columns`,
`head()`, `describe()`, `to_pandas()` / `to_polars()` / `to_arrow()`, and the
Arrow and dataframe protocols, so `pl.DataFrame(A)` or `duckdb.sql("select
* from A")` work as they are. `len(A)` counts the vertices with at least one
edge.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-edges.svg %}
</div><figcaption>
<strong>Views are slices.</strong> Edges are stored sorted by block — each chromosome's cis edges together, the trans edges last — so <code>A.chrom('chr2')</code>, <code>A.cis</code> and <code>A.trans</code> are two offsets into the same arrays, not copies, and <code>A.blocks</code> is that offset table. Neighbour lookups use one adjacency over every block, so a trans partner is never missed; <code>A.graph(backend=)</code> builds the engine's object from the arrays on demand.
</figcaption></figure>

{: .note }
> The core pipeline is small: **make → add_mcool → normalize → annotate →
> strength → prime_hubs**. Drawing lives in `genomeblocks.architecture_draw`
> (reached as `A.draw`) so the graph object stays free of matplotlib, and graph
> algorithms run through the graph backend so the object stays free of
> graph-tool.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-pipeline.svg %}
</div><figcaption>
<strong>Six calls, each writing one column.</strong> <code>make</code> builds the vertex and edge tables; <code>add_mcool</code> fills <code>ep.w</code>; <code>normalize</code> adds <code>ep.d</code> and <code>ep.n</code>; <code>annotate</code> writes <code>vp.annot</code> and <code>vp.gene</code>; <code>strength</code> sums <code>ep.n</code> per vertex into <code>vp.strength</code>; <code>prime_hubs</code> reads that column back and returns genes. Every method returns the object, so the chain reads top to bottom.
</figcaption></figure>

---

## Build from loops
{: .sec-purple }

`Architecture.make` assumes you already have a CRE catalogue. Each loop
anchor is mapped to the CREs within `±r` bp of its midpoint; every
(CRE at anchor 1, CRE at anchor 2) pair becomes one edge. Duplicate pairs
collapse into one edge.

```python
import genomeblocks as gb

cre = gb.Loci.make("cre.bed")
A = gb.Architecture.make(cre, "loops.bedpe", r=100)        # r: bp around each anchor midpoint (default 2500)
# -> [INFO] 4 loops | 4 mapped (100.0%) | loci=6, links=4 (1 trans)
A
# -> Architecture(name='Skeleton', loci=6, links=4 [3 cis · 1 trans], edge_props=[w], vertex_props=[])
```

| arg | meaning |
|---|---|
| `bedpe` | a BEDPE path, a [`Pairs`]({{ '/guide/bedpe/' | relative_url }}), or a pandas / polars / arrow frame with BEDPE columns |
| `r` | radius (bp) around each anchor midpoint to catch CREs (default 2500) |
| `dmax` | skip **cis** loops whose anchor midpoints are farther apart than this (default 1e9) |
| `trans` | keep inter-chromosomal loops (default `True`) |
| `backend` | interval engine for the anchor-to-CRE mapping (`'genomeblocks'`, `'cgranges'`, `'ncls'`, `'bioframe'`, `'pyranges'`, `'bedtools'`); every engine gives the same edges |
| `name` | the graph's name, shown in `repr` and plot titles (default `'Skeleton'`) |

The Loci must be in genome order; `Loci.make` sorts by default. Unsorted input
is sorted into a copy, with a warning, and vertex numbers then refer to
`A.loci`, not to the input.

Right after `make`, `A.ep.w` is all zeros. Read the edges and the blocks:

```python
A.head()
# ->    src  tgt               uid1                 uid2 chrom1 chrom2    cis    w
# -> 0    0    2   chr1:900-1100(+)    chr1:4900-5100(+)   chr1   chr1   True  0.0
# -> 1    1    4  chr1:1900-2100(-)  chr1:10900-11100(-)   chr1   chr1   True  0.0
# -> 2    5    6    chr2:500-600(+)    chr2:5000-5100(+)   chr2   chr2   True  0.0
# -> 3    0    5   chr1:900-1100(+)      chr2:500-600(+)   chr1   chr2  False  0.0
A.blocks
# -> {'chr1': (0, 2), 'chr2': (2, 3), 'trans': (3, 4)}
list(A)                      # (src, tgt) row pairs
# -> [(0, 2), (1, 4), (5, 6), (0, 5)]
A.degree                     # one value per Loci row
# -> array([2, 1, 1, 0, 1, 2, 1])
A.n_links, A.n_trans, A.n_loci
# -> (4, 1, 6)
```

### Other constructors

```python
# from vertex rows: src[k]–tgt[k] is edge k; keyword arrays become edge columns
B = gb.Architecture.from_edges(cre, [0, 0, 1], [2, 5, 4], w=[5.0, 2.0, 6.0])

# from an edge table (pandas / polars / arrow) with src / tgt rows or uid1 / uid2
B = gb.Architecture.from_frame(cre, A.edges_frame())
B = gb.Architecture.from_frame(cre, A.edges_frame()[["uid1", "uid2", "w"]])

# from a square sparse matrix over the Loci rows (upper triangle read)
B = gb.Architecture.from_scipy(cre, A.to_scipy("w"), weight="w")
```

`from_frame` reads what `edges_frame()` writes: every other numeric column
becomes an edge column. A uid that is not in the Loci raises.

---

## Adding Hi-C weights
{: .sec-purple }

```python
A.add_mcool("hic.mcool", resolution=5000)        # .mcool: pick a resolution
A.add_mcool("hic.cool")                          # .cool: single resolution
# -> [INFO] Set distributed weights for 4/4 edges from cooler. [w]
A.ep.w
# -> array([6., 3., 4., 9.])
```

Each CRE is assigned to the bin holding its start; the count of the (bin, bin)
pixel the edge spans is shared equally among the edges that fall in the same
pixel, and lands in `ep.w` (or the `name=` you pass — call it again with
another file and name to stack weights from several experiments). Pixels are
read one chromosome block at a time, so memory stays bounded and trans pixels
are found in the same pass.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-weights.svg %}
</div><figcaption>
<strong>From a loop to an O/E weight in three steps.</strong> <code>make</code> links every CRE under one anchor (midpoint ± <code>r</code>) to every CRE under the other, never two CREs under the same anchor. <code>add_mcool</code> reads the pixel each edge spans and shares its count equally among the edges in that pixel, so two edges in one 12-count pixel get <code>ep.w = 6</code> each. <code>normalize</code> fits <code>E(d) = C / d^α</code> over the cis edges and stores <code>w / E(d)</code> in <code>ep.n</code>; a trans edge has <code>d = ∞</code> and is divided by the mean trans weight instead.
</figcaption></figure>

The file type is checked: a multi-resolution file without `resolution=`
raises and lists the resolutions it holds; a single-resolution file with
`resolution=` raises too.

---

## Normalizing to O/E, and pruning
{: .sec-purple }

```python
A.normalize(source="w", name="n")
# -> [INFO] Power-law fit: alpha=0.802, C=4.143e+03 on 3 cis edges; 1 trans edges use the mean trans weight → ep.n
A.ep.d
# -> array([4000., 9000., 4500.,   inf])
```

`normalize` fits a power law `w ≈ C · d^-α` over the cis edges with positive
distance and weight, divides each cis weight by its expectation, and stores
the ratio in `ep.n`. The fitted `alpha` and `C` are kept in `A.fit`. Trans
edges have no distance: `ep.d` is `inf` for them and their expectation is the
mean trans weight, so a trans `ep.n` reads as "times the average trans
contact". The fit needs at least three cis edges with a positive weight
and distance, at two distinct distances or more. With fewer, `normalize`
issues a `UserWarning` that names the usable cis-edge count, `A.fit` holds
`nan` and the cis `ep.n` is 0.

```python
A.prune()          # drop zero-distance (co-located) cis edges; trans edges stay
# -> [INFO] prune: removed 0 zero-distance edges → 4 edges.
```

{: .warning }
> `prune` and `annotate` need `normalize` first: `prune` reads `ep.d`, and
> `annotate` ranks promoter contacts by `ep.n`. Both raise a `ValueError`
> naming the missing column when it is not there.

---

## Annotating vertices
{: .sec-purple }

```python
genes = gb.Genes.make("genes.gtf")
A.annotate(genes, key="n", name="gene")
# -> [INFO] Annotated 6 loci: 3 promoter CREs | 2/3 non-promoter CREs assigned to a top-'n' promoter gene → vp.gene.
A.vp.annot
# -> array(['Promoter-TSS', 'Promoter-TSS', '3UTR', 'Promoter-TSS', 'Promoter-TSS', 'Intergenic', 'Intergenic'], dtype=object)
A.vp.gene
# -> array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_A', ''], dtype=object)
```

Two-stage gene assignment:

| Vertex column | Value |
|---|---|
| `vp.annot` | region class from the [`Genes`]({{ '/guide/genes/' | relative_url }}) tables: `Promoter-TSS`, `5UTR`, `3UTR`, `Exonic`, `Intronic`, `Intergenic` |
| `vp.gene` | promoter CREs → the gene of their nearest TSS; every other CRE → the gene of its highest-`key` promoter neighbour, cis or trans (`''` if none) |

Both columns have one value per Loci row, including rows without edges. Pass a
distinct `name=` per `key` to keep several assignments side by side.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-hubs.svg %}
</div><figcaption>
<strong>Enhancers take their gene from their strongest promoter contact; hubs are cut at the knee.</strong> <code>annotate</code> gives a promoter CRE the gene of its nearest TSS and every other CRE the <code>vp.gene</code> of the promoter neighbour with the highest <code>ep[key]</code> (O/E by default), cis or trans. <code>prime_hubs</code> ranks vertices by <code>vp.strength</code>, cuts the [0, 1]-scaled curve where its slope reaches 1, and collects the hubs' genes — split into promoter and enhancer hubs by <code>vp.annot</code>.
</figcaption></figure>

---

## Hub genes
{: .sec-purple }

### Node strength — `strength`

```python
A.strength(key="n", name="strength")     # vp.strength = sum of ep.n per vertex
# -> [INFO] Summed ep.n → vp.strength (node strength).
```

Two `bincount`s over `src` and `tgt`. No normalization is applied.

### Cutoff — `elbow`

```python
cutoff, uids = A.elbow("strength")
# -> [INFO] Slope-1 on vp.strength (value≈1.12 at cutoff): cutoff at 3/6 (50.0%)
hub_uids = uids[:cutoff]
```

Sorts vertices by `vp.strength` descending and finds the **slope-1 knee** of
the normalized curve (the same rule `genomeblocks.se` uses for
super-enhancers). Returns the cutoff and the uids in descending order.

### The one-liner — `prime_hubs`

```python
res = A.prime_hubs(key="n")
# -> [INFO] Slope-1 on vp.strength (value≈1.12 at cutoff): cutoff at 3/6 (50.0%)
# -> [INFO] Prime hubs: 3 hubs (1 promoters, 2 enhancers)
# -> [INFO] Prime genes: 1 = 1 promoter + 1 enhancer (overlap: 1)
res["prime_genes"]       # genes of all hub CREs
# -> {'GENE_A'}
res["promoter_genes"]    # genes of hub CREs that are promoters
res["enhancer_genes"]    # genes of hub CREs that are not promoters
res["hub_uids"]          # all hub CRE uids, strongest first
# -> ['chr1:900-1100(+)', 'chr2:500-600(+)', 'chr1:4900-5100(+)']
res["promoter_uids"], res["enhancer_uids"]
# -> (['chr1:900-1100(+)'], ['chr2:500-600(+)', 'chr1:4900-5100(+)'])
res["cutoff"]            # number of hubs
# -> 3
```

`prime_hubs` runs `strength` (if `vp.strength` is missing) → `elbow` → splits
the hubs by `vp.annot`, collecting each hub's `vp.gene`. It needs `vp.annot`
and `vp.gene` from `annotate`.

---

## CREs around a gene — `support`
{: .sec-purple }

```python
A.support(genes, r=5000)
# -> {'GENE_A': ['chr1:900-1100(+)', 'chr1:1900-2100(-)', 'chr1:4900-5100(+)'], 'GENE_B': ['chr1:10900-11100(-)']}
A.support(genes, r=5000, rows=True)["GENE_A"]
# -> array([0, 1, 2])
```

For every gene at once: the linked CREs within TSS ± `r`. The TSS is the gene's
5′-most base (`start` on `+`, `end - 1` on `-`). `mode="overlap"` (default)
keeps CREs that intersect the window, `mode="center"` only those whose midpoint
is inside it. `rows=True` returns row arrays instead of uids; `linked=False`
includes CREs without edges. The single-window form is `A.near("chr1:0-3,000")`
(a Loci) or `A.near_rows(...)`.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-support.svg %}
</div><figcaption>
<strong><code>support()</code> keeps the linked CREs inside <code>TSS ± r</code>.</strong> With <code>mode='overlap'</code> (the default) a CRE counts when its interval touches the window, so <code>a</code> and <code>b</code> both count — <code>b</code> straddles the window's edge. With <code>mode='center'</code> only the midpoint decides, and <code>b</code>'s midpoint lies outside, so <code>a</code> alone remains. <code>c</code> is inside the window but has no edge in the Architecture, so it is left out unless <code>linked=False</code>; <code>d</code> is outside; <code>e</code> is the partner of <code>a</code> and <code>b</code>. <code>rows=True</code> returns row numbers in place of uids.
</figcaption></figure>

---

## Views: chromosome, cis, trans, region, subgraph
{: .sec-purple }

Every view is an Architecture over the same Loci and the same vertex columns;
slices of the edge blocks cost nothing.

```python
A.chrom("chr1")                        # cis edges of one chromosome
# -> Architecture(name='Skeleton:chr1', loci=4, links=2 [2 cis · 0 trans], ...)
A.cis                                  # every cis edge
A.trans                                # every trans edge
A.region("chr1:0-6 kb")                # both ends inside the region
A.region("chr1:0-6 kb", both=False)    # one end inside: trans partners too
A.subgraph(rows=[0, 2, 5])             # both ends among these rows
A.subgraph(vp="annot", values="Promoter-TSS")
A.subgraph(mask=A.vp.strength > 0.5)
for chrom, view in A.chroms():         # fan out per chromosome
    ...
```

Regions are written as `'chr1:1,000-2,000'`, `'chr8:127.5-128.5 Mb'`,
`(chrom, start, end)` or a `Locus`.

---

## Neighbours
{: .sec-purple }

```python
A.neighbor_rows(0)                 # (partner rows, edge ids), cis and trans
# -> (array([2, 5], dtype=int32), array([0, 3]))
A.neighbors(0)                     # the same as a DataFrame with the edge columns
# ->    row                uid chrom    cis    w       d         n
# -> 0    2  chr1:4900-5100(+)  chr1   True  6.0  4000.0  1.124534
# -> 1    5    chr2:500-600(+)  chr2  False  9.0     inf  1.000000
A["chr1:900-1100(+)"]              # {partner uid: w}
# -> {'chr1:4900-5100(+)': 6.0, 'chr2:500-600(+)': 9.0}
A["chr1:900-1100(+)", "chr1:4900-5100(+)"]   # one edge's columns (None if absent)
# -> {'w': 6.0, 'd': 4000.0, 'n': 1.1245336516005273}
"chr1:900-1100(+)" in A            # True when the CRE has at least one edge
```

A CRE is named by its row, its uid, or a `Locus` read from the table
(`cre[0]`). Neighbours come from one adjacency index over all edges, so trans
partners are always included.

---

## Graph algorithms and graph backends
{: .sec-purple }

Graph algorithms run through the `graph` backend family: graph-tool when
installed, else scipy (identical results); igraph and networkx on request.
Results are plain numpy arrays aligned to the Loci rows and are also stored as
vertex columns.

```python
A.components()                     # component per row, -1 without edges -> vp.component
# -> array([ 0,  1,  0, -1,  1,  0,  0])
A.components(backend="networkx")   # the same labels from any engine
A.pagerank("n")                    # weighted by an edge column -> vp.pagerank
A.pagerank(backend="igraph", damping=0.85)
```

Components are numbered by their first row, so every engine agrees. Switch
the engine for one call with `backend=`, or for a block of code:

```python
with gb.use_backend(graph="networkx"):
    A.components()
    A.pagerank("n")
```

A backend that is requested but not installed raises — there is no silent
switch:

```python
A.components(backend="graph-tool")                  # without graph-tool
# -> ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool
```

Get the engine's own graph object, with vertex *i* = Loci row *i*, every edge
(cis and trans), and the edge columns as edge attributes:

```python
g = A.graph()                       # the default engine's graph (graph-tool, else a scipy CSR)
g = A.graph(backend="networkx", vprops=["strength"])

A.to_networkx()                     # networkx.Graph; node attrs chrom / start / end + vprops
A.to_igraph(vprops=["gene"])        # igraph.Graph
A.to_scipy("n")                     # symmetric CSR over the Loci rows; weight=None -> 1 per edge
A.to_anndata()                      # obs = vertices (+ vp), obsp['w'], obsp['n'], ... adjacencies
A.to_graph_tool()                   # needs graph-tool
```

The graph-tool graph is built once and cached; its edge properties are
re-synced from `A.ep` on every `graph()` call. networkx holds one edge per
pair, so parallel edges (possible after `from_edges`) are merged and their
numeric columns summed — the same sum `to_scipy` puts in the adjacency.

See [Backends]({{ '/backends/' | relative_url }}) for the full list and
[Interoperability]({{ '/interoperability/' | relative_url }}) for what each
export carries.

---

## Tables out
{: .sec-purple }

```python
A.edges_frame()        # one row per edge: src, tgt, uid1, uid2, chrom1, chrom2, cis, edge columns
A.to_pandas()          # the same; to_polars() / to_arrow() too
A.vertices_frame()     # one row per linked vertex: coordinates, uid, vertex columns (index = row)
A.vertices_frame(linked=False)
A.block_counts()       # chrom x chrom edge counts, cis on the diagonal
# ->       chr1  chr2
# -> chr1     2     1
# -> chr2     1     1
A.describe()           # loci, linked loci, edges (cis / trans / blocks), columns
A.shape, A.columns     # of the edge table
```

---

## Drawing a region
{: .sec-purple }

`A.draw` plots the CREs of a region as nodes and the edges among them as
lines, with plain matplotlib. Positions come from the graph backend's layout.

```python
ax = A.draw("chr8:127.5-128.5 Mb")                 # spring layout, default engine
ax = A.draw(("chr1", 0, 12_000),
            layout="genomic",            # 'spring' | 'circular' | 'genomic' (nodes on a line, edges as arcs)
            merge_distance=1500,         # collapse CREs within 1.5 kb into one node
            vertex_size_by="strength",   # a vertex column or an array per Loci row
            edge_width_by="n",           # an edge column
            vertex_color="annot",        # a colour, a vertex column, or an array per row
            label_prop="gene",           # labels: unique values joined by '/', else coordinates
            backend="networkx")          # engine for the spring layout
ax.get_title()
# -> 'Skeleton: chr1:0-12,000\n3 nodes (4 CREs merged at 1,500 bp), 2 edges'
ax.figure.savefig("region.pdf")
```

Numeric `vertex_color` columns are mapped through `cmap` (default
`'viridis'`); string columns get one colour per value and a legend. Edge
weights and size values of merged CREs are summed. `linked_only=True`
(default) skips CREs without edges; a region with nothing to draw raises.
`ax=` draws into an existing Axes.

---

## Saving and loading
{: .sec-purple }

```python
A.save("A_dir")                       # vertices.parquet (Loci + vp), edges.parquet (src, tgt, ep), meta.json
B = gb.Architecture.load("A_dir")
B.loci.equals(A.loci), list(B) == list(A)
# -> (True, True)
```

Everything round-trips: vertex and edge columns, the name, the edge order.
`pickle` works too.

---

## Set operations and copies
{: .sec-purple }

```python
union     = A | A.chrom("chr1")    # edge union; edge columns from the first operand where both have the edge
intersect = A & A.chrom("chr1")    # common edges (a view; vertex columns are dropped)
copy      = A.copy()               # independent arrays, vertex columns included
```

Both operands must be built on the same Loci object.

---

## End-to-end recipe
{: .sec-purple }

```python
import genomeblocks as gb

cre   = gb.Loci.make("cre.bed")
genes = gb.Genes.make("genes.gtf", promoter_r=1000)

A = (gb.Architecture.make(cre, "loops.bedpe", r=2500)
       .add_mcool("hic.mcool", resolution=5000)
       .normalize()
       .prune()
       .annotate(genes)
       .strength())

res = A.prime_hubs()
print(f"{len(res['prime_genes'])} prime genes, {len(res['hub_uids'])} hub CREs")

A.save("architecture")                      # parquet tables
A.draw("chr8:127.5-128.5 Mb", layout="genomic", vertex_size_by="strength", vertex_color="annot")
```

For the viewers that take an Architecture — `gb.browser`, `gb.igv_html` and
`gb.View` — see [Browser]({{ '/guide/browser/' | relative_url }}).
