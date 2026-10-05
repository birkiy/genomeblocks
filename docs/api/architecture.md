---
title: architecture
parent: API Reference
layout: default
nav_order: 5
---

# `genomeblocks.architecture`
{: .no_toc }

A CRE interaction graph stored as two tables: the vertices are the rows of a
`Loci` (vertex `i` = row `i`), the edges are a sorted `src` / `tgt` table with
edge columns in `A.ep` (`w` weight, `d` distance, `n` observed / expected, ...)
and vertex columns in `A.vp` (`annot`, `gene`, `strength`, ...). See the
[Architecture guide]({{ '/guide/architecture/' | relative_url }}) for the
workflow and [Design: architecture]({{ '/design/architecture/' | relative_url }})
for the layout of the edge table.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Architecture`
{: .sec-navy }

```python
import genomeblocks as gb
from genomeblocks import Architecture

Architecture(loci, src=(), tgt=(), *, ep=None, vp=None, name="Architecture")
```

The constructor takes vertex rows `src[k]` / `tgt[k]` and dicts of edge and
vertex columns; `loci` must be genome-sorted (what `Loci.make` returns, or
call `.sort()`). The classmethods below are the usual way in.

Edges are kept in canonical order: `src < tgt`, cis edges grouped by
chromosome in genome order, every trans edge in one final block. That order
is what `blocks`, `chrom()`, `cis` and `trans` slice.

```python
cre   = gb.Loci.make("peaks.bed", keep=True)
genes = gb.Genes.make("genes.gtf")
A = Architecture.make(cre, "loops.bedpe", r=100)
# -> [INFO] 4 loops | 4 mapped (100.0%) | loci=6, links=4 (1 trans)
A
# -> Architecture(name='Skeleton', loci=6, links=4 [3 cis · 1 trans], edge_props=[w], vertex_props=[])
list(A)                      # the edges as (src row, tgt row)
# -> [(0, 2), (1, 4), (5, 6), (0, 5)]
A.blocks
# -> {'chr1': (0, 2), 'chr2': (2, 3), 'trans': (3, 4)}
```

As a table an `Architecture` *is* its edge table: `shape`, `columns`,
`head()`, `describe()`, `to_pandas()` / `to_polars()` / `to_arrow()` and the
Arrow / dataframe-interchange / narwhals protocols all describe the edges
(see [Tables out](#tables-out)). `len(A)` counts the vertices with at least
one edge.

---

## Construction
{: .sec-green }

| Method | Signature | One line |
|---|---|---|
| `Architecture.make` | `make(loci, bedpe, *, name="Skeleton", r=2500, dmax=1e9, trans=True, verbose=True, backend=None)` | Edges between every pair of CREs within `r` bp of the two anchor midpoints of each loop. |
| `Architecture.from_edges` | `from_edges(loci, src, tgt, *, name="Architecture", vp=None, **ep)` | From vertex rows; keyword arrays become edge columns. |
| `Architecture.from_frame` | `from_frame(loci, df, *, src="src", tgt="tgt", name="Architecture")` | From an edge table (pandas / polars / arrow) with rows in `src` / `tgt` or uids in `uid1` / `uid2`. |
| `Architecture.from_scipy` | `from_scipy(loci, M, *, weight="w", name="Architecture")` | From a square sparse matrix over the loci rows (upper triangle read). |
| `Architecture.load` | `load(path, *, genome=None)` | Read back what [`save`](#copy-save-and-load) wrote. |

### `make`

`loci` is anything [`as_loci`]({{ '/api/interop/' | relative_url }}) takes;
`bedpe` is a BEDPE path, a [`Pairs`]({{ '/api/bedpe/' | relative_url }}) or a
frame with BEDPE columns. Cis loops longer than `dmax` are skipped (`dmax`
does not apply to trans loops); `trans=False` drops trans loops. `backend`
picks the interval engine for the anchor-to-CRE mapping (any
[intervals backend]({{ '/backends/' | relative_url }}) gives the same edges).
The new graph has one edge column, `w`, set to 0 — `add_mcool` fills it.

```python
A = Architecture.make(cre, gb.Pairs.make("loops.bedpe"), r=100, dmax=5000, trans=False, verbose=False)
A.n_links, A.n_trans
# -> (2, 0)
```

If `loci` is not in genome order `make` prints a warning and works on a
sorted copy, exposed as `A.loci` — align vertex columns to `A.loci`, not to
the input.

### `from_edges`

```python
B = Architecture.from_edges(cre, [0, 0, 1], [2, 5, 4], w=[5.0, 1.0, 3.0])
B
# -> Architecture(name='Architecture', loci=5, links=3 [2 cis · 1 trans], edge_props=[w], vertex_props=[])
list(B), B.ep["w"].tolist()          # re-sorted into canonical order, columns follow
# -> ([(0, 2), (1, 4), (0, 5)], [5.0, 3.0, 1.0])
```

### `from_frame` and `from_scipy`

`from_frame` reads what `edges_frame()` writes: every numeric column other
than `src`, `tgt`, `uid1`, `uid2`, `chrom1`, `chrom2`, `cis` becomes an edge
column. A table with neither rows nor uids raises a `ValueError` naming the
accepted columns.

```python
df = A.edges_frame()
Architecture.from_frame(cre, df).ep.keys()
# -> dict_keys(['w', 'd', 'n'])
Architecture.from_frame(cre, df[["uid1", "uid2", "w"]]).n_links
# -> 4
Architecture.from_scipy(cre, A.to_scipy("w")).ep["w"].tolist()
# -> [6.0, 3.0, 4.0, 9.0]
Architecture.from_frame(cre, df[["w"]])
# -> ValueError: an edge table needs vertex rows in 'src' / 'tgt' or uids in 'uid1' / 'uid2'
#    (what edges_frame() writes); got columns ['w']
```

---

## Sizes and blocks
{: .sec-green }

| Name | Type | Meaning |
|---|---|---|
| `A.loci` | `Loci` | The vertex table (vertex `i` = row `i`). |
| `A.src`, `A.tgt` | `int32` arrays | Edge endpoints, `src < tgt`, canonical order. |
| `A.ep`, `A.vp` | `Props` | Edge / vertex columns; a dict with attribute access (`A.ep.w` is `A.ep["w"]`). |
| `A.name` | `str` | Shown in `repr` and plot titles. |
| `A.n_links` | `int` | Number of edges. |
| `A.n_loci`, `len(A)` | `int` | Vertices with at least one edge. |
| `A.n_trans` | `int` | Number of trans edges. |
| `A.degree` | `int` array, one per Loci row | Edges per vertex (0 for unlinked rows). |
| `A.is_cis` | `bool` array, one per edge | Both ends on one chromosome. |
| `A.blocks` | `dict` | `{chrom: (lo, hi)}` edge ranges of each cis block plus `'trans'`. |
| `A.block_counts()` | `DataFrame` | chrom x chrom edge counts (cis on the diagonal). |
| `A.fit` | `dict` | `{'alpha', 'C'}` of the last `normalize()` power-law fit. |

```python
A.n_links, A.n_trans, A.n_loci, A.degree.tolist()
# -> (4, 1, 6, [2, 1, 1, 0, 1, 2, 1])
A.is_cis.tolist()
# -> [True, True, True, False]
A.block_counts()
# ->       chr1  chr2
#    chr1     2     1
#    chr2     1     1
```

---

## Views
{: .sec-green }

A view shares the `Loci` and the vertex columns with its parent and selects
edges; slices (`chrom`, `cis`, `trans`) copy nothing.

| Method | Signature | One line |
|---|---|---|
| `A.chrom` | `chrom(chrom)` | Cis edges of one chromosome. |
| `A.cis` | property | All cis edges. |
| `A.trans` | property | All trans edges. |
| `A.chroms` | `chroms()` | Iterate `(chrom, view)` over the cis blocks. |
| `A.region` | `region(region, start=None, end=None, *, both=True)` | Edges inside a region; `both=False` keeps edges with one end inside (trans partners too). |
| `A.subgraph` | `subgraph(rows=None, *, mask=None, vp=None, values=None, name=None)` | Edges whose two ends are both selected: by rows, a boolean mask, or `vp[vp] in values`. |

```python
A.chrom("chr1").n_links, A.cis.n_links, A.trans.n_links
# -> (2, 3, 1)
for chrom, view in A.chroms():
    print(chrom, view.n_links)
# -> chr1 2
#    chr2 1
A.region("chr1:0-6,000")
# -> Architecture(name='Skeleton:chr1:0-6,000', loci=2, links=1 [1 cis · 0 trans], ...)
A.region("chr1", 0, 6000, both=False)
# -> Architecture(name='Skeleton:chr1:0-6000', loci=5, links=3 [2 cis · 1 trans], ...)
A.subgraph(rows=[0, 2]).n_links
# -> 1
A.subgraph(vp="annot", values="Promoter-TSS")       # after annotate()
# -> Architecture(name='Skeleton_sub', loci=2, links=1 [1 cis · 0 trans], ...)
```

Region strings accept thousands separators and units: `'chr1:1,000-2,000'`,
`'chr8:127.7-128.1 Mb'`, `'chr2:5kb-12kb'`; a `(chrom, start, end)` tuple or
a `Locus` works too.

---

## Lookups
{: .sec-green }

A vertex is named by its row number, its uid (`'chr1:900-1100(+)'`) or a
`Locus` read from the table; all three work wherever a vertex is expected.

| Method | Signature | One line |
|---|---|---|
| `A.neighbor_rows` | `neighbor_rows(x)` | `(partner rows, edge ids)` of one CRE, cis and trans, as arrays. |
| `A.neighbors` | `neighbors(x)` | The same as a DataFrame: `row`, `uid`, `chrom`, `cis` + every edge column. |
| `A[uid]` | `__getitem__` | `{neighbour uid: w}`. |
| `A[uid1, uid2]` | `__getitem__` | That edge's columns as a dict, or `None` when there is no edge. |
| `uid in A` | `__contains__` | True when the uid is a vertex with at least one edge. |
| `A.near_rows` | `near_rows(key, start=None, end=None, *, r=0, mode="overlap", linked=True)` | Row numbers of the CREs near a window (`key` is a Locus, a region string or a chromosome with `start` / `end`); `r` widens the window on each side. |
| `A.near` | `near(key, start=None, end=None, *, r=0, mode="overlap")` | The same CREs as a `Loci`. |

`mode='overlap'` keeps CREs that intersect the window, `'center'` only those
whose midpoint is in it.

```python
A.neighbor_rows(0)
# -> (array([2, 5], dtype=int32), array([0, 3]))
A.neighbors(0)
# ->    row                uid chrom    cis    w       d         n
#    0    2  chr1:4900-5100(+)  chr1   True  6.0  4000.0  1.124534
#    1    5    chr2:500-600(+)  chr2  False  9.0     inf  1.000000
A[cre.uid[0]]
# -> {'chr1:4900-5100(+)': 6.0, 'chr2:500-600(+)': 9.0}
A[cre.uid[0], cre.uid[2]]
# -> {'w': 6.0, 'd': 4000.0, 'n': 1.1245336516005273}
A[cre.uid[0], cre.uid[1]]
# -> None
cre.uid[0] in A, cre.uid[3] in A
# -> (True, False)
A.near_rows("chr1:0-3,000")
# -> array([0, 1])
A.near(genes["GENE_A"].tss, r=2000)           # TSS ± 2 kb
# -> Loci(n=2, chroms=1, cols=[name, score])
```

---

## Weights: `add_mcool`
{: .sec-navy }

```python
A.add_mcool(mcool, *, resolution=None, name="w", verbose=True) -> Architecture
```

Edge weight = the Hi-C count of the (bin, bin) pixel that holds the two
CREs, shared equally among the edges that fall in the same pixel. Pixels are
read one chromosome block at a time, so memory stays bounded and trans
pixels are found in the same pass. Pass `resolution=` for a multi-resolution
`.mcool`; a single-resolution `.cool` takes none. Both mistakes raise a
`ValueError` that lists the resolutions available.

```python
A.add_mcool("hic.cool")
# -> [INFO] Set distributed weights for 4/4 edges from cooler. [w]
A.ep["w"].tolist()
# -> [6.0, 3.0, 4.0, 9.0]
A.add_mcool("hic.cool", resolution=1000)
# -> ValueError: hic.cool is a single-resolution cooler: call add_mcool without resolution=
```

---

## Normalize, prune, annotate
{: .sec-navy }

| Method | Signature | One line |
|---|---|---|
| `A.normalize` | `normalize(*, source="w", name="n", verbose=True)` | Observed / expected: a power law on distance for cis edges, the mean trans weight for trans edges; writes `ep.d` and `ep[name]`, stores the fit in `A.fit`. |
| `A.prune` | `prune(*, dist_prop="d", verbose=True)` | Drop zero-distance (co-located) cis edges; trans edges (`d = inf`) stay. |
| `A.annotate` | `annotate(genes, *, key="n", name="gene", verbose=True)` | `vp.annot` = region label per CRE; `vp[name]` = the nearest-TSS gene for promoter CREs, else the gene of the highest-`key` promoter neighbour (cis or trans). |

Every one of these returns `self`, so they chain.

```python
A.normalize()
# -> [INFO] Power-law fit: alpha=0.802, C=4.143e+03 on 3 cis edges; 1 trans edges use the mean trans weight → ep.n
list(A.ep), A.ep["d"].tolist()
# -> (['w', 'd', 'n'], [4000.0, 9000.0, 4500.0, inf])
A.prune()
# -> [INFO] prune: removed 0 zero-distance edges → 4 edges.
A.annotate(genes)
# -> [INFO] Annotated 6 loci: 3 promoter CREs | 2/3 non-promoter CREs assigned to a top-'n' promoter gene → vp.gene.
A.vp["annot"].tolist()
# -> ['Promoter-TSS', 'Promoter-TSS', '3UTR', 'Promoter-TSS', 'Promoter-TSS', 'Intergenic', 'Intergenic']
A.vp["gene"].tolist()
# -> ['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_A', '']
```

{: .warning }
> `prune` reads `ep.d` and `annotate` ranks promoter contacts by `ep.n`, so
> both need `normalize()` first; they raise a `ValueError` naming the missing
> column otherwise. A fit needs at least three cis edges at two distinct
> distances; with fewer, `A.fit` holds `nan` and the cis `n` values are 0.

Vertex columns are one value per `Loci` row, so unlinked CREs get a label
too (`annot`) and an empty gene (`''`). `genes` is a
[`Genes`]({{ '/api/genes/' | relative_url }}) table; a TSS is the gene's
5'-most base, `[t, t+1)`.

---

## Hubs: strength, elbow, prime_hubs, support
{: .sec-navy }

| Method | Signature | One line |
|---|---|---|
| `A.strength` | `strength(key="n", name="strength", *, verbose=True)` | `vp[name]` = sum of incident `ep[key]` per vertex (node strength). |
| `A.elbow` | `elbow(key, *, verbose=True)` | Slope-1 knee on the sorted `vp[key]` curve; returns `(cutoff, uids sorted by value, descending)`. |
| `A.prime_hubs` | `prime_hubs(key="n", gene="gene", *, verbose=True)` | Hubs above the knee of node strength and the genes they carry. |
| `A.support` | `support(genes, *, r=5000, mode="overlap", uids=True, rows=False, linked=True, backend=None)` | `{gene_name: [uid, ...]}` of the linked CREs within TSS ± `r` of each gene. |

`prime_hubs` needs `vp.annot` and `vp[gene]` (run `annotate` first); it
computes the strength column when it is missing. Its dict has the keys
`prime_genes`, `promoter_genes`, `enhancer_genes`, `hub_uids`, `cutoff`,
`promoter_uids`, `enhancer_uids`. The knee is the same slope-1 rule as
[`se.knee`]({{ '/api/atlas/' | relative_url }}#knee), so super-enhancers and
hubs are cut the same way.

```python
A.strength()
# -> [INFO] Summed ep.n → vp.strength (node strength).
A.elbow("strength")
# -> [INFO] Slope-1 on vp.strength (value≈1.12 at cutoff): cutoff at 3/6 (50.0%)
#    (3, ['chr1:900-1100(+)', 'chr2:500-600(+)', 'chr1:4900-5100(+)', 'chr1:1900-2100(-)', ...])
hubs = A.prime_hubs(verbose=False)
hubs["hub_uids"], hubs["prime_genes"]
# -> (['chr1:900-1100(+)', 'chr2:500-600(+)', 'chr1:4900-5100(+)'], {'GENE_A'})
```

`support` computes the TSS ± `r` window for every gene at once;
`rows=True` returns CRE row arrays instead of uids (what array work wants),
`mode='center'` keeps only CREs whose midpoint is in the window,
`linked=False` includes CREs without edges.

```python
A.support(genes, r=5000)
# -> {'GENE_A': ['chr1:900-1100(+)', 'chr1:1900-2100(-)', 'chr1:4900-5100(+)'], 'GENE_B': ['chr1:10900-11100(-)']}
A.support(genes, r=5000, rows=True)
# -> {'GENE_A': array([0, 1, 2]), 'GENE_B': array([4])}
```

---

## Graph backends
{: .sec-navy }

Graph algorithms run through the graph backend family (`graph-tool` when
installed, else `scipy`; `igraph` and `networkx` on request — see
[Backends]({{ '/backends/' | relative_url }})). Results are numpy arrays
aligned to the `Loci` rows and are identical across engines. A requested
engine that is not installed raises `ImportError` with the install command.

| Method | Signature | One line |
|---|---|---|
| `A.graph` | `graph(vprops=(), *, backend=None)` | The graph as the engine's own object (vertex `i` = row `i`, every edge, edge columns as properties, `vprops` copied too). |
| `A.components` | `components(name="component", *, backend=None)` | Connected component per vertex (`-1` for unlinked rows), numbered by first row; stored as `vp[name]`. |
| `A.pagerank` | `pagerank(weight=None, *, damping=0.85, name="pagerank", backend=None)` | PageRank per vertex (`weight` = an edge column or `None`); stored as `vp[name]`. |
| `A.to_graph_tool` | `to_graph_tool(vprops=())` | `graph(vprops, backend="graph-tool")`. |
| `A.to_networkx` | `to_networkx(vprops=(), *, loci_cols=True)` | `networkx.Graph`; `loci_cols` adds `chrom` / `start` / `end` node attributes. |
| `A.to_igraph` | `to_igraph(vprops=(), *, loci_cols=True)` | `igraph.Graph`, same attributes. |
| `A.to_scipy` | `to_scipy(weight="w")` | Symmetric CSR adjacency over the Loci rows (`weight=None`: 1 per edge). |
| `A.to_anndata` | `to_anndata(*, weights=None)` | AnnData with the vertices as `obs` and one sparse adjacency per edge column in `obsp` (the layout scanpy's graph tools read). |

The graph-tool graph is built once and cached; its edge columns are re-synced
on every `graph()` call. `to_networkx` and `to_scipy` sum parallel edges.

```python
A.components().tolist()
# -> [0, 1, 0, -1, 1, 0, 0]
A.pagerank("w").round(4).tolist()
# -> [0.2438, 0.1626, 0.1073, 0.0244, 0.1626, 0.2179, 0.0814]
A.components(backend="igraph").tolist() == A.components(backend="scipy").tolist()
# -> True
A.components(backend="graph-tool")          # on a machine without it
# -> ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool

g = A.to_networkx()
g.number_of_edges(), g[0][2], g.nodes[0]
# -> (4, {'w': 6.0, 'd': 4000.0, 'n': 1.1245336516005273}, {'chrom': 'chr1', 'start': 900, 'end': 1100})
ig = A.to_igraph(vprops=["strength"])
ig.vs.attributes(), ig.es.attributes()
# -> (['strength', 'chrom', 'start', 'end'], ['w', 'd', 'n'])
A.to_scipy("w")[0, 2]
# -> 6.0
A.to_anndata()
# -> AnnData object with n_obs × n_vars = 7 × 0
#        obs: 'chrom', 'start', 'end', 'strand', 'name', 'score', 'annot', 'gene', 'strength', 'component', 'pagerank'
#        uns: 'genomeblocks'
#        obsp: 'w', 'd', 'n'
```

---

## Tables out
{: .sec-green }

| Method | Returns | One line |
|---|---|---|
| `A.edges_frame()` | `DataFrame` | One row per edge: `src`, `tgt`, `uid1`, `uid2`, `chrom1`, `chrom2`, `cis` + every edge column. |
| `A.to_pandas()` | `DataFrame` | Same as `edges_frame()`. |
| `A.to_polars()` | `polars.DataFrame` | The edge table through Arrow. |
| `A.to_arrow()` | `pyarrow.Table` | The edge table. |
| `A.columns` | `list` | Column names of the edge table. |
| `A.shape` | `(n_links, n_columns)` | |
| `A.head(n=5)`, `A.tail(n=5)` | `DataFrame` | First / last edges. |
| `A.describe()` (alias `summary()`) | `DataFrame` | Name, loci, linked loci, edges (cis / trans / blocks), edge and vertex columns. |
| `A.vertices_frame(linked=True)` | `DataFrame` | One row per vertex with links (`linked=False`: every row): coordinates, uid and every vertex column, indexed by row. |
| `A.to_frame()` | `DataFrame` | Vertex table for linked vertices: `uid` + vertex columns. |
| `A.block_counts()` | `DataFrame` | chrom x chrom edge counts. |
| `iter(A)` | | Yields `(src row, tgt row)` per edge. |
| `A._repr_html_()` | | Notebook summary: the cis blocks and the trans block. |

The Arrow C stream, `__dataframe__` and narwhals protocols come from the same
`to_arrow()`, so `pl.DataFrame(A)`, `duckdb.sql("select * from A")`,
seaborn, plotly and altair take an `Architecture` as it is.

```python
A.shape, A.columns
# -> ((4, 10), ['src', 'tgt', 'uid1', 'uid2', 'chrom1', 'chrom2', 'cis', 'w', 'd', 'n'])
A.to_polars()
# -> shape: (4, 10)
#    ┌─────┬─────┬───────────────────┬─────────────────────┬───┬───────┬─────┬────────┬──────────┐
#    │ src ┆ tgt ┆ uid1              ┆ uid2                ┆ … ┆ cis   ┆ w   ┆ d      ┆ n        │
#    ╞═════╪═════╪═══════════════════╪═════════════════════╪═══╪═══════╪═════╪════════╪══════════╡
#    │ 0   ┆ 2   ┆ chr1:900-1100(+)  ┆ chr1:4900-5100(+)   ┆ … ┆ true  ┆ 6.0 ┆ 4000.0 ┆ 1.124534 │
#    │ 1   ┆ 4   ┆ chr1:1900-2100(-) ┆ chr1:10900-11100(-) ┆ … ┆ true  ┆ 3.0 ┆ 9000.0 ┆ 1.077749 │
#    │ 5   ┆ 6   ┆ chr2:500-600(+)   ┆ chr2:5000-5100(+)   ┆ … ┆ true  ┆ 4.0 ┆ 4500.0 ┆ 0.823993 │
#    │ 0   ┆ 5   ┆ chr1:900-1100(+)  ┆ chr2:500-600(+)     ┆ … ┆ false ┆ 9.0 ┆ inf    ┆ 1.0      │
#    └─────┴─────┴───────────────────┴─────────────────────┴───┴───────┴─────┴────────┴──────────┘
A.describe()
# ->                                                                value
#    name                                                        Skeleton
#    loci (vertices)                                                    6
#    loci with links                                                    6
#    edges                                                              4
#    cis edges                                                          3
#    trans edges                                                        1
#    cis blocks (chromosomes)                                           2
#    edge columns                                                 w, d, n
#    vertex columns            annot, gene, strength, component, pagerank
A.to_frame()
# ->                    uid         annot    gene  strength  component  pagerank
#    0     chr1:900-1100(+)  Promoter-TSS  GENE_A  2.124534          0  0.243820
#    1    chr1:1900-2100(-)  Promoter-TSS  GENE_A  1.077749          1  0.162602
#    2    chr1:4900-5100(+)          3UTR  GENE_A  1.124534          0  0.107289
#    ...
```

---

## Set operations
{: .sec-green }

Both graphs must be built on the *same* `Loci` object.

| Expression | Result |
|---|---|
| `A \| B` | Union of the edge sets; an edge column present in one graph only is 0 on the other's edges. |
| `A & B` | The edges of `A` that are also in `B`, with `A`'s edge columns (a view; vertex columns are not carried). |

```python
U = A | B
list(U), U.ep["w"].tolist()
# -> ([(0, 2), (1, 4), (5, 6), (0, 5)], [6.0, 3.0, 4.0, 9.0])
list(A & B)
# -> [(0, 2), (1, 4), (0, 5)]
```

---

## Copy, save and load
{: .sec-green }

| Method | Signature | One line |
|---|---|---|
| `A.copy()` | | Deep copy of edges, edge columns and vertex columns (the `Loci` is shared). |
| `A.save` | `save(path)` | Writes a directory: `vertices.parquet` (Loci + `vp:` columns), `edges.parquet` (`src`, `tgt`, `ep:` columns), `meta.json`. |
| `Architecture.load` | `load(path, *, genome=None)` | Reads it back; `genome` re-codes the chromosomes into an existing `Genome`. |

```python
A.save("arch")
sorted(os.listdir("arch"))
# -> ['edges.parquet', 'meta.json', 'vertices.parquet']
back = Architecture.load("arch")
list(back) == list(A), back.loci.equals(A.loci)
# -> (True, True)
```

---

## Drawing
{: .sec-purple }

`A.draw(region, **kw)` calls `genomeblocks.architecture_draw.draw` with the
graph. Plain matplotlib on top of the tables: the CREs of a region become
nodes (merged when closer than `merge_distance`), the edges among them
become lines, and the positions come from the graph backend's layout.
Nothing here needs graph-tool.

```python
from genomeblocks.architecture_draw import draw

draw(A, region, *,
     layout="spring",                 # 'spring' | 'circular' | 'genomic'
     backend=None,                    # graph backend for the spring layout
     merge_distance=None,             # merge CREs whose centres are within this many bp
     vertex_size_by=None,             # vertex column name or array (one value per Loci row)
     edge_width_by="w",               # edge column that scales line widths
     vertex_size_range=(40, 400),
     edge_width_range=(0.6, 4.0),
     vertex_color=None,               # a colour, a vertex column name, or an array
     cmap="viridis",                  # for numeric vertex colours
     edge_color="#b8b8b8",
     figsize=(8, 6),
     show_labels=True,
     label_prop="gene",               # vertex column used for node labels
     font_size=8,
     linked_only=True,                # skip CREs without edges
     seed=0,
     title=True,
     ax=None) -> matplotlib.axes.Axes
```

| Argument | Meaning |
|---|---|
| `region` | `'chr1:1,000-2,000'`, `'chr8:127.5-128.5 Mb'`, `(chrom, start, end)` or a `Locus`. |
| `layout` | `'spring'`: the graph backend's force-directed layout (spectral for scipy); `'circular'`; `'genomic'`: nodes on a line at their genomic position, edges as arcs — a browser-like view. |
| `merge_distance` | CREs whose centres are within this many bp become one node; edge weights and size values are summed. |
| `vertex_size_by` | Default: the number of CREs in the node. |
| `vertex_color` | Numbers are mapped through `cmap`; strings get one colour per value and a legend. |
| `label_prop` | Unique non-empty values of that vertex column, joined by `/`; nodes without one are labelled with their coordinates. |

Raises `ValueError` when the region holds no (linked) CREs.

```python
ax = A.draw("chr1:0-12 kb")
ax.get_title()
# -> 'Skeleton: chr1:0-12,000\n4 nodes, 2 edges'
ax = A.draw(("chr1", 0, 12_000), layout="genomic", merge_distance=1500,
            vertex_size_by="strength", vertex_color="annot")
ax.get_title()
# -> 'Skeleton: chr1:0-12,000\n3 nodes (4 CREs merged at 1,500 bp), 2 edges'
A.draw("chr1:0-12kb", layout="circular", backend="networkx", edge_width_by="n", label_prop="gene")
A.draw("chr2:7,000-8,000")
# -> ValueError: no CREs with links in chr2:7,000-8,000
```

For a region view with signal tracks and gene models see
[`browser`]({{ '/api/browser/' | relative_url }}); for a shareable
interactive page of the whole graph see
[`View`]({{ '/api/browser/' | relative_url }}#view) and
[`igv_html`]({{ '/api/browser/' | relative_url }}#igv_html).
