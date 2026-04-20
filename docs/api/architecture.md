---
title: architecture
parent: API Reference
layout: default
nav_order: 5
---

# `genomeblocks.architecture`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Architecture(graph_tool.Graph)`

Undirected by default (some builders accept `directed=True`). Every vertex has `vp.uid`; every edge has `ep.w` (weight), `ep.n` (normalized), `ep.d` (distance).

### Constructor

```python
Architecture(name: str | None = None)
```

### Properties & shortcuts

| Name | Description |
|---|---|
| `n_loci` | `num_vertices()`. |
| `n_links` | `num_edges()`. |
| `index` | `{uid → Vertex}`. |

### Lookup

```python
arch[uid]                    # {neighbor_uid: ep.w[edge]}
arch[(uid1, uid2)]           # {ep_name: value} for the connecting edge
uid in arch                  # bool
len(arch)                    # vertex count
```

---

### Factory methods

```python
Architecture.make(loci, bedpe, *, name="Skeleton",
                  r=2500, dmax=1e9, verbose=True) -> Architecture
# Build from BEDPE loops; r = ±radius around anchor midpoints for CRE mapping.

Architecture.make_clique(loci, *, name="Clique", verbose=True) -> Architecture
# Fully connected graph over loci; ep.d pre-filled.

Architecture.make_spread(source_loci, bedpe, *, name="Spread",
                         r=2500, dmax=1e9, hops=1,
                         directed=False, verbose=True) -> Architecture
# Spreading constructor. Discovers loop-anchor centers reachable from source loci in N hops.
```

---

### Weight assignment & normalization

```python
arch.add_mcool(loci, mcool, *, resolution=None,
               name="w", verbose=True) -> Architecture
# Assign Hi-C pixel sums; distributes across overlapping edges per bin pair.

arch.normalize(loci, *, source="w", name="n",
               verbose=True) -> Architecture
# Fit w ≈ C·d^-α power law; write O/E to ep.n (or `name`).
```

---

### Annotation

```python
arch.annotate(loci, genes, *, verbose=True) -> Architecture
# Populates vp.annot, vp.gene, vp.genes (multi), vp.transcripts (multi).
```

---

### Network mining

```python
arch.aggregate(key="n", name="agg", *, verbose=True) -> Architecture
# Node strength: vp.agg = sum of incident ep[key]; vp.nagg = normalized.

arch.cluster(key="n", name="lc", *, verbose=True) -> Architecture
# Weighted local clustering coefficient per vertex.

arch.elbow(key, *, transform="double_exp", verbose=True) -> (cutoff, sorted_uids)
# Kneedle cutoff over vp[key] descending.

arch.focus(sources, key="n", *, verbose=True) -> dict
# dict keys: 'genes' (set), 'records' (list[{source, n_candidates, gene_weights, focus_genes}])

arch.prime_hubs(key="n", *, verbose=True) -> dict
# dict keys: 'prime_genes', 'promoter_genes', 'focus_genes',
#            'hub_uids', 'cutoff', 'promoter_uids', 'enhancer_uids'
```

---

### Subsetting & set operations

```python
arch.copy() -> Architecture                           # deep-copy all props
arch.subgraph(filter_func=None, vp_name=None, vp_values=None,
              uids=None, name=None) -> Architecture   # one-of filter modes

arch | other                                           # vertex+edge union
arch & other                                           # common vertices + common edges
```

---

### Visualization

```python
arch.draw(loci, region, *,
          merge_distance=None,     # collapse nearby loci into single nodes
          vertex_size_by=None,
          edge_width_by='w',
          vertex_size_range=(10, 50),
          edge_width_range=(1, 10),
          vertex_color=None,       # str hex | vp_name | None
          edge_color='#CCCCCC',
          figsize=(12, 8),
          layout='spring',         # 'spring' | 'circular' | 'kamada_kawai'
          show_labels=True,
          label_prop='gene',
          ax=None,
          **kwargs)
```

Returns the matplotlib axis.

---

### Serialization

```python
arch.to_frame() -> pandas.DataFrame    # one row per vertex, one col per vp
pickle.dumps(arch)                      # full property round-trip via __getstate__/__setstate__
```
