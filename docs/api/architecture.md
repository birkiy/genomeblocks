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

### Construction

```python
Architecture.make(loci, bedpe, *, name="Skeleton",
                  r=2500, dmax=1e9, verbose=True) -> Architecture
# Build from BEDPE loops (parsed via bedpe.read_bedpe); r = ±radius around
# anchor midpoints for CRE mapping, dmax drops far-apart loops.
```

---

### Weight assignment & normalization

```python
arch.add_mcool(loci, mcool, *, resolution=None,
               name="w", verbose=True) -> Architecture
# Assign Hi-C pixel sums; distributes across overlapping edges per bin pair.

arch.normalize(loci, *, source="w", name="n",
               verbose=True) -> Architecture
# Fit w ≈ C·d^-α power law on cis edges; write O/E to ep.n (or `name`),
# distances to ep.d. Trans edges: d = inf, expectation = mean trans weight.

arch.prune(*, dist_prop="d", verbose=True) -> Architecture
# Drop zero-distance (co-located) edges. Call once, after normalize.
```

---

### Annotation

```python
arch.annotate(loci, genes, *, key="n", name="gene", verbose=True) -> Architecture
# vp.annot = region class; vp[name] = nearest gene (promoters) or the gene of the
# top-`key`-weight promoter neighbour (other CREs). Run normalize first.
```

---

### Hub genes

```python
arch.strength(key="n", name="strength", *, verbose=True) -> Architecture
# Node strength: vp[name] = sum of incident ep[key] (no normalization).

arch.elbow(key, *, verbose=True) -> (cutoff, sorted_uids)
# Slope-1 knee over vp[key] descending; cutoff = number of hubs.

arch.prime_hubs(key="n", gene="gene", *, verbose=True) -> dict
# dict keys: 'prime_genes', 'promoter_genes', 'enhancer_genes',
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

Drawing lives in `genomeblocks.architecture_draw` (separate module, matplotlib):

```python
from genomeblocks.architecture_draw import draw

draw(arch, loci, region, *,
     merge_distance=None,     # collapse nearby loci into single nodes
     vertex_size_by=None,
     edge_width_by='w',
     vertex_color=None,       # str hex | vp_name | PropertyMap | None
     edge_color='#CCCCCC',
     figsize=(12, 8),
     layout='spring',         # 'spring' | 'circular' | 'kamada_kawai'
     show_labels=True,
     label_prop='gene',
     ax=None,
     **kwargs) -> matplotlib axis
```

---

### Serialization

```python
arch.to_frame() -> pandas.DataFrame    # one row per vertex, one col per vp
pickle.dumps(arch)                      # full property round-trip via __getstate__/__setstate__
```
