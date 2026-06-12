---
title: Architecture
parent: User Guide
layout: default
nav_order: 4
---

# Architecture
{: .no_toc }

Chromatin-contact networks in a few lines: build from HiChIP/loop files, overlay
Hi-C from mcool, normalize to O/E, annotate vertices with gene context, and mine
hub genes.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is

`Architecture` is a thin subclass of
[`graph_tool.Graph`](https://graph-tool.skewed.de/) (Peixoto 2014 — see
[Credits](../credits#graph-tool)). Each vertex is a CRE keyed by its `Locus`
UID; edges are contacts, with built-in edge properties:

| Edge property | Meaning |
|---|---|
| `ep.w` | raw contact weight (loop support or Hi-C sum) |
| `ep.n` | O/E-normalized weight (power-law expectation) |
| `ep.d` | genomic distance between endpoints |

Because it **is** a `graph_tool.Graph`, every graph-tool algorithm (centrality,
SBM, layouts, community detection) works out of the box.

{: .note }
> The core pipeline is intentionally small: **make → add_mcool → normalize →
> annotate → strength → prime_hubs**. Drawing lives in a separate
> `architecture_draw` module so the graph object stays dependency-light.

---

## Build from loops

`Architecture.make` assumes you already have a CRE catalogue. Each BEDPE loop
anchor is mapped to any CREs within `±r` bp of its midpoint; one edge per
(CRE₁, CRE₂) pair. BEDPE parsing is delegated to the package's single reader
(`bedpe.read_bedpe`) — you just pass a path.

```python
from genomeblocks import Architecture, Loci

cre  = Loci.make("cre.bed")
arch = Architecture.make(cre, "HiChIP_loops.bedpe", r=2500, dmax=1e9)
```

| arg | meaning |
|---|---|
| `r` | radius (bp) around each loop-anchor midpoint to catch CREs |
| `dmax` | drop loops whose anchors are farther apart than this |

---

## Adding Hi-C weights

```python
arch = arch.add_mcool(cre, "cohesin.mcool", resolution=5000, name="w")
```

Reads the `.mcool` at the requested resolution, assigns each CRE to its nearest
bin, sums `pixels` for every bin pair an edge spans, and distributes that sum
across the edges sharing the bin pair. The result lands in `ep.w` (or the `name`
you pass — call it repeatedly to stack weights from several experiments).

---

## Normalizing to O/E, and pruning

```python
arch.normalize(cre, source="w", name="n")   # power-law O/E -> ep.n, distances -> ep.d
arch.prune()                                 # drop zero-distance (co-located) edges
```

`normalize` fits a power law `w ≈ C · d^-α` over all edges with positive distance
and weight, then divides the raw weight by the fitted expectation; the
observed-over-expected ratio lands in `ep.n` and the printed `α`/`C` describe the
decay. Zero-distance edges are a power-law singularity — run `prune()` **once at
the end** to remove them (removing edges mid-pipeline desyncs graph-tool's
edge-index range).

---

## Annotating vertices

```python
arch.annotate(cre, genes, key="n")
```

Two-stage gene assignment (requires the edge weights from `normalize`, so it can
rank promoter contacts):

| Property | Value |
|---|---|
| `vp.annot` | region class (`Promoter-TSS`, `5UTR`, `3UTR`, `Exonic`, `Intronic`, `Intergenic`) |
| `vp.gene` | promoters → their nearest gene; other CREs → the gene of their **highest-`key`-weight promoter neighbour** (`''` if none) |

So every CRE is tied to a gene either by being a promoter or by its strongest
promoter contact. Pass a distinct `name=` per `key` to keep multiple assignments
side by side.

---

## Hub genes

### Node strength — `strength`

```python
arch.strength(key="n", name="strength")   # vp.strength = sum of ep.n per vertex
```

Sums the incident edge weights per vertex. No normalization is applied — divide
yourself if you want a fraction.

### Cutoff — `elbow`

```python
cutoff, sorted_uids = arch.elbow("strength")
hub_uids = sorted_uids[:cutoff]
```

Sorts vertices by `vp.strength` descending and finds the **slope-1 knee** (the
45° point on the [0,1]-normalized ascending curve, measured on a smoothed
curve). Hubs are everything above the cutoff.

### The one-liner — `prime_hubs`

```python
result = arch.prime_hubs(key="n")
result["prime_genes"]     # genes of all hub CREs
result["promoter_genes"]  # genes of hub CREs that are promoters
result["enhancer_genes"]  # genes of hub CREs that are not promoters
result["hub_uids"]        # all hub CRE uids
result["cutoff"]          # number of hubs
```

`prime_hubs` runs `strength` (if needed) → `elbow` → splits the hub CREs by
whether they are promoters, collecting each hub's `vp.gene`. Requires
`vp.annot` and `vp.gene` from `annotate`.

---

## Subsetting, set operations, serialization

```python
sub  = arch.subgraph(filter_func=lambda v: arch.vp.strength[v] > 0.5)
sub  = arch.subgraph(vp_name="gene", vp_values=["MYC", "TP53"])
sub  = arch.subgraph(uids=["chr8:127000-128000(.)", ...])
copy = arch.copy()

union     = arch_a | arch_b   # vertex + edge union (edge props from the first operand)
intersect = arch_a & arch_b   # common vertices + common edges

import pickle
pickle.dump(arch, open("arch.pkl", "wb"))   # all vertex/edge properties round-trip
```

All subgraphs and copies preserve every vertex and edge property.

---

## Drawing a region

Drawing lives in `genomeblocks.architecture_draw` (kept out of the core so
`Architecture` carries no matplotlib dependency):

```python
from genomeblocks.architecture_draw import draw

ax = draw(
    arch, loci=cre,
    region=("chr8", 127_000_000, 130_000_000),
    merge_distance=1000,        # collapse loci within 1 kb into single nodes
    vertex_size_by="strength",
    edge_width_by="n",
    vertex_color="annot",       # color by a vertex property
    label_prop="gene",
    layout="spring",            # 'spring' | 'circular' | 'kamada_kawai'
)
```

---

## End-to-end recipe

```python
from genomeblocks import Loci, Genes, Architecture

cre   = Loci.make("cre.bed")
genes = Genes.make("gencode.v49.annotation.gtf", promoter_r=1000)

arch = (Architecture.make(cre, "HiChIP.bedpe", r=2500)
                    .add_mcool(cre, "cohesin.mcool", resolution=5000, name="w")
                    .normalize(cre, source="w", name="n")
                    .annotate(cre, genes)
                    .strength(key="n", name="strength"))
arch.prune()

result = arch.prime_hubs(key="n")
print(f"{len(result['prime_genes'])} prime genes, "
      f"{len(result['hub_uids'])} hub CREs")
```
