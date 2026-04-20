---
title: Architecture
parent: User Guide
layout: default
nav_order: 4
---

# Architecture
{: .no_toc }

Chromatin-contact networks in a few lines: build from HiChIP/loop files, overlay Hi-C from mcool, normalize to O/E, annotate vertices with gene context, and mine hubs and focus genes.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is

`Architecture` is a thin subclass of `graph_tool.Graph`. Each vertex represents a CRE (or loop anchor) and is keyed by a `Locus` UID. Edges represent contacts, with built-in edge properties:

| Edge property | Meaning |
|---|---|
| `ep.w` | raw contact weight (loop support or Hi-C sum) |
| `ep.n` | O/E-normalized weight (power-law expectation) |
| `ep.d` | genomic distance between endpoints |

Because it **is** a `graph_tool.Graph`, every graph-tool algorithm (centrality, SBM, layouts, community detection) works out of the box.

---

## Three ways to build a graph

### `Architecture.make(loci, bedpe, r=...)`

Assumes you already have a CRE catalogue. Each BEDPE loop anchor is mapped to any CREs within `±r` bp of its midpoint; one edge per (CRE₁, CRE₂) pair.

```python
from genomeblocks import Architecture, Loci

cre = Loci.make("cre.bed")
arch = Architecture.make(cre, "HiChIP_loops.bedpe", r=2500, dmax=1e9)
```

### `Architecture.make_spread(source_loci, bedpe, hops=...)`

Start from source anchors (e.g. gene promoters) and discover the network by hopping through loops. You do **not** need a pre-defined CRE list — vertices are loop-anchor centers.

```python
promoters = Loci(list(genes.annot["prom"]))
arch = Architecture.make_spread(
    source_loci=promoters,
    bedpe="loops.bedpe",
    hops=2,              # 1 = direct, 2 = neighbors of neighbors
    r=2500,
    directed=True,       # regulatory flow
)
```

At `hops=2` the network grows: promoters → direct contacts → contacts of contacts.

### `Architecture.make_clique(loci)`

A fully connected graph over a `Loci` set — useful as a null model for contact enrichment tests, and as a lightweight container when edges come from somewhere else (e.g. co-expression).

```python
clique = Architecture.make_clique(cre, name="null")
```

---

## Adding Hi-C weights

```python
arch = arch.add_mcool(cre, "cohesin.mcool", resolution=5000, name="w")
```

Reads the `.mcool` at the requested resolution, assigns each CRE to its nearest bin, sums `pixels` for every bin pair an edge spans, and distributes the sum evenly across all edges sharing that bin pair. The result lands in `ep.w` (or the `name` you pass).

{: .tip }
> Call `add_mcool(..., name="...")` multiple times to stack weights from different Hi-C experiments on the same graph.

---

## Normalizing to O/E

```python
arch.normalize(cre, source="w", name="n")
```

Fits a power law `w ≈ C · d^-α` over all edges with both `ep.d > 0` and `ep.w > 0`, then divides the raw weight by the fitted expectation. Observed-over-Expected lands in `ep.n`.

Typical alpha values from Hi-C: `0.8–1.2` intra-TAD, `1.5–2.0` inter-TAD. `arch.normalize()` prints the fitted `α` and `C`.

---

## Annotating vertices

```python
arch.annotate(cre, genes, verbose=True)
```

Populates, per vertex:

| Property | Value |
|---|---|
| `vp.annot` | region class (`Promoter-TSS`, `5UTR`, `Intronic`, `Intergenic`, …) |
| `vp.gene` | single nearest gene (back-compat convenience) |
| `vp.genes` | `;`-joined list of **all** candidate gene names with overlapping promoters |
| `vp.transcripts` | `,`-joined list of candidate transcript IDs (alt-promoter aware) |

Multi-candidate entries matter for bidirectional promoters (TP53 / WRAP53) and dense TSS clusters. Downstream `focus()` and `prime_hubs()` aggregate edge weight across candidates rather than forcing an arbitrary pick.

---

## Finding hubs and focus genes

### Node strength (`aggregate`)

```python
arch.aggregate(key="n", name="agg")
# vp.agg   = sum of ep.n per vertex
# vp.nagg  = normalized so it sums to 1
```

### Elbow (`elbow`)

```python
cutoff, sorted_uids = arch.elbow("agg", transform="double_exp")
hub_uids = sorted_uids[:cutoff]
```

Sorts by `vp.agg` descending and applies the Kneedle algorithm (with a double-exponential transform that amplifies curvature) to auto-pick the hub cutoff.

### Focus genes per source (`focus`)

For each source CRE, tally edge weight to each candidate promoter gene, then take the per-source elbow over the sorted gene-weight vector.

```python
result = arch.focus(source_uids=enhancer_hubs, key="n")
result["genes"]    # union of focus genes
result["records"]  # per-source gene_weights, focus_genes, n_focus
```

### The one-liner: `prime_hubs`

Runs the full pipeline end-to-end:

```python
result = arch.prime_hubs(key="n")
result["prime_genes"]    # hub-promoter genes ∪ focus genes of hub enhancers
result["promoter_genes"] # hub promoter CREs → their genes
result["focus_genes"]    # hub enhancer CREs → their focus genes
result["hub_uids"]       # all hub CREs
result["cutoff"]         # elbow index
```

---

## Subsetting

```python
# By custom predicate
sub = arch.subgraph(filter_func=lambda v: arch.vp.agg[v] > 0.5)

# By vertex property value
sub = arch.subgraph(vp_name="gene", vp_values=["MYC", "TP53"])

# By explicit UIDs
sub = arch.subgraph(uids=["chr8:127000-128000(.)", ...])

# Deep copy
copy = arch.copy()
```

All subgraphs preserve every vertex and edge property.

---

## Set operations on graphs

```python
union     = arch_a | arch_b   # vertex + edge union
intersect = arch_a & arch_b   # common vertices + common edges
```

Edge properties in the output are copied from the first operand.

---

## Drawing a region

```python
fig, ax = plt.subplots(figsize=(14, 10))

arch.draw(
    loci=cre,
    region=("chr8", 127_000_000, 130_000_000),
    merge_distance=1000,          # collapse loci within 1 kb into single nodes
    vertex_size_by="agg",
    edge_width_by="n",
    vertex_color="annot",         # color by vp.annot
    label_prop="gene",
    layout="spring",
    ax=ax,
)
fig.savefig("myc_locus.pdf")
```

Uses graph-tool's SFDP / Kamada-Kawai / circular layouts.

---

## Serialization

```python
import pickle
pickle.dump(arch, open("arch.pkl", "wb"))
arch2 = pickle.load(open("arch.pkl", "rb"))
```

`Architecture` overrides `__getstate__` / `__setstate__` to extract vertex and edge property maps into pickleable dicts and rebuild the graph on load. All properties round-trip.

---

## End-to-end recipe

```python
from genomeblocks import Loci, Genes, Architecture

cre   = Loci.make("cre.bed")
genes = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)

arch = (Architecture.make(cre, "HiChIP.bedpe", r=2500)
                    .add_mcool(cre, "cohesin.mcool", resolution=5000, name="w")
                    .normalize(cre, source="w", name="n")
                    .annotate(cre, genes)
                    .aggregate(key="n", name="agg"))

result = arch.prime_hubs(key="n")
print(f"{len(result['prime_genes'])} prime genes, "
      f"{len(result['hub_uids'])} hub CREs")
```
