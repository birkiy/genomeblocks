---
title: Concepts
layout: default
nav_order: 4
---

# Concepts
{: .no_toc }

The mental model behind `genomeblocks`. Read this before diving into the module guides — it explains *why* the API looks the way it does.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## UIDs are the lingua franca

Every interval carries a deterministic UID:

```
"{chrom}:{start}-{end}({strand})"
# e.g. "chr1:9800-10100(.)" or "chr6:122,707,489-122,714,633(+)"
```

The UID is how:

- `Loci` maps names to positions: `loci["chr1:9800-10100(.)"]`.
- `Architecture` names graph vertices — so intersecting two graphs or subsetting by locus is just set algebra on UID strings.
- `Genes.annotations()` and `Genes.nearest_genes()` return DataFrames keyed on UIDs.

{: .tip }
> UIDs are stable across copies, so pickle a `Loci` and the downstream
> `Architecture` built from it still line up.

---

## Lazy, cached indexes

Most expensive structures are built once and memoized:

| Attribute | Purpose | Built on first access |
|---|---|---|
| `Loci.cgr` | `cgranges` interval tree for O(log n) overlap queries | yes |
| `Loci.uids` | `{uid → position}` dict for `loci[uid]` lookup | yes |
| `Genes.annot` | Pre-computed promoter / exon / UTR `Loci` unions | yes |

The pickling layer (`__getstate__` / `__setstate__`) intentionally drops these so you never serialize a stale index.

---

## Chainable pipelines

All core operations return either the same type (`Loci`, `Architecture`) or a new instance — never `None` and never "mutate in place and also return." That means you can always chain:

```python
cre = (Loci.make("peaks.bed")
           .slop(100)
           .sort()
           .merge())

arch = (Architecture.make(cre, "loops.bedpe", r=2500)
                    .add_mcool(cre, "m.mcool", resolution=5000)
                    .normalize(cre)
                    .annotate(cre, genes)
                    .strength(key="n"))
```

This mirrors the `dplyr` / pandas method-chain style and keeps intermediate state out of your namespace.

---

## Set algebra everywhere

Python operators map to set operations on both `Loci` and `Architecture`:

| Operator | `Loci` | `Architecture` |
|---|---|---|
| `a & b` | overlap intersection | vertex + edge intersection |
| `a \| b` | union (concat) | vertex + edge union |
| `a - b` | loci in `a` not overlapping `b` | — |
| `a ^ b` | symmetric difference | — |
| `a + b` | concatenation | — |

Because `&`, `-`, `+` operate on UID membership, you compose CRE sets directly —
no separate query layer:

```python
active_distal = (atac & h3k27ac) - promoters
```

---

## Lazy imports at the top level

`import genomeblocks` is deliberately cheap — public names are loaded from `__init__.py` via an `__getattr__` shim only when first referenced. Practical consequences:

- You can import `Loci`, `Genes`, `browser` without paying the graph-tool import cost.
- `Architecture` only requires graph-tool at the moment you touch it.
- No circular-import headaches between `loci`, `signal`, and `motifs` — each attaches its own methods to `Loci` on import.

See `genomeblocks/__init__.py` for the full `_EXPORTS` map.

---

## Method attachment

You'll notice that many methods are defined at module scope and then bound to the class at the bottom of the file:

```python
# in loci.py
def intersect(s, o): ...
def slop(s, n): ...

Loci.intersect = intersect
Loci.slop = slop
```

Two reasons:

1. **Avoids circular imports** — `Loci.signal` lives in `signal.py`, which imports `Loci`. Defining `signal()` as a free function and attaching it lets both modules coexist.
2. **Keeps related code together** — BedPE utilities stay in `bedpe.py`, signal utilities in `signal.py`, even though they hang off `Loci`.

---

## Graph-tool under the hood

`Architecture` **is** a `graph_tool.Graph`. Every vertex / edge property you set with `g.new_vp(...)` or `g.ep[...]` is available through the native graph-tool API — filtering, layout, centrality, community detection all just work.

`genomeblocks` layers domain semantics on top:

- `vp.uid`, `vp.annot`, `vp.gene`, `vp.genes`, `vp.transcripts` — set by `annotate()`.
- `ep.w`, `ep.n`, `ep.d` — raw weight, O/E-normalized, genomic distance.
- Custom weights via `add_mcool(name=...)`.

If you need an algorithm `genomeblocks` doesn't wrap, reach straight into `graph_tool.all as gt` and the graph is yours.

---

## Backends and fallbacks

| Feature | Preferred | Fallback |
|---|---|---|
| BigWig reads | `pybigtools` (Rust, releases GIL) | pure-Python reader in `genomeblocks/bigwig.py` |
| Overlap queries | `cgranges` | — |
| `nearest` | `pyranges` | — |
| mcool I/O | `cooler` | — |

The signal extractor auto-detects `pybigtools`; pass `backend="bigwig"` to force the pure-Python reader (useful for debugging or when compiled deps are unavailable).
