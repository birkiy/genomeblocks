---
title: Atlas
parent: User Guide
layout: default
nav_order: 8
---

# Atlas
{: .no_toc }

GIGGLE-style enrichment of a region set against a large collection of BED
tracks (e.g. all of ChIP-Atlas). One sparse `bin × track` index turns each
query into a single sparse matrix-vector product, so thousands of tracks are
scored in one vectorised pass.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The idea

Tile the genome at a fixed resolution (default 1 kb) and store one bit per
`(bin, track)` cell in a single CSR matrix. A query set is reduced to its
covered bins; `M[query_bins].sum(0)` then returns the per-track shared-bin
count for **every track at once**. From those overlap counts the `Atlas`
computes a Fisher 2×2 enrichment (GIGGLE score) or an empirical Monte-Carlo
null.

At ChIP-Atlas scale (~25k tracks, 1 kb bins) the index is a few GB resident and
each query is sub-second.

> Layer et al., *GIGGLE: a search engine for large-scale integrated genome
> analysis.* Nat. Methods 15, 123–126 (2018).

---

## Building an index

```python
from genomeblocks import Atlas

atlas = Atlas.make(
    "chipatlas/*.bed.gz",              # glob, directory, or list of paths
    chromsizes="hg38.chrom.sizes",    # path, dict, or a cooler.Cooler
    bin_size=1000,
    workers=None,                     # None → cpu_count() - 1 (one BED per task)
    name_pattern=r"^[^.]+",           # keep only the SRX accession from the filename
    meta="chipatlas_meta.tsv",        # optional per-track metadata (see below)
    meta_id_col="srx",
)
print(atlas)                          # Atlas(tracks=..., bins=..., bin_size=1000, nnz=...)
```

`paths` accepts a glob string, a directory, or an explicit list. `chromsizes`
can be a UCSC `.chrom.sizes` path, a `{chrom: length}` dict, or anything with a
`.chromsizes` mapping (like a `cooler.Cooler`). Intervals on chromosomes not in
`chromsizes` are silently dropped.

Persist and reload with the compressed `.npz` format (never pickles — safe to
share):

```python
atlas.save("chipatlas_hg38_1kb.npz")
atlas = Atlas.load("chipatlas_hg38_1kb.npz")
```

---

## Enrichment search — `atlas.search(query)`

Fisher 2×2 of each track against the genome null, in bin units:

```python
cre = Loci.make("my_peaks.narrowPeak")
res = atlas.search(cre)                       # DataFrame sorted by giggle_score
res.head()[["name", "overlaps", "log2_odds", "p", "giggle_score"]]
```

Columns: `overlaps` (query bins hit), `log2_odds`, `p` (exact hypergeometric),
and `giggle_score = -log10(p) · log2(OR)` — negative for depletion. If you
attached metadata, its columns are merged in.

Pass a **reference set** to ask "more enriched in query than in ref?" — the
c/d cells of the 2×2 come from `ref` instead of the genome:

```python
res = atlas.search(query=up_cre, ref=all_cre)
```

The same call is available fluently on any `Loci`:

```python
res = cre.enrich(atlas)                # == atlas.search(cre)
res = up.enrich(atlas, ref=all_cre)
```

---

## Monte-Carlo null — `atlas.bootstrap(query, n=...)`

When you want an empirical null instead of the analytic Fisher test:

```python
res = atlas.bootstrap(cre, n=1000, seed=0)    # position-shuffle null (chrom-aware)
res[["name", "observed", "expected", "z", "p_emp"]].head()
```

Three interacting knobs control the null:

| Argument | Null model |
|---|---|
| `pool=None` (default) | Per-interval position shuffle; `keep_chrom=True` preserves each interval's chromosome. |
| `pool=<Loci>` | Draw the null from a curated CRE universe (LOLA / regioneR style) — controls for the universe's own bias. |
| `sample=<int>` | Per iteration, subsample every query group **and** the pool to `sample` regions; the pool draw is shared across groups so comparisons are paired. |

`query` may be a single `Loci` **or** a `dict[str, Loci]`; with a dict the
result is long-format with a `group` column, so several query sets share one
null. Fluent form: `cre.enrich_mc(atlas, n=1000)`.

---

## Attaching metadata

Track names alone are rarely enough — you usually want antigen / cell-type
labels alongside the scores. Attach a table keyed on the track id:

```python
atlas.attach_meta(
    "chipatlas_meta.tsv",
    id_col="srx",                     # column that joins against track_names
    columns=["srx", "antigen", "cell"],   # supply headers for a header-less file
)
```

Tracks with no metadata row get `NaN`; metadata rows for unknown tracks are
dropped. Every `search` / `bootstrap` result then carries the metadata columns.

---

## See also

- [API → Atlas]({{ '/api/atlas/' | relative_url }}) for full signatures.
- [Enrichment walkthrough]({{ '/walkthrough/enrichment/' | relative_url }}) for a worked ChIP-Atlas example.
