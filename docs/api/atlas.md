---
title: atlas
parent: API Reference
layout: default
nav_order: 9
---

# `genomeblocks.atlas`
{: .no_toc }

GIGGLE-style enrichment over a sparse `bin × track` index. See the
[Atlas guide](../guide/atlas) for narrative usage.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Atlas.make(paths, *, chromsizes, ...) -> Atlas`

```python
Atlas.make(
    paths: str | Iterable[str],       # glob, directory, or list of BED-family paths
    *,
    chromsizes: str | dict | object,  # .chrom.sizes path, {chrom: len}, or a cooler.Cooler
    names: Sequence[str] | None = None,
    bin_size: int = 1000,
    workers: int | None = None,       # None → cpu_count() - 1
    verbose: bool = True,
    meta: str | pd.DataFrame | None = None,
    meta_columns: Sequence[str] | None = None,
    meta_id_col: str | None = None,
    meta_sep: str = "\t",
    name_pattern: str | None = None,  # re.search on the basename → track id
) -> Atlas
```

Reads each BED/narrowPeak/broadPeak (optionally `.gz`) file into its set of
covered bins and assembles a single CSR `M` of shape `(n_bins, n_tracks)`.
Accepts `.bed`, `.narrowPeak`, `.broadPeak` and their `.gz` variants.

---

## `atlas.search(query, *, ref=None, alternative="two-sided") -> pd.DataFrame`

Per-track Fisher 2×2 enrichment, sorted by `giggle_score` descending.

- `ref=None` — c/d cells come from the genome null.
- `ref=<Loci>` — c/d cells come from the reference set ("more enriched in
  query than in ref").
- `alternative` — `"two-sided"` (default), `"greater"`, or `"less"`.

Result columns: `name`, `n_query_bins`, `n_track_bins`, `track_n_peaks`,
`overlaps`, `log2_odds`, `p`, `giggle_score`, plus any attached metadata.

Fluent alias: `Loci.enrich(atlas, ref=None, **kw)`.

---

## `atlas.bootstrap(query, *, n=10, ...) -> pd.DataFrame`

```python
atlas.bootstrap(
    query: Loci | dict[str, Loci],
    *,
    n: int = 10,
    pool: Loci | None = None,         # None → position-shuffle null
    sample: int | None = None,        # subsample query + pool per iteration
    replace: bool = False,
    keep_chrom: bool = True,          # shuffle within the original chromosome
    seed: int | None = None,
    verbose: bool = True,
) -> pd.DataFrame
```

Empirical Monte-Carlo enrichment. Result columns: `name`, `observed`,
`expected`, `obs_std`, `null_std`, `log2fc`, `z`, `p_emp`, `track_n_bins`
(plus `group` when `query` is a dict, and any attached metadata).

Fluent alias: `Loci.enrich_mc(atlas, n=10, **kw)`.

---

## `atlas.attach_meta(meta, *, id_col=None, columns=None, sep="\t") -> Atlas`

Attach per-track metadata aligned to `track_names`. `meta` is a TSV/CSV path or
a `DataFrame`; `id_col` selects the join column (defaults to the first);
`columns` supplies headers for a header-less file. Returns `self`.

---

## Persistence

| Method | Purpose |
|---|---|
| `atlas.save(path)` | Write a compressed `.npz` (no pickles — `allow_pickle=False`). |
| `Atlas.load(path)` | Reconstruct an `Atlas`, including metadata. |

---

## Attributes

| Attribute | Meaning |
|---|---|
| `bin_size` | Genome tiling resolution (bp). |
| `n_bins` | Total number of genome bins across all chroms. |
| `track_names` | List of track ids (order matches `M`'s columns). |
| `track_n_peaks` / `track_n_bins` | Per-track raw peak count / covered-bin count. |
| `M` | `scipy.sparse.csr_matrix` of shape `(n_bins, n_tracks)`. |
| `meta` | Optional `DataFrame` aligned to `track_names`. |
