---
title: motifs
parent: API Reference
layout: default
nav_order: 9
---

# `genomeblocks.motifs`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `make_genome(path) -> dict[str, str]`

Parse a FASTA (optionally `.gz`) into an in-memory `{chrom: sequence}` dict.

```python
genome = make_genome("hg38.fa.gz")
```

---

## `scan_motifs(loci, genome, motif_path, motif_format='jaspar', r=250, threshold=13.0, norm=True, verbose=True) -> dict[str, float]`

Scan each `Locus` window (±`r` bp around `.center`) for every motif in `motif_path`.

| Arg | Meaning |
|---|---|
| `loci` | `Loci` to scan. |
| `genome` | `dict[chrom → seq]` or path to FASTA (auto-loaded). |
| `motif_path` | Path to a motif collection readable by `lightmotif.load`. |
| `motif_format` | `'jaspar'`, `'meme'`, etc. |
| `r` | Half-window around `Locus.center`. |
| `threshold` | Log-odds score threshold. |
| `norm` | Divide counts by motif width. |
| `verbose` | Show `tqdm` progress over motifs. |

Also attached as `Loci.scan_motifs(genome, motif_path, ...)`.

Returns a dict `{motif_name: count_or_normalized}`.

{: .note }
> `lightmotif`'s `'jaspar'` format is the **raw 4-line count** layout, not the
> bracketed JASPAR-2016 form — use `motif_format='jaspar16'` for the latter.

---

## `scan_motifs_matrix(loci, genome, motif_path, *, r=250, threshold=13.0, norm=True, workers=None, ...) -> DataFrame`

Full `(n_loci × n_motifs)` count matrix (rows = locus uids, columns = motif
names). Windows are extracted once; motifs are scanned across a process pool.
Attached as `Loci.scan_motifs_matrix(...)`.

## `scan_motifs_matrix_masked(loci, genome, motif_path, anchors, *, window=10, skip_anchors=True, seed=None, ...) -> DataFrame`

Same matrix, but every match to an `anchors` motif (case-insensitive substring)
is masked with random bases first. Attached as
`Loci.scan_motifs_matrix_masked(...)`.

---

## Differential / enrichment

| Function | Purpose |
|---|---|
| `compare_motifs(mat_a, mat_b, *, pseudo=0.1, alternative='two-sided')` | Per-motif Mann-Whitney U + LFC + BH FDR between two matrices. |
| `compare_motifs_to_ref(query, ref, ...)` | One or more query matrices vs a reference pool; `query` may be a `dict` of groups. |
| `bootstrap_enrichment(groups, ref, *, boot=100, sample=500, ...)` | Resampled mean motif counts → LFC point estimates. |

---

## Archetypes (`genomeblocks.motifs`)

| Function | Purpose |
|---|---|
| `pwm_distance_matrix(motif_path, ...) -> (D, names, pwms)` | Pairwise Sandelin-Wasserman distance matrix. |
| `cluster_motifs(D, *, cutoff=0.3, linkage_method='average') -> (labels, Z)` | Hierarchical clustering. |
| `archetype(pwms, ...) -> pfm` | Consensus PFM for one cluster. |
| `archetype_from_names(motif_path, names, ...) -> dict` | Consensus from a named subset. |
| `build_archetypes(motif_path, *, cutoff=0.3, ...) -> dict` | End-to-end load → distance → cluster → consensus per cluster. |
| `write_meme(archetypes, path)` | Write archetypes as MEME so they can be re-scanned. |

## Logos (`genomeblocks.motifs_draw`)

`plot_archetype(pfm, ...)`, `plot_archetypes(archetypes, ...)`,
`plot_cluster_members(...)`, `plot_dendrogram(Z, ...)` — sequence-logo helpers
(lazy `logomaker` import).
