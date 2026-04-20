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
