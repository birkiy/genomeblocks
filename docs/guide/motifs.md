---
title: Motifs
parent: User Guide
layout: default
nav_order: 8
---

# Motifs
{: .no_toc }

Scan TF binding motifs (JASPAR / MEME / custom PSSM) across a `Loci` set using the fast [lightmotif](https://github.com/althonos/lightmotif) backend (Martin Larralde, EMBL) — see [Credits](../credits#lightmotif) for the citation.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Preparing the genome

```python
from genomeblocks import make_genome

genome = make_genome("hg38.fa.gz")
# → {"chr1": "ACGT...", "chr2": "ACGT...", ...}
```

Any FASTA (bgzipped or plain) works; a small helper parses it into an in-memory `{chrom: str}` dict. For very large genomes and one-off scans, consider slicing to just the chromosomes present in your `Loci`.

{: .warning }
> **`lightmotif`'s `jaspar` format is the raw 4-line count format**, *not* the
> bracketed JASPAR-2016 layout (`A [ … ]`). If your PWM library is in the
> bracketed form, convert it (one PWM = a header line then four space-separated
> count rows in `A C G T` order) or load it with `motif_format='jaspar16'`.
> A ready-to-use, lightmotif-format JASPAR **H14CORE** file lives at
> `…/motif-db/H14CORE_jaspar_format_lightmotif.txt`.

---

## Scanning motifs

```python
from genomeblocks import scan_motifs

counts = scan_motifs(
    loci,
    genome,
    motif_path="JASPAR2024_CORE.jaspar",
    motif_format="jaspar",    # or 'meme'
    r=250,                    # ±r bp window around each Locus center
    threshold=13.0,           # log-odds threshold
    norm=True,                # normalize hit count by motif length
    verbose=True,
)
# → {motif_name: hit_count_per_bp or raw_count}
```

Also attached as a method on `Loci`:

```python
counts = loci.scan_motifs(genome, motif_path="motifs.jaspar")
```

---

## What the scan does

For every motif in the file:

1. Build a PSSM (`counts.normalize(0.1).log_odds()`).
2. For every `Locus`:
   - Extract a ±`r` bp window around `.center` from `genome`.
   - Stripe the sequence (required by `lightmotif` for vectorized scanning).
   - Count all positions with a log-odds score above `threshold`.
3. If `norm=True`, divide by the motif width (so long motifs don't dominate).

Returns a `{motif_name: score}` dict that you can turn into a DataFrame:

```python
import pandas as pd
pd.Series(counts).sort_values(ascending=False).head(20)
```

---

## Per-locus × per-motif matrix

When you need a full `(n_loci × n_motifs)` matrix rather than aggregate counts,
use `scan_motifs_matrix` — it extracts each window once and distributes the
motifs across a process pool:

```python
mat = loci.scan_motifs_matrix(
    genome,
    motif_path="H14CORE_jaspar_format_lightmotif.txt",
    r=250, threshold=13.0, norm=True,
    workers=None,             # None → cpu_count() // 2
)
# → pandas.DataFrame, rows = locus uids, columns = motif names
```

This is the building block for the differential and enrichment helpers below.

---

## Differential enrichment between two sets

Given two matrices (e.g. condition-A vs condition-B CREs), `compare_motifs`
runs a per-motif Mann-Whitney U test plus a log2 fold change of the means, with
Benjamini-Hochberg FDR:

```python
from genomeblocks.motifs import compare_motifs, compare_motifs_to_ref

mat_a = a_cre.scan_motifs_matrix(genome, motif_path)
mat_b = b_cre.scan_motifs_matrix(genome, motif_path)

diff = compare_motifs(mat_a, mat_b, pseudo=0.1)     # Factor, mean_A, mean_B, LFC, U, p, p_adj
diff.head()
```

To ask "which motifs are enriched above a background pool?", use
`compare_motifs_to_ref` (accepts a single matrix or a `dict` of groups sharing
one reference). `bootstrap_enrichment` gives the resampled point-estimate
variant.

---

## Masking an anchor motif

To ask "which co-factors enrich *independent of* CTCF (or any anchor)?", mask
every anchor match before scanning:

```python
masked = loci.scan_motifs_matrix_masked(
    genome, motif_path,
    anchors=["CTCF"],        # case-insensitive substring against motif names
    window=10,               # ±bp replaced with random bases around each hit
    seed=0,
)
```

Anchor motifs are excluded from the output by default. **Mask the reference the
same way** before a differential test, or query sets look depleted for reasons
unrelated to biology.

---

## Motif archetypes (clustering → consensus PWM)

Collapse a redundant PWM library into consensus archetypes (Sandelin-Wasserman
similarity → hierarchical clustering → per-cluster consensus), then render
sequence logos:

```python
from genomeblocks.motifs import build_archetypes, write_meme
from genomeblocks.motifs_draw import plot_archetypes, plot_dendrogram

res = build_archetypes(motif_path, cutoff=0.3, workers=None)
plot_archetypes(res["archetypes"], members=res["members"])   # logo grid
write_meme(res["archetypes"], "archetypes.meme")             # feed back into scanning
```

`archetype_from_names(motif_path, names=[...])` builds a single consensus from a
named subset (e.g. the top hits of a differential test).

---

## Tips

- `threshold=13.0` is a reasonable default for JASPAR (log-odds in bits). Pick an empirical value by scanning a random-shuffle control.
- Motif counts scale roughly linearly with `r`; use the same `r` you plan to use in every downstream analysis to stay comparable.
- `Locus.sequence(genome, r=...)` does the windowed slice — if you already have the index built, reuse it rather than reopening the FASTA.
