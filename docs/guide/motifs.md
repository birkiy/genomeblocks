---
title: Motifs
parent: User Guide
layout: default
nav_order: 8
---

# Motifs
{: .no_toc }

Scan TF binding motifs (JASPAR / MEME / custom PSSM) across a `Loci` set using the fast `lightmotif` backend.
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

## Per-locus hits

`scan_motifs` currently returns *aggregated* counts per motif. If you need per-locus per-motif matrices, iterate explicitly:

```python
import lightmotif, numpy as np

motifs = list(lightmotif.load("motifs.jaspar", format="jaspar"))
n, m = len(loci), len(motifs)
mat = np.zeros((n, m), dtype=np.int32)

for mi, mdl in enumerate(motifs):
    pssm = mdl.counts.normalize(0.1).log_odds()
    for li, l in enumerate(loci):
        seq = l.sequence(genome, r=250).upper()
        if len(seq) != 500: continue
        seq = lightmotif.stripe(seq)
        mat[li, mi] = sum(1 for _ in lightmotif.scan(pssm, seq, threshold=13.0))
```

---

## Tips

- `threshold=13.0` is a reasonable default for JASPAR (log-odds in bits). Pick an empirical value by scanning a random-shuffle control.
- Motif counts scale roughly linearly with `r`; use the same `r` you plan to use in every downstream analysis to stay comparable.
- `Locus.sequence(genome, r=...)` does the windowed slice — if you already have the index built, reuse it rather than reopening the FASTA.
