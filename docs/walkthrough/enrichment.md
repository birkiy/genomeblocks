---
title: "4. Motif & ChIP-Atlas enrichment"
parent: "Example: AR & FOXA1"
layout: default
nav_order: 4
---

# Motif & ChIP-Atlas enrichment
{: .no_toc }

Two complementary "what's different?" analyses: **sequence motifs** enriched in
each AR set, and **published TF datasets** (ChIP-Atlas) enriched in each set.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

# Motif enrichment

## The idea

If AR+F and AR−F are bound by different cofactors, their underlying DNA should
carry different transcription-factor **motifs**. We scan a JASPAR motif library
across each set and ask which motifs are enriched relative to a background —
here all accessible chromatin, so we control for "open chromatin" sequence
composition and find what is specific to each AR class.

## Scanning motifs into a matrix

Pass the FASTA path straight to the scanner — the windows are read through the
indexed FASTA reader, so the genome is never loaded whole — and the motif file:

```python
from genomeblocks.motifs import scan_motifs_matrix, bootstrap_enrichment

pool = accessible                           # background = ATAC union

mat_pf   = scan_motifs_matrix(ARpF, GENOME_FA, MOTIF_DB, r=250, workers=8)
mat_mf   = scan_motifs_matrix(ARmF, GENOME_FA, MOTIF_DB, r=250, workers=8)
mat_pool = scan_motifs_matrix(pool, GENOME_FA, MOTIF_DB, r=250, workers=8)
```

`scan_motifs_matrix` returns a `pandas.DataFrame` shaped **(loci x motifs)**:
rows are the loci in their order (index = `uid`), columns are motif names, and
each cell is the motif's hit count in that locus's **centre ± `r`** window
(normalised by motif width by default). `MOTIF_DB` is a JASPAR-format file
(`format='jaspar'`, the raw four-row count layout, is the default; pass
`format='jaspar16'` for the bracketed layout, or `'meme'`). The scan runs on
lightmotif; `backend='moods'` or `'biopython'` give the same hits.

{: .note }
> Scanning the full accessible pool can be the slow step — it is the largest set.
> `workers` parallelises the scan over motifs; if it is still heavy, subsample
> `pool` to a few thousand regions for the background. If you would rather
> hold the genome in memory, `gb.read_fasta(GENOME_FA)` returns a
> `{chrom: sequence}` dict that every scanner accepts in place of the path.

## Bootstrapped log-fold-change

`bootstrap_enrichment` compares the per-motif mean hit rate of each group to the
reference pool, with bootstrap resampling for stability:

```python
enr = bootstrap_enrichment({"AR+F": mat_pf, "AR-F": mat_mf},
                           ref=mat_pool, boot=100, sample=500, seed=0)
```

It returns one row per motif with:

- `mean_ref` — mean rate in the background pool;
- `mean_<g>`, `LFC_<g>` — mean rate and `log2((mean_g + pseudo)/(mean_ref + pseudo))`
  for each group;
- `LFC` — when exactly two groups are given, the **difference** `LFC_AR+F −
  LFC_AR−F`, i.e. the motif's preference for one set over the other.

`boot` / `sample` are the number of bootstrap iterations and the per-iteration
subsample size; `seed` makes it reproducible.

Sorting by `LFC` puts AR+F-preferential motifs on top and AR−F-preferential ones
at the bottom:

```python
ranked = enr.sort_values("LFC", ascending=False)
ranked.head(10)[["Factor", "LFC", "LFC_AR+F", "LFC_AR-F"]]   # top in AR+F
ranked.tail(10)[["Factor", "LFC", "LFC_AR+F", "LFC_AR-F"]]   # top in AR-F
```

Biologically you expect **forkhead (FOX) motifs** to surface among the
AR+F-preferential set — the sequence signature of the FOXA1 pioneering that
defines those sites.

---

# ChIP-Atlas enrichment

## The idea

Motifs tell you about *sequence*; ChIP-Atlas tells you about *measured binding*.
[ChIP-Atlas](https://chip-atlas.org) aggregates tens of thousands of public
ChIP-seq / ATAC experiments. `genomeblocks`' `Atlas` indexes them so you can ask:
*which published datasets does my region set overlap more than expected?* Using
AR−F as the reference for AR+F (and vice versa) gives a **differential** answer.

## Building the index

`Atlas.make` builds a sparse **bin x track** matrix: every chromosome is tiled
into `bin_size` (1 kb) bins, and each ChIP-Atlas BED becomes a column marking the
bins it covers. A query then reduces to one sparse mat-vec — sub-second per
query after the one-time build.

```python
from genomeblocks import Atlas

atlas = Atlas.make(
    GIGGLE_DIR, chromsizes=CHROMSIZES,
    meta=GIGGLE_META,
    name_pattern=r"([^.]+)",                                # SRX...20.bed.gz -> SRX...
    meta_columns=["id", "antigen", "class", "cell_line"],   # meta TSV is header-less
    meta_id_col="id",
)
atlas.save("chipatlas_hg38_1kb.npz")       # reuse across sessions: Atlas.load(...)
```

`chromsizes` takes a `.chrom.sizes` path as here, a `{chrom: length}` dict, a
`gb.Genome.from_fasta(GENOME_FA)`, or a cooler. The result is a table of its
tracks (`atlas.head()`, `atlas.to_pandas()`), with the metadata columns joined.

Two arguments are essential for the metadata to line up — and are a common
first-run trip-up:

{: .warning }
> **`name_pattern`** extracts the track id from each filename. ChIP-Atlas files
> look like `SRX23002840.20.bed.gz`; `r"([^.]+)"` keeps just the accession
> (`SRX23002840`) so it matches the metadata id. Without it the `.20` suffix
> never joins and every metadata column comes back `NaN`.
>
> **`meta_columns` + `meta_id_col`** handle a **header-less** metadata TSV:
> passing `meta_columns` tells the reader there is no header row and names the
> columns; `meta_id_col` picks the one to join on. Omit these and the first data
> row is mistaken for the header.

## Differential search

`atlas.search(query, ref=...)` runs a per-track Fisher 2x2 of the query against
the reference (in bin units) and returns a `DataFrame` sorted by `giggle_score`:

```python
enr_pf = atlas.search(ARpF, ref=ARmF)    # tracks enriched in AR+F over AR-F
enr_mf = atlas.search(ARmF, ref=ARpF)    # tracks enriched in AR-F over AR+F

cols = ["name", "antigen", "class", "cell_line", "overlaps", "log2_odds", "giggle_score"]
enr_pf.head(15)[cols]
```

Columns: `name` (track id), `overlaps` (query bins hit), `log2_odds` and `p`
(the Fisher test), `giggle_score` (the signed significance used for ranking), and
the joined `antigen / class / cell_line` metadata. The query and the reference
are anything `as_loci` takes, and `ARpF.enrich(atlas, ref=ARmF)` is the same
call as a method.

## Reading the result

This is the validation step. On this data the top AR+F (FOXA1-dependent) tracks
resolve to **FOXA1 in prostate / LNCaP** — independent experiments confirming
those AR sites are FOXA1 territory — while the top AR−F (FOXA1-independent)
tracks shift to other factors such as **NR3C1**, a different nuclear receptor.
The `antigen`/`cell_line` columns (now populated thanks to the metadata args)
are what make that readable.

{: .note }
> `search(query)` with no `ref` tests against a genome-wide null instead, and
> `bootstrap(query, n)` gives a shuffled-position null — use those when you want
> "enriched vs the genome" rather than "enriched vs the other set." See the
> [Atlas guide]({{ '/guide/atlas/' | relative_url }}).

Next: **[Genome browser →]({{ '/walkthrough/browser/' | relative_url }})**
