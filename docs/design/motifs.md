---
title: Motifs
parent: Design
layout: default
nav_order: 4
---

# Motifs
{: .no_toc }

`scan_motifs_matrix` counts position-weight-matrix hits in a fixed window
around every region, for every motif in a library, and returns a
regions × motifs table that the enrichment helpers compare between sets.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## Scanning many short windows
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/motifs-block.svg %}
</div><figcaption>
<strong>One scan per motif, not per window.</strong> A 500-bp window is so short that calling the SIMD scanner once per window costs more in call overhead than in scoring. The windows are joined and striped once, and each motif is scanned over the whole block in a single <code>lightmotif.scan()</code> call. Each hit is assigned to its window by binary search on the window offsets. A hit that would span two windows is dropped, so the counts equal scanning each window on its own.
</figcaption></figure>

- **Windows.** Each region contributes `sequence(genome, r)`: its centre
  ± `r`. Windows that run off a contig (shorter than `2r`) are skipped, so the
  result's index lists the regions actually scanned. Windows containing a
  character the scanner cannot encode count zero for every motif.
- **Scores.** Each motif's counts become a log-odds matrix (pseudocount 0.1,
  uniform background), and a hit is any position scoring at least
  `threshold` (13 by default). `norm=True` divides counts by motif width.
- **Parallelism.** Motifs are independent, so they are split into batches
  across a process pool. Each worker builds the striped block once in its
  initializer. Jobs with fewer than 16 motifs or 200 windows run serially,
  because process start-up would dominate.

## Masking an anchor motif
{: .sec-navy }

`scan_motifs_matrix_masked(..., anchors=["CTCF"])` first finds every hit of the
anchor motifs and overwrites ± `window` bp around each with random bases. It
then scans the masked windows. Co-factors are counted without the anchor's
own sites (and the anchor motifs are left out of the result). The masked bases
are random A/C/G/T, seeded by `seed=`.

## From matrices to enrichment
{: .sec-navy }

- `bootstrap_enrichment(groups, ref)` resamples rows of each matrix and reports
  per-motif log2 fold change against a reference set (and the difference when
  exactly two groups are given).
- `compare_motifs(mat_a, mat_b)` runs a per-motif Mann-Whitney U test on the
  counts of two matrices and reports the log2 fold change of their means.
- `compare_motifs_to_ref` does the same against a reference matrix.

## Motif archetypes
{: .sec-navy }

Large libraries repeat themselves (dozens of near-identical FOX or ETS
matrices). `pwm_distance_matrix` scores every pair of PWMs, `cluster_motifs`
cuts a hierarchical clustering at `cutoff`, and `build_archetypes` aligns each
cluster into one information-weighted archetype PWM. `write_meme` exports them
for MEME tools. `motifs_draw` draws archetype logos, cluster members and the
dendrogram (`logomaker` is optional).

Measured numbers are on the [Benchmarks]({{ '/benchmarks/#motifs' | relative_url }}) page.
