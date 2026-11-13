---
title: Landscape
parent: Benchmarks
layout: default
nav_order: 1
permalink: /benchmarks/landscape/
---

# Landscape: genomeblocks against the field
{: .no_toc }

Every building block of genomeblocks does a job other tools also do. This page
lists those tools theme by theme, runs the ones that install here on the same
data and the same task, and scores each of them the way Saelens *et al.*
(*Nat. Biotechnol.* 2019) scored trajectory-inference methods: **accuracy**,
**scalability**, **stability** and **usability**, aggregated into one overall
score. The figure follows that paper's Fig. 2: tool properties and capability
icons on the left, the aggregated scores on the right.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## The figure
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/landscape-figure.svg %}
</div><figcaption>
<strong>Tools by theme, graded per theme.</strong> Rows are tools, grouped
by the theme they serve (tools doing three or more core jobs in three or more
themes come first). Panel <strong>a</strong> gives each tool's platform,
interface, package index, first publication year, whether it saw a release or
commit in the last twelve months, and one dot per theme: its grade on that
theme's ten core jobs (see <a href="#the-grade">The grade</a>) — a filled dot
sized best / good / moderate / low, a hollow ring when the tool does some core
jobs but was not graded, a grey dot when it does none. Panel <strong>b</strong>
gives the four aspect scores and the overall score, averaged over the themes a
tool was benchmarked in; "not benchmarked here" rows show only usability.
</figcaption></figure>

Download: [SVG]({{ '/assets/images/landscape.svg' | relative_url }}) · [PNG]({{ '/assets/images/landscape.png' | relative_url }}) · [PDF]({{ '/assets/images/landscape.pdf' | relative_url }})

## The grade
{: .sec-navy #the-grade }

Each theme has ten **core jobs** (and ten auxiliary ones), chosen on
2026-10-23 from what the field's main tools are built around, before looking at
what genomeblocks does (`benchmarks/landscape/jobs.json`). A tool's grade in a
theme combines how much of the theme it covers with how well it does what it
was timed on:

- **C, coverage**: the share of the ten core jobs the tool does, among those
  whose status is known. A job counts when the tool has its own function for
  it (measured here, or documented in the inventory); a recipe composed from
  other functions does not count. genomeblocks' status comes from its source
  code, and its "partial" jobs count as no. A tool is graded only when at least
  seven of the ten statuses are known; unknown is never read as no.
- **Q, quality**: over the core jobs the tool was timed on *against at least one
  other tool* (a task only one tool ran gives no comparison and is left out),
  the mean of speed x accuracy, where speed is `1 − log10(t / t_fastest) / 3`
  and accuracy the task's agreement with its reference or simulated truth.
- **S = C x Q**, and the tier is S relative to the best S in the theme: best
  (≥ 0.9 of it), good (≥ 0.6), moderate (≥ 0.3), low.

The grade does not average across themes, and usability (a property of the
whole tool) is not part of it.

## The metrics
{: .sec-navy }

| Aspect | What is measured | How it is scored |
|---|---|---|
| **Accuracy** | Does the tool return the reference answer for the task? The reference is an independent tool or an exact computation: bedtools for interval sets and nearest distances, GenomicFeatures/TxDb for region classes, per-base bigWig means for signal matrices, an exact numpy sliding-window scan for log-odds motif hits and MEME FIMO for p-value hits (genomeblocks' own scan used to be the reference until 2026-10-23), the exact interval-level Fisher test for enrichment rankings, bedtools `pairtobed` for loop overlaps, R tximport / edgeR / DESeq2 for expression; and, for the jobs added on 2026-10-23, a simulated truth wherever one exists — simulated reads for BAM counting and coverage, simulated fold changes, dispersions and batch effects for the DE family, planted motifs, a Hi-C map with planted TADs, loops, compartments and bias. | 1 when the normalised answer is identical; otherwise the agreement — Pearson *r* for matrices and vectors, Spearman *ρ* for rankings, the fraction of rows agreeing for labels and distances — averaged over the tasks the tool ran in the theme. |
| **Scalability** | Wall time of the task at its largest size (1 M peaks, 100 k regions, 10 k windows × 200 bins, 1,000 windows × 100 motifs, 500 tracks, 50 k loops, 60 k genes × 24 samples), medians of repeated runs after a warm-up, one tool at a time. | `1 − log10(t / t_fastest) / 3`, floored at 0: each decade behind the fastest tool costs a third of the score; averaged over the tool's tasks. The slopes across sizes are in the tables. |
| **Stability** | The same data presented differently — rows shuffled, chromosomes in another order, samples or genes reordered — must give the identical answer (after undoing the permutation). | The fraction of perturbation tests passed. Tools that were not perturbed show "n/a". |
| **Usability** | A transparent ten-item checklist from the inventory, in the spirit of the paper's quality-control sheet: packaged on an index (pip / conda / Bioconductor / CRAN / npm), a documentation site, a tutorial, a function-level API reference, automated tests, continuous integration, a release or commit in the last twelve months, a peer-reviewed paper, an open license, a programmatic (library) interface. | The share of the known items a tool fulfils, equal weights. |
| **Overall** | | The geometric mean of the aspects a tool has, so a tool cannot buy back a failed aspect with the others. |

Timings of the main suite (bedtools, pyranges 0.x, bioframe, intervaltree,
cgranges, ncls, deepTools, pyBigWig, pybigtools, FIMO, MOODS, lightmotif,
Biopython, cooler, graph-tool, igraph, networkx, scipy, tximport, edgeR,
DESeq2) come from the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page's
results; the tools added for this page (the Bioconductor stack, pyranges 1.x,
polars-bio, BiocPy, BEDOPS, HOMER, UROPA, EnrichedHeatmap, genomation, UCSC
tools, motifmatchr, universalmotif, Biostrings, LOLA, regioneR, GIGGLE,
GenomicInteractions, igraph, pytximport, PyDESeq2, RNAnorm) were run by
`benchmarks/bench_landscape.py` on the same data, as one SLURM job per theme on
cluster nodes with the same CPU model as the main suite's machine (Intel Xeon
Gold 6130, 6 CPUs per job; the parts are single-threaded except where noted).


## Interval operations
{: .sec-navy #intervals }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Overlap join, all A x B pairs (timed); **C2** Overlap filter / count (-u, -c, -v) (timed); **C3** Nearest / closest (ties, k nearest, strand) (timed); **C4** Coverage of A over B (counts, bp, fraction) (timed); **C5** Genome-wide depth track (timed); **C6** Merge / cluster with distance and column aggregation (timed); **C7** Subtract / difference (timed); **C8** Complement against a genome (timed); **C9** Parse + sort large BED / narrowPeak / GFF (timed); **C10** Sequence extraction, strand-aware (timed).

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| pyranges (0.x) | best | 0.82 | 10/10 | 0.82 | C2 0.76, C4 0.74, C6 0.92, C7 0.78, C9 0.88 |
| GenomicRanges + IRanges | best | 0.81 | 9/10 | 0.90 | C2 0.85, C3 0.89, C5 0.89, C6 0.95, C7 0.89, C8 0.92 |
| **genomeblocks** | best | 0.80 | 8/10 | 1.00 | C2 1.00, C4 1.00, C5 1.00, C7 1.00, C8 1.00, C9 1.00 |
| plyranges | best | 0.77 | 9/10 | 0.86 | C2 0.84, C3 0.90, C5 0.87, C6 0.96, C7 0.85, C8 0.73 |
| valr | best | 0.76 | 9/10 | 0.85 | C2 0.76, C4 0.86, C6 0.91, C7 0.91, C8 0.79 |
| pyranges (1.x) | good | 0.69 | 9/9 (1 unknown) | 0.69 | C2 0.62, C3 0.81, C5 0.70, C6 0.66, C7 0.77, C8 0.58 |
| BEDOPS | good | 0.67 | 9/10 | 0.74 | C2 0.63, C3 0.88, C4 0.69, C6 0.74, C7 0.81, C8 0.70 |
| bioframe | good | 0.67 | 10/10 | 0.67 | C2 0.64, C4 0.67, C6 0.71, C7 0.67, C8 0.64 |
| polars-bio | good | 0.66 | 9/10 | 0.73 | C2 0.63, C3 0.90, C4 0.70, C6 0.69, C7 0.83, C8 0.65 |
| bedtools | good | 0.61 | 10/10 | 0.61 | C2 0.56, C3 0.85, C4 0.54, C5 0.33, C6 0.72, C7 0.65, C8 0.60 |
| BiocPy genomicranges / iranges | good | 0.52 | 9/10 | 0.58 | C2 0.55, C3 0.89, C6 0.45, C8 0.42 |
| superintervals | moderate | 0.26 | 5/10 | 0.52 | C2 0.52 |
| ailist (Python) | low | 0.24 | 5/10 | 0.48 | C2 0.48 |
| data.table | low | 0.19 | 3/10 | 0.63 | C2 0.63 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| **genomeblocks** | 0.97 | 1.00 | 1.00 | 1.00 | 0.90 | 12 | 6 |
| GenomicRanges + IRanges | 0.97 | 1.00 | 0.89 | 1.00 | 1.00 | 7 | 6 |
| plyranges | 0.96 | 1.00 | 0.86 | 1.00 | 1.00 | 7 | 6 |
| valr | 0.95 | 1.00 | 0.83 | 1.00 | 1.00 | 7 | 6 |
| polars-bio | 0.92 | 1.00 | 0.72 | 1.00 | 1.00 | 7 | 6 |
| pyranges (0.x) | 0.90 | 1.00 | 0.73 | – | 1.00 | 7 | 0 |
| pyranges (1.x) | 0.88 | 1.00 | 0.68 | 1.00 | 0.90 | 7 | 6 |
| data.table | 0.87 | 1.00 | 0.63 | 1.00 | 0.90 | 1 | 2 |
| bedtools | 0.86 | 1.00 | 0.60 | 1.00 | 0.90 | 8 | 4 |
| intervaltree | 0.84 | 1.00 | 0.98 | – | 0.60 | 2 | 0 |
| BiocPy genomicranges / iranges | 0.83 | 0.93 | 0.57 | 1.00 | 0.90 | 6 | 6 |
| bioframe | 0.83 | 1.00 | 0.58 | – | 1.00 | 7 | 0 |
| ailist (Python) | 0.78 | 1.00 | 0.48 | 1.00 | 0.75 | 1 | 2 |
| superintervals | 0.75 | 1.00 | 0.52 | 1.00 | 0.60 | 1 | 2 |
| BEDOPS | 0.66 | 1.00 | 0.73 | 0.33 | 0.80 | 7 | 3 |

### `complement`
{: .no_toc }

Scored at n = 1M; fastest 74 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Loci.complement** | genomeblocks | 1M | 74 ms | 1.00 | identical (670,312) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 128 ms | 1.00 | identical (670,312) | landscape |
| valr | valr | 1M | 314 ms | 1.00 | identical (670,312) | landscape |
| plyranges | plyranges | 1M | 484 ms | 1.00 | identical (670,312) | landscape |
| BEDOPS bedops --difference genome (CLI) | BEDOPS | 1M | 585 ms | 1.00 | identical (670,312) | landscape |
| polars-bio 0.36.2 | polars-bio | 1M | 823 ms | 1.00 | identical (670,312) | landscape |
| bioframe 0.8.0 | bioframe | 1M | 908 ms | 1.00 | identical (670,312) | landscape |
| bedtools complement (CLI) | bedtools | 1M | 1.1 s | 1.00 | reference | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 1.3 s | 1.00 | identical (670,312) | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 1M | 4.0 s | 1.00 | identical (670,312) | landscape |
| **genomeblocks Loci.complement** | genomeblocks | 100k | 8.6 ms | 1.00 | identical (95,545) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 24 ms | 1.00 | identical (95,545) | landscape |
| valr | valr | 100k | 54 ms | 1.00 | identical (95,545) | landscape |
| BEDOPS bedops --difference genome (CLI) | BEDOPS | 100k | 86 ms | 1.00 | identical (95,545) | landscape |
| polars-bio 0.36.2 | polars-bio | 100k | 114 ms | 1.00 | identical (95,545) | landscape |
| plyranges | plyranges | 100k | 139 ms | 1.00 | identical (95,545) | landscape |
| bioframe 0.8.0 | bioframe | 100k | 149 ms | 1.00 | identical (95,545) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 160 ms | 1.00 | identical (95,545) | landscape |
| bedtools complement (CLI) | bedtools | 100k | 163 ms | 1.00 | reference | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 100k | 503 ms | 1.00 | identical (95,545) | landscape |
| **genomeblocks Loci.complement** | genomeblocks | 10k | 1.6 ms | 1.00 | identical (9,967) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 10k | 15 ms | 1.00 | identical (9,967) | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 21 ms | 1.00 | identical (9,967) | landscape |
| valr | valr | 10k | 22 ms | 1.00 | identical (9,967) | landscape |
| BEDOPS bedops --difference genome (CLI) | BEDOPS | 10k | 24 ms | 1.00 | identical (9,967) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 24 ms | 1.00 | identical (9,967) | landscape |
| bedtools complement (CLI) | bedtools | 10k | 34 ms | 1.00 | reference | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 10k | 63 ms | 1.00 | identical (9,967) | landscape |
| bioframe 0.8.0 | bioframe | 10k | 76 ms | 1.00 | identical (9,967) | landscape |
| plyranges | plyranges | 10k | 89 ms | 1.00 | identical (9,967) | landscape |

### `coverage`
{: .no_toc }

Scored at n = 1M; fastest 181 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Loci.coverage** | genomeblocks | 1M | 181 ms | 1.00 | 999,990/999,990 intervals agree | landscape |
| valr | valr | 1M | 478 ms | 1.00 | 999,990/999,990 intervals agree | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 1M | 1.1 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| polars-bio 0.36.2 | polars-bio | 1M | 1.4 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| BEDOPS bedmap --bases-uniq (CLI) | BEDOPS | 1M | 1.6 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| bioframe 0.8.0 | bioframe | 1M | 1.8 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| bedtools coverage (CLI) | bedtools | 1M | 4.3 s | 1.00 | reference | landscape |
| **genomeblocks Loci.coverage** | genomeblocks | 100k | 18 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| valr | valr | 100k | 48 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| polars-bio 0.36.2 | polars-bio | 100k | 149 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| BEDOPS bedmap --bases-uniq (CLI) | BEDOPS | 100k | 161 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| bioframe 0.8.0 | bioframe | 100k | 209 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 100k | 229 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| bedtools coverage (CLI) | bedtools | 100k | 367 ms | 1.00 | reference | landscape |
| **genomeblocks Loci.coverage** | genomeblocks | 10k | 2.2 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| valr | valr | 10k | 13 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 26 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| BEDOPS bedmap --bases-uniq (CLI) | BEDOPS | 10k | 34 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| bedtools coverage (CLI) | bedtools | 10k | 52 ms | 1.00 | reference | landscape |
| bioframe 0.8.0 | bioframe | 10k | 57 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 10k | 142 ms | 1.00 | 10,000/10,000 intervals agree | landscape |

### `genome_depth`
{: .no_toc }

Scored at n = 1M; fastest 243 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks hichip.coverage** | genomeblocks | 1M | 243 ms | 1.00 | identical (1,326,511) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 531 ms | 1.00 | identical (1,326,511) | landscape |
| plyranges | plyranges | 1M | 576 ms | 1.00 | identical (1,326,511) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 2.0 s | 1.00 | identical (1,326,511) | landscape |
| bedtools genomecov -bg (CLI) | bedtools | 1M | 24 s | 1.00 | reference | landscape |
| **genomeblocks hichip.coverage** | genomeblocks | 100k | 18 ms | 1.00 | identical (104,427) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 71 ms | 1.00 | identical (104,427) | landscape |
| plyranges | plyranges | 100k | 92 ms | 1.00 | identical (104,427) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 162 ms | 1.00 | identical (104,427) | landscape |
| bedtools genomecov -bg (CLI) | bedtools | 100k | 20 s | 1.00 | reference | landscape |
| **genomeblocks hichip.coverage** | genomeblocks | 10k | 2.2 ms | 1.00 | identical (10,054) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 24 ms | 1.00 | identical (10,054) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 10k | 33 ms | 1.00 | identical (10,054) | landscape |
| plyranges | plyranges | 10k | 41 ms | 1.00 | identical (10,054) | landscape |
| bedtools genomecov -bg (CLI) | bedtools | 10k | 20 s | 1.00 | reference | landscape |

### `getfasta`
{: .no_toc }

Scored at n = 10k; fastest 5.4 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **memory (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 5.4 ms | 1.00 |  | backends |
| **genomeblocks (indexed FASTA reader)** | genomeblocks | 10k | 37 ms | 1.00 |  | backends |
| **pysam (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 40 ms | 1.00 |  | backends |
| **pyfaidx (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 135 ms | 1.00 |  | backends |
| **biopython (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 5.0 s | 1.00 |  | backends |

### `intersect`
{: .no_toc }

Scored at n = 1M; fastest 87 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 1M | 87 ms | 1.00 |  | loci |
| **genomeblocks** | genomeblocks | 1M | 95 ms | 1.00 | identical (729,413) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 276 ms | 1.00 | identical (729,413) | landscape |
| plyranges | plyranges | 1M | 311 ms | 1.00 | identical (729,413) | landscape |
| pyranges | pyranges (0.x) | 1M | 336 ms | 1.00 |  | loci |
| valr | valr | 1M | 558 ms | 1.00 | differs: 729,696 vs 729,413 (Jaccard 0.9996) | landscape |
| bioframe | bioframe | 1M | 1.0 s | 1.00 |  | loci |
| data.table::foverlaps | data.table | 1M | 1.1 s | 1.00 | identical (729,413) | landscape |
| BEDOPS bedops --element-of (CLI) | BEDOPS | 1M | 1.1 s | 1.00 | identical (729,413) | landscape |
| polars-bio 0.36.2 | polars-bio | 1M | 1.3 s | 1.00 | identical (729,413) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 1.3 s | 1.00 | identical (729,413) | landscape |
| **genomeblocks (cgranges)** | genomeblocks | 1M | 1.6 s | 1.00 |  | loci |
| bedtools (CLI) | bedtools | 1M | 1.8 s | 1.00 |  | loci |
| bedtools intersect -u (CLI) | bedtools | 1M | 2.1 s | 1.00 | reference | landscape |
| superintervals 1.0.3 | superintervals | 1M | 2.4 s | 1.00 | identical (729,413) | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 1M | 3.0 s | 1.00 | identical (729,413) | landscape |
| AIList (ailist) 2.2.2 | ailist (Python) | 1M | 3.1 s | 1.00 | identical (729,413) | landscape |
| **genomeblocks** | genomeblocks | 100k | 7.4 ms | 1.00 |  | loci |
| **genomeblocks** | genomeblocks | 100k | 9.2 ms | 1.00 | identical (46,096) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 36 ms | 1.00 | identical (46,096) | landscape |
| plyranges | plyranges | 100k | 37 ms | 1.00 | identical (46,096) | landscape |
| valr | valr | 100k | 39 ms | 1.00 | differs: 46,105 vs 46,096 (Jaccard 0.9998) | landscape |
| bioframe | bioframe | 100k | 61 ms | 1.00 |  | loci |
| pyranges | pyranges (0.x) | 100k | 62 ms | 1.00 |  | loci |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 83 ms | 1.00 | identical (46,096) | landscape |
| data.table::foverlaps | data.table | 100k | 90 ms | 1.00 | identical (46,096) | landscape |
| **genomeblocks (cgranges)** | genomeblocks | 100k | 99 ms | 1.00 |  | loci |
| polars-bio 0.36.2 | polars-bio | 100k | 113 ms | 1.00 | identical (46,096) | landscape |
| BEDOPS bedops --element-of (CLI) | BEDOPS | 100k | 120 ms | 1.00 | identical (46,096) | landscape |
| bedtools (CLI) | bedtools | 100k | 169 ms | 1.00 |  | loci |
| bedtools intersect -u (CLI) | bedtools | 100k | 182 ms | 1.00 | reference | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 100k | 201 ms | 1.00 | identical (46,096) | landscape |
| superintervals 1.0.3 | superintervals | 100k | 228 ms | 1.00 | identical (46,096) | landscape |
| AIList (ailist) 2.2.2 | ailist (Python) | 100k | 288 ms | 1.00 | identical (46,096) | landscape |
| intervaltree | intervaltree | 100k | 22 s | 1.00 |  | loci |
| **genomeblocks** | genomeblocks | 10k | 901 µs | 1.00 |  | loci |
| **genomeblocks** | genomeblocks | 10k | 1.2 ms | 1.00 | identical (4,069) | landscape |
| **genomeblocks (cgranges)** | genomeblocks | 10k | 9.8 ms | 1.00 |  | loci |
| bioframe | bioframe | 10k | 12 ms | 1.00 |  | loci |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 12 ms | 1.00 | identical (4,069) | landscape |
| valr | valr | 10k | 12 ms | 1.00 | identical (4,069) | landscape |
| data.table::foverlaps | data.table | 10k | 15 ms | 1.00 | identical (4,069) | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 18 ms | 1.00 | identical (4,069) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 10k | 18 ms | 1.00 | identical (4,069) | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 10k | 20 ms | 1.00 | identical (4,069) | landscape |
| plyranges | plyranges | 10k | 20 ms | 1.00 | identical (4,069) | landscape |
| BEDOPS bedops --element-of (CLI) | BEDOPS | 10k | 28 ms | 1.00 | identical (4,069) | landscape |
| bedtools intersect -u (CLI) | bedtools | 10k | 34 ms | 1.00 | reference | landscape |
| bedtools (CLI) | bedtools | 10k | 36 ms | 1.00 |  | loci |
| superintervals 1.0.3 | superintervals | 10k | 37 ms | 1.00 | identical (4,069) | landscape |
| AIList (ailist) 2.2.2 | ailist (Python) | 10k | 43 ms | 1.00 | identical (4,069) | landscape |
| pyranges | pyranges (0.x) | 10k | 44 ms | 1.00 |  | loci |
| intervaltree | intervaltree | 10k | 464 ms | 1.00 |  | loci |
| **genomeblocks** | genomeblocks | 1k | 268 µs | 1.00 |  | loci |
| **genomeblocks (cgranges)** | genomeblocks | 1k | 1.1 ms | 1.00 |  | loci |
| bioframe | bioframe | 1k | 7.3 ms | 1.00 |  | loci |
| bedtools (CLI) | bedtools | 1k | 16 ms | 1.00 |  | loci |
| intervaltree | intervaltree | 1k | 24 ms | 1.00 |  | loci |
| pyranges | pyranges (0.x) | 1k | 49 ms | 1.00 |  | loci |

### `merge`
{: .no_toc }

Scored at n = 1M; fastest 95 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 1M | 95 ms | 1.00 | identical (670,293) | landscape |
| plyranges | plyranges | 1M | 129 ms | 1.00 | identical (670,293) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 130 ms | 1.00 | identical (670,293) | landscape |
| **genomeblocks** | genomeblocks | 1M | 135 ms | – |  | loci |
| pyranges | pyranges (0.x) | 1M | 164 ms | – |  | loci |
| valr | valr | 1M | 175 ms | 1.00 | identical (670,293) | landscape |
| BEDOPS bedops --merge (CLI) | BEDOPS | 1M | 554 ms | 1.00 | identical (670,293) | landscape |
| bedtools merge (CLI) | bedtools | 1M | 638 ms | 1.00 | reference | landscape |
| bioframe | bioframe | 1M | 699 ms | – |  | loci |
| polars-bio 0.36.2 | polars-bio | 1M | 810 ms | 1.00 | differs: 670,804 vs 670,293 (Jaccard 0.9977) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 980 ms | 1.00 | differs: 670,804 vs 670,293 (Jaccard 0.9977) | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 1M | 4.3 s | 1.00 | identical (670,293) | landscape |
| bedtools (CLI) | bedtools | 1M | 5.1 s | – |  | loci |
| **genomeblocks** | genomeblocks | 100k | 10 ms | 1.00 | identical (95,525) | landscape |
| **genomeblocks** | genomeblocks | 100k | 18 ms | – |  | loci |
| plyranges | plyranges | 100k | 25 ms | 1.00 | identical (95,525) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 25 ms | 1.00 | identical (95,525) | landscape |
| valr | valr | 100k | 27 ms | 1.00 | identical (95,525) | landscape |
| pyranges | pyranges (0.x) | 100k | 52 ms | – |  | loci |
| BEDOPS bedops --merge (CLI) | BEDOPS | 100k | 81 ms | 1.00 | identical (95,525) | landscape |
| bedtools merge (CLI) | bedtools | 100k | 92 ms | 1.00 | reference | landscape |
| bioframe | bioframe | 100k | 102 ms | – |  | loci |
| polars-bio 0.36.2 | polars-bio | 100k | 120 ms | 1.00 | differs: 95,535 vs 95,525 (Jaccard 0.9997) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 129 ms | 1.00 | differs: 95,535 vs 95,525 (Jaccard 0.9997) | landscape |
| bedtools (CLI) | bedtools | 100k | 543 ms | – |  | loci |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 100k | 561 ms | 1.00 | identical (95,525) | landscape |
| **genomeblocks** | genomeblocks | 10k | 1.4 ms | 1.00 | identical (9,946) | landscape |
| **genomeblocks** | genomeblocks | 10k | 2.0 ms | – |  | loci |
| GenomicRanges | GenomicRanges + IRanges | 10k | 11 ms | 1.00 | identical (9,946) | landscape |
| plyranges | plyranges | 10k | 13 ms | 1.00 | identical (9,946) | landscape |
| valr | valr | 10k | 14 ms | 1.00 | identical (9,946) | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 17 ms | 1.00 | identical (9,946) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 18 ms | 1.00 | identical (9,946) | landscape |
| BEDOPS bedops --merge (CLI) | BEDOPS | 10k | 23 ms | 1.00 | identical (9,946) | landscape |
| bioframe | bioframe | 10k | 23 ms | – |  | loci |
| bedtools merge (CLI) | bedtools | 10k | 24 ms | 1.00 | reference | landscape |
| pyranges | pyranges (0.x) | 10k | 30 ms | – |  | loci |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 10k | 65 ms | 1.00 | identical (9,946) | landscape |
| bedtools (CLI) | bedtools | 10k | 92 ms | – |  | loci |
| **genomeblocks** | genomeblocks | 1k | 366 µs | – |  | loci |
| bioframe | bioframe | 1k | 16 ms | – |  | loci |
| bedtools (CLI) | bedtools | 1k | 27 ms | – |  | loci |
| pyranges | pyranges (0.x) | 1k | 28 ms | – |  | loci |

### `nearest`
{: .no_toc }

Scored at n = 1M; fastest 722 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 1M | 722 ms | 1.00 | 1,000,000/1,000,000 rows agree with bedtools closest | landscape |
| plyranges | plyranges | 1M | 1.4 s | 1.00 | 1,000,000/1,000,000 rows agree | landscape |
| polars-bio 0.36.2 | polars-bio | 1M | 1.5 s | 1.00 | 1,000,000/1,000,000 rows agree | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 1.6 s | 1.00 | 1,000,000/1,000,000 rows agree | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 1M | 1.6 s | 1.00 | 1,000,000/1,000,000 rows agree | landscape |
| BEDOPS closest-features (CLI) | BEDOPS | 1M | 1.6 s | 1.00 | 999,934/1,000,000 rows agree | landscape |
| bedtools closest (CLI) | bedtools | 1M | 2.0 s | 1.00 | reference | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 2.7 s | 1.00 | 1,828,167/1,828,167 rows agree | landscape |
| **genomeblocks** | genomeblocks | 100k | 30 ms | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 100k | 37 ms | 1.00 | 100,000/100,000 rows agree with bedtools closest | landscape |
| valr | valr | 100k | 75 ms | 1.00 | 100,000/100,000 rows agree | landscape |
| polars-bio 0.36.2 | polars-bio | 100k | 147 ms | 1.00 | 100,000/100,000 rows agree | landscape |
| BEDOPS closest-features (CLI) | BEDOPS | 100k | 166 ms | 1.00 | 100,000/100,000 rows agree | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 100k | 175 ms | 1.00 | 100,000/100,000 rows agree | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 183 ms | 1.00 | 105,435/105,435 rows agree | landscape |
| bedtools closest (CLI) | bedtools | 100k | 209 ms | 1.00 | reference | landscape |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 228 ms | 1.00 |  | backends |
| plyranges | plyranges | 100k | 262 ms | 1.00 | 100,000/100,000 rows agree | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 265 ms | 1.00 | 100,000/100,000 rows agree | landscape |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 392 ms | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 1.8 s | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 10k | 3.2 ms | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 10k | 4.1 ms | 1.00 | 10,000/10,000 rows agree with bedtools closest | landscape |
| valr | valr | 10k | 16 ms | 1.00 | 10,000/10,000 rows agree | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 10k | 19 ms | 1.00 | 10,000/10,000 rows agree | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 19 ms | 1.00 | 10,000/10,000 rows agree | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 22 ms | 1.00 | 10,046/10,046 rows agree | landscape |
| BEDOPS closest-features (CLI) | BEDOPS | 10k | 33 ms | 1.00 | 10,000/10,000 rows agree | landscape |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 37 ms | 1.00 |  | backends |
| bedtools closest (CLI) | bedtools | 10k | 39 ms | 1.00 | reference | landscape |
| GenomicRanges | GenomicRanges + IRanges | 10k | 74 ms | 1.00 | 10,000/10,000 rows agree | landscape |
| plyranges | plyranges | 10k | 92 ms | 1.00 | 10,000/10,000 rows agree | landscape |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 214 ms | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 290 ms | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 1k | 695 µs | 1.00 |  | backends |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 20 ms | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 126 ms | 1.00 |  | backends |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 190 ms | 1.00 |  | backends |

### `overlap_count`
{: .no_toc }

Scored at n = 1M; fastest 115 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Loci.count_overlaps** | genomeblocks | 1M | 115 ms | 1.00 | 999,990/999,990 intervals agree | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 289 ms | 1.00 | 999,990/999,990 intervals agree | landscape |
| plyranges | plyranges | 1M | 289 ms | 1.00 | 999,990/999,990 intervals agree | landscape |
| valr | valr | 1M | 489 ms | 1.00 | 998,030/999,990 intervals agree | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 1M | 779 ms | 1.00 | 999,990/999,990 intervals agree | landscape |
| polars-bio 0.36.2 | polars-bio | 1M | 1.4 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| bioframe 0.8.0 | bioframe | 1M | 1.4 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 1.5 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| BEDOPS bedmap --count (CLI) | BEDOPS | 1M | 1.5 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 1M | 1.6 s | 1.00 | 999,990/999,990 intervals agree | landscape |
| bedtools intersect -c (CLI) | bedtools | 1M | 2.3 s | 1.00 | reference | landscape |
| **genomeblocks Loci.count_overlaps** | genomeblocks | 100k | 9.5 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 37 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| plyranges | plyranges | 100k | 38 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| valr | valr | 100k | 48 ms | 1.00 | 99,984/100,000 intervals agree | landscape |
| bioframe 0.8.0 | bioframe | 100k | 116 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 100k | 133 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| polars-bio 0.36.2 | polars-bio | 100k | 142 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 143 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 100k | 155 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| BEDOPS bedmap --count (CLI) | BEDOPS | 100k | 157 ms | 1.00 | 100,000/100,000 intervals agree | landscape |
| bedtools intersect -c (CLI) | bedtools | 100k | 237 ms | 1.00 | reference | landscape |
| **genomeblocks Loci.count_overlaps** | genomeblocks | 10k | 1.1 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| valr | valr | 10k | 15 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| GenomicRanges | GenomicRanges + IRanges | 10k | 17 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 10k | 19 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 20 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| plyranges | plyranges | 10k | 20 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| bioframe 0.8.0 | bioframe | 10k | 22 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 32 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| BEDOPS bedmap --count (CLI) | BEDOPS | 10k | 32 ms | 1.00 | 10,000/10,000 intervals agree | landscape |
| bedtools intersect -c (CLI) | bedtools | 10k | 36 ms | 1.00 | reference | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 10k | 78 ms | 1.00 | 10,000/10,000 intervals agree | landscape |

### `overlap_pairs`
{: .no_toc }

Scored at n = 1M; fastest 288 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **ncls (through genomeblocks&#x27; backend)** | genomeblocks | 1M | 288 ms | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 1M | 291 ms | 1.00 |  | backends |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 1M | 966 ms | 1.00 |  | backends |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 1M | 1.3 s | 1.00 |  | backends |
| **cgranges (through genomeblocks&#x27; backend)** | genomeblocks | 1M | 1.4 s | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 1M | 9.3 s | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 100k | 11 ms | 1.00 |  | backends |
| **ncls (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 20 ms | 1.00 |  | backends |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 83 ms | 1.00 |  | backends |
| **cgranges (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 97 ms | 1.00 |  | backends |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 131 ms | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 100k | 670 ms | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 10k | 1.2 ms | 1.00 |  | backends |
| **ncls (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 2.3 ms | 1.00 |  | backends |
| **cgranges (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 9.3 ms | 1.00 |  | backends |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 13 ms | 1.00 |  | backends |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 57 ms | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 102 ms | 1.00 |  | backends |
| **genomeblocks** | genomeblocks | 1k | 347 µs | 1.00 |  | backends |
| **ncls (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 897 µs | 1.00 |  | backends |
| **cgranges (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 1.2 ms | 1.00 |  | backends |
| **bioframe (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 7.1 ms | 1.00 |  | backends |
| **bedtools (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 30 ms | 1.00 |  | backends |
| **pyranges (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 57 ms | 1.00 |  | backends |

### `point_lookup`
{: .no_toc }

Scored at n = 100k; fastest 9 µs.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks (cgranges)** | genomeblocks | 100k | 9 µs | – |  | loci |
| **genomeblocks (ncls)** | genomeblocks | 100k | 10 µs | – |  | loci |
| intervaltree | intervaltree | 100k | 10 µs | – |  | loci |
| **genomeblocks (numpy point index)** | genomeblocks | 100k | 14 µs | – |  | loci |
| pyranges | pyranges (0.x) | 100k | 1.5 ms | – |  | loci |
| bioframe | bioframe | 100k | 5.4 ms | – |  | loci |

### `read_bed`
{: .no_toc }

Scored at n = 1M; fastest 394 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 1M | 394 ms | – |  | loci |
| pyranges | pyranges (0.x) | 1M | 928 ms | – |  | loci |
| **genomeblocks** | genomeblocks | 100k | 56 ms | – |  | loci |
| pyranges | pyranges (0.x) | 100k | 91 ms | – |  | loci |
| **genomeblocks** | genomeblocks | 10k | 6.6 ms | – |  | loci |
| pyranges | pyranges (0.x) | 10k | 13 ms | – |  | loci |
| **genomeblocks** | genomeblocks | 1k | 2.2 ms | – |  | loci |
| pyranges | pyranges (0.x) | 1k | 5.7 ms | – |  | loci |

### `subtract`
{: .no_toc }

Scored at n = 1M; fastest 259 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Loci.subtract** | genomeblocks | 1M | 259 ms | 1.00 | identical (508,459) | landscape |
| valr | valr | 1M | 499 ms | 1.00 | identical (508,459) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 1M | 548 ms | 1.00 | identical (508,459) | landscape |
| plyranges | plyranges | 1M | 713 ms | 1.00 | identical (508,459) | landscape |
| polars-bio 0.36.2 | polars-bio | 1M | 815 ms | 1.00 | identical (508,459) | landscape |
| BEDOPS bedops --difference (CLI) | BEDOPS | 1M | 966 ms | 1.00 | identical (508,459) | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 1M | 1.2 s | 1.00 | identical (508,459) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 1M | 1.3 s | 1.00 | identical (508,459) | landscape |
| bioframe 0.8.0 | bioframe | 1M | 2.6 s | 1.00 | identical (508,459) | landscape |
| bedtools subtract (CLI) | bedtools | 1M | 3.0 s | 1.00 | reference | landscape |
| **genomeblocks Loci.subtract** | genomeblocks | 100k | 22 ms | 1.00 | identical (81,436) | landscape |
| valr | valr | 100k | 68 ms | 1.00 | identical (81,436) | landscape |
| polars-bio 0.36.2 | polars-bio | 100k | 109 ms | 1.00 | identical (81,436) | landscape |
| BEDOPS bedops --difference (CLI) | BEDOPS | 100k | 122 ms | 1.00 | identical (81,436) | landscape |
| GenomicRanges | GenomicRanges + IRanges | 100k | 129 ms | 1.00 | identical (81,436) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 100k | 136 ms | 1.00 | identical (81,436) | landscape |
| plyranges | plyranges | 100k | 163 ms | 1.00 | identical (81,436) | landscape |
| bedtools subtract (CLI) | bedtools | 100k | 260 ms | 1.00 | reference | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 100k | 297 ms | 1.00 | identical (81,436) | landscape |
| bioframe 0.8.0 | bioframe | 100k | 333 ms | 1.00 | identical (81,436) | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 100k | 3.0 min | 0.42 | differs: 95,519 vs 81,436 (Jaccard 0.4199) | landscape |
| **genomeblocks Loci.subtract** | genomeblocks | 10k | 2.5 ms | 1.00 | identical (8,674) | landscape |
| polars-bio 0.36.2 | polars-bio | 10k | 15 ms | 1.00 | identical (8,674) | landscape |
| pyranges 1.5.0 | pyranges (1.x) | 10k | 19 ms | 1.00 | identical (8,674) | landscape |
| BEDOPS bedops --difference (CLI) | BEDOPS | 10k | 29 ms | 1.00 | identical (8,674) | landscape |
| valr | valr | 10k | 30 ms | 1.00 | identical (8,674) | landscape |
| bedtools subtract (CLI) | bedtools | 10k | 46 ms | 1.00 | reference | landscape |
| GenomicRanges | GenomicRanges + IRanges | 10k | 84 ms | 1.00 | identical (8,674) | landscape |
| plyranges | plyranges | 10k | 86 ms | 1.00 | identical (8,674) | landscape |
| bioframe 0.8.0 | bioframe | 10k | 114 ms | 1.00 | identical (8,674) | landscape |
| pyranges 0.1.4 | pyranges (0.x) | 10k | 202 ms | 1.00 | identical (8,674) | landscape |
| BiocPy genomicranges 0.9.0 | BiocPy genomicranges / iranges | 10k | 14 s | 0.47 | differs: 9,944 vs 8,674 (Jaccard 0.4662) | landscape |

### Not run or did not finish
{: .no_toc }

| Engine | task | why |
|---|---|---|
| BiocPy genomicranges 0.9.0 | genome_depth | not run: GenomicRanges.coverage() returns a dense per-base float64 array per chromosome (about 25 GB for hg38) |
| BiocPy genomicranges 0.9.0 | genome_depth | not run: GenomicRanges.coverage() returns a dense per-base float64 array per chromosome (about 25 GB for hg38) |
| valr | nearest | not run at n = 1000000: bed_closest did not finish at 1 M (4.5 h, 0.2% CPU); timed at 100k |
| BiocPy genomicranges 0.9.0 | subtract | not run at n = 1,000,000: 3,000 s per run at 1 M on 2026-10-23 (superlinear; Jaccard 0.24); timed at 10k and 100k |
| BiocPy genomicranges 0.9.0 | genome_depth | not run: GenomicRanges.coverage() returns a dense per-base float64 array per chromosome (about 25 GB for hg38) |

### Perturbation tests
{: .no_toc }

| Tool | Engine | task | variant | identical answer |
|---|---|---|---|---|
| genomeblocks | genomeblocks | intersect | chrom-order | yes |
| genomeblocks | genomeblocks | intersect | shuffled | yes |
| genomeblocks | genomeblocks | merge | chrom-order | yes |
| genomeblocks | genomeblocks | merge | shuffled | yes |
| genomeblocks | genomeblocks | nearest | chrom-order | yes |
| genomeblocks | genomeblocks | nearest | shuffled | yes |
| BEDOPS | BEDOPS (CLI) | intersect | chrom-order | **no** |
| BEDOPS | BEDOPS (CLI) | merge | chrom-order | yes |
| BEDOPS | BEDOPS (CLI) | merge | shuffled | **no** |
| BiocPy genomicranges / iranges | BiocPy genomicranges 0.9.0 | intersect | chrom-order | yes |
| BiocPy genomicranges / iranges | BiocPy genomicranges 0.9.0 | intersect | shuffled | yes |
| BiocPy genomicranges / iranges | BiocPy genomicranges 0.9.0 | merge | chrom-order | yes |
| BiocPy genomicranges / iranges | BiocPy genomicranges 0.9.0 | merge | shuffled | yes |
| BiocPy genomicranges / iranges | BiocPy genomicranges 0.9.0 | nearest | chrom-order | yes |
| BiocPy genomicranges / iranges | BiocPy genomicranges 0.9.0 | nearest | shuffled | yes |
| GenomicRanges + IRanges | GenomicRanges | intersect | chrom-order | yes |
| GenomicRanges + IRanges | GenomicRanges | intersect | shuffled | yes |
| GenomicRanges + IRanges | GenomicRanges | merge | chrom-order | yes |
| GenomicRanges + IRanges | GenomicRanges | merge | shuffled | yes |
| GenomicRanges + IRanges | GenomicRanges | nearest | chrom-order | yes |
| GenomicRanges + IRanges | GenomicRanges | nearest | shuffled | yes |
| ailist (Python) | AIList (ailist) 2.2.2 | intersect | chrom-order | yes |
| ailist (Python) | AIList (ailist) 2.2.2 | intersect | shuffled | yes |
| bedtools | bedtools (CLI) | intersect | chrom-order | yes |
| bedtools | bedtools (CLI) | intersect | shuffled | yes |
| bedtools | bedtools (CLI) | merge | chrom-order | yes |
| bedtools | bedtools (CLI) | nearest | chrom-order | yes |
| data.table | data.table::foverlaps | intersect | chrom-order | yes |
| data.table | data.table::foverlaps | intersect | shuffled | yes |
| plyranges | plyranges | intersect | chrom-order | yes |
| plyranges | plyranges | intersect | shuffled | yes |
| plyranges | plyranges | merge | chrom-order | yes |
| plyranges | plyranges | merge | shuffled | yes |
| plyranges | plyranges | nearest | chrom-order | yes |
| plyranges | plyranges | nearest | shuffled | yes |
| polars-bio | polars-bio 0.36.2 | intersect | chrom-order | yes |
| polars-bio | polars-bio 0.36.2 | intersect | shuffled | yes |
| polars-bio | polars-bio 0.36.2 | merge | chrom-order | yes |
| polars-bio | polars-bio 0.36.2 | merge | shuffled | yes |
| polars-bio | polars-bio 0.36.2 | nearest | chrom-order | yes |
| polars-bio | polars-bio 0.36.2 | nearest | shuffled | yes |
| pyranges (1.x) | pyranges 1.5.0 | intersect | chrom-order | yes |
| pyranges (1.x) | pyranges 1.5.0 | intersect | shuffled | yes |
| pyranges (1.x) | pyranges 1.5.0 | merge | chrom-order | yes |
| pyranges (1.x) | pyranges 1.5.0 | merge | shuffled | yes |
| pyranges (1.x) | pyranges 1.5.0 | nearest | chrom-order | yes |
| pyranges (1.x) | pyranges 1.5.0 | nearest | shuffled | yes |
| superintervals | superintervals 1.0.3 | intersect | chrom-order | yes |
| superintervals | superintervals 1.0.3 | intersect | shuffled | yes |
| valr | valr | intersect | chrom-order | yes |
| valr | valr | intersect | shuffled | yes |
| valr | valr | merge | chrom-order | yes |
| valr | valr | merge | shuffled | yes |
| valr | valr | nearest | chrom-order | yes |
| valr | valr | nearest | shuffled | yes |

### Every tool in this theme
{: .no_toc }

Capability checklist: `read_bed`, `overlap`, `set_algebra`, `merge`, `slop`, `sort`, `nearest`, `point_lookup_index`, `tile`, `liftover`, `getfasta`, `dataframe_interop`, `coverage`, `jaccard_fisher`, `out_of_core`, `parallel`, `columnar_storage`, `shuffle_regions`, `multi_set_intersect`, `window_search`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 75% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| pybedtools | Python | library | pip · conda | 2011 | 0.12.1 2026-08 | yes | 85% | 1.00 | Requires the bedtools binary; IntervalFile.all_hits gives an in-memory index; to_dataframe/from_dataframe for pandas; pybedtools.parallel for multi-process shuffles. |
| pyranges (0.x) | Python | library | pip · conda | 2020 | 0.1.6 2026-10 | yes | 80% | 1.00 | Legacy line still released (0.1.6, 2026-10) but README says v1 &#x27;will soon replace&#x27; it; nb_cpu multiprocessing via ray; stats module (jaccard/forbes/fisher); features.tss()/tes()/introns() from GTF. |
| bedtools | C/C++ | cli | conda | 2010 | v2.31.1 2023-11 | yes | 70% | 0.90 | Reference CLI for genome arithmetic; set ops limited to subtract/complement (union = cat|sort|merge); &#x27;closest&#x27; against a TSS BED is the annotation use; last tagged release 2023 but repo pushed 2026-06; genomecov -bg -scale is the classic BAM-&gt;bedGraph route; GitHub push 2026-06 but no release since 2023-11; pybedtools wraps it; Primitives only: fisher (2x2 on one pair), jaccard, reldist, shuffle (-incl/-excl backgrounds) and pairtobed -type either/both for loop anchors; collection ranking must be scripted. |
| pyranges (1.x) | Python | mixed | pip | 2025 | 1.5.0 2026-10 | yes | 65% | 0.90 | Complete rewrite published as package &#x27;pyranges1&#x27; (requires Python&gt;=3.12); set_num_threads for the Rust engine; read_chain exists but liftover not verified; also a &#x27;pyranger&#x27; CLI; preprint DOI copied verbatim from the README. |
| valr | R | library | CRAN | 2017 | 0.10.0 2026-07 | yes | 65% | 1.00 | bedtools semantics in the tidyverse (bed_intersect/closest/merge/slop/makewindows/jaccard/fisher/projection); create_tss/create_introns/create_utrs from BED12; read_gtf is deprecated. |
| BEDOPS | C/C++ | cli | conda | 2012 | v2.4.42 2025-04 | no | 60% | 0.80 | First-class union/difference/symmetric-difference (-u/-d/-s) over any number of sorted files; sort-bed --max-mem for out-of-core; gtf2bed/gff2bed converters; last push 2025-04. |
| GenomicRanges + IRanges | R | library | Bioc | 2013 | 1.64.0 (IRanges 2.46.0) 2026-04 | yes | 60% | 1.00 | The de-facto R standard; reading BED, liftOver and getSeq live in rtracklayer/BSgenome; GNCList() gives a reusable NCList index; parallelism via BiocParallel outside the package. |
| plyranges | R | library | Bioc | 2019 | 1.32.0 2026-04 | yes | 60% | 1.00 | Grammar (join_overlap_*, join_nearest, reduce_ranges, stretch/anchor, tile_ranges, compute_coverage) on top of GenomicRanges; now maintained under tidyomics. |
| BiocPy genomicranges / iranges | Python | library | pip | – | 0.9.0 2026-09 | yes | 55% | 0.90 | Genentech&#x27;s port of the Bioconductor API (intra/inter-range, set ops, tiles, coverage, find_overlaps/nearest, promoters); read_gtf/read_ucsc; from_pandas/from_polars; BioC2023 abstract, no paper. |
| bioframe | Python | library | pip · conda | 2024 | 0.8.0 2025-04 | yes | 55% | 1.00 | Pure-pandas ops (overlap/closest/cluster/merge/complement/subtract/setdiff/coverage/expand/binnify); getfasta via load_fasta (pysam) + frac_gc; no union op, no index, single-threaded. |
| polars-bio | mixed | library | pip | 2025 | 0.36.2 2026-10 | yes | 45% | 1.00 | overlap/nearest/count_overlaps/coverage/merge/cluster/complement/subtract with lazy scans, streaming sinks and cloud (S3/GCS) reads; readers for BED/GFF3/GTF/VCF/BAM/FASTQ etc.; no union/liftover. |
| bx-python | Python | library | pip · conda | – | 0.15.1 2026-08 | yes | 40% | 0.80 | Galaxy-era toolbox (IntervalTree.before/after/upstream_of_interval, binned bitsets, bed_* scripts, 2bit/bigWig access); still released but mostly maintenance. |
| granges | Rust | mixed | src | – | 0.2.2 2024-03 | no | 35% | 0.60 | Vince Buffalo&#x27;s Rust library + CLI (adjust, filter, flank, map_*, windows); claims 10-30% faster than bedtools on some ops; inactive since 2024-05. |
| ailist (Python) | Python | library | pip | – | 2.2.2 2026-08 | yes | 30% | 0.75 | Single-chromosome AIList with intersect/subtract/common/merge/coverage/wps; numpy in, pandas Series out; no nearest. |
| bedr | R | library | CRAN | 2016 | 1.1.5 2026-01 | yes | 30% | 0.90 | Thin R wrapper that shells out to bedtools/BEDOPS/tabix (must be installed); jaccard/reldist wrappers; low activity apart from CRAN maintenance. |
| bedtk | C/C++ | cli | conda | – | v1.2 2025-08 | no | 30% | 0.30 | Minimal, very fast subset of bedtools (isec, flt, cov, sort, merge); README warns the CLI may change; no sorting required; Engine only (fast overlap), no statistics; v1.2 tagged 2025-08-24, just before the activity window. |
| gtars | Rust | mixed | pip | – | 0.10.0 2026-09 | yes | 30% | 0.90 | Sheffield-lab Rust workspace (overlaprs/igd, uniwig, tokenizers for geniml, lola enrichment, refget, bbcache, fragsplit); ML-oriented rather than a bedtools replacement. |
| HTSeq GenomicArray | Python | library | pip · conda | 2015 | 2.1.2 2026-02 | yes | 30% | 1.00 | Interval-indexed arrays for counting/coverage (htseq-count) plus a streaming GFF_Reader; feature classes must be built by the user; HTSeq 2.0 paper 2022 (btac166). |
| rtracklayer | R | library | Bioc | 2009 | 1.72.0 2026-04 | yes | 30% | 1.00 | The I/O layer of Bioconductor: import(GTF/GFF/BED), export.bed (BED12 via asBED), export.gff, liftOver(chain), indexed &#x27;which=&#x27; queries on bigWig/tabix; no arithmetic of its own; The R bigWig I/O layer every Bioconductor profile package builds on; summary() gives per-region means. |
| data.table | R | library | CRAN | – | 1.18.6.1 2026-08 | yes | 25% | 0.90 | General-purpose table library; foverlaps() is a single overlap-join primitive, not a genomics toolkit; data.table itself is OpenMP-parallel but foverlaps parallelism is unverified. |
| GenomeKit | mixed | library | conda | – | 7.7.0 2026-10 | yes | 25% | 0.90 | Prebuilt GENCODE/RefSeq/Ensembl/MANE annotation binaries queried by interval (find_overlapping etc.), dna(), tracks, VCF/BAM access; capability flags for overlap/slop/parse_gtf from memory of the API, docs page did not enumerate them. |
| GIGGLE | C/C++ | mixed | conda | 2018 | – | yes | 25% | 0.70 | Search engine over thousands of interval files with odds-ratio/Fisher/combo-score enrichment; Travis CI config only; no tagged releases; Origin of the Fisher x odds-ratio &#x27;combo score&#x27;; metadata index/filter (-m/-u); no GitHub releases or bioconda recipe found (conda env file only); python/go bindings are community-maintained and dated; last commit 2026-03. |
| IGD | C/C++ | cli | src | 2021 | – | no | 25% | 0.22 | Precursor of gtars-igd; single/batch queries and Seqpare similarity; R wrapper IGDr; last push 2021; GIGGLE-like index with overlap counts/Seqpare similarity but no p-values; R and Python bindings; last commit 2021-09. |
| superintervals | mixed | library | pip | – | 1.0.3 2026-10 | yes | 25% | 0.60 | Single-chromosome interval index with SIMD counting; merge_overlaps/union/intersection/difference/symmetric_difference/gaps/span; R binding in a separate repo; README-level docs only. |
| tabix / htslib | C/C++ | mixed | conda | 2011 | 1.24 2026-07 | yes | 25% | 0.90 | Random-access region index, not an arithmetic toolkit; faidx sequence fetch is the same library; used by pysam, bedr, giggle. |
| AIList | C/C++ | mixed | src | 2019 | – | no | 20% | 0.40 | Reference C implementation + CLI from the Sheffield lab; last push 2022; the Python &#x27;ailist&#x27; package below is an independent reimplementation. |
| biobear | mixed | library | pip | – | 0.23.7 2025-01 | no | 20% | 0.78 | File-format-to-Arrow/SQL reader (BED, GFF, VCF, BAM, FASTA...) rather than an interval engine; overlaps only via SQL joins; repo quiet since 2025-03. |
| cgranges | C/C++ | library | src | – | – | no | 20% | 0.30 | Heng Li&#x27;s overlap-query library (C API + IITree.h); in-repo Cython wrapper but no PyPI package; underlies bedtk; last push 2024-05, no tagged releases. |
| intervaltree | Python | library | pip · conda | – | 3.2.1 2025-12 | yes | 20% | 0.60 | Generic (non-genomic) per-interval objects; set ops are set-of-intervals semantics (union/difference of trees), merge_overlaps, slice; no chromosome notion. |
| ncls | Python | library | pip · conda | 2020 | 0.0.70 2025-07 | no | 20% | 0.50 | Index-only library (all_overlaps_both on numpy arrays) that powers pyranges 0.x; last push 2025-07. |
| coitrees | Rust | library | src | – | 0.4.0 2024-01 | no | 15% | 0.60 | Overlap-query engine used by polars-bio and granges; bed-intersect example binary; no Python binding. |
| pybedlite | Python | library | pip · conda | – | 1.1.0 2025-09 | yes | 15% | 0.80 | Deliberately lightweight BED read/write + overlap detection for pipelines; no arithmetic beyond overlaps. |
| pyfaidx | Python | library | pip · conda | 2017 | 0.9 2025 | yes | 15% | 0.90 | Random access to FASTA by .fai index; its `faidx` CLI extracts BED regions (bedtools getfasta&#x27;s library-side counterpart). |


## Gene annotation
{: .sec-navy #annotation }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Gene models from GTF / GFF3 (timed); **C2** Derived regions: promoters, introns, UTRs, intergenic; **C3** Region -&gt; feature class with priority rules (timed); **C4** Nearest TSS / gene with signed distance (timed); **C5** All genes within a window (multi-gene assignment) (timed); **C6** Region-gene-set enrichment (regulatory domains); **C7** Custom promoter / regulatory windows (timed); **C8** Transcript- vs gene-level annotation; **C9** Summary plots (feature mix, distance to TSS); **C10** GO / pathway enrichment of the assigned genes.

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| **genomeblocks** | best | 0.50 | 5/10 | 1.00 | C3 1.00, C4 1.00, C7 1.00 |
| GenomicFeatures + txdbmaker | good | 0.45 | 5/10 | 0.90 | C1 0.79, C3 0.95, C7 0.95 |
| ChIPseeker | good | 0.34 | 8/10 | 0.42 | C3 0.49, C4 0.29, C7 0.49 |
| pyranges (0.x) | good | 0.32 | 4/9 (1 unknown) | 0.73 | C5 0.73 |
| pyranges (1.x) | good | 0.31 | 3/8 (2 unknown) | 0.83 | C5 0.83 |
| valr | moderate | 0.27 | 3/10 | 0.91 | C5 0.91 |
| ChIPpeakAnno | moderate | 0.18 | 10/10 | 0.18 | C4 0.14, C5 0.22 |
| bedtools | low | 0.14 | 2/10 | 0.72 | C4 0.75, C5 0.69 |
| rtracklayer | low | 0.10 | 1/10 | 1.00 | C1 1.00 |
| HOMER | low | 0.01 | 9/9 (1 unknown) | 0.01 | C3 0.01, C4 0.00, C7 0.01 |
| UROPA | low | 0.00 | 6/10 | 0.00 | C4 0.00 |
| GenomicRanges + IRanges | untested | – | 2/2 (8 unknown) | 0.76 | C4 0.53, C5 1.00 |
| plyranges | untested | – | 1/1 (9 unknown) | 0.94 | C5 0.94 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| plyranges | 0.98 | 1.00 | 0.94 | – | 1.00 | 1 | 0 |
| **genomeblocks** | 0.97 | 1.00 | 1.00 | 1.00 | 0.90 | 4 | 2 |
| valr | 0.97 | 1.00 | 0.91 | – | 1.00 | 1 | 0 |
| GenomicFeatures + txdbmaker | 0.97 | 1.00 | 0.87 | 1.00 | 1.00 | 2 | 2 |
| GenomicRanges + IRanges | 0.93 | 1.00 | 0.76 | 1.00 | 1.00 | 2 | 2 |
| pyranges (1.x) | 0.91 | 1.00 | 0.83 | – | 0.90 | 1 | 0 |
| pyranges (0.x) | 0.90 | 1.00 | 0.73 | – | 1.00 | 1 | 0 |
| bedtools | 0.86 | 1.00 | 0.72 | – | 0.90 | 2 | 0 |
| ChIPseeker | 0.79 | 1.00 | 0.39 | 1.00 | 1.00 | 2 | 1 |
| ChIPpeakAnno | 0.65 | 1.00 | 0.18 | 1.00 | 1.00 | 2 | 1 |
| HOMER | 0.18 | 0.96 | 0.01 | – | 0.88 | 2 | 0 |
| UROPA | 0.08 | 1.00 | 0.00 | – | 0.50 | 1 | 0 |
| rtracklayer | – | – | 1.00 | – | 1.00 | 1 | 0 |

### `genes_in_window`
{: .no_toc }

Scored at n = 100k; fastest 39 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| GenomicRanges findOverlaps (maxgap) | GenomicRanges + IRanges | 100k | 39 ms | 1.00 | identical (66,501) | landscape |
| plyranges join_overlap_inner (maxgap) | plyranges | 100k | 59 ms | 1.00 | identical (66,501) | landscape |
| valr bed_window (both) | valr | 100k | 71 ms | 1.00 | differs: 66,500 vs 66,501 (Jaccard 0.9999) | landscape |
| pyranges 1.5.0 join_overlaps (slack) | pyranges (1.x) | 100k | 124 ms | 1.00 | identical (66,501) | landscape |
| pyranges 0.1.4 join (slack) | pyranges (0.x) | 100k | 254 ms | 1.00 | identical (66,501) | landscape |
| bedtools window -w 50 kb (CLI) | bedtools | 100k | 334 ms | 1.00 | reference (66,501 peak-gene pairs) | landscape |
| ChIPpeakAnno annotatePeakInBatch (overlapping, maxgap) | ChIPpeakAnno | 100k | 8.3 s | 1.00 | identical (66,501) | landscape |

### `isoform_selection`
{: .no_toc }

Scored at n = 70k; fastest 40 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks select_isoforms (peaks)** | genomeblocks | 70k | 40 ms | – |  | genes |
| **genomeblocks select_isoforms (peaks + bigWig)** | genomeblocks | 70k | 292 ms | – |  | genes |
| **genomeblocks select_isoforms (bigWig only)** | genomeblocks | 70k | 1.3 s | – |  | genes |

### `labels`
{: .no_toc }

Scored at n = 100k; fastest 98 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks labels** | genomeblocks | 100k | 98 ms | – |  | genes |
| GenomicFeatures DIY labels (transcript-level TSS) | GenomicFeatures + txdbmaker | 100k | 140 ms | 1.00 | reference for the classes (transcript-level TSS windows) | landscape |
| GenomicFeatures DIY labels (gene-level TSS) | GenomicFeatures + txdbmaker | 100k | 146 ms | 1.00 | reference for the classes (gene-level TSS windows) | landscape |
| **genomeblocks labels (ncls backend)** | genomeblocks | 100k | 149 ms | – |  | genes |
| **genomeblocks labels** | genomeblocks | 100k | 152 ms | 1.00 | 99,845/100,000 rows: same class as GenomicFeatures (gene-level TSS windows, same precedence) | landscape |
| **genomeblocks labels (bioframe backend)** | genomeblocks | 100k | 416 ms | – |  | genes |
| **genomeblocks labels (cgranges backend)** | genomeblocks | 100k | 461 ms | – |  | genes |
| **genomeblocks labels (pyranges backend)** | genomeblocks | 100k | 649 ms | – |  | genes |
| **genomeblocks labels (bedtools backend)** | genomeblocks | 100k | 2.6 s | – |  | genes |
| ChIPseeker annotatePeak | ChIPseeker | 100k | 3.3 s | 1.00 | 100,000/100,000 rows: same class as GenomicFeatures (transcript-level TSS windows) | landscape |
| HOMER annotatePeaks.pl -gtf (CLI) | HOMER | 100k | 89 s | 0.92 | 92,082/100,000 rows: same class as GenomicFeatures (transcript-level; HOMER&#x27;s promoter is -1 kb..+100 bp) | landscape |
| **genomeblocks labels** | genomeblocks | 10k | 17 ms | – |  | genes |
| **genomeblocks labels (ncls backend)** | genomeblocks | 10k | 52 ms | – |  | genes |
| **genomeblocks labels (cgranges backend)** | genomeblocks | 10k | 168 ms | – |  | genes |
| **genomeblocks labels (bioframe backend)** | genomeblocks | 10k | 194 ms | – |  | genes |
| **genomeblocks labels (pyranges backend)** | genomeblocks | 10k | 426 ms | – |  | genes |
| **genomeblocks labels (bedtools backend)** | genomeblocks | 10k | 1.4 s | – |  | genes |

### `nearest_tss`
{: .no_toc }

Scored at n = 100k; fastest 25 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks nearest_tss** | genomeblocks | 100k | 25 ms | 1.00 | 100,000/100,000 rows agree with bedtools closest (gene-level TSS) | landscape |
| **genomeblocks nearest_genes** | genomeblocks | 100k | 92 ms | – |  | genes |
| bedtools closest (gene-level TSS windows, CLI) | bedtools | 100k | 145 ms | 1.00 | reference (20,000 windows) | landscape |
| bedtools closest (transcript-level TSS windows, CLI) | bedtools | 100k | 208 ms | 1.00 | reference (70,028 windows) | landscape |
| GenomicRanges distanceToNearest (transcript-level TSS windows) | GenomicRanges + IRanges | 100k | 661 ms | 1.00 | 99,896/100,000 rows agree (transcript-level) | landscape |
| GenomicRanges distanceToNearest (gene-level TSS windows) | GenomicRanges + IRanges | 100k | 690 ms | 1.00 | 99,943/100,000 rows agree (gene-level) | landscape |
| ChIPseeker annotatePeak | ChIPseeker | 100k | 3.3 s | 1.00 | 99,896/100,000 rows: same distance to the nearest TSS window (transcript-level) | landscape |
| ChIPpeakAnno annotatePeakInBatch | ChIPpeakAnno | 100k | 9.6 s | 1.00 | 99,914/100,018 rows: same distance to the nearest TSS window (transcript-level) | landscape |
| HOMER annotatePeaks.pl -gtf (CLI) | HOMER | 100k | 89 s | 1.00 | 100,000/100,000 rows: same distance to the nearest TSS window (transcript-level) | landscape |
| UROPA (CLI) | UROPA | 100k | 4.3 min | 1.00 | 99,999/100,000 rows: same distance to the nearest TSS window (transcript-level) | landscape |
| **genomeblocks nearest_genes** | genomeblocks | 10k | 16 ms | – |  | genes |

### `parse_gtf`
{: .no_toc }

Scored at n = 100k; fastest 3.6 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Genes.make** | genomeblocks | 877k | 609 ms | – |  | genes |
| **genomeblocks Genes.make** | genomeblocks | 877k | 1.0 s | – |  | landscape |
| rtracklayer::import | rtracklayer | 100k | 3.6 s | – |  | landscape |
| txdbmaker::makeTxDbFromGFF | GenomicFeatures + txdbmaker | 100k | 16 s | – |  | landscape |

### Perturbation tests
{: .no_toc }

| Tool | Engine | task | variant | identical answer |
|---|---|---|---|---|
| genomeblocks | genomeblocks | labels | shuffled | yes |
| genomeblocks | genomeblocks | nearest_tss | shuffled | yes |
| ChIPpeakAnno | ChIPpeakAnno annotatePeakInBatch | nearest_tss | shuffled | yes |
| ChIPseeker | ChIPseeker annotatePeak | labels | shuffled | yes |
| GenomicFeatures + txdbmaker | GenomicFeatures DIY labels (gene-level TSS) | labels | shuffled | yes |
| GenomicFeatures + txdbmaker | GenomicFeatures DIY labels (transcript-level TSS) | labels | shuffled | yes |
| GenomicRanges + IRanges | GenomicRanges distanceToNearest (gene-level TSS windows) | nearest_tss | shuffled | yes |
| GenomicRanges + IRanges | GenomicRanges distanceToNearest (transcript-level TSS windows) | nearest_tss | shuffled | yes |

### Every tool in this theme
{: .no_toc }

Capability checklist: `parse_gtf`, `parse_genepred`, `feature_labels`, `nearest_tss`, `nearest_gene`, `custom_promoter_window`, `isoform_selection`, `tss_table`, `export_bed12`, `export_gtf`, `region_gene_enrichment`, `distance_plots`, `gtf_parse_speed_claim`, `uses_txdb`, `dataframe_output`, `derived_gene_regions`, `genes_in_window`, `transcript_level_annotation`, `gene_set_function_enrichment`, `gene_id_mapping`, `gene_filters`, `custom_annotation_sets`, `chrom_naming_conversion`, `sequence_context_gc_cpg`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 54% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| ChIPpeakAnno | R | library | Bioc | 2010 | 3.46.1 2026-04 | yes | 67% | 1.00 | annotatePeakInBatch (nearest TSS/gene, bindingRegion windows), assignChromosomeRegion / genomicElementDistribution pies, getEnrichedGO/PATH, binOverFeature; toGRanges reads BED/GFF/narrowPeak and TxDb/EnsDb. |
| ChIPseeker | R | library | Bioc | 2015 | 1.48.0 2026-04 | yes | 54% | 1.00 | Most-cited R peak annotator: tssRegion-configurable promoter class, plotAnnoPie/Bar, plotDistToTSS, covplot, seq2gene + clusterProfiler enrichment; GTF must first become a TxDb (txdbmaker); 2022 Current Protocols update (10.1002/cpz1.585); Profiles/heatmaps are of peak occupancy (optionally weighted), not of bigWig signal; primarily a peak annotator; Adjacent: covplot draws genome-wide peak coverage per chromosome and TSS profiles; peak annotation, not a locus track viewer. |
| HOMER | Perl | cli | conda | 2010 | 5.1 | yes | 54% | 0.88 | Fixed promoter-TSS class (-1 kb/+100 bp), nearest-TSS distance, -gtf custom models, -go/-genomeOntology enrichment, -hist/-annStats tables (no plotting); version 5.1 taken from the bioconda recipe (built 2026-04) because homer.ucsd.edu refused connections on 2026-10-22; No public VCS (homer.ucsd.edu, HTTP-only; version taken from bioconda); -ghist writes a matrix but plotting is external; motif thresholds are log-odds, GC-matched background enrichment. |
| GenomicFeatures + txdbmaker | R | library | Bioc | 2013 | 1.64.0 (txdbmaker 1.8.0) 2026-04 | yes | 46% | 1.00 | The gene-model substrate most R annotators depend on; makeTxDbFromGFF/UCSC/Ensembl/Biomart moved to txdbmaker in Bioc 3.19; genePred only via live UCSC table queries; BED12/GTF export through rtracklayer&#x27;s asBED/asGFF. |
| genomation | R | library | Bioc | 2015 | 1.44.0 2026-04 | no | 42% | 0.90 | annotateWithGeneParts / getAssociationWithTSS / plotTargetAnnotation working from BED12 gene models (no TxDb needed) plus ScoreMatrix heatmaps; maintenance-only since 2025-04; Last GitHub push 2025-04; the 2026-04 date is the Bioconductor 3.23 version bump. |
| pygtftk | Python | mixed | pip · conda | 2019 | 1.6.2 2022-03 | yes | 42% | 1.00 | get_5p_3p_coords/tss_dist/closest_genes/intronic/get_tx_seq/convert plus OLOGRAM overlap-enrichment statistics; PyPI release 2022 but repo pushed 2026-01; Travis badge stale. |
| rGREAT | R | mixed | Bioc | 2023 | 2.14.0 2026-04 | yes | 42% | 1.00 | Local GREAT (binomial over basal+extension domains) for any organism/gene-set plus the web-service client; plotRegionGeneAssociations gives distance-to-TSS plots; Shiny viewer; Adjacent (annotation theme): functional enrichment of one region set; last commit 2026-01. |
| annotatr | R | library | Bioc | 2017 | 1.38.3 2026-04 | yes | 38% | 1.00 | Overlap-based multi-label annotation (promoters/5&#x27;UTR/exons/introns/3&#x27;UTR/intergenic, CpG, enhancers, chromatin states, MANE/canonical) with plot_annotation/plot_coannotations; no distance-to-TSS; promoters fixed at 1 kb unless built via build_txdb_annotations. |
| ensembldb | R | library | Bioc | 2019 | 2.36.1 2026-04 | yes | 38% | 1.00 | Ensembl-flavoured alternative to TxDb (ensDbFromGtf/ensDbFromGff, AnnotationHub EnsDb packages, tx support level / canonical filters, genome-transcript-protein coordinate mapping); used by Signac and ChIPpeakAnno. |
| GenomicDistributions | R | library | Bioc | 2022 | 1.20.0 2026-04 | no | 38% | 0.90 | QC-style summaries: calcFeatureDist (TSS distance), calcPartitions (promoterCore/Prox, exon, intron, UTRs, intergenic), chromosome bins, GC; plots for each; repo last pushed 2025-04; Adjacent: QC/characterisation of region sets, no enrichment test; last push 2025-04. |
| pyranges (0.x) | Python | library | pip · conda | 2020 | 0.1.6 2026-10 | yes | 38% | 1.00 | Legacy line still released (0.1.6, 2026-10) but README says v1 &#x27;will soon replace&#x27; it; nb_cpu multiprocessing via ray; stats module (jaccard/forbes/fisher); features.tss()/tes()/introns() from GTF. |
| ArchR | R | library | src | 2021 | v1.0.3 2024-11 | no | 33% | 0.78 | Single-cell ATAC framework; peaks get nearestGene/distToGeneStart/peakType (Promoter/Exonic/Intronic/Distal)/distToTSS/nearestTSS/GC with promoterRegion user-set; GitHub-only install; last push 2025-02; Single-cell framework; hypergeometric motif enrichment in marker peaks; last GitHub push 2025-02. |
| goldmine | R | library | src | 2016 | – | no | 33% | 0.60 | Context annotation (promoter/5&#x27;UTR/exon/intron/3&#x27;UTR/intergenic) + genomic-feature enrichment with data.table output; unmaintained since 2019 (GitHub only, never on CRAN/Bioc); the NAR 2016 DOI was not re-verified so left null. |
| greatpy | Python | library | pip | – | 0.0.1 2022-08 | no | 33% | 0.80 | Only Python GREAT reimplementation found (basal+extension/OneClosest/TwoClosest, binomial + hypergeometric, dotplots); 0.0.1 and no commits since 2022-10. |
| pyranges (1.x) | Python | mixed | pip | 2025 | 1.5.0 2026-10 | yes | 33% | 0.90 | Complete rewrite published as package &#x27;pyranges1&#x27; (requires Python&gt;=3.12); set_num_threads for the Rust engine; read_chain exists but liftover not verified; also a &#x27;pyranger&#x27; CLI; preprint DOI copied verbatim from the README. |
| UROPA | Python | cli | pip · conda | 2017 | 4.0.3 2023-10 | no | 33% | 0.50 | Any GTF feature type with user anchors/distances and prioritised multi-query annotation (allhits/finalhits/besthits TSV + summary PDF); last release 2023, repo last pushed 2023-12. |
| chipenrich | R | library | Bioc | 2014 | 2.36.0 2026-04 | no | 29% | 0.90 | Gene-set enrichment with length/mappability correction (chipenrich, polyenrich, broadenrich, hybridenrich, proxReg); 10 genomes via chipenrich.data; plot_dist_to_tss; repo last pushed 2025-06. |
| GREAT (web) | Web | web | web | 2010 | 4.0.4 | ? | 29% | 0.50 | Web-only, four assemblies; version 4.0.4 is the last documented release (per rGREAT vignettes); great.stanford.edu refused connections on 2026-10-22 so current status is unverified; 2022 update paper Tanigawa et al. PLoS Comput Biol. |
| BiocPy genomicranges / iranges | Python | library | pip | – | 0.9.0 2026-09 | yes | 25% | 0.90 | Genentech&#x27;s port of the Bioconductor API (intra/inter-range, set ops, tiles, coverage, find_overlaps/nearest, promoters); read_gtf/read_ucsc; from_pandas/from_polars; BioC2023 abstract, no paper. |
| PAVIS | Web | web | web | 2013 | – | no | 25% | 0.20 | NIEHS web tool; its own update log ends in 2018 and the server returned HTTP 503 on 2026-10-22, so treat as effectively retired. |
| valr | R | library | CRAN | 2017 | 0.10.0 2026-07 | yes | 25% | 1.00 | bedtools semantics in the tidyverse (bed_intersect/closest/merge/slop/makewindows/jaccard/fisher/projection); create_tss/create_introns/create_utrs from BED12; read_gtf is deprecated. |
| AGAT | Perl | cli | conda | – | v1.7.0 2026-04 | yes | 21% | 0.80 | Annotation-file sanitiser/converter (gff2bed, gxf2gxf, add_introns, manage_UTRs, extract_sequences, statistics); no region annotation; README says the paper is not yet published. |
| bedtools | C/C++ | cli | conda | 2010 | v2.31.1 2023-11 | yes | 21% | 0.90 | Reference CLI for genome arithmetic; set ops limited to subtract/complement (union = cat|sort|merge); &#x27;closest&#x27; against a TSS BED is the annotation use; last tagged release 2023 but repo pushed 2026-06; genomecov -bg -scale is the classic BAM-&gt;bedGraph route; GitHub push 2026-06 but no release since 2023-11; pybedtools wraps it; Primitives only: fisher (2x2 on one pair), jaccard, reldist, shuffle (-incl/-excl backgrounds) and pairtobed -type either/both for loop anchors; collection ranking must be scripted. |
| gffread | C/C++ | cli | conda | 2020 | v0.12.9 2026-04 | yes | 21% | 0.60 | StringTie/Cufflinks companion: GFF&lt;-&gt;GTF (-T), --bed BED12, -w/-x/-y sequence extraction, region filtering; no region annotation. |
| gffutils | Python | library | pip · conda | – | 0.14 2026-03 | yes | 21% | 0.90 | Hierarchical (children/parents) model with region() lookups, bed12() export, GFF&lt;-&gt;GTF conversion and pybedtools integration; DB creation is the slow step. |
| Signac | R | library | CRAN | 2021 | 1.17.1 2026-04 | yes | 21% | 1.00 | Single-cell framework; peak annotation = ClosestFeature (nearest annotation feature with type exon/cds/utr/gap, deprecated in 1.17.0), GetTSSPositions, LinkPeaks; gene models come from ensembldb EnsDb; FindMotifs: hypergeometric test against GC-matched background peaks; JASPAR via TFBSTools getMatrixSet. |
| GenomeKit | mixed | library | conda | – | 7.7.0 2026-10 | yes | 17% | 0.90 | Prebuilt GENCODE/RefSeq/Ensembl/MANE annotation binaries queried by interval (find_overlapping etc.), dna(), tracks, VCF/BAM access; capability flags for overlap/slop/parse_gtf from memory of the API, docs page did not enumerate them. |
| GTFtools | Python | cli | pip · conda | 2018 | 0.9.0 2022-06 | no | 17% | 0.40 | Source hosted at genemine.org (no canonical GitHub found; a RacconC/gtftools mirror is listed by BioContainers); -t TSS-flank BED with -w up/down window; dormant since 2022. |
| rtracklayer | R | library | Bioc | 2009 | 1.72.0 2026-04 | yes | 17% | 1.00 | The I/O layer of Bioconductor: import(GTF/GFF/BED), export.bed (BED12 via asBED), export.gff, liftOver(chain), indexed &#x27;which=&#x27; queries on bigWig/tabix; no arithmetic of its own; The R bigWig I/O layer every Bioconductor profile package builds on; summary() gives per-region means. |
| polars-bio | mixed | library | pip | 2025 | 0.36.2 2026-10 | yes | 12% | 1.00 | overlap/nearest/count_overlaps/coverage/merge/cluster/complement/subtract with lazy scans, streaming sinks and cloud (S3/GCS) reads; readers for BED/GFF3/GTF/VCF/BAM/FASTQ etc.; no union/liftover. |
| pyensembl | Python | mixed | pip · conda | – | 2.31.1 2026-10 | yes | 12% | 0.88 | Locus/ID/name lookups (genes_at_locus etc.) and sequences for Ensembl releases or custom GTF/FASTA; overlap lookups only, no nearest or feature classes. |
| BCBio.GFF | Python | library | pip · conda | – | 0.7.1 2024-01 | no | 8% | 0.56 | Biopython-ecosystem GFF reader/writer (GFF.parse/GFF.write); no tables, no coordinates ops; lives inside the bcbb monorepo. |
| biobear | mixed | library | pip | – | 0.23.7 2025-01 | no | 8% | 0.78 | File-format-to-Arrow/SQL reader (BED, GFF, VCF, BAM, FASTA...) rather than an interval engine; overlaps only via SQL joins; repo quiet since 2025-03. |
| gffpandas | Python | library | pip · conda | – | 1.2.0 2019-03 | no | 8% | 0.80 | GFF3 only (no GTF), writes GFF3/CSV/TSV; dormant since 2019 release / 2024 push. |
| gtfparse | Python | library | pip · conda | – | 3.0.5 2026-10 | yes | 8% | 0.60 | Parser only (read_gtf), pyarrow/pandas based with a polars extra; the GTF engine behind pyensembl. |
| HTSeq GenomicArray | Python | library | pip · conda | 2015 | 2.1.2 2026-02 | yes | 4% | 1.00 | Interval-indexed arrays for counting/coverage (htseq-count) plus a streaming GFF_Reader; feature classes must be built by the user; HTSeq 2.0 paper 2022 (btac166). |


## Signal (bigWig / BAM)
{: .sec-navy #signal }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** BAM -&gt; normalised coverage bigWig (CPM / RPGC, read extension) (timed); **C2** Region x bin matrix (reference-point / scale-regions) (timed); **C3** Per-region summary over many tracks (timed); **C4** Read counts per region across many BAMs (timed); **C5** Profile + heatmap plots (sorting, k-means); **C6** Sample correlation / PCA / fingerprint QC; **C7** Track comparison (log2 ratio / subtract vs input) (timed); **C8** Normalisation across samples (scale factors, TMM, spike-in) (timed); **C9** Write bigWig (from bedGraph / arrays) (timed); **C10** Random-access bigWig reading, zoom levels (timed).

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| WiggleTools | best | 0.60 | 6/10 | 1.00 | C7 1.00 |
| megadepth | good | 0.50 | 5/10 | 1.00 | C1 1.00 |
| deepTools | good | 0.49 | 10/10 | 0.49 | C1 0.71, C4 0.26, C7 0.51 |
| **genomeblocks** | good | 0.40 | 4/10 | 1.00 | C2 1.00, C3 1.00, C8 1.00 |
| pybigtools | good | 0.39 | 4/10 | 0.97 | C9 0.97 |
| csaw | good | 0.36 | 4/9 (1 unknown) | 0.82 | C8 0.82 |
| pyBigWig | moderate | 0.33 | 4/10 | 0.81 | C9 0.81 |
| UCSC kent utilities | moderate | 0.32 | 4/10 | 0.79 | C3 0.66, C9 0.92 |
| bedtools | moderate | 0.30 | 4/10 | 0.76 | C1 0.84, C4 0.67 |
| rtracklayer | moderate | 0.30 | 3/10 | 1.00 | C9 1.00 |
| genomation | moderate | 0.29 | 4/9 (1 unknown) | 0.66 | C2 0.66 |
| pysam | moderate | 0.28 | 3/10 | 0.93 | C4 0.93 |
| EnrichedHeatmap | low | 0.15 | 3/10 | 0.50 | C2 0.50 |
| edgeR | untested | – | 1/1 (9 unknown) | 0.77 | C8 0.77 |
| GenomicAlignments | untested | – | 2/4 (6 unknown) | 0.98 | C4 0.98 |
| RNAnorm | untested | – | 1/1 (9 unknown) | 0.49 | C8 0.49 |
| Rsubread featureCounts | untested | – | 1/1 (9 unknown) | 1.00 | C4 1.00 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| rtracklayer | 1.00 | 1.00 | 1.00 | – | 1.00 | 2 | 0 |
| megadepth | 1.00 | 1.00 | 1.00 | – | 1.00 | 1 | 0 |
| Rsubread featureCounts | 1.00 | 1.00 | 1.00 | – | 1.00 | 1 | 0 |
| GenomicAlignments | 0.99 | 1.00 | 0.98 | – | 1.00 | 1 | 0 |
| pybigtools | 0.99 | 1.00 | 0.97 | – | 1.00 | 2 | 0 |
| **genomeblocks** | 0.97 | 1.00 | 1.00 | 1.00 | 0.90 | 5 | 1 |
| pysam | 0.94 | 1.00 | 0.93 | – | 0.90 | 1 | 0 |
| csaw | 0.94 | 0.99 | 0.83 | – | 1.00 | 1 | 0 |
| edgeR | 0.92 | 1.00 | 0.77 | – | 1.00 | 1 | 0 |
| genomation | 0.88 | 1.00 | 0.66 | 1.00 | 0.90 | 1 | 1 |
| WiggleTools | 0.87 | 1.00 | 1.00 | – | 0.67 | 1 | 0 |
| bedtools | 0.87 | 1.00 | 0.74 | – | 0.90 | 3 | 0 |
| EnrichedHeatmap | 0.84 | 1.00 | 0.50 | 1.00 | 1.00 | 1 | 1 |
| pyBigWig | 0.83 | 1.00 | 0.81 | – | 0.70 | 2 | 0 |
| deepTools | 0.78 | 1.00 | 0.47 | – | 1.00 | 5 | 0 |
| UCSC kent utilities | 0.73 | 1.00 | 0.79 | – | 0.50 | 2 | 0 |
| RNAnorm | 0.73 | 1.00 | 0.49 | – | 0.80 | 1 | 0 |

### `bam_to_bigwig`
{: .no_toc }

Scored at n = 10M; fastest 28 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| megadepth --bigwig (CLI) | megadepth | 10M | 28 s | 1.00 | max |Δ| 0.00e+00 | landscape |
| bedtools genomecov -bg + UCSC bedGraphToBigWig (CLI) | bedtools | 10M | 83 s | 1.00 | max |Δ| 0.00e+00 | landscape |
| deepTools bamCoverage --binSize 50 (CLI) | deepTools | 10M | 3.4 min | 1.00 | max |Δ| 0.00e+00 | landscape |
| deepTools bamCoverage --binSize 1 (CLI) | deepTools | 10M | 44.8 min | 1.00 | max |Δ| 0.00e+00 | landscape |
| megadepth --bigwig (CLI) | megadepth | 1M | 6.9 s | 1.00 | max |Δ| 0.00e+00 | landscape |
| bedtools genomecov -bg + UCSC bedGraphToBigWig (CLI) | bedtools | 1M | 25 s | 1.00 | max |Δ| 0.00e+00 | landscape |
| deepTools bamCoverage --binSize 50 (CLI) | deepTools | 1M | 90 s | 1.00 | max |Δ| 0.00e+00 | landscape |
| deepTools bamCoverage --binSize 1 (CLI) | deepTools | 1M | 42.5 min | 1.00 | max |Δ| 0.00e+00 | landscape |

### `bigwig_read`
{: .no_toc }

Scored at n = 10k; fastest 2.0 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| rtracklayer::import (bigWig, which=windows) | rtracklayer | 10k | 2.0 s | – |  | landscape |

### `hichip_shortrange_track`
{: .no_toc }

Scored at n = 5.1M; fastest 17 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks hichip** | genomeblocks | 5.1M | 17 s | 1.00 |  | shortrange |
| awk | sort | bedtools genomecov (shell) | bedtools | 5.1M | 2.3 min | 1.00 |  | shortrange |

### `multi_bam_counts`
{: .no_toc }

Scored at n = 10M; fastest 18 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| Subread featureCounts -F SAF -O (CLI) | Rsubread featureCounts | 10M | 18 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| GenomicAlignments summarizeOverlaps (Union, inter.feature = FALSE) | GenomicAlignments | 10M | 21 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| Rsubread featureCounts (SAF, allowMultiOverlap) | Rsubread featureCounts | 10M | 26 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| pysam 0.23.3 AlignmentFile.count | pysam | 10M | 29 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| bedtools multicov (CLI) | bedtools | 10M | 2.8 min | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| deepTools multiBamSummary BED-file (CLI) | deepTools | 10M | 48.9 min | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| Subread featureCounts -F SAF -O (CLI) | Rsubread featureCounts | 1M | 2.1 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| GenomicAlignments summarizeOverlaps (Union, inter.feature = FALSE) | GenomicAlignments | 1M | 2.9 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| Rsubread featureCounts (SAF, allowMultiOverlap) | Rsubread featureCounts | 1M | 3.4 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| pysam 0.23.3 AlignmentFile.count | pysam | 1M | 20 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| bedtools multicov (CLI) | bedtools | 1M | 59 s | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |
| deepTools multiBamSummary BED-file (CLI) | deepTools | 1M | 71.6 min | 1.00 | 100,000/100,000 regions agree in all 4 BAMs | landscape |

### `region_matrix`
{: .no_toc }

Scored at n = 10k; fastest 455 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **pybigtools (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 455 ms | 1.00 |  | backends |
| **genomeblocks signal (pybigtools)** | genomeblocks | 10k | 642 ms | 1.00 | max |Δ| 7.38e-06 vs per-base pyBigWig means | landscape |
| **python (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 3.1 s | 1.00 |  | backends |
| genomation::ScoreMatrixBin | genomation | 10k | 4.8 s | 1.00 | Pearson 0.99796, max |Δ| 9.01 | landscape |
| rtracklayer + EnrichedHeatmap::normalizeToMatrix | EnrichedHeatmap | 10k | 14 s | 1.00 | max |Δ| 6.39e-14 | landscape |
| **pybigwig (through genomeblocks&#x27; backend)** | genomeblocks | 10k | 4.9 min | 1.00 |  | backends |
| **genomeblocks · pybigtools (exact)** | genomeblocks | 5k | 372 ms | 1.00 | max |Δ| 0, Pearson 1.0000 | signal |
| pybigtools values() + numpy binning | pybigtools | 5k | 951 ms | 1.00 | max |Δ| 0, Pearson 1.0000 | signal |
| pyBigWig values() + numpy binning | pyBigWig | 5k | 1.3 s | 1.00 | max |Δ| 5.7e-06, Pearson 1.0000 | signal |
| deepTools computeMatrix (-p 1) | deepTools | 5k | 37 s | 1.00 | max |Δ| 3.8, Pearson 0.9980 | signal |

### `region_matrix_16_tracks_8_workers`
{: .no_toc }

Scored at n = 5k; fastest 1.4 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks · processes** | genomeblocks | 5k | 1.4 s | – |  | signal |
| deepTools computeMatrix | deepTools | 5k | 96 s | – |  | signal |

### `region_mean`
{: .no_toc }

Scored at n = 10k; fastest 602 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks signal (span=True, 1 bin)** | genomeblocks | 10k | 602 ms | 1.00 | max |Δ| 3.48e-06 (within 0.0001 relative) vs per-base means | landscape |
| UCSC bigWigAverageOverBed (CLI) | UCSC kent utilities | 10k | 6.1 s | 1.00 | max |Δ| 4.87e-05 (within 0.0001 relative) vs per-base means (its mean column: covered bases only) | landscape |

### `track_comparison`
{: .no_toc }

Scored at n = 249M; fastest 4.6 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| WiggleTools log 2 ratio of 50 bp bin means (CLI) | WiggleTools | 249M | 4.6 s | 1.00 | max |Δ| 6.43e-07 | landscape |
| deepTools bigwigCompare --binSize 50 (CLI) | deepTools | 249M | 2.3 min | 1.00 | max |Δ| 4.89e-06 | landscape |

### `track_tmm`
{: .no_toc }

Scored at n = 10M; fastest 79 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks signal.tmm** | genomeblocks | 10M | 79 ms | 1.00 | max |Δ| 5.26e-13 (within 1e-06 relative) | landscape |
| csaw normFactors (TMM) | csaw | 10M | 254 ms | 0.97 | Pearson 0.97109, max |Δ| 44.1 | landscape |
| edgeR calcNormFactors + cpm | edgeR | 10M | 372 ms | 1.00 | reference | landscape |
| rnanorm 2.2.0 TMM | RNAnorm | 10M | 2.7 s | 1.00 | max |Δ| 5.40e-13 (within 1e-06 relative) | landscape |
| **genomeblocks signal.tmm** | genomeblocks | 1M | 40 ms | 1.00 | max |Δ| 5.12e-13 (within 1e-06 relative) | landscape |
| csaw normFactors (TMM) | csaw | 1M | 106 ms | 1.00 | Pearson 1.00000, max |Δ| 0.585 | landscape |
| edgeR calcNormFactors + cpm | edgeR | 1M | 170 ms | 1.00 | reference | landscape |
| rnanorm 2.2.0 TMM | RNAnorm | 1M | 2.3 s | 1.00 | max |Δ| 5.12e-13 (within 1e-06 relative) | landscape |

### `write_bigwig`
{: .no_toc }

Scored at n = 10M; fastest 10 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| rtracklayer export.bw | rtracklayer | 10M | 10 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| pybigtools 0.2.5 write | pybigtools | 10M | 13 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| UCSC bedGraphToBigWig (CLI) | UCSC kent utilities | 10M | 17 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| pyBigWig 0.3.24 addEntries | pyBigWig | 10M | 37 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| rtracklayer export.bw | rtracklayer | 1000k | 1.1 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| UCSC bedGraphToBigWig (CLI) | UCSC kent utilities | 1000k | 1.5 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| pybigtools 0.2.5 write | pybigtools | 1000k | 1.9 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |
| pyBigWig 0.3.24 addEntries | pyBigWig | 1000k | 6.3 s | 1.00 | max |Δ| 2.98e-08 (within 1e-06 relative) | landscape |

### Perturbation tests
{: .no_toc }

| Tool | Engine | task | variant | identical answer |
|---|---|---|---|---|
| genomeblocks | genomeblocks | region_matrix | shuffled | yes |
| EnrichedHeatmap | rtracklayer + EnrichedHeatmap::normalizeToMatrix | region_matrix | shuffled | yes |
| genomation | genomation::ScoreMatrixBin | region_matrix | shuffled | yes |

### Every tool in this theme
{: .no_toc }

Capability checklist: `bigwig_read`, `bigwig_write`, `region_matrix`, `region_mean`, `profile_plot`, `heatmap_plot`, `multi_track`, `track_normalisation`, `bam_coverage`, `bam_pileup`, `bedgraph_to_bigwig`, `super_enhancer`, `hichip_shortrange_track`, `parallel`, `strand_aware`, `in_memory_api`, `gpu`, `multi_bam_region_counts`, `sample_correlation_qc`, `track_comparison`, `paired_end_fragments`, `blacklist_masking`, `peak_set_qc`, `signal_smoothing`, `bigwig_merge`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 60% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| deepTools | Python | mixed | pip · conda | 2014 | 4.0.0 2026-09 | yes | 80% | 1.00 | De-facto standard; reference-point and scale-regions matrices, --scaleFactor / bigwigCompare for between-track scaling; 4.0.0 requires Python&gt;=3.12. |
| DiffBind | R | library | Bioc | 2012 | 3.22.2 2026-04 | yes | 52% | 1.00 | dba.normalize: library-size, RLE, TMM, spike-in, background bins; dba.plotProfile delegates matrices/heatmaps to profileplyr. |
| HOMER | Perl | cli | conda | 2010 | 5.1 | yes | 44% | 0.88 | Fixed promoter-TSS class (-1 kb/+100 bp), nearest-TSS distance, -gtf custom models, -go/-genomeOntology enrichment, -hist/-annStats tables (no plotting); version 5.1 taken from the bioconda recipe (built 2026-04) because homer.ucsd.edu refused connections on 2026-10-22; No public VCS (homer.ucsd.edu, HTTP-only; version taken from bioconda); -ghist writes a matrix but plotting is external; motif thresholds are log-odds, GC-matched background enrichment. |
| csaw | R | library | Bioc | 2016 | 1.46.0 2026-04 | yes | 40% | 1.00 | TMM / efficiency / composition normalisation factors between samples; profileSites gives aggregate coverage profiles. |
| WiggleTools | C/C++ | cli | conda | 2014 | 1.2.11 (bioconda) 2025-09 | yes | 40% | 0.67 | mean/sum/scale across tracks, apply over BED regions, profile command (numbers only); last commit 2026-05. |
| genomation | R | library | Bioc | 2015 | 1.44.0 2026-04 | no | 36% | 0.90 | annotateWithGeneParts / getAssociationWithTSS / plotTargetAnnotation working from BED12 gene models (no TxDb needed) plus ScoreMatrix heatmaps; maintenance-only since 2025-04; Last GitHub push 2025-04; the 2026-04 date is the Bioconductor 3.23 version bump. |
| ngs.plot | mixed | cli | src | 2014 | 2.63 (no tagged release) | no | 36% | 0.40 | BAM-only, bundled genome annotation databases; last commit 2023-05; non-commercial licence. |
| TOBIAS | Python | cli | pip · conda | 2020 | 0.17.5 2026-08 | yes | 32% | 0.80 | ATAC footprinting suite; TFBScan uses MOODS with p-value thresholds, BINDetect compares conditions, ClusterMotifs builds motif similarity clusters. |
| trackplot | R | library | src | 2024 | – | no | 32% | 0.50 | Minimal-dependency R script layered on bwtool; IGV-style locus plots plus profiles/heatmaps; GitHub push 2024-10; DOI inferred from the article id vbae031; Minimal-dependency IGV-style locus plots, profile plots and heatmaps from bigWig (needs bwtool); gene models from UCSC/GTF; last commit 2024-10, no tagged releases. Distinct from the Python sashimi tool of the same name. |
| bamsignals | R | library | Bioc | – | 1.44.1 2026-04 | yes | 28% | 0.90 | bamProfile / bamCount / bamCoverage straight from indexed BAM, strand-specific; no plotting; GitHub push 2026-05. |
| bedtools | C/C++ | cli | conda | 2010 | v2.31.1 2023-11 | yes | 28% | 0.90 | Reference CLI for genome arithmetic; set ops limited to subtract/complement (union = cat|sort|merge); &#x27;closest&#x27; against a TSS BED is the annotation use; last tagged release 2023 but repo pushed 2026-06; genomecov -bg -scale is the classic BAM-&gt;bedGraph route; GitHub push 2026-06 but no release since 2023-11; pybedtools wraps it; Primitives only: fisher (2x2 on one pair), jaccard, reldist, shuffle (-incl/-excl backgrounds) and pairtobed -type either/both for loop anchors; collection ranking must be scripted. |
| EnrichedHeatmap | R | library | Bioc | 2018 | 1.42.0 2026-04 | yes | 28% | 1.00 | normalizeToMatrix takes a GRanges with a value column (read bigWig via rtracklayer); heatmap lists concatenate tracks; GitHub push 2026-01. |
| megadepth | C/C++ | mixed | conda · Bioc | 2021 | 1.2.0 (C++) / 1.22.0 (R) 2022-02 | yes | 28% | 1.00 | Fast per-region quantification over thousands of BAMs/bigWigs (recount3 engine); GitHub push 2026-09, no release since 2022. |
| metagene2 | R | library | Bioc | – | 1.28.0 2026-04 | yes | 28% | 0.90 | Input/control normalisation and design groups; activity verified only through the Bioconductor bump. |
| profileplyr | R | library | Bioc | – | 1.28.3 2026-04 | yes | 28% | 0.90 | Imports deepTools computeMatrix output or builds matrices via soGGi, plots with EnrichedHeatmap; range annotation and grouping. |
| seqplots | R | mixed | src | 2016 | 1.27.0 (Bioconductor 3.12, last before removal) | no | 28% | 0.67 | Removed from Bioconductor; last commit 2022-05; also plots motif-density tracks alongside signal. |
| soGGi | R | library | Bioc | – | 1.44.0 2025-10 | no | 28% | 0.67 | Deprecated and removed from Bioconductor 3.23; last version in 3.22; backend of profileplyr. |
| bigwig-loader | Python | library | pip · conda | 2024 | 0.3.2 2025-12 | yes | 24% | 0.80 | GPU (CUDA) bigWig decompression for ML dataloading; the only GPU region-cube reader found; conda-forge availability not verified. |
| bwtool | C/C++ | cli | src | 2014 | – | no | 24% | 0.50 | matrix / aggregate / summary subcommands output numbers only (no plots); last commit 2017-09; used as the engine by trackplot. |
| CoverageView | R | library | Bioc | – | 1.50.0 2026-04 | yes | 24% | 0.89 | Small single-track package; activity is the Bioconductor version bump only (no public GitHub). |
| pybigtools | Rust | mixed | pip · conda | 2024 | pybigtools 0.3.0 / bigtools v0.5.8 2026-06 | yes | 24% | 1.00 | Rust reimplementation of UCSC kent bigWig/bigBed I/O with multithreaded writing; average_over_bed in the Python API. |
| samtools | C/C++ | cli | conda | 2009 | 1.24 2026-07 | yes | 24% | 0.80 | mpileup is the canonical per-base allele/base-count tool; no bigWig support. |
| UCSC kent utilities | C/C++ | cli | conda | 2010 | 502 (bioconda ucsc-bigwigaverageoverbed) 2026-10 | yes | 24% | 0.50 | Reference implementation of the bigWig format; one-region-one-number tools, no matrices or plots. |
| ChIPseeker | R | library | Bioc | 2015 | 1.48.0 2026-04 | yes | 20% | 1.00 | Most-cited R peak annotator: tssRegion-configurable promoter class, plotAnnoPie/Bar, plotDistToTSS, covplot, seq2gene + clusterProfiler enrichment; GTF must first become a TxDb (txdbmaker); 2022 Current Protocols update (10.1002/cpz1.585); Profiles/heatmaps are of peak occupancy (optionally weighted), not of bigWig signal; primarily a peak annotator; Adjacent: covplot draws genome-wide peak coverage per chromosome and TSS profiles; peak annotation, not a locus track viewer. |
| pyBigWig | Python | library | pip · conda | – | 0.3.26 2026-09 | yes | 20% | 0.70 | Low-level reader/writer used by deepTools; remote http access; threaded calls serialise on the GIL. |
| pysam | Python | library | pip · conda | – | 0.24.1 2026-09 | yes | 20% | 0.90 | The Python BAM layer: count_coverage gives per-base ACGT counts in memory; no bigWig. |
| rtracklayer | R | library | Bioc | 2009 | 1.72.0 2026-04 | yes | 20% | 1.00 | The I/O layer of Bioconductor: import(GTF/GFF/BED), export.bed (BED12 via asBED), export.gff, liftOver(chain), indexed &#x27;which=&#x27; queries on bigWig/tabix; no arithmetic of its own; The R bigWig I/O layer every Bioconductor profile package builds on; summary() gives per-region means. |
| mosdepth | C/C++ | cli | conda | 2018 | v0.3.14 2026-04 | yes | 16% | 0.70 | Written in Nim; --by BED gives mean depth per region; bedGraph output needs bedGraphToBigWig for bigWig. |
| GenomicAlignments | R | library | conda · Bioc | 2013 | – | yes | 12% | 1.00 | Reads BAM alignments into R; summarizeOverlaps counts reads per feature (Union / IntersectStrict / IntersectNotEmpty modes). Added 2026-10-23 as an engine the landscape benchmark runs. |
| LILY | R | cli | src | 2017 | – | no | 12% | 0.22 | ROSE-style stitching with explicit copy-number correction for cancer samples; depends on HMCan; last commit 2018-08. |
| MAnorm2 | R | library | CRAN | 2021 | 1.2.2 2022-10 | no | 12% | 0.67 | Between-sample MA normalisation using common peaks; counting is done by the separate MAnorm2_utils Python package; last commit 2022-11. |
| ROSE | Python | cli | src | 2013 | 1.3.2 2024-05 | yes | 12% | 0.50 | Stitching (12.5 kb) + TSS exclusion + input-subtracted rpm ranking + tangent-line knee; GitHub push 2025-11; Docker image provided; original Young-lab code was Python 2; Adjacent: SE calling only; stjude mirror of the Young-lab Bitbucket code (README commit 2025-11), Python 2 era. |
| ROSE2 | Python | cli | pip | 2013 | 1.1.0 2019-10 | no | 12% | 0.43 | Requires bamliquidator; last commit 2019-10; licence not stated on PyPI or GitHub. |
| CREAM | R | library | CRAN | 2019 | 1.1.1 2018-06 | no | 8% | 0.67 | Peak-only alternative to ROSE (no signal ranking, data-driven stitching window); last commit 2021-02. |
| FitHiChIP | mixed | cli | src | 2019 | 11.0 2022-12 | no | 4% | 0.50 | PeakInferHiChIP.sh builds 1D coverage from HiC-Pro valid pairs (and DE/SC/RE reads when present) and calls peaks with MACS2; last commit 2024-04; v11 (2022-12) added .hic/.cool input; last commit 2024-04; no bioconda package found. |
| HiC-Pro | mixed | cli | src | 2015 | v3.1.0 2021-08 | no | 4% | 0.43 | Produces the allValidPairs that hichipper / HiChIP-Peaks / FitHiChIP consume; no 1D track itself; last commit 2022-06; allValidPairs is the de-facto HiChIP input for FitHiChIP/hichipper; last commit 2022-06; no bioconda package found under hic-pro/hicpro. |
| HiChIP-Peaks | Python | cli | pip | 2020 | 0.1.2 2020-04 | no | 4% | 0.50 | Restriction-site-centred model of short-range HiChIP reads; last commit 2021-01. |
| hichipper | Python | cli | pip | 2018 | 0.7.7 2019-06 | no | 4% | 0.62 | Calls peaks directly from HiChIP reads with a restriction-site background model; last commit 2022-07; Unmaintained since 2022-07. |
| imPROSE | Python | cli | src | 2019 | – | no | 4% | 0.22 | Predicts SEs from features rather than stitching signal; last commit 2018-03; publication DOI from memory (low confidence). |
| SEanalysis 2.0 | Web | web | web | 2019 | 2.0 | ? | 0% | 0.38 | Annotates user regions against precomputed ROSE-called SEs from SEdb (2.0 paper NAR 2023, PMC10320134); does not call SEs from user signal. |


## Motifs
{: .sec-navy #motifs }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Genome-scale PWM scan with a p-value threshold (timed); **C2** Library-scale scan (1,000+ motifs x 100k regions) (timed); **C3** Known-motif enrichment vs background (timed); **C4** Background construction (GC/CpG-matched, shuffled); **C5** De novo motif discovery (timed); **C6** Comparison to motif databases (timed); **C7** Clustering into non-redundant archetypes; **C8** Differential enrichment between two sets (timed); **C9** Per-sample motif activity (deviation scores) (timed); **C10** Positional distribution / centrality around summits.

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| **genomeblocks** | best | 0.60 | 6/10 | 1.00 | C1 1.00, C2 1.00, C3 1.00, C8 1.00 |
| MEME suite | good | 0.51 | 8/10 | 0.64 | C1 0.31, C2 0.22, C3 0.79, C5 0.83, C6 0.90, C8 0.79 |
| HOMER | good | 0.46 | 9/10 | 0.51 | C1 0.33, C3 0.62, C5 0.48, C8 0.62 |
| universalmotif | good | 0.39 | 9/10 | 0.43 | C1 0.36, C6 0.50 |
| chromVAR | moderate | 0.29 | 6/10 | 0.49 | C9 0.49 |
| MOODS | moderate | 0.27 | 3/10 | 0.91 | C1 0.91, C2 0.90 |
| motifmatchr | low | 0.15 | 3/10 | 0.50 | C1 0.50 |
| Biopython | low | 0.05 | 3/10 | 0.18 | C1 0.18 |
| Biostrings | low | 0.00 | 1/10 | 0.00 | C1 0.00 |
| pychromVAR | untested | – | 2/3 (7 unknown) | 1.00 | C9 1.00 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| pychromVAR | 1.00 | 1.00 | 1.00 | – | 1.00 | 1 | 0 |
| **genomeblocks** | 0.97 | 1.00 | 1.00 | – | 0.90 | 5 | 0 |
| MOODS | 0.80 | 1.00 | 0.91 | – | 0.56 | 2 | 0 |
| motifmatchr | 0.79 | 0.69 | 0.73 | – | 1.00 | 1 | 0 |
| chromVAR | 0.79 | 0.50 | 0.97 | – | 1.00 | 1 | 0 |
| universalmotif | 0.78 | 0.80 | 0.60 | – | 1.00 | 2 | 0 |
| HOMER | 0.78 | 0.89 | 0.61 | – | 0.88 | 4 | 0 |
| MEME suite | 0.75 | 0.95 | 0.57 | – | 0.78 | 7 | 0 |
| Biopython | 0.56 | 1.00 | 0.18 | – | 1.00 | 1 | 0 |
| Biostrings | 0.10 | 1.00 | 0.00 | – | 0.90 | 1 | 0 |

### `db_comparison`
{: .no_toc }

Scored at n = 20; fastest 22 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| MEME Tomtom (CLI) | MEME suite | 20 | 22 s | 0.90 | 18/20 queries have their source motif in the top 3 | landscape |
| universalmotif compare_motifs (PCC) | universalmotif | 20 | 66 s | 0.60 | 12/20 queries have their source motif in the top 3 | landscape |

### `de_novo`
{: .no_toc }

Scored at n = 1k; fastest 15 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| MEME STREME (CLI) | MEME suite | 1k | 15 s | 0.75 | best match to X 0.722, to Z 0.778 (5 motifs) | landscape |
| HOMER homer2 denovo (CLI) | HOMER | 1k | 38 s | 0.55 | best match to X 0.556, to Z 0.542 (5 motifs) | landscape |
| MEME meme -mod zoops (CLI) | MEME suite | 1k | 15.4 min | 0.92 | best match to X 0.944, to Z 0.889 (3 motifs) | landscape |

### `differential_enrichment`
{: .no_toc }

Scored at n = 1k; fastest 410 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks scan_motifs_matrix + compare_motifs** | genomeblocks | 1k | 410 ms | 1.00 | ranks: MA0149.1 #1, MA0729.1 #10 (top 5 wanted for MA0149.1; not wanted for MA0729.1) | landscape |
| MEME AME (CLI) | MEME suite | 1k | 1.8 s | 1.00 | ranks: MA0149.1 #1, MA0729.1 #None (top 5 wanted for MA0149.1; not wanted for MA0729.1) | landscape |
| MEME SEA (CLI) | MEME suite | 1k | 5.5 s | 1.00 | ranks: MA0149.1 #1, MA0729.1 #None (top 5 wanted for MA0149.1; not wanted for MA0729.1) | landscape |
| HOMER homer2 known (CLI) | HOMER | 1k | 5.8 s | 1.00 | ranks: MA0149.1 #1, MA0729.1 #17 (top 5 wanted for MA0149.1; not wanted for MA0729.1) | landscape |

### `known_enrichment`
{: .no_toc }

Scored at n = 1k; fastest 414 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks scan_motifs_matrix + compare_motifs_to_ref** | genomeblocks | 1k | 414 ms | 1.00 | ranks: MA0149.1 #1, MA0729.1 #3 (top 5 wanted for MA0149.1, MA0729.1) | landscape |
| MEME AME (CLI) | MEME suite | 1k | 1.8 s | 1.00 | ranks: MA0149.1 #1, MA0729.1 #2 (top 5 wanted for MA0149.1, MA0729.1) | landscape |
| MEME SEA (CLI) | MEME suite | 1k | 5.4 s | 1.00 | ranks: MA0149.1 #1, MA0729.1 #2 (top 5 wanted for MA0149.1, MA0729.1) | landscape |
| HOMER homer2 known (CLI) | HOMER | 1k | 5.9 s | 1.00 | ranks: MA0149.1 #1, MA0729.1 #3 (top 5 wanted for MA0149.1, MA0729.1) | landscape |

### `motif_activity`
{: .no_toc }

Scored at n = 2k; fastest 3.9 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| pychromVAR 0.0.4 | pychromVAR | 2k | 3.9 s | 1.00 | X (MA0149.1) #3 of 100 by B - A t; Y (MA0731.1) #2 from the bottom (top / bottom 3 wanted) | landscape |
| chromVAR matchMotifs + computeDeviations | chromVAR | 2k | 4.6 s | 0.50 | X (MA0149.1) #4 of 100 by B - A t; Y (MA0731.1) #2 from the bottom (top / bottom 3 wanted) | landscape |

### `scan_library_1019`
{: .no_toc }

Scored at n = 1k; fastest 368 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks scan_motifs_matrix (lightmotif)** | genomeblocks | 1k | 368 ms | – |  | motifs |
| **genomeblocks scan_motifs_matrix (MOODS)** | genomeblocks | 1k | 565 ms | – |  | motifs |
| MOODS (C++, all motifs per pass) | MOODS | 1k | 721 ms | – |  | motifs |
| MEME FIMO --text (CLI) | MEME suite | 1k | 78 s | – |  | motifs |

### `scan_logodds`
{: .no_toc }

Scored at n = 1k; fastest 37 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **lightmotif (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 37 ms | 1.00 |  | backends |
| **genomeblocks scan_motifs_matrix (lightmotif)** | genomeblocks | 1k | 42 ms | 1.00 | 993 hits vs 993, Pearson 1.0000 | motifs |
| **moods (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 53 ms | 1.00 |  | backends |
| **genomeblocks scan_motifs_matrix (lightmotif)** | genomeblocks | 1k | 55 ms | 1.00 | identical: 993 hits vs the exact numpy scan | landscape |
| **genomeblocks scan_motifs_matrix (MOODS)** | genomeblocks | 1k | 58 ms | 1.00 | 993 hits vs 993, Pearson 1.0000 | motifs |
| MOODS (C++, all motifs per pass) | MOODS | 1k | 69 ms | 1.00 | 993 hits vs 993, Pearson 1.0000 | motifs |
| universalmotif::scan_sequences | universalmotif | 1k | 3.2 s | 1.00 | 992 vs 993 hits, Pearson 1.0000 over motifs | landscape |
| **biopython (through genomeblocks&#x27; backend)** | genomeblocks | 1k | 3.3 s | 1.00 |  | backends |
| HOMER homer2 find (CLI) | HOMER | 1k | 3.9 s | 1.00 | 997 vs 993 hits, Pearson 0.9998 over motifs | landscape |
| MEME FIMO --text (CLI) | MEME suite | 1k | 7.6 s | 1.00 | 996 hits vs 993, Pearson 0.9999 | motifs |
| Biopython PSSM.search | Biopython | 1k | 11 s | 1.00 | 993 hits vs 993, Pearson 1.0000 | motifs |
| Biostrings::countPWM | Biostrings | 1k | 62 s | 1.00 | identical: 993 hits | landscape |

### `scan_pvalue`
{: .no_toc }

Scored at n = 1k; fastest 221 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks scan_motifs_matrix (lightmotif, both strands)** | genomeblocks | 1k | 221 ms | 1.00 | 11,355 vs 11,296 hits, Pearson 0.9992 over motifs vs the exact p-value scan | landscape |
| motifmatchr::matchMotifs (p-value cutoff) | motifmatchr | 1k | 1.4 s | 0.69 | 12,379 vs 11,296 hits, Pearson 0.6862 over motifs | landscape |
| MEME FIMO --text, both strands (CLI) | MEME suite | 1k | 16 s | 1.00 | 11,255 vs 11,296 hits, Pearson 0.9999 over motifs vs the exact p-value scan | landscape |

### Every tool in this theme
{: .no_toc }

Capability checklist: `parse_jaspar`, `parse_meme`, `parse_transfac`, `pwm_scan`, `pvalue_threshold`, `per_position_profile`, `enrichment_vs_background`, `differential_enrichment`, `motif_comparison`, `motif_clustering`, `logo_plot`, `de_novo_discovery`, `multi_engine`, `parallel`, `gpu`, `in_memory_api`, `web_only`, `library_scale_scan`, `background_construction`, `per_sample_motif_activity`, `motif_format_export`, `motif_model_conversion`, `motif_hits_genomic`, `motif_cooccurrence_spacing`, `variant_motif_effect`, `masked_scanning`, `attribution_motifs`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 63% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| GimmeMotifs | Python | mixed | pip · conda | 2011 | 0.18.4 2026-06 | yes | 74% | 1.00 | Most complete single Python package: FPR thresholds from genomic background, maelstrom differential across many sets, ensemble de novo (MEME/HOMER/BioProspector...), non-redundant clustered database. |
| universalmotif | R | library | Bioc | 2024 | 1.30.1 2026-04 | yes | 70% | 1.00 | Import/export of JASPAR, MEME, TRANSFAC, HOMER, CIS-BP, UniPROBE; scan_sequences with p-value thresholds, enrich_motifs vs shuffled/background, compare_motifs + motif trees; de novo only via run_meme wrapper. |
| MEME suite | C/C++ | mixed | conda | 2009 | 5.5.9 (bioconda) 2026-05 | yes | 67% | 0.78 | Reference toolkit (FIMO scan, SEA/AME enrichment with control sets, CentriMo positional, Tomtom comparison, STREME/MEME de novo); no public VCS, meme-suite.org TLS failed on 2026-10-22 so version is from bioconda. |
| RSAT | Perl | mixed | conda | 2008 | 2026.10.08 (bioconda rsat-core) 2026-10 | yes | 67% | 0.88 | matrix-clustering is the engine behind JASPAR&#x27;s motif clusters; Markov background models for p-values; 2022 update NAR 10.1093/nar/gkac312. |
| memes | R | library | Bioc | 2021 | 1.20.0 2026-04 | yes | 56% | 1.00 | Needs a local MEME install; wraps FIMO, AME (with control sets), Tomtom, DREME/STREME/MEME but not CentriMo or SEA. |
| monaLisa | R | library | Bioc | 2022 | 1.18.0 2026-04 | yes | 52% | 1.00 | Binned enrichment (HOMER-style known-motif test, Biostrings or HOMER engine) plus randomised-lasso stability selection; GC-matched background bins. |
| HOMER | Perl | cli | conda | 2010 | 5.1 | yes | 48% | 0.88 | Fixed promoter-TSS class (-1 kb/+100 bp), nearest-TSS distance, -gtf custom models, -go/-genomeOntology enrichment, -hist/-annStats tables (no plotting); version 5.1 taken from the bioconda recipe (built 2026-04) because homer.ucsd.edu refused connections on 2026-10-22; No public VCS (homer.ucsd.edu, HTTP-only; version taken from bioconda); -ghist writes a matrix but plotting is external; motif thresholds are log-odds, GC-matched background enrichment. |
| motifTestR | R | library | Bioc | – | 1.8.1 2026-04 | yes | 44% | 0.90 | R-native CentriMo alternative (binomial per bin + harmonic-mean p); enrichment with iterative, quasi-Poisson or hypergeometric models; parsing delegated to universalmotif; no paper found. |
| PWMEnrich | R | library | Bioc | – | 4.48.0 2026-04 | yes | 44% | 0.89 | Threshold-free lognormal affinity model with precomputed genome backgrounds; motifDiffEnrichment for two sets; activity is the Bioconductor bump only. |
| ArchR | R | library | src | 2021 | v1.0.3 2024-11 | no | 41% | 0.78 | Single-cell ATAC framework; peaks get nearestGene/distToGeneStart/peakType (Promoter/Exonic/Intronic/Distal)/distToTSS/nearestTSS/GC with promoterRegion user-set; GitHub-only install; last push 2025-02; Single-cell framework; hypergeometric motif enrichment in marker peaks; last GitHub push 2025-02. |
| Biopython | Python | library | pip · conda | 2009 | 1.88 2026-08 | yes | 41% | 1.00 | Broadest format parser (jaspar, meme, transfac, pfm, sites, xms, clusterbuster...); ScoreDistribution gives FPR/patser thresholds; logo via the WebLogo web service; pure-Python scanning is slow. |
| JASPAR web tools | Web | web | web | 2026 | JASPAR 2026 (11th release) | yes | 41% | 0.60 | 2026 site lists no per-sequence Scan tool; pwm_scan refers to precomputed genome-wide TFBS tracks; enrichment is LOLA-based (universe or set1 vs set2); data licence from the site&#x27;s general terms (not re-verified). |
| Signac | R | library | CRAN | 2021 | 1.17.1 2026-04 | yes | 41% | 1.00 | Single-cell framework; peak annotation = ClosestFeature (nearest annotation feature with type exon/cds/utr/gap, deprecated in 1.17.0), GetTSSPositions, LinkPeaks; gene models come from ensembldb EnsDb; FindMotifs: hypergeometric test against GC-matched background peaks; JASPAR via TFBSTools getMatrixSet. |
| TOBIAS | Python | cli | pip · conda | 2020 | 0.17.5 2026-08 | yes | 41% | 0.80 | ATAC footprinting suite; TFBScan uses MOODS with p-value thresholds, BINDetect compares conditions, ClusterMotifs builds motif similarity clusters. |
| chromVAR | R | library | Bioc | 2017 | 1.34.1 2026-04 | yes | 37% | 1.00 | GC/accessibility-matched background peak sets; scanning via motifmatchr; the single-cell accessibility standard reused by ArchR/Signac. |
| TFBSTools | R | library | Bioc | 2016 | 1.50.0 2026-04 | yes | 37% | 1.00 | Matrix classes used across Bioconductor; searchSeq thresholds are relative scores (min.score), de novo only via the runMEME wrapper. |
| gReLU | Python | library | pip | 2025 | 1.1.0 2026-03 | yes | 33% | 1.00 | Adjacent deep-learning framework; motif scanning delegated to memesuite-lite FIMO, de novo via TF-MoDISco on attributions. |
| tangermeme | Python | library | pip | 2026 | 1.6.0 2026-10 | yes | 33% | 1.00 | Deep-learning interpretation toolkit (adjacent); FIMO/Tomtom CLI moved to memesuite-lite, so scanning keys left null; preprint bioRxiv 10.1101/2025.08.08.669296. |
| lightmotif | Rust | library | pip | – | 0.10.1 2026-03 | yes | 30% | 0.90 | AVX2/SSE2/NEON vectorised scanning (parallel = SIMD, single-threaded); TFM-PVALUE reimplementation for p-value thresholds; no paper or CITATION.cff found. |
| MOODS | C/C++ | library | pip · conda | 2009 | 1.9.4.1 2019-05 | no | 26% | 0.56 | Lookahead-filtration scanner used as the backend of motifmatchr, TOBIAS and the genomeblocks bioconda build; last commit 2019-05. |
| motifStack | R | library | Bioc | 2018 | 1.56.0 2026-04 | yes | 26% | 1.00 | importMatrix reads pfm/pcm/meme/transfac/jaspar/cisbp/homer; motif trees and signatures for visualising clustered libraries. |
| TF-MoDISco-lite | Python | mixed | pip | 2018 | 3.0.0 2026-10 | yes | 26% | 0.62 | Adjacent: de novo motifs from deep-learning attributions, not from sequence sets; arXiv DOI follows the standard arXiv DOI scheme. |
| memesuite-lite | Python | mixed | pip | – | 0.5.0 2026-09 | yes | 22% | 0.62 | Pure-Python drop-in for FIMO and Tomtom without the MEME licence; used by gReLU and tangermeme. |
| motifmatchr | R | library | Bioc | 2017 | 1.34.0 2026-04 | yes | 22% | 1.00 | MOODS C++ backend behind chromVAR, ArchR and Signac; parsing delegated to TFBSTools/universalmotif. |
| oPOSSUM 3 | Web | web | web | 2012 | – | no | 22% | 0.33 | opossum.cisreg.ca did not resolve on 2026-10-22 (apparently offline); a standalone Perl distribution once existed. |
| pycistarget | Python | library | src | 2023 | – | yes | 22% | 0.60 | Not on PyPI (pip from GitHub, Rust extension); last commit 2026-10-22; databases built separately with Cluster-Buster. |
| SnapATAC2 | Python | library | pip · conda | 2024 | 2.10.0 2026-07 | yes | 22% | 1.00 | Rust-backed single-cell framework; hypergeometric (no background) or binomial (with background) motif enrichment over region groups. |
| STAMP | Web | web | web | 2007 | – | no | 22% | 0.22 | Server still answers (2026-10-22) but its databases stop at JASPAR 2010 / TRANSFAC 11.3. |
| MotifDb | R | library | Bioc | – | 1.54.0 2026-04 | yes | 19% | 0.90 | Resource package (JASPAR, HOCOMOCO, cisbp, UniPROBE...) feeding PWMEnrich/universalmotif; activity is the Bioconductor bump only. |
| PWMScan | C/C++ | mixed | web | 2018 | 2021-04 | no | 19% | 0.57 | Web server moved to epd.expasy.org/pwmtools (&#x27;last update April 2021&#x27;); score&lt;-&gt;p-value conversion with background models; standalone C package not verified. |
| pymemesuite | Python | library | pip | – | 0.1.0a4 2025-05 | ? | 19% | 0.67 | Alpha-stage in-process FIMO; GitHub activity after 2025-05 not checked. |
| Biostrings | R | library | conda · Bioc | – | 2.78 2026-04 | yes | 15% | 0.90 | matchPWM / countPWM: one PWM against one sequence at an absolute score; the primitive under motifmatchr and TFBSTools. |
| Pscan | Web | web | web | 2009 | – | ? | 15% | 0.38 | HTTP-only server at 159.149.160.88/pscan; unreachable over HTTPS on 2026-10-22, status unverified. |
| pychromVAR | Python | library | pip | – | 0.0.4 | ? | 15% | 1.00 | Python port of chromVAR (GC-matched background peaks, motif deviation z-scores). Added 2026-10-23 as an engine the landscape benchmark runs. |
| RcisTarget | R | library | Bioc | 2017 | 1.32.0 2026-04 | yes | 15% | 1.00 | Deprecated in Bioconductor 3.23 (removal planned for 3.24); superseded by pycistarget; no scanning, relies on precomputed cisTarget databases. |
| Vierstra motif archetypes | Python | mixed | src | 2020 | 2.1-beta (README) | no | 15% | 0.44 | Resource/notebook workflow rather than a package: Tomtom-similarity clustering of CIS-BP + JASPAR2022 (5193 models -&gt; 693 clusters) and precomputed MOODS scans; last commit 2023-04. |
| Logomaker | Python | library | pip · conda | 2020 | 0.8.7 2025-03 | ? | 11% | 1.00 | Logos only (also saliency / attribution logos); GitHub activity after the 2025-03 release not checked. |
| MACRO-APE / PERFECTOS-APE | Java | mixed | web | 2013 | 3.0.6 (jar per README) | no | 11% | 0.57 | Companion of the HOCOMOCO collection (v14 site live; HOCOMOCO 2024 NAR 10.1093/nar/gkad1077); exact p-value thresholds for PWMs; last commit 2022-07. |
| ggseqlogo | R | library | CRAN | 2017 | 0.2.2 2025-12 | yes | 7% | 1.00 | Logos only; the logo backend of Signac&#x27;s MotifPlot. |
| pyJASPAR | Python | library | pip · conda | – | 4.0.0 2025-10 | yes | 7% | 0.88 | Offline JASPAR access (through JASPAR 2024/2026 releases); scanning is left to Biopython. |
| seqLogo | R | library | Bioc | – | 1.78.0 2026-04 | yes | 7% | 0.89 | Oldest Bioconductor logo package (since BioC 2.0); activity is the Bioconductor bump only. |


## Co-occupancy / enrichment
{: .sec-navy #cooccupancy }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Index a large track collection (timed); **C2** Query-vs-collection enrichment ranking (timed); **C3** Exact overlap statistics (Fisher) (timed); **C4** Permutation null with genome-aware randomisation (timed); **C5** Custom universe / restricted background (timed); **C6** Covariate-matched background (GC, distance, signal); **C7** Pairwise co-occurrence across many sets (timed); **C8** Metadata search (TF, cell type, antigen); **C9** Scale to 10k+ tracks; **C10** Multiple-testing correction across tracks.

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| **genomeblocks** | best | 0.42 | 6/10 | 0.70 | C1 1.00, C2 0.90, C3 0.90, C5 0.00 |
| GIGGLE | moderate | 0.21 | 5/10 | 0.42 | C1 0.70, C2 0.28, C3 0.28 |
| LOLA | moderate | 0.20 | 6/9 (1 unknown) | 0.30 | C1 0.49, C2 0.00, C3 0.00, C5 0.70 |
| bedtools | moderate | 0.13 | 4/10 | 0.33 | C2 0.00, C3 0.00, C7 1.00 |
| regioneR | low | 0.11 | 2/10 | 0.56 | C5 0.56 |
| valr | untested | – | 1/1 (9 unknown) | 0.99 | C7 0.99 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| valr | 1.00 | 1.00 | 0.99 | – | 1.00 | 1 | 0 |
| regioneR | 0.84 | 0.98 | 0.78 | – | 0.78 | 2 | 0 |
| bedtools | 0.77 | 1.00 | 0.50 | – | 0.90 | 2 | 0 |
| **genomeblocks** | 0.74 | 0.45 | 1.00 | – | 0.90 | 4 | 0 |
| LOLA | 0.71 | 0.99 | 0.40 | – | 0.90 | 3 | 0 |
| GIGGLE | 0.68 | 0.88 | 0.51 | – | 0.70 | 2 | 0 |
| pyranges (0.x) | – | – | 0.44 | – | 1.00 | 1 | 0 |

### `enrichment_search`
{: .no_toc }

Scored at n = 500; fastest 15 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Atlas.search** | genomeblocks | 500 | 15 ms | 0.90 | Spearman 0.9049 vs exact Fisher; vs designed truth 0.8575 | landscape |
| GIGGLE search -s (CLI) | GIGGLE | 500 | 1.6 s | 0.88 | Spearman 0.8773 vs exact Fisher; vs designed truth 0.8705 | landscape |
| LOLA::runLOLA (1 kb-tile universe) | LOLA | 500 | 19 s | 0.97 | Spearman 0.9737 vs exact Fisher; vs designed truth 0.7126 | landscape |
| bedtools fisher ×500 (CLI) | bedtools | 500 | 26 s | 1.00 | Spearman 0.9961 vs exact Fisher; vs designed truth 0.6466 | landscape |

### `index_collection`
{: .no_toc }

Scored at n = 500; fastest 4.1 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Atlas.make (1 kb bins, 4 workers)** | genomeblocks | 500 | 4.1 s | – |  | landscape |
| GIGGLE index (CLI) | GIGGLE | 500 | 33 s | – |  | landscape |
| LOLA::loadRegionDB (500 tracks) | LOLA | 500 | 2.3 min | – |  | landscape |

### `overlap_counts_500_tracks`
{: .no_toc }

Scored at n = 500; fastest 442 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **per-track overlap_any loop** | genomeblocks | 500 | 442 ms | – |  | atlas |
| pyranges overlap(), per track | pyranges (0.x) | 500 | 21 s | – |  | atlas |

### `pairwise_cooccurrence`
{: .no_toc }

Scored at n = 20; fastest 5.3 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| bedtools jaccard x 190 pairs (CLI) | bedtools | 20 | 5.3 s | 1.00 | Pearson 1.00000, max |Δ| 4.77e-08 | landscape |
| valr bed_jaccard x 190 pairs | valr | 20 | 5.6 s | 1.00 | max |Δ| 9.37e-17 | landscape |

### `permutation_null`
{: .no_toc }

Scored at n = 500; fastest 11.5 min.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| regioneR::overlapPermTest (100 perms, per track) | regioneR | 500 | 11.5 min | 0.96 | Spearman 0.9567 over 10 tracks | landscape |

### `universe_enrichment`
{: .no_toc }

Scored at n = 500; fastest 1.5 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Atlas.bootstrap (pool = universe, 100 draws)** | genomeblocks | 500 | 1.5 s | 0.00 | Spearman -0.4141 vs the exact universe Fisher test | landscape |
| LOLA::runLOLA (userUniverse) | LOLA | 500 | 12 s | 1.00 | Spearman 1.0000 over 500 tracks vs the exact universe Fisher test | landscape |
| regioneR::permTest (resampleRegions in the universe, 100 perms, per track) | regioneR | 500 | 32 s | 1.00 | Spearman 1.0000 over 10 tracks vs the exact universe Fisher test | landscape |

### Not run or did not finish
{: .no_toc }

| Engine | task | why |
|---|---|---|
| intervene pairwise --compute jaccard (CLI) | pairwise_cooccurrence | RuntimeError: intervene: ImportError: cannot import name &#x27;Iterable&#x27; from &#x27;collections&#x27; (gb_bench_cli/lib/python3.11/collections/__init__.py) |

### Every tool in this theme
{: .no_toc }

Capability checklist: `index_collection`, `enrichment_search`, `fisher_test`, `permutation_null`, `reference_background`, `metadata_search`, `pairwise_cooccupancy`, `giggle_score`, `scale_10k_tracks`, `local_install`, `web_only`, `parallel`, `in_memory_api`, `bootstrap_resampling_of_query`, `matched_background`, `multiple_testing_correction`, `local_zscore_shift`, `block_bootstrap`, `effect_sizes`, `set_similarity_search`, `overlap_set_visualisation`, `enrichment_plots`, `aggregate_by_metadata`, `public_atlas_bundles`, `per_track_overlap_counts`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 60% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| ChIP-Atlas Enrichment Analysis | Web | web | web | 2018 | – | yes | 56% | 0.62 | Background = random permutation of the query (x1/x10/x100) or a user BED; filters by experiment type and cell-type class; separate Colocalization tool gives antigen x antigen partners (TSV/GML); 2022 NAR update 10.1093/nar/gkac199; Fisher test per the paper (page does not name it). |
| LOLA | R | library | Bioc | 2015 | 1.42.0 | no | 56% | 0.90 | Fisher&#x27;s exact test against a user universe; LOLACore/LOLAExt region DBs with cell-type/antibody metadata; GitHub last commit 2020-09, Bioconductor version bumps only. |
| LOLAweb | Web | web | web | 2018 | – | no | 48% | 0.50 | Containerised LOLA; lolaweb.databio.org refused connections on 2026-10-22; last commit 2022-04. |
| TFEA.ChIP | R | library | Bioc | 2019 | 1.32.0 | no | 48% | 0.75 | Adjacent: input is a gene list (e.g. DE genes), regions enter only via the prebuilt TF-gene DB; last commit 2025-07. |
| UniBind Enrichment | mixed | mixed | web | 2021 | – | ? | 48% | 0.67 | unibind.uio.no/enrichment up on 2026-10-22; CLI version on Bitbucket (activity not checked); LOLA DBs on Zenodo (2020-11). |
| regioneReloaded | R | library | Bioc | 2023 | 1.14.0 | no | 44% | 0.78 | crosswisePermTest: z-score association matrix across many region sets with clustering/plots; last commit 2023-06. |
| Cistrome DB Toolkit | Web | web | web | 2018 | – | ? | 40% | 0.33 | dbtoolkit.cistrome.org reachable 2026-10-22; peak-set mode capped at 50,000 peaks; GIGGLE scoring per the 2019 NAR paper (page itself does not state method). |
| GIGGLE | C/C++ | mixed | conda | 2018 | – | yes | 40% | 0.70 | Search engine over thousands of interval files with odds-ratio/Fisher/combo-score enrichment; Travis CI config only; no tagged releases; Origin of the Fisher x odds-ratio &#x27;combo score&#x27;; metadata index/filter (-m/-u); no GitHub releases or bioconda recipe found (conda env file only); python/go bindings are community-maintained and dated; last commit 2026-03. |
| OLOGRAM (pygtftk) | Python | mixed | pip · conda | 2019 | 1.6.2 2022-03 | yes | 40% | 1.00 | gtftk ologram with -k threads and --bed-incl/--bed-excl; OLOGRAM-MODL mines multi-set combinations; last commit 2026-01, last PyPI release 2022-03. |
| ReMapEnrich | R | library | src | 2017 | – | no | 40% | 0.56 | Introduced with ReMap 2018 (issued online 2017); shuffles query against catalogue with optional universe; README mentions &#x27;single core&#x27; default; last commit 2021-11; install via devtools only. |
| GAT | Python | cli | pip · conda | 2013 | 1.3.6 2019-06 | no | 36% | 0.70 | --num-threads, --workspace-file, --isochore-file; fold/l2fold/empirical p/q per annotation; unmaintained since 2019. |
| regioneR | R | library | Bioc | 2015 | 1.44.0 | no | 36% | 0.78 | Circular/random/resample-from-universe permutation tests for one pair of region sets; mc.cores parallelism; last GitHub commit 2025-05. |
| RELI | C/C++ | cli | src | 2018 | v0.90 2019-01 | yes | 36% | 0.44 | SNP/GWAS-centric: null draws from background variant sets with LD expansion; last commit 2026-06 (README). |
| Coloc-stats | Web | web | web | 2018 | – | ? | 32% | 0.33 | Page returned only a title on 2026-10-22 (JS app); method list from the 2018 paper. |
| geniml / BEDbase | Python | mixed | pip | – | 0.8.4 2026-01 | yes | 32% | 1.00 | Ranking is embedding similarity, not a statistical enrichment test; bedshift module provides perturbation nulls; geniml/BEDbase journal citation not found via Crossref (docs point to a how-to-cite page). |
| Genomic HyperBrowser | Python | web | web | 2010 | – | no | 28% | 0.43 | GSuite HyperBrowser (GigaScience 2017, 10.1093/gigascience/gix032) adds collection-level ranking; hyperbrowser.uio.no returned 503 on 2026-10-22; last commit 2023-09. |
| pyranges | Python | library | pip · conda | 2019 | 0.1.6 2026-10 | yes | 28% | 1.00 | Engine with per-pair statistics only; pyranges v0 line (0.1.x) still released 2026-10 alongside pyranges 1.x. |
| IGD | C/C++ | cli | src | 2021 | – | no | 24% | 0.22 | Precursor of gtars-igd; single/batch queries and Seqpare similarity; R wrapper IGDr; last push 2021; GIGGLE-like index with overlap counts/Seqpare similarity but no p-values; R and Python bindings; last commit 2021-09. |
| Intervene | Python | mixed | pip · conda | 2017 | 0.6.5 2020-11 | no | 24% | 0.75 | Track x track matrix tool (also a Shiny app); last commit 2023-10. |
| nullranges | R | library | Bioc | 2023 | 1.18.0 | yes | 24% | 1.00 | Block bootstrap of the reference/feature ranges with segmentation and exclusion regions; the test statistic is left to the user (plyranges); last commit 2026-08. |
| rGREAT | R | mixed | Bioc | 2023 | 2.14.0 2026-04 | yes | 24% | 1.00 | Local GREAT (binomial over basal+extension domains) for any organism/gene-set plus the web-service client; plotRegionGeneAssociations gives distance-to-TSS plots; Shiny viewer; Adjacent (annotation theme): functional enrichment of one region set; last commit 2026-01. |
| TF-COMB | Python | library | pip | 2022 | 1.1.1 2025-06 | no | 24% | 0.88 | Pairwise co-occupancy of TFBS/peaks (not query-vs-collection); last commit/release 2025-06. |
| bedtools | C/C++ | cli | conda | 2010 | v2.31.1 2023-11 | yes | 20% | 0.90 | Reference CLI for genome arithmetic; set ops limited to subtract/complement (union = cat|sort|merge); &#x27;closest&#x27; against a TSS BED is the annotation use; last tagged release 2023 but repo pushed 2026-06; genomecov -bg -scale is the classic BAM-&gt;bedGraph route; GitHub push 2026-06 but no release since 2023-11; pybedtools wraps it; Primitives only: fisher (2x2 on one pair), jaccard, reldist, shuffle (-incl/-excl backgrounds) and pairtobed -type either/both for loop anchors; collection ranking must be scripted. |
| GenometriCorr | R | library | src | 2012 | 1.1.24 2021-03 | no | 20% | 0.50 | Version string from last commit message (2021-03), no tags; not on CRAN/Bioconductor. |
| BART2 | Python | cli | pip | 2018 | 2.0 2023-12 | yes | 16% | 0.71 | Requires multi-GB prebuilt library download; the bioconda package &#x27;bart&#x27; is an unrelated tool; last commit 2026-03. |
| Bedshift | Python | mixed | pip | 2021 | 1.1.1 2021-04 | yes | 16% | 1.00 | Null-set generator, not a test; repository archived 2026-03 and functionality folded into geniml. |
| poverlap | Python | cli | src | – | – | no | 16% | 0.33 | No paper; the PyPI package named &#x27;poverlap&#x27; is an unrelated 2D-geometry package; last commit 2020-07. |
| ColoWeb | Web | web | web | 2015 | – | no | 12% | 0.14 | Host projects.insilico.us.com no longer resolves (2026-10-22); historical entry. |
| RegRegSEA | Web | web | web | – | – | ? | 12% | 0.33 | Saarland CCB web server; citation on site is a 2026 placeholder (unpublished); GitHub link is an issue tracker only. |
| GenomicDistributions | R | library | Bioc | 2022 | 1.20.0 2026-04 | no | 8% | 0.90 | QC-style summaries: calcFeatureDist (TSS distance), calcPartitions (promoterCore/Prox, exon, intron, UTRs, intergenic), chromosome bins, GC; plots for each; repo last pushed 2025-04; Adjacent: QC/characterisation of region sets, no enrichment test; last push 2025-04. |
| IntervalStats | C/C++ | cli | src | 2012 | – | ? | 8% | 0.17 | Analytic (non-permutation) proximity test; distributed from a Princeton server, no public repo found; wrapped by Coloc-stats. |
| bedtk | C/C++ | cli | conda | – | v1.2 2025-08 | no | 4% | 0.30 | Minimal, very fast subset of bedtools (isec, flt, cov, sort, merge); README warns the CLI may change; no sorting required; Engine only (fast overlap), no statistics; v1.2 tagged 2025-08-24, just before the activity window. |
| Cistrome-GO | Web | web | web | 2019 | – | ? | 4% | 0.33 | Adjacent (annotation theme): functional enrichment of one peak set, not co-occupancy against a track collection. |


## 3D genome & contact graphs
{: .sec-navy #graph3d }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Parse + deduplicate pairs from alignments; **C2** Pairs -&gt; binned multi-resolution matrix (.cool / .mcool / .hic) (timed); **C3** Matrix balancing (ICE / KR) (timed); **C4** Expected / observed-over-expected (distance decay) (timed); **C5** Read .cool / .mcool / .hic (timed); **C6** Loop / interaction calling (timed); **C7** TAD / insulation boundaries (timed); **C8** A/B compartments (eigenvectors, saddle) (timed); **C9** Differential contacts / loops; **C10** Loop-feature overlap and anchor-gene assignment (timed).

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| HiCExplorer | best | 0.43 | 10/10 | 0.43 | C3 0.96, C6 0.01, C7 0.64, C8 0.11 |
| Juicer / Juicer Tools | best | 0.42 | 9/10 | 0.47 | C6 0.47 |
| cooltools | good | 0.34 | 5/10 | 0.67 | C6 0.46, C7 0.58, C8 0.97 |
| hictk | good | 0.29 | 3/10 | 0.96 | C3 0.96 |
| cooler | good | 0.27 | 3/10 | 0.91 | C2 0.85, C3 0.98 |
| MAPS | moderate | 0.20 | 3/10 | 0.66 | C6 0.66 |
| Mustache | moderate | 0.20 | 4/10 | 0.49 | C6 0.49 |
| **genomeblocks** | moderate | 0.19 | 2/10 | 0.95 | C6 0.91, C10 1.00 |
| FitHiChIP | low | 0.12 | 4/10 | 0.29 | C6 0.29 |
| GenomicInteractions | low | 0.09 | 1/10 | 0.86 | C10 0.86 |
| bedtools | low | 0.05 | 1/10 | 0.54 | C10 0.54 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| igraph | 1.00 | 1.00 | 1.00 | – | 1.00 | 3 | 0 |
| hictk | 0.99 | 0.98 | 0.98 | – | 1.00 | 1 | 0 |
| cooler | 0.97 | 0.98 | 0.93 | – | 1.00 | 2 | 0 |
| **genomeblocks** | 0.95 | 0.95 | 1.00 | – | 0.90 | 17 | 0 |
| GenomicInteractions | 0.84 | 1.00 | 0.80 | – | 0.75 | 3 | 0 |
| cooltools | 0.83 | 0.73 | 0.86 | – | 0.90 | 3 | 0 |
| bedtools | 0.79 | 1.00 | 0.54 | – | 0.90 | 2 | 0 |
| HiCExplorer | 0.67 | 0.45 | 0.82 | – | 0.80 | 4 | 0 |
| Mustache | 0.62 | 0.74 | 0.66 | – | 0.50 | 2 | 0 |
| Juicer / Juicer Tools | 0.61 | 0.82 | 0.57 | – | 0.50 | 2 | 0 |
| MAPS | 0.57 | 1.00 | 0.66 | – | 0.29 | 1 | 0 |
| FitHiChIP | 0.53 | 1.00 | 0.29 | – | 0.50 | 1 | 0 |

### `anchor_gene_annotation`
{: .no_toc }

Scored at n = 73k; fastest 35 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Architecture.annotate** | genomeblocks | 73k | 35 ms | – |  | architecture |

### `bin_pairs_to_matrix`
{: .no_toc }

Scored at n = 5M; fastest 3.3 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks count_pairs_2d (500 kb)** | genomeblocks | 5M | 3.3 s | – |  | pairs |
| cooler cload pairs (500 kb, CLI) | cooler | 5M | 9.2 s | – |  | pairs |

### `build_graph`
{: .no_toc }

Scored at n = 50k; fastest 187 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| igraph (R) | igraph | 50k | 187 ms | – |  | landscape |

### `compartments`
{: .no_toc }

Scored at n = 8k; fastest 16 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| cooltools eigs-cis (CLI) | cooltools | 8k | 16 s | 0.97 | |Pearson| 0.9738 with the simulated A/B compartments, mean over 2 chromosomes (0.975, 0.973) | landscape |
| HiCExplorer hicPCA (CLI) | HiCExplorer | 8k | 103 s | 0.16 | |Pearson| 0.1560 with the simulated A/B compartments, mean over 2 chromosomes (0.140, 0.172) | landscape |

### `components`
{: .no_toc }

Scored at n = 114k; fastest 14 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **graph-tool (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 14 ms | 1.00 |  | backends |
| **genomeblocks (graph-tool)** | genomeblocks | 114k | 16 ms | 1.00 | reference: 6,638 components | landscape |
| **scipy (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 18 ms | 1.00 |  | backends |
| **igraph (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 59 ms | 1.00 |  | backends |
| **networkx (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 447 ms | 1.00 |  | backends |
| igraph (R) | igraph | 50k | 19 ms | 1.00 | 6,638 components of 2+ vertices | landscape |

### `contact_graph`
{: .no_toc }

Scored at n = 114k; fastest 155 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Architecture.make** | genomeblocks | 114k | 155 ms | – |  | architecture |

### `hub_mining`
{: .no_toc }

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Architecture.prime_hubs** | genomeblocks | – | 8.2 ms | – |  | architecture |

### `loop_calling`
{: .no_toc }

Scored at n = 8k; fastest 3.6 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks mustache() (port of Mustache 1.3.3)** | genomeblocks | 8k | 3.6 s | 0.48 | F1 0.478 (precision 0.696, recall 0.364; 23 called, 44 planted, within 2 bins) | landscape |
| **genomeblocks hiccups() (port of Juicer HiCCUPS, cooler weights)** | genomeblocks | 8k | 6.2 s | 0.62 | F1 0.615 (precision 0.952, recall 0.455; 21 called, 44 planted, within 2 bins) | landscape |
| cooltools expected-cis + dots (CLI) | cooltools | 8k | 25 s | 0.65 | F1 0.647 (precision 0.917, recall 0.500; 24 called, 44 planted, within 2 bins) | landscape |
| Mustache 1.3.3 (CLI) | Mustache | 8k | 44 s | 0.48 | F1 0.478 (precision 0.696, recall 0.364; 23 called, 44 planted, within 2 bins) | landscape |
| HiCExplorer hicDetectLoops (CLI) | HiCExplorer | 8k | 60 s | 0.02 | F1 0.017 (precision 0.013, recall 0.023; 76 called, 44 planted, within 2 bins) | landscape |
| Juicer juicer_tools 2.13.06 hiccups --cpu (CLI) | Juicer / Juicer Tools | 8k | 90 s | 0.64 | F1 0.636 (precision 0.955, recall 0.477; 22 called, 44 planted, within 2 bins) | landscape |

### `loop_overlap_both`
{: .no_toc }

Scored at n = 50k; fastest 24 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Pairs.overlapping** | genomeblocks | 50k | 24 ms | 1.00 | identical (19,717) | landscape |
| InteractionSet::overlapsAny (both anchors) | GenomicInteractions | 50k | 86 ms | 1.00 | identical (19,717) | landscape |
| bedtools pairtobed (CLI) | bedtools | 50k | 539 ms | 1.00 | reference | landscape |

### `loop_overlap_either`
{: .no_toc }

Scored at n = 50k; fastest 23 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Pairs.overlapping** | genomeblocks | 50k | 23 ms | 1.00 | identical (42,500) | landscape |
| InteractionSet::overlapsAny (either anchor) | GenomicInteractions | 50k | 45 ms | 1.00 | identical (42,500) | landscape |
| bedtools pairtobed (CLI) | bedtools | 50k | 565 ms | 1.00 | reference | landscape |

### `loops_fithichip`
{: .no_toc }

Scored at n = 53k; fastest 2.5 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks fithichip()** | genomeblocks | 53k | 2.5 s | 1.00 | 52,797 of the original&#x27;s 52,797 loops, same input (FitHiChIP 11.0 peak-to-all) | loops |
| FitHiChIP | FitHiChIP | 53k | 5.4 min | 1.00 | reference: the original&#x27;s 52,797 loops (FitHiChIP 11.0 peak-to-all) | loops |

### `loops_hiccups`
{: .no_toc }

Scored at n = 182; fastest 10 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks hiccups()** | genomeblocks | 182 | 10 s | 1.00 | 182 of the original&#x27;s 182 loops, same input (Juicer HiCCUPS, 8 threads each, the port on Juicer&#x27;s SCALE vectors) | hiccups |
| juicer_tools 2.13.06 hiccups --cpu | Juicer / Juicer Tools | 182 | 2.7 min | 1.00 | reference: the original&#x27;s 182 loops (Juicer HiCCUPS, 8 threads each, the port on Juicer&#x27;s SCALE vectors) | hiccups |

### `loops_maps`
{: .no_toc }

Scored at n = 24k; fastest 4.7 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks maps()** | genomeblocks | 24k | 4.7 s | 1.00 | 23,941 of the original&#x27;s 23,941 loops, same input (MAPS (MAPS.py + regression; feather&#x27;s file split not counted)) | maps |
| MAPS | MAPS | 24k | 49 s | 1.00 | reference: the original&#x27;s 23,941 loops (MAPS (MAPS.py + regression; feather&#x27;s file split not counted)) | maps |

### `loops_mustache`
{: .no_toc }

Scored at n = 660; fastest 8.0 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks mustache()** | genomeblocks | 660 | 8.0 s | 1.00 | 660 of the original&#x27;s 660 loops, same input (Mustache 1.3.3, 4 processes each) | mustache |
| Mustache 1.3.3 | Mustache | 660 | 73 s | 1.00 | reference: the original&#x27;s 660 loops (Mustache 1.3.3, 4 processes each) | mustache |

### `matrix_balancing`
{: .no_toc }

Scored at n = 8k; fastest 61 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| cooler balance (CLI) | cooler | 8k | 61 s | 0.98 | Pearson 0.9774 between log weight and minus log simulated bias (7,812 of 8,000 bins) | landscape |
| hictk balance ice (CLI) | hictk | 8k | 70 s | 0.98 | Pearson 0.9774 between log weight and minus log simulated bias (7,812 of 8,000 bins) | landscape |
| HiCExplorer hicCorrectMatrix correct --correctionMethod ICE (CLI) | HiCExplorer | 8k | 72 s | 0.98 | Pearson 0.9795 between log weight and minus log simulated bias (7,817 of 8,000 bins) | landscape |

### `oe_normalisation`
{: .no_toc }

Scored at n = 114k; fastest 9.6 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Architecture.normalize** | genomeblocks | 114k | 9.6 ms | – |  | architecture |

### `pagerank`
{: .no_toc }

Scored at n = 114k; fastest 106 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **scipy (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 106 ms | 1.00 |  | backends |
| **genomeblocks (graph-tool)** | genomeblocks | 114k | 108 ms | 1.00 | reference | landscape |
| **graph-tool (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 112 ms | 1.00 |  | backends |
| **igraph (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 148 ms | 1.00 |  | backends |
| **networkx (through genomeblocks&#x27; backend)** | genomeblocks | 114k | 686 ms | 1.00 |  | backends |
| igraph (R) | igraph | 50k | 160 ms | 1.00 | max |Δ| 3.24e-14 | landscape |

### `pairs_streaming`
{: .no_toc }

Scored at n = 5M; fastest 3.2 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks count_pairs (50 kb x partner chrom)** | genomeblocks | 5M | 3.2 s | – |  | pairs |

### `read_bedpe`
{: .no_toc }

Scored at n = 50k; fastest 26 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Pairs.make** | genomeblocks | 50k | 26 ms | – |  | landscape |
| GenomicInteractions::makeGenomicInteractionsFromFile | GenomicInteractions | 50k | 233 ms | – |  | landscape |

### `read_cool_mcool`
{: .no_toc }

Scored at n = 114k; fastest 4.1 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Architecture.add_mcool** | genomeblocks | 114k | 4.1 s | – |  | architecture |

### `tad_boundaries`
{: .no_toc }

Scored at n = 8k; fastest 3.7 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| cooltools insulation 100 kb (CLI) | cooltools | 8k | 3.7 s | 0.58 | F1 0.579 (precision 0.410, recall 0.988; 205 called, 85 planted, within 2 bins) | landscape |
| HiCExplorer hicFindTADs (CLI) | HiCExplorer | 8k | 4.5 s | 0.65 | F1 0.653 (precision 0.497, recall 0.953; 163 called, 85 planted, within 2 bins) | landscape |

### Every tool in this theme
{: .no_toc }

Capability checklist: `read_bedpe`, `loop_region_overlap`, `pairs_streaming`, `bin_pairs_to_matrix`, `read_cool_mcool`, `read_hic_juicer`, `oe_normalisation`, `loop_calling`, `contact_graph`, `graph_algorithms`, `hub_mining`, `anchor_gene_annotation`, `hichip_support`, `graph_export`, `arc_visualisation`, `parallel`, `in_memory_api`, `pairs_parse_dedup`, `matrix_balancing`, `tad_calling`, `compartment_calling`, `differential_contacts`, `aggregate_pileups`, `virtual_4c`, `contact_qc`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 64% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| cLoops2 | Python | cli | pip | 2021 | 0.0.5 2023-12 | no | 76% | 0.43 | Closest end-to-end analogue of a loop-&gt;graph-&gt;hub workflow (networkx, HITS, .sif export, closest-TSS anchor classes); cLoops (2020, 10.1093/bioinformatics/btz651) is the predecessor; last commit 2025-06. |
| HiCExplorer | Python | cli | pip · conda | 2018 | 3.7.6 2024-11 | no | 68% | 0.80 | Bins from BAM (not .pairs); .hic via hicConvertFormat; Galaxy HiCExplorer 3 server (NAR 2020, 10.1093/nar/gkaa220); PyPI stale at 2.2.1 (2019), bioconda current; last commit 2024-11. |
| FAN-C | Python | mixed | pip · conda | 2020 | 0.9.28 2024-03 | no | 64% | 0.88 | Last commit/release 2024-03. |
| HiCExperiment / HiContacts | R | library | Bioc | 2024 | HiCExperiment 1.12.0 / HiContacts 1.14.0 | yes | 56% | 1.00 | Binning via HiCool (hicstuff wrapper); loops called externally (chromosight) and imported; arcs via plotgardener; last commits 2026-10. |
| Juicer / Juicer Tools | Java | cli | src | 2016 | 1.6 2020-06 | no | 56% | 0.50 | HiCCUPS (loops, Rao 2014 Cell 10.1016/j.cell.2014.11.021, GPU) and Arrowhead (TADs) ship in juicer_tools; dump oe; last commit 2025-08; no bioconda package found. |
| cooltools | Python | mixed | pip · conda | 2024 | 0.7.1 2024-07 | no | 44% | 0.90 | expected-cis + obs/exp, dots (HiCCUPS-style) loop caller; last commit 2025-07, just before the window. |
| GENOVA | R | library | src | 2021 | v1.0.0-TheCartographer 2021-11 | no | 40% | 0.75 | Loops consumed as BEDPE for APA; no caller; last commit 2024-03; GitHub install only; Hi-C analysis + pyramid/hic.matrixplot figures with loop/TAD overlays (APA, saddle, insulation); last commit 2024-03. |
| FitHiChIP | mixed | cli | src | 2019 | 11.0 2022-12 | no | 36% | 0.50 | PeakInferHiChIP.sh builds 1D coverage from HiC-Pro valid pairs (and DE/SC/RE reads when present) and calls peaks with MACS2; last commit 2024-04; v11 (2022-12) added .hic/.cool input; last commit 2024-04; no bioconda package found. |
| GenomicInteractions | R | library | conda · Bioc | 2015 | 1.46.0 | no | 32% | 0.75 | Closest R analogue of a loop-&gt;graph workflow (export.igraph) but graph analysis is delegated to igraph; last GitHub commit 2022-04; Data structures underpinning diffHic/HiCExperiment; last GitHub commit 2025-07. |
| HiC-DC+ | R | library | Bioc | 2021 | 1.20.0 | no | 32% | 0.75 | Last GitHub commit 2022-06. |
| hicstuff | Python | mixed | pip · conda | 2020 | v3.2.5 2026-05 | yes | 32% | 0.88 | Zenodo-only citation; backs Bioconductor HiCool; last commit 2026-09. |
| chromosight | Python | cli | pip · conda | 2020 | 1.6.3 2022-10 | yes | 28% | 0.88 | quantify scores user BED2D/BEDPE loops; last commit 2026-09, last release 2022-10. |
| HiC-Pro | mixed | cli | src | 2015 | v3.1.0 2021-08 | no | 28% | 0.43 | Produces the allValidPairs that hichipper / HiChIP-Peaks / FitHiChIP consume; no 1D track itself; last commit 2022-06; allValidPairs is the de-facto HiChIP input for FitHiChIP/hichipper; last commit 2022-06; no bioconda package found under hic-pro/hicpro. |
| hictk | C/C++ | mixed | conda | 2024 | v2.2.0 2025-11 | yes | 28% | 1.00 | Format layer (convert/dump/load/merge/zoomify/balance), no expected/o-e or loop calling; last commit 2026-10. |
| QuIN | Web | web | web | 2016 | – | no | 28% | 0.50 | Earliest dedicated contact-network hub tool; quin.jax.org refused connections 2026-10-22; last commit 2020-12. |
| cooler | Python | mixed | pip · conda | 2019 | 0.10.4 2025-07 | yes | 24% | 1.00 | Reference .cool format; ICE balancing (no o/e, that is cooltools); last commit 2026-01. |
| coolpup.py | Python | mixed | pip · conda | 2020 | 1.1.0 2023-07 | no | 24% | 0.88 | Last commit 2024-10. |
| diffHic | R | library | Bioc | 2015 | 1.44.0 | yes | 24% | 0.88 | Only commit in window is the automated Bioconductor version bump (2026-04). |
| Mustache | Python | cli | pip · conda | 2020 | 1.3.3 2024-05 | no | 24% | 0.50 | Last commit 2025-01. |
| Activity-by-Contact (ABC) | Python | cli | src | 2019 | ABCv1.1.2 2024-03 | yes | 20% | 0.57 | hub_mining marked as adjacent (E-G scoring, not graph hubs); last commit 2026-08. |
| CHiCAGO | R | mixed | Bioc | 2016 | 1.40.0 | ? | 20% | 1.00 | Capture Hi-C specific; Bitbucket activity not checked. |
| Cicero | R | library | Bioc | 2018 | 1.30.0 | yes | 20% | 1.00 | Not 3C data: edges are co-accessibility, communities via igraph; last commit 2026-03. |
| graph-tool | C/C++ | library | conda | 2014 | release-3.9 2026-09 | yes | 20% | 0.88 | Generic graph engine; figshare DOI not Crossref-verifiable (DataCite); not on PyPI. |
| HiCHub | Python | cli | pip | 2022 | 1.2.0 2022-09 | no | 20% | 0.38 | Differential hubs between two conditions; no journal version found in Crossref; last commit 2025-04. |
| pairtools | Python | cli | pip · conda | 2024 | 1.1.2 2024-12 | yes | 20% | 0.90 | Pairs-level only; binning delegated to cooler; last commit 2026-07. |
| dcHiC | mixed | cli | src | 2022 | – | no | 16% | 0.29 | Adjacent (compartments); last commit 2024-01. |
| Hi-Cociety | R | mixed | CRAN | 2025 | – | yes | 16% | 0.60 | CRAN package archived 2026-03-21 (uncorrected check issues); GitHub last commit 2026-02. |
| Hub-enhancer / enhancer-community analyses (method only) | mixed | mixed | src | 2018 | – | ? | 16% | 0.33 | No packaged tool found: Huang 2018 (Nat Commun) and Madsen 2020 (Nat Genet, 10.1038/s41588-020-0709-z, enhancer-capture Hi-C) are analyses; Hnisz 2013 (Cell, 10.1016/j.cell.2013.09.053) super-enhancers were called with ROSE, not a graph method. |
| igraph | mixed | library | pip · conda · CRAN | 2006 | 1.0.0 2025-10 | yes | 16% | 1.00 | Generic graph engine (used by GenomicInteractions, Cicero, HiCHub); 2006 paper has no DOI, a 2026 PLOS One paper (10.1371/journal.pone.0355108) is the modern citation. |
| MAPS | mixed | cli | src | 2019 | – | no | 16% | 0.29 | Python+R pipeline (feather/bwa), last commit 2023-06. |
| MINE | Python | cli | src | 2023 | – | yes | 16% | 0.44 | PyTorch/CUDA research code (Python 3.6 env); last commit 2025-11. |
| multiHiCcompare | R | library | Bioc | 2019 | 1.30.0 | yes | 16% | 0.88 | Differential-contact tool; cool import via HiCcompare helpers; last commit 2026-07. |
| networkx | Python | library | pip · conda | 2008 | 3.7 2026-09 | yes | 16% | 1.00 | Generic graph engine (used by cLoops2). |
| Peakachu | Python | cli | pip | 2020 | 2.3 2024-03 | no | 16% | 0.50 | The bioconda package &#x27;peakachu&#x27; is an unrelated CLIP-seq tool; last commit 2024-03. |
| plotgardener | R | library | Bioc | 2022 | 1.18.0 | yes | 16% | 1.00 | Viewer used by the OHCA/HiContacts stack for arcs; GitHub activity not checked; Sushi&#x27;s successor: absolute-coordinate multi-panel pages (plotHicTriangle/Square, plotPairsArches, plotSignal, plotGenes from TxDb/OrgDb) for publication figures. |
| HiCcompare | R | library | Bioc | 2018 | 1.34.0 | yes | 12% | 0.88 | Last commit 2026-09. |
| hichipper | Python | cli | pip | 2018 | 0.7.7 2019-06 | no | 12% | 0.62 | Calls peaks directly from HiChIP reads with a restriction-site background model; last commit 2022-07; Unmaintained since 2022-07. |
| hictkpy | Python | library | pip · conda | 2024 | v1.4.0 2025-12 | yes | 12% | 1.00 | Writer takes already-binned pixels (binning of raw pairs is hictk CLI&#x27;s job). |
| Juicebox | Java | gui | web | 2016 | v2.20.00 2022-10 | yes | 12% | 0.62 | Viewer; last commit 2026-09; Hi-C-centric viewer (.hic) with 1D bigWig/gene tracks and 2D loop/domain annotations, image export, shareable sessions; desktop commits to 2026-09, juicebox.js web. |
| pgltools | Python | mixed | src | 2017 | – | no | 12% | 0.43 | bedtools-like for loops; not on PyPI/bioconda (404 on both); last commit 2022-11. |
| pyGenomeTracks | Python | cli | pip · conda | 2020 | 3.9 2024-07 | yes | 12% | 0.90 | Viewer; last commit 2026-08; De-facto standard static track figure (bigwig, bedgraph, GTF/bed12 genes, links arcs, cool/h5 Hi-C triangles, TADs); no BAM track; last PyPI release 2024-07 but commits to 2026-08. |
| straw | mixed | library | pip · CRAN | – | 1.3.1 2022-05 | no | 12% | 0.50 | No standalone paper (described with Juicebox); C++/Python/R/Java/JS bindings; last commit 2025-04. |
| bedtools | C/C++ | cli | conda | 2010 | v2.31.1 2023-11 | yes | 8% | 0.90 | Reference CLI for genome arithmetic; set ops limited to subtract/complement (union = cat|sort|merge); &#x27;closest&#x27; against a TSS BED is the annotation use; last tagged release 2023 but repo pushed 2026-06; genomecov -bg -scale is the classic BAM-&gt;bedGraph route; GitHub push 2026-06 but no release since 2023-11; pybedtools wraps it; Primitives only: fisher (2x2 on one pair), jaccard, reldist, shuffle (-incl/-excl backgrounds) and pairtobed -type either/both for loop anchors; collection ranking must be scripted. |
| HiC1Dmetrics | Python | mixed | pip | 2021 | 0.2.12 2025-08 | no | 8% | 0.75 | Adjacent (1D metrics); last commit/release 2025-08. |
| HiCseg | R | library | src | 2014 | – | no | 8% | 0.43 | Adjacent (TADs); removed from CRAN and archived 2024-03-24. |
| HiGlass | mixed | web | pip | 2018 | v2.3.5 2026-06 | yes | 8% | 0.90 | Viewer; higlass-python 1.4.0 (2026-03) for notebooks; Multiscale interactive Hi-C + track browser with linked views; Jupyter via higlass-python; needs a tile server (not a single-file export); SVG export of views. |
| iced | Python | library | pip | 2015 | 0.6.0 2025-05 | no | 8% | 0.86 | Normalisation function only (balancing, not o/e); KR equivalents live inside Juicer/cooler/hictk. |
| SciPy (csgraph) | Python | library | pip · conda | 2020 | 1.17 2026 | yes | 8% | 1.00 | connected_components and sparse linear algebra on an adjacency matrix: genomeblocks&#x27; pip-default graph backend. |
| 3D-GNOME | Web | web | web | 2016 | – | ? | 4% | 0.50 | Adjacent (3D modelling); 2.0 in NAR 2020 (10.1093/nar/gkaa388); 3dgnome.cent.uw.edu.pl refused connections 2026-10-22. |
| C3D | mixed | cli | src | 2018 | – | no | 4% | 0.29 | Predicts interactions from accessibility correlation (bash + R); last commit 2018-08. |
| ROSE | Python | cli | src | 2013 | 1.3.2 2024-05 | yes | 4% | 0.50 | Stitching (12.5 kb) + TSS exclusion + input-subtracted rpm ranking + tangent-line knee; GitHub push 2025-11; Docker image provided; original Young-lab code was Python 2; Adjacent: SE calling only; stjude mirror of the Young-lab Bitbucket code (README commit 2025-11), Python 2 era. |


## Expression / DE
{: .sec-navy #expression }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Count reads per gene from BAM (timed); **C2** Import Salmon / kallisto / RSEM (timed); **C3** Transcript -&gt; gene with length offsets (timed); **C4** Library-size normalisation (size factors, TMM, upper quartile) (timed); **C5** Low-count filtering (timed); **C6** Dispersion estimation (NB, empirical Bayes) (timed); **C7** DE test with a design matrix (Wald, LRT, QL F, voom) (timed); **C8** LFC shrinkage (timed); **C9** Variance-stabilising transforms (VST, rlog, voom) (timed); **C10** Batch correction / covariates (timed).

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|
| InMoose | best | 0.46 | 6/9 (1 unknown) | 0.69 | C6 0.90, C7 0.82, C9 0.73, C10 0.33 |
| DESeq2 | good | 0.40 | 5/10 | 0.80 | C4 0.64, C6 0.79, C7 0.67, C8 0.93, C9 1.00 |
| **genomeblocks** | good | 0.40 | 4/10 | 1.00 | C2 1.00, C3 1.00, C4 0.98, C5 1.00 |
| edgeR | good | 0.40 | 5/10 | 0.79 | C4 0.67, C5 0.90, C6 0.84, C7 0.76 |
| limma | good | 0.39 | 4/10 | 0.98 | C7 0.98, C10 0.98 |
| PyDESeq2 | good | 0.31 | 5/10 | 0.62 | C4 1.00, C6 0.59, C7 0.46, C8 0.66, C9 0.38 |
| tximport | moderate | 0.18 | 2/10 | 0.88 | C2 0.90, C3 0.86 |
| pytximport | moderate | 0.14 | 2/10 | 0.72 | C3 0.72 |
| Rsubread featureCounts | low | 0.10 | 1/10 | 1.00 | C1 1.00 |
| RNAnorm | low | 0.08 | 1/10 | 0.79 | C4 0.79 |
| sva / ComBat-seq | low | 0.05 | 1/10 | 0.49 | C10 0.49 |
| HTSeq | low | 0.04 | 1/10 | 0.44 | C1 0.44 |
| GenomicAlignments | untested | – | 1/1 (9 unknown) | 0.88 | C1 0.88 |

### Scores
{: .no_toc }

| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |
|---|---|---|---|---|---|---|---|
| Rsubread featureCounts | 1.00 | 1.00 | 1.00 | – | 1.00 | 1 | 0 |
| limma | 0.99 | 0.98 | 1.00 | – | 1.00 | 2 | 0 |
| **genomeblocks** | 0.96 | 1.00 | 0.95 | 1.00 | 0.90 | 8 | 2 |
| GenomicAlignments | 0.96 | 1.00 | 0.89 | – | 1.00 | 1 | 0 |
| tximport | 0.96 | 1.00 | 0.88 | – | 1.00 | 2 | 0 |
| edgeR | 0.92 | 0.97 | 0.82 | – | 1.00 | 5 | 0 |
| DESeq2 | 0.91 | 0.95 | 0.80 | – | 1.00 | 5 | 0 |
| pytximport | 0.89 | 1.00 | 0.72 | – | 1.00 | 1 | 0 |
| PyDESeq2 | 0.89 | 0.96 | 0.64 | 1.00 | 1.00 | 5 | 2 |
| InMoose | 0.88 | 0.93 | 0.74 | – | 1.00 | 4 | 0 |
| RNAnorm | 0.79 | 1.00 | 0.49 | 1.00 | 0.80 | 2 | 2 |
| sva / ComBat-seq | 0.79 | 0.80 | 0.62 | – | 1.00 | 1 | 0 |
| HTSeq | 0.76 | 1.00 | 0.44 | – | 1.00 | 1 | 0 |
| data.table | – | – | 1.00 | – | 0.90 | 1 | 0 |

### `batch_correction`
{: .no_toc }

Scored at n = 60k; fastest 3.5 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| limma removeBatchEffect | limma | 60k | 3.5 s | 1.00 | batch removed 1.000 (residual batch effect vs the simulated one: r -0.051, 0.929 before) x group effect kept 1.000 (r with the simulated log2 FC 0.961, 0.920 before) | landscape |
| sva ComBat_seq | sva / ComBat-seq | 60k | 48 s | 0.74 | batch removed 0.736 (residual batch effect vs the simulated one: r 0.245, 0.929 before) x group effect kept 1.000 (r with the simulated log2 FC 0.920, 0.920 before) | landscape |
| InMoose 0.9.1 pycombat_seq | InMoose | 60k | 3.5 min | 0.74 | batch removed 0.736 (residual batch effect vs the simulated one: r 0.245, 0.929 before) x group effect kept 1.000 (r with the simulated log2 FC 0.920, 0.920 before) | landscape |
| limma removeBatchEffect | limma | 20k | 303 ms | 0.96 | batch removed 0.964 (residual batch effect vs the simulated one: r 0.032, 0.879 before) x group effect kept 0.995 (r with the simulated log2 FC 0.867, 0.871 before) | landscape |
| sva ComBat_seq | sva / ComBat-seq | 20k | 5.4 s | 0.86 | batch removed 0.861 (residual batch effect vs the simulated one: r 0.122, 0.879 before) x group effect kept 0.995 (r with the simulated log2 FC 0.867, 0.871 before) | landscape |
| InMoose 0.9.1 pycombat_seq | InMoose | 20k | 44 s | 0.86 | batch removed 0.861 (residual batch effect vs the simulated one: r 0.122, 0.879 before) x group effect kept 0.995 (r with the simulated log2 FC 0.867, 0.871 before) | landscape |

### `cpm`
{: .no_toc }

Scored at n = 60k; fastest 11 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Expression.cpm** | genomeblocks | 60k | 11 ms | 1.00 | reference | landscape |
| rnanorm 2.2.0 | RNAnorm | 60k | 2.9 s | 1.00 | max |Δ| 1.75e-10 | landscape |

### `de_test`
{: .no_toc }

Scored at n = 60k; fastest 9.0 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| limma-voom | limma | 60k | 9.0 s | 0.99 | AUROC 0.9936 against the simulated DE genes; at BH 5 %: 5,867 called, 95.5% true, 95.7% of the true DE genes found | landscape |
| InMoose 0.9.1 deseq2 | InMoose | 60k | 29 s | 0.99 | AUROC 0.9946 against the simulated DE genes; at BH 5 %: 6,016 called, 94.3% true, 96.9% of the true DE genes found | landscape |
| edgeR QL | edgeR | 60k | 45 s | 1.00 | AUROC 0.9961 against the simulated DE genes; at BH 5 %: 6,076 called, 94.0% true, 97.6% of the true DE genes found | landscape |
| DESeq2 | DESeq2 | 60k | 81 s | 0.99 | AUROC 0.9946 against the simulated DE genes; at BH 5 %: 6,016 called, 94.3% true, 96.9% of the true DE genes found | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 60k | 6.1 min | 0.99 | AUROC 0.9946 against the simulated DE genes; at BH 5 %: 6,015 called, 94.4% true, 96.9% of the true DE genes found | landscape |
| InMoose 0.9.1 edgepy | InMoose | 60k | 33.0 min | 0.99 | AUROC 0.9925 against the simulated DE genes; at BH 5 %: 7,648 called, 74.5% true, 97.3% of the true DE genes found | landscape |
| limma-voom | limma | 20k | 1.8 s | 0.97 | AUROC 0.9716 against the simulated DE genes; at BH 5 %: 1,638 called, 96.4% true, 80.4% of the true DE genes found | landscape |
| InMoose 0.9.1 deseq2 | InMoose | 20k | 3.4 s | 0.97 | AUROC 0.9728 against the simulated DE genes; at BH 5 %: 1,801 called, 92.6% true, 84.9% of the true DE genes found | landscape |
| edgeR QL | edgeR | 20k | 6.9 s | 0.98 | AUROC 0.9775 against the simulated DE genes; at BH 5 %: 1,733 called, 95.1% true, 83.9% of the true DE genes found | landscape |
| DESeq2 | DESeq2 | 20k | 12 s | 0.97 | AUROC 0.9728 against the simulated DE genes; at BH 5 %: 1,801 called, 92.6% true, 84.9% of the true DE genes found | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 20k | 81 s | 0.97 | AUROC 0.9728 against the simulated DE genes; at BH 5 %: 1,802 called, 92.6% true, 84.9% of the true DE genes found | landscape |
| InMoose 0.9.1 edgepy | InMoose | 20k | 5.4 min | 0.97 | AUROC 0.9718 against the simulated DE genes; at BH 5 %: 1,721 called, 93.7% true, 82.1% of the true DE genes found | landscape |

### `dispersion`
{: .no_toc }

Scored at n = 60k; fastest 23 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| InMoose 0.9.1 deseq2 | InMoose | 60k | 23 s | 0.95 | Spearman 0.9538 with the simulated dispersions (60,000 genes) | landscape |
| edgeR estimateDisp | edgeR | 60k | 38 s | 0.95 | Spearman 0.9537 with the simulated dispersions (60,000 genes) | landscape |
| DESeq2 | DESeq2 | 60k | 56 s | 0.95 | Spearman 0.9538 with the simulated dispersions (60,000 genes) | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 60k | 4.3 min | 0.95 | Spearman 0.9538 with the simulated dispersions (60,000 genes) | landscape |
| InMoose 0.9.1 deseq2 | InMoose | 20k | 2.6 s | 0.86 | Spearman 0.8561 with the simulated dispersions (19,997 genes) | landscape |
| DESeq2 | DESeq2 | 20k | 5.1 s | 0.86 | Spearman 0.8561 with the simulated dispersions (19,997 genes) | landscape |
| edgeR estimateDisp | edgeR | 20k | 5.7 s | 0.85 | Spearman 0.8537 with the simulated dispersions (20,000 genes) | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 20k | 63 s | 0.86 | Spearman 0.8561 with the simulated dispersions (19,997 genes) | landscape |

### `filter_by_expr`
{: .no_toc }

Scored at n = 60k; fastest 55 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 60k | 55 ms | – |  | expression |
| R edgeR filterByExpr | edgeR | 60k | 107 ms | 1.00 |  | expression |

### `lfc_shrinkage`
{: .no_toc }

Scored at n = 60k; fastest 4.4 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| DESeq2 lfcShrink (ashr) | DESeq2 | 60k | 4.4 s | 0.98 | Pearson 0.9843 with the simulated log2 FC; RMSE 0.086 | landscape |
| DESeq2 lfcShrink (normal) | DESeq2 | 60k | 24 s | 0.95 | Pearson 0.9465 with the simulated log2 FC; RMSE 0.159 | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 60k | 32 s | 0.96 | Pearson 0.9583 with the simulated log2 FC; RMSE 0.141 | landscape |
| DESeq2 lfcShrink (apeglm) | DESeq2 | 60k | 33 s | 0.96 | Pearson 0.9584 with the simulated log2 FC; RMSE 0.141 | landscape |
| DESeq2 lfcShrink (ashr) | DESeq2 | 20k | 1.5 s | 0.93 | Pearson 0.9328 with the simulated log2 FC; RMSE 0.175 | landscape |
| DESeq2 lfcShrink (normal) | DESeq2 | 20k | 6.5 s | 0.85 | Pearson 0.8498 with the simulated log2 FC; RMSE 0.266 | landscape |
| DESeq2 lfcShrink (apeglm) | DESeq2 | 20k | 8.6 s | 0.89 | Pearson 0.8897 with the simulated log2 FC; RMSE 0.229 | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 20k | 11 s | 0.89 | Pearson 0.8897 with the simulated log2 FC; RMSE 0.229 | landscape |

### `log_cpm`
{: .no_toc }

Scored at n = 60k; fastest 60 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 60k | 60 ms | – |  | expression |
| R edgeR cpm(log=TRUE) | edgeR | 60k | 207 ms | 1.00 |  | expression |

### `read_counts_bam`
{: .no_toc }

Scored at n = 10M; fastest 6.0 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| Subread featureCounts (CLI) | Rsubread featureCounts | 10M | 6.0 s | 1.00 | 20,000/20,000 genes equal the simulated counts | landscape |
| Rsubread featureCounts (GTF) | Rsubread featureCounts | 10M | 8.0 s | 1.00 | 20,000/20,000 genes equal the simulated counts | landscape |
| GenomicAlignments summarizeOverlaps (Union) | GenomicAlignments | 10M | 13 s | 1.00 | 19,936/20,000 genes equal the simulated counts | landscape |
| HTSeq htseq-count (CLI) | HTSeq | 10M | 4.7 min | 1.00 | 20,000/20,000 genes equal the simulated counts | landscape |
| GenomicAlignments summarizeOverlaps (Union) | GenomicAlignments | 1M | 1.7 s | 1.00 | 19,964/20,000 genes equal the simulated counts | landscape |
| Subread featureCounts (CLI) | Rsubread featureCounts | 1M | 1.9 s | 1.00 | 20,000/20,000 genes equal the simulated counts | landscape |
| Rsubread featureCounts (GTF) | Rsubread featureCounts | 1M | 2.3 s | 1.00 | 20,000/20,000 genes equal the simulated counts | landscape |
| HTSeq htseq-count (CLI) | HTSeq | 1M | 50 s | 1.00 | 20,000/20,000 genes equal the simulated counts | landscape |

### `read_featurecounts`
{: .no_toc }

Scored at n = 60k; fastest 78 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| R data.table::fread | data.table | 60k | 78 ms | – |  | expression |
| **genomeblocks (polars)** | genomeblocks | 60k | 178 ms | 1.00 |  | expression |
| **genomeblocks (pandas)** | genomeblocks | 60k | 644 ms | 1.00 |  | expression |

### `read_salmon`
{: .no_toc }

Scored at n = 60k; fastest 13 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks (polars)** | genomeblocks | 60k | 13 s | 1.00 |  | expression |
| R tximport (txOut) | tximport | 60k | 25 s | 1.00 |  | expression |
| **genomeblocks (pandas)** | genomeblocks | 60k | 52 s | 1.00 |  | expression |

### `size_factors`
{: .no_toc }

Scored at n = 60k; fastest 55 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| pydeseq2 0.5.4 | PyDESeq2 | 60k | 55 ms | 1.00 | max |Δ| 1.33e-15 | landscape |
| **genomeblocks Expression.size_factors** | genomeblocks | 60k | 72 ms | 1.00 | reference (equals DESeq2 estimateSizeFactors, see the Expression bench) | landscape |
| **genomeblocks** | genomeblocks | 60k | 205 ms | – |  | expression |
| R DESeq2 | DESeq2 | 60k | 681 ms | 1.00 |  | expression |

### `tmm`
{: .no_toc }

Scored at n = 60k; fastest 530 ms.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks Expression.tmm** | genomeblocks | 60k | 530 ms | 1.00 | reference (equals edgeR calcNormFactors to 4e-16, see the Expression bench) | landscape |
| **genomeblocks** | genomeblocks | 60k | 1.7 s | – |  | expression |
| rnanorm 2.2.0 | RNAnorm | 60k | 2.3 s | 1.00 | max |Δ| 4.44e-16 | landscape |
| R edgeR calcNormFactors | edgeR | 60k | 5.3 s | 1.00 |  | expression |

### `tx_to_gene`
{: .no_toc }

Scored at n = 60k; fastest 1.1 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| **genomeblocks** | genomeblocks | 60k | 1.1 s | – |  | expression |
| R tximport::summarizeToGene | tximport | 60k | 2.8 s | 1.00 |  | expression |
| **genomeblocks Expression.make + to_genes (length_scaled_tpm)** | genomeblocks | 60k | 3.9 s | 1.00 | reference (equals R tximport to 1e-9, see the Expression bench) | landscape |
| pytximport 0.13.0 | pytximport | 60k | 7.7 s | 1.00 | max |Δ| 9.02e-10 | landscape |

### `vst`
{: .no_toc }

Scored at n = 60k; fastest 1.9 s.

| Engine (as run) | Tool | n | time | accuracy | detail | source |
|---|---|---|---|---|---|---|
| DESeq2 vst | DESeq2 | 60k | 1.9 s | 1.00 | reference | landscape |
| InMoose 0.9.1 deseq2 | InMoose | 60k | 12 s | 1.00 | Pearson 0.99997, max |Δ| 0.151 | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 60k | 2.3 min | 1.00 | Pearson 0.99997, max |Δ| 0.154 | landscape |
| DESeq2 vst | DESeq2 | 20k | 913 ms | 1.00 | reference | landscape |
| InMoose 0.9.1 deseq2 | InMoose | 20k | 1.3 s | 1.00 | Pearson 0.99999, max |Δ| 0.0644 | landscape |
| PyDESeq2 0.5.4 | PyDESeq2 | 20k | 37 s | 1.00 | Pearson 1.00000, max |Δ| 0.0619 | landscape |

### Perturbation tests
{: .no_toc }

| Tool | Engine | task | variant | identical answer |
|---|---|---|---|---|
| genomeblocks | genomeblocks | size_factors | samples shuffled | yes |
| genomeblocks | genomeblocks | tmm | samples shuffled | yes |
| PyDESeq2 | pydeseq2 0.5.4 | size_factors | genes shuffled | yes |
| PyDESeq2 | pydeseq2 0.5.4 | size_factors | samples shuffled | yes |
| RNAnorm | rnanorm 2.2.0 | tmm | genes shuffled | yes |
| RNAnorm | rnanorm 2.2.0 | tmm | samples shuffled | yes |

### Every tool in this theme
{: .no_toc }

Capability checklist: `read_salmon_kallisto`, `read_rsem`, `read_featurecounts`, `read_star_htseq`, `tx_to_gene`, `tmm`, `size_factors`, `cpm_logcpm`, `filter_by_expr`, `de_test`, `lfc_shrinkage`, `aligned_to_gene_models`, `anndata_or_sce`, `batch_correction`, `python_native`, `r_native`, `multi_format_sniffing`, `parallel`, `count_reads_from_bam`, `dispersion_estimation`, `variance_stabilising_transform`, `sample_qc_pca`, `ma_volcano_plots`, `independent_filtering_fdr`, `gene_set_tests`, `pseudobulk`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 54% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| edgeR | R | library | conda · Bioc | 2010 | 4.10.7 2026-04 | yes | 62% | 1.00 | TMM/RLE/upperquartile in calcNormFactors, cpm/logCPM, filterByExpr, QL F-test; catchSalmon/catchKallisto/catchRSEM read quant dirs (transcript-level, with overdispersion), featureCounts2DGEList and readDGE for count files; no tximport-style gene summarising. |
| RNAlysis | Python | mixed | pip | 2023 | 4.3.0 2026-08 | yes | 62% | 1.00 | GUI-first (no-code) with a Python API; TMM / median-of-ratios / RPM normalisation and low-read filters native, DESeq2 and limma-voom run through an R installation. |
| DESeq2 | R | library | conda · Bioc | 2014 | 1.52.0 2026-04 | yes | 50% | 1.00 | Median-of-ratios size factors, NB GLM Wald/LRT, lfcShrink (apeglm/ashr/normal), vst/rlog; DESeqDataSetFromHTSeqCount reads htseq-count files, salmon/kallisto come in via tximport; BiocParallel. |
| iDEP | R | web | web | 2018 | v2.4.0 2025-11 | yes | 50% | 0.62 | Web (Shiny) front end over edgeR/DESeq2/limma with pathway analysis; last commit 2026-10; GitHub About shows no license field. |
| DEBrowser | R | web | Bioc | 2019 | 1.40.0 2026-04 | yes | 46% | 0.90 | Interactive DESeq2/edgeR/limma front end with ComBat; active only via Bioconductor version bumps (upstream last commit 2024-02). |
| nf-core/rnaseq | mixed | cli | src | 2020 | 3.27.0 2026-09 | yes | 42% | 0.80 | Pipeline, not a library: runs salmon/kallisto/RSEM/STAR, tximport summarising and DESeq2 QC against the supplied GTF; DE lives in nf-core/differentialabundance. |
| InMoose | Python | library | pip | 2023 | 0.9.1 2026-01 | yes | 38% | 1.00 | Ports of edgeR (partial), DESeq2 and limma plus pycombat_seq (ComBat-seq) and cohort QC in one Python package; the only Python package covering TMM, size factors, QL/NB DE and batch correction together. |
| OmicVerse | Python | library | pip · conda | 2024 | 2.3.2 2026-09 | yes | 38% | 1.00 | Umbrella bulk+single-cell framework; ov.bulk.pyDEG runs DESeq2 (PyDESeq2), t-test or Wilcoxon, plus ComBat batch correction and gene-ID mapping. |
| PyDESeq2 | Python | library | pip · conda | 2023 | 0.5.4 2026-01 | yes | 38% | 1.00 | Pure-Python DESeq2 (Wald test, multi-factor designs via formulaic, apeglm-style lfc_shrink); takes a pandas/AnnData count matrix, no quantifier readers; Python &gt;=3.11. |
| sleuth | R | library | conda | 2017 | v0.30.2 2025-05 | no | 38% | 0.89 | Reads kallisto HDF5 bootstraps (salmon via wasabi), gene_mode aggregation, measurement-error model with Wald/LRT and a Shiny explorer; GitHub-only R package, last commit 2025-05. |
| Degust | mixed | web | web | 2019 | v4.1.1 2019-10 | yes | 35% | 0.50 | Monash web tool running edgeR QL / limma-voom with interactive MA/parallel-coordinate views; commits to 2026-05 but no tagged release since 2019. |
| limma | R | library | conda · Bioc | 2015 | 3.68.5 2026-04 | yes | 35% | 1.00 | voom (2014, 10.1186/gb-2014-15-2-r29) + empirical-Bayes linear models and limma-trend; removeBatchEffect; relies on edgeR for TMM/CPM/filtering of counts. |
| DESeq (legacy) | R | library | src | 2010 | 1.42.0 | no | 31% | 0.67 | Origin of median-of-ratios size factors; deprecated and removed from Bioconductor after release 3.12 (last version 1.42.0), superseded by DESeq2. |
| pertpy | Python | library | pip · conda | 2025 | 1.4.0 2026-10 | yes | 31% | 1.00 | Adjacent (perturbation / single-cell): unified DE interface wrapping PyDESeq2, edgeR (via rpy2) and statsmodels on pseudobulks. |
| scanpy | Python | library | pip · conda | 2018 | 1.12.4 2026-08 | yes | 31% | 1.00 | Adjacent (single-cell): normalize_total/log1p is CPM-like, rank_genes_groups is t-test/Wilcoxon/logreg rather than NB GLM, sc.pp.combat for batch. |
| fishpond (swish) | R | library | Bioc | 2019 | 2.18.0 2026-04 | yes | 27% | 1.00 | Non-parametric DE (Wilcoxon / Mann-Whitney over salmon Gibbs/bootstrap replicates); scaleInfReps does median-ratio scaling, labelKeep filters; reads come via tximeta. |
| glmGamPoi | R | library | Bioc | 2020 | 1.24.0 2026-04 | yes | 27% | 1.00 | Fast Gamma-Poisson GLM with quasi-likelihood F-tests and pseudobulk helper; works on-disk / sparse matrices, used by Seurat sctransform. |
| NOISeq | R | library | Bioc | 2011 | 2.56.0 2026-04 | yes | 27% | 1.00 | Non-parametric DE (NOISeq/NOISeqBIO), QC reports, filtered.data (CPM/Wilcoxon/proportion filters), ARSyNseq batch removal; rpkm/uqua/tmm normalisation. |
| baySeq | R | library | Bioc | 2010 | 2.46.0 2026-04 | yes | 23% | 1.00 | Empirical-Bayes posterior likelihoods of DE patterns across groups; snow-cluster parallelism; getLibsizes offers quantile/edgeR scaling. |
| decoupler | Python | library | pip · conda | 2022 | 2.2.0 2026-07 | yes | 23% | 1.00 | Adjacent: TF/pathway enrichment; ships a pseudobulk helper and an edgeR filter_by_expr port, DE itself is delegated to PyDESeq2 in its tutorials. |
| diffxpy | Python | library | pip | – | 0.7.4 2020-03 | no | 23% | 0.75 | NB GLM Wald/LRT/t-test/Wilcoxon for single-cell; README says unmaintained and recommends PyDESeq2; repo archived 2025-10-21. |
| EBSeq | R | library | Bioc | 2013 | 2.10.0 2026-04 | yes | 23% | 1.00 | Empirical-Bayes isoform/gene DE (RSEM&#x27;s companion); MedianNorm size factors and PostFC posterior fold changes; Bioc version bumps only. |
| edgePython | Python | library | pip | – | 0.5.0 2026-10 | yes | 23% | 0.83 | Pachter-lab Python port of edgeR (TMM, filter_by_expr, glm_ql_fit / glm_ql_ftest, NEBULA-LN) released 2026-10-22; too new to have docs site or journal paper. |
| GEO2R | Web | web | web | 2013 | – | ? | 23% | 0.50 | NCBI web tool: limma for arrays, DESeq2 (beta) on NCBI-generated RNA-seq counts, exports the R script; help page is behind a captcha so details are from search snippets. |
| tximeta | R | library | Bioc | 2020 | 1.30.0 2026-04 | yes | 23% | 1.00 | tximport + automatic transcriptome identification by index checksum; attaches GRanges/TxDb gene models and summarizeToGene(); the closest R analogue of GTF-aligned expression tables. |
| dearseq | R | library | Bioc | 2020 | 1.24.0 2026-04 | yes | 19% | 1.00 | Variance-component score / permutation tests for gene and gene-set DE (robust to small n and longitudinal designs). |
| pytximport | Python | mixed | pip · conda | 2024 | 0.13.0 2026-02 | yes | 19% | 1.00 | Python tximport (salmon, kallisto, RSEM, piscem/oarfish, generic tsv) with CLI and library; outputs AnnData or BiocPy SummarizedExperiment; helpers to build tx2gene maps. |
| bioinfokit | Python | library | pip | 2020 | 2.1.4 2024-05 | no | 15% | 0.62 | General toolkit with CPM/RPKM/TPM helpers and volcano/MA plots; no DE model; last commit 2024-05. |
| DEGseq | R | library | Bioc | 2010 | 1.66.0 2026-04 | yes | 15% | 1.00 | Legacy MA-plot based tests (Fisher/LRT) without a replicate dispersion model; kept alive by Bioconductor version bumps. |
| HTSeq | Python | mixed | pip · conda | 2015 | 2.1.2 2026-02 | yes | 15% | 1.00 | Adjacent: produces the htseq-count tables others read (rows = GTF genes); HTSeq 2.0 paper 10.1093/bioinformatics/btac166; -n threads. |
| Rsubread featureCounts | R | library | conda · Bioc | 2014 | 2.26.0 2026-04 | yes | 15% | 1.00 | Adjacent counting step (produces the featureCounts tables); Rsubread paper 10.1093/nar/gkz114; nthreads. |
| tximport | R | library | conda · Bioc | 2015 | 1.40.0 2026-04 | yes | 15% | 1.00 | Reference implementation of lengthScaledTPM / scaledTPM transcript-to-gene summarising; quantifier type must be given explicitly (type=); also reads sailfish/StringTie/alevin. |
| BiocPy summarizedexperiment | Python | library | pip | – | 0.8.0 2026-09 | yes | 12% | 0.90 | Python port of the Bioconductor container (BiocPy; 2023 BioC conference abstract, no journal paper found); converts to/from AnnData; pytximport can emit it. |
| conorm | Python | library | pip | – | 1.2.0 2021-10 | no | 12% | 0.29 | Small TMM / MRN (median-of-ratios) / GeTMM normaliser; no license metadata, dormant since 2021. |
| RNAnorm | Python | mixed | pip · conda | – | 2.2.0 2025-01 | no | 12% | 0.80 | CPM, FPKM, TPM, UQ, CUF, TMM, CTF as sklearn transformers plus a CLI; normalisation only, last commit 2025-04. |
| RUVSeq | R | library | Bioc | 2014 | 1.46.0 2026-04 | yes | 12% | 1.00 | Removes unwanted variation (RUVg/RUVs/RUVr) as covariates for edgeR/DESeq2; upstream GitHub last commit 2022-11, Bioc bumps only. |
| SummarizedExperiment | R | library | Bioc | 2015 | 1.42.0 2026-04 | yes | 12% | 1.00 | The container every Bioconductor DE tool consumes; rowRanges aligns assays to gene models, which is the design genomeblocks&#x27; gene-row expression table mirrors. |
| apeglm | R | library | Bioc | 2019 | 1.34.0 2026-04 | yes | 8% | 1.00 | Approximate posterior estimation for GLM coefficients (heavy-tailed prior); default shrinker for DESeq2 and the model PyDESeq2 re-implements. |
| ashr | R | library | CRAN | 2017 | 2.2-63 2023-08 | yes | 8% | 1.00 | Empirical-Bayes adaptive shrinkage used by DESeq2::lfcShrink(type=&#x27;ashr&#x27;); CRAN release 2023-08, GitHub commits to 2026-05. |
| sva / ComBat-seq | R | library | Bioc | 2012 | 3.60.0 2026-04 | yes | 8% | 1.00 | Surrogate variable analysis, ComBat and ComBat-seq (NB batch adjustment on counts); the reference implementation that InMoose ports. |
| Wilcoxon rank-sum DE (Li et al. 2022) | R | library | src | 2022 | – | ? | 8% | 0.25 | A method recommendation (Wilcoxon beats DESeq2/edgeR FDR control at large n), not a package; include only as a benchmark baseline row. |
| GenomicAlignments | R | library | conda · Bioc | 2013 | – | yes | 4% | 1.00 | Reads BAM alignments into R; summarizeOverlaps counts reads per feature (Union / IntersectStrict / IntersectNotEmpty modes). Added 2026-10-23 as an engine the landscape benchmark runs. |


## Viewers
{: .sec-navy #viewers }

### Core jobs and grades
{: .no_toc }

Core jobs: **C1** Multi-track static locus figure (signal, genes, peaks); **C2** bigWig signal tracks; **C3** Gene model tracks; **C4** BAM read / coverage tracks; **C5** Loop arcs; **C6** Hi-C contact-matrix track; **C7** Interactive zoom / pan browsing; **C8** Notebook-inline display; **C9** Many regions / multi-panel layouts; **C10** Programmatic figure API.

| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |
|---|---|---|---|---|---|

### Every tool in this theme
{: .no_toc }

Capability checklist: `static_region_plot`, `bigwig_tracks`, `bam_tracks`, `gene_tracks`, `loop_arcs`, `hic_matrix_track`, `interactive_html_export`, `igv_session_or_igvjs`, `notebook_inline`, `config_file_driven`, `cli`, `programmatic_api`, `multi_region_layout`, `graph_network_view`, `web_app`, `interactive_browsing`, `variant_tracks`, `region_highlights`, `linked_views`, `publication_export`.

| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **genomeblocks** | Python | library | pip · conda | 2025 | 2.0.1 2026-10 | yes | 70% | 0.90 | Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1. |
| JBrowse 2 | JavaScript | mixed | npm | 2023 | v4.3.0 (npm latest; v5.0.0-beta.17 tagged 2026-10) 2026-05 | yes | 85% | 1.00 | Web + desktop + embedded React browser; @jbrowse/cli manages configs, @jbrowse/img (jb2export) renders static SVG, .hic and bedpe arcs supported; v4.3.0 release date year inferred from the release ordering (page shows day/month only). |
| Gosling / gos | mixed | library | pip · npm | 2022 | gosling.js 1.0.7 npm latest (2.0.0-alpha.10 tag 2026-03); gos 0.3.0 2026-03 | yes | 80% | 1.00 | Grammar-of-graphics for genomics (HiGlass-backed rendering): linked multi-view, circular layouts, bigWig/BAM/matrix/bedpe; gos Python emits standalone HTML; gos last PyPI release 2025-05. |
| HiGlass | mixed | web | pip | 2018 | v2.3.5 2026-06 | yes | 75% | 0.90 | Viewer; higlass-python 1.4.0 (2026-03) for notebooks; Multiscale interactive Hi-C + track browser with linked views; Jupyter via higlass-python; needs a tile server (not a single-file export); SVG export of views. |
| IGV desktop | Java | gui | conda | 2011 | v3.0.0-beta.7 (stable line 2.19.x) 2026-09 | yes | 70% | 0.80 | Interactive desktop browser; batch scripts / port commands give snapshots (PNG/SVG), multi-locus view, bedpe arcs and .hic contact-map tracks; latest tag is a 3.0 beta. |
| igv.js | JavaScript | library | npm | 2023 | 3.8.10 2026-10 | yes | 70% | 1.00 | Embeddable browser component behind igv-web, igv-notebook and igv-reports; JSON sessions, multi-locus, interact/bedpe arcs, SVG export; Hi-C matrices live in the sibling juicebox.js. |
| UCSC Genome Browser | Web | web | web | 2002 | – | yes | 70% | 0.62 | Reference web browser: PDF/SVG view export, interact/bigInteract arcs, .hic tracks, multi-region mode, REST API and track hubs; kent repo commits to 2026-10 (no tagged releases). |
| igv-notebook | Python | library | pip | – | 0.6.2 2025-02 | yes | 65% | 0.67 | Thin Python wrapper that embeds igv.js (any igv.js track, session JSON load, SVG save) in Jupyter and Colab; commits to 2026-10, release 2025-02. |
| trackViewer | R | library | Bioc | 2019 | 1.48.0 2026-04 | yes | 65% | 1.00 | Lollipop/dandelion plots, interaction heatmaps and loopBouquetPlot (a 3D-ish contact-network layout), browseTracks htmlwidget for interactive HTML. |
| WashU Epigenome Browser | JavaScript | web | npm | 2011 | – | yes | 65% | 0.71 | Rewritten 3rd-generation codebase (eg3, 2025 NAR update); longrange/bedpe arcs, .hic/cool heatmaps, gene-set multi-region view, SVG screenshots; no tagged releases, commits to 2026-08. |
| CoolBox | Python | mixed | pip · conda | 2021 | 0.4.0 2025-11 | yes | 60% | 1.00 | pyGenomeTracks-derived, Jupyter-first: Frame + Track composition with an interactive widget browser; .cool/.hic, arcs, BAM, bigWig, GTF. |
| GenomeSpy | JavaScript | mixed | npm | 2024 | v1.2.0 2026-10 | yes | 55% | 1.00 | Grammar-based, GPU-accelerated cohort/sample-faceted genome views (cancer genomics focus); embeds in any page; PNG export rather than figure-grade static output. |
| ggcoverage | R | library | src | 2023 | v1.3.0 2024-04 | no | 55% | 0.70 | ggplot2 grammar for BAM/bigWig/bedGraph coverage + GTF genes + peaks + bedpe links + Hi-C matrix; removed from CRAN 2023-11, GitHub-only, last commit 2025-01. |
| Juicebox | Java | gui | web | 2016 | v2.20.00 2022-10 | yes | 55% | 0.62 | Viewer; last commit 2026-09; Hi-C-centric viewer (.hic) with 1D bigWig/gene tracks and 2D loop/domain annotations, image export, shareable sessions; desktop commits to 2026-09, juicebox.js web. |
| pyGenomeTracks | Python | cli | pip · conda | 2020 | 3.9 2024-07 | yes | 55% | 0.90 | Viewer; last commit 2026-08; De-facto standard static track figure (bigwig, bedgraph, GTF/bed12 genes, links arcs, cool/h5 Hi-C triangles, TADs); no BAM track; last PyPI release 2024-07 but commits to 2026-08. |
| ggbio | R | library | Bioc | 2012 | 1.60.0 2026-04 | yes | 50% | 1.00 | autoplot() for GRanges/TxDb/BamFile, geom_arch links, tracks() stacking and circular layouts; maintenance mode, Bioc bumps. |
| karyoploteR | R | library | Bioc | 2017 | 1.38.0 2026-04 | yes | 50% | 1.00 | Base-graphics ideograms with zoom; kpPlotBigWig, kpPlotBAMDensity, kpPlotGenes (TxDb), kpPlotLinks arcs; whole-genome or region. |
| plotgardener | R | library | Bioc | 2022 | 1.18.0 | yes | 50% | 1.00 | Viewer used by the OHCA/HiContacts stack for arcs; GitHub activity not checked; Sushi&#x27;s successor: absolute-coordinate multi-panel pages (plotHicTriangle/Square, plotPairsArches, plotSignal, plotGenes from TxDb/OrgDb) for publication figures. |
| pyGenomeViz | Python | mixed | pip · conda | – | 1.7.0 2026-06 | yes | 50% | 0.90 | Prokaryote / comparative-genomics feature maps with synteny links (MUMmer/BLAST), HTML export and a Streamlit app; not an NGS signal-track viewer. |
| trackc | Python | mixed | pip | – | 0.0.28 2026-06 | yes | 50% | 0.90 | Newer matplotlib-native 3D-genome track plotter aimed at rearranged genomes (multi-region/zoom-in layouts, Virtual4C); 0.0.x versions, no paper found. |
| Circos | Perl | cli | conda | 2009 | 0.69-9 (bioconda 0.69.9) | no | 45% | 0.62 | Adjacent circular ideogram/link plots (contacts as ribbons = a network-like view); upstream release date not verified (circos.ca TLS failure), bioconda package updated 2025-04. |
| FAN-C (fancplot) | Python | mixed | pip · conda | 2020 | 0.9.28 2024-03 | no | 45% | 0.70 | fancplot CLI and fanc.plotting API: triangular/square Hi-C, bigWig line/bar panels, GTF gene panel, loop/TAD annotations; last commit 2024-03. |
| LocusZoom | JavaScript | library | CRAN · npm | 2010 | v0.14.0 2022-03 | no | 45% | 0.90 | Adjacent (variant association plots with a gene track); LocusZoom.js 2021 Bioinformatics 10.1093/bioinformatics/btab186; R port locuszoomr on CRAN makes static ggplot versions; JS commits to 2025-05. |
| Gviz | R | library | Bioc | 2016 | 1.56.0 2026-04 | yes | 40% | 1.00 | Classic R track plotting (DataTrack from bigWig, AlignmentsTrack from BAM, GeneRegionTrack from TxDb/biomaRt, ideogram); no native arcs or Hi-C. |
| HiCExplorer (hicPlotMatrix / hicPlotTADs) | Python | cli | pip · conda | 2018 | 3.7.6 2024-11 | no | 40% | 0.80 | Hi-C processing suite whose hicPlotTADs delegates to pyGenomeTracks and hicPlotMatrix draws square/triangular heatmaps; PyPI listing is stale (2019), bioconda current; last commit 2024-11. |
| HiCPlotter | Python | cli | src | 2015 | – | no | 40% | 0.30 | Python 2.7 script for side-by-side Hi-C heatmaps with bedGraph/bed/arc tracks; last commit 2020-07, no releases. |
| jbrowse-jupyter | Python | library | pip | – | 1.5.9 2023-12 | yes | 40% | 0.88 | Deprecated 2026-07 in favour of jbrowse-anywidget (per its last commit); listed for completeness as the JBrowse notebook path. |
| rGenomeTracks | R | library | Bioc | – | 1.18.0 2026-04 | yes | 40% | 0.89 | R wrapper around pyGenomeTracks (needs the Python tool installed); inherits its track set and single-region output. |
| Sushi | R | library | src | 2014 | 1.32.0 | no | 40% | 0.56 | Deprecated (README points to plotgardener); removed from Bioconductor after release 3.14; plotBedpe arcs / plotHic from data.frames, bedGraph not bigWig. |
| genomeview | Python | library | pip | 2018 | 1.1.4 2023-08 | no | 35% | 0.80 | Lightweight SVG renderer strong on BAM (short/long reads, quick consensus, grouping by tag); bigWig via pyBigWig, annotations as BED; dormant since 2023. |
| gghic | R | library | src | 2024 | – | yes | 35% | 0.70 | ggplot2 extension for Hi-C-like data with loop/TAD/gene/signal annotation layers; GitHub-only install (CRAN-SUBMISSION file present), last commit 2026-01. |
| igv-reports | Python | cli | pip · conda | 2026 | 1.17.0 2026-09 | yes | 35% | 0.78 | create_report builds one standalone HTML per site list (VCF/BED/MAF/BEDPE) with BAM/CRAM/VCF/BED/GTF/WIG/bedGraph tracks, no dependency on input files; bigWig not listed as a track input. |
| svist4get | Python | mixed | pip · conda | 2019 | 1.3.1.1 2024-07 | no | 35% | 0.88 | Publication-style locus images from bedGraph tracks + GTF + FASTA (sequence shown at base resolution); last commit 2024-07. |
| Biopython GenomeDiagram | Python | library | pip · conda | 2006 | 1.88 2026-08 | yes | 30% | 1.00 | Feature maps (linear/circular) with cross-links from GenBank records; graph tracks from arbitrary data; no NGS file readers. |
| GENOVA | R | library | src | 2021 | v1.0.0-TheCartographer 2021-11 | no | 30% | 0.75 | Loops consumed as BEDPE for APA; no caller; last commit 2024-03; GitHub install only; Hi-C analysis + pyramid/hic.matrixplot figures with loop/TAD overlays (APA, saddle, insulation); last commit 2024-03. |
| trackplot | R | library | src | 2024 | – | no | 30% | 0.50 | Minimal-dependency R script layered on bwtool; IGV-style locus plots plus profiles/heatmaps; GitHub push 2024-10; DOI inferred from the article id vbae031; Minimal-dependency IGV-style locus plots, profile plots and heatmaps from bigWig (needs bwtool); gene models from UCSC/GTF; last commit 2024-10, no tagged releases. Distinct from the Python sashimi tool of the same name. |
| tracktables | R | library | Bioc | – | 1.46.0 2026-04 | yes | 25% | 0.89 | Builds IGV sessions and HTML region tables that open loci in a running IGV; export-to-IGV helper rather than a renderer. |
| ChIPseeker | R | library | Bioc | 2015 | 1.48.0 2026-04 | yes | 20% | 1.00 | Most-cited R peak annotator: tssRegion-configurable promoter class, plotAnnoPie/Bar, plotDistToTSS, covplot, seq2gene + clusterProfiler enrichment; GTF must first become a TxDb (txdbmaker); 2022 Current Protocols update (10.1002/cpz1.585); Profiles/heatmaps are of peak occupancy (optionally weighted), not of bigWig signal; primarily a peak annotator; Adjacent: covplot draws genome-wide peak coverage per chromosome and TSS profiles; peak annotation, not a locus track viewer. |


## Usability checklist
{: .sec-navy #usability }

| Tool | packaged | docs site | tutorial | API reference | tests | CI | active (12 mo) | peer-reviewed | open license | library API | score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| genomeblocks | yes | yes | yes | yes | yes | yes | yes | no | yes | yes | 0.90 |
| ailist (Python) | yes | ? | no | ? | yes | yes | yes | no | yes | yes | 0.75 |
| BEDOPS | yes | yes | yes | yes | yes | yes | no | yes | yes | no | 0.80 |
| bedtools | yes | yes | yes | yes | yes | yes | yes | yes | yes | no | 0.90 |
| BiocPy genomicranges / iranges | yes | yes | yes | yes | yes | yes | yes | no | yes | yes | 0.90 |
| bioframe | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| Biopython | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| Biostrings | yes | yes | yes | yes | yes | yes | yes | no | yes | yes | 0.90 |
| ChIPpeakAnno | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| ChIPseeker | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| chromVAR | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| cooler | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| cooltools | yes | yes | yes | yes | yes | yes | no | yes | yes | yes | 0.90 |
| csaw | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| data.table | yes | yes | yes | yes | yes | yes | yes | no | yes | yes | 0.90 |
| deepTools | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| DESeq2 | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| edgeR | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| EnrichedHeatmap | yes | yes | yes | yes | ? | yes | yes | yes | yes | yes | 1.00 |
| FitHiChIP | no | yes | yes | no | ? | ? | no | yes | yes | no | 0.50 |
| genomation | yes | yes | yes | yes | yes | yes | no | yes | yes | yes | 0.90 |
| GenomicAlignments | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| GenomicFeatures + txdbmaker | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| GenomicInteractions | yes | no | yes | yes | ? | ? | no | yes | yes | yes | 0.75 |
| GenomicRanges + IRanges | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| GIGGLE | yes | no | no | no | yes | yes | yes | yes | yes | yes | 0.70 |
| HiCExplorer | yes | yes | yes | yes | yes | yes | no | yes | yes | no | 0.80 |
| hictk | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| HOMER | yes | yes | yes | yes | ? | ? | yes | yes | yes | no | 0.88 |
| HTSeq | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| igraph | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| InMoose | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| intervaltree | yes | no | no | no | yes | yes | yes | no | yes | yes | 0.60 |
| Juicer / Juicer Tools | no | yes | yes | no | ? | ? | no | yes | yes | no | 0.50 |
| limma | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| LOLA | yes | yes | yes | yes | yes | yes | no | yes | yes | yes | 0.90 |
| MAPS | no | no | yes | no | ? | ? | no | yes | ? | no | 0.29 |
| megadepth | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| MEME suite | yes | yes | yes | no | yes | ? | yes | yes | no | yes | 0.78 |
| MOODS | yes | no | yes | no | ? | no | no | yes | yes | yes | 0.56 |
| motifmatchr | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| Mustache | yes | no | yes | no | ? | ? | no | yes | yes | no | 0.50 |
| plyranges | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| polars-bio | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| pybigtools | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| pyBigWig | yes | no | no | yes | yes | yes | yes | no | yes | yes | 0.70 |
| pychromVAR | yes | yes | yes | yes | ? | ? | ? | ? | yes | yes | 1.00 |
| PyDESeq2 | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| pyranges (0.x) | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| pyranges (1.x) | yes | yes | yes | yes | yes | yes | yes | no | yes | yes | 0.90 |
| pysam | yes | yes | yes | yes | yes | yes | yes | no | yes | yes | 0.90 |
| pytximport | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| regioneR | yes | no | yes | yes | yes | ? | no | yes | yes | yes | 0.78 |
| RNAnorm | yes | yes | yes | yes | yes | yes | no | no | yes | yes | 0.80 |
| Rsubread featureCounts | yes | yes | yes | yes | ? | yes | yes | yes | yes | yes | 1.00 |
| rtracklayer | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| superintervals | yes | no | no | no | yes | yes | yes | no | yes | yes | 0.60 |
| sva / ComBat-seq | yes | yes | yes | yes | ? | yes | yes | yes | yes | yes | 1.00 |
| tximport | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| UCSC kent utilities | yes | yes | no | no | ? | ? | yes | yes | no | no | 0.50 |
| universalmotif | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| UROPA | yes | yes | yes | no | no | no | no | yes | yes | no | 0.50 |
| valr | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | 1.00 |
| WiggleTools | yes | no | yes | no | yes | ? | yes | yes | yes | no | 0.67 |

## Reproduce
{: .sec-navy #reproduce }

```bash
cd benchmarks
python make_data.py                       # the synthetic hg38-shaped data
GB_LANDSCAPE_R=/path/to/Rscript \
GB_LANDSCAPE_PY=/path/to/python \
GB_LANDSCAPE_CLI=/path/to/cli-env/bin GB_LANDSCAPE_HIC=/path/to/hic-env \
GIGGLE=/path/to/giggle HOMER_BIN=/path/to/homer/bin \
MUSTACHE=/path/to/mustache JUICER_TOOLS=/path/to/juicer_tools_2.13.06.jar \
python bench_landscape.py                 # results/landscape.json
python build_landscape.py                 # this page and the figure
```

One part per call (`GB_BENCH_PARTS=intervals`) keeps the other parts' rows;
with `GB_BENCH_OUT=landscape_<part>` each part writes its own file, so the
parts can run as separate cluster jobs, and `landscape/merge_results.py`
folds them back into `results/landscape.json`.

| Variable | What it points at |
|---|---|
| `GB_LANDSCAPE_R` | an Rscript with GenomicRanges, plyranges, valr, data.table, txdbmaker, ChIPseeker, ChIPpeakAnno, rtracklayer, EnrichedHeatmap, genomation, universalmotif, Biostrings, motifmatchr, TFBSTools, chromVAR, LOLA, regioneR, GenomicInteractions, InteractionSet, igraph, Rsubread, GenomicAlignments, csaw, edgeR, limma, DESeq2, apeglm, ashr, sva and tximport |
| `GB_LANDSCAPE_PY` | a Python with pyranges1, polars-bio, genomicranges, superintervals, ailist, pytximport, pydeseq2, inmoose, rnanorm, pychromvar and uropa, plus BEDOPS and `bigWigAverageOverBed` in its `bin/` |
| `GB_LANDSCAPE_CLI` | a `bin/` with featureCounts, htseq-count, megadepth, wiggletools, intervene and bedGraphToBigWig |
| `GB_LANDSCAPE_HIC` | an environment with cooler, cooltools, hictk and HiCExplorer (the `hic` part) |
| `GIGGLE`, `HOMER_BIN` | the GIGGLE binary (built from source) and HOMER's `bin/` |
| `MUSTACHE`, `JUICER_TOOLS`, `JAVA` | Mustache 1.3.3's command, juicer_tools 2.13.06 and the Java to run it (the `hic_loops` part; genomeblocks' ports run without them) |

The four loop callers genomeblocks has ports of are also compared with the
originals on real HiChIP and Hi-C samples by `bench_loops.py`, `bench_maps.py`,
`bench_mustache.py` and `bench_hiccups.py` (see [Benchmarks]({{ '/benchmarks/' | relative_url }})); their
synthetic runs feed the `loops_*` tasks here. The inventory in
`landscape/tools.json` was compiled on 2026-10-22 from each tool's repository
and package index; corrections are welcome as pull requests.
