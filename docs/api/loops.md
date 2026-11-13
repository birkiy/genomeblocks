---
title: loops
parent: API Reference
layout: default
nav_order: 7.5
---

# Loop calling
{: .no_toc }

Four loop callers, each a numpy / scipy reimplementation that reproduces the
original tool's calls on the same input, and is faster: FitHiChIP 11.0
(`fithichip`, HiChIP, from a cooler and anchor peaks), MAPS (`maps`, HiChIP /
PLAC-seq, from a pairs file), Mustache 1.3.3 (`mustache`, Hi-C / HiChIP, from a
balanced cooler) and Juicer's HiCCUPS (`hiccups`, Hi-C, from an .mcool). The
[benchmarks]({{ '/benchmarks/' | relative_url }}#loops) compare each with the
original tool.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
r = gb.fithichip(cool, peaks)          # gb.maps, gb.mustache, gb.hiccups
r.loops                                # a Pairs table of loops
r.write("out", "sample")               # the original tool's output files
```

1. TOC
{:toc}

## `fithichip`

```python
fithichip(cool, peaks, *, resolution=None, chromsizes=None, low=5_000, high=10_000_000, qvalue=0.01,
          nbins=200, bias_range=(0.2, 5.0), neigh=2, verbose=True) -> FitHiChIP
```

| argument | meaning |
|---|---|
| `cool` | a `.cool`, a cooler URI (`file.mcool::resolutions/5000`) or an `.mcool` with `resolution=`. Raw counts are used; balancing weights are ignored (coverage-bias mode). |
| `peaks` | the anchors: a narrowPeak / BED path, a `Loci`, a frame — anything `as_loci` reads. |
| `chromsizes` | the chromosomes to use (FitHiChIP's `ChrSizeFile`); default: the cooler's. |
| `low`, `high` | distance range of the tested bin pairs (bin starts), inclusive. |
| `qvalue` | FDR of a loop (`Q < qvalue`). |
| `nbins` | equal-occupancy groups of the distance model. |
| `bias_range` | biases of the pairs that fit the bias regression. |
| `neigh` | merging neighbourhood, in bins. |

```python
r = gb.fithichip("hichip.mcool", "anchors.narrowPeak", resolution=5000)
# [fithichip] 33,197,566 peak-to-all pairs, 1,222,113 loops (q < 0.01), 439,048 merged | read 7.1s, ...
r.loops                    # Pairs: bin pairs with count, bias1/2, P, Q, ...
r.merged.to_bedpe("loops.bedpe")
A = gb.Architecture.make(cre, r.merged)          # loops straight into a CRE graph
r.write("out", "sample")   # FitHiChIP's two IGV bedpe files (byte for byte) + the tables
```

### The model

1. **Coverage** of a bin: contacts of the pixels whose first bin it is (cis and
   trans) plus the cis pixels whose second bin it is.
2. **Peak bins**: bins touching a peak (closed intervals).
3. **Bias** = coverage / mean nonzero coverage of the bins of the same class
   (peak bins, other bins).
4. **Tested pairs**: cis pixels with a count in [`low`, `high`] and at least
   one peak bin.
5. **Distance model**: equal-occupancy groups along distance; prior =
   contacts / possible peak-to-all pairs / all contacts; R's
   `smooth.spline(cv = TRUE)` (ported: `genomeblocks._smoothspline`), evaluated
   on the bin grid and fitted again.
6. **Bias regression** per group: `log10 count ~ log10 bias1 + log10 bias2`,
   the coefficients smoothed over distance the same way.
7. **Test**: binomial upper tail against the expected count, Benjamini–Hochberg.
8. **Merging**: 8-connected loops form components; in each, loops are taken by
   q (then count) unless both bins lie within `neigh` bins of a taken one.

FitHiChIP's own definitions are kept where they differ from the obvious
choice: the diagonal counts twice in the coverage, bins shorter than 2 bp have
none, and the first bin of the coverage table sorted by name (e.g.
`chr1:0-5000`) is never tested. The [benchmarks]({{ '/benchmarks/' | relative_url }}#loops)
compare the result with FitHiChIP on a synthetic cooler and on HiChIP samples:
the same tested pairs, loops and merged loops; P and q agree to 1e-7–1e-5
(relative), the order of FitHiChIP's own differences between two machines.

## `FitHiChIP`

The result.

| attribute | |
|---|---|
| `tested` | `Pairs` of every tested bin pair, sorted by distance: `count`, `bias1`, `bias2`, `is_peak1`, `is_peak2`, `distance`, `p` (distance prior), `expected`, `p_bias`, `dbinom`, `P`, `Q` |
| `loops` | the rows of `tested` with `Q < qvalue` |
| `merged` | merged loops: `count`, `P`, `Q`, the component's span (`span1_start`, `span1_end`, `span2_start`, `span2_end`), `sum_count`, `strong_conn` (fraction of the span's bin pairs that are loops) |
| `groups` | the equal-occupancy table (FitHiChIP's `Bin_Info.log`) |
| `coefficients` | the bias regression per group |
| `timings` | seconds per step |

`write(outdir, prefix, *, tested=False)` writes
`{prefix}.interactions_FitHiC_Q{q}_IGV.bedpe` and
`Merge_Nearby_Interactions/{prefix}.interactions_FitHiC_Q{q}_MergeNearContacts_IGV.bedpe`
in FitHiChIP's format, the loop tables as TSV and, with `tested=True`, every
tested pair as parquet; it returns the paths.

## `maps`

```python
maps(pairs_or_cool, peaks, *, features, resolution=5000, binning_range=1_000_000, chroms=None, short=None,
     short_cutoff=1000, format="auto", columns=None, strand_columns=None, count_cutoff=12, ratio_cutoff=2.0,
     fdr_cutoff=2.0, gap=15_000, singleton_fdr=1e-4, verbose=True) -> MAPS
```

MAPS (Juric et al. 2019, PLoS Comput Biol 15:e1006982; github.com/ijuric/MAPS) from a pairs file: feather's
split into short (≤ 1 kb, opposite strands) and long pairs, 1D peak bins, AND / XOR bin pairs at 2 bins to
`binning_range`, covariates from MAPS's genomic-features file (effective fragment length, GC, mappability) and
the short-pair count of each bin, a zero-truncated Poisson regression per chromosome and set fitted as VGAM
fits it, MAPS's two-pass refit, its p-value P(Y > c) / P(Y > 0) with the truncated mean, BH per chromosome and
set, the calls (count ≥ 12, observed / expected ≥ 2, fdr < 0.01) and their 15 kb clusters. MAPS's
`*.sig3Dinteractions.bedpe` comes out of `write()`.

| argument | meaning |
|---|---|
| `pairs_or_cool` | HiC-Pro allValidPairs, 4DN `.pairs`, or `id strand1 chr1 pos1 frag1 strand2 chr2 pos2 frag2` (1-based 5' ends; auto-detected, or `columns=` + `strand_columns=`); or a cooler for the long-pair counts, with `short=` |
| `peaks` | the 1D peaks (ChIP-seq or HiChIP-derived), anything `as_loci` reads |
| `features` | MAPS's features table for the genome / enzyme / resolution, e.g. `MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_5Kb_el.GRCh38.txt` |
| `chroms` | default chr1–22, chrX (MAPS's hg38 run) |
| `short` | short-pair counts per bin (`{chrom: array}`) or a pairs file, instead of the input's own |

```python
r = gb.maps("sample.allValidPairs", "peaks.narrowPeak",
            features="F_GC_M_MboI_5Kb_el.GRCh38.txt")
# [maps] 6,804,225 bin pairs tested, 69,814 significant, 61,894 reported | ...
r.loops                       # Pairs: count, expected, fdr, ClusterLabel, ClusterType, ClusterSummit, ...
r.tested                      # every bin pair of the regressions, with covariates, expected, P, fdr
r.write("out", "sample.5k")   # MAPS's sample.5k.2.sig3Dinteractions.bedpe
```

On the same pairs and peaks MAPS tests the same bin pairs with the same counts and calls the same loops;
expected, P and fdr agree to ~1e-13 except where VGAM stops one iteration early (≤ 1e-7, see the
benchmark). Not ported: feather's mapping and filtering, the negative-binomial model, the blacklist file.

## `hiccups`

```python
hiccups(cool, *, resolutions=(5000, 10000, 25000), weights="weight", divisive=False, expected=None, chroms=None,
        fdr=0.1, peak_width=None, window=None, cluster_radius=None, thresholds=(0.02, 1.5, 1.75, 2.0),
        restrict=True, matrix_size=512, workers=1, verbose=True) -> HiCCUPS
```

Juicer's HiCCUPS (Rao et al. 2014; `juicer_tools hiccups --cpu`, 2.13.06 as in the ENCODE pipeline) on an
`.mcool`: expected by distance, the lower-left / donut / horizontal / vertical backgrounds, lambda-chunk Poisson
FDR thresholds, clustering and the 5 / 10 / 25 kb merge, with Juicer's defaults per resolution. Window sums run on
float64 summed-area tables; the pixels whose float32 sum could differ from Juicer's, and every reported pixel,
are re-summed in Juicer's cell order, so with Juicer's own SCALE vectors (`weights=`, `divisive=True`, from
`juicer_tools dump norm`) its four output files come out byte for byte. By default the cooler's balancing
weights are used. `workers=` forks processes over the search tiles; the result does not depend on it.

| argument | meaning |
|---|---|
| `cool` | an `.mcool` with the resolutions, or `{resolution: cooler URI}` |
| `weights` | a bins column, an array or `{resolution: ...}`; cooler convention, or Juicer's (`divisive=True`) |
| `expected` | `{resolution: {chrom: array}}` (Juicer's `dump expected`); default: computed as Juicer does |
| `restrict` | the CPU path's search limit (tiles up to ~8 Mb apart; at 25 kb only diagonal tiles) |

```python
r = gb.hiccups("sample.mcool", workers=16)
r.loops            # merged_loops: count, expected*, fdr*, numCollapsed, centroid1/2, radius, resolution
r.pixels[5000]     # postprocessed_pixels_5000
r.write("out")     # merged_loops.bedpe, postprocessed / enriched pixels, fdr_thresholds, Juicer's format
```

## `mustache`

```python
mustache(cool, *, resolution=None, chromosomes=None, distance=None, pt=0.2, st=0.88, sigma0=1.6,
         octaves=2, weight="weight", threads=4, verbose=True) -> Mustache
```
Mustache's multi-scale loop caller (Roayaei Ardakany et al. 2020) on a balanced cooler: per-diagonal local
z-scores, a difference-of-Gaussian scale space in overlapping 2000-bin blocks, an exponential P per scale,
BH within each block, Mustache's sparsity and diagonal-enrichment filters. It calls the loops Mustache 1.3.3
calls with its defaults (`mustache -f cool -r res -o out.tsv`), FDRs within 1e-11, 9–15× faster on 4 cores.
`weight` is the bin-table column of the balancing weights (Mustache's `-norm`), an array with one weight per
bin, or `None` for raw counts; `threads` scans blocks in parallel. The result's `loops` are bin pairs with
`count`, `balanced`, `P`, `FDR` and `scale` (`DETECTION_SCALE`); `r.write("out", "sample")` writes
`out/sample.tsv` in Mustache's format.

```python
r = gb.mustache("hic.mcool", resolution=5000)
# [mustache] 24 chromosomes, 394 blocks, 760 loops (FDR < 0.2) | read 4.2s, normalize 8.4s, scale space 121.6s, tables 0.1s
r.loops.to_pandas().sort_values("FDR").head()
r.write("out", "sample")
```
