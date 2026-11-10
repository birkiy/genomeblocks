---
title: loops
parent: API Reference
layout: default
nav_order: 7.5
---

# `genomeblocks.loops`
{: .no_toc }

HiChIP loop calling with FitHiChIP's model — peak-to-all bin pairs, a
coverage-bias regression, a binomial test and the merging of adjacent loops —
on the pixels of a cooler, in numpy. It returns FitHiChIP 11.0's tested pairs,
loops and merged loops (`IntType=3`, `BiasType=1`, `UseP2PBackgrnd=0`,
`MergeInt=1`, its defaults otherwise) in seconds to minutes instead of hours.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
from genomeblocks.loops import fithichip, FitHiChIP
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
