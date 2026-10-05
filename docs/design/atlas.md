---
title: Atlas
parent: Design
layout: default
nav_order: 6
---

# Atlas
{: .no_toc }

`Atlas` answers "which of these thousands of ChIP-seq tracks overlap my regions
more than expected?" with one sparse matrix and one vectorised test. The
query can be any interval set `as_loci` takes.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## One bit per bin and track
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/atlas-index.svg %}
</div><figcaption>
<strong>A query is a column sum.</strong> The genome is tiled into fixed bins, and every track becomes one column of a sparse bins × tracks matrix with a 1 wherever the track has a peak. A query's intervals are turned into the set of bins they touch, and <code>M[query_bins].sum(axis=0)</code> gives the overlap count for every track in one sparse operation. That fills the 2×2 table of a Fisher test for all tracks at once.
</figcaption></figure>

| Attribute | What it holds |
|---|---|
| `M` | `scipy.sparse.csr_matrix`, shape `(n_bins, n_tracks)`, `uint8` ones |
| `bin_size`, `chrom_names`, `chrom_sizes`, `chrom_offsets`, `n_bins` | the tiling: chromosomes laid end to end, `chrom_offsets[c]` = first bin of `c` |
| `track_names`, `track_n_peaks`, `track_n_bins` | one entry per column |
| `meta` | optional DataFrame aligned to `track_names` (cell type, antibody, ...) |

```python
import genomeblocks as gb

atlas = gb.Atlas.make("atlas/", chromsizes="genome.chrom.sizes", bin_size=500, verbose=False)
atlas
# -> Atlas(tracks=3, bins=56, bin_size=500, nnz=21)
type(atlas.M).__name__, atlas.M.shape, atlas.M.dtype, atlas.M.nnz
# -> ('csr_matrix', (56, 3), dtype('uint8'), 21)
atlas.chrom_offsets, atlas.n_bins
# -> ({'chr1': 0, 'chr2': 40}, 56)
```

## Building the index
{: .sec-navy }

`Atlas.make(paths, chromsizes=, bin_size=1000, workers=None)` takes a
directory, a glob, one path or a list of BED-like files, and chromosome sizes
as a `.chrom.sizes` path, a dict, a Series, a `Genome` with sizes or a cooler
(`genome.read_sizes`). It tiles each chromosome at `bin_size`, lays the
chromosomes end to end with fixed offsets, and then:

1. **reads each file into bin ids** in a worker process (`workers=None` is
   half the cores): pandas reads the three coordinate columns, rows on
   unknown chromosomes are dropped, each interval becomes the bin range
   `[start // bin, (end − 1) // bin + 1)` clipped to the chromosome, and the
   ranges are exploded into one sorted, unique bin-id array per track;
2. **assembles one CSC matrix** in a single pass (tracks are columns, the bin
   ids are the row indices, every value is 1) and stores it as CSR, which
   makes `M[query_bins]` a row slice;
3. **attaches metadata** when `meta=` is given (`attach_meta` joins a TSV /
   CSV or a DataFrame on the track name; unknown tracks get NaN).

`save()` / `load()` round-trip the whole index, including the metadata,
through one compressed `.npz`.

## Searching
{: .sec-navy }

Every query goes through `_intervals_to_bin_ranges`: `as_loci` on the input,
then the same bin arithmetic as the build, in numpy over all rows at once
(intervals on chromosomes the atlas lacks are dropped). The union of those
ranges is the query's bin set, and `_overlap_counts` is the column sum.

```python
cre = gb.Loci.make("peaks.bed")
atlas._intervals_to_bin_ranges(cre).tolist()
# -> [[1, 3], [3, 5], [9, 11], [19, 21], [21, 23], [41, 42], [50, 51]]
```

| Mode | Cells of the 2×2 table |
|---|---|
| `search(query)` | a = query bins in the track, b = rest of the query, c = rest of the track, d = rest of the genome |
| `search(query, ref=other)` | a / b from the query, c / d from the reference set: "more in query than in ref" |
| `bootstrap(query, n)` | an empirical null from `n` position shuffles of the query (each interval stays on its chromosome by default) |
| `bootstrap(query, n, pool=universe)` | the null is sampled from a curated set of regions instead of the genome |
| `bootstrap({g1: ..., g2: ...}, n, sample=k)` | several groups, each subsampled to `k` regions per iteration against one shared pool draw, so groups are compared on a paired null |

`_fisher_vec` runs the exact hypergeometric test on all tracks as arrays and
returns `log2_odds` (with 0.5 added to each cell), `p` and the GIGGLE-style
score `−log10(p) · log2(OR)`, whose sign flips for depletion. The result is a
DataFrame per track, joined with the metadata, sorted by score.

```python
atlas.search(cre)
# ->   name  n_query_bins  track_n_bins  track_n_peaks  overlaps  log2_odds         p  giggle_score
# -> 0   t2            12             7              7         3   1.729352  0.320905      0.853650
# -> 1   t1            12             7              7         2   0.773960  0.939026      0.021147
# -> 2   t0            12             7              7         0  -2.321928  0.330467     -1.116549

atlas.bootstrap(cre, n=5, seed=0, verbose=False)[["name", "observed", "expected", "z", "p_emp"]]
# ->   name  observed  expected         z     p_emp
# -> 0   t2       3.0       2.4  0.395628  0.666667
# -> 1   t1       2.0       1.6  0.350823  0.666667
# -> 2   t0       0.0       2.2 -2.629503  1.000000
```

`Loci.enrich(atlas)` and `Loci.enrich_mc(atlas, n=)` call `search` and
`bootstrap` on a `Loci` directly, and a frame, a BED path or a dict of
columns works as the query just the same (a dict of groups is told apart
from a dict of columns by its keys).

Working in bins rather than base pairs is what makes a query one sparse
product. The cost is resolution: two peaks in the same bin count as one bin.
Choose `bin_size` near your peak width.

## The track table
{: .sec-navy }

As a table, an `Atlas` is its track table: `len`, `shape`, `columns`,
`head()`, `describe()`, `to_pandas()` (name, `n_peaks`, `n_bins`, then the
metadata columns), `atlas["SRX..."]`, and the Arrow / dataframe protocols from
`TableMixin`, so the catalogue of a 25,000-track atlas can be filtered in
polars or DuckDB before any query.

```python
import pandas as pd, polars as pl
atlas.attach_meta(pd.DataFrame({"id": ["t0", "t1"], "factor": ["AR", "FOXA1"]}), id_col="id")
atlas.columns
# -> ['name', 'n_peaks', 'n_bins', 'factor']
pl.DataFrame(atlas).shape
# -> (3, 4)
atlas["t1"]
# -> {'name': 't1', 'n_peaks': 7, 'n_bins': 7, 'factor': 'FOXA1'}
atlas.save("atlas.npz"); gb.Atlas.load("atlas.npz").meta.loc["t0", "factor"]
# -> 'AR'
```

## Costs
{: .sec-navy }

- **Build:** one BED read per track in parallel, then one `csc_matrix`
  construction; peak memory is about four bytes per set bit. ChIP-Atlas scale
  (~25,000 mouse tracks at 1 kb) is a few GB resident.
- **Query:** one `as_loci`, one vectorised bin-range computation, one
  `unique`, one sparse row slice and column sum, one vectorised Fisher test:
  sub-second whatever the number of tracks.
- **Bootstrap:** `n` repetitions of the query step on shuffled or sampled bin
  sets; the shuffle itself is a few numpy draws per iteration.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
