---
title: atlas
parent: API Reference
layout: default
nav_order: 9
---

# `genomeblocks.atlas`
{: .no_toc }

GIGGLE-style enrichment over a sparse `bin x track` index, plus the two
modules that grow out of peak sets: `genomeblocks.se` (ROSE-style
super-enhancers) and `genomeblocks.hichip` (the ChIP-like short-range track
hidden in HiChIP pairs). See the [Atlas guide]({{ '/guide/atlas/' | relative_url }})
for the workflow.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Atlas`
{: .sec-navy }

```python
import genomeblocks as gb
from genomeblocks import Atlas
```

The genome is tiled at `bin_size` and one bit per (bin, track) cell is kept
in a single CSR matrix `M` of shape `(n_bins, n_tracks)`. A query becomes one
sparse matrix-vector product — `M[query_bins].sum(0)` — that returns per-track
shared-bin counts for thousands of tracks at once. Sized for
ChIP-Atlas-scale collections (~25k tracks: ~2.5 GB resident at 1 kb,
sub-second per query). Every query set is anything
[`as_loci`]({{ '/api/interop/' | relative_url }}) takes.

| Attribute | Meaning |
|---|---|
| `bin_size` | Tiling resolution (bp). |
| `n_bins` | Bins over all chromosomes. |
| `chrom_names`, `chrom_sizes`, `chrom_offsets` | Chromosome order, lengths and the first bin of each. |
| `track_names` | Track ids, in the order of `M`'s columns. |
| `track_n_peaks`, `track_n_bins` | Per-track interval count / covered-bin count. |
| `M` | `scipy.sparse.csr_matrix` `(n_bins, n_tracks)` of `uint8` ones. |
| `meta` | Optional `DataFrame` aligned to `track_names` (`None` until `attach_meta`). |

### `Atlas.make`

```python
Atlas.make(paths, *,
           chromsizes,                        # .chrom.sizes path, {chrom: length}, a Series, a Genome with sizes, or a cooler
           names: Sequence[str] | None = None,
           bin_size: int = 1000,
           workers: int | None = None,        # None = half the cores
           verbose: bool = True,
           meta=None,                         # TSV / CSV path or DataFrame -> attach_meta
           meta_columns: Sequence[str] | None = None,
           meta_id_col: str | None = None,
           meta_sep: str = "\t",
           name_pattern: str | None = None)   # re.search on the basename -> track id
    -> Atlas
```

`paths` is a glob, a directory (every `.bed`, `.narrowPeak`, `.broadPeak`,
plain or `.gz`), one path or a list of paths. Track ids default to the
basename without the BED-family extension; `name_pattern` is applied with
`re.search` and the first group (or the whole match) is kept — e.g.
`r'^[^.]+'` turns `SRX10895570.05.bed.gz` into `SRX10895570`. Files are read
in parallel over `workers` processes and assembled into one matrix. No file,
a missing file, or no interval on the given chromosomes raise with the
reason.

```python
CHROM_SIZES = {"chr1": 20_000, "chr2": 8_000}
atlas = Atlas.make(["t0.bed", "t1.bed", "t2.bed"], chromsizes=CHROM_SIZES, bin_size=500, workers=1)
atlas
# -> Atlas(tracks=3, bins=56, bin_size=500, nnz=21)
atlas.track_names, atlas.M.shape, atlas.chrom_offsets
# -> (['t0', 't1', 't2'], (56, 3), {'chr1': 0, 'chr2': 40})
Atlas.make("tracks/", chromsizes="genome.chrom.sizes", bin_size=500, verbose=False)      # a directory
Atlas.make("t*.bed", chromsizes=gb.Genome.from_sizes(CHROM_SIZES), bin_size=500, names=["x", "y", "z"], verbose=False)
Atlas.make("t0.bed", chromsizes=CHROM_SIZES, bin_size=500, name_pattern=r"^t(\d)", verbose=False).track_names
# -> ['0']
Atlas.make("none/*.bed", chromsizes=CHROM_SIZES)
# -> ValueError: No BED files matched: 'none/*.bed'
```

### `attach_meta`

```python
atlas.attach_meta(meta, *, id_col=None, columns=None, sep="\t") -> Atlas
```

Per-track metadata aligned to `track_names`: a TSV / CSV path or a
`DataFrame`. `columns` gives headers to a header-less file; `id_col` is the
column joined against `track_names` (default: the first). Tracks without a
row get `NaN`; rows for unknown tracks are dropped. The columns are appended
to every `search` / `bootstrap` result and to the track table. Returns
`self`.

```python
import pandas as pd
atlas.attach_meta(pd.DataFrame({"id": ["t0", "t1", "t9"], "factor": ["AR", "FOXA1", "X"]}), id_col="id")
atlas.meta
# ->    factor
#    id
#    t0     AR
#    t1  FOXA1
#    t2    NaN
atlas.attach_meta("meta.tsv", columns=["id", "factor"])            # header-less file
atlas.columns
# -> ['name', 'n_peaks', 'n_bins', 'factor']
```

### `search`

```python
atlas.search(query, *, ref=None, alternative="two-sided") -> pandas.DataFrame
```

Per-track Fisher 2x2 enrichment of the query bins, sorted by `giggle_score`
descending. With `ref=None` the other two cells come from the genome null
(`track bins - overlaps`, `genome bins - query bins - ...`); with
`ref=<intervals>` they come from the reference set ("more enriched in the
query than in the reference"). `p` is the exact hypergeometric probability;
`log2_odds` uses a 0.5 pseudocount; `giggle_score = -log10(p) * log2_odds`,
so depletion comes out negative. `alternative` is `'two-sided'`, `'greater'`
or `'less'`.

Columns: `name`, `n_query_bins`, `n_ref_bins` (with `ref` only),
`track_n_bins`, `track_n_peaks`, `overlaps`, `log2_odds`, `p`,
`giggle_score`, then the metadata columns. Also `Loci.enrich(atlas, **kw)`.
A query with no bins on the atlas's chromosomes raises `ValueError`.

```python
cre = gb.Loci.make("peaks.bed", keep=True)
atlas.search(cre)[["name", "overlaps", "giggle_score", "factor"]]
# ->   name  overlaps  giggle_score factor
#    0   t2         3      0.853650    NaN
#    1   t1         2      0.021147  FOXA1
#    2   t0         0     -1.116549     AR
atlas.search(cre.to_pandas(), ref="t1.bed", alternative="greater").columns.tolist()
# -> ['name', 'n_query_bins', 'n_ref_bins', 'track_n_bins', 'track_n_peaks', 'overlaps', 'log2_odds', 'p', 'giggle_score', 'factor']
cre.enrich(atlas).shape
# -> (3, 9)
atlas.search(gb.as_loci([("chrZ", 1, 2)]))
# -> ValueError: Query has no bins on the atlas's chromosomes.
```

### `bootstrap`

```python
atlas.bootstrap(query, *,
                n: int = 10,                 # iterations
                pool=None,                   # None: position-shuffle null; intervals: sample the null from this universe
                sample: int | None = None,   # per iteration, subsample every group AND the pool to this many regions
                replace: bool = False,
                keep_chrom: bool = True,     # shuffle within the original chromosome
                seed: int | None = None,
                verbose: bool = True) -> pandas.DataFrame
```

Observed overlaps against a resampled null. `query` is one interval set or a
`{group: intervals}` dict (a dict with `chrom` / `start` / `end` keys is one
set of columns, not groups). Three knobs interact:

- `pool=None` — the null shuffles each interval's position (within its
  chromosome when `keep_chrom`); each group shuffles independently.
- `pool=<intervals>` — the null is drawn from a curated universe (LOLA /
  regioneR style); use it when the query is a subset of the pool.
- `sample=<int>` — every group and the pool are subsampled to `sample`
  regions per iteration; the pool draw is shared across groups within an
  iteration, so groups of different size are paired against the same null.

Columns: `name`, `observed`, `expected`, `obs_std`, `null_std`, `log2fc`
(`log2((obs + 1) / (null + 1))`), `z` (paired z of obs - null over the
iterations), `p_emp` (`(1 + #(obs <= null)) / (n + 1)`), `track_n_bins`,
the metadata columns, and `group` first when `query` is a dict. Sorted by
`z` descending within each group. Also `Loci.enrich_mc(atlas, *, n=10, **kw)`.

```python
atlas.bootstrap(cre, n=20, seed=0, verbose=False)[["name", "observed", "expected", "log2fc", "z", "p_emp"]]
# ->   name  observed  expected    log2fc         z     p_emp
#    0   t2       3.0      1.85  0.489038  0.807351  0.333333
#    1   t1       2.0      1.70  0.152003  0.290953  0.619048
#    2   t0       0.0      1.55 -1.350497 -1.476102  1.000000
groups = atlas.bootstrap({"a": cre.head(4), "b": "t1.bed"}, n=5, seed=0, verbose=False)
groups.columns.tolist()
# -> ['group', 'name', 'observed', 'expected', 'obs_std', 'null_std', 'log2fc', 'z', 'p_emp', 'track_n_bins', 'factor']
atlas.bootstrap(cre, pool=cre.to_polars(), n=5, sample=3, seed=1, verbose=False).shape
# -> (3, 10)
```

### `save` and `load`

```python
atlas.save(path) -> None        # a compressed .npz, no pickles
Atlas.load(path) -> Atlas
```

Everything is stored, including the metadata (as strings; a missing value
comes back as `NaN`, as `attach_meta` stores it), so `save` → `load` is an
exact round trip.

```python
atlas.save("atlas.npz")
back = Atlas.load("atlas.npz")
back.track_names == atlas.track_names, (back.M != atlas.M).nnz == 0, back.meta.loc["t0", "factor"]
# -> (True, True, 'AR')
```

### The track table

As a table an `Atlas` is its track table: `name`, `n_peaks`, `n_bins` and
the metadata columns.

| Member | One line |
|---|---|
| `len(atlas)`, `iter(atlas)` | Number of tracks; iterate the track names. |
| `atlas.shape`, `atlas.columns` | Like a DataFrame. |
| `atlas["SRX..."]`, `atlas[i]` | That track's row as a dict. |
| `atlas.to_pandas()`, `to_polars()`, `to_arrow()` | The track table. |
| `atlas.head(n=5)`, `tail(n=5)` | First / last tracks. |
| `atlas.describe()` (alias `summary()`) | Tracks, chromosomes, bin size, bins, nonzero cells, density, bins per track, metadata columns. |
| `atlas._repr_html_()` | Notebook summary. |

The Arrow C stream, `__dataframe__` and narwhals protocols are implemented,
so `pl.DataFrame(atlas)` and `duckdb.sql(...)` take it directly.

```python
atlas.shape, atlas.columns, list(atlas)
# -> ((3, 4), ['name', 'n_peaks', 'n_bins', 'factor'], ['t0', 't1', 't2'])
atlas["t1"]
# -> {'name': 't1', 'n_peaks': 7, 'n_bins': 7, 'factor': 'FOXA1'}
atlas.describe()
# ->                        value
#    tracks                     3
#    chromosomes                2
#    bin size                 500
#    bins                      56
#    nonzero cells             21
#    density                0.125
#    bins per track min         7
#    bins per track median    7.0
#    bins per track max         7
#    metadata columns      factor
```

---

## `genomeblocks.se`
{: .sec-navy }

Super-enhancers, ROSE-style, on columnar tables: stitch peaks within
`stitch` bp into regions, score each region as mean signal x width, and cut
the ranked scores at the slope-1 knee.

```python
from genomeblocks import se
```

| Function | Signature | Returns |
|---|---|---|
| `knee` | `knee(values)` | `(cutoff, order)`: `order` indexes the positive values sorted descending; the first `cutoff` of them are above the knee. |
| `stitch_peaks` | `stitch_peaks(peaks, stitch=12_500, *, backend=None)` | A `Loci` of stitched regions in genome order with an `n_peaks` column. |
| `call_se` | `call_se(peaks, bigwigs, *, stitch=12_500, workers=1, verbose=False, return_all=False, backend=None, **signal_kw)` | The SEs as a `Loci` with `n_peaks`, `score`, `rank`; with `return_all=True`, `(se, all_regions)`. |
| `nearest_gene_within` | `nearest_gene_within(regions, tss, maxd)` | Per region, the row of the closest TSS, or -1 when it is farther than `maxd`. |

### `knee`

The same slope-1 rule as `Architecture.elbow`, so super-enhancers and prime
hubs are cut the same way. The curve is normalised to the unit square and
smoothed for long inputs; the knee is where its gradient first reaches 1.

```python
se.knee([10, 9, 8, 1, 0.5, 0.2, 0])
# -> (4, array([0, 1, 2, 3, 4, 5]))
```

### `stitch_peaks`

Peaks whose gap is at most `stitch` bp form one region running from the
first peak's start to the last peak's end. `peaks` is anything `as_loci`
takes; `backend` picks the interval engine.

```python
peaks = gb.as_loci([("chr1", 1000, 1500), ("chr1", 1600, 2000), ("chr1", 2500, 3000), ("chr1", 10_000, 11_000),
                    ("chr1", 15_000, 15_500), ("chr2", 100, 200)])
st = se.stitch_peaks(peaks, stitch=1000)
st.to_records(), st.cols["n_peaks"].tolist()
# -> ([('chr1', 1000, 3000, '.'), ('chr1', 10000, 11000, '.'), ('chr1', 15000, 15500, '.'), ('chr2', 100, 200, '.')], [3, 1, 1, 1])
se.stitch_peaks(peaks.to_pandas()).to_records()        # the ROSE default, 12.5 kb
# -> [('chr1', 1000, 15500, '.'), ('chr2', 100, 200, '.')]
```

### `call_se`

`bigwigs` is one or several bigWig paths / open handles / a `{name: path}`
dict — what [`Loci.signal`]({{ '/api/signal/' | relative_url }}) takes;
the score averages the tracks. Extra keywords (e.g. `exact=False`) go to
`signal`; `backend` picks the bigWig engine. Regions with a zero score get
rank 0. Also `Loci.call_se(bigwigs, **kw)`.

```python
ses, every = se.call_se(peaks, "signal.bw", stitch=1000, return_all=True, verbose=True)
# -> [INFO] 6 peaks → 4 stitched regions → 2 SEs (median 1,500 bp)
every.to_pandas()
# ->   chrom  start    end strand  n_peaks         score  rank
#    0  chr1   1000   3000      .        3  11472.912788     1
#    1  chr1  10000  11000      .        1   5850.211143     2
#    2  chr1  15000  15500      .        1   1742.084980     3
#    3  chr2    100    200      .        1      0.000000     0
ses
# -> Loci(n=2, chroms=1, cols=[n_peaks, score, rank])
peaks.call_se(["signal.bw", "signal2.bw"], stitch=1000, backend="python")
```

### `nearest_gene_within`

`tss` is a 1-bp TSS `Loci`, e.g. `genes.get_tss()`; the distance is 0 when
the TSS lies inside the region, else the distance to the nearest edge. The
two tables may carry different `Genome`s.

```python
genes = gb.Genes.make("genes.gtf")
regions = gb.as_loci([("chr3", 5, 10), ("chr1", 0, 2000), ("chr1", 10_500, 10_600)])
se.nearest_gene_within(regions, genes.get_tss(), 1000).tolist()
# -> [-1, 0, 1]
genes.genes.cols["gene_name"].tolist()
# -> ['GENE_A', 'GENE_B', 'GENE_C']
```

---

## `genomeblocks.hichip`
{: .sec-navy }

Ligation pairs closer than ~1 kb are mostly undigested ChIP fragments, so
their 5' read ends behave like ChIP-seq reads. This module does, in one
process, what the usual shell recipe does with awk, sort, bedtools and a
bedGraph converter: cis pairs within `max_dist`, both 5' ends (stranded), a
BED6 for MACS3, ends extended to fragments, coverage, bigWig. Peak calling
stays with MACS3.

```python
from genomeblocks import hichip
```

| Function | Signature | Returns |
|---|---|---|
| `shortrange_ends` | `shortrange_ends(pairs, max_dist=1000, *, columns=None, genome=None)` | Stranded 1-bp `Loci` of both 5' ends of every cis pair with `\|pos2 - pos1\| <= max_dist`, genome-sorted. |
| `write_bed` | `write_bed(ends, path)` | BED6 of the ends (the MACS3 input); returns `path`. |
| `fragments` | `fragments(ends, extsize=147, chrom_sizes=None)` | Each end extended `extsize` bp in its read direction, clipped to the chromosome when `chrom_sizes` is given. |
| `coverage` | `coverage(L)` | Generator of `(chrom, starts, ends, depth)` runs with depth > 0, adjacent equal depths merged (the records of `bedtools genomecov -bg`). |
| `to_bigwig` | `to_bigwig(L, path, chrom_sizes)` | Coverage of `L` as a bigWig; returns `path`. |
| `to_bedgraph` | `to_bedgraph(L, path)` | Coverage as bedGraph; returns `path`. |
| `macs3` | `macs3(bed, name, outdir, *, gsize="hs", extsize=147, q=0.01, exe="macs3")` | Runs `macs3 callpeak --nomodel --extsize ... --keep-dup all`; returns the narrowPeak path. |
| `shortrange_track` | `shortrange_track(pairs, out_prefix, chrom_sizes, *, max_dist=1000, extsize=147, genome=None)` | The whole recipe minus peak calling: `{'ends': Loci, 'bed': path, 'bigwig': path}`. |

`pairs` is a HiC-Pro `allValidPairs` file (plain or `.gz`); `columns` maps
`chr1` / `pos1` / `strand1` / `chr2` / `pos2` / `strand2` to 0-based column
numbers for other layouts (4DN `.pairs`:
`{"chr1": 1, "pos1": 2, "chr2": 3, "pos2": 4, "strand1": 5, "strand2": 6}`).
The file is scanned lazily with polars when installed, else in pandas chunks
— the same rows either way. `chrom_sizes` is a dict, a `.chrom.sizes` path or
a `Genome` with sizes. Every other input is anything `as_loci` takes.

```python
ends = hichip.shortrange_ends("hichip.allValidPairs", 1000)
ends
# -> Loci(n=40, chroms=1, sorted)
ends.head(3).to_pandas()
# ->   chrom  start  end strand
#    0  chr1     99  100      +
#    1  chr1    109  110      +
#    2  chr1    119  120      +
hichip.write_bed(ends, "ends.bed")
# -> 'ends.bed'
frags = hichip.fragments(ends, 147, {"chr1": 20_000, "chr2": 8_000})
frags.lengths.max()
# -> 147
chrom, starts, stops, depth = next(hichip.coverage(frags))
chrom, starts[:3], stops[:3], depth[:3]
# -> ('chr1', array([ 99, 109, 119]), array([109, 119, 129]), array([1, 2, 3]))
hichip.to_bigwig(frags, "cov.bw", "genome.chrom.sizes")
# -> 'cov.bw'
hichip.to_bedgraph(frags, "cov.bedGraph")
out = hichip.shortrange_track("hichip.allValidPairs", "sample", {"chr1": 20_000, "chr2": 8_000})
out["bed"], out["bigwig"]
# -> ('sample_shortrange_ends.bed', 'sample_shortrange.bw')
hichip.macs3(out["bed"], "sample", "macs_out")          # on a machine without MACS3
# -> RuntimeError: 'macs3' not found on PATH: pip install macs3 (or conda install -c bioconda macs3), or pass exe=/path/to/macs3
```
