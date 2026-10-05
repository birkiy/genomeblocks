---
title: Signal
parent: User Guide
layout: default
nav_order: 5
---

# Signal
{: .no_toc }

bigWig signal over a set of loci as one `(rows × tracks × bins)` cube, TMM
normalization, and heatmaps and profiles of that cube. Reading goes through
the bigwig backend: [pybigtools](https://github.com/jackh726/bigtools) by
default (see [Credits]({{ '/credits/#pybigtools' | relative_url }})), pyBigWig
or the pure-Python reader on request.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `loci.signal()` — pull a `(rows × tracks × bins)` cube
{: .sec-green }

```python
import genomeblocks as gb

cre = gb.Loci.make("cre.bed")
S = cre.signal(["atac.bw", "h3k27ac.bw"], n_bins=200, flank=3_000)
# -> [INFO] Extracting 2 bigwigs for 7 loci into 200 bins (span=False, agg='mean', backend='pybigtools', exact=True, workers=1).
S.shape, S.dtype
# -> ((7, 2, 200), dtype('float32'))
```

Row *i* of the cube is row *i* of the loci, so it lines up with every other
table built from them (vertex columns of an
[Architecture]({{ '/guide/architecture/' | relative_url }}), motif matrices,
annotations). The result is a plain `numpy.ndarray`.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-cube.svg %}
</div><figcaption>
<strong>One row of the cube per locus, one grid per track.</strong> Each row becomes the window centre ± <code>flank</code> (or the interval itself with <code>span=True</code>); the bigwig backend returns its <code>n_bins</code> summaries in one call, and they land in <code>cube[i, t, :]</code>. Row <em>i</em> is still locus <em>i</em>, so a mask over the loci selects the same rows of the cube, and the tracks keep the order you passed.
</figcaption></figure>

| arg | default | meaning |
|---|---|---|
| `bigwigs` | | one bigWig or several: paths, open pyBigWig / pybigtools handles, or a `{name: path}` dict (tracks in its order) |
| `n_bins` | 200 | bins per row |
| `flank` | 3000 | bp on each side of the row's centre (ignored with `span=True`) |
| `span` | `False` | bin the whole interval instead of centre ± flank |
| `agg` | `'mean'` | per-bin statistic: `mean`, `min`, `max`, `sum`, `std`, `coverage` |
| `dtype` | `float32` | cube dtype |
| `workers` | 1 | 1 = sequential; > 1 = that many processes; `None` = half the cores |
| `backend` | `None` | `'pybigtools'`, `'pybigwig'` or `'python'`; `None` = the default |
| `exact` | `True` | base-accurate binning; `False` lets pybigtools / pyBigWig use zoom levels (approximate, ~3× faster) |
| `progress`, `verbose` | `True` | progress bar and the `[INFO]` line |

The function form takes anything [`as_loci`]({{ '/interoperability/' | relative_url }})
reads — a frame, a BED path, a list of regions:

```python
from genomeblocks.signal import signal
S = signal(cre.to_pandas(), "atac.bw", n_bins=200, flank=3_000)
S = signal("cre.bed", {"ATAC": "atac.bw", "H3K27ac": "h3k27ac.bw"}, n_bins=100, flank=1_000)
```

### Inputs: paths, handles, dicts

Paths and open handles mix freely in one call. An open handle is read with its
own library, so `backend=` cannot apply to it — pass the path instead:

```python
import pyBigWig
h = pyBigWig.open("atac.bw")
S = cre.signal([h, "h3k27ac.bw"], n_bins=50, flank=1_000)        # fine
S = cre.signal([h], n_bins=50, flank=1_000, backend="python", verbose=False)
# -> ValueError: the track is an open pybigwig handle, so backend='python' cannot apply: pass the file path instead, or drop backend=
```

### Windows, chromosome edges, missing chromosomes

With `span=False` the window of a row is `centre - flank` to `centre + flank`,
where `centre = (start + end) // 2`. With `span=True` it is the row itself.
Windows are not clipped: a window that leaves the chromosome is binned over
the full grid, bases outside have no data, and the bins they fall in read 0.
A row on a chromosome the file lacks stays 0 for that track.

```python
edges = gb.as_loci([("chr1", 20, 80), ("chr1", 19_950, 20_000), ("chr3", 10, 20)])
edges.signal("atac.bw", n_bins=5, flank=100, verbose=False)[:, 0].round(2)
# -> array([[0.  , 6.73, 6.73, 3.43, 3.43],      # window starts before base 0
# ->        [0.  , 0.  , 0.  , 0.  , 0.  ],      # no data in the last 100 bp of chr1
# ->        [0.  , 0.  , 0.  , 0.  , 0.  ]], dtype=float32)     # chr3 is not in the file
```

### Aggregators

`mean`, `min`, `max`, `sum`, `std` (population) and `coverage` (fraction of
the bin's bases with data). A bin without data is 0. Anything else raises:

```python
cre.signal("atac.bw", n_bins=4, flank=100, agg="median", verbose=False)
# -> ValueError: unknown stat 'median'; use mean, min, max, sum, std or coverage
```

---

## bigWig backends and the one binning rule
{: .sec-green }

| Backend | Library | Notes |
|---|---|---|
| `pybigtools` (default) | Rust, `pip install pybigtools` | one native call per region does I/O, decompression and binning |
| `pybigwig` | C (libBigWig), `pip install pyBigWig` | |
| `python` | pure Python in `genomeblocks.backends._bbi` (numpy + mmap) | no compiled dependency; always available |

All three give the same numbers. Bin `b` of a window of `n` bases covers
`[floor(n·b/n_bins), floor(n·(b+1)/n_bins))` for every engine — the integer
edges libBigWig and pybigtools use — so fractional bin widths bin alike, and
bins outside the chromosome are 0.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-bins-rule.svg %}
</div><figcaption>
<strong>One binning rule for every engine.</strong> A 10-base window with <code>n_bins = 4</code> has edges <code>floor(10·b/4) = 0, 2, 5, 7, 10</code>, so the bins are 2, 3, 2 and 3 bases wide and the per-bin means are 2, 3.33, 4 and 4.67; pybigtools, pyBigWig and the pure-Python reader place the edges identically and report the same numbers. A window that hangs off the start of the chromosome keeps the same edge rule; its first bin lies entirely before base 0, has no data, and reads 0 (the <code>missing</code> value).
</figcaption></figure>

```python
for name in ("pybigtools", "pybigwig", "python"):
    Sb = cre.signal("atac.bw", n_bins=20, flank=500, backend=name, verbose=False, progress=False)

with gb.use_backend(bigwig="python"):             # for a block of code
    Sb = cre.signal("atac.bw", n_bins=20, flank=500)
```

A requested backend that is not installed raises an `ImportError` with the
install command; nothing is swapped in silently. The same handles serve the
[browser]({{ '/guide/browser/' | relative_url }}) tracks, so `backend=` means
the same thing there.

---

## Execution model and sizing
{: .sec-green }

- **`workers=1` (default) is sequential**: one native call per region and
  track, which is the fastest choice for typical heatmap-sized jobs.
- **`workers > 1` uses processes, not threads.** Children write into one
  shared-memory cube; work is split by (track chunk, loci chunk) and each child
  opens its own handles. This needs **paths** — open handles cannot cross a
  process boundary:

```python
tiles = gb.Loci.tile_genome("genome.chrom.sizes", 500)
S_tiles = tiles.signal(["atac.bw", "h3k27ac.bw"], n_bins=8, flank=400, workers=2)   # same cube as workers=1
cre.signal([pyBigWig.open("atac.bw"), pyBigWig.open("h3k27ac.bw")], n_bins=8, workers=2, verbose=False)
# -> ValueError: workers > 1 needs bigWig paths (open handles cannot be shared between processes); pass the file paths, or workers=1.
```

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-parallel.svg %}
</div><figcaption>
<strong>Workers fill one shared cube in place.</strong> With <code>workers &gt; 1</code> the cube lives in shared memory; each process opens its own handles for its slab of tracks and loci and writes straight into it, so nothing is copied back through pipes. That is why the call needs file paths rather than open handles, and why the cube is identical to the one <code>workers=1</code> produces.
</figcaption></figure>

- An explicit `workers` is capped at `min(workers, n_tracks · ⌈n_loci / 1000⌉, cpu_count)`,
  so a small job (one track, under 1000 loci) runs sequentially even with
  `workers=2`, handles and all. `workers=None` asks for half the cores.
- The cube takes `n_loci × n_tracks × n_bins × itemsize` bytes. A cube over
  half of the available RAM is refused before anything is read:

```python
gb.Loci.tile_genome({"chr1": 10**9}, 10).signal("atac.bw", n_bins=200)
# -> MemoryError: The signal cube needs ~80.0 GB, over half of the available RAM; extract in chunks of loci (L[a:b].signal(...)).
```

Chunk by rows (`cre[a:b].signal(...)`) and stream to disk when you need more.

---

## The cube and its converters
{: .sec-green }

The cube is a numpy array; `genomeblocks.interop` labels it with the loci:

```python
from genomeblocks.interop import cube_to_xarray, cube_to_anndata, cube_to_pandas

S = cre.signal(["atac.bw", "h3k27ac.bw"], n_bins=20, flank=500, verbose=False)
da = cube_to_xarray(S, cre, ["ATAC", "H3K27ac"], flank=500)
da.dims, da.shape
# -> (('region', 'track', 'bin'), (7, 2, 20))
da.coords["bin"].values[:3], da.coords["region"].values[:2]
# -> (array([-475., -425., -375.]), array(['chr1:900-1100', 'chr1:1900-2100'], dtype=object))

ad = cube_to_anndata(S, cre, ["ATAC", "H3K27ac"])        # obs = tracks, var = loci
ad
# -> AnnData object with n_obs × n_vars = 2 × 7
# ->     var: 'chrom', 'start', 'end', 'strand'
# ->     varm: 'bins:ATAC', 'bins:H3K27ac'

df = cube_to_pandas(S, cre, ["ATAC", "H3K27ac"])         # (track, bin) columns, indexed by region name
df.shape
# -> (7, 40)
```

`cube_to_anndata` puts the per-row aggregate over bins (`agg='mean'`, `'sum'`
or `'max'`) in `X` and keeps each track's full profile in `varm['bins:<track>']`.
With a single bin, `cube_to_pandas` gives one column per track — the natural
"signal per region" table:

```python
cube_to_pandas(cre.signal("atac.bw", n_bins=1, span=True, verbose=False), cre, ["ATAC"])
# ->                       ATAC
# -> region
# -> chr1:900-1100     8.879168
# -> chr1:1900-2100    3.029531
# -> chr1:4900-5100    8.288841
# -> chr1:9950-10050   5.319891
# -> chr1:10900-11100  7.058707
# -> chr2:500-600      0.000000
# -> chr2:5000-5100    0.000000
```

---

## TMM normalization
{: .sec-green }

```python
N = gb.tmm(S)
```

Per-track TMM factors (Robinson & Oshlack 2010, edgeR's rule) computed on the
per-row means, then library-size scaling to counts per million, so tracks become
comparable. The algorithm is vendored in `genomeblocks.signal`; there is no
external normalization dependency.

A track with no signal in any row is left at 0 with a warning; a cube with no
signal at all raises:

```python
gb.tmm(S_with_an_empty_track)
# -> RuntimeWarning: tmm: track(s) [1] have no signal and are left at 0
gb.tmm(np.zeros((5, 2, 3)))
# -> ValueError: tmm: every track is empty (no signal in any region)
```

---

## `plot_heatmap`
{: .sec-green }

A deeptools-style grid: one column per track, one block of rows per group,
and the mean profile of each group above. `cre.plot_heatmap(S, ...)` and
`genomeblocks.signal_draw.plot_heatmap(cre, S, ...)` are the same function.

```python
fig = cre.plot_heatmap(S, samples=["ATAC", "H3K27ac"], vmax=[10, 6], height=3_000)
fig.savefig("heatmap.pdf")
```

### Groups as Loci, masks or rows

`groups` is a `{name: rows}` dict. Each value may be a boolean mask over the
plotted loci, an array of row numbers, or an interval set — a Loci or anything
`as_loci` reads, matched to the plotted loci by chromosome, start and end:

```python
genes  = gb.Genes.make("genes.gtf")
prom   = genes.get_tss().slop(1_000)                     # TSS ± 1 kb
groups = {
    "promoter": cre.intersect(prom),                      # a Loci
    "distal":   cre.difference(prom).to_pandas(),         # a frame
    "top 3":    np.arange(3),                             # row numbers
    "chr2":     cre.chroms == "chr2",                     # a mask
}
fig = cre.plot_heatmap(S, groups=groups, sets=["promoter", "distal"],   # sets: block order
                       samples=["ATAC", "H3K27ac"], height=3_000)
```

A group name in `sets` that is not in `groups` raises. Omit `groups` to plot
all rows as one block.

### Scales, colours, order

| arg | meaning |
|---|---|
| `samples` | column titles, one per track (a wrong count raises) |
| `vmax` | colour scale top: a scalar, or one value per track |
| `ymax`, `ymin` | profile y-limits: scalar or one per track |
| `cmap` | colormap: a name, or one per track (`["Blues", "Reds"]`) |
| `colors` | `{group: colour}` for the profile lines (default tab10) |
| `height` | bp from the centre to the window edge, for the x tick labels (`-3.0kb`, `center`, `+3.0kb`) |
| `profile` | draw the mean profiles above the heatmaps (default `True`) |
| `sort` | `'group'` (by mean signal within each block, default), `'global'`, or `None` (input order) |
| `ylabel`, `dpi` | profile axis label (default `'CPM signal'`), figure dpi |

```python
fig = cre.plot_heatmap(S, groups=groups, samples=["ATAC", "H3K27ac"],
                       cmap=["Blues", "Reds"], vmax=[8, 8], ymax=[8, 8], sort="global")
```

---

## `plot_profiles`
{: .sec-green }

The mean profile of one track per group, side by side:

```python
fig = cre.plot_profiles(S, groups=groups, track=1, ylim=5, height=3_000)
```

`track` picks the cube column (default 0); `ylim` fixes the y-axis (default:
the 99th percentile of the profile).

---

## `compare_heatmap` — A-only / shared / B-only
{: .sec-green }

Compare two interval sets across several marks in one figure:

```python
fig, union, S, groups = gb.compare_heatmap(
    a=cre_mesc, b=cre_hesc,                      # anything as_loci reads
    bigwigs=["ATAC_mESC.bw", "ATAC_hESC.bw", "H3K27ac_mESC.bw", "H3K27ac_hESC.bw"],
    a_name="mESC", b_name="hESC", common_name="shared",
    n_bins=200, flank=3_000,
    normalize=True,                               # tmm() before plotting
    samples={"ATAC": [0, 1], "H3K27ac": [2, 3]},  # average tracks 0 and 1 into one column
    cmap=["Blues", "Reds"], vmax=[10, 6],         # one per plotted column (after merging)
)
fig.savefig("compare.pdf")
```

The cube is extracted for the **union** of the two sets in the order a-only,
shared, b-only, and `groups` maps each block name to a boolean mask over
`union`. `samples` may be a list of titles (no merging) or a
`{column: [track indices]}` dict that averages replicate tracks into one
column; per-column `cmap`, `vmax` and `ymax` lists then have one entry per
merged column. `sets` controls block order (default `[a_name, common_name, b_name]`);
`sort`, `cmap`, `vmax`, `ymax`, `colors`, `profile` are passed through to
`plot_heatmap`. Extra extraction arguments go in `signal_kw`
(`signal_kw={"backend": "python", "workers": 4}`).

The function returns `(fig, union, S, groups)`. Pass the cube back as `S=` to
re-plot with other parameters without reading the files again:

```python
fig, union, S, groups = gb.compare_heatmap(a, b, bigwigs, n_bins=200, flank=3_000)
fig2, *_ = gb.compare_heatmap(a, b, [], S=S, flank=3_000, vmax=4, sort=None)
```

---

## Motif cubes — `plot_motif_heatmap`
{: .sec-green }

`scan_motifs_profile` (see [Motifs]({{ '/guide/motifs/' | relative_url }}))
returns a `(rows × motifs × bins)` cube of hit positions. `gb.plot_motif_heatmap`
draws it with the same machinery, one column per motif:

```python
M, names = cre.scan_motifs_profile("genome.fa", "motifs.jaspar", r=500, n_bins=20)
fig = gb.plot_motif_heatmap(cre, M, names, r=500, groups=groups)
```

`r` is the scan half-window (x tick labels). `vmax` and `ymax` default to a
per-motif scale — the 98th percentile of the non-zero cells and the top of the
group mean profiles — and the default colormap is `'Purples'` with the y label
`'motif hits/bin'`. Any other `plot_heatmap` argument passes through.

---

## Example: histone marks around CRE centres
{: .sec-green }

```python
import genomeblocks as gb

cre = gb.Loci.make("cre.bed")
marks = {"H3K4me1": "H3K4me1.bw", "H3K4me3": "H3K4me3.bw",
         "H3K27ac": "H3K27ac.bw", "H3K27me3": "H3K27me3.bw"}

S = gb.tmm(cre.signal(marks, n_bins=200, flank=3_000))
genes = gb.Genes.make("genes.gtf")
prom = genes.get_tss().slop(1_000)

fig = cre.plot_heatmap(S, groups={"promoter": cre.intersect(prom), "distal": cre.difference(prom)},
                       samples=list(marks), vmax=6, height=3_000)
fig.savefig("marks.pdf")
```
