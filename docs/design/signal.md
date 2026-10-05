---
title: Signal
parent: Design
layout: default
nav_order: 3
---

# Signal
{: .no_toc }

`signal()` reads bigWig tracks into one dense array, the *signal cube*, with
one row per interval, one plane per track and one column per bin. Heatmaps,
profiles, TMM normalisation and the converters to xarray / AnnData / pandas
all read that array. The files are read through the `bigwig` backend.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## One row per interval, one call per window
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-cube.svg %}
</div><figcaption>
<strong>The cube is filled row by row.</strong> Each interval becomes a window around its centre (± <code>flank</code>), or its full span with <code>span=True</code>. The window is split into <code>n_bins</code>, and the backend handle returns all <code>n_bins</code> summaries in one call. Those values go straight into <code>cube[i, t, :]</code>. No per-base array reaches Python for <code>mean</code>, <code>min</code> or <code>max</code>.
</figcaption></figure>

```python
import genomeblocks as gb

cre = gb.Loci.make("peaks.bed")
S = cre.signal(["signal.bw", "signal2.bw"], n_bins=10, flank=500)
# -> [INFO] Extracting 2 bigwigs for 7 loci into 10 bins (span=False, agg='mean', backend='pybigtools', exact=True, workers=1).
S.shape, S.dtype
# -> ((7, 2, 10), dtype('float32'))
```

- **Inputs.** `loci` is anything `as_loci` takes (a frame works);
  `bigwigs` is one path, a list, a `{name: path}` dict (tracks in its order),
  or open pyBigWig / pybigtools handles. Row `i` of the cube is row `i` of
  the loci, whatever the input was.
- **Chromosomes the file lacks** keep zeros: the two chr2 rows above are 0 in
  the first track, which has no chr2, and non-zero in the second.
- **Edges.** A window that runs off a chromosome end is not clipped before
  binning: the handle bins the full window grid and the bases outside have no
  data, so edge bins keep the signal they do have and the centre bin stays
  the centre bin for every row.
- **Statistics.** `mean`, `min`, `max` are the engines' native summaries;
  `sum`, `std` (population) and `coverage` (fraction of a bin's bases with
  data) are reduced from per-base values in numpy where the engine has no
  native call.
- **Memory guard.** `signal()` refuses a cube larger than half the available
  RAM instead of swapping; extract big jobs as `L[a:b].signal(...)`.

## Three readers behind one handle
{: .sec-green }

`backends.bigwig.open_bigwig(src, backend=None)` wraps a path or an open
handle in one small interface: `chroms()`, `stats_array(chrom, start, end,
n_bins=, stat=, exact=, missing=)`, `values(chrom, start, end)` and
`close()`. The extraction loop only knows that interface.

| Backend | When | How it reads |
|---|---|---|
| `pybigtools` | default (`pip install genomeblocks`) | Rust: block I/O, decompression and binning in one FFI call |
| `pybigwig` | `backend="pybigwig"` | libBigWig (C) `stats(..., nBins=, exact=)` |
| `python` | `backend="python"`, or when neither is installed | `backends/_bbi.py`: `mmap` + `zlib` + numpy; the R-tree is cached per handle and the last decompressed block is kept |

The three engines give the same numbers, because the normalisation lives in
the shared handle class rather than in each engine:

- bin `b` of a window of `n` bases is `[floor(n·b/n_bins), floor(n·(b+1)/n_bins))`,
  the integer edges libBigWig and pybigtools use, so fractional bins land in
  the same places for every engine (the pure-Python reader takes that path
  whenever `n` is not a multiple of `n_bins`);
- a window inside the chromosome with `n_bins ≤ n` goes to the engine's native
  binning; anything else (a window that leaves the chromosome, a chromosome
  the file lacks, a window shorter than its bin count) is read per base,
  padded with NaN outside the file, and reduced by one numpy routine;
- a bin with no data is `missing` (0 in a cube), a base with no data is NaN.

```python
from genomeblocks.backends.bigwig import open_bigwig

for b in ("pybigtools", "pybigwig", "python"):
    h = open_bigwig("signal.bw", backend=b)
    print(b, h.stats_array("chr1", 1000, 1300, n_bins=3, stat="mean").round(3),
          h.stats_array("chr1", -50, 100, n_bins=3, stat="mean").round(3))   # off the chromosome start
    h.close()
# -> pybigtools [8.343 1.025 8.717] [0.    6.733 0.   ]
# -> pybigwig [8.343 1.025 8.717] [0.    6.733 0.   ]
# -> python [8.343 1.025 8.717] [0.    6.733 0.   ]
```

An already-open handle keeps its engine: `cre.signal([pyBigWig.open(p)])`
reads through pyBigWig, and asking for another backend on it is an error
rather than a silent re-open.

```python
import pyBigWig
h = pyBigWig.open("signal.bw")
open_bigwig(h).backend
# -> 'pybigwig'
open_bigwig(h, backend="python")
# -> ValueError: the track is an open pybigwig handle, so backend='python' cannot apply: pass the file path instead, or drop backend=
```

`exact=False` lets pybigtools and pyBigWig answer from the zoom levels
(about 3× faster, approximate); the default is base-accurate.

## Parallel extraction
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-parallel.svg %}
</div><figcaption>
<strong>Processes, shared memory, no pickling of results.</strong> With <code>workers &gt; 1</code> the cube is allocated in <code>SharedMemory</code>. The work is split into (track range × locus range) slabs, and each worker process attaches to the cube and opens its own bigWig handles through the same backend. It writes its slab in place, and the main process copies the finished cube out once.
</figcaption></figure>

Processes rather than threads: the compiled readers serialise concurrent calls
from Python threads, so threads do not scale. `workers=1` (the default) is a
single sequential pass in chunks of 512 loci, which is fastest for typical
heatmap sizes. `workers=n` runs `n` processes, capped at the core count and at
the available work (one unit per track and per 1,000 loci). `workers=None`
uses half the cores. Open handles cannot cross a process boundary, so
`workers > 1` needs paths.

```python
L = gb.Loci.tile_genome({"chr1": 20_000, "chr2": 8_000}, 500)
a = L.signal(["signal.bw", "signal2.bw"], n_bins=8, flank=400, verbose=False, progress=False)
b = L.signal(["signal.bw", "signal2.bw"], n_bins=8, flank=400, workers=2, verbose=False, progress=False)
a.shape, bool((a == b).all())
# -> ((56, 2, 8), True)
```

## After extraction
{: .sec-green }

- **`gb.tmm(cube)`** scales tracks by edgeR's trimmed-mean-of-M-values factors
  computed on the per-region means, then to counts per million, so tracks of
  different depth become comparable. A track with no signal anywhere stays 0
  with a warning; a cube with no signal at all is an error.
- **Converters** (`genomeblocks.interop`): `cube_to_xarray(S, L, tracks,
  flank=)` labels the three axes with region names, chromosome / start / end
  and bin centres in bp; `cube_to_anndata` puts tracks in `obs`, loci in `var`
  and per-track profiles in `varm`; `cube_to_pandas` gives a region × (track,
  bin) frame.
- **Drawing** lives in `signal_draw`, so extraction never imports matplotlib.
  `plot_heatmap(loci, S, groups=...)` / `plot_profiles` / `compare_heatmap`
  take groups as a plain `{name: rows}` dict where `rows` is a `Loci` (or
  anything `as_loci` reads, matched by coordinates), a boolean mask, or row
  numbers.

```python
from genomeblocks.interop import cube_to_xarray
xa = cube_to_xarray(S, cre, ["ATAC", "H3K27ac"], flank=500)
xa.dims, xa.shape, xa.coords["bin"].values[:3]
# -> (('region', 'track', 'bin'), (7, 2, 10), array([-450., -350., -250.]))

from genomeblocks import signal_draw as sd
fig = sd.plot_heatmap(cre, S, groups={"chr1": cre["chrom"] == "chr1", "chr2": cre.tail(2)},
                      samples=["ATAC", "H3K27ac"])
```

## Costs
{: .sec-green }

- **Per window and track:** one backend call that returns `n_bins` floats;
  the Python loop costs the call overhead only. `exact=False` trades accuracy
  for the zoom levels.
- **Edges and small windows:** one per-base read plus a numpy reduction
  (cumulative sums over the window), a few microseconds more than the native
  path.
- **Workers:** each process opens its own handles; the cube crosses no pipe.
- **Memory:** `n_loci × n_tracks × n_bins × 4` bytes (`float32`), checked
  before the first read.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
