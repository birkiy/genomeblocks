---
title: Signal
parent: Design
layout: default
nav_order: 3
---

# Signal
{: .no_toc }

`signal()` reads bigWig tracks into one dense array, the *signal cube*, with
one row per region, one plane per track and one column per bin. Heatmaps,
profiles and TMM normalisation all read that array.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## One region, one call, one row
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-cube.svg %}
</div><figcaption>
<strong>The cube is filled row by row.</strong> Each region becomes a window around its centre (± <code>flank</code>), or its full span with <code>span=True</code>. The window is split into <code>n_bins</code>, and the bigWig returns all <code>n_bins</code> summaries in a single native call. Those values go straight into <code>cube[i, t, :]</code>. No per-base array reaches Python for <code>mean</code>, <code>min</code> or <code>max</code>.
</figcaption></figure>

- **Edges.** A window that runs off a chromosome end is clipped, and its
  missing bins stay zero, so every row has the same length and the same
  centre bin.
- **Statistics.** `mean`, `min` and `max` are native pybigtools summaries
  (`exact=True` uses base-accurate binning, `exact=False` the zoom levels).
  `sum`, `std` and `coverage` are reduced from per-base values in numpy.
- **Memory guard.** `signal()` refuses a cube larger than half the available
  RAM instead of swapping. Extract big jobs in chunks.

## Two readers behind one interface
{: .sec-green }

| Backend | When | How it reads |
|---|---|---|
| `pybigtools` | installed (default) | Rust: block I/O, decompression and binning in one call |
| `bigwig.py` | `backend="bigwig"`, or pybigtools missing | pure Python: `mmap` + `zlib`, R-tree cached per handle, last decompressed block cached |

Both expose the same `stats_array()` call, so the extraction loop does not care
which one it has. pybigtools 0.3 renamed `missing=` to `fillna=`. The adapter
picks the right keyword for the installed version, so neither version warns.

## Parallel extraction
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/signal-parallel.svg %}
</div><figcaption>
<strong>Processes, shared memory, no pickling of results.</strong> With <code>workers &gt; 1</code> the cube is allocated in <code>SharedMemory</code>. The work is split into (track range × locus range) slabs, and each worker process attaches to the cube and opens its own bigWig handles. It writes its slab in place, and the main process copies the finished cube out once.
</figcaption></figure>

Processes rather than threads: pybigtools serialises concurrent calls from
Python threads, so threads do not scale. `workers=1` (the default) is a
single sequential pass, which is fastest for typical heatmap sizes.
`workers=n` runs `n` processes, capped at the core count and at the
available work. `workers=None` uses half the cores, which suits shared
machines.

## After extraction
{: .sec-green }

- `tmm(cube)` scales tracks by edgeR's trimmed-mean-of-M-values factors
  (computed on per-region means), so tracks with different depths become
  comparable.
- `signal_draw.plot_heatmap` / `plot_profiles` / `compare_heatmap` draw the
  cube; groups are a plain `{name: Loci}` dictionary.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/#signal' | relative_url }}) page.
