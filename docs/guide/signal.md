---
title: Signal
parent: User Guide
layout: default
nav_order: 5
---

# Signal
{: .no_toc }

Threaded bigWig extraction, TMM normalization, and comparative heatmaps for a `Loci` set. Powered by [pybigtools](https://github.com/jackh726/bigtools) (Huey, 2023) on the fast path — see [Credits](../credits#pybigtools) for the full citation.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `loci.signal()` — pull a `(loci × tracks × bins)` cube

```python
cube = loci.signal(
    bigwigs=["atac.bw", "h3k4me1.bw", "h3k27ac.bw"],
    n_bins=200,         # bins per region
    flank=3_000,        # ±flank around the region center (ignored if span=True)
    agg="mean",         # per-bin aggregator: mean / max / min / std / sum / coverage
    span=False,         # True → use the full Locus span rather than center±flank
    dtype=np.float32,
    progress=True,
    workers=1,          # 1 = sequential (fastest for typical jobs); >1 = multiprocessing; None = half the cores
    backend=None,       # None → auto-detect pybigtools; 'bigwig' forces pure-Python
    exact=True,         # pybigtools base-accurate binning (False = ~3× faster, approximate)
)
print(cube.shape)       # (n_loci, n_tracks, n_bins)
```

Returns a dense `numpy.ndarray` — immediately suitable for `tmm()`, plotting, or clustering.

### Execution model

- **Default (`workers=1`) is sequential.** A single native-Rust pass with
  `pybigtools` reaches tens of thousands of region-tracks/s — fastest for the
  typical heatmap / browser / per-locus workload.
- **`workers > 1` uses multiprocessing, not threading.** `pybigtools`
  serialises concurrent Python threads, so scaling comes from a
  `ProcessPoolExecutor` whose children write directly into one shared-memory
  cube. Work is split by (track-chunk, loci-chunk); each child opens its own
  bigWig handles.
- An explicit `workers` is honoured, capped at `min(workers, n_tracks·⌈n_loci/1000⌉, cpu_count())`.
  `workers=None` asks for half the cores, leaving the rest for the OS and each
  worker's own Rust pool — the polite choice on a shared machine.

### Backends

| Backend | Source | Notes |
|---|---|---|
| `pybigtools` (default) | Rust, via `pybigtools` | Fastest; releases the GIL; exact base-pair reads. |
| `bigwig` | Pure Python in `genomeblocks.bigwig` | Zero compiled deps; falls back when `pybigtools` isn't importable. |

---

## TMM normalization

```python
from genomeblocks import tmm
cube_n = tmm(cube)
```

Per-track TMM normalization factors computed over per-region means, then scaled to library size in per-million. Useful when pooling biological replicates or comparing cell types. The edgeR TMM algorithm (Robinson & Oshlack, 2010) is vendored directly in `genomeblocks.signal` — no external normalization dependency.

---

## `compare_heatmap` — A-specific / shared / B-specific grid

Typical task: compare enhancer sets between two conditions across multiple marks.

```python
from genomeblocks import compare_heatmap

fig, union, S, groups = compare_heatmap(
    a=cre_mesc,
    b=cre_hesc,
    bigwigs=["ATAC_mESC.bw", "ATAC_hESC.bw",
             "H3K27ac_mESC.bw", "H3K27ac_hESC.bw"],
    a_name="mESC",
    b_name="hESC",
    common_name="shared",
    n_bins=200,
    flank=3_000,
    normalize=True,          # run tmm() before plotting
    cmap=["Blues", "Blues", "Reds", "Reds"],
    vmax=[10, 10, 6, 6],
    samples={"ATAC":   [0, 1],    # merge bigwigs 0 and 1 into one column
             "H3K27ac": [2, 3]},
)
fig.savefig("cre_compare.pdf")
```

- `sets` controls row order (default `[a_name, common_name, b_name]`).
- `samples` can be a list (no merging) or a dict `{column_name: [bigwig_indices]}` to average replicates inline.
- `sort="group"` orders rows within each group by mean signal; `"global"` orders across all rows; `None` keeps input order.
- Pass a pre-computed `S` to skip extraction entirely (useful for iterating on plot params).

Returns `(fig, union_loci, S, groups)` so you can re-plot with different params.

---

## `plot_heatmap` / `plot_profiles`

The same plotting machinery as `compare_heatmap`, decoupled from Loci comparison.
Row groups are a plain `dict[str, Loci]` — there is no separate annotation
object:

```python
groups = {
    "promoter":  promoter_cre,
    "enhancer":  enhancer_cre,
    "quiescent": quiescent_cre,
}

fig = loci.plot_heatmap(cube, groups=groups,
                        sets=["promoter", "enhancer", "quiescent"],
                        cmap="Blues", vmax=8)
fig.savefig("heatmap.pdf")

fig = loci.plot_profiles(cube, groups=groups, ylim=5)
fig.savefig("profiles.pdf")
```

Omit `groups` to treat all loci as one group. `plot_heatmap` / `plot_profiles`
live in `genomeblocks.signal_draw` (and are attached as `Loci` methods); the
processing functions `signal` / `tmm` stay in `genomeblocks.signal`.

---

## Sizing considerations

Memory for the cube is `n_loci × n_tracks × n_bins × dtype_size`. The extractor refuses to allocate more than half of **available RAM** as a safety check (raising `MemoryError`). If you need to scan a million CREs × 50 tracks × 200 bins, chunk by loci and stream to disk (`np.save` per chunk).

---

## Example: ChIP-seq enrichment profile around CRE centers

```python
from genomeblocks import Loci

cre = Loci.make("cre.bed")
bigwigs = [
    "H3K4me1.bw", "H3K4me3.bw", "H3K27ac.bw",
    "H3K27me3.bw", "H3K9me3.bw",
]

cube = cre.signal(bigwigs, n_bins=200, flank=3_000, agg="mean")
cube = gb.tmm(cube)

fig = cre.plot_profiles(cube, ylim=6)
fig.savefig("marks.pdf")
```
