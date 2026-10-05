---
title: signal
parent: API Reference
layout: default
nav_order: 6
---

# `genomeblocks.signal`
{: .no_toc }

bigWig signal over loci as a cube of `(rows x tracks x bins)`, TMM
normalisation, and the plotting modules that draw a cube:
`genomeblocks.signal_draw` (heatmaps and profiles) and
`genomeblocks.motifs_draw` (the motif cube and sequence logos). See the
[Signal guide]({{ '/guide/signal/' | relative_url }}) for the workflow.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `signal(loci, bigwigs, ...)`
{: #signal .sec-navy }

```python
from genomeblocks.signal import signal

signal(loci, bigwigs, *,
       n_bins: int = 200,
       flank: int = 3_000,
       agg: str = "mean",            # 'mean' | 'min' | 'max' | 'sum' | 'std' | 'coverage'
       dtype=np.float32,
       progress: bool = True,
       workers: int | None = 1,      # 1 = sequential; >1 = that many processes; None = half the cores
       span: bool = False,           # True: bin the whole locus, not centre ± flank
       verbose: bool = True,
       backend: str | None = None,   # 'pybigtools' (default) | 'pybigwig' | 'python'
       exact: bool = True) -> np.ndarray      # shape (n_loci, n_tracks, n_bins)
```

Also available as `Loci.signal(bigwigs, **kw)`.

| Argument | Meaning |
|---|---|
| `loci` | Anything [`as_loci`]({{ '/api/interop/' | relative_url }}) takes: a `Loci`, a frame, a BED path, a list of regions. |
| `bigwigs` | One bigWig or several: paths, open pyBigWig / pybigtools handles, or a `{name: path}` dict (tracks in its order). |
| `n_bins` | Bins per locus. |
| `flank` | bp each side of the locus centre (ignored with `span=True`). |
| `agg` | Per-bin statistic. |
| `span` | Bin the whole locus `[start, end)` instead of centre ± flank. |
| `workers` | `1` reads sequentially, one native call per region and track. `>1` uses that many processes writing into a shared-memory cube (needs paths, not handles); the count is capped at `n_tracks x ceil(n_loci / 1000)` and at the core count. `None` = half the cores. |
| `backend` | The [bigwig backend]({{ '/backends/' | relative_url }}); every engine gives the same numbers. |
| `exact` | Base-accurate binning. `False` lets pybigtools / pyBigWig use zoom levels (about 3x faster, approximate). |

Row `i` of the cube is row `i` of the loci, so the cube lines up with every
other table built from them. Rows on chromosomes a bigWig lacks stay 0; a
window that leaves the chromosome keeps the signal of the bases it has.
Raises `MemoryError` before reading when the cube would take more than half
of the available RAM (extract in chunks of loci), and `ValueError` for an
empty input or an unknown `agg`.

```python
import genomeblocks as gb
cre = gb.Loci.make("peaks.bed", keep=True)

S = signal(cre, ["signal.bw", "signal2.bw"], n_bins=20, flank=500)
# -> [INFO] Extracting 2 bigwigs for 7 loci into 20 bins (span=False, agg='mean', backend='pybigtools', exact=True, workers=1).
S.shape, S.dtype
# -> ((7, 2, 20), dtype('float32'))

S2 = cre.signal({"ATAC": "signal.bw", "H3K27ac": "signal2.bw"}, n_bins=20, flank=500, verbose=False, progress=False)
np.array_equal(S, S2)
# -> True
cre.signal("signal.bw", span=True, n_bins=1, agg="max", verbose=False, progress=False).shape   # one max per locus
# -> (7, 1, 1)

import pyBigWig
h = pyBigWig.open("signal.bw")
cre.signal([h], n_bins=20, flank=500, verbose=False, progress=False)        # an open handle is a track too
cre.signal("signal.bw", agg="median", verbose=False, progress=False)
# -> ValueError: unknown stat 'median'; use mean, min, max, sum, std or coverage
cre.signal([h], workers=2, verbose=False, progress=False)
# -> ValueError: workers > 1 needs bigWig paths (open handles cannot be shared between processes); pass the file paths, or workers=1.
```

The cube converts with `genomeblocks.interop`:
`cube_to_xarray(S, L, tracks=None, *, flank=None, name="signal")` (dims
`region`, `track`, `bin`, with `chrom` / `start` / `end` coordinates),
`cube_to_anndata(S, L, tracks=None, *, agg="mean")` and
`cube_to_pandas(S, L, tracks=None)` — see
[interop]({{ '/api/interop/' | relative_url }}).

---

## `tmm(cube)`
{: #tmm .sec-navy }

```python
from genomeblocks.signal import tmm      # also gb.tmm
tmm(cube: np.ndarray) -> np.ndarray
```

TMM-normalise a `(regions x tracks x bins)` cube: edgeR's trimmed mean of
M-values per track (Robinson & Oshlack 2010) plus library-size scaling to
counts per million, so tracks become comparable. Input and output have the
same shape. A track with no signal anywhere is left at 0 with a
`RuntimeWarning`; a cube with no signal at all raises `ValueError`.

```python
N = gb.tmm(S)
N.shape, np.isfinite(N).all()
# -> ((7, 2, 20), True)
S0 = S.copy(); S0[:, 1] = 0
gb.tmm(S0)
# -> RuntimeWarning: tmm: track(s) [1] have no signal and are left at 0
gb.tmm(np.zeros((5, 2, 3)))
# -> ValueError: tmm: every track is empty (no signal in any region)
```

---

## `genomeblocks.signal_draw`
{: .sec-purple }

Heatmaps and average profiles of a cube. Row groups are a plain
`{name: rows}` dict, where `rows` is a `Loci` (or anything `as_loci` reads —
matched to the plotted loci by chromosome, start and end), a boolean mask over
the rows, or row numbers.

| Function | Returns | One line |
|---|---|---|
| `group_mask(loci, rows)` | `bool` array | One group as a mask over `loci`. |
| `plot_heatmap(loci, S, ...)` | `Figure` | deeptools-style heatmap, one column per track, a mean profile on top. |
| `plot_profiles(loci, S, ...)` | `Figure` | Mean profile of one track per group, side by side. |
| `compare_heatmap(a, b, bigwigs, ...)` | `(fig, union, S, groups)` | Two interval sets as one grid: a-only, common, b-only. |

### `group_mask`

```python
from genomeblocks import signal_draw as sd
sd.group_mask(loci, rows) -> np.ndarray
```

```python
np.flatnonzero(sd.group_mask(cre, [0, 2]))                      # row numbers
# -> array([0, 2])
np.flatnonzero(sd.group_mask(cre, cre.head(2).to_pandas()))     # intervals, matched by coordinates
# -> array([0, 1])
np.flatnonzero(sd.group_mask(cre, cre["score"] > 30))           # a mask
# -> array([3, 4, 5, 6])
sd.group_mask(cre, np.array([True, False]))
# -> ValueError: a group mask has 2 values for 7 rows
```

### `plot_heatmap`

```python
sd.plot_heatmap(loci, S, *,
                groups: dict[str, object] | None = None,   # {name: Loci | mask | rows}; None = one block
                sets: list[str] | None = None,             # block order (default: the keys of groups)
                samples: list[str] | None = None,          # column titles, one per track
                colors: dict[str, tuple] | None = None,    # profile colour per block
                ymax=10, ymin=0,                           # profile y-limits: scalar or one per track
                height: int = 3000,                        # bp from the centre to the edge (x tick labels)
                cmap="Blues", vmax=10,                     # heatmap colours and ceiling: scalar or one per track
                profile: bool = True,                      # mean profile of each block above the heatmaps
                sort: str | None = "group",                # 'group' | 'global' | None
                dpi: int = 100,
                ylabel: str = "CPM signal") -> matplotlib.figure.Figure
```

Also available as `Loci.plot_heatmap(S, **kw)`. `sort='group'` orders rows
by mean signal within each block, `'global'` across blocks, `None` keeps the
input order. `S` must have as many rows as `loci`; `samples`, `cmap`,
`vmax`, `ymax` and `ymin` must be scalars or have one entry per track.

```python
fig = sd.plot_heatmap(cre, S, samples=["ATAC", "H3K27ac"], vmax=[8, 8], ymax=8, height=500)
fig = cre.plot_heatmap(S, groups={"strong": cre["score"] > 30, "weak": cre["score"] <= 30},
                       samples=["ATAC", "H3K27ac"], height=500)
fig = sd.plot_heatmap(cre.to_pandas(), S, groups={"up": np.arange(3), "down": cre.tail(4)},
                      samples=["ATAC", "H3K27ac"], cmap=["Blues", "Reds"], sort="global", profile=False, height=500)
sd.plot_heatmap(cre, S, samples=["only one"])
# -> ValueError: samples must have length 2 (one label per track), got 1
sd.plot_heatmap(cre, S[:3])
# -> ValueError: S must be (rows, tracks, bins) with 7 rows, got shape (3, 2, 20)
```

### `plot_profiles`

```python
sd.plot_profiles(loci, S, *,
                 groups: dict[str, object] | None = None,
                 sets: list[str] | None = None,
                 colors: dict[str, tuple] | None = None,
                 ylim: float | None = None,       # default: the 99th percentile of the mean profile
                 dpi: int = 100,
                 height: int = 3000,
                 track: int = 0) -> matplotlib.figure.Figure
```

Also available as `Loci.plot_profiles(S, **kw)`. One panel per group, all
for the same `track`.

```python
fig = sd.plot_profiles(cre, S, groups={"strong": cre["score"] > 30, "weak": cre["score"] <= 30}, track=1, height=500)
len(fig.axes)
# -> 2
```

### `compare_heatmap`

```python
sd.compare_heatmap(a, b, bigwigs, *,
                   a_name: str = "A", b_name: str = "B", common_name: str = "common",
                   sets: list[str] | None = None,
                   samples=None,                   # list of column titles, or {column: [track indices]} to average replicates
                   n_bins: int = 200, flank: int = 3_000, agg: str = "mean",
                   normalize: bool = True,         # tmm() before plotting
                   cmap="Blues", vmax=10, ymax=10, ymin=0,
                   profile: bool = True, sort: str | None = "group",
                   colors=None, dpi: int = 100,
                   S: np.ndarray | None = None,    # a cube already extracted for the union
                   signal_kw: dict | None = None)  # extra keywords for signal()
    -> (fig, union_loci, S, groups)
```

Re-exported as `genomeblocks.compare_heatmap`. `a` and `b` are anything
`as_loci` takes. The cube is extracted for their union in the order a-only,
common, b-only; `groups` maps each block name to a boolean mask over
`union_loci`, so the returned cube and loci can be re-plotted or saved.

```python
a, b = cre.head(4), cre.tail(5)
fig, union, S_u, groups = sd.compare_heatmap(a, b, ["signal.bw", "signal2.bw"], n_bins=20, flank=500, normalize=False)
len(union), S_u.shape, {k: int(v.sum()) for k, v in groups.items()}
# -> (7, (7, 2, 20), {'A': 2, 'common': 2, 'B': 3})
fig, *_ = gb.compare_heatmap(a.to_pandas(), b, ["signal.bw", "signal2.bw"], n_bins=20, flank=500,
                             samples={"mean": [0, 1]})                      # two replicates, one column
```

---

## `genomeblocks.motifs_draw`
{: .sec-purple }

The motif cube of [`scan_motifs_profile`]({{ '/api/motifs/' | relative_url }}#scan_motifs_profile)
drawn like a signal heatmap, and sequence logos for PFMs and archetypes.
Logos need `logomaker` (`pip install genomeblocks[viz]`); it is imported
inside each function.

| Function | Returns | One line |
|---|---|---|
| `plot_motif_heatmap(loci, M, names, ...)` | `Figure` | Heatmap + profile of a `(rows x motifs x bins)` cube, one column per motif. |
| `plot_archetype(pfm, ...)` | `Axes` | One `(4, W)` PFM as an information-content logo (bits). |
| `plot_archetypes(archetypes, ...)` | `Figure` | A grid of logos, one per archetype. |
| `plot_cluster_members(archetype_pfm, member_pwms, member_names, ...)` | `Figure` | The archetype on top, each member underneath (QC). |
| `plot_dendrogram(Z, names, ...)` | `Axes` | The clustering tree of `cluster_motifs`, with the cutoff line. |

### `plot_motif_heatmap`

```python
gb.plot_motif_heatmap(loci, M, names=None, *,
                      r: int = 500,            # half-window used for the scan (x tick labels)
                      groups=None,             # {name: Loci | mask | rows}, as plot_heatmap
                      vmax=None, ymax=None,    # default: a per-motif scale (98th percentile of non-zero cells; top of the profiles)
                      cmap="Purples",
                      ylabel="motif hits/bin",
                      **kw)                    # passed to signal_draw.plot_heatmap
    -> matplotlib.figure.Figure
```

Re-exported as `genomeblocks.plot_motif_heatmap`. `M` must be
three-dimensional; `names` default to `motif_0`, `motif_1`, ...

```python
M, names = cre.scan_motifs_profile("genome.fa", "motifs.jaspar", r=100, n_bins=20, threshold=7.0, verbose=False)
M.shape, names
# -> ((7, 2, 20), ['M1', 'M2'])
fig = gb.plot_motif_heatmap(cre, M, names, r=100)
fig = gb.plot_motif_heatmap(cre, M, names, r=100, groups={"strong": cre["score"] > 30, "weak": cre["score"] <= 30}, cmap="Greens")
gb.plot_motif_heatmap(cre, M[:, 0])
# -> ValueError: M must be (rows, motifs, bins), got shape (7, 20)
```

### Logos

```python
from genomeblocks import motifs_draw as md

md.plot_archetype(pfm, ax=None, *, title=None, alphabet="ACGT", show_xticks=True, ylim=(0, 2)) -> Axes
md.plot_archetypes(archetypes, *, ncols=4, figsize_per=(2.8, 1.3), members=None, sort_by_size=True) -> Figure
md.plot_cluster_members(archetype_pfm, member_pwms, member_names, *, ncols=3, figsize_per=(2.8, 1.3),
                        archetype_title="ARCHETYPE") -> Figure
md.plot_dendrogram(Z, names=None, *, cutoff=None, ax=None, color_threshold=None, leaf_font_size=6,
                   no_labels=False) -> Axes
```

`archetypes` is the `{name: pfm}` dict of
[`build_archetypes`]({{ '/api/motifs/' | relative_url }}#build_archetypes);
with `members=` (its `{name: [motif names]}` dict) each title carries the
cluster size and the grid is sorted by it. `plot_dendrogram` takes the
linkage matrix `Z` of `cluster_motifs` / `build_archetypes` and draws a red
line at `cutoff`.

```python
from genomeblocks import motifs as gm
res = gm.build_archetypes("motifs.jaspar", cutoff=0.5, min_overlap=3, verbose=False)
ax  = md.plot_archetype(res["archetypes"]["ARCH_01"], title="ARCH_01")
fig = md.plot_archetypes(res["archetypes"], members=res["members"], ncols=2)
one = gm.archetype_from_names("motifs.jaspar", ["M1", "M2"], min_overlap=3, verbose=False)
fig = md.plot_cluster_members(one["archetype"], one["pwms"], one["members"], archetype_title="M1+M2")
ax  = md.plot_dendrogram(res["Z"], res["names"], cutoff=0.5)
ax.get_ylabel()
# -> 'SW distance'
```
