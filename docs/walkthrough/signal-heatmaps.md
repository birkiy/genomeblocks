---
title: "2. Signal heatmaps"
parent: "Example: AR & FOXA1"
layout: default
nav_order: 2
---

# Signal heatmaps
{: .no_toc }

Turning bigWig coverage into a region × track × bin **signal cube**, averaging
replicates, and drawing grouped heatmaps of AR+F vs AR−F.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Peaks vs signal

Peaks (last page) are discrete calls — *where* a factor binds. **bigWig** files
hold the continuous read coverage genome-wide — *how strong* the signal is at
every base. Heatmaps and browser tracks are built from bigWigs.

## The signal cube

`signal()` (also available as `Loci.signal`) extracts binned coverage for a set
of loci across a list of bigWigs:

```python
regions = ARpF + ARmF                       # rows: AR+F peaks then AR-F peaks
bw_order = ["ATAC_0h_r1", "ATAC_0h_r2", "ATAC_4h_r1", "ATAC_4h_r2",
            "AR_0h", "AR_4h", "FOXA1_0h", "FOXA1_4h"]
bigwigs = [BW[k] for k in bw_order]

S = np.nan_to_num(regions.signal(bigwigs, n_bins=200, flank=2000, workers=4))
```

The return value `S` is a NumPy array of shape **`(n_loci, n_tracks, n_bins)`**:

- **`n_loci`** rows — here every AR+F peak followed by every AR−F peak.
- **`n_tracks`** — one slab per bigWig, in the order you passed them.
- **`n_bins`** columns — each region is split into `n_bins` equal bins.

Key arguments:

| arg | meaning |
|---|---|
| `n_bins=200` | bins per region (heatmap column resolution) |
| `flank=2000` | window is the locus **center ± `flank`** bp (set `span=True` to use the full interval instead) |
| `workers=4` | parallel extraction across processes; the default `1` is a fast single Rust pass |

{: .note }
> bigWigs can have gaps (no coverage); `np.nan_to_num` turns those `NaN`s into 0
> so the heatmap and any averaging behave.

## Grouping replicates by averaging columns

The two ATAC replicates are two adjacent track slabs. Averaging them is just a
mean over the track axis — exactly what the browser does internally for a list
of bigWigs. We collapse 8 tracks into 6 conditions:

```python
S6 = np.stack([
    S[:, [0, 1], :].mean(1),   # ATAC 0h  (rep1 + rep2)
    S[:, [2, 3], :].mean(1),   # ATAC 4h  (rep1 + rep2)
    S[:, 4, :], S[:, 5, :],    # AR 0h, AR 4h
    S[:, 6, :], S[:, 7, :],    # FOXA1 0h, FOXA1 4h
], axis=1)
samples = ["ATAC 0h", "ATAC 4h", "AR 0h", "AR 4h", "FOXA1 0h", "FOXA1 4h"]
```

`S[:, [0, 1], :].mean(1)` averages the two ATAC-0h slabs over axis 1, leaving
`(n_loci, n_bins)`; `np.stack(..., axis=1)` reassembles the six condition slabs
into a `(n_loci, 6, n_bins)` cube.

## Drawing the heatmap

`plot_heatmap` (in `genomeblocks.signal_draw`, also `Loci.plot_heatmap`) renders
one heatmap column per track, optionally split into row groups:

```python
from genomeblocks.signal_draw import plot_heatmap

vmax = float(np.percentile(S6, 99))         # robust color cap
fig = plot_heatmap(regions, S6,
                   groups={"AR+F": ARpF, "AR-F": ARmF},
                   sets=["AR+F", "AR-F"],
                   samples=samples,
                   vmax=vmax, ymax=vmax)
```

How the arguments shape the figure:

- **`groups`** — a plain `dict[str, Loci]`. Each group becomes a block of rows;
  membership is by `uid`, so `regions` is split into its AR+F and AR−F halves.
  (There is no separate "tags" object — grouping is just a dict of `Loci`.)
- **`sets`** — the order (top → bottom) the groups are drawn in.
- **`samples`** — column labels, one per track in `S6`.
- **`vmax` / `ymax`** — color scale cap and profile-track y-max. Using the 99th
  percentile keeps a few hot bins from washing out the rest.

{: .tip }
> To make tracks comparable across the figure, TMM-normalise the cube first with
> `from genomeblocks import tmm; S = tmm(S)` before grouping. The notebook leaves
> this commented out so the raw coverage is shown; flip it on when comparing
> libraries of different depth.

## Reading the result

Rows are AR peaks, split into the FOXA1-dependent (AR+F) and FOXA1-independent
(AR−F) blocks; columns are the six conditions. The interesting contrasts:

- **ATAC 0 h vs 4 h** — is the site already open before androgen, or does it open
  on stimulation? FOXA1-pioneered (AR+F) sites tend to be pre-accessible.
- **AR 0 h vs 4 h** — confirms the DHT-induced AR gain that defines the cistrome.
- **FOXA1 0 h vs 4 h** — by construction, strong in AR+F and weak in AR−F.

Next: **[Genomic annotation →]({{ '/walkthrough/annotation/' | relative_url }})**
