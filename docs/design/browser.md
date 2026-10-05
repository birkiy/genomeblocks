---
title: Browser
parent: Design
layout: default
nav_order: 8
---

# Browser
{: .no_toc }

`browser(region, tracks)` draws an IGV-like figure of one region with
matplotlib: one axis per track, all sharing the same x range.
{: .fs-5 .fw-300 }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/browser-tracks.svg %}
</div><figcaption>
<strong>The track decides the drawer.</strong> <code>tracks</code> is an ordered <code>{name: source}</code> dictionary. Each source's type is detected from its file extension or Python type, and the matching drawer fills that track's axis. Only data inside the region is read: bigWigs return binned summaries, and BAMs a per-base pileup of the window.
</figcaption></figure>

| Source | Drawn as |
|---|---|
| `.bw` / `.bigwig` | area of binned means (one native call per track) |
| list of bigWig paths | one track: bin means averaged across the replicates |
| `.bam` | per-base coverage; bases that differ from `reference` are coloured (needs pysam) |
| `.bed` / `.narrowPeak` / `Loci` | rectangles |
| `.bedpe` / `list[Pair]` | half-sine arcs between anchor midpoints |
| `Genes` | stacked gene models; CDS taller than UTR exons |

- **Layout.** Track heights come from per-type defaults (bedpe and genes are
  taller), and `track_heights=` overrides them. A coordinate ruler sits at the
  bottom, with a reference-sequence strip above it when `reference` is given
  (letters drawn up to 5 kb).
- **Scales.** bigWig and BAM tracks can share a y-limit (`bw_share`,
  `bam_share`) or take fixed values (`bw_ymax`, `bam_ymax`), so samples are
  comparable at a glance.
- **Output.** Every element is a vector patch or line, nothing is rasterised,
  so the SVG export stays clean for figures.
