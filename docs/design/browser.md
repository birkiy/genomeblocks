---
title: Browser
parent: Design
layout: default
nav_order: 8
---

# Browser
{: .no_toc }

`browser(region, tracks)` draws an IGV-like figure of one region with
matplotlib: one axis per track, all sharing the same x range. Each track is
read through the same boundary and backends as the rest of the package, and
only the data inside the region is read.
{: .fs-5 .fw-300 }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/browser-tracks.svg %}
</div><figcaption>
<strong>The track decides the drawer.</strong> <code>tracks</code> is an ordered <code>{name: source}</code> dictionary. Each source's type is detected from its Python type or file extension, and the matching drawer fills that track's axis. Intervals go through <code>as_loci</code>, loops through <code>as_pairs</code>, bigWigs through the bigwig backend as binned summaries of the region, and BAMs as a per-base pileup of the window.
</figcaption></figure>

| Source | Detected as | Drawn as |
|---|---|---|
| `.bw` / `.bigwig` path, an open pyBigWig / pybigtools handle | `bw` | area of binned means: one `stats_array` call on the region |
| a list of bigWigs | `bw` | one track: bin means averaged across the replicates |
| `.bam` | `bam` | per-base coverage; bases that differ from `reference` are coloured (needs pysam) |
| `Pairs`, `.bedpe` | `bedpe` | half-sine arcs between anchor midpoints, width by `score` |
| `Genes` | `genes` | stacked gene models from the three tables; CDS taller than exons |
| anything else: `Loci`, a frame, `.bed` / `.narrowPeak`, a list of regions, a region string | `bed` | rectangles of the rows `overlap_rows(region)` returns |

```python
import genomeblocks as gb
from genomeblocks.browserview import _detect_track_type

cre, genes, pairs = gb.Loci.make("peaks.bed"), gb.Genes.make("genes.gtf"), gb.Pairs.make("loops.bedpe")
[_detect_track_type(t) for t in (genes, pairs, "signal.bw", ["signal.bw", "signal2.bw"], cre, cre.to_pandas(), "chr1:1-2")]
# -> ['genes', 'bedpe', 'bw', 'bw', 'bed', 'bed', 'bed']

fig, axes = gb.browser("chr1:0-12 kb", {
    "ATAC": "signal.bw",
    "ATAC x2": ["signal.bw", "signal2.bw"],
    "CRE": cre,
    "frame": cre.to_pandas(),
    "loops": pairs,
    "genes": genes,
    "regions": ["chr1:500-700", "chr1:3000-3500"],
}, bw_share=[["ATAC", "ATAC x2"]])
list(axes)
# -> ['ATAC', 'ATAC x2', 'CRE', 'frame', 'loops', 'genes', 'regions', '_axis']
axes["ATAC"].get_ylim() == axes["ATAC x2"].get_ylim()
# -> True
```

- **How each drawer reads.** Intervals: `as_loci(track)` once, then
  `overlap_rows(chrom, start, end)` through the cached lookup index, and one
  `PatchCollection` of the clipped rows, no Python loop per rectangle.
  bigWigs: `open_bigwig(track, backend=)` and one `stats_array(chrom, start,
  end, n_bins=bw_n_bins, stat="mean")` per file; the three bigwig engines
  bin alike, so `backend="python"` draws the same track. Loops: `as_pairs`,
  the cis pairs on the shown chromosome whose arc crosses the view, drawn as
  one `LineCollection`. Genes: the gene rows from `genes.genes.overlap_rows`,
  their transcripts by `transcripts['gene']`, their exons and CDS by
  `features['transcript']`, stacked greedily into rows;
  `genes_max_transcripts=1` keeps the canonical (else longest) isoform.
  BAM: `bam.pileup_counts` (pysam) gives a `4 × n` base-count array, the
  reference bases come from the FASTA, and mismatching reads above
  `bam_allele_freq` are stacked in IGV's nucleotide colours.
- **Layout.** Track heights come from per-type defaults (loops and genes are
  taller), and `track_heights=` overrides them. A coordinate ruler sits at
  the bottom (`axes['_axis']`), with a reference-sequence strip above it
  when `reference` is given (`axes['_sequence']`: letters up to 200 bp, a
  colour strip up to 5 kb).
- **Scales.** bigWig and BAM tracks can share a y-limit (`bw_share`,
  `bam_share`: groups of track names that take the group's region maximum)
  or take fixed values (`bw_ymax`, `bam_ymax`, scalar or per track), so
  samples are comparable at a glance.
- **Regions** are anything `parse_region` reads: `"chr8:127.7-128.1 Mb"`,
  `"chr1:1,000-2,000"`, `(chrom, start, end)` or a `Locus`.
- **Output.** Every element is a vector patch or line, nothing is rasterised,
  so the SVG export stays clean for figures.

```python
fig, axes = gb.browser(("chr1", 0, 12_000), {"genes": genes, "ATAC": "signal.bw"},
                       genes_max_transcripts=1, bw_ymax=5, colors={"ATAC": "#aa3377"}, backend="python")
tuple(map(float, axes["ATAC"].get_ylim()))
# -> (0.0, 5.0)
with gb.use_backend(bigwig="pybigwig"):
    fig, axes = gb.browser("chr1:0.5-1.5 kb", {"ATAC": "signal.bw"})
axes["_axis"].get_xlabel()
# -> 'chr1:500-1,500'
fig, axes = gb.browser("chr1:950-1150", {"reads": "reads.bam", "CRE": cre}, reference="genome.fa")
list(axes)
# -> ['reads', 'CRE', '_sequence', '_axis']
```

{: .note }
`gb.View` (one-file interactive view) and `gb.igv_html` (an igv.js page)
share the same inputs: `Loci` or anything `as_loci` takes, `Pairs` /
BEDPE, `Genes` (drawn from `to_bed12` / `representative`, 0-based) and
bigWigs through the bigwig backend. They write their data into a single HTML
file and need no server.

## Costs

- **Per interval track:** one lookup into the cached point index and one
  patch collection.
- **Per bigWig:** one backend call returning `bw_n_bins` means (1,000 by
  default); the file is opened and closed per draw, so a long session can
  pass open handles instead.
- **Per BAM:** one pysam pileup over the window; this is the only drawer that
  reads per base, which is why it needs a region of browser size.
- **Genes:** one `overlap_rows` on the genes table and two `isin` selections
  on the linked tables.
