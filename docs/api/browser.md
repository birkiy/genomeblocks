---
title: browser
parent: API Reference
layout: default
nav_order: 7
---

# `genomeblocks.browserview`
{: .no_toc }

Three ways to look at a region: `browser()` draws it with matplotlib,
`igv_html()` writes a one-file IGV page, and `View` builds a one-file
interactive viewer of genomeblocks tables. `genomeblocks.bam` holds the
pileup helpers the browser uses for `.bam` tracks. See the
[Browser guide]({{ '/guide/browser/' | relative_url }}) for the workflow.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
gb.browser(...)      # genomeblocks.browserview.browser
gb.igv_html(...)     # genomeblocks.igv.igv_html
gb.View(...)         # genomeblocks.view.View
gb.coverage(...)     # genomeblocks.bam.coverage
```

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `browser(region, tracks, ...)`
{: #browser .sec-purple }

```python
gb.browser(region, tracks, *,
           figsize: tuple[float, float | None] = (10, None),   # height None: derived from the track heights
           dpi: int = 100,
           track_heights: dict[str, float] | None = None,      # per-track height overrides (inches)
           colors: dict[str, str] | None = None,               # per-track colour overrides
           bw_n_bins: int = 1000,                              # bins requested per bigWig track
           bw_ymax: float | dict[str, float] | None = None,    # bigWig y-maximum: one for all, or per track
           bw_share: Sequence[Sequence[str]] | None = None,    # groups of bigWig tracks that share one y-scale
           reference: str | None = None,                       # indexed FASTA (.fa + .fai); needed for .bam tracks
           show_sequence: bool = True,                         # with reference: a sequence track at the bottom
           bam_min_baseq: int = 15,
           bam_allele_freq: float = 0.2,
           bam_ymax: float | dict[str, float] | None = None,
           bam_share: Sequence[Sequence[str]] | None = None,
           label_fontsize: int = 8,
           hspace: float = 0.15,
           genes_max_transcripts: int | None = None,           # isoforms per gene (1 = canonical / longest only)
           backend: str | None = None)                         # bigwig backend: 'pybigtools' | 'pybigwig' | 'python'
    -> (matplotlib.figure.Figure, dict[str, matplotlib.axes.Axes])
```

Returns `(fig, axes_by_name)`. The axes dict has one entry per track, plus
`'_axis'` for the coordinate ruler and `'_sequence'` when a reference
sequence track is drawn. Every panel shares the x axis; all output is
vector (no rasterised patches), and only binned values are read from bigWig
files, so wide regions stay cheap.

### Region

`region` is a `Locus`, a `(chrom, start, end)` tuple or a string:
`'chr6:122,600,000-122,800,000'`, `'chr8:127.7-128.1 Mb'`, `'chr1:0-12 kb'`.
`end` must be greater than `start`.

### Track kinds

`tracks` is an ordered `{name: source}` dict; the kind is detected from the
source.

| Source | Kind | Drawn as |
|---|---|---|
| bigWig path (`.bw` / `.bigwig`, `.gz` too) or an open pyBigWig / pybigtools handle | `bw` | Binned coverage: a filled step line. |
| A list of bigWigs | `bw` | The per-bin mean of the list — one track per replicate group. |
| `.bam` path | `bam` | Per-base depth in gray; positions where the non-reference allele fraction reaches `bam_allele_freq` get nucleotide-coloured bars (IGV colours). Needs `reference` and a `.bai` index alongside. |
| `.bedpe` path or a [`Pairs`]({{ '/api/bedpe/' | relative_url }}) | `bedpe` | Half-sine arcs between the anchor midpoints of cis pairs; line width follows the `score` column. |
| A [`Genes`]({{ '/api/genes/' | relative_url }}) | `genes` | Stacked gene models from the three tables: body line, thin exon boxes, taller CDS boxes, a label with the strand arrow per gene. |
| `.narrowPeak` path | `narrowPeak` | A row of rectangles. |
| Anything else [`as_loci`]({{ '/api/interop/' | relative_url }}) takes: a `Loci`, a BED path, a frame, a list of region strings | `bed` | A row of rectangles. |

Per-kind defaults:

| Kind | Height | Colour |
|---|---|---|
| `bed` / `narrowPeak` | 0.2 | `#444444` |
| `bw` | 0.5 | `#4c78a8` |
| `bam` | 0.6 | `#a6a6a6` |
| `bedpe` | 1.5 | `#888888` |
| `genes` | 1.5 | `#000000` |
| sequence | 0.2 | IGV nucleotide colours |

The sequence track shows coloured letters at 200 bp or less, a coloured strip
up to 5 kb, and stays blank beyond; its height is set with
`track_heights={'_sequence': ...}`.

```python
cre   = gb.Loci.make("peaks.bed", keep=True)
genes = gb.Genes.make("genes.gtf")
pairs = gb.Pairs.make("loops.bedpe")

fig, axes = gb.browser("chr1:0-12 kb", {
    "ATAC": "signal.bw",
    "H3K27ac (2 reps)": ["signal.bw", "signal2.bw"],
    "CREs": cre,
    "peaks": "peaks.bed",
    "loops": pairs,
    "regions": ["chr1:500-700", "chr1:3000-3500"],
    "genes": genes,
}, bw_share=[["ATAC", "H3K27ac (2 reps)"]])
list(axes)
# -> ['ATAC', 'H3K27ac (2 reps)', 'CREs', 'peaks', 'loops', 'regions', 'genes', '_axis']
axes["ATAC"].get_ylim() == axes["H3K27ac (2 reps)"].get_ylim()      # shared scale
# -> True

fig, axes = gb.browser(("chr1", 0, 12_000), {"genes": genes, "ATAC": "signal.bw"},
                       genes_max_transcripts=1, bw_ymax=5, colors={"ATAC": "#d1495b"},
                       track_heights={"genes": 1.0}, backend="python")
fig, axes = gb.browser("chr1:950-1150", {"reads": "reads.bam", "CREs": cre}, reference="ref.fa")
list(axes)
# -> ['reads', 'CREs', '_sequence', '_axis']

gb.browser("chr1:100-50", {"CRE": cre})
# -> ValueError: Invalid region: end (50) must be > start (100)
gb.browser("chr1:900-1200", {"bam": "reads.bam"})
# -> ValueError: BAM tracks need a `reference` FASTA for mismatch colouring; pass reference='genome.fa'.
```

{: .tip }
> `bw_share=[["AR 0h", "AR 4h"]]` gives the named bigWig tracks one y-scale
> (the group's region maximum); an explicit `bw_ymax` wins. `bam_share` and
> `bam_ymax` do the same for BAM tracks.

---

## `genomeblocks.bam`
{: .sec-navy }

Pileup extraction for the browser, matplotlib-free. `pysam` is imported on
first use (`pip install genomeblocks[bam]`); a `chr1` / `1` naming mismatch
between the request and the file is resolved automatically.

| Function | Signature | Returns |
|---|---|---|
| `pileup_counts` | `pileup_counts(bam, chrom, start, end, *, min_baseq=15)` | `int64` array `(4, end - start)`: A, C, G, T depth per base. |
| `reference_seq` | `reference_seq(fasta, chrom, start, end)` | Upper-case bases over `[start, end)`, padded with `N` past a contig edge. |
| `coverage` | `coverage(bam, region, *, min_baseq=15)` | Total depth per base over a region (`Locus`, tuple or string); also `gb.coverage`. |

`bam` and `fasta` are paths (opened and closed inside) or open
`pysam.AlignmentFile` / `pysam.FastaFile` handles. Bases below `min_baseq`
and unmapped / secondary / QC-fail / duplicate reads are skipped — the reads
IGV drops from its coverage track.

```python
from genomeblocks.bam import pileup_counts, reference_seq, coverage

c = pileup_counts("reads.bam", "chr1", 1000, 1100)
c.shape, c.sum(axis=0).max()
# -> ((4, 100), 14)
reference_seq("ref.fa", "chr1", -4, 8)
# -> 'NNNNACGTACGT'
gb.coverage("reads.bam", "chr1:1,000-1,100")[:12].tolist()
# -> [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
```

---

## `igv_html(path, *, regions, ...)`
{: #igv_html .sec-purple }

```python
gb.igv_html(path, *,
            regions: Sequence,                        # loci offered as one-click buttons
            loci: dict | None = None,                 # {track name: intervals}, embedded genome-wide
            genes=None,                               # a Genes: one BED12 model per gene, genome-wide
            signal: dict[str, str] | None = None,     # {track name: bigWig path or open handle}, around the regions only
            architecture=None,                        # an Architecture: loops touching the regions, as arcs
            score: str = "n",                         # edge column shown as arc height
            chrom_sizes: dict[str, int] | None = None,
            genome_id: str | None = None,             # e.g. 'hg38': igv.js's hosted genome (needs internet)
            flank: int = 250_000,                     # signal and loops are embedded for region ± flank
            bin_size: int = 50,                       # signal bin (bp)
            title: str = "Genome browser",
            notes: dict[str, str] | None = None,      # {region: one-line note} next to each button
            igv_js: str = "cdn",                      # 'cdn' | path to a local igv.min.js (embedded) | any URL
            standalone: bool = True,                  # False: page content only, for a host that adds <html>/<head>
            subtitle: str | None = None,
            backend: str | None = None) -> dict       # bytes per embedded track, plus 'total'
```

One self-contained HTML file: every track is a gzipped data URI and igv.js
draws it, so the reader needs only a web browser. A region string may hold
two loci (`'chr1:0-5 kb chr2:0-6 kb'`) to open a split view — how a trans
loop is shown. Without `genome_id`, the page carries chromosome sizes only
(from `chrom_sizes`, else the extent of the embedded loci plus `flank`).

```python
A = gb.Architecture.make(cre, pairs, r=100, verbose=False)
A.ep["w"][:] = [5.0, 3.0, 2.0, 1.0]            # or A.add_mcool("hic.mcool", resolution=5000)
A.normalize(verbose=False)
sizes = gb.igv_html("share.html",
                    regions=["chr1:0-12 kb", ("chr1", 0, 5000), "chr1:0-5 kb chr2:0-6 kb"],
                    loci={"CREs": cre}, genes=genes, signal={"ATAC": "signal.bw"}, architecture=A,
                    title="Example", notes={"chr1:0-12 kb": "GENE_A and GENE_B"})
sizes
# -> {'ATAC': 2005, 'CREs': 197, 'loops (n)': 149, 'genes': 185, 'total': 6985}
```

Coordinates in the page are 1-based closed, as igv.js expects
(`chr1:0-12 kb` becomes `chr1:1-12000`); the embedded BED tracks keep
genomeblocks' 0-based starts.

---

## `View`
{: .sec-purple }

```python
gb.View(architecture=None, *,
        cre=None,                   # a Loci (or anything as_loci takes) when there is no Architecture
        genes=None,                 # a Genes
        samples=None,               # ['LNCaP', ...] or {'LNCaP': '#46a8e4', ...}
        title="genomeblocks view", subtitle="",
        anchor_r: int = 5000,       # TSS ± r picks a gene's anchor CREs
        anchor_mode: str = "center",
        gene_half: int = 500_000,   # half-width of the window a gene button opens
        chrom_sizes=None)           # dict, .chrom.sizes path or Genome; default: the tables' extent
```

Our own interactive browser for genomeblocks tables, saved as one HTML file
(no install, no server, works offline). Data travels as the columns
genomeblocks keeps in memory (gzipped typed arrays) and a small canvas viewer
draws them. Because the viewer knows the tables it can search a gene and
light up its anchor CREs, draw the anchor's contact profile per sample, list
partners with their O/E and open a trans partner side by side. Pass an
`Architecture` or `cre=`; passing neither raises `ValueError`.

Build it like a figure, one track per call; tracks appear in call order and
every method returns the `View`.

| Method | Signature | Track |
|---|---|---|
| `v.anchor_profile` | `anchor_profile(name="contacts of the anchor (O/E)", smooth=1.5, height=None)` | Per-sample sum of edge weight from the selected gene / CRE to every partner. |
| `v.cre_values` | `cre_values(name, values, smooth=1.0, height=None)` | Per-CRE numbers aligned to the CRE rows, as a profile; `{sample: array}` overlays. |
| `v.signal` | `signal(name, bigwigs, *, bin_size=50, flank=250_000, genome_bin=20_000, height=None)` | bigWig signal: fine bins around the listed regions, coarse bins genome-wide; `{sample: path}` or one path. |
| `v.intervals` | `intervals(name, sets, color=None, height=None)` | Interval rows: `{sample or label: intervals / boolean mask over the CRE rows}`, or one set. |
| `v.points` | `points(name, tables, *, chrom="chromosome", start="start", end="end", value="log2", segments=None, seg_cols=("chrom", "start", "end", "label"), ylim=None, height=None)` | Values along the genome (e.g. CNVkit `.cnr` frames or paths), with optional labelled segments. |
| `v.cres` | `cres(name="CREs", height=None)` | The CRE rows (promoters darker); the anchor and its partners are highlighted. |
| `v.loops` | `loops(score="n", name="loops", height=None)` | Edges as arcs (trans edges as labelled stubs); `score` is an edge column or `{sample: column}`. Needs an `Architecture`. |
| `v.genes` | `genes(name="genes", height=None)` | Gene models (the canonical isoform when selected, else the longest). |

Navigation and output:

| Method | Signature | One line |
|---|---|---|
| `v.region` | `region(label, locus=None, *, gene=None, note="", anchor=None)` | A one-click button. `locus` may hold two loci for a split view; `gene` anchors that gene (and sets the window when `locus` is None); `anchor` a CRE row. |
| `v.mark` | `mark(label, chrom, pos)` | A labelled position marker. |
| `v.hubs` | `hubs(values, n=12, label=None)` | The "Top hubs" list: the `n` CRE rows with the largest `values` (default: `vp.strength` when present). |
| `v.to_html` | `to_html(standalone=True)` | `(page, bytes per part)`. |
| `v.save` | `save(path, standalone=True)` | Write the page; returns bytes per part. |
| `v._repr_html_()` | | The page in an iframe — display `v` in a notebook. |

Without any track call, `save` / `to_html` build a default figure: anchor
profile, node strength (when `vp.strength` exists), CREs, loops, genes.

```python
A.annotate(genes, verbose=False).strength(verbose=False)
v = gb.View(A, genes=genes, samples={"LNCaP": "#46a8e4"}, title="MYC hubs")
v.anchor_profile("graph O/E (anchor)")
v.cre_values("node strength", {"LNCaP": A.vp["strength"]})
v.signal("ATAC", {"LNCaP": "signal.bw"}, bin_size=50, flank=2000)
v.intervals("peaks", {"LNCaP": cre})
v.intervals("strong", cre["score"] > 30, color="#d1495b")     # a mask over the CRE rows
v.cres(); v.loops("n"); v.genes()
v.region("GENE_B", gene="GENE_B", note="the lncRNA")
v.region("split", "chr1:0-5 kb chr2:0-6 kb")
v.mark("E1", "chr1", 1000)
v.hubs(A.vp["strength"], n=3, label="strength")
v.save("view.html")
# -> {'CREs': 224, 'edges': 136, 'genes': 366, 'track: node strength': 60, 'track: ATAC': 1712,
#     'track: peaks': 144, 'track: strong': 120, 'viewer (js + css)': 43118, 'total': 49726}

gb.View(cre=cre).loops()
# -> ValueError: loops need an Architecture
```

### `view_html`

```python
from genomeblocks.view import view_html
view_html(path, *, architecture=None, cre=None, genes=None, signal=None, regions=(), notes=None,
          title="genomeblocks view", subtitle="", standalone=True, **kw) -> dict
```

The one-call version with the default figure: anchor profile, node strength,
one `signal` track per `{name: bigWig}` entry, CREs, loops, genes, and one
button per region. `**kw` goes to `View(...)`.

```python
view_html("view3.html", architecture=A, genes=genes, signal={"ATAC": "signal.bw"},
          regions=["chr1:0-12 kb"], title="one call")
# -> {'CREs': 224, 'edges': 136, 'genes': 366, 'track: node strength': 60, 'track: ATAC': 1228,
#     'viewer (js + css)': 43118, 'total': 48733}
```

{: .note }
> For BAM / VCF review or huge remote files use `igv_html`; the `View` embeds
> genomeblocks tables and bigWig bins only.
