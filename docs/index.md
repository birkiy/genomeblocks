---
title: Home
layout: default
nav_order: 1
description: "Building blocks for regulatory genomics in Python, as columnar tables that work with pandas, polars, Arrow and the rest of the ecosystem."
permalink: /
---

<div class="gb-eyebrow-row">
  <span class="gb-eyebrow">v2.0 · MIT</span>
  <button type="button" class="gb-install-chip" onclick="gbCopyInstall(this)" title="Copy to clipboard">
    <span class="chip-text">pip install genomeblocks</span>
    <span class="copy-icon" aria-hidden="true">⧉</span>
  </button>
</div>

# genomeblocks
{: .fs-9 }

Building blocks for regulatory genomics — peaks, genes, loops, contact networks, signal and motifs as columnar tables that line up by row and hand themselves to pandas, polars, Arrow and friends.
{: .fs-6 .fw-300 }

<div class="gb-works-with" aria-label="Supported data types">
  <span class="label">Works with</span>
  <span class="chip">ChIP-seq · ATAC-seq</span>
  <span class="chip">Hi-C · HiChIP</span>
  <span class="chip">GTF / GFF3 / genePred</span>
  <span class="chip">bigWig · BAM · FASTA</span>
  <span class="chip">JASPAR · MEME motifs</span>
  <span class="chip">pandas · polars · Arrow</span>
</div>

[Quickstart →]({{ '/quickstart/' | relative_url }}){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/birkiy/genomeblocks){: .btn .fs-5 .mb-4 .mb-md-0 }


## Why genomeblocks?
{: .sec-navy }

A regulatory-genomics analysis is usually a cocktail of bedtools, a GTF parser, cooler, a bigWig reader, a motif scanner, a graph library and one-off heatmap code, glued together with DataFrames whose rows no longer line up. `genomeblocks` keeps every object as a **table of numpy columns**: the row number is the join key, so the signal cube, the motif matrix, the annotation labels and the contact graph computed from one set of peaks all share its row order. Everything comes in through one boundary (`gb.as_loci`) and goes out through explicit converters (`to_pandas`, `to_polars`, `to_arrow`, …); the heavy work runs through backends you can swap per call.

### Tables
{: .no_toc }

| Table | Rows | What it holds | Guide |
|---|---|---|---|
| `Loci` | one interval | `codes` / `starts` / `ends` / `strands` columns plus any extra columns (`name`, `score`, …); set algebra, `slop`, `merge`, `nearest`, `signal`, `scan_motifs`, `enrich` | [Loci]({{ '/guide/loci/' | relative_url }}) |
| `Locus` | — | one interval; `L[i]` is a `LocusView` that reads the table in place; the uid `chrom:start-end(strand)` | [Locus API]({{ '/api/locus/' | relative_url }}) |
| `Genome` | — | chromosome names ↔ integer codes, plus sizes; every table carries one | [Concepts]({{ '/concepts/' | relative_url }}) |
| `Genes` | one gene | three linked `Loci` — `genes`, `transcripts`, `features` — 0-based, from GTF / GFF3 / genePred; `annotations`, `nearest_genes`, `get_tss`, `select_isoforms` | [Genes]({{ '/guide/genes/' | relative_url }}) |
| `Pairs` | one loop | BEDPE as two row-aligned anchor `Loci` (`a`, `b`) plus columns; `overlapping`, `filter`, `distance` | [BEDPE]({{ '/guide/bedpe/' | relative_url }}) |
| `Architecture` | one edge | vertices = a `Loci`, edge table `src` / `tgt` + `ep` columns (`w`, `d`, `n`), vertex columns `vp`; `add_mcool`, `normalize`, `annotate`, `strength`, `prime_hubs`, `components`, `draw` | [Architecture]({{ '/guide/architecture/' | relative_url }}) |
| `Atlas` | one track | a sparse bins × tracks index over thousands of BED files; `search`, `bootstrap` | [Atlas]({{ '/guide/atlas/' | relative_url }}) |
| `Library` | one motif | count matrices from JASPAR / TRANSFAC / MEME / uniprobe / Biopython / arrays (`gb.load_motifs`) | [Motifs]({{ '/guide/motifs/' | relative_url }}) |

Every table has `shape`, `columns`, `head()`, `tail()`, `describe()`, a notebook `_repr_html_`, and the Arrow C stream, dataframe-interchange and narwhals protocols — `pl.DataFrame(cre)`, `duckdb.sql("select * from cre")`, `sns.histplot(data=cre, x="score")`, `alt.Chart(cre)` take a `Loci` as it is.

### Functions
{: .no_toc }

| Module | Entry points | Guide |
|---|---|---|
| `signal` | `loci.signal(bigwigs)` → a `(rows, tracks, bins)` cube; `gb.tmm`; `plot_heatmap`, `plot_profiles`, `gb.compare_heatmap` | [Signal]({{ '/guide/signal/' | relative_url }}) |
| `motifs` | `scan_motifs`, `scan_motifs_matrix`, `scan_motifs_matrix_masked`, `scan_motifs_profile`, `bootstrap_enrichment`, `compare_motifs`, `pwm_distance_matrix`, `cluster_motifs`, archetypes; `gb.plot_motif_heatmap` | [Motifs]({{ '/guide/motifs/' | relative_url }}) |
| `browser` · `igv_html` · `View` | a matplotlib region view; a one-file IGV page; a one-file interactive view of an `Architecture` | [Browser]({{ '/guide/browser/' | relative_url }}) |
| `bedpe` | `count_pairs`, `count_pairs_2d`, `read_pairs_chunks` over HiC-Pro / 4DN / Juicer pairs files | [BEDPE]({{ '/guide/bedpe/' | relative_url }}) |
| `se` · `hichip` · `coverage` | ROSE-style super-enhancers; HiChIP short-range tracks; BAM pileups | [Signal]({{ '/guide/signal/' | relative_url }}) |

### The boundary
{: .no_toc }

`gb.as_loci(x)` turns anything interval-like into a `Loci`, and every public function calls it on its inputs — so you never convert first.

| In | Out |
|---|---|
| a `Loci`, a `Locus`, a region string `'chr1:1,000-2,000'`, a list of regions or `(chrom, start, end)` tuples | `to_records()`, `to_numpy()`, `uid` |
| a path: BED / narrowPeak / broadPeak, CSV / TSV with a header, parquet (`.gz` too) | `to_bed(path)`, `save(path)` (parquet) |
| pandas, polars (eager or lazy), pyarrow, bioframe, PyRanges, pybedtools frames; dicts of columns; structured arrays | `to_pandas()`, `to_polars()`, `to_arrow()`, `to_bioframe()`, `to_pyranges()`, `to_bedtool()`, `to_cgranges()` |
| an AnnData (its `var` regions) | `to_anndata(X)` |
| a FASTA path, dict or open handle wherever sequence is read | `sequences(fasta)`, `to_seqrecords(fasta)`, `to_fasta(path, fasta)` |

Column names are lenient on the way in (`chrom` / `chr` / `Chromosome` / `seqnames`, `start` / `Start` / `chromStart`, …) and fixed inside. See [Interoperability]({{ '/interoperability/' | relative_url }}).

### Backends
{: .no_toc }

| Family | Default (installed by `pip install genomeblocks`) | Others on request |
|---|---|---|
| `intervals` — overlap, nearest, merge, point lookups | genomeblocks (numpy) | cgranges, ncls, bioframe, pyranges, bedtools |
| `bigwig` — reading signal | pybigtools | pybigwig, python |
| `motifs` — scoring matrices along sequences | lightmotif | moods, biopython |
| `fasta` — fetching sequence | genomeblocks (indexed `.fai` reader) | pysam, pyfaidx, memory, biopython |
| `tables` — parsing BED / GTF / pairs text | polars, else pandas | pandas |
| `graph` — algorithms on an `Architecture` | graph-tool, else scipy | igraph, networkx |

`gb.backends()` lists them; `backend=` on a call or `with gb.use_backend(intervals="bioframe"):` picks one; a requested backend that is not installed raises an `ImportError` with the install command — nothing falls back silently. See [Backends]({{ '/backends/' | relative_url }}).

---

## 60-second example
{: .sec-green }

Three stages; each result is a table the next stage takes as it is:

```python
import genomeblocks as gb

# 1. CREs: ATAC peaks, widened ±100 bp, merged — a table of intervals
cre = gb.Loci.make("atac_peaks.narrowPeak").slop(100).merge()

# 2. Contact graph: loops anchored on the CREs, Hi-C weights, O/E, genes, hubs
genes = gb.Genes.make("gencode.v38.annotation.gtf")
A = (gb.Architecture.make(cre, "loops.bedpe", r=2500)
       .add_mcool("hic.mcool", resolution=5000)
       .normalize()
       .annotate(genes)
       .strength())
hubs = A.prime_hubs()                 # {'prime_genes': {...}, 'hub_uids': [...], ...}

# 3. Signal, labels and the hand-off — the rows of `cre` are the join key
S = cre.signal(["ATAC.bw", "H3K27ac.bw"], n_bins=200, flank=3000)   # (rows, tracks, bins)
labels = genes.annotations(cre)                                      # one row per CRE
df = cre.to_pandas()                                                 # or to_polars(), to_arrow(), ...
```

<div class="gb-stage-rail">
  <div class="gb-stage green">
    <span class="gb-stage-num">1</span>
    <span><strong>Build intervals.</strong> Peaks → CREs as a <code>Loci</code> table.</span>
  </div>
  <div class="gb-stage purple">
    <span class="gb-stage-num">2</span>
    <span><strong>Wire the network.</strong> Loops + contacts as an <code>Architecture</code> over those rows.</span>
  </div>
  <div class="gb-stage navy">
    <span class="gb-stage-num">3</span>
    <span><strong>Annotate and hand off.</strong> Signal cube, gene labels and a DataFrame, all in the same row order.</span>
  </div>
</div>

---

## Documentation map
{: .sec-navy }

| Section | When to read it |
|---|---|
| [Installation]({{ '/installation/' | relative_url }}) | `pip install genomeblocks`, the extras, and the three conda-only engines. |
| [Quickstart]({{ '/quickstart/' | relative_url }}) | A 10-minute tour end to end, with the outputs. |
| [Concepts]({{ '/concepts/' | relative_url }}) | The mental model — tables, the row as join key, Genome codes, 0-based coordinates, the boundary, backends. |
| [Interoperability]({{ '/interoperability/' | relative_url }}) | Every input `as_loci` takes and every `to_*` exit, with round trips. |
| [Backends]({{ '/backends/' | relative_url }}) | The six backend families, how to pick one, and what each engine supports. |
| [Example: AR & FOXA1]({{ '/walkthrough/' | relative_url }}) | A complete real-data walkthrough, concept by concept. |
| [User Guide → Loci]({{ '/guide/loci/' | relative_url }}) | Interval algebra in depth. |
| [User Guide → Genes]({{ '/guide/genes/' | relative_url }}) | GTF / GFF3 / genePred tables, annotation and isoform selection. |
| [User Guide → Architecture]({{ '/guide/architecture/' | relative_url }}) | Chromatin-contact networks. |
| [User Guide → Signal]({{ '/guide/signal/' | relative_url }}) | bigWig extraction, TMM and heatmaps. |
| [User Guide → Browser]({{ '/guide/browser/' | relative_url }}) | Region plots, IGV pages and the interactive view. |
| [User Guide → BEDPE]({{ '/guide/bedpe/' | relative_url }}) | `Pairs` and contact counting from pairs files. |
| [User Guide → Motifs]({{ '/guide/motifs/' | relative_url }}) | Motif libraries, scanning and enrichment. |
| [User Guide → Atlas]({{ '/guide/atlas/' | relative_url }}) | GIGGLE-style enrichment against BED collections. |
| [API Reference]({{ '/api/' | relative_url }}) | Full signatures. |
| [Design]({{ '/design/' | relative_url }}) | How each module works inside, with diagrams. |
| [Benchmarks]({{ '/benchmarks/' | relative_url }}) | Every backend of every family timed through the genomeblocks API. |
| [Release Notes]({{ '/release-notes/' | relative_url }}) | What changed in 2.0 and the 1.x → 2.0 migration table. |
| [Credits]({{ '/credits/' | relative_url }}) | Upstream libraries and how to cite them. |

---

## Citing
{: .sec-navy }

If `genomeblocks` is useful in your work, please cite the GitHub repository:

```
Altintas, U. B. (2024). genomeblocks: Fluent building blocks for regulatory genomics.
https://github.com/birkiy/genomeblocks
```

`genomeblocks` stands on top of several excellent upstream libraries — please also cite the ones whose backend did the work in your analysis. The [Credits]({{ '/credits/' | relative_url }}) page lists them all with BibTeX.

---

## License

MIT. See [LICENSE](https://github.com/birkiy/genomeblocks/blob/main/LICENSE).

<script>
function gbCopyInstall(btn) {
  var text = btn.querySelector('.chip-text').textContent;
  if (navigator.clipboard) { navigator.clipboard.writeText(text); }
  var icon = btn.querySelector('.copy-icon');
  var original = icon.textContent;
  icon.textContent = '✓';
  btn.classList.add('copied');
  setTimeout(function() {
    icon.textContent = original;
    btn.classList.remove('copied');
  }, 1400);
}
</script>
