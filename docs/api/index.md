---
title: API Reference
layout: default
nav_order: 6
has_children: true
permalink: /api/
---

# API Reference

Signatures and short descriptions for every public name, grouped by module,
with examples that run on 2.0. Use Ctrl-F / the search box to jump around.
Everything you need day to day is one import away:

```python
import genomeblocks as gb
```

| Top-level name | What it is | Page |
|---|---|---|
| `gb.Loci`, `gb.Locus`, `gb.Genome` | Intervals as a table; one interval; chromosome names ↔ codes. | [loci]({{ '/api/loci/' | relative_url }}), [locus]({{ '/api/locus/' | relative_url }}) |
| `gb.Genes` | Genes / transcripts / features, three linked tables. | [genes]({{ '/api/genes/' | relative_url }}) |
| `gb.Pairs` | BEDPE as a table (anchors `a` and `b`). | [bedpe]({{ '/api/bedpe/' | relative_url }}) |
| `gb.Architecture` | A CRE interaction graph: vertices = a `Loci`, edges = a table. | [architecture]({{ '/api/architecture/' | relative_url }}) |
| `gb.Atlas` | A sparse bins x tracks index for enrichment. | [atlas]({{ '/api/atlas/' | relative_url }}) |
| `gb.as_loci` | Anything interval-like → `Loci`; what every function calls on its inputs. | [interop]({{ '/api/interop/' | relative_url }}) |
| `gb.read_fasta`, `gb.load_motifs` | A FASTA as `{chrom: seq}`; a motif file as a `Library`. | [backends]({{ '/api/backends/' | relative_url }}) |
| `gb.backends()`, `gb.use_backend` | List the engines; pick one for a block of code. | [backends]({{ '/api/backends/' | relative_url }}) |
| `gb.tmm`, `gb.compare_heatmap`, `gb.plot_motif_heatmap` | Signal normalisation and the two cross-set plots. | [signal]({{ '/api/signal/' | relative_url }}) |
| `gb.browser`, `gb.View`, `gb.igv_html`, `gb.coverage` | The matplotlib browser, the one-file HTML viewer, igv.js pages, BAM coverage. | [browser]({{ '/api/browser/' | relative_url }}) |

Names are imported on first use, so `import genomeblocks` stays light.

## Modules

| Module | Public names | Page |
|---|---|---|
| `genomeblocks.locus` | `Locus`, `parse_region`, `parse_regions` (and `loci.LocusView`, the row view) | [locus]({{ '/api/locus/' | relative_url }}) |
| `genomeblocks.genome` | `Genome` (`from_sizes`, `from_fasta`, `encode` / `decode`, `rank`, `sizes`), `read_sizes` | [loci]({{ '/api/loci/' | relative_url }}#genome), [Concepts]({{ '/concepts/' | relative_url }}) |
| `genomeblocks.loci` | `Loci`, `LocusView`, `STRANDS`, `SCODE` | [loci]({{ '/api/loci/' | relative_url }}) |
| `genomeblocks.genes` | `Genes`, `GeneView`, `tss_base`, `LABELS`, `KINDS` | [genes]({{ '/api/genes/' | relative_url }}) |
| `genomeblocks.bedpe` | `Pairs`, `read_bedpe`, `as_pairs`, `read_pairs_chunks`, `PAIRS_FORMAT_COLUMNS`, `count_pairs`, `count_pairs_2d`, `pair_2d_block`, `pair_2d_to_frame` | [bedpe]({{ '/api/bedpe/' | relative_url }}) |
| `genomeblocks.architecture` | `Architecture`, `Props` | [architecture]({{ '/api/architecture/' | relative_url }}) |
| `genomeblocks.architecture_draw` | `draw` | [architecture]({{ '/api/architecture/' | relative_url }}) |
| `genomeblocks.signal` | `signal`, `tmm` | [signal]({{ '/api/signal/' | relative_url }}) |
| `genomeblocks.signal_draw` | `group_mask`, `plot_heatmap`, `plot_profiles`, `compare_heatmap` | [signal]({{ '/api/signal/' | relative_url }}) |
| `genomeblocks.motifs_draw` | `plot_motif_heatmap`, `plot_archetype`, `plot_archetypes`, `plot_cluster_members`, `plot_dendrogram` | [signal]({{ '/api/signal/' | relative_url }}) |
| `genomeblocks.browserview` | `browser` (re-exported as `gb.browser`) | [browser]({{ '/api/browser/' | relative_url }}) |
| `genomeblocks.bam` | `coverage`, `pileup_counts`, `reference_seq`, `BASES` | [browser]({{ '/api/browser/' | relative_url }}) |
| `genomeblocks.igv` | `igv_html`, `IGV_VERSION`, `IGV_CDN` | [browser]({{ '/api/browser/' | relative_url }}) |
| `genomeblocks.view` | `View`, `view_html` | [browser]({{ '/api/browser/' | relative_url }}) |
| `genomeblocks.motifs` | `scan_motifs`, `scan_motifs_matrix`, `scan_motifs_matrix_masked`, `scan_motifs_profile`, `bootstrap_enrichment`, `compare_motifs`, `compare_motifs_to_ref`, `pwm_distance_matrix`, `cluster_motifs`, `archetype`, `archetype_from_names`, `build_archetypes`, `write_meme` (re-exports `Library`, `load_motifs`) | [motifs]({{ '/api/motifs/' | relative_url }}) |
| `genomeblocks.atlas` | `Atlas` | [atlas]({{ '/api/atlas/' | relative_url }}) |
| `genomeblocks.se` | `knee`, `stitch_peaks`, `call_se`, `nearest_gene_within` | [atlas]({{ '/api/atlas/' | relative_url }}) |
| `genomeblocks.hichip` | `shortrange_ends`, `write_bed`, `fragments`, `coverage`, `to_bigwig`, `to_bedgraph`, `macs3`, `shortrange_track` | [atlas]({{ '/api/atlas/' | relative_url }}) |
| `genomeblocks.interop` | `as_loci`, `frame`, `loci_from_frame` / `loci_from_pyranges` / `loci_from_bedtool` / `loci_from_anndata`, `loci_to_pandas` / `loci_to_polars` / `loci_to_arrow` / `loci_to_bioframe` / `loci_to_pyranges` / `loci_to_bedtool` / `loci_to_cgranges` / `loci_to_anndata`, `cube_to_xarray` / `cube_to_anndata` / `cube_to_pandas`, `sequences`, `windows`, `revcomp`, `to_seqrecords`, `write_fasta`, `liftover` | [interop]({{ '/api/interop/' | relative_url }}) |
| `genomeblocks.backends` | `families`, `installed`, `resolve`, `use_backend`, `backends`, `unsupported`, `AUTO` | [backends]({{ '/api/backends/' | relative_url }}) |
| `genomeblocks.backends.intervals` | `overlap_pairs`, `overlap_any`, `nearest`, `merge`, `point_rows` | [backends]({{ '/api/backends/' | relative_url }}#backendsintervals) |
| `genomeblocks.backends.bigwig` | `open_bigwig`, `handle_backend`, `PyBigToolsHandle`, `PyBigWigHandle`, `PythonHandle` | [backends]({{ '/api/backends/' | relative_url }}#backendsbigwig) |
| `genomeblocks.backends.motifs` | `Library`, `load_motifs`, `logodds_matrix`, `threshold_from_pvalue`, `write_meme`, `Block`, `FORMATS`, `PSEUDOCOUNT` | [backends]({{ '/api/backends/' | relative_url }}#backendsmotifs) |
| `genomeblocks.backends.fasta` | `read_fasta`, `open_fasta`, `IndexedSource`, `PysamSource`, `PyfaidxSource`, `MemorySource`, `BiopythonSource` | [backends]({{ '/api/backends/' | relative_url }}#backendsfasta) |
| `genomeblocks.backends.tables` | `sniff`, `read_columns` | [backends]({{ '/api/backends/' | relative_url }}#backendstables) |
| `genomeblocks.backends.graph` | `to_scipy`, `native`, `components`, `pagerank`, `layout`, `layout_edges` | [backends]({{ '/api/backends/' | relative_url }}#backendsgraph) |

## Conventions

- **Coordinates** are 0-based, half-open in every table (`[start, end)`;
  a GTF's 1-based start becomes `start - 1` once, while parsing). A gene's
  TSS is the 1-bp interval `[t, t + 1)` with `t = start` on `+` and
  `t = end - 1` on `-`. See [Concepts]({{ '/concepts/' | relative_url }}).
- **Rows are the join key.** Anything computed from a `Loci` — a signal
  cube, a motif matrix, labels, `Architecture` vertex columns — has one row
  per locus in the table's order.
- **Inputs** go through `as_loci`: every function that takes intervals takes
  a `Loci`, a path, a region string, a `Locus`, a pandas / polars / pyarrow /
  bioframe / pyranges / pybedtools frame, a dict of columns, an AnnData or a
  list of regions.
- **Outputs**: every container has `shape`, `columns`, `head()`, `tail()`,
  `describe()`, `to_pandas()` / `to_polars()` / `to_arrow()`, `_repr_html_()`
  and the Arrow C stream, dataframe-interchange and narwhals protocols, so
  polars, pyarrow, DuckDB, seaborn, plotly and altair take it as is.
- **Backends**: `backend=` on a call or `with gb.use_backend(family=name):`
  picks the engine; a requested engine that is missing raises `ImportError`
  with the install command, never a silent switch. There is no global
  `Genome` and no global backend setting.
- **Signatures** on these pages are copied from the source; keyword-only
  arguments follow the `*`.
