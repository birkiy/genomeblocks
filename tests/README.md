# Tests

Fast, synthetic tests for genomeblocks 2.0 — no external genome, bigWig or
ChIP-Atlas data. Every fixture is a tiny file written into a temporary
directory or a small in-memory table (see [`conftest.py`](conftest.py)).

## Running

```bash
pip install -e ".[test]"      # pytest and every pip-installable backend
pytest                         # from the repo root
```

The parity tests run over the backends installed in the environment
(`conftest.installed_backends(family)`); graph-tool, cgranges and bedtools
are conda-only and skipped when absent. The default engines need nothing
beyond `pip install genomeblocks`.

## Coverage

| file | what |
| --- | --- |
| `test_api.py` | public surface, lazy exports, import cost, the backends registry, `use_backend` |
| `test_locus.py` | `Locus`, `parse_region` / `parse_regions` (units, split views) |
| `test_loci.py` | `Loci.make` (headers, kept columns, bad files), set algebra vs brute force on every interval backend, zero-length rule, merge / nearest, table ergonomics, protocols (polars, Arrow, duckdb, narwhals), round trips (pandas, polars, Arrow, bioframe, pyranges, parquet, BED, AnnData), sequences |
| `test_interop.py` | `as_loci` on every input, lenient column names, CSV / TSV readers, file order, nulls, cube converters, liftover |
| `test_genes.py` | 0-based GTF / GFF3 / genePred tables, both parsers agree, Ensembl spellings, duplicate gene ids, UCSC layouts, `from_frame`, annotation on every backend, `select_isoforms`, `to_gtf` / `to_bed12` / save / load, table protocols |
| `test_bedpe.py` | `Pairs` (headers, frames, filters, overlaps, protocols), pairs-file format detection, `count_pairs` vs brute force |
| `test_architecture.py` | `make` on every interval backend, normalize / prune / annotate / strength / hubs / support, views, graph backends agree, exports, `add_mcool`, `draw` |
| `test_signal.py` | the three bigWig backends bin alike (every stat, fractional bins, chromosome edges), `signal()` inputs and handles, multiprocess parity, `tmm` |
| `test_signal_draw.py` | `plot_heatmap` / `plot_profiles` / `compare_heatmap` / `plot_motif_heatmap` headless |
| `test_motifs.py` | motif formats, `Library` as a table, every scanning engine finds the planted sites, masked / profile scans, statistics, FASTA backends agree |
| `test_atlas.py` | `Atlas.make` inputs, columnar bin ranges, search / bootstrap on any interval input, metadata, save / load, table protocols |
| `test_browser.py` | `browser()` with every track kind (incl. BAM), `igv_html`, `View` payload (0-based genes, TSS) |
| `test_bam.py` | BAM pileups, `se.stitch_peaks` / `call_se`, HiChIP short-range tracks |
| `test_columnar.py` | the table contract on every container; seaborn / plotly / altair take a Loci; no hidden global state; backend differences do not leak |
| `test_functional.py` | end to end: peaks → genes → architecture → signal → motifs → enrichment → views |
