# genomeblocks examples

Each example is a self-contained directory: a `download_data.sh` to fetch its
data and one focused notebook that runs a single coherent analysis end to end.

| example | what it shows |
| --- | --- |
| [`columnar_prototype/`](columnar_prototype/) | Prototype tour of `genomeblocks.columnar` (branch `columnar-prototype`): CREs, genes and the Architecture as row-aligned tables; per-chromosome / cis / trans views of one graph; graph-tool algorithms; save/reload; hand-off to pandas, polars, Arrow. Runs on the synthetic benchmark data. |
| [`ar_foxa1_lncap/`](ar_foxa1_lncap/) | AR / FOXA1 cistrome dependence in LNCaP (±DHT): accessibility, set logic, signal heatmaps, annotation, motif + ChIP-Atlas enrichment, browser view. Covers `Loci` set ops, `signal`/`plot_heatmap`, `Genes.annotations`, `scan_motifs_matrix`/`bootstrap_enrichment`, `Atlas`, and `browser`. |
