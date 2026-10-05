# genomeblocks examples

Each example is a self-contained directory: a `download_data.sh` to fetch its
data and one focused notebook that runs a single coherent analysis end to end.

| example | what it shows |
| --- | --- |
| [`case_study_se/`](case_study_se/) | Case study: H3K27ac HiChIP (BEDPE, .mcool, allValidPairs) + ATAC + GTF/GFF3 → isoform-fixed genes, Architecture, HiChIP short-range track + MACS3, ROSE super-enhancers, prime vs SE genes, and a shareable genomeblocks view. Runs on the synthetic benchmark data (`benchmarks/make_data.py`). |
| [`ar_foxa1_lncap/`](ar_foxa1_lncap/) | AR / FOXA1 cistrome dependence in LNCaP (±DHT): accessibility, set logic, signal heatmaps, annotation, motif + ChIP-Atlas enrichment, browser view. Covers `Loci` set ops, `signal`/`plot_heatmap`, `Genes.annotations`, `scan_motifs_matrix`/`bootstrap_enrichment`, `Atlas`, and `browser`. |
