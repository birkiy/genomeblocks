# AR & FOXA1 in LNCaP (±DHT)

A single, self-contained genomeblocks analysis: how the androgen receptor (AR)
cistrome depends on the pioneer factor FOXA1 after androgen (DHT) stimulation in
LNCaP prostate cancer cells.

## The workflow ([`ar_foxa1_lncap.ipynb`](ar_foxa1_lncap.ipynb))

1. **Accessible chromatin** — union of ATAC-seq peaks (0h + 4h).
2. **Inside vs outside** — fraction of AR / FOXA1 binding in accessible chromatin.
3. **Filter** to accessible peaks.
4. **Venn** — FOXA1 (0h, 4h) vs AR (4h), defining:
   - **AR+F** — AR 4h peaks that *are* FOXA1-bound (FOXA1-dependent AR)
   - **AR−F** — AR 4h peaks that are *not* FOXA1-bound (FOXA1-independent AR)
5. **Heatmaps** of AR+F / AR−F across all 6 conditions (ATAC reps grouped).
6. **Genomic annotation** pie charts (AR+F vs AR−F).
7. **Motif enrichment** — AR+F vs AR−F, ATAC union as background pool.
8. **ChIP-Atlas (atlas) enrichment** — which TFs differ between AR+F and AR−F.
9. **Browser view** of `chr19:50,792,009-50,923,669`.

## Data

```bash
./download_data.sh        # → ./data/*.bed, ./data/*.bw  (ChIP-Atlas, hg38)
```

The ChIP-Atlas accessions are documented at the top of `download_data.sh`.

### hg38 reference files (edit paths at the top of the notebook)

| what | default path |
| --- | --- |
| genome FASTA | `/groups/lackgrp/genome_annotations/hg38/hg38.fa` |
| gene GTF | `/groups/lackgrp/genome_annotations/hg38/gencode.v49.annotation_protein_coding.gtf` |
| JASPAR motifs | `/groups/lackgrp/databases/motifs/motif-db/H14CORE_jaspar_format_lightmotif.txt` |
| chrom sizes | `/groups/lackgrp/genome_annotations/hg38/hg38_chr.chrom.sizes` |
| ChIP-Atlas meta | `/groups/lackgrp/databases/giggle/giggle_hg38/hg38_tfs_meta.tsv` |
| ChIP-Atlas index dir | `/groups/lackgrp/databases/giggle/giggle_hg38/ChIP-Atlas-ALL` |
