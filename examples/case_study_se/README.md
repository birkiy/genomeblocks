# Case study: HiChIP + ATAC → super-enhancers, prime hubs, shareable view

Runs the whole recipe in the columnar prototype (`genomeblocks.columnar`):

```
chrom sizes ─▶ Genome
ATAC peaks ─▶ CREs        GTF/GFF3 + ATAC bw ─▶ Genes (isoform-fixed)
HiChIP bedpe + mcool ─▶ Architecture ─▶ O/E · genes · prime hubs
HiChIP allValidPairs ─▶ short-range ends ─▶ MACS3 peaks + coverage bigWig ─▶ SEs
                     ─▶ shared / SE-only / prime-only genes ─▶ View (one HTML)
```

| file | what |
| --- | --- |
| `case_study_se.ipynb` | executed notebook (synthetic hg38-shaped stand-ins); first code cell has the mESC/mm10 paths |
| `make_notebook.py` | rebuilds and re-executes it |

Inputs: ATAC peaks + bigWig, H3K27ac HiChIP loops (BEDPE), .mcool, HiC-Pro
allValidPairs, GTF or GFF3. Needs MACS3 for the short-range peaks. Writes
`work/` (short-range bigWig, peaks, `mesc_view.html`).

The short-range step (`hichip.shortrange_ends` → `fragments` → `to_bigwig`)
gives the same ends and a byte-identical bedGraph as the awk | sort |
`bedtools genomecov` recipe (`benchmarks/bench_shortrange.py`).
