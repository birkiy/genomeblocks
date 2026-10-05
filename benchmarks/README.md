# genomeblocks benchmarks

How fast each building block is against the tools people would otherwise use.
The current results, with charts and every table, are on the docs site:
**[birkiy.github.io/genomeblocks/benchmarks](https://birkiy.github.io/genomeblocks/benchmarks/)**.
That page is generated from `results/*.json` by `build_site.py`.

Every comparison first checks that the engines return the same answer
(same intervals, same matrix, same hit counts, same graph), then times them:
medians of 3–5 runs after a warm-up, one bench at a time.

## Reproduce

```bash
cd benchmarks
python make_data.py            # synthetic hg38-shaped data, ~1.5 GB, seeded
PY=python ./run_all.sh         # every bench, then build_site.py
```

`make_data.py` needs `bioframe` (chromosome sizes) and `pyjaspar` (JASPAR
CORE motifs); nothing is downloaded at run time. External baselines are found
on `PATH`: bedtools, deepTools `computeMatrix`, cooler. Optional:
`FIMO=/path/to/fimo` (MEME suite) and `GIGGLE=/path/to/giggle` (built from
github.com/ryanlayer/giggle; skipped when missing). `Architecture` needs
graph-tool (conda-forge).

## Files

| file | what |
|---|---|
| `make_data.py` | peaks (1k–1M), genome-wide bigWigs, 500 Atlas tracks, a GENCODE-shaped GTF, 5M Hi-C pairs, loops (cis + trans), coolers, HiChIP pairs, JASPAR motifs |
| `bench_loci.py`, `bench_loci_columnar.py` | intersect, merge, BED parsing, single-query latency, index caching |
| `bench_signal.py` | bigWig → matrix engines, scaling, processes vs threads vs deepTools |
| `bench_atlas.py` | Atlas build, query, load, accuracy, bootstrap (+ GIGGLE when available) |
| `bench_motifs.py` | genomeblocks vs FIMO, MOODS, Biopython, numpy; library size; workers |
| `bench_pairs.py` | `count_pairs` / `count_pairs_2d` vs cooler cload and a per-pair loop |
| `bench_genes.py` | GTF parsing, annotation index, labelling, nearest gene, isoform choice |
| `bench_architecture.py` | the classic Architecture pipeline step by step |
| `bench_prototype.py` | `genomeblocks.columnar` vs the classic modules: steps, views, trans loops, scaling, fresh-process end-to-end |
| `bench_shortrange.py` | HiChIP short-range track vs the awk / sort / bedtools recipe |
| `bench_import.py` | import cost of each entry point |
| `build_site.py` | writes `docs/benchmarks/index.md` |

## Design reports (historical)

`report/*.html`, `figures/`, `prototypes/` and the `build_report.py`,
`build_proposal.py`, `build_explainer.py`, `build_prototype_report.py`,
`plot.py` scripts are the design studies measured on 1.0.1. They found the
slow paths fixed in 1.1 and led to `genomeblocks.columnar`. `bench_fixes.py`,
`bench_backends.py`, `bench_columnar.py`, `bench_genes_make.py` and
`bench_motif_backends.py` belong to those studies.
