# genomeblocks benchmarks

How fast each building block is against the tools people would otherwise use,
and what each swappable backend costs. The current results, with charts and
every table, are on the docs site:
**[birkiy.github.io/genomeblocks/benchmarks](https://birkiy.github.io/genomeblocks/benchmarks/)**.
That page is generated from `results/*.json` by `build_site.py`.

Every comparison first checks that the engines return the same answer
(same intervals, same matrix, same hit counts, same components), then times
them: medians of 3–5 runs after a warm-up, one bench at a time.

## Reproduce

```bash
cd benchmarks
python make_data.py            # synthetic hg38-shaped data, ~1.5 GB, seeded
PY=python ./run_all.sh         # every bench, then build_site.py
PY=python ./run_all.sh backends loci     # just these
```

`make_data.py` needs `bioframe` (chromosome sizes) and `pyjaspar` (JASPAR
CORE motifs); nothing is downloaded at run time. External baselines are found
on `PATH`: bedtools, deepTools `computeMatrix`, cooler. Optional:
`FIMO=/path/to/fimo` (MEME suite) and `GIGGLE=/path/to/giggle` (built from
github.com/ryanlayer/giggle; skipped when missing). Backends that are not
installed are skipped and listed in the results (`graph-tool`, `cgranges`
and `bedtools` are conda-only). `GB_BENCH_SIZES=1000,10000` limits the
interval sizes for a quick run.

## Files

| file | what |
|---|---|
| `make_data.py` | peaks (1k–1M), genome-wide bigWigs, 500 Atlas tracks, a GENCODE-shaped GTF, 5M Hi-C pairs, loops (cis + trans), coolers, HiChIP pairs, JASPAR motifs |
| `bench_backends.py` | **every backend of every family through the genomeblocks API**: intervals (A & B, merge, nearest), bigWig signal, motif scanning, FASTA fetching, table parsing, graph algorithms — with the agreement check |
| `bench_loci.py` | intersect, merge, BED parsing, single-query latency, index caching vs pyranges, bioframe, bedtools, intervaltree |
| `bench_signal.py` | bigWig → matrix engines, scaling, processes vs threads vs deepTools |
| `bench_atlas.py` | Atlas build, query, load, accuracy, bootstrap (+ GIGGLE when available) |
| `bench_motifs.py` | genomeblocks vs FIMO, MOODS, Biopython, numpy; library size; workers |
| `bench_pairs.py` | `count_pairs` / `count_pairs_2d` vs cooler cload and a per-pair loop |
| `bench_genes.py` | GTF parsing, annotation index, labelling, nearest gene, isoform choice |
| `bench_architecture.py` | the Architecture pipeline step by step, graph backends, save / load |
| `bench_shortrange.py` | HiChIP short-range track vs the awk / sort / bedtools recipe |
| `bench_import.py` | import cost of each entry point |
| `build_site.py` | writes `docs/benchmarks/index.md` |

## Design studies (historical)

`archive/` holds the design studies measured on 1.0.1 and 1.1 that led to the
columnar 2.0 core: `bench_prototype.py` (columnar vs classic), the backend and
column prototypes (`prototypes/`), their reports (`report/`, `figures/`) and
builders (`build_report.py`, `build_proposal.py`, `build_explainer.py`,
`build_prototype_report.py`, `plot.py`), plus `bench_fixes.py`,
`bench_columnar.py`, `bench_genes_make.py`, `bench_motif_backends.py` and
`bench_loci_columnar.py` with their results. They import the 1.x API and are
kept for the record, not to run.
