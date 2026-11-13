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
interval sizes for a quick run. `bench_expression.py` runs R through
`GB_RSCRIPT=/path/Rscript:/other/Rscript` (the first Rscript with the step's
packages is used; steps without one are skipped and listed) and
`GB_BENCH_EXPR=20000x6,60000x24` picks its sizes. `bench_loops.py` runs
FitHiChIP when `FITHICHIP_HOME` names its directory (`FITHICHIP_EXEC`, e.g.
`singularity exec fithichip.sif`, prefixes the call; without it genomeblocks is
timed alone) and compares real samples with existing FitHiChIP runs listed in
`GB_BENCH_LOOPS_REAL=samples.tsv`. `bench_maps.py` runs MAPS when `MAPS_HOME` points at a clone of
github.com/ijuric/MAPS (`MAPS_PYTHON` / `MAPS_RSCRIPT`: a Python with pandas and an Rscript with VGAM, MASS,
data.table); without it genomeblocks alone is timed on a made-up features table. `GB_BENCH_MAPS_REAL=samples.tsv`
compares real samples with existing MAPS runs. `bench_hiccups.py` runs Juicer when `JUICER_TOOLS` names
juicer_tools_2.13.06.jar (`JAVA` for the java), with `GB_BENCH_THREADS` threads / workers; real samples (Juicer
output + dumped vectors) via `GB_BENCH_HICCUPS_REAL=samples.tsv`. `bench_mustache.py` runs Mustache when `MUSTACHE` names its
command (`MUSTACHE_EXEC` prefixes it; `GB_BENCH_THREADS`, default 4, sets both tools' cores) and compares real
samples with existing Mustache runs listed in `GB_BENCH_MUSTACHE_REAL=samples.tsv` (name, cool, mustache_tsv,
mustache_seconds, note).

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
| `bench_loops.py` | HiChIP loop calling: `fithichip` vs FitHiChIP 11.0 on the same cooler and peaks — tested pairs, loops, merged loops, P / q values and FitHiChIP's IGV files compared, wall time and memory; synthetic cooler, plus real samples on request |
| `bench_mcool_split.py` | `add_mcool`'s pixel split: each sample's own graph against a union graph under eight sharing rules (equal split, no split, each pixel once per CRE, centre binning, ATAC-weighted, peak mask, with and without the self-contact pixel) — super-hubs lost, CREs diluted, contact conservation; seeded synthetic cooler, plus the T-cell HiChIP numbers recorded as measured |
| `bench_maps.py` | HiChIP loop calling: `maps` vs MAPS (ijuric/MAPS) on the same pairs and peaks — tested bin pairs, counts, expected / P / fdr, calls and MAPS's sig3Dinteractions.bedpe compared, wall time and memory; synthetic HiChIP pairs with loops, plus real samples on request |
| `bench_mustache.py` | Hi-C / HiChIP loop calling: `mustache()` vs Mustache 1.3.3 on the same balanced cooler — loop sets, FDR / detection scale, Mustache's output file, wall time and memory with the same cores; synthetic cooler (built on first run), plus real samples on request |
| `bench_hiccups.py` | Hi-C loop calling: `hiccups` vs juicer_tools 2.13.06 `hiccups --cpu` on the same matrix and SCALE vectors — merged loops, postprocessed / enriched pixels, thresholds compared file by file, wall time and memory; synthetic Hi-C with planted loops, plus real samples on request |
| `bench_import.py` | import cost of each entry point |
| `bench_expression.py` | reading salmon / featureCounts tables with both table backends, tximport-style summarising and TMM / size factors / CPM / filterByExpr against R's tximport, edgeR and DESeq2 — exactness (max \|Δ\|) and time; generates its own GTF-scale data |
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
