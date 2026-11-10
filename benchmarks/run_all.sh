#!/usr/bin/env bash
# Run the whole benchmark suite, one bench at a time (several use every core,
# so running them concurrently would skew each other's timings).
#
#   PY=/path/to/python ./run_all.sh            # data + every bench + the docs page
#   PY=/path/to/python ./run_all.sh loci atlas  # just these benches
#
# External baselines are looked up on PATH (or set the variables):
#   bedtools, computeMatrix (deepTools), cooler, bgzip,
#   GIGGLE=/path/to/giggle  (github.com/ryanlayer/giggle, built from source)
#   FIMO=/path/to/fimo      (MEME suite: micromamba create -c bioconda meme)
#   GB_RSCRIPT=/path/Rscript:/other/Rscript  (R with tximport / DESeq2 / edgeR, for bench_expression)
#   FITHICHIP_HOME=/path/to/FitHiChIP [FITHICHIP_EXEC="singularity exec img.sif"]  (bench_loops; else gb only)
#   GB_BENCH_LOOPS_REAL=samples.tsv  (bench_loops on real samples against existing FitHiChIP runs)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python}

[[ -f data/hg38.chrom.sizes ]] || "$PY" make_data.py

BENCHES=("$@")
[[ ${#BENCHES[@]} -gt 0 ]] || BENCHES=(import backends loci signal atlas motifs pairs genes architecture shortrange loops expression)

for b in "${BENCHES[@]}"; do
    echo "=== bench_$b ==="
    "$PY" -W ignore "bench_$b.py" 2>&1 | tee "results/$b.log"
done
"$PY" build_site.py           # the docs Benchmarks page
