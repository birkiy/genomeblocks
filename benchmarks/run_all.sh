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
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python}

[[ -f data/hg38.chrom.sizes ]] || "$PY" make_data.py

BENCHES=("$@")
[[ ${#BENCHES[@]} -gt 0 ]] || BENCHES=(import loci loci_columnar signal atlas motifs pairs genes architecture shortrange prototype)

for b in "${BENCHES[@]}"; do
    echo "=== bench_$b ==="
    "$PY" -W ignore "bench_$b.py" 2>&1 | tee "results/$b.log"
done
"$PY" build_site.py           # the docs Benchmarks page
