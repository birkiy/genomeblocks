#!/usr/bin/env bash
# Run the whole benchmark suite, one bench at a time (several use every core,
# so running them concurrently would skew each other's timings).
#
#   PY=/path/to/python ./run_all.sh            # data + every bench + figures
#   PY=/path/to/python ./run_all.sh loci atlas  # just these benches
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python}

[[ -f data/hg38.chrom.sizes ]] || "$PY" make_data.py

BENCHES=("$@")
[[ ${#BENCHES[@]} -gt 0 ]] || BENCHES=(import loci signal atlas motifs pairs genes architecture fixes)

for b in "${BENCHES[@]}"; do
    echo "=== bench_$b ==="
    "$PY" -W ignore "bench_$b.py" 2>&1 | tee "results/$b.log"
done
"$PY" plot.py
