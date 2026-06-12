#!/usr/bin/env bash
# Download the ChIP-Atlas data for the AR / FOXA1 LNCaP example into ./data
#
# LNCaP, ±4h DHT (dihydrotestosterone). SRX accessions (ChIP-Atlas, hg38):
#
#   SRX23002839  FOXA1  0h   (peak + bigwig)
#   SRX23002841  FOXA1  4h   (peak + bigwig)
#   SRX23002834  AR     0h   (peak + bigwig)
#   SRX23002836  AR     4h   (peak + bigwig)
#   SRX23002894  ATAC   0h rep1   (peak + bigwig)
#   SRX23002895  ATAC   0h rep2   (peak + bigwig)
#   SRX23002898  ATAC   4h rep1   (peak + bigwig)
#   SRX23002899  ATAC   4h rep2   (peak + bigwig)
#
# Peaks are bed05 (q < 1e-5). Re-running skips files already present.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data

BASE="https://chip-atlas.dbcls.jp/data/hg38/eachData"
SRX=(SRX23002839 SRX23002841 SRX23002834 SRX23002836 \
     SRX23002894 SRX23002895 SRX23002898 SRX23002899)

get() {  # url -> dest (skip if non-empty)
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then echo "  have $(basename "$dest")"; return; fi
    echo "  get  $(basename "$dest")"
    curl -fL --retry 3 --retry-delay 2 -o "$dest" "$url"
}

echo "[1/2] peak files (bed05)"
for s in "${SRX[@]}"; do get "$BASE/bed05/${s}.05.bed" "data/${s}.05.bed"; done

echo "[2/2] signal files (bigwig)"
for s in "${SRX[@]}"; do get "$BASE/bw/${s}.bw" "data/${s}.bw"; done

echo "done -> $(pwd)/data"
