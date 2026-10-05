#!/usr/bin/env bash
# Build and test the conda package from this checkout, with only conda-forge
# and bioconda enabled — what Bioconda's CI will see — before the release is
# on PyPI. The recipe is the release one with its source pointed here.
#
#   conda-recipe/build-local.sh [output-dir]       # needs conda-build on PATH
#
# Output defaults to ../genomeblocks-conda-build next to the checkout (it must
# not live inside the checkout: the whole tree is copied as the source).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/.." && pwd)"
out="${1:-$(dirname "$repo")/genomeblocks-conda-build}"
mkdir -p "$out/recipe"
sed -e "s|^  url: .*|  path: $repo|" -e "/^  sha256:/d" -e "/^  # filled from/d" \
    "$here/genomeblocks/meta.yaml" > "$out/recipe/meta.yaml"
cp "$here/genomeblocks/run_test.py" "$out/recipe/"
conda-build "$out/recipe" -c conda-forge -c bioconda --override-channels \
    --croot "$out/croot" --no-anaconda-upload "${@:2}"
echo "built: $(ls "$out"/croot/noarch/genomeblocks-*.conda "$out"/croot/noarch/genomeblocks-*.tar.bz2 2>/dev/null)"
