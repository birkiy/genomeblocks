# Conda packaging

Staging recipe for [Bioconda](https://bioconda.github.io). genomeblocks is a
bioinformatics package (its deps — `cooler`, `pyranges`, `pybigtools`,
`cgranges`, `graph-tool` — live on Bioconda / conda-forge), so **Bioconda** is
the right channel.

```
conda-recipe/
  genomeblocks/meta.yaml   # the package recipe (PyPI source)
```

## Dependency status

Bioconda can only depend on packages already on conda channels. Checked:

| dependency | channel | status |
| --- | --- | --- |
| graph-tool | conda-forge | ✅ available |
| cooler, pyranges, pybigtools, cgranges | bioconda | ✅ available |
| numpy, pandas, scipy, matplotlib-base, tqdm | conda-forge | ✅ available |
| **lightmotif** | — | ❌ PyPI-only → recipe lives at `~/apps/lightmotif-conda`, submit first |

~~conorm~~ was **dropped** in v1: it is unlicensed (un-packageable) and was only
used by `tmm()`, so the TMM algorithm is now vendored in `genomeblocks/signal.py`
(`_tmm_norm_factors`, numerically identical to conorm). No dependency remains.

### lightmotif

`lightmotif` (used by motif scanning, lazy-imported) is a compiled Rust/PyO3
package built with `maturin`. Its conda recipe is kept **outside this repo** at
`~/apps/lightmotif-conda/` — submit it to Bioconda (or conda-forge) first so the
dependency resolves. Alternatively make it optional (drop from `run:` here and
document it as a pip extra) since `scan_motifs*` only need it on demand.

## Publishing checklist

The recipe sources the **PyPI sdist** (Bioconda's CI fetches it), so:

1. **Publish to PyPI** (genomeblocks is not there yet):
   ```bash
   python -m build --sdist --wheel
   twine upload dist/*
   ```
2. **Fill the sha256** in `genomeblocks/meta.yaml`:
   ```bash
   openssl dgst -sha256 dist/genomeblocks-1.0.0.tar.gz
   # or, from the live PyPI release:
   #   curl -sL <pypi-sdist-url> | openssl dgst -sha256
   ```
3. **Submit to Bioconda**: fork
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes), copy
   `genomeblocks/` here into `recipes/genomeblocks/`, and open a PR. Submit
   `lightmotif` (from `~/apps/lightmotif-conda`) first/with it.

## Local sanity check

```bash
conda install -n base -c conda-forge conda-build
conda build conda-recipe/genomeblocks -c conda-forge -c bioconda
```
