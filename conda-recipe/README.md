# Conda packaging

Staging recipes for [Bioconda](https://bioconda.github.io). genomeblocks is a
bioinformatics package (its deps — `cooler`, `pyranges`, `pybigtools`,
`cgranges`, `graph-tool` — live on Bioconda / conda-forge), so **Bioconda** is
the right channel, not conda-forge.

```
conda-recipe/
  genomeblocks/meta.yaml   # the package recipe (PyPI source)
  conorm/meta.yaml         # a missing dependency, packaged to unblock
```

## Dependency status

Bioconda can only depend on packages already on conda channels. Checked:

| dependency | channel | status |
| --- | --- | --- |
| graph-tool | conda-forge | ✅ available |
| cooler, pyranges, pybigtools, cgranges | bioconda | ✅ available |
| numpy, pandas, scipy, matplotlib-base, tqdm | conda-forge | ✅ available |
| **conorm** | — | ❌ PyPI-only → recipe in `conorm/` (submit first) |
| **lightmotif** | — | ❌ PyPI-only → **blocker**, see below |

### lightmotif

`lightmotif` (used by the motif-scanning functions) is a compiled (Rust/PyO3)
package and is not on any conda channel. Options, in order of preference:

1. **Package it for Bioconda/conda-forge** (it has wheels, but conda needs a
   build recipe with the Rust toolchain). This is the clean fix.
2. **Make it an optional runtime dependency** — `scan_motifs*` already
   lazy-imports it, so the conda package works for everything except motif
   scanning until lightmotif is available. If you take this route, drop it from
   `run:` in `genomeblocks/meta.yaml` and note it as a pip extra.

## Publishing checklist

The recipe sources the **PyPI sdist** (Bioconda's CI fetches it), so:

1. **Publish to PyPI** (genomeblocks is not there yet):
   ```bash
   python -m build --sdist --wheel
   twine upload dist/*
   ```
2. **Fill the sha256** in `genomeblocks/meta.yaml` (and `conorm/meta.yaml`):
   ```bash
   openssl dgst -sha256 dist/genomeblocks-1.0.0.tar.gz
   # or, from the live PyPI release:
   #   curl -sL <pypi-sdist-url> | openssl dgst -sha256
   ```
3. **Submit to Bioconda**: fork
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes), copy each
   directory here into `recipes/<name>/`, and open a PR. Submit `conorm` (and
   `lightmotif`, if packaging it) before/with `genomeblocks` so the dependency
   graph resolves.

## Local sanity check

```bash
conda install -n base -c conda-forge conda-build
conda build conda-recipe/conorm        # build the dep first
conda build conda-recipe/genomeblocks -c conda-forge -c bioconda
```
