# Conda packaging (Bioconda)

`genomeblocks/meta.yaml` is a ready Bioconda recipe (`noarch: python`), and
`genomeblocks/run_test.py` is its smoke test. Bioconda is the right channel:
the bioinformatics dependencies (`cooler`, `pybigtools`, `moods`, and the
optional `pysam`, `cgranges`, `bedtools`, ...) live there, everything else on
conda-forge.

## Dependencies (checked 2026-10-05)

| run dependency | channel |
| --- | --- |
| numpy, pandas >=2, pyarrow, narwhals, scipy, matplotlib-base, tqdm | conda-forge |
| cooler, pybigtools | bioconda |
| **moods** — the motif engine | bioconda |

**No lightmotif.** lightmotif is on PyPI only, so the conda package scans
motifs with MOODS (C++), which reports the same hits: genomeblocks feeds every
engine the same log-odds matrices, parses motif files itself (JASPAR,
jaspar16, TRANSFAC, uniprobe, MEME) and computes p-value cutoffs itself, so
nothing needs lightmotif. In a conda environment `pip install lightmotif`
adds it back as an alternative engine (`backend="lightmotif"`).

Optional engines users can add from conda: `graph-tool` (conda-forge),
`polars`, `pysam`, `pyfaidx`, `pybigwig`, `cgranges`, `ncls`, `bedtools` +
`pybedtools`, `bioframe`, `pyranges`, `anndata`, `xarray`, `biopython`,
`networkx`, `python-igraph`, `logomaker`.

MOODS is GPL-3.0-or-later (or the Biopython licence); genomeblocks stays MIT
and only depends on it at run time.

## Verified locally

```bash
conda-recipe/build-local.sh /path/to/output     # needs conda-build on PATH
```

builds the recipe from this checkout (the release recipe with its source
pointed here) using only `conda-forge` and `bioconda`, the way Bioconda's CI
will. The package test imports genomeblocks, checks that MOODS is the motif
engine with no lightmotif present, runs interval set algebra and scans a
planted motif with a score and with a p-value cutoff. The full test suite
also passes against the installed package in a conda environment with every
optional engine (Python 3.12, pandas 3.0, numpy 2.5).

## From a release to Bioconda

1. **Release on PyPI** (PyPI has 1.0.0 and 1.0.1 only; 1.1.0 was never
   uploaded):
   ```bash
   python -m build              # sdist + wheel in dist/
   twine upload dist/*
   ```
2. **Fill the sha256** of the sdist in `genomeblocks/meta.yaml`:
   ```bash
   curl -sL https://pypi.org/packages/source/g/genomeblocks/genomeblocks-2.0.0.tar.gz | sha256sum
   ```
3. **Submit**: fork [bioconda-recipes](https://github.com/bioconda/bioconda-recipes),
   copy `genomeblocks/` (meta.yaml and run_test.py) to `recipes/genomeblocks/`,
   open a pull request. Bioconda's CI lints, builds and tests it; a Bioconda
   member reviews and merges. After that, BiocondaBot opens the version-bump
   pull requests itself whenever a new release lands on PyPI.
4. **Install**:
   ```bash
   conda install -c conda-forge -c bioconda genomeblocks     # or mamba / micromamba
   ```

## pip, for comparison

`pip install genomeblocks` has no motif engine (MOODS has no wheels on PyPI,
it compiles): `pip install "genomeblocks[motifs]"` builds MOODS-python, and
`pip install "genomeblocks[lightmotif]"` takes lightmotif's prebuilt wheels.
Both report the same hits.
