---
title: Installation
layout: default
nav_order: 2
---

# Installation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Requirements

- **Python** ≥ 3.10
- Everything else is declared in `pyproject.toml` and installed by pip.

`pip install genomeblocks` brings the **default backend of every family except motif scanning** — the motif engines are compiled, so they are extras (below) — and everything else works after a plain pip install:

| Family | Default engine | Pulled in by pip |
|---|---|---|
| intervals (overlap, nearest, merge) | genomeblocks' numpy kernel | — (built in) |
| bigWig reading | [pybigtools](https://github.com/jackh726/bigtools) | `pybigtools` |
| motif scanning | [MOODS](https://github.com/jhkorhonen/MOODS), else [lightmotif](https://github.com/althonos/lightmotif) | — (an extra: `genomeblocks[motifs]` or `genomeblocks[lightmotif]`; conda: `moods`) |
| FASTA | genomeblocks' indexed (`.fai`) reader | — (built in) |
| text tables (BED, GTF, pairs) | pandas (polars when installed) | `pandas`, `pyarrow` |
| graph algorithms | scipy (graph-tool when installed) | `scipy` |
| Hi-C `.cool` / `.mcool` | [cooler](https://github.com/open2c/cooler) | `cooler` |
| the boundary and the protocols | [narwhals](https://github.com/narwhals-dev/narwhals), [pyarrow](https://arrow.apache.org/) | `narwhals`, `pyarrow` |

TMM normalisation (`gb.tmm`) is vendored — the edgeR algorithm ships inside `genomeblocks.signal`. Heavy modules are imported lazily, so `import genomeblocks` itself takes a few milliseconds. See the [Credits page]({{ '/credits/' | relative_url }}) for the citation of every upstream tool.

---

## pip

```bash
pip install genomeblocks
```

That is the default install: everything but a motif engine. Optional extras add the motif engines, other backends and converter targets:

| Extra | Adds | Use it for |
|---|---|---|
| `motifs` | `MOODS-python` | the default motif engine, [MOODS](https://github.com/jhkorhonen/MOODS) (C++, compiled at install) |
| `lightmotif` | `lightmotif` | the alternative motif engine (prebuilt wheels, no compiler); same hits; used automatically only when MOODS is absent |
| `fast` | `polars` | faster GTF / GFF3 / BED / pairs parsing (same columns as pandas) |
| `bam` | `pysam` | BAM pileups in `browser` / `gb.coverage`, FASTA through pysam |
| `viz` | `logomaker` | motif logos (`genomeblocks.motifs_draw`) |
| `interop` | `polars`, `bioframe`, `pyranges`, `pybedtools`, `anndata`, `xarray`, `biopython`, `networkx`, `igraph`, `MOODS-python`, `pyBigWig`, `ncls`, `pyfaidx`, `pyliftover` | every other backend and every `from_*` / `to_*` target |
| `all` | `motifs` + `fast` + `bam` + `viz` + `interop` | everything that pip can install |
| `test` | the above plus `lightmotif`, `pytest`, `duckdb`, `seaborn`, `plotly`, `altair` | running the test suite |

```bash
pip install "genomeblocks[motifs]"     # + MOODS, the motif engine (compiled at install)
pip install "genomeblocks[lightmotif]" # or lightmotif (prebuilt wheels); same hits
pip install "genomeblocks[all]"        # every pip-installable backend
pip install "genomeblocks[fast,bam]"   # pick the ones you need
```

{: .note }
> Apart from `motifs` / `lightmotif`, extras only add *choices*: the default engine of every other family is already there, and results are identical across engines (the test suite checks it), so install an extra for speed on your data or to hand results to a library you already use. The motif functions are the one exception — without an engine they raise `ImportError` naming both extras (see below).

---

## conda: the three conda-only engines

Three engines cannot come from PyPI. They are optional backends; `genomeblocks` never needs them:

| Engine | Family | Install | What it adds |
|---|---|---|---|
| [graph-tool](https://graph-tool.skewed.de/) | graph | `conda install -c conda-forge graph-tool` | the default graph engine when installed (scipy otherwise; identical components, centralities within tolerance); `A.to_graph_tool()` |
| [cgranges](https://github.com/lh3/cgranges) | intervals | `conda install -c bioconda cgranges` (or `pip install git+https://github.com/lh3/cgranges`) | a C interval index; `L.to_cgranges()` |
| [bedtools](https://bedtools.readthedocs.io/) | intervals | `conda install -c bioconda bedtools` + `pip install pybedtools` | the `bedtools` backend (the binary must be on `PATH`); `L.to_bedtool()` |

The conda package of `genomeblocks` (Bioconda) depends on `moods`, so a conda install has the motif engine from the start; `conda install -c bioconda moods` adds it to a pip environment too.

The repository ships an `environment.yml` with every backend, including these three:

```bash
git clone https://github.com/birkiy/genomeblocks.git
cd genomeblocks
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate genomeblocks
```

{: .tip }
> `Architecture` works on a plain pip install. Graph algorithms run on scipy when graph-tool is absent, and `gb.backends()` tells you which engine is in use.

---

## Verifying the install

`gb.backends()` is the one check. It returns a DataFrame with every family and backend, whether it is installed, which one is the automatic default here, which one a call would use right now, and the install command for the missing ones:

```python
import genomeblocks as gb
gb.__version__
# -> '2.0.0'

gb.backends()
#        family       backend  installed  default  in use   install
# 0   intervals  genomeblocks       True     True    True
# 1   intervals      cgranges      False    False   False   conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)
# 2   intervals          ncls       True    False   False
# 3   intervals      bioframe       True    False   False
# 4   intervals      pyranges       True    False   False
# 5   intervals      bedtools      False    False   False   pip install pybedtools; conda install -c bioconda bedtools  (the bedtools binary must be on PATH)
# 6      bigwig    pybigtools       True     True    True
# ...
# 19      graph    graph-tool      False    False   False   conda install -c conda-forge graph-tool
# 20      graph         scipy       True     True    True
# 21      graph        igraph       True    False   False
# 22      graph      networkx       True    False   False

gb.backends.families()
# -> {'intervals': ['genomeblocks', 'cgranges', 'ncls', 'bioframe', 'pyranges', 'bedtools'],
#     'bigwig': ['pybigtools', 'pybigwig', 'python'],
#     'motifs': ['moods', 'lightmotif', 'biopython'],
#     'fasta': ['genomeblocks', 'pysam', 'pyfaidx', 'memory', 'biopython'],
#     'tables': ['polars', 'pandas'],
#     'graph': ['graph-tool', 'scipy', 'igraph', 'networkx']}
```

A no-data smoke test:

```python
from genomeblocks import Loci

loci = Loci.from_records([("chr1", 100, 500), ("chr1", 300, 700)])
print(loci)                  # -> Loci(n=2, chroms=1)
print(loci.merge())          # -> Loci(n=1, chroms=1, sorted)
print(loci.slop(100))        # -> Loci(n=2, chroms=1)
print(loci.merge().to_pandas())
#   chrom  start  end strand
# 0  chr1    100  700      .
```

The public surface is the set of names on the package:

```python
import genomeblocks as gb
gb.__all__
# -> ['Architecture', 'Atlas', 'Genes', 'Genome', 'Loci', 'Locus', 'Pairs', 'View',
#     'as_loci', 'backends', 'browser', 'compare_heatmap', 'coverage', 'igv_html',
#     'load_motifs', 'plot_motif_heatmap', 'read_fasta', 'tmm', 'use_backend']
```

---

## When a backend is missing

Asking for an engine that is not installed is an error that names the fix — never a silent switch to another engine:

```python
loci.merge(backend="cgranges")
# ImportError: the 'cgranges' intervals backend is not installed:
#   conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)

with gb.use_backend(graph="graph-tool"):     # checked on entry, before any work
    ...
# ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool

loci.merge(backend="nope")
# ValueError: unknown intervals backend 'nope'; choose from: genomeblocks, cgranges, ncls, bioframe, pyranges, bedtools

from genomeblocks import motifs as gm        # a bare pip install: no motif engine
gm.scan_motifs_matrix(loci, "genome.fa", "motifs.jaspar")
# ImportError: no motifs backend is installed: conda install -c bioconda moods  (or pip install 'genomeblocks[motifs]', which compiles MOODS-python) or pip install 'genomeblocks[lightmotif]'  (prebuilt wheels)
```

{: .warning }
> `pybedtools` writes temporary files. Point `TMPDIR` (or `pybedtools.set_tempdir()`) at a project directory before using the `bedtools` backend on large sets.

---

## Running the tests

The suite is synthetic (no genome, bigWig or ChIP-Atlas download) and runs every installed backend of every family against the default:

```bash
pip install -e ".[test]"
pytest
```

Backends that are not installed are skipped, so the suite passes on a plain pip install as well as on the full `environment.yml`.
