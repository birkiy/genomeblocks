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
- **[graph-tool](https://graph-tool.skewed.de/)** — chromatin-architecture graphs (Peixoto 2014; not pip-installable, use conda/mamba)
- Standard scientific stack: `numpy`, `pandas`, `scipy`, `matplotlib`
- Genomics I/O:
  - **[pybigtools](https://github.com/jackh726/bigtools)** — threaded bigWig reader (Huey 2023)
  - **[pyranges](https://github.com/pyranges/pyranges)** — genomic interval DataFrames (Stovner & Sætrom 2020)
  - **[cgranges](https://github.com/lh3/cgranges)** — C-level interval overlap index (Heng Li)
  - `cooler` — Hi-C `.mcool` I/O
- Optional:
  - **[lightmotif](https://github.com/althonos/lightmotif)** — SIMD-accelerated PSSM scanning (Larralde 2023)
  - `conorm` — TMM normalization

Because `graph-tool` is a compiled C++/Boost library, the recommended path is conda. See the [Credits page](credits) for full citations of every upstream tool.

---

## Recommended: conda environment

```bash
git clone https://github.com/birkiy/genomeblocks.git
cd genomeblocks
conda env create -f environment.yml
conda activate genomeblocks
pip install -e .
```

The shipped `environment.yml` pins tested versions of `graph-tool`, `cooler`, and `pybigtools`.

---

## Alternative: mamba (faster)

```bash
mamba env create -f environment.yml
mamba activate genomeblocks
pip install -e .
```

---

## pip-only (no graph-tool → no Architecture)

If you only need `Loci`, `Genes`, `Signal`, and `Browser`, you can skip graph-tool:

```bash
pip install genomeblocks
```

The lazy-import in `genomeblocks/__init__.py` means `Architecture` will only fail the moment you touch it. All other subsystems work.

{: .note }
> Browser, signal extraction, and Genes don't need graph-tool. You'll get a clean import error only if you reference `Architecture`.

---

## Verifying the install

```python
import genomeblocks as gb
print(gb.__all__)
# ['Architecture', 'Atlas', 'CDS', 'Exon', 'Gene', 'Genes', 'Loci', 'Locus',
#  'Transcript', 'UTR', 'browser', 'compare_heatmap', 'make_genome',
#  'scan_motifs', 'tmm']
```

Try a no-data smoke test:

```python
from genomeblocks import Loci, Locus
loci = Loci([Locus("chr1", 100, 500), Locus("chr1", 300, 700)])
print(loci.merge())          # Loci(n=1)
print(loci.slop(100))        # Loci(n=2)  with ±100 bp
```

---

## Optional extras

| Feature | Extra dependency |
|---|---|
| `Architecture.*` | `graph-tool` (conda) |
| `scan_motifs()` | `lightmotif` |
| `tmm()` | `conorm` |
| `Architecture.add_mcool()` | `cooler` |
| BigWig signal (fast path) | `pybigtools` (falls back to pure Python) |
| `pyranges`-backed ops (`Loci.nearest`, `Genes.nearest_genes`) | `pyranges` |

All are declared in `pyproject.toml`; the only one that absolutely needs conda is `graph-tool`.
