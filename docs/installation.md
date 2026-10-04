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
- Genomics I/O (all pip-installable):
  - **[pybigtools](https://github.com/jackh726/bigtools)** — threaded bigWig reader (Huey 2023)
  - **[pyranges](https://github.com/pyranges/pyranges)** — genomic interval DataFrames (Stovner & Sætrom 2020)
  - **[lightmotif](https://github.com/althonos/lightmotif)** — SIMD-accelerated PSSM scanning (Larralde 2023)
  - `cooler` — Hi-C `.mcool` I/O
- **Optional speedup** — **[cgranges](https://github.com/lh3/cgranges)** (Heng Li), a C interval-overlap index. It is **not on PyPI**; genomeblocks uses a pure-Python fallback when it's absent, so overlap ops work without it. Install it (conda, or `pip install git+https://github.com/lh3/cgranges`) for the fast path on large sets.

Everything except `graph-tool` and the optional `cgranges` is pip-installable
and declared in `pyproject.toml`. TMM normalization (`tmm()`) is **vendored** —
the edgeR algorithm ships inside `genomeblocks.signal`, so there is no external
normalization dependency. Because `graph-tool` is a compiled C++/Boost library,
the recommended path is conda. See the [Credits page]({{ '/credits/' | relative_url }}) for full
citations of every upstream tool.

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

## pip-only install → **no `Architecture`**

{: .warning }
> **`graph-tool` is not on PyPI.** `pip install genomeblocks` gives you every
> subsystem **except `Architecture`** (the chromatin-contact graph), which
> imports `graph-tool` at first use. If you need `Architecture`, you must use a
> conda/mamba environment. This is the single most common install surprise.

```bash
pip install genomeblocks
```

Thanks to the lazy imports in `genomeblocks/__init__.py`, the missing dependency
only surfaces the moment you touch `Architecture` (or `architecture_draw`) — you
get a clean `ImportError` for `graph_tool`, not a broken package.

| Works pip-only | Needs conda (`graph-tool`) |
|---|---|
| `Loci`, `Locus`, `Genes` | `Architecture` |
| `signal` / `tmm` / heatmaps | `architecture_draw.draw` |
| `browser`, `Atlas`, `scan_motifs`, `bedpe` | — |

---

## Verifying the install

```python
import genomeblocks as gb
print(gb.__all__)
# ['Architecture', 'Atlas', 'CDS', 'Exon', 'Gene', 'Genes', 'Loci', 'Locus',
#  'Transcript', 'UTR', 'browser', 'compare_heatmap', 'coverage',
#  'make_genome', 'scan_motifs', 'tmm']
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

| Feature | Dependency | On PyPI? |
|---|---|---|
| `Architecture.*` | `graph-tool` | ❌ conda only |
| `scan_motifs()` / motif matrices | `lightmotif` | ✅ (auto) |
| `tmm()` | — (edgeR TMM vendored) | — |
| `Architecture.add_mcool()` | `cooler` | ✅ (auto) |
| BigWig signal (fast path) | `pybigtools` | ✅ (auto; falls back to pure Python) |
| Interval overlap (fast path) | `cgranges` | ❌ conda / from source; falls back to pure Python |
| `pyranges`-backed ops (`Loci.nearest`, `Genes.nearest_genes`) | `pyranges` | ✅ (auto) |
| BAM tracks in `browser` / `coverage()` | `pysam` | ✅ `pip install genomeblocks[bam]` |
| Motif logos (`motifs_draw`) | `logomaker` | ✅ (install separately) |

Everything marked ✅ is declared in `pyproject.toml` and installed by pip. The
`graph-tool` and `cgranges` fast paths need conda/source — but only
`graph-tool` is truly required (for `Architecture`); `cgranges` is a pure
speedup with a built-in pure-Python fallback.
