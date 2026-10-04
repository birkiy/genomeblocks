---
title: Credits
layout: default
nav_order: 90
description: "Third-party tools that genomeblocks builds on — and how to cite them."
permalink: /credits/
---

# Credits & references
{: .no_toc }

`genomeblocks` is a thin, opinionated surface on top of excellent open-source libraries. If a particular module was load-bearing in your analysis, please cite the upstream tool in addition to `genomeblocks`. The list below maps each module to its dependencies and the canonical citation.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## cgranges — interval overlap index
{: #cgranges }

Used by `Loci` for every overlap query (`&`, `-`, `^`, `overlaps`, `intersect`, `map`) and transitively by `Architecture.make` when anchoring BEDPE pairs to CREs.

- **Author:** Heng Li (Broad Institute)
- **Repository:** [github.com/lh3/cgranges](https://github.com/lh3/cgranges)
- **License:** MIT

> Li, H. (2021). *cgranges: a C library for computing genomic interval overlaps.* GitHub: lh3/cgranges.

```bibtex
@misc{cgranges,
  author = {Li, Heng},
  title  = {cgranges: A C library for computing interval overlaps on a generic genome},
  year   = {2021},
  url    = {https://github.com/lh3/cgranges}
}
```

---

## PyRanges — fast genomic interval DataFrames
{: #pyranges }

Used by `Loci.nearest()`, `Loci.to_pyranges()`, `Genes.nearest_genes()`, and anywhere a DataFrame-shaped overlap report is more natural than a raw `Loci`.

- **Authors:** Endre Bakken Stovner and Pål Sætrom (NTNU)
- **Repository:** [github.com/pyranges/pyranges](https://github.com/pyranges/pyranges)
- **License:** MIT

> Stovner, E. B. and Sætrom, P. (2020). *PyRanges: efficient comparison of genomic intervals in Python.* Bioinformatics 36(3): 918–919. doi: [10.1093/bioinformatics/btz615](https://doi.org/10.1093/bioinformatics/btz615)

```bibtex
@article{stovner2020pyranges,
  author  = {Stovner, Endre Bakken and S{\ae}trom, P{\aa}l},
  title   = {PyRanges: efficient comparison of genomic intervals in Python},
  journal = {Bioinformatics},
  volume  = {36},
  number  = {3},
  pages   = {918--919},
  year    = {2020},
  doi     = {10.1093/bioinformatics/btz615}
}
```

---

## graph-tool — chromatin-contact networks
{: #graph-tool }

`Architecture` is a subclass of `graph_tool.Graph`. Every network operation — spreading, clique construction, O/E normalization, centrality, block models — is a graph-tool call under the hood. graph-tool is also the only genomeblocks dependency that must come from conda (it's a compiled Boost/C++ library).

- **Author:** Tiago P. Peixoto
- **Site:** [graph-tool.skewed.de](https://graph-tool.skewed.de/)
- **License:** LGPL-3.0

> Peixoto, T. P. (2014). *The graph-tool python library.* figshare. doi: [10.6084/m9.figshare.1164194](https://doi.org/10.6084/m9.figshare.1164194)

```bibtex
@article{peixoto_graph-tool_2014,
  author  = {Peixoto, Tiago P.},
  title   = {The graph-tool python library},
  journal = {figshare},
  year    = {2014},
  doi     = {10.6084/m9.figshare.1164194},
  url     = {https://graph-tool.skewed.de/}
}
```

{: .note }
> If you use block-model inference, community detection, or SBM-based layouts (e.g. `arch.partition()`, nested SBM plots), please also cite the specific method paper — Peixoto maintains [a citation guide](https://graph-tool.skewed.de/static/doc/index.html) for each algorithm.

---

## lightmotif — Rust-backed PSSM scanning
{: #lightmotif }

Powers `scan_motifs()` and the JASPAR-to-counts pipeline. lightmotif's SIMD-accelerated Rust core makes genome-wide motif scans fast enough to run inside a notebook.

- **Author:** Martin Larralde (EMBL)
- **Repository:** [github.com/althonos/lightmotif](https://github.com/althonos/lightmotif)
- **License:** MIT

> Larralde, M. (2023). *lightmotif: PSSM scoring with SIMD in Python and Rust.* GitHub: althonos/lightmotif.

```bibtex
@misc{lightmotif,
  author = {Larralde, Martin},
  title  = {lightmotif: A lightweight library for PSSM scoring, with SIMD backends},
  year   = {2023},
  url    = {https://github.com/althonos/lightmotif}
}
```

---

## pybigtools — threaded bigWig reader
{: #pybigtools }

The fast path for `loci.signal()` and every module that extracts continuous coverage from `.bw` / `.bigWig` files. When pybigtools isn't available, `signal` falls back to a pure-Python reader — but the threaded Rust backend is the one you want for real work.

- **Author:** Jack Huey (the bigtools project; related to the ENCODE DCC tooling)
- **Repository:** [github.com/jackh726/bigtools](https://github.com/jackh726/bigtools)
- **License:** MIT

> Huey, J. (2023). *bigtools: A high-performance BigWig and BigBed library in Rust, with Python bindings.* GitHub: jackh726/bigtools.

```bibtex
@misc{pybigtools,
  author = {Huey, Jack},
  title  = {bigtools: A high-performance BigWig/BigBed library in Rust with Python bindings (pybigtools)},
  year   = {2023},
  url    = {https://github.com/jackh726/bigtools}
}
```

---

## At-a-glance map

| Module | Upstream tool(s) |
|---|---|
| [Loci]({{ '/guide/loci/' | relative_url }}) — overlap / set algebra | cgranges, pyranges |
| [Genes]({{ '/guide/genes/' | relative_url }}) — enhancer-to-gene | pyranges |
| [Architecture]({{ '/guide/architecture/' | relative_url }}) — contact networks | graph-tool |
| [Signal]({{ '/guide/signal/' | relative_url }}) — bigWig extraction | pybigtools |
| [Motifs]({{ '/guide/motifs/' | relative_url }}) — PSSM scanning | lightmotif |

---

## Citing `genomeblocks`

If `genomeblocks` itself is useful, please also cite the repository:

```
Altintas, U. B. (2024). genomeblocks: Fluent building blocks for regulatory genomics.
https://github.com/birkiy/genomeblocks
```

```bibtex
@misc{genomeblocks,
  author = {Altintas, Umut Berkay},
  title  = {genomeblocks: Fluent building blocks for regulatory genomics},
  year   = {2024},
  url    = {https://github.com/birkiy/genomeblocks}
}
```
