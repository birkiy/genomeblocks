---
title: Credits
layout: default
nav_order: 90
description: "Third-party tools that genomeblocks builds on — and how to cite them."
permalink: /credits/
---

# Credits & references
{: .no_toc }

`genomeblocks` is a thin, opinionated surface on top of excellent open-source libraries. If a particular backend or converter was load-bearing in your analysis, please cite the upstream tool in addition to `genomeblocks`. The list below maps each library to where it is used and gives the canonical citation. `gb.backends()` tells you which engines ran.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The boundary and the protocols

### narwhals — one frame API over pandas, polars and Arrow
{: #narwhals }

Used by `gb.as_loci` / `genomeblocks.interop.frame` to read any pandas, polars (eager or lazy), pyarrow or interchange-protocol table once at the boundary, and by the `__narwhals_dataframe__` protocol that lets altair, plotly and other narwhals-based tools take every genomeblocks table directly.

- **Authors:** Marco Gorelli and the narwhals contributors
- **Repository:** [github.com/narwhals-dev/narwhals](https://github.com/narwhals-dev/narwhals)
- **License:** MIT

> Gorelli, M. et al. (2024). *Narwhals: an extremely lightweight compatibility layer between dataframe libraries.* GitHub: narwhals-dev/narwhals.

{% raw %}
```bibtex
@misc{narwhals,
  author = {Gorelli, Marco and the Narwhals contributors},
  title  = {Narwhals: Extremely lightweight compatibility layer between dataframe libraries},
  year   = {2024},
  url    = {https://github.com/narwhals-dev/narwhals}
}
```
{% endraw %}

### Apache Arrow (pyarrow) — the exchange format
{: #pyarrow }

Every container's `to_arrow()` is the basis of `to_polars()`, the Arrow C stream protocol, the dataframe interchange protocol, and parquet `save` / `load`. Numeric columns cross to polars and DuckDB without a copy.

- **Authors:** The Apache Arrow developers
- **Site:** [arrow.apache.org](https://arrow.apache.org/)
- **License:** Apache-2.0

> Apache Arrow Developers (2016–). *Apache Arrow: A cross-language development platform for in-memory data.* https://arrow.apache.org/

{% raw %}
```bibtex
@misc{arrow,
  author = {{Apache Arrow Developers}},
  title  = {Apache Arrow: A cross-language development platform for in-memory data},
  year   = {2016},
  url    = {https://arrow.apache.org/}
}
```
{% endraw %}

### polars — table parsing and a hand-off target
{: #polars }

The default `tables` backend when installed (`pip install "genomeblocks[fast]"`): GTF / GFF3 / BED / BEDPE parsing in `Loci.make`, `Genes.make` and `Pairs.make` (pandas gives the same columns otherwise). Also `to_polars()` on every container and polars frames as inputs everywhere.

- **Authors:** Ritchie Vink and the polars contributors
- **Repository:** [github.com/pola-rs/polars](https://github.com/pola-rs/polars)
- **License:** MIT

> Vink, R. et al. (2023). *Polars: Blazingly fast DataFrames in Rust, Python, Node.js, R and SQL.* GitHub: pola-rs/polars.

{% raw %}
```bibtex
@misc{polars,
  author = {Vink, Ritchie and the Polars contributors},
  title  = {Polars: Blazingly fast DataFrames in Rust, Python, Node.js, R and SQL},
  year   = {2023},
  url    = {https://github.com/pola-rs/polars}
}
```
{% endraw %}

---

## Interval engines

### cgranges — C interval overlap index
{: #cgranges }

An optional `intervals` backend (`backend="cgranges"`) for overlap and point lookups; `Loci.to_cgranges()` builds the index for your own code. conda-only (`conda install -c bioconda cgranges`).

- **Author:** Heng Li (Broad Institute)
- **Repository:** [github.com/lh3/cgranges](https://github.com/lh3/cgranges)
- **License:** MIT

> Li, H. (2021). *cgranges: a C library for computing genomic interval overlaps.* GitHub: lh3/cgranges.

{% raw %}
```bibtex
@misc{cgranges,
  author = {Li, Heng},
  title  = {cgranges: A C library for computing interval overlaps on a generic genome},
  year   = {2021},
  url    = {https://github.com/lh3/cgranges}
}
```
{% endraw %}

### PyRanges — genomic interval DataFrames
{: #pyranges }

An optional `intervals` backend (`backend="pyranges"`: overlap, nearest, merge, point lookups) and a converter target: `Loci.from_pyranges` / `to_pyranges`, `Genes.from_frame(pyranges.read_gtf(...))`.

- **Authors:** Endre Bakken Stovner and Pål Sætrom (NTNU)
- **Repository:** [github.com/pyranges/pyranges](https://github.com/pyranges/pyranges)
- **License:** MIT

> Stovner, E. B. and Sætrom, P. (2020). *PyRanges: efficient comparison of genomic intervals in Python.* Bioinformatics 36(3): 918–919. doi: [10.1093/bioinformatics/btz615](https://doi.org/10.1093/bioinformatics/btz615)

{% raw %}
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
{% endraw %}

### NCLS — nested containment lists
{: #ncls }

An optional `intervals` backend (`backend="ncls"`) for overlap and point lookups, through the `ncls` Python package from the PyRanges project.

- **Authors:** Alexander V. Alekseyenko and Christopher J. Lee (algorithm); Endre Bakken Stovner (Python package)
- **Repository:** [github.com/pyranges/ncls](https://github.com/pyranges/ncls)
- **License:** BSD-3-Clause

> Alekseyenko, A. V. and Lee, C. J. (2007). *Nested Containment List (NCList): a new algorithm for accelerating interval query of genome alignment and interval databases.* Bioinformatics 23(11): 1386–1393. doi: [10.1093/bioinformatics/btl647](https://doi.org/10.1093/bioinformatics/btl647)

{% raw %}
```bibtex
@article{alekseyenko2007ncls,
  author  = {Alekseyenko, Alexander V. and Lee, Christopher J.},
  title   = {Nested Containment List (NCList): a new algorithm for accelerating interval query of genome alignment and interval databases},
  journal = {Bioinformatics},
  volume  = {23},
  number  = {11},
  pages   = {1386--1393},
  year    = {2007},
  doi     = {10.1093/bioinformatics/btl647}
}
```
{% endraw %}

### bioframe — genomic intervals in pandas
{: #bioframe }

An optional `intervals` backend (`backend="bioframe"`: overlap, nearest, merge, point lookups) and the column convention of `Loci.to_pandas()` / `to_bioframe()` (`chrom`, `start`, `end`).

- **Authors:** Open2C (Nezar Abdennur et al.)
- **Repository:** [github.com/open2c/bioframe](https://github.com/open2c/bioframe)
- **License:** MIT

> Open2C, Abdennur, N., Fudenberg, G., Flyamer, I. M., Galitsyna, A. A., Goloborodko, A., Imakaev, M. and Venev, S. V. (2024). *Bioframe: operations on genomic intervals in Pandas dataframes.* Bioinformatics 40(2): btae088. doi: [10.1093/bioinformatics/btae088](https://doi.org/10.1093/bioinformatics/btae088)

{% raw %}
```bibtex
@article{open2c2024bioframe,
  author  = {{Open2C} and Abdennur, Nezar and Fudenberg, Geoffrey and Flyamer, Ilya M. and Galitsyna, Aleksandra A. and Goloborodko, Anton and Imakaev, Maxim and Venev, Sergey V.},
  title   = {Bioframe: operations on genomic intervals in Pandas dataframes},
  journal = {Bioinformatics},
  volume  = {40},
  number  = {2},
  pages   = {btae088},
  year    = {2024},
  doi     = {10.1093/bioinformatics/btae088}
}
```
{% endraw %}

### pybedtools / BEDTools
{: #pybedtools }

An optional `intervals` backend (`backend="bedtools"`: overlap, nearest, merge, point lookups through the `bedtools` binary) and a converter target: `Loci.from_bedtool` / `to_bedtool`. The binary is conda-only (`conda install -c bioconda bedtools`).

- **Authors:** Ryan K. Dale, Brent S. Pedersen and Aaron R. Quinlan (pybedtools); Aaron R. Quinlan and Ira M. Hall (BEDTools)
- **Repository:** [github.com/daler/pybedtools](https://github.com/daler/pybedtools), [github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)
- **License:** MIT

> Dale, R. K., Pedersen, B. S. and Quinlan, A. R. (2011). *Pybedtools: a flexible Python library for manipulating genomic datasets and annotations.* Bioinformatics 27(24): 3423–3424. doi: [10.1093/bioinformatics/btr539](https://doi.org/10.1093/bioinformatics/btr539)
>
> Quinlan, A. R. and Hall, I. M. (2010). *BEDTools: a flexible suite of utilities for comparing genomic features.* Bioinformatics 26(6): 841–842. doi: [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

{% raw %}
```bibtex
@article{dale2011pybedtools,
  author  = {Dale, Ryan K. and Pedersen, Brent S. and Quinlan, Aaron R.},
  title   = {Pybedtools: a flexible Python library for manipulating genomic datasets and annotations},
  journal = {Bioinformatics},
  volume  = {27},
  number  = {24},
  pages   = {3423--3424},
  year    = {2011},
  doi     = {10.1093/bioinformatics/btr539}
}
@article{quinlan2010bedtools,
  author  = {Quinlan, Aaron R. and Hall, Ira M.},
  title   = {BEDTools: a flexible suite of utilities for comparing genomic features},
  journal = {Bioinformatics},
  volume  = {26},
  number  = {6},
  pages   = {841--842},
  year    = {2010},
  doi     = {10.1093/bioinformatics/btq033}
}
```
{% endraw %}

---

## Graph engines

### graph-tool — chromatin-contact networks
{: #graph-tool }

The default `graph` backend when installed (scipy otherwise): `Architecture.components`, `pagerank`, the spring layout of `A.draw`, and `A.graph()` / `A.to_graph_tool()` build a `graph_tool.Graph` from the edge table on demand. graph-tool is a compiled Boost/C++ library and comes from conda (`conda install -c conda-forge graph-tool`).

- **Author:** Tiago P. Peixoto
- **Site:** [graph-tool.skewed.de](https://graph-tool.skewed.de/)
- **License:** LGPL-3.0

> Peixoto, T. P. (2014). *The graph-tool python library.* figshare. doi: [10.6084/m9.figshare.1164194](https://doi.org/10.6084/m9.figshare.1164194)

{% raw %}
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
{% endraw %}

{: .note }
> If you run block-model inference, community detection or SBM layouts on `A.to_graph_tool()`, please also cite the specific method paper — Peixoto maintains [a citation guide](https://graph-tool.skewed.de/static/doc/index.html) for each algorithm.

### igraph
{: #igraph }

An optional `graph` backend (`backend="igraph"`) and a converter target: `Architecture.to_igraph()`.

- **Authors:** Gábor Csárdi and Tamás Nepusz
- **Site:** [igraph.org](https://igraph.org/)
- **License:** GPL-2.0-or-later

> Csárdi, G. and Nepusz, T. (2006). *The igraph software package for complex network research.* InterJournal, Complex Systems 1695.

{% raw %}
```bibtex
@article{csardi2006igraph,
  author  = {Cs{\'a}rdi, G{\'a}bor and Nepusz, Tam{\'a}s},
  title   = {The igraph software package for complex network research},
  journal = {InterJournal},
  volume  = {Complex Systems},
  pages   = {1695},
  year    = {2006},
  url     = {https://igraph.org}
}
```
{% endraw %}

### NetworkX
{: #networkx }

An optional `graph` backend (`backend="networkx"`) and a converter target: `Architecture.to_networkx()`.

- **Authors:** Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart
- **Site:** [networkx.org](https://networkx.org/)
- **License:** BSD-3-Clause

> Hagberg, A. A., Schult, D. A. and Swart, P. J. (2008). *Exploring network structure, dynamics, and function using NetworkX.* Proceedings of the 7th Python in Science Conference (SciPy 2008), 11–15.

{% raw %}
```bibtex
@inproceedings{hagberg2008networkx,
  author    = {Hagberg, Aric A. and Schult, Daniel A. and Swart, Pieter J.},
  title     = {Exploring network structure, dynamics, and function using {NetworkX}},
  booktitle = {Proceedings of the 7th Python in Science Conference (SciPy 2008)},
  pages     = {11--15},
  year      = {2008}
}
```
{% endraw %}

---

## Sequence and motif engines

### MOODS — motif occurrence detection
{: #moods }

The default `motifs` backend: `scan_motifs`, `scan_motifs_matrix`, `scan_motifs_matrix_masked` and `scan_motifs_profile` hand each log-odds matrix — a whole batch of them for a library — to a `MOODS.scan.Scanner`, and `Library.to_moods()` exports the matrices. The conda package depends on bioconda's `moods`; on pip it is `pip install 'genomeblocks[motifs]'` (MOODS-python, compiled at install).

- **Authors:** Janne H. Korhonen, Petri Martinmäki, Cinzia Pizzi, Pasi Rastas, Esko Ukkonen
- **Repository:** [github.com/jhkorhonen/MOODS](https://github.com/jhkorhonen/MOODS)
- **License:** GPL-3.0 / Biopython License (dual)

> Korhonen, J., Martinmäki, P., Pizzi, C., Rastas, P. and Ukkonen, E. (2009). *MOODS: fast search for position weight matrix matches in DNA sequences.* Bioinformatics 25(23): 3181–3182. doi: [10.1093/bioinformatics/btp554](https://doi.org/10.1093/bioinformatics/btp554)
>
> Korhonen, J. H., Palin, K., Taipale, J. and Ukkonen, E. (2017). *Fast motif matching revisited: high-order PWMs, SNPs and indels.* Bioinformatics 33(4): 514–521. doi: [10.1093/bioinformatics/btw683](https://doi.org/10.1093/bioinformatics/btw683)

{% raw %}
```bibtex
@article{korhonen2009moods,
  author  = {Korhonen, Janne and Martinm{\"a}ki, Petri and Pizzi, Cinzia and Rastas, Pasi and Ukkonen, Esko},
  title   = {MOODS: fast search for position weight matrix matches in DNA sequences},
  journal = {Bioinformatics},
  volume  = {25},
  number  = {23},
  pages   = {3181--3182},
  year    = {2009},
  doi     = {10.1093/bioinformatics/btp554}
}
```
{% endraw %}

### lightmotif — Rust-backed PSSM scanning
{: #lightmotif }

The alternative `motifs` backend (`backend="lightmotif"`, `pip install 'genomeblocks[lightmotif]'`, prebuilt wheels), picked automatically only when MOODS is absent: its SIMD core scans the striped block one matrix at a time and reports the same hits. genomeblocks' `logodds_matrix` reproduces lightmotif's float32 log-odds arithmetic bit for bit, and `load_motifs` accepts lightmotif motif objects.

- **Author:** Martin Larralde (EMBL)
- **Repository:** [github.com/althonos/lightmotif](https://github.com/althonos/lightmotif)
- **License:** MIT

> Larralde, M. (2023). *lightmotif: PSSM scoring with SIMD in Python and Rust.* GitHub: althonos/lightmotif.

{% raw %}
```bibtex
@misc{lightmotif,
  author = {Larralde, Martin},
  title  = {lightmotif: A lightweight library for PSSM scoring, with SIMD backends},
  year   = {2023},
  url    = {https://github.com/althonos/lightmotif}
}
```
{% endraw %}

### pyfaidx — indexed FASTA access
{: #pyfaidx }

An optional `fasta` backend (`backend="pyfaidx"`) for `Loci.sequences`, motif scanning windows and the browser's reference track.

- **Authors:** Matthew D. Shirley, Zhaorong Ma, Brent S. Pedersen, Sarah J. Wheelan
- **Repository:** [github.com/mdshw5/pyfaidx](https://github.com/mdshw5/pyfaidx)
- **License:** BSD-3-Clause

> Shirley, M. D., Ma, Z., Pedersen, B. S. and Wheelan, S. J. (2015). *Efficient "pythonic" access to FASTA files using pyfaidx.* PeerJ PrePrints 3: e1196. doi: [10.7287/peerj.preprints.970v1](https://doi.org/10.7287/peerj.preprints.970v1)

{% raw %}
```bibtex
@article{shirley2015pyfaidx,
  author  = {Shirley, Matthew D. and Ma, Zhaorong and Pedersen, Brent S. and Wheelan, Sarah J.},
  title   = {Efficient ``pythonic'' access to FASTA files using pyfaidx},
  journal = {PeerJ PrePrints},
  volume  = {3},
  pages   = {e1196},
  year    = {2015},
  doi     = {10.7287/peerj.preprints.970v1}
}
```
{% endraw %}

### pyliftover — assembly lift-over
{: #pyliftover }

`Loci.liftover(chain_file)` reads UCSC chain files through pyliftover (`pip install pyliftover`); genomeblocks adds the inward walk that implements UCSC's base-fraction rule.

- **Author:** Konstantin Tretyakov
- **Repository:** [github.com/konstantint/pyliftover](https://github.com/konstantint/pyliftover)
- **License:** MIT

> Tretyakov, K. (2013). *pyliftover: pure-Python implementation of UCSC liftOver genome coordinate conversion.* GitHub: konstantint/pyliftover.

{% raw %}
```bibtex
@misc{pyliftover,
  author = {Tretyakov, Konstantin},
  title  = {pyliftover: Pure-python implementation of UCSC liftOver genome coordinate conversion},
  year   = {2013},
  url    = {https://github.com/konstantint/pyliftover}
}
```
{% endraw %}

---

## Signal engines

### pybigtools — threaded bigWig reader
{: #pybigtools }

The default `bigwig` backend: `loci.signal()`, `compare_heatmap`, `select_isoforms(bw=...)`, the browser's bigWig tracks and `igv_html`'s embedded signal all read `.bw` / `.bigWig` files through pybigtools' Rust core. `backend="pybigwig"` and the pure-Python `backend="python"` reader give the same bins.

- **Author:** Jack Huey (the bigtools project)
- **Repository:** [github.com/jackh726/bigtools](https://github.com/jackh726/bigtools)
- **License:** MIT

> Huey, J. (2023). *bigtools: A high-performance BigWig and BigBed library in Rust, with Python bindings.* GitHub: jackh726/bigtools.

{% raw %}
```bibtex
@misc{pybigtools,
  author = {Huey, Jack},
  title  = {bigtools: A high-performance BigWig/BigBed library in Rust with Python bindings (pybigtools)},
  year   = {2023},
  url    = {https://github.com/jackh726/bigtools}
}
```
{% endraw %}

---

## Data containers

### AnnData
{: #anndata }

`Loci.from_anndata(adata)` reads scATAC peaks from an AnnData's `var` (columns or `chr1:100-200` names) in its order; `Loci.to_anndata(X)` and `Architecture.to_anndata()` (adjacencies in `obsp`) build one; `interop.cube_to_anndata` stores a signal cube.

- **Authors:** Isaac Virshup, Sergei Rybakov, Fabian J. Theis, Philipp Angerer, F. Alexander Wolf
- **Repository:** [github.com/scverse/anndata](https://github.com/scverse/anndata)
- **License:** BSD-3-Clause

> Virshup, I., Rybakov, S., Theis, F. J., Angerer, P. and Wolf, F. A. (2024). *anndata: Access and store annotated data matrices.* Journal of Open Source Software 9(101): 4371. doi: [10.21105/joss.04371](https://doi.org/10.21105/joss.04371)

{% raw %}
```bibtex
@article{virshup2024anndata,
  author  = {Virshup, Isaac and Rybakov, Sergei and Theis, Fabian J. and Angerer, Philipp and Wolf, F. Alexander},
  title   = {anndata: Access and store annotated data matrices},
  journal = {Journal of Open Source Software},
  volume  = {9},
  number  = {101},
  pages   = {4371},
  year    = {2024},
  doi     = {10.21105/joss.04371}
}
```
{% endraw %}

### xarray
{: #xarray }

`interop.cube_to_xarray(S, loci, tracks)` labels a `(rows, tracks, bins)` signal cube with region, track and bin-offset coordinates.

- **Authors:** Stephan Hoyer, Joe Hamman and the xarray developers
- **Repository:** [github.com/pydata/xarray](https://github.com/pydata/xarray)
- **License:** Apache-2.0

> Hoyer, S. and Hamman, J. (2017). *xarray: N-D labeled arrays and datasets in Python.* Journal of Open Research Software 5(1): 10. doi: [10.5334/jors.148](https://doi.org/10.5334/jors.148)

{% raw %}
```bibtex
@article{hoyer2017xarray,
  author  = {Hoyer, Stephan and Hamman, Joe},
  title   = {xarray: {N-D} labeled arrays and datasets in {Python}},
  journal = {Journal of Open Research Software},
  volume  = {5},
  number  = {1},
  pages   = {10},
  year    = {2017},
  doi     = {10.5334/jors.148}
}
```
{% endraw %}

---

## Also used

| Library | Where | Citation |
|---|---|---|
| [cooler](https://github.com/open2c/cooler) (MIT) | `Architecture.add_mcool`; chromosome sizes from a cooler in `Loci.tile` / `Genome.from_sizes` | Abdennur, N. and Mirny, L. A. (2020). Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 36(1): 311–316. doi: [10.1093/bioinformatics/btz540](https://doi.org/10.1093/bioinformatics/btz540) |
| [pyBigWig](https://github.com/deeptools/pyBigWig) (MIT) | the `pybigwig` bigwig backend | Ramírez, F. et al. (2016). deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Research* 44(W1): W160–W165. doi: [10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257) |
| [pysam](https://github.com/pysam-developers/pysam) / htslib (MIT) | the `pysam` fasta backend; BAM pileups in `browser` and `gb.coverage` | Bonfield, J. K. et al. (2021). HTSlib: C library for reading/writing high-throughput sequencing data. *GigaScience* 10(2): giab007. doi: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007) |
| [Biopython](https://biopython.org/) (Biopython License) | the `biopython` motifs and fasta backends; `Library.to_biopython`, `load_motifs(Bio.motifs)`, `Loci.to_seqrecords` | Cock, P. J. A. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11): 1422–1423. doi: [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) |
| [SciPy](https://scipy.org/) (BSD-3) | the default `graph` backend when graph-tool is absent (sparse adjacency, components, PageRank); the power-law fit of `normalize`; Fisher tests in `Atlas` | Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17: 261–272. doi: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) |
| [pandas](https://pandas.pydata.org/) (BSD-3) | the fallback `tables` backend; `to_pandas()`; `gb.backends()` | McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61. doi: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a) |
| [NumPy](https://numpy.org/) (BSD-3) | every table | Harris, C. R. et al. (2020). Array programming with NumPy. *Nature* 585: 357–362. doi: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| [edgeR](https://bioconductor.org/packages/edgeR/) (GPL) | the TMM algorithm vendored in `gb.tmm` | Robinson, M. D. and Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology* 11: R25. doi: [10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25) |
| [igv.js](https://github.com/igvteam/igv.js) (MIT) | `gb.igv_html` loads igv.js from a CDN or a local file | Robinson, J. T. et al. (2023). igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). *Bioinformatics* 39(1): btac830. doi: [10.1093/bioinformatics/btac830](https://doi.org/10.1093/bioinformatics/btac830) |

---

## At-a-glance map

| Family / module | Default | Other engines and targets |
|---|---|---|
| [Loci]({{ '/guide/loci/' | relative_url }}) — intervals | genomeblocks (numpy) | cgranges, ncls, bioframe, pyranges, bedtools; pandas, polars, Arrow, AnnData, pyliftover |
| [Genes]({{ '/guide/genes/' | relative_url }}) — annotation tables | polars / pandas parsing | the interval engines above; PyRanges frames in |
| [Architecture]({{ '/guide/architecture/' | relative_url }}) — contact networks | graph-tool, else scipy | igraph, networkx, AnnData; cooler for Hi-C weights |
| [Signal]({{ '/guide/signal/' | relative_url }}) — bigWig cubes | pybigtools | pyBigWig, the pure-Python reader; xarray, AnnData out; edgeR TMM |
| [Motifs]({{ '/guide/motifs/' | relative_url }}) — PSSM scanning | MOODS, else lightmotif (`genomeblocks[motifs]` / `[lightmotif]`) | Biopython; FASTA via genomeblocks, pysam, pyfaidx, Biopython |
| [Browser]({{ '/guide/browser/' | relative_url }}) — views | matplotlib | pysam for BAM tracks; igv.js for `igv_html` |
| [Atlas]({{ '/guide/atlas/' | relative_url }}) — enrichment | genomeblocks (scipy sparse) | — |

---

## Citing `genomeblocks`

If `genomeblocks` itself is useful, please also cite the repository:

```
Altintas, U. B. (2024). genomeblocks: Fluent building blocks for regulatory genomics.
https://github.com/birkiy/genomeblocks
```

{% raw %}
```bibtex
@misc{genomeblocks,
  author = {Altintas, Umut Berkay},
  title  = {genomeblocks: Fluent building blocks for regulatory genomics},
  year   = {2024},
  url    = {https://github.com/birkiy/genomeblocks}
}
```
{% endraw %}
