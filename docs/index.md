---
title: Home
layout: default
nav_order: 1
description: "Fluent building blocks for regulatory genomics in Python."
permalink: /
---

# genomeblocks
{: .fs-9 }

Fluent building blocks for regulatory genomics — from peaks to chromatin networks in a handful of expressive chained calls.
{: .fs-6 .fw-300 }

[Quickstart](quickstart){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/birkiy/genomeblocks){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Why genomeblocks?

Regulatory-genomics analyses usually end up as a cocktail of bedtools, PyRanges, cooler, pyBigWig, GTF parsing boilerplate, graph libraries, and one-off heatmap code. `genomeblocks` unifies those pieces behind a small set of composable objects:

- **`Loci`** — a list-of-intervals container with set algebra (`&`, `|`, `-`, `^`), `slop`, `sort`, `merge`, `nearest`, indexed overlap queries, and signal extraction.
- **`Locus`** — a single interval with a canonical UID (`chrom:start-end(strand)`) used everywhere as a stable key.
- **`Tags`** — an in-memory annotation store for a `Loci` set, queryable by boolean expression (`l.atac & (l.h3k27ac > 1.5)`).
- **`Genes`** — GENCODE / GTF / UCSC parser yielding `Gene` → `Transcript` → `Exon`/`CDS`/`UTR` hierarchies, with enhancer-to-gene assignment in one call.
- **`Architecture`** — a chromatin-contact graph (graph-tool) built from BEDPE loops or mcool matrices; supports spreading, clique construction, O/E normalization, and hub/focus-gene discovery.
- **`signal`** — threaded bigWig extraction (pybigtools backend, pure-Python fallback), TMM normalization, and comparative heatmaps.
- **`browser`** — an IGV-like, SVG-clean multi-track region viewer built on matplotlib.
- **`bedpe`** / **`motifs`** — BEDPE parsing + pair-to-bed intersection; JASPAR motif scanning over a FASTA genome.

Everything is **chainable**: the output of one stage is always a first-class object accepted by the next.

---

## 60-second example

```python
from genomeblocks import Architecture, Genes, Loci

# CREs: ATAC peaks, extended ±100 bp, sorted, merged
cre = (Loci.make("atac_peaks.narrowPeak")
           .slop(100)
           .sort()
           .merge())

# Super-enhancers that overlap CREs
se = cre.intersect(Loci.make("H3K27ac_SE.bed"))

# Contact graph: CRE-resolved HiChIP loops, Hi-C weights, O/E normalization
arch = (Architecture.make(cre, "RNAP_loops.bedpe", r=2500)
                    .add_mcool(cre, "RNAP.mcool", resolution=5000)
                    .normalize(cre))

# Gene annotation: which annotation class does each super-enhancer fall into?
genes = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)
counts = genes.annotations(se & cre).groupby("annotation").size()
```

---

## Documentation map

| Section | When to read it |
|---|---|
| [Installation](installation) | Setup with conda or pip. |
| [Quickstart](quickstart) | A 10-minute tour end-to-end. |
| [Concepts](concepts) | The mental model — UIDs, lazy indexes, chainable APIs. |
| [User Guide → Loci](guide/loci) | Interval algebra in depth. |
| [User Guide → Tags](guide/tags) | Boolean queries over CREs. |
| [User Guide → Genes](guide/genes) | GTF parsing & enhancer-to-gene. |
| [User Guide → Architecture](guide/architecture) | Chromatin-contact networks. |
| [User Guide → Signal](guide/signal) | BigWig extraction & heatmaps. |
| [User Guide → Browser](guide/browser) | Multi-track region plots. |
| [User Guide → BEDPE](guide/bedpe) | Loops & paired intervals. |
| [User Guide → Motifs](guide/motifs) | TF motif scanning. |
| [API Reference](api/) | Full method signatures. |
| [Examples](examples) | End-to-end scripts. |

---

## Citing

If `genomeblocks` is useful in your work, please cite the GitHub repository:

```
Altintas, U. B. (2024). genomeblocks: Fluent building blocks for regulatory genomics.
https://github.com/birkiy/genomeblocks
```

---

## License

MIT. See [LICENSE](https://github.com/birkiy/genomeblocks/blob/main/LICENSE).
