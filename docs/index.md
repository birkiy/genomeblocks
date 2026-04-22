---
title: Home
layout: default
nav_order: 1
description: "Fluent building blocks for regulatory genomics in Python."
permalink: /
---

<span class="gb-eyebrow">v0.9 · MIT · <code>conda env create -f environment.yml</code></span>

# genomeblocks
{: .fs-9 }

Fluent building blocks for regulatory genomics; from peaks to chromatin networks in a handful of expressive chained calls.
{: .fs-6 .fw-300 }

[Quickstart](quickstart){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/birkiy/genomeblocks){: .btn .fs-5 .mb-4 .mb-md-0 }

<img src="{{ '/assets/images/browser_Nanog.svg' | relative_url }}" alt="Genomeblocks browser — IGV-like region viewer" class="gb-browser-hero" />

---

## Why genomeblocks?

Regulatory-genomics analyses usually end up as a cocktail of bedtools, PyRanges, cooler, pyBigWig, GTF parsing boilerplate, graph libraries, and one-off heatmap code. `genomeblocks` unifies those pieces behind a small set of composable objects:

<ul class="gb-modules">
  <li><strong>Loci</strong><br><small>A list-of-intervals container with set algebra (<code>&</code>, <code>|</code>, <code>-</code>, <code>^</code>), <code>slop</code>, <code>sort</code>, <code>merge</code>, <code>nearest</code>, indexed overlap queries, and signal extraction.</small></li>
  <li><strong>Locus</strong><br><small>A single interval with a canonical UID (<code>chrom:start-end(strand)</code>) used everywhere as a stable key.</small></li>
  <li><strong>Tags</strong><br><small>An in-memory annotation store for a <code>Loci</code> set, queryable by boolean expression (<code>l.atac & (l.h3k27ac > 1.5)</code>).</small></li>
  <li><strong>Genes</strong><br><small>GENCODE / GTF / UCSC parser yielding <code>Gene</code> → <code>Transcript</code> → <code>Exon</code>/<code>CDS</code>/<code>UTR</code> hierarchies, with enhancer-to-gene assignment in one call.</small></li>
  <li><strong>Architecture</strong><br><small>A chromatin-contact graph (graph-tool) built from BEDPE loops or mcool matrices; spreading, cliques, O/E, hub discovery.</small></li>
  <li><strong>signal</strong><br><small>Threaded bigWig extraction (pybigtools backend, pure-Python fallback), TMM normalization, and comparative heatmaps.</small></li>
  <li><strong>browser</strong><br><small>An IGV-like, SVG-clean multi-track region viewer built on matplotlib.</small></li>
  <li><strong>bedpe / motifs</strong><br><small>BEDPE parsing + pair-to-bed intersection; JASPAR motif scanning over a FASTA genome.</small></li>
</ul>

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
| [Credits](credits) | Upstream tools & citations. |

---

## Citing

If `genomeblocks` is useful in your work, please cite the GitHub repository:

```
Altintas, U. B. (2024). genomeblocks: Fluent building blocks for regulatory genomics.
https://github.com/birkiy/genomeblocks
```

`genomeblocks` stands on top of several excellent upstream libraries — please also cite the ones whose module was load-bearing in your analysis. The [Credits](credits) page lists them all with BibTeX.

---

## License

MIT. See [LICENSE](https://github.com/birkiy/genomeblocks/blob/main/LICENSE).
