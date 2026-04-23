---
title: Home
layout: default
nav_order: 1
description: "Fluent building blocks for regulatory genomics in Python."
permalink: /
---

<div class="gb-eyebrow-row">
  <span class="gb-eyebrow">v0.9 · MIT</span>
  <button type="button" class="gb-install-chip" onclick="gbCopyInstall(this)" title="Copy to clipboard">
    <span class="chip-text">pip install genomeblocks</span>
    <span class="copy-icon" aria-hidden="true">⧉</span>
  </button>
</div>

# genomeblocks
{: .fs-9 }

Fluent building blocks for regulatory genomics — from peaks to chromatin networks in a handful of expressive chained calls.
{: .fs-6 .fw-300 }

<div class="gb-works-with" aria-label="Supported data types">
  <span class="label">Works with</span>
  <span class="chip">ChIP-seq</span>
  <span class="chip">ATAC-seq</span>
  <span class="chip">Hi-C · HiChIP</span>
  <span class="chip">GTF / GENCODE</span>
  <span class="chip">BigWig</span>
  <span class="chip">JASPAR motifs</span>
</div>

[Quickstart →](quickstart){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/birkiy/genomeblocks){: .btn .fs-5 .mb-4 .mb-md-0 }


## Why genomeblocks?
{: .sec-navy }

Regulatory-genomics analyses usually end up as a cocktail of bedtools, PyRanges, cooler, pyBigWig, GTF parsing boilerplate, graph libraries, and one-off heatmap code. `genomeblocks` unifies those pieces behind a small set of composable objects, grouped into three color-coded families:

<ul class="gb-modules">
  <li class="green">
    <span class="mod-name"><code>Loci</code><span class="mod-tag">Intervals</span></span>
    <small>A list-of-intervals container with set algebra (<code>&amp;</code>, <code>|</code>, <code>-</code>, <code>^</code>), <code>slop</code>, <code>sort</code>, <code>merge</code>, <code>nearest</code>, indexed overlap queries, and signal extraction.</small>
  </li>
  <li class="green">
    <span class="mod-name"><code>Locus</code><span class="mod-tag">Intervals</span></span>
    <small>A single interval with a canonical UID (<code>chrom:start-end(strand)</code>) used everywhere as a stable key.</small>
  </li>
  <li class="green">
    <span class="mod-name"><code>signal</code><span class="mod-tag">Intervals</span></span>
    <small>Threaded bigWig extraction (pybigtools backend, pure-Python fallback), TMM normalization, and comparative heatmaps.</small>
  </li>
  <li class="navy">
    <span class="mod-name"><code>Tags</code><span class="mod-tag">Annotation</span></span>
    <small>An in-memory annotation store for a <code>Loci</code> set, queryable by boolean expression (<code>l.atac &amp; (l.h3k27ac &gt; 1.5)</code>).</small>
  </li>
  <li class="navy">
    <span class="mod-name"><code>Genes</code><span class="mod-tag">Annotation</span></span>
    <small>GENCODE / GTF / UCSC parser yielding <code>Gene</code> → <code>Transcript</code> → <code>Exon</code>/<code>CDS</code>/<code>UTR</code> hierarchies, with enhancer-to-gene assignment in one call.</small>
  </li>
  <li class="navy">
    <span class="mod-name"><code>motifs</code><span class="mod-tag">Annotation</span></span>
    <small>JASPAR motif scanning over a FASTA genome; TF family lookups and enrichment helpers.</small>
  </li>
  <li class="purple">
    <span class="mod-name"><code>Architecture</code><span class="mod-tag">Networks</span></span>
    <small>A chromatin-contact graph (graph-tool) built from BEDPE loops or mcool matrices; spreading, cliques, O/E, hub discovery.</small>
  </li>
  <li class="purple">
    <span class="mod-name"><code>bedpe</code><span class="mod-tag">Networks</span></span>
    <small>BEDPE parsing + pair-to-bed intersection for loops and pairwise-interval data.</small>
  </li>
  <li class="purple">
    <span class="mod-name"><code>browser</code><span class="mod-tag">Networks</span></span>
    <small>An IGV-like, SVG-clean multi-track region viewer built on matplotlib.</small>
  </li>
</ul>

Everything is **chainable**: the output of one stage is always a first-class object accepted by the next.

---

## 60-second example
{: .sec-green }

Three stages, each one a pure object you can hand to the next:

```python
from genomeblocks import Architecture, Genes, Loci

# ① CREs: ATAC peaks, extended ±100 bp, sorted, merged
cre = (Loci.make("atac_peaks.narrowPeak")
           .slop(100)
           .sort()
           .merge())

# Super-enhancers that overlap CREs
se = cre.intersect(Loci.make("H3K27ac_SE.bed"))

# ② Contact graph: CRE-resolved HiChIP loops, Hi-C weights, O/E normalization
arch = (Architecture.make(cre, "RNAP_loops.bedpe", r=2500)
                    .add_mcool(cre, "RNAP.mcool", resolution=5000)
                    .normalize(cre))

# ③ Gene annotation: which annotation class does each super-enhancer fall into?
genes = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)
counts = genes.annotations(se & cre).groupby("annotation").size()
```

<div class="gb-stage-rail">
  <div class="gb-stage green">
    <span class="gb-stage-num">1</span>
    <span><strong>Build intervals.</strong> Peaks → CREs via <code>Loci</code>.</span>
  </div>
  <div class="gb-stage purple">
    <span class="gb-stage-num">2</span>
    <span><strong>Wire the network.</strong> Loops + contacts via <code>Architecture</code>.</span>
  </div>
  <div class="gb-stage navy">
    <span class="gb-stage-num">3</span>
    <span><strong>Annotate.</strong> Gene classes for each CRE via <code>Genes</code>.</span>
  </div>
</div>

---

## Documentation map
{: .sec-navy }

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
| [Tools → Venn diagram](tools/venn) | Drop 2–3 BED / narrowPeak files in the browser and get a live overlap Venn. |
| [Credits](credits) | Upstream tools & citations. |

---

## Citing
{: .sec-navy }

If `genomeblocks` is useful in your work, please cite the GitHub repository:

```
Altintas, U. B. (2024). genomeblocks: Fluent building blocks for regulatory genomics.
https://github.com/birkiy/genomeblocks
```

`genomeblocks` stands on top of several excellent upstream libraries — please also cite the ones whose module was load-bearing in your analysis. The [Credits](credits) page lists them all with BibTeX.

---

## License

MIT. See [LICENSE](https://github.com/birkiy/genomeblocks/blob/main/LICENSE).

<script>
function gbCopyInstall(btn) {
  var text = btn.querySelector('.chip-text').textContent;
  if (navigator.clipboard) { navigator.clipboard.writeText(text); }
  var icon = btn.querySelector('.copy-icon');
  var original = icon.textContent;
  icon.textContent = '✓';
  btn.classList.add('copied');
  setTimeout(function() {
    icon.textContent = original;
    btn.classList.remove('copied');
  }, 1400);
}
</script>
