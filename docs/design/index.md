---
title: Design
layout: default
nav_order: 8
has_children: true
permalink: /design/
---

# Design
{: .no_toc }

How genomeblocks works inside: what each module stores, the path data takes
through it, and why each step costs what it does. Every page draws the
mechanism. The [User Guide]({{ '/guide/' | relative_url }}) and [API Reference]({{ '/api/' | relative_url }}) cover how to
call it, and the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page measures it.
{: .fs-5 .fw-300 }

## The package at a glance
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/overview.svg %}
</div><figcaption>
<strong>Four layers.</strong> <code>Locus</code> and <code>Loci</code> are the core every module builds on. Each domain module handles one data type and imports its heavy dependency only when it is called. Plotting lives in separate <code>*_draw</code> modules, so processing never pulls in matplotlib. <code>genomeblocks.columnar</code> is the table-first foundation added in 1.1: its rows are still <code>Locus</code> objects, and classic functions accept its tables.
</figcaption></figure>

## Principles
{: .sec-navy }

**One key for an interval.** `Locus.uid` (`chr1:100-200(+)`) is how every
table, graph and dictionary refers to a region. Two loci are equal when their
uids are equal, so set operations and lookups never compare floats or
objects.

**Indexes are built once, on first use.** A `Loci` builds its uid map and its
interval index the first time an operation needs them, then keeps them. Repeated
queries against the same set cost only the lookups.

**Methods attach to `Loci`.** `signal`, `bedpe`, `motifs` and `atlas` add their
entry points to `Loci` when imported (`loci.signal(...)`,
`loci.count_pairs(...)`, `loci.scan_motifs(...)`, `loci.enrich(...)`), and
`loci.py` imports them at its end so the methods are always there.

**Heavy dependencies load on first call.** pandas, scipy, matplotlib, cooler,
graph-tool, lightmotif and pybigtools are imported inside the functions that
use them. `import genomeblocks` is a few milliseconds and
`from genomeblocks import Loci` loads numpy and little else.

**Processing and plotting are separate modules.** `signal_draw`,
`motifs_draw`, `architecture_draw` and `browserview` read results; nothing in
a processing module imports matplotlib.

**Chainable results.** Every step returns a first-class object (a `Loci`, a
`Genes`, an `Architecture`, a numpy array or a DataFrame), and the next step
accepts it unchanged.

## Module pages
{: .sec-navy }

| Page | What it explains |
|---|---|
| [Loci & Locus]({{ '/design/loci/' | relative_url }}) | set algebra semantics and the overlap index |
| [Genes]({{ '/design/genes/' | relative_url }}) | GTF → objects, region labels, ATAC-supported isoforms |
| [Signal]({{ '/design/signal/' | relative_url }}) | bigWig → signal cube, backends, process-parallel extraction |
| [Motifs]({{ '/design/motifs/' | relative_url }}) | block scanning of many windows × many motifs |
| [BEDPE & pairs]({{ '/design/bedpe/' | relative_url }}) | the BEDPE reader and streaming Hi-C pair counts |
| [Atlas]({{ '/design/atlas/' | relative_url }}) | one sparse bins × tracks index and vectorised Fisher tests |
| [Architecture]({{ '/design/architecture/' | relative_url }}) | the contact-graph pipeline, weights, genes and hubs |
| [Browser]({{ '/design/browser/' | relative_url }}) | one axis per track, drawer by track type |
| [Columnar]({{ '/design/columnar/' | relative_url }}) | tables aligned by row, the sorted edge table, short-range HiChIP |
