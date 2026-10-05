---
title: Design
layout: default
nav_order: 8
has_children: true
has_toc: false
permalink: /design/
---

# Design
{: .no_toc }

How genomeblocks works inside: what each table stores, the path data takes
through it, and why each step costs what it does. Every page draws the
mechanism. The [User Guide]({{ '/guide/' | relative_url }}) and [API Reference]({{ '/api/' | relative_url }}) cover how to
call it, and the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page measures it.
{: .fs-5 .fw-300 }

## The package at a glance
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/overview.svg %}
</div><figcaption>
<strong>Five bands, four layers.</strong> The <em>core</em> tables are the package: <code>Genome</code> maps names to codes, <code>Loci</code> holds intervals as numpy columns, and a <code>Locus</code> is one row read from them. The <em>domain</em> tables (<code>Genes</code>, <code>Pairs</code>, <code>Architecture</code>, <code>Atlas</code>) and the compute modules (<code>signal</code>, <code>motifs</code>) are built on those columns and return arrays aligned to the rows. Every input enters through the <em>boundary</em> (<code>genomeblocks.interop</code>), so the tables never see pandas or polars. Heavy work runs through one of six <em>backend</em> families, invisible unless you ask for a particular engine. <em>Plot</em> modules read the tables and are never imported by them, so processing never pulls in matplotlib.
</figcaption></figure>

```python
import genomeblocks as gb                          # 3 ms; imports no submodule

cre   = gb.Loci.make("atac.narrowPeak", keep=True) # tables
genes = gb.Genes.make("gencode.gtf")
A = (gb.Architecture.make(cre, "loops.bedpe")      # boundary: a path, a Pairs or a frame
       .add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes))
A.components(backend="networkx")                   # backends: same labels on any engine
A.draw("chr8:127.5-128.5 Mb")                      # drawing
```

## Principles
{: .sec-navy }

**The row number is the join key.** Row *i* of a `Loci` is row *i* of the
signal cube, the motif matrix, the label vector and every `Architecture`
vertex column computed from it; `transcripts['gene'][k]` is the row of
transcript *k*'s gene; an edge is `(src row, tgt row)`. Joining a result back
to its intervals is array indexing. The uid (`chr1:100-200(+)`) remains as a
name for lookups and exports, not as the key the tables rely on.

**One Genome per table, re-coded on contact.** Chromosomes are stored as small
integer codes into the table's `Genome`. Tables built from one another share
it; tables read separately each get their own, and an operation between two
of them re-codes the second onto the first through a small lookup
(`Loci._check`). There is no session-wide dictionary and no required setup.

**Boundary first.** Every public function calls `as_loci` on its interval
inputs and the tables convert back out through explicit `to_*` methods and
the Arrow / dataframe protocols. Column names are lenient on input, strict
internally. Nothing deep in the package branches on an input type.

**Backends are invisible.** A call runs on the default engine; `backend=` or
`gb.use_backend()` picks another; results are normalised to one set of rules
so the engine never shows in them. Automatic choice happens only where the
candidates are known to give identical answers. A requested engine that is
missing raises with its install command rather than switching silently.

**Protocols, not adapters.** Every container implements
`__arrow_c_stream__`, `__dataframe__` and `__narwhals_dataframe__` (and `Loci`
`__array__`), so polars, DuckDB, seaborn, plotly and altair take it as it is.
`shape`, `columns`, `head()`, `describe()` and `_repr_html_` make it behave
like a frame in a notebook.

**Lazy imports.** `import genomeblocks` imports nothing but the names table
and resolves `gb.Loci`, `gb.backends` and the rest on first use. pandas,
pyarrow, narwhals, scipy, matplotlib and every engine are imported inside the
function that needs them, so the import costs milliseconds and a pipeline
pays only for what it calls.

**Processing and drawing are separate modules.** `signal_draw`,
`motifs_draw`, `architecture_draw`, `browserview`, `view` and `igv` read
results; nothing in a processing module imports matplotlib.

**Chainable results.** Every step returns a table, a numpy array or a frame,
and the next step accepts it unchanged.

## Module pages
{: .sec-navy }

| Page | What it explains |
|---|---|
| [Tables]({{ '/design/tables/' | relative_url }}) | columns and codes, the shared Genome, the sorted edge table and its views, parquet, short-range HiChIP and super-enhancers |
| [Loci & Locus]({{ '/design/loci/' | relative_url }}) | set algebra on the genome-wide axis, point lookups, the interval backends |
| [Genes]({{ '/design/genes/' | relative_url }}) | GTF / GFF3 → three linked tables, 0-based coordinates, region labels, supported isoforms |
| [Signal]({{ '/design/signal/' | relative_url }}) | bigWig → signal cube, the bigWig backends, process-parallel extraction |
| [Motifs]({{ '/design/motifs/' | relative_url }}) | one library, block scanning of many windows × many motifs, three engines |
| [BEDPE & pairs]({{ '/design/bedpe/' | relative_url }}) | the `Pairs` table and streaming Hi-C pair counts |
| [Atlas]({{ '/design/atlas/' | relative_url }}) | one sparse bins × tracks index and vectorised Fisher tests |
| [Architecture]({{ '/design/architecture/' | relative_url }}) | the contact-graph pipeline, weights, genes, hubs and the graph backends |
| [Browser]({{ '/design/browser/' | relative_url }}) | one axis per track, a drawer per track type |

The boundary and the engines have their own pages outside this section:
[Interoperability]({{ '/interoperability/' | relative_url }}) and
[Backends]({{ '/backends/' | relative_url }}).
