---
title: Concepts
layout: default
nav_order: 4
---

# Concepts
{: .no_toc }

The mental model behind `genomeblocks` 2.0. Read this before the module guides — it explains *why* the API looks the way it does.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Everything is a table; the row is the join key

A `Loci` is four numpy columns on a `Genome`, plus a dict of extra columns:

```
codes    int32   chromosome code into the table's Genome
starts   int64   0-based start
ends     int64   end, exclusive
strands  int8    0 '.', 1 '+', 2 '-'
cols     dict    name, score, signalValue, ... — one array per column, aligned to the rows
```

Everything computed *from* a `Loci` has one row per locus, in the same order: the signal cube `S[i]`, the motif matrix `M.iloc[i]`, the annotation label `labels[i]`, the vertex `i` of an `Architecture`. There is no uid → position dictionary to keep in sync; a boolean mask or a row array selects the same rows everywhere:

```python
import numpy as np
from genomeblocks import Loci

L = Loci.make("atac_peaks.narrowPeak", keep=["signalValue"])
S = L.signal("ATAC.bw", n_bins=10, flank=500)          # (7, 1, 10)
open_rows = S[:, 0, :].mean(axis=1) > 5                 # a mask over the rows
L[open_rows]                                            # -> Loci(n=4, chroms=1, sorted, cols=[signalValue])
L["mean_atac"] = S[:, 0, :].mean(axis=1)                # a new column, one value per row
L[0].mean_atac                                          # -> 6.55...
```

The other tables are built from `Loci`:

| Table | Rows | Linked by |
|---|---|---|
| `Genes` | `genes`, `transcripts`, `features` — three `Loci` | `transcripts['gene'][k]` is the gene row of transcript `k`; `features['transcript'][j]` the transcript row of feature `j` |
| `Pairs` | anchors `a` and `b` — two row-aligned `Loci` | pair `i` is `a[i]` and `b[i]` |
| `Architecture` | vertices = the `Loci`; edges = `src`, `tgt` + `ep` columns | `vp` columns are aligned to the `Loci` rows, `ep` columns to the edges |
| `Atlas` | a sparse bins × tracks matrix | tracks are rows of the track table |

`L[i]` hands out a `LocusView` — a real `Locus` whose fields read (and write through to) the columns; `genes['TP53']` a `GeneView` with `.tss`, `.transcripts`, `.exons`. Whole-set work never iterates over views.

---

## One `Genome` per table, re-coded on contact

Chromosomes are stored as small integer codes relative to the table's `Genome` (names ↔ codes, plus sizes when known). Tables made from one another — `take`, `slop`, `merge`, an `Architecture`'s vertices, a `Genes`' three tables — share their `Genome`. Tables read separately each get their own, and an operation between two of them re-codes the right-hand side onto the left through a small lookup (`O(n)`), adding unseen names to the left `Genome`:

```python
from genomeblocks import as_loci, Genome

a = as_loci([("chr2", 10, 20), ("chr1", 5, 9)])
b = as_loci([("chr1", 1, 6), ("chrX", 1, 6)])
a.genome.names, b.genome.names        # -> ['chr2', 'chr1'], ['chr1', 'chrX']
a & b                                 # -> Loci(n=1, chroms=1)
a.genome.names                        # -> ['chr2', 'chr1', 'chrX']   b was re-coded onto a's Genome
```

Pass one `Genome` to several constructors to skip even that:

```python
g = Genome.from_sizes("hg38.chrom.sizes")      # or Genome.from_fasta(path), Genome(sizes={...})
x = Loci.make("atac_peaks.bed", genome=g)
y = as_loci("chr1:0-100", genome=g)
x.genome is y.genome                           # -> True
```

Codes are append-only and never change. Sorting uses a separate natural-order rank (chr1, chr2, …, chr10, …, chrX), so `sort()` gives the same order whatever the order in which names were first seen. On a sorted table each chromosome is one contiguous block: `L.chrom_offsets` → `{'chr1': (0, 5), 'chr2': (5, 7)}`, and `L.by_chrom('chr1')` is a slice.

{: .note }
> There is no session-wide `Genome`. Nothing has to be set up before the first call, and two analyses in one process cannot interfere.

---

## Coordinates: 0-based, half-open, everywhere

Every table stores `[start, end)` in BED convention — including `Genes`. A GTF or GFF3 1-based start becomes `start - 1` while parsing; UCSC genePred tables are already 0-based; frames from pandas / polars / bioframe / PyRanges are read as 0-based (`Genes.from_frame(one_based=...)` says otherwise). So a gene and a peak compare directly, and the first base of a chromosome is `[0, 1)`.

```python
from genomeblocks import Genes

genes = Genes.make("genes.gtf")                 # GTF: chr1 ... gene 1001 5000 . + ...
genes.to_pandas()[["gene_name", "chrom", "start", "end", "strand"]]
#   gene_name chrom  start    end strand
# 0    GENE_A  chr1   1000   5000      +
# 1    GENE_B  chr1  10000  11000      -
```

**The TSS rule.** A gene's TSS is the 1-bp interval `[t, t+1)` on the gene's strand, with `t = start` on `+` and `t = end - 1` on `-`:

```python
genes["GENE_A"].tss        # -> Locus(chrom='chr1', start=1000,  end=1001,  strand='+')
genes["GENE_B"].tss        # -> Locus(chrom='chr1', start=10999, end=11000, strand='-')
genes.get_tss()            # one row per gene, with gene_name / gene_id / gene (row) columns
```

Promoter windows are `get_tss().slop(r)`; `annotations()` labels Promoter-TSS > 5UTR > 3UTR > Exonic > Intronic > Intergenic from that same index. Two intervals overlap when `s1 < e2 and s2 < e1`; a zero-length interval `[p, p)` is a point between two bases and overlaps `[s, e)` when `s < p < e`. Every interval backend is held to these rules.

---

## The boundary: `as_loci`

`gb.as_loci(x)` is the one entry point, and every public function calls it on its inputs — `signal(peaks_df, ...)`, `genes.annotations("peaks.bed")`, `Architecture.make(cre, loops_frame)` all work without converting first. It accepts:

| Input | Example |
|---|---|
| a `Loci` (returned as is), a `Locus` | `as_loci(cre)` |
| a region string, a list of regions or tuples, a `(chroms, starts, ends)` tuple of arrays | `as_loci("chr1:1,000-2,000")`, `as_loci([("chr1", 10, 20), "chr2:5-9"])` |
| a path: BED / narrowPeak / broadPeak / bedGraph, CSV / TSV with a header, parquet (`.gz` too) | `as_loci("peaks.narrowPeak")` — rows stay in file order |
| any table: pandas, polars (eager or lazy), pyarrow, bioframe, PyRanges, pybedtools, dicts of columns, structured arrays, anything with `__arrow_c_stream__` / `__dataframe__` | `as_loci(pl.DataFrame({"chr": [...], "start": [...], "end": [...]}))` |
| an AnnData (its `var` regions, from columns or from names like `chr1:100-200`) | `as_loci(adata)` |

Columns are found by name — `chrom` / `chr` / `Chromosome` / `seqnames` / `#chrom`, `start` / `Start` / `chromStart`, `end` / `End` / `chromEnd`, `strand` — and positionally only for a header-less text + int + int frame. Anything else raises an error that names the accepted spellings and the fix (`chrom=`, `start=`, `end=`):

```python
import pandas as pd
as_loci(pd.DataFrame({"a": [1], "b": [2], "c": [3]}))
# ValueError: expected columns for the chromosome, start and end, e.g. chrom/start/end (bioframe),
#   Chromosome/Start/End (pyranges), seqnames/start/end or chr/chromStart/chromEnd — or pass chrom=, start=, end= to name them
```

Normalising once at the boundary means no function branches on input type inside. See [Interoperability]({{ '/interoperability/' | relative_url }}).

---

## The exits: `to_*`

Every container has explicit converters, and they say when they copy:

| Method | Returns | Copies? |
|---|---|---|
| `to_pandas()` | DataFrame — `chrom` categorical in genome order, `start`, `end`, `strand`, then the extra columns | the integer columns are shared |
| `to_arrow()` / `to_polars()` | pyarrow Table / polars DataFrame | numeric columns zero-copy through Arrow |
| `to_bioframe()`, `to_pyranges()`, `to_bedtool()`, `to_cgranges()` | that library's object | yes (their layout) |
| `to_numpy()` / `np.asarray(L)` | a structured array `(chrom, start, end, strand, ...)` | yes |
| `to_records()` | `(chrom, start, end, strand)` tuples | yes |
| `to_bed(path)`, `save(path)` (parquet), `to_fasta(path, fasta)` | a file | — |
| `to_anndata(X)` | AnnData with the loci as `var` | — |

`Genes` adds `to_pandas(table)`, `to_gtf`, `to_bed12`; `Pairs` adds `to_bedpe`; `Architecture` adds `to_networkx`, `to_igraph`, `to_graph_tool`, `to_scipy`, `to_anndata`, `edges_frame`, `vertices_frame`; a signal cube leaves through `interop.cube_to_xarray` / `cube_to_pandas` / `cube_to_anndata`. Round trips are exact: `Loci.from_frame(L.to_polars()).equals(L)` is `True` for every supported format.

---

## Protocols: other libraries take the tables as they are

Every container (`Loci`, `Genes`, `Pairs`, `Architecture`, `Atlas`, `Library`) implements, on top of its `to_arrow()`:

| Protocol | Who reads it |
|---|---|
| `__arrow_c_stream__` (Arrow PyCapsule) | polars, pyarrow, DuckDB, pandas ≥ 2.2 |
| `__dataframe__` (dataframe interchange) | seaborn, pandas |
| `__narwhals_dataframe__` | altair, plotly, marimo and other narwhals-based tools |
| `__len__`, `__iter__`, `__getitem__`, `__array__`, `shape`, `columns` | Python, numpy |

```python
import polars as pl, duckdb, seaborn as sns, altair as alt, plotly.express as px

pl.DataFrame(L)                                        # a polars frame
duckdb.sql("select count(*) from L").fetchall()        # -> [(7,)]
sns.histplot(data=L, x="mean_atac")
alt.Chart(L).mark_bar().encode(x="chrom:N", y="count()")
px.histogram(L, x="mean_atac")
```

For notebooks: a short `repr`, an HTML head / summary, `head()`, `tail()`, `describe()`.

---

## Backends: explicit choice, no silent fallback

The heavy work on the tables runs through a *backend*, one per family of operations:

| Family | Default | Others | Operations |
|---|---|---|---|
| `intervals` | genomeblocks (numpy) | cgranges, ncls, bioframe, pyranges, bedtools | overlap, nearest, merge, point lookups (cgranges and ncls: overlap and point lookups only) |
| `bigwig` | pybigtools | pybigwig, python | reading signal |
| `motifs` | lightmotif | moods, biopython | scoring matrices along sequences |
| `fasta` | genomeblocks (indexed) | pysam, pyfaidx, memory, biopython | fetching sequence |
| `tables` | polars, else pandas | pandas | parsing BED / GTF / pairs text |
| `graph` | graph-tool, else scipy | igraph, networkx | components, PageRank, layouts |

Three rules:

1. **Automatic only where the answers are identical.** polars and pandas parse to the same columns; graph-tool and scipy give the same components; pybigtools and the pure-Python reader bin alike. Anything else is a request.
2. **A request is per call or per block**, never global state:
   ```python
   import genomeblocks as gb
   L.merge(backend="pyranges")                         # this call
   with gb.use_backend(intervals="bioframe", graph="networkx"):
       ...                                             # this block; names are checked on entry
   ```
3. **A missing or unknown backend raises** with the install command — never a silent switch:
   ```python
   L.merge(backend="cgranges")
   # ImportError: the 'cgranges' intervals backend is not installed: conda install -c bioconda cgranges  (...)
   L.merge(backend="ncls")
   # NotImplementedError: the 'ncls' intervals backend has no merge; use one of: genomeblocks, bioframe, pyranges, bedtools
   ```

Backend differences do not leak into results: row order, nulls and dtypes are normalised to genomeblocks' rules whatever engine ran, and the test suite runs every public function against every installed backend. `gb.backends()` shows what is installed, the automatic default, and what a call would use right now. See [Backends]({{ '/backends/' | relative_url }}) and [Benchmarks]({{ '/benchmarks/' | relative_url }}).

---

## No hidden global state

Nothing in `genomeblocks` has to be configured before the first call, and nothing a call does changes a later call elsewhere:

- no session-wide `Genome` — each table carries its own (see above);
- no global backend switch — `backend=` and `use_backend` are scoped;
- indexes (point-lookup trees, the uid → row dict, chromosome offsets, a `Genes`' annotation index) are built on first use, cached on the object, and dropped when its rows change; pickling drops them too, so nothing stale is ever serialised.

---

## Lazy imports

`import genomeblocks` takes a few milliseconds. Public names are loaded from `__init__.py` through `__getattr__` on first use, and heavy optional dependencies are imported inside the functions that need them:

```python
import sys, genomeblocks
from genomeblocks import Loci
"pandas" in sys.modules, "matplotlib" in sys.modules, "cooler" in sys.modules
# -> (False, False, False)
```

So a script that only does interval algebra never pays for matplotlib or cooler, and an optional library is required only at the moment its backend or converter is asked for.

---

## Chainable pipelines and set algebra

Operations return a table (the same type, or the natural next one) — never `None`, never "mutate and also return something else" — so stages chain, and intermediate state stays out of your namespace:

```python
import genomeblocks as gb

cre = Loci.make("peaks.narrowPeak").slop(100).merge()
A = (gb.Architecture.make(cre, "loops.bedpe", r=2500)
       .add_mcool("hic.mcool", resolution=5000)
       .normalize()
       .annotate(genes)
       .strength())
```

The in-place methods of `Architecture` (`add_mcool`, `normalize`, `annotate`, `strength`, `prune`) add columns to the same graph and return it. Python operators map to set operations:

| Operator | `Loci` | `Architecture` |
|---|---|---|
| `a & b` | rows of `a` overlapping `b` (`intersect`) | vertex + edge intersection |
| `a \| b` | concatenation (`union`, `+`) | vertex + edge union |
| `a - b` | rows of `a` not overlapping `b` (`difference`) | — |
| `a ^ b` | rows of either not overlapping the other | — |

Because they take any interval input, a CRE set composes directly: `active_distal = (atac & h3k27ac) - genes.get_tss().slop(1000)`.
