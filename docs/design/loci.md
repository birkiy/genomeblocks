---
title: Loci & Locus
parent: Design
layout: default
nav_order: 1
---

# Loci & Locus
{: .no_toc }

`Loci` is a table of intervals: four numpy columns plus any extra columns, all
aligned by row. `Locus` is what one row looks like. Everything else in the
package is built on these two, and every whole-set operation on a `Loci` runs
through the `intervals` backend.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## What a Loci stores
{: .sec-green }

| Column | dtype | Meaning |
|---|---|---|
| `codes` | `int32` | chromosome code into the table's `Genome` |
| `starts` | `int64` | 0-based start (BED convention) |
| `ends` | `int64` | end, exclusive |
| `strands` | `int8` | `0` = `.`, `1` = `+`, `2` = `-` |
| `cols` | `dict[str, ndarray]` | any other column (`name`, `score`, `gene_name`, ...), one value per row |

```python
import genomeblocks as gb

L = gb.Loci.make("peaks.bed", keep=True)
L
# -> Loci(n=7, chroms=2, sorted, cols=[name, score])
L.codes.dtype, L.starts.dtype, L.ends.dtype, L.strands.dtype
# -> (dtype('int32'), dtype('int64'), dtype('int64'), dtype('int8'))
L.codes, L.strands
# -> (array([0, 0, 0, 0, 0, 1, 1], dtype=int32), array([1, 2, 1, 0, 2, 1, 1], dtype=int8))
L.genome
# -> Genome(2 chroms: chr1, chr2)
L.shape, L.columns
# -> ((7, 6), ['chrom', 'start', 'end', 'strand', 'name', 'score'])
```

**The row number is the join key.** A signal cube, a motif matrix, an
annotation vector or an `Architecture` vertex column computed from a `Loci`
has one row per interval, in this order. Joining a result back to its
intervals is array indexing; there is no dictionary in between.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-rows.svg %}
</div><figcaption>
<strong>Row <em>i</em> everywhere.</strong> Row <em>i</em> of the CREs is row <em>i</em> of every label, signal and graph column, and an Architecture's edge table stores row numbers, not uids. Genes are three tables linked the same way: transcripts point to their gene's row, features to their transcript's row.
</figcaption></figure>

**Chromosomes are codes.** Each table owns a `Genome`, a `names ↔ codes`
dictionary (plus sizes when known). Codes are append-only, so a code never
changes once given, and sorting uses a separate natural order (chr1, chr2, ...,
chr10, ..., chrX) that does not depend on the order names were seen in. There
is no session-wide genome: two tables read separately have two `Genome`s, and
an operation between them re-codes the second onto the first through a small
lookup table (`Loci._check`, O(n)). Passing one `genome=` to several
constructors skips even that.

```python
a = gb.Loci.make("peaks.bed")
b = gb.as_loci("chr2:400-700")
a.genome.names, b.genome.names
# -> (['chr1', 'chr2'], ['chr2'])
(a & b).to_records()
# -> [('chr2', 500, 600, '+')]
```

**Sorted rows are chromosome blocks.** `Loci.make` sorts into genome order, so
each chromosome is one contiguous run of rows; `chrom_offsets` gives the
`(first, last + 1)` of each block and `by_chrom` is a slice. `take`,
`from_frame` and `as_loci` keep the input order instead, because a frame's
rows, an AnnData's `var` and the resulting Loci must line up one to one.

```python
L.is_sorted, L.chrom_offsets
# -> (True, {'chr1': (0, 5), 'chr2': (5, 7)})
```

## Locus and LocusView
{: .sec-green }

A `Locus` is a small dataclass: `chrom`, `start`, `end` (half-open) and
`strand` (`.` by default). `L[i]` does not copy anything: it returns a
`LocusView`, a real `Locus` whose fields read and write the columns, with the
extra columns one attribute away.

| Attribute | Value | Used for |
|---|---|---|
| `uid` | `"chr1:100-200(+)"` | equality, hashing, `L[uid]`, `Architecture[uid]`, index of motif matrices |
| `center` | `(start + end) // 2` | signal windows, motif windows, loop distances |
| `row` | position in the table (`LocusView` only) | going back to the columns |
| ordering | by `chrom`, then `start` | `sorted()` on plain `Locus` objects |

```python
L[1]
# -> Locus[1](chr1:1900-2100(-))
L[1].score, L[1].row
# -> (20.0, 1)
```

`Locus.parse("chr8:127.7-128.1 Mb")` and `Locus.from_uid` read the two string
spellings; `L.uid` and `L.uids` (`uid → row`) are built once per table and
dropped when a coordinate changes.

## The boundary: anything in, one representation inside
{: .sec-green }

Every public function calls `as_loci()` on its interval inputs, so a pandas or
polars frame, a pyarrow table, a bioframe / pyranges / pybedtools object, a
dict of columns, a structured array, an AnnData, a BED / CSV / parquet path, a
region string, a `Locus` or a list of any of these is accepted wherever a
`Loci` is. Normalisation happens once, at the boundary, through narwhals
(`interop.frame`), and the code inside never branches on the input type.

```python
import pandas as pd

df = pd.DataFrame({"Chromosome": ["chr1", "chr2"], "Start": [1000, 0], "End": [1050, 10_000]})
hits = L & df                      # a frame is fine where a Loci is expected
hits.to_records()
# -> [('chr1', 900, 1100, '+'), ('chr2', 500, 600, '+'), ('chr2', 5000, 5100, '+')]
hits.cols["name"]                  # rows of L come back whole, with their columns
# -> array(['p1', 'p6', 'p7'], dtype=object)
```

`Loci.from_frame` finds the coordinate columns by name, case-insensitively
(`chrom` / `chromosome` / `chr` / `seqnames` / `Chromosome` ..., `start` /
`chromStart` / `Start`, `end` / `chromEnd` / `End`, `strand`), falls back to
position only for a header-less text + int + int frame, and otherwise raises
naming the accepted spellings. Nulls in a coordinate column are an error, not
a dropped row. Every other column is kept (`keep=True`) unless a list says
which.

{: .note }
Column names are lenient on the way in and strict inside: once a frame is a
`Loci`, the columns are `codes`, `starts`, `ends`, `strands` and nothing else
is consulted.

## Set algebra keeps whole rows
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-setops.svg %}
</div><figcaption>
<strong>Operators filter rows; they do not cut intervals.</strong> <code>A &amp; B</code> keeps every row of A that touches B (like <code>bedtools intersect -u</code>), <code>A − B</code> keeps the rest, and <code>merge()</code> fuses overlapping or book-ended intervals after sorting. A row comes back with its own columns, so the result still indexes any table keyed on A.
</figcaption></figure>

| Expression | Result |
|---|---|
| `A & B`, `A.intersect(B)` | rows of A overlapping at least one row of B |
| `A - B`, `A / B`, `A.difference(B)` | rows of A overlapping nothing in B |
| `A + B`, `A | B`, `A.union(B)` | concatenation, columns both sides have kept (follow with `.sort().merge()` to fuse) |
| `A ^ B` | `(A - B) + (B - A)` |
| `A.overlap_any(B)` | boolean mask over A |
| `A.overlap_pairs(B)` | `(rows of A, rows of B)` for every overlapping pair, sorted |
| `A.nearest(B)` | `(row in B, distance)` per row of A; overlaps are 0, book-ended are 0, nothing on the chromosome is `(-1, -1)` |
| `A.slop(n)` | every row widened by `n` bp (start clamped at 0) |
| `A.sort()`, `A.merge()` | genome order; fused overlapping / book-ended rows, strand of the first |
| `A.overlaps(region)`, `A.overlap_rows(region)` | one window: the rows (or row numbers) it touches |

Two intervals overlap when `s1 < e2 and s2 < e1`. A zero-length interval
`[p, p)` is the point between two bases and overlaps `[s, e)` exactly when
`s < p < e`; every backend's answer is filtered to this rule, so the engines
agree on it.

```python
other = gb.as_loci([("chr1", 1000, 1050), "chr2:0-10,000"])
L.overlap_pairs(other)
# -> (array([0, 5, 6]), array([0, 1, 1]))
L.nearest(other)
# -> (array([0, 0, 0, 0, 0, 1, 1]), array([   0,  850, 3850, 8900, 9850,    0,    0]))
```

## The numpy kernels
{: .sec-green }

The default engine is `genomeblocks._intervals`: a handful of numpy calls over
the whole genome, with no per-chromosome or per-interval Python loop. Every
chromosome is laid end to end on one axis, `code · 2⁴⁰ + position`
(`STRIDE = 1 << 40`, larger than any chromosome), so two tables that share a
`Genome` compare directly on that axis.

| Kernel | How |
|---|---|
| `overlaps_any` | sort the reference by start, keep a running maximum of ends; query `[s, e)` overlaps iff the largest end among references starting before `e` exceeds `s`. Works for references of any length, O((n + m) log m). |
| `overlap_pairs` | references sorted once; a reference can only reach `[s, e)` if it starts in `(s − maxlen, e)`, so each query's candidates are one `searchsorted` range, expanded with `repeat` and filtered by `end > s`. |
| `merge` | sort by start, running maximum of ends; a new block starts where a start exceeds the running max. Returns the first input row of each block too. |
| `nearest` | left neighbour from the running max and its argmax, right neighbour from the first start at or after `e`, overlaps from `overlap_pairs` (lowest start wins); the closer side wins, the left one on a tie. |

The `maxlen` window in `overlap_pairs` is tight for peak-like data. A set
that mixes peaks with megabase domains widens every window to the longest
domain; `overlaps_any`, `merge` and `nearest` do not have that dependence,
and for pair enumeration on such a set `backend="ncls"` or `"bioframe"` is the
alternative.

## Single-window lookups and the index cache
{: .sec-green }

Browsers, `Architecture.near`, `Genes` lookups and `L.overlaps(region)` ask
for one window at a time. That goes through `overlap_rows`, which keeps a
lookup index on the `Loci`: built on first use, one per backend, cached until
the rows change.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-index.svg %}
</div><figcaption>
<strong>Two <code>searchsorted</code> calls bound the scan.</strong> The default <code>PointIndex</code> keeps the rows sorted by start together with the running maximum of their ends. No row before the first position whose running max passes <code>qs</code> can reach the query, and no row starting at or after <code>qe</code> can either, so only <code>lo:hi</code> is tested for <code>end &gt; qs</code>. A lookup costs O(log n + k). The other engines answer the same window with their own structure (a cgranges or NCLS tree cached on the Loci, or a one-row <code>overlap_pairs</code> for bioframe / pyranges / bedtools) and come back as the same sorted row numbers.
</figcaption></figure>

```python
L.overlap_rows("chr1:1,000-2,000")
# -> array([0, 1])
list(L._indexes), type(L._indexes["genomeblocks"]).__name__
# -> (['genomeblocks'], 'PointIndex')
L.overlaps("chr1:1,000-2,000", backend="ncls").to_records()
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '-')]
list(L._indexes)                   # one index per engine
# -> ['genomeblocks', 'ncls']
```

Adding or replacing a column (`L["score"] = ...`) keeps the indexes: only a
coordinate write through a `LocusView`, or a new `Loci`, invalidates them
(`_dirty`), together with the cached `uid` array and `uids` map.

## The backend seam
{: .sec-green }

Whole-set operations dispatch to `genomeblocks.backends.intervals`, which
takes two `Loci` on one `Genome` and returns row numbers. The caller never
sees which engine ran, and the engines' answers are normalised to one set of
rules before they come back:

- pairs are filtered to the half-open rule and sorted by `(query row,
  reference row)`;
- `nearest` distances are recomputed as the gap in bases from the engine's
  chosen row, overlaps are re-ranked by the genomeblocks tie rule (lowest
  start, then lowest row), and -1 marks a chromosome with nothing;
- `merge` results are re-sorted and the first input row of every block is
  recovered with one `searchsorted`;
- chromosome codes, not names, are handed to the engines as labels, so both
  sides agree without decoding and names never clash.

| Operation | genomeblocks | cgranges | ncls | bioframe | pyranges | bedtools |
|---|---|---|---|---|---|---|
| overlap (`&`, `-`, `overlap_pairs`) | yes | yes | yes | yes | yes | yes |
| `nearest` | yes | – | – | yes | yes | yes |
| `merge` | yes | – | – | yes | yes | yes |
| point lookup | yes | yes | yes | yes | yes | yes |

The default is always `genomeblocks`; there is no automatic switch to another
engine, because the others exist for comparison and for users who already
depend on them. `backend=` on a call, or `with gb.use_backend(intervals=...)`
for a block, chooses one; a missing engine raises with its install command,
and an operation an engine lacks raises `NotImplementedError` naming the ones
that have it.

```python
for b in ("genomeblocks", "ncls", "bioframe", "pyranges"):
    qi, ri = L.overlap_pairs(other, backend=b)
    print(b, qi.tolist(), ri.tolist())
# -> genomeblocks [0, 5, 6] [0, 1, 1]
# -> ncls [0, 5, 6] [0, 1, 1]
# -> bioframe [0, 5, 6] [0, 1, 1]
# -> pyranges [0, 5, 6] [0, 1, 1]

with gb.use_backend(intervals="bioframe"):
    len(L & other), L.nearest(other)[1].tolist()
# -> (3, [0, 850, 3850, 8900, 9850, 0, 0])

L.overlap_pairs(other, backend="cgranges")
# -> ImportError: the 'cgranges' intervals backend is not installed: conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)
L.nearest(other, backend="ncls")
# -> NotImplementedError: the 'ncls' intervals backend has no nearest; use one of: genomeblocks, bioframe, pyranges, bedtools
```

{: .warning }
bedtools widens zero-length intervals on its side before intersecting, so it
can still miss a pair that needs the widened base; every other engine matches
the default exactly. The parity tests in `tests/test_loci.py` run every
installed engine against brute force.

## Out: converters and protocols
{: .sec-green }

Every `Loci` converts back out, and the protocols other libraries look for are
there, so polars, pyarrow, DuckDB, pandas, seaborn, plotly and altair take the
object as it is.

| Call | Gives | Copies? |
|---|---|---|
| `to_arrow()` | pyarrow Table: `chrom` and `strand` as dictionary arrays, `start` / `end` int64, then `cols` | numeric columns are wrapped, not copied; object columns are converted |
| `to_polars()` | polars DataFrame, through Arrow | numeric columns not copied |
| `to_pandas(uid=False)` | DataFrame with categorical `chrom` / `strand` (bioframe's column names) | builds the categoricals; numeric columns are shared |
| `to_bioframe()` | plain `chrom` / `start` / `end` / `strand` strings | decodes names |
| `to_pyranges()`, `to_bedtool()`, `to_cgranges()` | the library's own object | yes (their formats) |
| `to_anndata(X, obs=)` | AnnData with the loci as `var` (names `chrom:start-end`) | `var` built from `to_pandas` |
| `to_records()`, `to_numpy()`, `to_bed(path)` | tuples, a structured array, BED6 text | yes |
| `save(path)` / `Loci.load(path)` | parquet via Arrow, read back on any `Genome` | file |

```python
import polars as pl, pyarrow as pa, duckdb

pl.DataFrame(L).shape                 # __arrow_c_stream__
# -> (7, 6)
pa.table(L).schema.field("chrom").type
# -> DictionaryType(dictionary<values=string, indices=int32, ordered=0>)
duckdb.sql("select chrom, count(*) n from L group by chrom order by chrom").fetchall()
# -> [('chr1', 5), ('chr2', 2)]
L.save("peaks.parquet"); gb.Loci.load("peaks.parquet").equals(L, cols=True)
# -> True
```

`_table.TableMixin` derives `__arrow_c_stream__`, `__dataframe__`,
`__narwhals_dataframe__`, `to_polars` and `shape` from `to_arrow()` alone, so
`Genes`, `Pairs`, `Architecture`, `Atlas` and the motif `Library` get the same
protocols by implementing one method. `__narwhals_dataframe__` is what makes
altair accept the object, since altair cannot consume interchange-only
objects.

## Costs
{: .sec-green }

- **Reading:** `Loci.make` parses the file through the `tables` backend
  (polars when installed, else pandas, identical columns), encodes chromosome
  names once per distinct name, and sorts once.
- **`A & B`:** one sort of B plus two `searchsorted` passes over A, all in
  numpy; nothing is copied for the mask, and `take` copies only the kept rows.
- **Point lookups:** one index build per `Loci` and engine, then O(log n + k)
  per window.
- **Memory:** 21 bytes per interval for the four columns (`4 + 8 + 8 + 1`),
  plus the extra columns; `_repr_html_` shows the figure.
- **Converters:** Arrow and polars share the numeric buffers; pandas builds
  two categoricals; the library-specific exports copy.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
