---
title: loci
parent: API Reference
layout: default
nav_order: 2
---

# `genomeblocks.loci`
{: .no_toc }

Genomic intervals as a table of numpy columns on a shared
[`Genome`](#genome): `codes` (int32 chromosome code), `starts` (int64,
0-based), `ends` (int64, exclusive), `strands` (int8: 0 `.`, 1 `+`, 2 `-`) and
`cols`, a dict of extra columns aligned to the rows. The row number is the
join key — a signal cube, a motif matrix, an annotation vector or an
`Architecture` vertex column computed from a `Loci` has one row per locus in
this order. See the [Loci guide]({{ '/guide/loci/' | relative_url }}) for the
workflow and [Design: loci]({{ '/design/loci/' | relative_url }}) for the
kernels behind the interval operations.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
from genomeblocks import Loci, Genome, as_loci
```

Every method that takes another interval set (`overlap_any`, `nearest`,
`&`, ...) accepts anything [`as_loci`]({{ '/api/interop/' | relative_url }})
takes — a path, a frame, a region string, a list of tuples — and re-codes it
onto this table's `Genome`. Methods with a `backend=` keyword run on the
[intervals backend]({{ '/api/backends/' | relative_url }}) you name
(`'genomeblocks'` by default; `'cgranges'`, `'ncls'`, `'bioframe'`,
`'pyranges'`, `'bedtools'`), with identical results.

The examples use a seven-row BED6 (`peaks.bed`, chr1 and chr2, columns
`name` and `score`) and a 20 kb / 8 kb synthetic genome.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Loci`
{: .sec-navy }

```python
Loci(codes=(), starts=(), ends=(), strands=None, *, genome=None, cols=None, filename=None, is_sorted=None)
```

The constructor takes the four columns as arrays (`strands=None` means all
`.`) plus `cols={name: array}`; every array must have one value per row.
You rarely call it: the constructors below and `as_loci` build a `Loci` from
files, frames and records.

```python
Loci([0, 0], [5, 50], [9, 60], genome=Genome(["chr1"])).to_records()
# -> [('chr1', 5, 9, '.'), ('chr1', 50, 60, '.')]
cre = Loci.make("peaks.bed", keep=True)
cre
# -> Loci(n=7, chroms=2, sorted, cols=[name, score])
cre.columns, cre.shape
# -> (['chrom', 'start', 'end', 'strand', 'name', 'score'], (7, 6))
```

### Construction

| Method | Signature | One line |
|---|---|---|
| `Loci.make` | `make(filename, *, genome=None, sort=True, keep=False, backend=None)` | Read a BED / narrowPeak / broadPeak file (`.gz` too) through the tables backend; rows sorted into genome order unless `sort=False`. |
| `Loci.from_frame` | `from_frame(df, chrom=None, start=None, end=None, strand=None, *, keep=True, genome=None)` | From any table — pandas, polars (eager or lazy), pyarrow, bioframe, PyRanges, a dict of columns — columns found by name; row order kept. |
| `Loci.from_pandas`, `from_polars`, `from_arrow`, `from_bioframe` | same arguments as `from_frame` | Aliases of `from_frame`. |
| `Loci.from_pyranges` | `from_pyranges(gr, **kw)` | From a PyRanges (0.x or 1.x). |
| `Loci.from_bedtool` | `from_bedtool(bt, **kw)` | From a pybedtools `BedTool`. |
| `Loci.from_anndata` | `from_anndata(adata, axis="var", **kw)` | The regions of an AnnData axis, in its order: chrom/start/end columns if present, else parsed from names like `chr1:100-200`. |
| `Loci.from_records` | `from_records(items, *, genome=None)` | From `Locus` objects, `(chrom, start, end[, strand])` tuples, region strings or uids, in order. |
| `Loci.from_uids` | `from_uids(uids, *, genome=None)` | From `chrom:start-end(strand)` strings (what `uid` makes). |
| `Loci.tile` | `tile(chrom, size, chromsizes, *, genome=None)` | Uniform `size`-bp tiles over one chromosome (the last one clipped). |
| `Loci.tile_genome` | `tile_genome(chromsizes, size, chroms=None, *, genome=None)` | Tiles over every chromosome of `chromsizes` (or just `chroms`, in that order). |
| `Loci.load` | `load(path, *, genome=None)` | Read back a parquet file written by [`save`](#persistence). |

`make` always reads coordinates and strand (column 6 when present).
`keep=True` also keeps the other standard columns under their BED /
narrowPeak names (`name`, `score`, `signalValue`, `pValue`, `qValue`,
`peak`, ...); a list keeps just those. `.gtf` / `.gff` and `.bedpe` files
are refused with a pointer to `Genes.make` / `Pairs.make`.

```python
Loci.make("peaks.bed", keep=["name"]).columns
# -> ['chrom', 'start', 'end', 'strand', 'name']
Loci.make("peaks.bed", sort=False).is_sorted, Loci.make("peaks.bed").is_sorted
# -> (True, True)
Loci.make("genes.gtf")
# -> ValueError: this is a gene annotation: use Genes.make()
```

`from_frame` finds the coordinate columns case-insensitively among
`chrom / chromosome / chr / seqnames / Chromosome / #chrom`, `start /
chromStart / Start`, `end / chromEnd / End`, `strand / Strand`, or takes the
first three columns of a header-less text + int + int frame; otherwise it
raises naming the accepted spellings. `keep=True` keeps every other column.

```python
import pandas as pd
df = pd.DataFrame({"Chromosome": ["chr2", "chr1"], "Start": [5, 100], "End": [9, 200], "gene": ["a", "b"]})
L = Loci.from_frame(df)
L, L.to_records(), list(L.cols)
# -> (Loci(n=2, chroms=2, cols=[gene]), [('chr2', 5, 9, '.'), ('chr1', 100, 200, '.')], ['gene'])
Loci.from_frame(df, chrom="Chromosome", keep=False).columns
# -> ['chrom', 'start', 'end', 'strand']
Loci.from_frame(pd.DataFrame([["chr1", 5, 9, "x"]])).to_records()       # header-less: by position
# -> [('chr1', 5, 9, '.')]
```

```python
Loci.from_records([("chr1", 5, 15, "+"), "chr2:100-200", gb.Locus("chr1", 50, 60)]).uid
# -> array(['chr1:5-15(+)', 'chr2:100-200(.)', 'chr1:50-60(.)'], dtype=object)
Loci.from_uids(["chr1:5-15(+)", "chr2:100-200(.)"]).to_records()
# -> [('chr1', 5, 15, '+'), ('chr2', 100, 200, '.')]
Loci.tile("chr1", 5000, {"chr1": 20_000, "chr2": 8_000}).to_records()
# -> [('chr1', 0, 5000, '.'), ('chr1', 5000, 10000, '.'), ('chr1', 10000, 15000, '.'), ('chr1', 15000, 20000, '.')]
Loci.tile_genome("genome.chrom.sizes", 5000)
# -> Loci(n=6, chroms=2)
Loci.tile_genome("genome.chrom.sizes", 5000, chroms=["chr2"]).to_records()
# -> [('chr2', 0, 5000, '.'), ('chr2', 5000, 8000, '.')]
```

`chromsizes` is a length (for `tile`), a `{chrom: length}` dict, a
`.chrom.sizes` path, a pandas Series, a `Genome` with sizes, or a cooler.

### Copies and selections

| Method | Signature | One line |
|---|---|---|
| `take` | `take(idx) -> Loci` | Rows `idx` (int array, bool mask or slice) as a new Loci; slices and masks keep the sorted flag. |
| `copy` | `copy() -> Loci` | Deep copy of the arrays (same `Genome`). |
| `equals` | `equals(other, cols=False) -> bool` | Same rows (chromosome, start, end, strand) in the same order; `cols=True` also compares the extra columns (NaN equals NaN). |
| `head`, `tail` | `head(n=5)`, `tail(n=5) -> Loci` | The first / last rows. |

```python
cre.take([0, 1]).to_records()
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '-')]
cre.take(cre["score"] > 50).uid.tolist()
# -> ['chr2:500-600(+)', 'chr2:5000-5100(+)']
cre.take([0, 1]).equals(cre.head(2)), cre.copy().equals(cre, cols=True)
# -> (True, True)
```

---

## Like a table
{: .sec-green }

| Name | Returns | One line |
|---|---|---|
| `len(L)`, `bool(L)` | `int`, `bool` | Rows; an empty Loci is falsy. |
| `L.shape` | `(rows, columns)` | Like a DataFrame. |
| `L.columns` | `list` | `['chrom', 'start', 'end', 'strand']` + the extra columns. |
| `L.head(n)`, `L.tail(n)` | `Loci` | First / last `n` rows. |
| `L.describe()`, `L.summary()` | `DataFrame` | Rows, chromosomes, bases covered (merged), length quantiles, strand counts, sorted flag. |
| `L['start']`, `L['score']` | `ndarray` | A core column (`chrom` decoded, `strand` as symbols) or an extra column. |
| `L['score'] = values` | | Add or replace a column (one value per row, or a scalar broadcast). Coordinates cannot be set this way. |
| `L[i]` | `LocusView` | Row `i` (negative counts from the end) as a [`Locus`]({{ '/api/locus/' | relative_url }}) reading the columns. |
| `L['chr1:900-1100(+)']` | `LocusView` | The row with that uid. |
| `L[mask]`, `L[rows]`, `L[a:b]`, `L[[uid, ...]]` | `Loci` | A selection (`take`). |
| `uid in L`, `locus in L` | `bool` | Membership by uid. |
| `L.row(uid)` | `int` | Row number of a uid. |
| `for locus in L` | `LocusView` | Iterate the rows. |
| `L.to_numpy()`, `np.asarray(L)` | structured array | `chrom, start, end, strand` + extra columns. |
| `L._repr_html_()` | | Head / tail table in Jupyter and VS Code. |
| `__arrow_c_stream__`, `__dataframe__`, `__narwhals_dataframe__` | | `pl.DataFrame(L)`, `pa.table(L)`, `duckdb.sql("... from L")`, seaborn / plotly / altair take a Loci as is. |

```python
cre.head(2)
# -> Loci(n=2, chroms=1, sorted, cols=[name, score])
cre.describe()
# ->                             value
#    rows                             7
#    chromosomes                      2
#    bases covered (merged)        1100
#    length min                     100
#    length 25%                   100.0
#    length median                200.0
#    length 75%                   200.0
#    length max                     200
#    length mean             157.142857
#    strand + / - / .         4 / 2 / 1
#    genome-sorted                 True
cre["start"][:3], cre["name"][:3]
# -> (array([ 900, 1900, 4900]), array(['p1', 'p2', 'p3'], dtype=object))
cre["score2"] = cre["score"] * 2                  # a new column
cre["flag"] = True                                # a scalar is broadcast
cre.columns
# -> ['chrom', 'start', 'end', 'strand', 'name', 'score', 'score2', 'flag']
cre["start"] = np.zeros(7)
# -> KeyError: "'start' is a coordinate; build a new Loci instead"
```

```python
cre[0], cre["chr1:900-1100(+)"], cre.row("chr1:900-1100(+)")
# -> (Locus[0](chr1:900-1100(+)), Locus[0](chr1:900-1100(+)), 0)
cre[["chr2:500-600(+)", "chr1:900-1100(+)"]].uid.tolist()        # uids select rows, in that order
# -> ['chr2:500-600(+)', 'chr1:900-1100(+)']
cre[cre["score"] > 40].to_records()
# -> [('chr1', 10900, 11100, '-'), ('chr2', 500, 600, '+'), ('chr2', 5000, 5100, '+')]
len(cre[1:3]), len(cre[[0, 2]])
# -> (2, 2)
"chr1:900-1100(+)" in cre, gb.Locus("chr1", 900, 1100, "+") in cre
# -> (True, True)
for l in cre.head(2):
    print(l, l.name)
# -> Locus[0](chr1:900-1100(+)) p1
#    Locus[1](chr1:1900-2100(-)) p2
```

```python
np.asarray(cre)[:2]
# -> array([('chr1',  900, 1100, '+', 'p1', 10.), ('chr1', 1900, 2100, '-', 'p2', 20.)],
#          dtype=[('chrom', 'O'), ('start', '<i8'), ('end', '<i8'), ('strand', 'O'), ('name', 'O'), ('score', '<f8')])
import polars as pl, pyarrow as pa, duckdb
pl.DataFrame(cre).shape, pa.table(cre).num_rows
# -> ((7, 6), 7)
duckdb.sql("select chrom, count(*) n from cre group by chrom order by chrom").fetchall()
# -> [('chr1', 5), ('chr2', 2)]
```

{: .note }
A `LocusView` writes through: `cre[0].start = 850` changes `cre.starts[0]`
and drops the cached uids and indexes. It is meant for the odd fix-up — see
[`LocusView`]({{ '/api/locus/' | relative_url }}#locusviewlocus).

---

## Coordinates and derived columns
{: .sec-green }

The stored arrays are attributes; everything else is derived on access
(`uid` and `uids` are built once and cached until the rows change).

| Name | Type | One line |
|---|---|---|
| `codes`, `starts`, `ends`, `strands` | `int32`, `int64`, `int64`, `int8` arrays | The stored columns. |
| `cols` | `dict[str, ndarray]` | The extra columns. |
| `genome` | `Genome` | Chromosome names ↔ codes (+ sizes). |
| `filename` | `str` or `None` | Where the table was read from. |
| `chroms` | `object` array | Chromosome name per row (`genome.decode(codes)`). |
| `strand` | `object` array | `'.'`, `'+'`, `'-'` per row. |
| `centers` | `int64` array | `(start + end) // 2`. |
| `lengths` | `int64` array | `end - start`. |
| `uid` | `object` array | `chrom:start-end(strand)` per row. |
| `names` | `object` array | `chrom:start-end` per row — the region name AnnData / scATAC tools use. |
| `uids` | `dict` | `uid -> row`. |
| `is_sorted` | `bool` | Rows in genome order (natural chromosome order, then start, then end). |
| `chrom_offsets` | `dict` | `{chrom: (first_row, last_row + 1)}`; needs sorted rows. |
| `chrom_slice(chrom)` | `slice` | The row range of one chromosome (empty when absent). |
| `by_chrom(chrom)` | `Loci` | The rows on one chromosome. |

```python
cre.chroms[:3], cre.codes[:3], cre.starts[:3], cre.ends[:3], cre.strands[:3], cre.strand[:3]
# -> (array(['chr1', 'chr1', 'chr1'], dtype=object), array([0, 0, 0], dtype=int32),
#     array([ 900, 1900, 4900]), array([1100, 2100, 5100]), array([1, 2, 1], dtype=int8),
#     array(['+', '-', '+'], dtype=object))
cre.centers[:3], cre.lengths[:3]
# -> (array([1000, 2000, 5000]), array([200, 200, 200]))
cre.uid[:2].tolist(), cre.names[:2].tolist(), cre.uids["chr2:500-600(+)"]
# -> (['chr1:900-1100(+)', 'chr1:1900-2100(-)'], ['chr1:900-1100', 'chr1:1900-2100'], 5)
cre.genome, cre.is_sorted, cre.chrom_offsets, cre.chrom_slice("chr2")
# -> (Genome(2 chroms: chr1, chr2), True, {'chr1': (0, 5), 'chr2': (5, 7)}, slice(5, 7, None))
cre.by_chrom("chr2").uid.tolist()
# -> ['chr2:500-600(+)', 'chr2:5000-5100(+)']
```

### `Genome`

```python
Genome(names=(), sizes=None, name=None)
Genome.from_sizes(sizes, name=None)      # {chrom: length}, a .chrom.sizes path, a Series, a cooler
Genome.from_fasta(fasta, name=None)      # names and lengths from a FASTA (its .fai when present)
```

Chromosome names ↔ integer codes, plus sizes when known. Codes are
append-only (`encode` adds unseen names); `rank[code]` is the position in
natural sort order (chr1, chr2, ..., chr10, chrX), which is what `sort` uses,
so results never depend on the order names were first seen. Tables made from
one another share a `Genome`; tables read separately get their own, and an
operation between two of them re-codes the right-hand side onto the left
(see [Concepts]({{ '/concepts/' | relative_url }})). Pass `genome=g` to
several constructors to skip even that.

| Name | One line |
|---|---|
| `names`, `code` | `list` of names (index = code) and `{name: code}`. |
| `sizes`, `size(name)`, `set_sizes(sizes)` | `{chrom: length}` for the chromosomes whose size is known; add or update sizes. |
| `encode(values)`, `decode(codes)` | Names → int32 codes (adding unseen names) / codes → names. |
| `rank` | `rank[code]` = position in natural order. |
| `name in g`, `g[name]`, `len(g)` | Membership, code lookup, number of chromosomes. |

```python
g = Genome.from_sizes("genome.chrom.sizes")
g, g.sizes
# -> (Genome(2 chroms: chr1, chr2), {'chr1': 20000, 'chr2': 8000})
a = Loci.make("peaks.bed", genome=g); b = as_loci("chr1:1,000-2,000", genome=g)
a.genome is b.genome
# -> True
h = Genome(["chr2", "chr1"])
h.rank, h.encode(["chr1", "chrX"]), h.names
# -> (array([1, 0]), array([1, 2], dtype=int32), ['chr2', 'chr1', 'chrX'])
```

---

## Interval operations
{: .sec-green }

All of these take `backend=`. Intervals are half-open: `[s1, e1)` and
`[s2, e2)` overlap when `s1 < e2 and s2 < e1`, so book-ended intervals do not
overlap (but `merge` fuses them, like bedtools). Row numbers in the results
refer to `self` and to the other set *as passed* (after `as_loci`, which
keeps row order).

| Method | Signature | Returns |
|---|---|---|
| `overlap_rows` | `overlap_rows(key, start=None, end=None, *, backend=None)` | Row numbers overlapping one region: a region string, a `(chrom, start, end)` or `Locus`, or `chrom, start, end`. The lookup index is built once per backend and cached. |
| `overlaps` | `overlaps(key, start=None, end=None, *, backend=None)` | The same rows as a `Loci`. |
| `overlap_any` | `overlap_any(o, *, backend=None)` | Boolean mask: row overlaps anything in `o`. |
| `overlap_pairs` | `overlap_pairs(o, *, backend=None)` | `(rows of self, rows of o)` for every overlapping pair, sorted by (self row, o row). |
| `intersect` | `intersect(o, *, backend=None)` | Rows that overlap `o` (`A & B`, bedtools `intersect -u`). |
| `difference` | `difference(o, *, backend=None)` | Rows that do not overlap `o` (`A - B`, bedtools `intersect -v`). |
| `nearest` | `nearest(o, *, backend=None)` | `(row in o, distance)` per row: 0 for an overlap, else the gap in bases (book-ended = 0); `(-1, -1)` when `o` has nothing on that chromosome. |
| `merge` | `merge(*, backend=None)` | Overlapping or book-ended rows fused (strand of the first row of each block), sorted; extra columns dropped. |
| `slop` | `slop(n)` | Every interval widened by `n` bp on each side (clipped at 0); columns kept. |
| `sort` | `sort()` | Genome order: natural chromosome order, start, end; columns follow. |

```python
cre.overlap_rows("chr1:1,000-2,000"), cre.overlap_rows("chr1", 1000, 2000), cre.overlap_rows(gb.Locus("chr9", 0, 5))
# -> (array([0, 1]), array([0, 1]), array([], dtype=int64))
cre.overlaps("chr1:1,000-2,000").to_records()
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '-')]
other = as_loci([("chr1", 1000, 1050), ("chr2", 0, 10_000)])
cre.overlap_any(other)
# -> array([ True, False, False, False, False,  True,  True])
cre.overlap_pairs(other)
# -> (array([0, 5, 6]), array([0, 1, 1]))
cre.intersect(other).uid.tolist()
# -> ['chr1:900-1100(+)', 'chr2:500-600(+)', 'chr2:5000-5100(+)']
cre.slop(100).to_records()[:2]
# -> [('chr1', 800, 1200, '+'), ('chr1', 1800, 2200, '-')]
```

```python
q = as_loci([("chr1", 100, 110), ("chr1", 200, 210), ("chr3", 0, 10)])
r = as_loci([("chr1", 110, 120), ("chr1", 150, 160)])
q.nearest(r)                      # book-ended -> 0; 40 bp gap; nothing on chr3
# -> (array([ 0,  1, -1]), array([ 0, 40, -1]))
m = as_loci([("chr1", 10, 20, "+"), ("chr1", 20, 30, "-"), ("chr1", 40, 50, "-"), ("chr2", 0, 5)])
m.merge().to_records()
# -> [('chr1', 10, 30, '+'), ('chr1', 40, 50, '-'), ('chr2', 0, 5, '.')]
m.merge(backend="bioframe").to_records()          # any engine, the same blocks
# -> [('chr1', 10, 30, '+'), ('chr1', 40, 50, '-'), ('chr2', 0, 5, '.')]
as_loci([("chr10", 5, 6), ("chr2", 5, 6), ("chr1", 9, 10)]).sort().chroms
# -> array(['chr1', 'chr2', 'chr10'], dtype=object)
```

```python
with gb.use_backend(intervals="pyranges"):
    cre.intersect(other).to_records()
# -> [('chr1', 900, 1100, '+'), ('chr2', 500, 600, '+'), ('chr2', 5000, 5100, '+')]
cre.intersect(other, backend="cgranges")
# -> ImportError: the 'cgranges' intervals backend is not installed: conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)
```

{: .tip }
`ncls` and `cgranges` do overlaps and point lookups only; asking them for
`nearest` or `merge` raises `NotImplementedError` naming the engines that can
(`genomeblocks`, `bioframe`, `pyranges`, `bedtools`).

---

## Set operators
{: .sec-green }

| Operator | Method | Semantics |
|---|---|---|
| `a & b` | `a.intersect(b)` | Rows of `a` overlapping any row of `b`. |
| `a - b`, `a / b` | `a.difference(b)` | Rows of `a` overlapping no row of `b`. |
| `a + b`, `a \| b` | `a.union(b)` | Concatenation, `a`'s rows then `b`'s; only the extra columns both sides have are kept. |
| `a ^ b` | `(a - b) + (b - a)` | Rows of either set that do not overlap the other. |

The operators never clip intervals (no bedtools `intersect` without `-u`):
`a & b` returns whole rows of `a`. The right-hand side is anything
`as_loci` takes. Half-open rule: `[10, 20)` and `[20, 30)` share no base, so
`&` is empty, while `merge` fuses them.

```python
len(cre & other), len(cre - other), len(cre + other), len(cre ^ other), len(cre | other)
# -> (3, 4, 9, 4, 9)
(cre - other).uid.tolist()
# -> ['chr1:1900-2100(-)', 'chr1:4900-5100(+)', 'chr1:9950-10050(.)', 'chr1:10900-11100(-)']
(cre + other).columns                    # `other` has no name / score: nothing shared is kept
# -> ['chrom', 'start', 'end', 'strand']
a = as_loci([("chr1", 10, 20)]); b = as_loci([("chr1", 20, 30)])
len(a & b), (a + b).merge().to_records()
# -> (0, [('chr1', 10, 30, '.')])
```

---

## Converters
{: .sec-purple }

Every converter is one line in [`genomeblocks.interop`]({{ '/api/interop/' | relative_url }});
`Loci.from_*` are the reverse directions (above). Numeric columns go to
Arrow / polars without a copy; pandas, bioframe and pyranges frames are
copies.

| Method | Returns | Notes |
|---|---|---|
| `to_pandas(uid=False)` | `DataFrame` | `chrom` (categorical in genome order), `start`, `end`, `strand` (categorical), the extra columns, `uid` when asked — bioframe's column names. |
| `to_polars()` | `polars.DataFrame` | Through Arrow. |
| `to_arrow()` | `pyarrow.Table` | `chrom` / `strand` as dictionary arrays; object columns with NaN / None become Arrow nulls. |
| `to_bioframe()` | `DataFrame` | Plain `chrom` / `start` / `end` / `strand` string columns, ready for `bioframe.*`. |
| `to_pyranges()` | `PyRanges` | `Chromosome / Start / End / Strand` + extra columns (pyranges 0.x groups rows by chromosome: sort first if row order matters). |
| `to_bedtool()` | `BedTool` | BED6: `name` = `cols['name']` or the uid, `score` = `cols['score']` or 0. |
| `to_cgranges()` | `cgranges` index | Built, label = row (needs `cgranges`). |
| `to_anndata(X=None, *, obs=None, layers=None, **kw)` | `AnnData` | Loci as `var` (index `chrom:start-end`), `X` is `(n_obs, n_loci)`, `obs` a DataFrame or a list of names. |
| `to_records()` | `list` | `(chrom, start, end, strand)` tuples. |
| `to_numpy()` | structured array | Also `np.asarray(L)`. |

```python
cre = Loci.make("peaks.bed", keep=True)
cre.to_pandas().head(3)
# ->   chrom  start   end strand name  score
#    0  chr1    900  1100      +   p1   10.0
#    1  chr1   1900  2100      -   p2   20.0
#    2  chr1   4900  5100      +   p3   30.0
cre.to_pandas(uid=True).columns.tolist()
# -> ['chrom', 'start', 'end', 'strand', 'name', 'score', 'uid']
cre.to_polars().schema
# -> Schema([('chrom', Categorical), ('start', Int64), ('end', Int64), ('strand', Categorical), ('name', String), ('score', Float64)])
cre.to_arrow().schema.names
# -> ['chrom', 'start', 'end', 'strand', 'name', 'score']
cre.to_bioframe().dtypes.to_dict()
# -> {'chrom': dtype('O'), 'start': dtype('int64'), 'end': dtype('int64'), 'strand': dtype('O'), 'name': dtype('O'), 'score': dtype('float64')}
cre.to_pyranges()                                 # a PyRanges with 7 rows, Chromosome / Start / End / Strand / name / score
type(cre.to_bedtool()).__name__
# -> 'BedTool'
cre.to_records()[:2]
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '-')]
X = np.random.default_rng(0).random((3, len(cre)))
ad = cre.to_anndata(X, obs=["s1", "s2", "s3"])
ad
# -> AnnData object with n_obs × n_vars = 3 × 7
#        var: 'chrom', 'start', 'end', 'strand', 'name', 'score'
Loci.from_anndata(ad).equals(cre)
# -> True
```

---

## Persistence
{: .sec-purple }

| Method | Signature | One line |
|---|---|---|
| `save` | `save(path)` | Parquet of `to_arrow()` (every column, dictionary-encoded chrom / strand). |
| `Loci.load` | `load(path, *, genome=None)` | Reads it back onto any `Genome`; sets `filename`. |
| `to_bed` | `to_bed(path=None, *, name=None, score=None)` | BED6 text (returned, or written to `path`): name = `cols[name]` or the uid, score = `cols[score]` or 0. |
| `to_fasta` | `to_fasta(path, fasta, r=None, *, strand=False, **kw)` | The rows' sequences as FASTA, headers = uids (see [Sequences](#sequences-and-liftover)). |

```python
cre.save("cre.parquet")
Loci.load("cre.parquet").equals(cre, cols=True), Loci.load("cre.parquet").filename
# -> (True, 'cre.parquet')
as_loci("cre.parquet")                            # as_loci reads parquet too
# -> Loci(n=7, chroms=2, cols=[name, score])
cre.to_bed().splitlines()[:2]
# -> ['chr1\t900\t1100\tchr1:900-1100(+)\t0\t+', 'chr1\t1900\t2100\tchr1:1900-2100(-)\t0\t-']
cre.to_bed(name="name", score="score").splitlines()[:2]
# -> ['chr1\t900\t1100\tp1\t10.0\t+', 'chr1\t1900\t2100\tp2\t20.0\t-']
cre.to_bed("cre.bed", name="name")
Loci.make("cre.bed", keep=True).columns
# -> ['chrom', 'start', 'end', 'strand', 'name', 'score']
```

Pickling works and drops the cached uids and lookup indexes, so a pickle
never carries a stale tree.

```python
import pickle
pickle.loads(pickle.dumps(cre)).equals(cre)
# -> True
```

---

## Sequences and liftover
{: .sec-purple }

| Method | Signature | One line |
|---|---|---|
| `sequences` | `sequences(fasta, r=None, *, strand=False, upper=False, backend=None) -> list[str]` | Sequence of every row (or of `center ± r`) from a FASTA path, a `{chrom: str}` dict or an open pysam / pyfaidx / Biopython handle, through the [fasta backend]({{ '/api/backends/' | relative_url }}). Windows are clipped to the chromosome; `strand=True` reverse-complements `-` rows. |
| `to_seqrecords` | `to_seqrecords(fasta, r=None, *, strand=False, **kw)` | Biopython `SeqRecord` per row (id = uid). |
| `to_fasta` | `to_fasta(path, fasta, r=None, *, strand=False, **kw) -> path` | Write them as FASTA. |
| `liftover` | `liftover(chain_file, *, min_match=0.95, verbose=True) -> Loci` | Lift to another assembly with a UCSC chain file (needs `pyliftover`). Rows that do not lift are dropped; `cols['source_row']` says where each row came from. |

```python
cre.sequences("genome.fa", r=5)[:3]
# -> ['TCAGCACGTA', 'TGATGTTGTG', 'CGTTAGGAAG']
cre.sequences("genome.fa", r=5, strand=True)[1]           # row 1 is on '-': reverse complement
# -> 'CACAACATCA'
cre.to_seqrecords("genome.fa")[0].id
# -> 'chr1:900-1100(+)'
cre.to_fasta("cre.fa", "genome.fa", r=10)
# -> 'cre.fa'
```

`liftover` walks each end inward by up to `1 - min_match` of the length
when it falls in a chain gap (UCSC's base-fraction rule) and keeps a row only
when both ends land on one chromosome and strand.

```python
# a chain that shifts chr1 onto chrA by +10,000 and says nothing about chr2
out = cre.liftover("t.chain")
# -> [INFO] liftover: 7 → 5 (2 dropped, min_match=0.95)
out, out.to_records()[:2], out["source_row"]
# -> (Loci(n=5, chroms=1, cols=[source_row]), [('chrA', 10900, 11100, '+'), ('chrA', 11900, 12100, '-')], array([0, 1, 2, 3, 4]))
cre.liftover("t.chain")                                  # without pyliftover
# -> ImportError: Loci.liftover needs pyliftover: pip install pyliftover
```

---

## Analysis delegations
{: .sec-purple }

These methods call the module function of the same name with `self` as the
first argument, so the full signatures live on the module pages. Results line
up with the rows.

| Method | Calls | Returns |
|---|---|---|
| `signal(bigwigs, **kw)` | [`signal.signal`]({{ '/api/signal/' | relative_url }}) | `float32` cube `(rows, tracks, bins)`; `backend=` picks the bigwig engine. |
| `plot_heatmap(S, **kw)`, `plot_profiles(S, **kw)` | [`signal_draw`]({{ '/api/signal/' | relative_url }}) | matplotlib `Figure`. |
| `scan_motifs(fasta, motifs, **kw)` | [`motifs.scan_motifs`]({{ '/api/motifs/' | relative_url }}) | `{motif: total hits}`; `backend=` picks the motifs engine, `fasta_backend=` (matrix) the FASTA reader. |
| `scan_motifs_matrix(fasta, motifs, **kw)` | `motifs.scan_motifs_matrix` | `DataFrame` (rows x motifs), index = uid. |
| `scan_motifs_matrix_masked(fasta, motifs, anchors, **kw)` | `motifs.scan_motifs_matrix_masked` | The same after masking the anchor motifs' matches. |
| `scan_motifs_profile(fasta, motifs, select=None, **kw)` | `motifs.scan_motifs_profile` | `(cube (rows, motifs, bins), names)`. |
| `enrich(atlas, **kw)` | [`Atlas.search`]({{ '/api/atlas/' | relative_url }}) | Per-track enrichment `DataFrame`. |
| `enrich_mc(atlas, *, n=10, **kw)` | `Atlas.bootstrap` | Bootstrapped enrichment. |
| `count_pairs(pairs_file, **kw)` | [`bedpe.count_pairs`]({{ '/api/bedpe/' | relative_url }}) | Contacts per row by partner chromosome. |
| `count_pairs_2d(pairs_file, **kw)` | `bedpe.count_pairs_2d` | Sparse `(rows, rows)` contact matrix. |
| `call_se(bigwigs, **kw)` | [`se.call_se`]({{ '/api/atlas/' | relative_url }}) | ROSE-style super-enhancers as a `Loci` with `n_peaks`, `score`, `rank`; `backend=` picks the intervals engine for stitching. |
| `liftover`, `sequences`, `to_seqrecords`, `to_fasta` | `interop` | See above. |

```python
S = cre.signal(["signal.bw", "signal2.bw"], n_bins=20, flank=500, verbose=False, progress=False)
S.shape, S.dtype
# -> ((7, 2, 20), dtype('float32'))
cre.scan_motifs("genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
# -> {'M1': 4.0, 'M2': 4.0}
cre.scan_motifs_matrix("genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
# ->                      M1  M2
#    uid
#    chr1:900-1100(+)      3   1
#    chr1:1900-2100(-)     0   0
#    chr1:4900-5100(+)     0   3
#    chr1:9950-10050(.)    0   0
#    chr1:10900-11100(-)   0   0
#    chr2:500-600(+)       1   0
#    chr2:5000-5100(+)     0   0
prof, names = cre.scan_motifs_profile("genome.fa", "motifs.jaspar", r=100, n_bins=10, threshold=7.0, verbose=False)
prof.shape, prof.dtype, names
# -> ((7, 2, 10), dtype('float32'), ['M1', 'M2'])
cre.count_pairs("contacts.allValidPairs", verbose=False).head(3)
# ->   chrom  start   end                uid  chr1  chr2
#    0  chr1    900  1100   chr1:900-1100(+)     3     1
#    1  chr1   1900  2100  chr1:1900-2100(-)     7     7
#    2  chr1   4900  5100  chr1:4900-5100(+)     4     2
se = cre.call_se("signal.bw", stitch=2000, verbose=False)
se, se.columns
# -> (Loci(n=3, chroms=1, cols=[n_peaks, score, rank]), ['chrom', 'start', 'end', 'strand', 'n_peaks', 'score', 'rank'])
type(cre.plot_heatmap(S, samples=["ATAC", "H3K27ac"])).__name__
# -> 'Figure'
```

{: .warning }
The motif functions default to `threshold=13.0` log2-odds, a sensible cutoff
for real 10–20 bp motifs; the 4 bp toy motifs here never reach it (their best
score is 7.83), so the examples pass `threshold=7.0`. Use `pvalue=` for a
per-motif cutoff.
