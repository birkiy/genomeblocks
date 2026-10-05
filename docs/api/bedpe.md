---
title: bedpe
parent: API Reference
layout: default
nav_order: 7
---

# `genomeblocks.bedpe`
{: .no_toc }

BEDPE as a table — [`Pairs`](#pairs) holds the two anchors of every loop as
two row-aligned [`Loci`]({{ '/api/loci/' | relative_url }}) (`P.a`, `P.b`)
plus extra columns — and contact counting from pairs files
(`count_pairs`, `count_pairs_2d`), which stream HiC-Pro allValidPairs, 4DN
`.pairs` or Juicer files once. This module is the one BEDPE reader in
genomeblocks: `Architecture.make` and the browser call it under the hood.
See the [BEDPE guide]({{ '/guide/bedpe/' | relative_url }}) for the workflow
and [Design: bedpe]({{ '/design/bedpe/' | relative_url }}) for the counting
pass.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
from genomeblocks import Pairs
from genomeblocks.bedpe import (read_bedpe, as_pairs, read_pairs_chunks, PAIRS_FORMAT_COLUMNS,
                                count_pairs, count_pairs_2d, pair_2d_block, pair_2d_to_frame)
```

The examples use a four-loop BEDPE (`loops.bedpe`: three cis loops, one
chr1–chr2 trans loop, with `name` and `score` columns) and a 400-read
HiC-Pro allValidPairs file on chr1 / chr2.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Pairs`
{: .sec-navy }

```python
Pairs(a: Loci, b: Loci, cols=None, *, filename=None)
```

Anchor `a` and anchor `b` need one row per pair; `b` is re-coded onto `a`'s
`Genome` when they differ. `cols` is `{name: array}`. Like every table a
`Pairs` has `len`, `shape`, `columns`, `head()`, `tail()`, `describe()`,
`P['score']`, `to_pandas()` / `to_polars()` / `to_arrow()` and the Arrow /
dataframe-interchange / narwhals protocols.

| Column | Where it lives |
|---|---|
| `chrom1`, `start1`, `end1`, `strand1` | `P.a` (a `Loci`: `codes`, `starts`, `ends`, `strands`). |
| `chrom2`, `start2`, `end2`, `strand2` | `P.b`. |
| `name`, `score`, `col11`, ... | `P.cols`; `P.columns` lists them in BEDPE order (name, score, strands, then the rest). |

```python
P = Pairs.make("loops.bedpe")
P
# -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
P.columns, P.shape
# -> (['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2'], (4, 10))
P.a, P.b
# -> (Loci(n=4, chroms=2), Loci(n=4, chroms=2))
P.a.to_records()
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '+'), ('chr1', 900, 1100, '-'), ('chr2', 500, 600, '+')]
P.cols
# -> {'name': array(['l1', 'l2', 'l3', 'l4'], dtype=object), 'score': array([5., 3., 1., 2.])}
```

### Construction

| Method | Signature | One line |
|---|---|---|
| `Pairs.make` | `make(filename, *, genome=None, min_score=None, max_distance=None, backend=None)` | Read a BEDPE (`.gz` too): the six anchor columns, then `name`, `score`, `strand1`, `strand2` when present; further columns are kept as `col11`, `col12`, ... `backend` picks the tables parser. |
| `Pairs.from_frame` | `from_frame(df, *, genome=None)` | From a pandas / polars / pyarrow frame with BEDPE columns (`chrom1 ... end2` by name, case-insensitive, else the first six); every other column is kept. |
| `Pairs.load` | `load(path, *, genome=None)` | Read back a parquet written by `save`. |
| `read_bedpe` | `read_bedpe(filename, **kw) -> Pairs` | Module-level alias of `Pairs.make`. |
| `as_pairs` | `as_pairs(x, *, genome=None) -> Pairs` | A `Pairs` as is; a `.bedpe` path (`make`) or a `.parquet` path (`load`); any frame (`from_frame`). What every consumer calls on its input. |

`min_score` / `max_distance` filter while reading (see [`filter`](#selections)).
A file with fewer than six columns raises; a `#` header line is skipped.

```python
Pairs.make("loops.bedpe", min_score=2, max_distance=5000)
# -> Pairs(n=2, cis=2, trans=0, cols=[name, score])
Pairs.make("loops.bedpe", backend="pandas").a.equals(P.a)
# -> True
Pairs.from_frame(P.to_pandas()).a.equals(P.a)
# -> True
import pandas as pd
Pairs.from_frame(pd.DataFrame({"c1": ["chr1"], "s1": [1], "e1": [5], "c2": ["chr1"], "s2": [100], "e2": [105]}))
# -> Pairs(n=1, cis=1, trans=0)
read_bedpe("loops.bedpe"), as_pairs(P) is P, as_pairs("loops.parquet"), as_pairs(P.to_polars())
# -> (Pairs(n=4, cis=3, trans=1, cols=[name, score]), True, Pairs(n=4, cis=3, trans=1, cols=[name, score]), Pairs(n=4, cis=3, trans=1, cols=[name, score]))
Pairs.make("x.bedpe")                           # a three-column file
# -> ValueError: x.bedpe: a BEDPE needs at least 6 columns, found 3
```

### Rows, columns and derived arrays

| Name | Returns | One line |
|---|---|---|
| `len(P)`, `P.shape`, `P.columns` | | Like a DataFrame. |
| `P[i]` | `(LocusView, LocusView)` | The two anchors of pair `i`. |
| `P['score']` | `ndarray` | An extra column. |
| `P[mask]`, `P[rows]`, `P[a:b]` | `Pairs` | A selection (`take`). |
| `for a, b in P` | | Iterate the anchor pairs. |
| `take(idx)` | `Pairs` | Rows `idx` (int array, mask or slice); keeps `filename`. |
| `head(n=5)`, `tail(n=5)` | `Pairs` | First / last rows. |
| `genome` | `Genome` | `P.a.genome`. |
| `mids1`, `mids2` | `int64` arrays | Anchor midpoints (`a.centers`, `b.centers`). |
| `is_cis` | `bool` array | `a.codes == b.codes`. |
| `distance` | `float` array | `abs(mid2 - mid1)` for cis pairs, `inf` for trans. |
| `describe()`, `summary()` | `DataFrame` | Pairs, cis / trans, chromosomes, cis-distance quantiles, score range. |

```python
P[0]
# -> (Locus[0](chr1:900-1100(+)), Locus[0](chr1:4900-5100(-)))
P["score"], len(P[[0, 2]])
# -> (array([5., 3., 1., 2.]), 2)
for a, b in P.head(2):
    print(a.uid, b.uid)
# -> chr1:900-1100(+) chr1:4900-5100(-)
#    chr1:1900-2100(+) chr1:10900-11100(+)
P.genome, P.is_cis, P.distance
# -> (Genome(2 chroms: chr1, chr2), array([ True,  True, False,  True]), array([4000., 9000.,   inf, 4500.]))
P.mids1, P.mids2
# -> (array([1000, 2000, 1000,  550]), array([ 5000, 11000,   550,  5050]))
P.describe()
# ->                       value
#    pairs                     4
#    cis                       3
#    trans                     1
#    chromosomes               2
#    cis distance min       4000
#    cis distance median  4500.0
#    cis distance max       9000
#    score min               1.0
#    score median            2.5
#    score max               5.0
```

### Selections

| Method | Signature | One line |
|---|---|---|
| `filter` | `filter(*, min_score=None, max_distance=None) -> Pairs` | Keep pairs with `score >= min_score` (when a `score` column exists) and `distance <= max_distance` (trans pairs are `inf`, so a `max_distance` drops them). |
| `anchors_overlap` | `anchors_overlap(loci, *, r=0, backend=None) -> (mask_a, mask_b)` | Whether anchor `a` / anchor `b` of each pair overlaps `loci`, anchors widened by `r` bp first. |
| `overlapping` | `overlapping(loci, *, r=0, both=False, backend=None) -> Pairs` | Pairs with one anchor (or `both`) overlapping `loci` — bedtools `pairtobed`. |

`loci` is anything [`as_loci`]({{ '/api/interop/' | relative_url }})
takes; `backend` picks the intervals engine.

```python
P.filter(min_score=3), P.filter(max_distance=5000)
# -> (Pairs(n=2, cis=2, trans=0, cols=[name, score]), Pairs(n=2, cis=2, trans=0, cols=[name, score]))
P.anchors_overlap([("chr1", 900, 1100)])
# -> (array([ True, False,  True, False]), array([False, False, False, False]))
cre = gb.Loci.make("peaks.bed")
P.overlapping(cre), P.overlapping(cre, both=True)
# -> (Pairs(n=4, cis=3, trans=1, cols=[name, score]), Pairs(n=4, cis=3, trans=1, cols=[name, score]))
P.overlapping("chr2:0-1000"), P.overlapping("chr2:0-1000", both=True)
# -> (Pairs(n=2, cis=1, trans=1, cols=[name, score]), Pairs(n=0, cis=0, trans=0, cols=[name, score]))
```

### Export and persistence

| Method | Returns | One line |
|---|---|---|
| `to_pandas()` | `DataFrame` | `chrom1 start1 end1 chrom2 start2 end2 [name score] strand1 strand2 ...` — chromosomes and strands as plain strings. |
| `to_polars()`, `to_arrow()` | | Through pandas. |
| `to_numpy()` | structured array | The same columns. |
| `to_bedpe(path)` | `path` | Tab-separated BEDPE, no header; `name` / `score` filled with `.` / 0 when absent, so the strands land in columns 9–10. |
| `save(path)`, `Pairs.load(path)` | | Parquet round trip. |

```python
P.to_pandas()
# ->   chrom1  start1  end1 chrom2  start2   end2 name  score strand1 strand2
#    0   chr1     900  1100   chr1    4900   5100   l1    5.0       +       -
#    1   chr1    1900  2100   chr1   10900  11100   l2    3.0       +       +
#    2   chr1     900  1100   chr2     500    600   l3    1.0       -       +
#    3   chr2     500   600   chr2    5000   5100   l4    2.0       +       +
P.to_polars().schema
# -> Schema([('chrom1', String), ('start1', Int64), ('end1', Int64), ('chrom2', String), ('start2', Int64), ('end2', Int64), ('name', String), ('score', Float64), ('strand1', String), ('strand2', String)])
P.to_numpy()[:1]
# -> rec.array([('chr1', 900, 1100, 'chr1', 4900, 5100, 'l1', 5., '+', '-')], ...)
P.to_bedpe("out.bedpe"); open("out.bedpe").readline()
# -> 'chr1\t900\t1100\tchr1\t4900\t5100\tl1\t5.0\t+\t-\n'
P.save("loops.parquet"); Pairs.load("loops.parquet")
# -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
import polars as pl
pl.DataFrame(P).shape                           # the Arrow protocol
# -> (4, 10)
```

{: .note }
`Architecture.make(cre, P, r=100)` takes a `Pairs` directly (or anything
`as_pairs` takes) and maps every anchor to the CREs within `r` bp of its
midpoint — see [`Architecture`]({{ '/api/architecture/' | relative_url }}).

---

## Pairs files: contacts per window
{: .sec-green }

`count_pairs` and `count_pairs_2d` read a pairs file in chunks, locate each
end in the windows of a `Loci` by a sorted search (so the windows on one
chromosome must not overlap), and count. One pass, constant memory in the
file size.

### `read_pairs_chunks` and `PAIRS_FORMAT_COLUMNS`

```python
read_pairs_chunks(filename, *, format="auto", columns=None, chunksize=2_000_000) -> Iterator[DataFrame]
```

Streams the file as DataFrames of `chrom1, pos1, chrom2, pos2` (categorical
chromosomes, int64 positions). `format` is `'auto'` (sniffed from the first
data line), `'allvalidpairs'`, `'pairs'` or `'juicer'`; `columns=(c1, p1,
c2, p2)` gives the 0-based columns of any other layout. `#` lines are
skipped (the 4DN `.pairs` header).

| Format | Columns `(chrom1, pos1, chrom2, pos2)` | Layout |
|---|---|---|
| `allvalidpairs` | `(1, 2, 4, 5)` | HiC-Pro: `readID chr1 pos1 strand1 chr2 pos2 strand2 ...` |
| `pairs` | `(1, 2, 3, 4)` | pairtools / 4DN: `readID chrom1 pos1 chrom2 pos2 strand1 strand2 ...` |
| `juicer` | `(2, 3, 6, 7)` | Juicer medium: `readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 ...` |

```python
PAIRS_FORMAT_COLUMNS
# -> {'allvalidpairs': (1, 2, 4, 5), 'pairs': (1, 2, 3, 4), 'juicer': (2, 3, 6, 7)}
chunk = next(read_pairs_chunks("contacts.allValidPairs"))
chunk.head(3)
# ->   chrom1  pos1 chrom2  pos2
#    0   chr2  9658   chr1  9193
#    1   chr2  3196   chr2  2021
#    2   chr2  7985   chr2  3429
next(read_pairs_chunks("contacts.allValidPairs", chunksize=10)).shape
# -> (10, 4)
next(read_pairs_chunks("contacts.allValidPairs", columns=(1, 2, 4, 5))).shape
# -> (400, 4)
next(read_pairs_chunks("contacts.allValidPairs", format="nope"))
# -> ValueError: Unknown pairs format: 'nope'. Known: ['allvalidpairs', 'juicer', 'pairs'], or pass columns=(c1, p1, c2, p2).
```

### `count_pairs`

```python
count_pairs(loci, pairs_file, *, target_chrom=None, format="auto", columns=None,
            chunksize=2_000_000, verbose=True) -> DataFrame
```

Contacts landing in each window of `loci`, split by partner chromosome. For
each pair, each end inside a window counts once for that window under the
other end's chromosome (a cis pair with both ends in windows counts once per
end). The result is aligned to the rows of `loci`: `chrom`, `start`, `end`,
`uid`, then one column per partner chromosome with at least one contact, or
a single `count` column with `target_chrom`. Also `Loci.count_pairs(pairs_file, **kw)`.

```python
W = gb.as_loci([("chr1", s, s + 1000) for s in range(0, 10_000, 2000)] + [("chr2", 0, 5000)])
count_pairs(W, "contacts.allValidPairs")
# -> [INFO] processed 400 pairs total
#      chrom  start   end                uid  chr1  chr2
#    0  chr1      0  1000     chr1:0-1000(.)    19     5
#    1  chr1   2000  3000  chr1:2000-3000(.)    28    16
#    2  chr1   4000  5000  chr1:4000-5000(.)    24    11
#    3  chr1   6000  7000  chr1:6000-7000(.)    22    14
#    4  chr1   8000  9000  chr1:8000-9000(.)    19    12
#    5  chr2      0  5000     chr2:0-5000(.)    55   153
count_pairs(W, "contacts.allValidPairs", target_chrom="chr2", verbose=False)["count"].tolist()
# -> [5, 16, 11, 14, 12, 153]
W.count_pairs("contacts.allValidPairs", verbose=False).columns.tolist()
# -> ['chrom', 'start', 'end', 'uid', 'chr1', 'chr2']
```

### `count_pairs_2d`

```python
count_pairs_2d(loci_a, pairs_file, *, loci_b=None, format="auto", columns=None,
               chunksize=2_000_000, verbose=True) -> scipy.sparse.csr_matrix
```

Window-to-window contact counts, rows of `loci_a` x rows of `loci_b`
(`loci_b=None` uses `loci_a`, which makes the matrix symmetric). Both
orientations of every pair are counted, so a pair with both ends inside one
window adds 2 to its diagonal cell in the symmetric case. `int64` CSR. Also
`Loci.count_pairs_2d(pairs_file, **kw)`.

```python
M = count_pairs_2d(W, "contacts.allValidPairs", verbose=False)
type(M).__name__, M.shape, M.dtype
# -> ('csr_matrix', (6, 6), dtype('int64'))
M.toarray()
# -> array([[ 2,  1,  2,  1,  3,  1],
#           [ 1,  4,  1,  3,  0,  9],
#           [ 2,  1,  2,  4,  0,  5],
#           [ 1,  3,  4,  2,  0,  7],
#           [ 3,  0,  0,  0,  0,  3],
#           [ 1,  9,  5,  7,  3, 74]])
count_pairs_2d(W, "contacts.allValidPairs", loci_b=W.head(2), verbose=False).shape
# -> (6, 2)
```

### `pair_2d_block` and `pair_2d_to_frame`

```python
pair_2d_block(mat, loci_a, loci_b, chrom_a, chrom_b) -> (ndarray, Loci, Loci)
pair_2d_to_frame(mat, loci_a, loci_b) -> DataFrame
```

`pair_2d_block` cuts the dense block of one chromosome pair out of a
`count_pairs_2d` matrix and returns it with the windows of `loci_a` on
`chrom_a` and of `loci_b` on `chrom_b` (an empty block when a chromosome is
absent). `pair_2d_to_frame` lists the non-zero cells as a long table
`chrom1 start1 end1 chrom2 start2 end2 count`.

```python
blk, ra, cb = pair_2d_block(M, W, W, "chr1", "chr1")
blk, ra, cb
# -> (array([[2, 1, 2, 1, 3],
#            [1, 4, 1, 3, 0],
#            [2, 1, 2, 4, 0],
#            [1, 3, 4, 2, 0],
#            [3, 0, 0, 0, 0]]), Loci(n=5, chroms=1), Loci(n=5, chroms=1))
pair_2d_block(M, W, W, "chr1", "chr3")[0].shape
# -> (5, 0)
pair_2d_to_frame(M, W, W).head()
# ->   chrom1  start1  end1 chrom2  start2  end2  count
#    0   chr1       0  1000   chr1       0  1000      2
#    1   chr1       0  1000   chr1    2000  3000      1
#    2   chr1       0  1000   chr1    4000  5000      2
#    3   chr1       0  1000   chr1    6000  7000      1
#    4   chr1       0  1000   chr1    8000  9000      3
```

{: .tip }
For a heatmap of one chromosome, `pair_2d_block(M, W, W, "chr1", "chr1")[0]`
is the array to draw; `pair_2d_to_frame` is the shape BEDPE-like tools and
`Pairs.from_frame` (with a `count` column kept) expect.
