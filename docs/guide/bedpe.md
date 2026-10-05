---
title: BEDPE
parent: User Guide
layout: default
nav_order: 7
---

# BEDPE
{: .no_toc }

Paired intervals as one table: read loop files into `Pairs`, filter them,
intersect their anchors with `Loci` (bedtools `pairtobed`), and count raw
Hi-C / HiChIP contacts from `.pairs`, `allValidPairs` and Juicer files into
windows.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is
{: .sec-purple }

`gb.Pairs` holds a BEDPE as two row-aligned [`Loci`]({{ '/guide/loci/' | relative_url }})
plus extra columns. Pair *i* is anchor `P.a[i]` and anchor `P.b[i]`; both
tables share one `Genome`.

| Attribute | Content |
|---|---|
| `P.a` | a `Loci`: `chrom1`, `start1`, `end1`, `strand1` |
| `P.b` | a `Loci`: `chrom2`, `start2`, `end2`, `strand2` |
| `P.cols` | `name`, `score` when the file has them, then any further columns (`col11`, ...) |

```python
import genomeblocks as gb

P = gb.Pairs.make("loops.bedpe")
P
# -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
P.shape, P.columns
# -> ((4, 10), ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2'])
P.to_pandas()
# ->   chrom1  start1  end1 chrom2  start2   end2 name  score strand1 strand2
# -> 0   chr1     900  1100   chr1    4900   5100   l1    5.0       +       -
# -> 1   chr1    1900  2100   chr1   10900  11100   l2    3.0       +       +
# -> 2   chr1     900  1100   chr2     500    600   l3    1.0       -       +
# -> 3   chr2     500   600   chr2    5000   5100   l4    2.0       +       +
```

Like every table it has `len`, `shape`, `columns`, `head()`, `tail()`,
`describe()`, `P['score']`, an HTML repr, and the Arrow / dataframe protocols
(`pl.DataFrame(P)`, `duckdb.sql("select * from P")`). This module is the one
BEDPE reader in `genomeblocks`: `Architecture.make` and the browser call it
for you.

---

## Reading and building
{: .sec-purple }

### `Pairs.make`

```python
P = gb.Pairs.make("loops.bedpe", min_score=3.0, max_distance=2e6)
P
# -> Pairs(n=2, cis=2, trans=0, cols=[name, score])
```

| arg | default | meaning |
|---|---|---|
| `min_score` | `None` | drop pairs with `score` below this (needs a score column) |
| `max_distance` | `None` | drop cis pairs whose anchor midpoints are farther apart — and every trans pair, whose distance is `inf` |
| `genome` | `None` | a `Genome` to share with other tables |
| `backend` | `None` | the table parser: `'polars'` or `'pandas'` |

The six anchor columns are required (fewer raises); `name`, `score`,
`strand1`, `strand2` are read when present; columns beyond 10 are kept as
`col11`, `col12`, ... `.gz` files and leading `#` header lines are fine.
`genomeblocks.bedpe.read_bedpe(path, **kw)` is the same call.

### `Pairs.from_frame` and `as_pairs`

Any pandas / polars / pyarrow frame with `chrom1 start1 end1 chrom2 start2 end2`
columns (by name, case-insensitive, else the first six); `strand1` / `strand2`
are read by name and every other column is kept.

```python
import pandas as pd
df = pd.DataFrame({"chrom1": ["chr1"], "start1": [100], "end1": [200],
                   "chrom2": ["chr1"], "start2": [900], "end2": [1000], "fdr": [0.01]})
F = gb.Pairs.from_frame(df)
F, F.columns
# -> (Pairs(n=1, cis=1, trans=0, cols=[fdr]), ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'strand1', 'strand2', 'fdr'])
```

`genomeblocks.bedpe.as_pairs(x)` takes a `Pairs` (returned as is), a BEDPE or
`.parquet` path, or a frame — it is what `Architecture.make` and the browser
call on their `bedpe` argument, so all of those forms work there too.

A `.bedpe` handed to `Loci.make` raises *"use Pairs.make()"*.

---

## Rows, columns and anchors
{: .sec-purple }

```python
P = gb.Pairs.make("loops.bedpe")
P[0]                                   # the two anchors of pair 0, as Locus views
# -> (Locus[0](chr1:900-1100(+)), Locus[0](chr1:4900-5100(-)))
P["score"], P["name"]                  # a column
# -> (array([5., 3., 1., 2.]), array(['l1', 'l2', 'l3', 'l4'], dtype=object))
P[P["score"] >= 3]                     # a mask, row numbers or a slice: a new Pairs
# -> Pairs(n=2, cis=2, trans=0, cols=[name, score])
for a, b in P.head(2):                 # iteration yields (anchor a, anchor b)
    print(a.uid, "->", b.uid)
# -> chr1:900-1100(+) -> chr1:4900-5100(-)
# -> chr1:1900-2100(+) -> chr1:10900-11100(+)
```

Per-pair arrays:

| Property | Meaning |
|---|---|
| `P.is_cis` | both anchors on one chromosome |
| `P.distance` | `\|mid2 - mid1\|` for cis pairs, `inf` for trans |
| `P.mids1`, `P.mids2` | anchor midpoints (`P.a.centers`, `P.b.centers`) |

```python
P.is_cis, P.distance
# -> (array([ True,  True, False,  True]), array([4000., 9000.,   inf, 4500.]))
P.describe()
# ->                      value
# -> pairs                   4.0
# -> cis                     3.0
# -> trans                   1.0
# -> chromosomes             2.0
# -> cis distance min     4000.0
# -> cis distance median  4500.0
# -> cis distance max     9000.0
# -> score min               1.0
# -> score median            2.5
# -> score max               5.0
```

The anchors are full `Loci`, so everything on the
[Loci]({{ '/guide/loci/' | relative_url }}) page applies to them. The set of
distinct anchors is one line:

```python
anchors = (P.a + P.b).merge()
anchors
# -> Loci(n=6, chroms=2, sorted)
```

---

## Filtering
{: .sec-purple }

```python
P.filter(min_score=3)                  # needs a score column
# -> Pairs(n=2, cis=2, trans=0, cols=[name, score])
P.filter(max_distance=5000)            # trans pairs (distance inf) go too
# -> Pairs(n=2, cis=2, trans=0, cols=[name, score])
P[P.is_cis & (P.distance < 5000)]      # any mask works
# -> Pairs(n=2, cis=2, trans=0, cols=[name, score])
```

`filter` returns the same object when nothing is dropped, a new `Pairs`
otherwise.

---

## Intersecting anchors with loci: `overlapping` / `anchors_overlap`
{: .sec-purple }

`P.overlapping(loci)` keeps the pairs with at least one anchor overlapping
`loci` — bedtools `pairtobed`; `both=True` requires both anchors. `r` widens
the anchors by `r` bp on each side before the test. `loci` is anything
`as_loci` takes.

```python
cre  = gb.Loci.make("peaks.bed")
prom = gb.as_loci("chr1:800-1,200")
P.overlapping(prom)                    # either anchor
# -> Pairs(n=2, cis=1, trans=1, cols=[name, score])
P.overlapping(prom, both=True)         # both anchors
# -> Pairs(n=0, cis=0, trans=0, cols=[name, score])
P.overlapping(cre, both=True)          # CRE-CRE loops only
# -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
P.overlapping(gb.as_loci("chr1:1,200-1,300"), r=200)     # anchors widened by ±200 bp
# -> Pairs(n=2, cis=1, trans=1, cols=[name, score])
```

`anchors_overlap(loci, r=0)` is the primitive underneath: two boolean masks,
anchor `a` overlaps and anchor `b` overlaps, aligned to the pairs. Combine
them when the two anchors must satisfy different conditions:

```python
genes = gb.Genes.make("gencode.v44.annotation.gtf")
promoters = genes.annot["prom"]
enhancers = cre - promoters

pa, pb = P.anchors_overlap(promoters)          # promoter at a / at b
ea, eb = P.anchors_overlap(enhancers)          # enhancer at a / at b
ep = P[(pa & eb) | (pb & ea)]                  # enhancer-promoter loops
ep
# -> Pairs(n=2, cis=1, trans=1, cols=[name, score])
```

Both take `backend=` for the interval engine; every engine gives the same
masks.

---

## Export
{: .sec-purple }

```python
P.to_bedpe("filtered_loops.bedpe")     # ten tab-separated columns ('.' / 0 fill in name / score)
P.save("loops.parquet")                # parquet; Pairs.load / as_pairs read it back
gb.Pairs.load("loops.parquet")
# -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
P.to_pandas(), P.to_polars(), P.to_arrow(), P.to_numpy()   # frames, or a structured array
```

---

## Pairs files: streaming raw contacts
{: .sec-purple }

The functions below work on **contact files** — one read pair per line, as
written by HiC-Pro, pairtools / 4DN and Juicer — not on BEDPE loops. They
stream the file once in chunks and never build a full matrix.

`read_pairs_chunks(path, format='auto', columns=None, chunksize=2_000_000)`
yields pandas frames of `chrom1, pos1, chrom2, pos2` (categorical
chromosomes). The format is detected from the first data line:

| `format` | Layout | Columns used (0-based) |
|---|---|---|
| `'pairs'` | pairtools / 4DN `.pairs`: `readID chrom1 pos1 chrom2 pos2 strand1 strand2 ...` | 1, 2, 3, 4 |
| `'allvalidpairs'` | HiC-Pro: `readID chr1 pos1 strand1 chr2 pos2 strand2 ...` | 1, 2, 4, 5 |
| `'juicer'` | Juicer medium: `readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 ...`, strands as `+` / `-` or as SAM flags `0` / `16` | 2, 3, 6, 7 |

`#` header lines are skipped, `.gz` is read. For any other layout pass
`columns=(chrom1, pos1, chrom2, pos2)`; an unknown `format` name raises an
error listing the known ones.

```python
from genomeblocks import bedpe
chunk = next(bedpe.read_pairs_chunks("sample.allValidPairs"))
chunk.head(2)
# ->   chrom1  pos1 chrom2  pos2
# -> 0   chr1  8482   chr2  3933
# -> 1   chr2  6000   chr1  1727
bedpe.PAIRS_FORMAT_COLUMNS
# -> {'allvalidpairs': (1, 2, 4, 5), 'pairs': (1, 2, 3, 4), 'juicer': (2, 3, 6, 7)}
```

### Contacts per window: `count_pairs`

`count_pairs(loci, pairs_file)` counts, for every window of `loci`, the
contacts landing in it, split by the chromosome of the partner end. **Both
ends of every pair are counted**: each end inside a window adds one to that
window under the other end's chromosome, so a cis pair with both ends in
windows counts once per end. The result is a DataFrame aligned to the rows of
`loci` — `chrom, start, end, uid`, then one column per partner chromosome
that has at least one contact.

```python
W = gb.Loci.tile_genome({"chr1": 10_000, "chr2": 10_000}, 2500)
df = bedpe.count_pairs(W, "sample.allValidPairs")
# -> [INFO] processed 400 pairs total
df
# ->   chrom  start    end                 uid  chr1  chr2
# -> 0  chr1      0   2500      chr1:0-2500(.)    65    31
# -> 1  chr1   2500   5000   chr1:2500-5000(.)    73    23
# -> 2  chr1   5000   7500   chr1:5000-7500(.)    53    28
# -> 3  chr1   7500  10000  chr1:7500-10000(.)    69    29
# -> 4  chr2      0   2500      chr2:0-2500(.)    38    74
# -> 5  chr2   2500   5000   chr2:2500-5000(.)    33    91
# -> 6  chr2   5000   7500   chr2:5000-7500(.)    19    82
# -> 7  chr2   7500  10000  chr2:7500-10000(.)    21    71
df[["chr1", "chr2"]].to_numpy().sum()          # 400 pairs, both ends in a window: 800
# -> 800
```

`target_chrom="chr8"` keeps only the contacts whose partner is on that
chromosome and returns a single `count` column:

```python
W.count_pairs("sample.allValidPairs", target_chrom="chr2", verbose=False)
# ->   chrom  start    end                 uid  count
# -> 0  chr1      0   2500      chr1:0-2500(.)     31
# -> 1  chr1   2500   5000   chr1:2500-5000(.)     23
# -> 2  chr1   5000   7500   chr1:5000-7500(.)     28
# -> 3  chr1   7500  10000  chr1:7500-10000(.)     29
# -> 4  chr2      0   2500      chr2:0-2500(.)     74
# -> 5  chr2   2500   5000   chr2:2500-5000(.)     91
# -> 6  chr2   5000   7500   chr2:5000-7500(.)     82
# -> 7  chr2   7500  10000  chr2:7500-10000(.)     71
```

`Loci.count_pairs(file, **kw)` is the method form. `format`, `columns` and
`chunksize` are passed to the reader; `verbose=False` silences the progress
line.

{: .warning }
> Windows on one chromosome must not overlap — each position is placed in at
> most one window by binary search. Tiles from `Loci.tile_genome` or merged
> peaks (`L.merge()`) qualify; raw overlapping peaks do not.

### Window-to-window matrix: `count_pairs_2d`

`count_pairs_2d(loci_a, pairs_file, loci_b=None)` returns a
`scipy.sparse.csr_matrix` of shape (rows of `loci_a`) x (rows of `loci_b`).
With `loci_b=None` the matrix is over `loci_a` on both axes and **symmetric**:
both orientations of every pair are counted, so a pair with both ends in
window *i* adds 2 to cell `(i, i)`, and a pair between *i* and *j* adds 1 to
`(i, j)` and 1 to `(j, i)`.

```python
M = bedpe.count_pairs_2d(W, "sample.allValidPairs", verbose=False)
M.shape, M.nnz, M.sum(), (M != M.T).nnz
# -> ((8, 8), 64, 800, 0)
M.toarray()[:4, :4]                                    # the chr1 x chr1 block
# -> array([[18, 22,  8, 17],
# ->        [22, 18, 18, 15],
# ->        [ 8, 18, 10, 17],
# ->        [17, 15, 17, 20]])
```

Two helpers read the matrix back through the windows:

```python
bedpe.pair_2d_to_frame(M, W, W).head(3)                # non-zero cells as a long frame
# ->   chrom1  start1  end1 chrom2  start2  end2  count
# -> 0   chr1       0  2500   chr1       0  2500     18
# -> 1   chr1       0  2500   chr1    2500  5000     22
# -> 2   chr1       0  2500   chr1    5000  7500      8
block, wa, wb = bedpe.pair_2d_block(M, W, W, "chr1", "chr2")   # one dense chromosome pair
block.shape, wa, wb
# -> ((4, 4), Loci(n=4, chroms=1), Loci(n=4, chroms=1))
```

With a second set, the matrix is `loci_a` x `loci_b` and no longer symmetric:

```python
W.count_pairs_2d("sample.allValidPairs", loci_b=W.by_chrom("chr2"), verbose=False).shape
# -> (8, 4)
```

---

## Hand-off to Architecture
{: .sec-purple }

A `Pairs` is what [`Architecture.make`]({{ '/guide/architecture/' | relative_url }})
builds its edges from: every anchor is mapped to the CREs within `r` bp of its
midpoint and each (CRE, CRE) combination becomes an edge. Pass the `Pairs`, a
BEDPE path or a frame — `as_pairs` runs on it either way — so you can filter
first:

```python
cre = gb.Loci.make("peaks.bed")
loops = gb.Pairs.make("loops.bedpe").filter(min_score=2)
A = gb.Architecture.make(cre, loops, r=100)
# -> [INFO] 3 loops | 3 mapped (100.0%) | loci=5, links=3 (0 trans)
A
# -> Architecture(name='Skeleton', loci=5, links=3 [3 cis · 0 trans], edge_props=[w], vertex_props=[])
```

The browser draws a `Pairs` or a BEDPE path as loop arcs
(`gb.browser(..., tracks={"loops": loops})`) — see
[Browser]({{ '/guide/browser/' | relative_url }}).

---

## End-to-end: enhancer-promoter loops
{: .sec-purple }

```python
import genomeblocks as gb

genes = gb.Genes.make("gencode.v44.annotation.gtf", promoter_r=1000)
cre   = gb.Loci.make("atac.narrowPeak")
loops = gb.Pairs.make("loops.bedpe", min_score=1)

promoters = genes.annot["prom"]
enhancers = cre - promoters

pa, pb = loops.anchors_overlap(promoters, r=500)
ea, eb = loops.anchors_overlap(enhancers, r=500)
ep = loops[(pa & eb) | (pb & ea)]

ep["gene_a"], _ = genes.nearest_tss(ep.a)      # the gene at each anchor
ep["gene_b"], _ = genes.nearest_tss(ep.b)
ep.to_bedpe("ep_loops.bedpe")
ep.to_pandas()
```

`ep.a` and `ep.b` are `Loci`, so `nearest_tss` labels each anchor, and the
result is a column on the pairs.

---

## See also
{: .sec-purple }

- [BEDPE API]({{ '/api/bedpe/' | relative_url }}) — every function and argument.
- [Design: BEDPE]({{ '/design/bedpe/' | relative_url }}) — the streaming reader and the window index.
- [Architecture]({{ '/guide/architecture/' | relative_url }}) — loops to a contact graph.
- [Loci]({{ '/guide/loci/' | relative_url }}) — the anchor tables.
