---
title: BEDPE & pairs
parent: Design
layout: default
nav_order: 5
---

# BEDPE & pairs
{: .no_toc }

`Pairs` is the BEDPE table: two row-aligned `Loci`, anchor `a` and anchor
`b`, plus columns. `bedpe` is also the module that streams Hi-C read-pair
files into window counts without loading them into memory.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## Pairs: one table, two anchors
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/bedpe-reader.svg %}
</div><figcaption>
<strong>Every loop file goes through <code>Pairs.make</code>.</strong> The six anchor columns are parsed by the tables backend into two Loci on one Genome, row <em>k</em> of <code>P.a</code> and row <em>k</em> of <code>P.b</code> being the two ends of pair <em>k</em>; name, score and any further columns sit in <code>P.cols</code>. <code>Architecture.make</code>, the browser's arc track and <code>Pairs.overlapping</code> all consume that table through <code>as_pairs()</code> and never parse BEDPE themselves.
</figcaption></figure>

| Attribute | What it holds |
|---|---|
| `P.a`, `P.b` | the two anchors as `Loci` (codes / starts / ends / strands), same `Genome`, same length |
| `P.cols` | `name`, `score` when the file has them, then `col11`, `col12`, ... |
| `P.is_cis` | `a.codes == b.codes` |
| `P.distance` | `|mid2 − mid1|` for cis pairs, `inf` for trans |
| `P.mids1`, `P.mids2` | anchor centres |
| `P.genome` | the shared `Genome` |

```python
import genomeblocks as gb

P = gb.Pairs.make("loops.bedpe")
P
# -> Pairs(n=4, cis=3, trans=1, cols=[name, score])
P.a, P.b.to_records()
# -> (Loci(n=4, chroms=2), [('chr1', 4900, 5100, '-'), ('chr1', 10900, 11100, '+'), ('chr2', 500, 600, '+'), ('chr2', 5000, 5100, '+')])
P.is_cis, P.distance
# -> (array([ True,  True, False,  True]), array([4000., 9000.,   inf, 4500.]))
P[0]                                       # the two anchors of one pair, as LocusViews
# -> (Locus[0](chr1:900-1100(+)), Locus[0](chr1:4900-5100(-)))
P.shape, P.columns
# -> ((4, 10), ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2'])
```

**Reading.** `Pairs.make` sniffs the first data line for the column count
(leading `#`, `track` and `browser` lines are skipped), reads the columns it
needs through `backends.tables.read_columns` (polars when installed, else
pandas, same arrays), encodes both chromosome columns into one `Genome`, and
applies `min_score` / `max_distance` as a row filter. `.gz` files are read
directly.

**Other inputs.** `Pairs.from_frame` takes any table with the BEDPE columns
(by name, or the first six by position); `as_pairs(x)` dispatches a path
(BEDPE or parquet), a `Pairs` or a frame. `to_pandas` / `to_polars` /
`to_arrow` / `to_bedpe` / `save` go back out, and the `TableMixin` protocols
hand the table to polars, DuckDB, seaborn and friends as it is.

```python
P2 = gb.Pairs.from_frame(P.to_polars())
P2.a.equals(P.a) and P2.b.equals(P.b)
# -> True
from genomeblocks.bedpe import as_pairs
as_pairs("loops.bedpe").shape, as_pairs(P.to_pandas()).shape
# -> ((4, 10), (4, 10))
```

## Selections run on the anchors
{: .sec-purple }

Every selection is a boolean mask over the pairs, computed on the anchor
`Loci` with the `intervals` backend, so `backend=` applies here too:

| Call | Mask |
|---|---|
| `P.filter(min_score=, max_distance=)` | `score >= min_score`, `distance <= max_distance` (trans pairs are `inf`, so a distance cap drops them) |
| `P.anchors_overlap(loci, r=0)` | `(a.overlap_any(loci), b.overlap_any(loci))`, anchors widened by `r` |
| `P.overlapping(loci, both=False)` | one anchor (or both) touches `loci`: `bedtools pairtobed` |
| `P.take(idx)`, `P[mask]`, `P.head(n)` | row selections, columns follow |

```python
cre = gb.Loci.make("peaks.bed")
len(P.filter(min_score=3)), len(P.filter(max_distance=5000))
# -> (2, 2)
len(P.overlapping(cre, both=True)), len(P.overlapping("chr1:900-1100"))
# -> (4, 2)
P.anchors_overlap(cre, backend="bioframe")[0]
# -> array([ True,  True,  True,  True])
```

## Who consumes a Pairs
{: .sec-purple }

- **`Architecture.make(loci, bedpe)`** calls `as_pairs`, re-codes the anchors
  onto the CREs' `Genome`, widens each anchor midpoint by `r` and asks the
  `intervals` backend for the CREs under every anchor; see
  [Architecture]({{ '/design/architecture/' | relative_url }}).
- **`browser`, `View`, `igv_html`** draw cis pairs of the shown chromosome as
  arcs between `mids1` and `mids2`, line width following `score`.
- **`Pairs.overlapping`** is the `pairtobed` of the package.

```python
A = gb.Architecture.make(cre, P, r=100, verbose=False)
A2 = gb.Architecture.make(cre, "loops.bedpe", r=100, verbose=False)
list(A) == list(A2)
# -> True
```

## Counting read pairs in windows
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/bedpe-count.svg %}
</div><figcaption>
<strong>One sorted key, one binary search per read end.</strong> Every window is placed on a single genome-wide key, <code>chromosome code · 2⁴⁰ + start</code>, and the keys are sorted once. The pairs file is read in chunks with categorical chromosome columns, so chromosome names are turned into codes once per category rather than once per row. Each read end is then found with one <code>searchsorted</code>, and the counts are a <code>bincount</code>.
</figcaption></figure>

- **Formats.** HiC-Pro `allValidPairs`, 4DN `.pairs` and Juicer medium are
  recognised from the first data line (`format="auto"`), and `columns=(c1,
  p1, c2, p2)` reads any other layout. `read_pairs_chunks` yields pandas
  frames of `chrom1, pos1, chrom2, pos2` with categorical chromosomes,
  `chunksize` rows at a time, so memory stays flat whatever the file size.
- **`count_pairs(windows, pairs)`** returns one row per window (aligned to the
  input rows: `chrom, start, end, uid`) and one column per partner
  chromosome with at least one contact, or a single `count` column with
  `target_chrom=`. Each read end is counted once for its window under the
  other end's chromosome, so a cis pair with both ends in windows adds two
  counts.
- **`count_pairs_2d(windows, pairs, loci_b=None)`** returns a sparse CSR matrix
  of window × window counts, symmetric when `loci_b` is omitted. Both
  orientations of every pair are counted, so a pair inside one window adds 2
  to its diagonal cell. `pair_2d_block` slices one chromosome pair as a dense
  block with its windows, `pair_2d_to_frame` flattens the non-zero cells.
- **Windows** must not overlap within a chromosome (tiles from
  `Loci.tile_genome` never do). Gaps are fine; a read end in a gap is not
  counted.

```python
W = gb.as_loci([("chr1", s, s + 1000) for s in range(0, 10_000, 2000)] + [("chr2", 0, 5000)])
df = W.count_pairs("hic.allValidPairs", verbose=False)
df
# ->   chrom  start   end                uid  chr1  chr2
# -> 0  chr1      0  1000     chr1:0-1000(.)    33    19
# -> 1  chr1   2000  3000  chr1:2000-3000(.)    24    12
# -> 2  chr1   4000  5000  chr1:4000-5000(.)    23    14
# -> 3  chr1   6000  7000  chr1:6000-7000(.)    26    14
# -> 4  chr1   8000  9000  chr1:8000-9000(.)    26    10
# -> 5  chr2      0  5000     chr2:0-5000(.)    61   141

from genomeblocks import bedpe
M = bedpe.count_pairs_2d(W, "hic.allValidPairs", verbose=False)
type(M).__name__, M.shape, M.nnz, (M.toarray() == M.toarray().T).all()
# -> ('csr_matrix', (6, 6), 31, True)
```

## Costs
{: .sec-purple }

- **`Pairs.make`:** one pass of the table parser, two `encode` calls, one
  filter.
- **Selections:** two `overlap_any` calls (anchors against the set) per
  selection; the anchors' lookup indexes are cached like any `Loci`.
- **`count_pairs`:** per chunk, one categorical → code lookup per
  chromosome column, two `searchsorted` over the window keys (one per end)
  and one `bincount`; the pairs file is streamed once. `count_pairs_2d`
  replaces the `bincount` by a `unique` on `(row, col)` keys per chunk and
  sums the chunks' COO triplets at the end.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
