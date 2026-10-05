---
title: Loci
parent: User Guide
layout: default
nav_order: 1
---

# Loci
{: .no_toc }

`Loci` is the table at the heart of `genomeblocks`: genomic intervals as numpy
columns, with set algebra, merge / nearest / point lookups through a swappable
interval engine, and converters to and from every frame library in the field.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is
{: .sec-green }

A `Loci` holds one interval per row in four aligned numpy arrays plus any number
of extra columns:

| Column | Array | dtype | Meaning |
|---|---|---|---|
| `chrom` | `L.codes` | int32 | chromosome *code* into `L.genome` (`L.chroms` decodes to names) |
| `start` | `L.starts` | int64 | 0-based start |
| `end` | `L.ends` | int64 | end, exclusive |
| `strand` | `L.strands` | int8 | 0 `.`, 1 `+`, 2 `-` (`L.strand` gives the symbols) |
| anything else | `L.cols[name]` or `L[name]` | any | `name`, `score`, `signalValue`, a label you add, ... |

Every table is 0-based, half-open: `[start, end)` covers `end - start` bases,
and two intervals overlap when `s1 < e2 and s2 < e1`. A BED file is read as it
is; a GTF is shifted by one on the way in ([Genes]({{ '/guide/genes/' | relative_url }})).

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-halfopen.svg %}
</div><figcaption>
<strong>Half-open coordinates make lengths and overlaps plain arithmetic.</strong> <code>[3, 7)</code> covers bases 3–6, so its length is <code>end − start = 4</code>. <code>[7, 10)</code> starts exactly where it ends: 3 &lt; 10 holds but 7 &lt; 7 does not, so the two do not overlap. <code>[5, 9)</code> shares bases 5 and 6 with it because both <code>s1 &lt; e2</code> and <code>s2 &lt; e1</code> hold. <code>[5, 5)</code> has no bases — it is the point between bases 4 and 5 — and lies inside <code>[3, 7)</code> because 3 &lt; 5 &lt; 7. A GTF row's 1-based start 4 is stored as 3; its end stays 7.
</figcaption></figure>

**The row number is the join key.** A signal cube, a motif matrix, an
annotation vector or an Architecture vertex column computed from a `Loci` has
one row per locus, in this order. Selecting rows (`L[mask]`, `L.take(rows)`)
gives a new table whose rows line up with `mask` / `rows`, so you subset the
other arrays the same way.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/columnar-rows.svg %}
</div><figcaption>
<strong>The row number is the join key.</strong> Every array you build from a <code>Loci</code> — annotation labels, a signal cube, the vertex columns of an Architecture — has one entry per row, in this order, and an Architecture's edges store row numbers rather than uids. Subset them together: <code>L[mask]</code> and <code>cube[mask]</code> keep lining up, and <code>cube[i]</code> is always locus <code>L[i]</code>.
</figcaption></figure>

Each `Loci` has a [`Genome`]({{ '/concepts/' | relative_url }}) mapping
chromosome names to codes. Tables derived from one another share it; two tables
read separately each have their own, and any operation between them re-codes
the second onto the first (an O(n) lookup). Pass `genome=` to several
constructors to share one from the start.

```python
import genomeblocks as gb

peaks = gb.Loci.make("peaks.bed", keep=True)
peaks
# -> Loci(n=7, chroms=2, sorted, cols=[name, score])
peaks.shape, peaks.columns
# -> ((7, 6), ['chrom', 'start', 'end', 'strand', 'name', 'score'])
peaks.codes, peaks.starts[:3], peaks.strands[:3]
# -> (array([0, 0, 0, 0, 0, 1, 1], dtype=int32), array([ 900, 1900, 4900]), array([1, 2, 1], dtype=int8))
peaks.genome
# -> Genome(2 chroms: chr1, chr2)
```

A `Loci` behaves like a DataFrame — `len`, `shape`, `columns`, `head()`,
`tail()`, `describe()`, `L['score']` — and it speaks the Arrow C stream,
dataframe interchange and narwhals protocols, so `pl.DataFrame(L)`,
`pa.table(L)`, `duckdb.sql("select * from L")`, seaborn, plotly and altair
take it as it is. In a notebook it renders as an HTML table.

---

## Building one
{: .sec-green }

### From a BED-like file: `Loci.make`

```python
peaks = gb.Loci.make("peaks.bed")                 # chrom, start, end, strand only
peaks = gb.Loci.make("peaks.bed", keep=True)      # plus name, score (every standard column)
atac  = gb.Loci.make("atac.narrowPeak", keep=["signalValue", "peak"])
atac.columns
# -> ['chrom', 'start', 'end', 'strand', 'signalValue', 'peak']
```

| arg | default | meaning |
|---|---|---|
| `keep` | `False` | `True` keeps every standard column under its usual name (`name`, `score`, `signalValue`, `pValue`, `qValue`, `peak` for narrowPeak; `thickStart` ... `blockStarts` for BED12; `col11` ... beyond); a list keeps those |
| `sort` | `True` | sort into genome order (natural chromosome order, then start, end) |
| `genome` | `None` | a `Genome` to code chromosomes on (shared with other tables) |
| `backend` | `None` | the table parser: `'polars'` (default when installed) or `'pandas'`; both give the same columns |

`.gz` files are fine. Leading `#`, `track` and `browser` lines are skipped.
Coordinates and strand are always read; column 6 is the strand when present.
narrowPeak / broadPeak columns are typed (`peak` int, `signalValue` float).

```python
summits = gb.Loci(atac.codes, atac.starts + atac["peak"], atac.starts + atac["peak"] + 1,
                  genome=atac.genome)             # 1-bp summits, same rows
```

{: .note }
> `Loci.make` is for intervals. A `.gtf` / `.gff3` raises *"use Genes.make()"*
> and a `.bedpe` raises *"use Pairs.make()"* — see
> [Genes]({{ '/guide/genes/' | relative_url }}) and [BEDPE]({{ '/guide/bedpe/' | relative_url }}).

### From a frame: `Loci.from_frame`

Any pandas, polars (eager or lazy), pyarrow, bioframe or PyRanges frame, a dict
of columns or a structured array. Columns are found by name, case-insensitive:
`chrom` / `chromosome` / `chr` / `seqnames` / `Chromosome` / `#chrom`;
`start` / `Start` / `chromStart` / `begin`; `end` / `End` / `chromEnd` / `stop`;
`strand` / `Strand`. A header-less frame whose first three columns are text,
int, int is read by position. Row order is kept.

```python
import pandas as pd, polars as pl

df = pd.DataFrame({"Chromosome": ["chr1", "chr1", "chr2"], "Start": [100, 400, 50],
                   "End": [200, 450, 90], "score": [1.5, 2.0, 0.5]})
L = gb.Loci.from_frame(df)
L
# -> Loci(n=3, chroms=2, cols=[score])
gb.Loci.from_frame(pl.DataFrame({"chrom": ["chr1"], "start": [5], "end": [9], "strand": ["-"]}).lazy()).to_records()
# -> [('chr1', 5, 9, '-')]
```

`chrom=`, `start=`, `end=`, `strand=` name the columns when the names are
unusual; `keep=` works as in `make`. A frame without recognisable columns
raises an error that lists the accepted spellings. Missing values in `chrom`,
`start` or `end` raise too — drop or fill those rows first.

`from_pandas`, `from_polars`, `from_arrow`, `from_bioframe` are the same
function; `from_pyranges`, `from_bedtool` and `from_anndata` (the regions of
`adata.var`, from columns or from `chr1:100-200` names) complete the set.

### From anything: `gb.as_loci`

Every public function in `genomeblocks` calls `as_loci` on its interval inputs,
so wherever a `Loci` is expected you can hand over the raw thing:

| Input | Example |
|---|---|
| a `Loci` | returned as is |
| a path | `.bed`, `.narrowPeak`, `.broadPeak`, `.bedGraph`, `.csv` / `.tsv` / `.txt` (with or without a header), `.parquet`, `.gz` of any of these |
| a region string | `'chr1:1,000-2,000'`, `'chr8:127.7-128.1 Mb'`, `'chr2:5kb-12kb'` |
| a `Locus` | one row |
| a frame | pandas, polars, pyarrow, bioframe, PyRanges, BedTool, dict of columns, structured array, any Arrow / interchange object |
| an AnnData | its `var` regions |
| a tuple of arrays | `(chroms, starts, ends[, strands])` |
| a list of records | `Locus` objects, `(chrom, start, end[, strand])` tuples, region strings, uids |
| a list of any of the above | concatenated in order |

```python
gb.as_loci("chr1:1,000-2,000").to_records()
# -> [('chr1', 1000, 2000, '.')]
gb.as_loci([("chr1", 5, 15, "+"), "chr2:100-200", gb.Locus("chr1", 50, 60)]).uid
# -> array(['chr1:5-15(+)', 'chr2:100-200(.)', 'chr1:50-60(.)'], dtype=object)
gb.as_loci({"chrom": ["chr1", "chr1"], "start": [1, 7], "end": [5, 9]}).to_records()
# -> [('chr1', 1, 5, '.'), ('chr1', 7, 9, '.')]
gb.as_loci("peaks.bed")                           # file order, every column kept
# -> Loci(n=7, chroms=2, cols=[name, score])
```

{: .tip }
> `as_loci(path)` keeps the file's row order and all its columns;
> `Loci.make(path)` sorts and keeps only coordinates unless you ask. Use
> `make` when you are building a table to work with, `as_loci` when you are
> passing a file straight into a function.

`Loci.from_records(items)` and `Loci.from_uids(uids)` build from lists of
records or `chrom:start-end(strand)` strings; `Loci.load(path)` reads a parquet
written by `save`.

---

## Rows, columns and indexing
{: .sec-green }

`L[i]` is a `LocusView`: a real `Locus` (`chrom`, `start`, `end`, `strand`,
`uid`, `length`, `center`, `overlaps`, `distance_to`) whose fields read the
columns. Nothing is copied, and extra columns are one attribute away. Writing
to a view writes through to the table.

```python
peaks = gb.Loci.make("peaks.bed", keep=True)
peaks[0]
# -> Locus[0](chr1:900-1100(+))
peaks[0].chrom, peaks[0].start, peaks[0].name, peaks[0].score, peaks[0].row
# -> ('chr1', 900, 'p1', 10.0, 0)
for l in peaks.head(2):
    print(l.uid, l.length, l.center)
# -> chr1:900-1100(+) 200 1000
# -> chr1:1900-2100(-) 200 2000
```

`L[key]` dispatches on the key:

| Key | Returns |
|---|---|
| `L[3]` | the row as a `LocusView` |
| `L['start']`, `L['score']` | the column as a numpy array (`chrom` / `start` / `end` / `strand` or any extra column) |
| `L['chr1:900-1100(+)']` | the row with that uid |
| `L[['chr2:500-600(+)', 'chr1:900-1100(+)']]` | those rows, in that order, as a `Loci` |
| `L[mask]`, `L[[0, 2]]`, `L[1:3]` | a new `Loci` of those rows (`L.take(idx)` does the same) |

```python
peaks["start"][:3], peaks["name"][:3]
# -> (array([ 900, 1900, 4900]), array(['p1', 'p2', 'p3'], dtype=object))
peaks[peaks["score"] > 40].to_records()
# -> [('chr1', 10900, 11100, '-'), ('chr2', 500, 600, '+'), ('chr2', 5000, 5100, '+')]
peaks[["chr2:500-600(+)", "chr1:900-1100(+)"]].to_records()
# -> [('chr2', 500, 600, '+'), ('chr1', 900, 1100, '+')]
"chr1:900-1100(+)" in peaks, peaks.row("chr1:900-1100(+)")
# -> (True, 0)
```

`L['name'] = values` adds or replaces a column: one value per row, or a scalar
broadcast to every row. Coordinates cannot be assigned this way — build a new
`Loci` (for example with `slop`, or from the arrays) instead.

```python
import numpy as np
peaks["log_score"] = np.log10(peaks["score"])
peaks["set"] = "atac"
peaks.columns
# -> ['chrom', 'start', 'end', 'strand', 'name', 'score', 'log_score', 'set']
```

Derived columns are properties: `L.chroms` (names), `L.strand` (symbols),
`L.centers`, `L.lengths`, `L.uid` (`chrom:start-end(strand)`, built once),
`L.names` (`chrom:start-end`, the scATAC region name), `L.uids` (`{uid: row}`).

`head(n)`, `tail(n)` return a `Loci`; `describe()` a one-column summary:

```python
peaks.describe()
# ->                              value
# -> rows                             7
# -> chromosomes                      2
# -> bases covered (merged)        1100
# -> length min                     100
# -> length 25%                   100.0
# -> length median                200.0
# -> length 75%                   200.0
# -> length max                     200
# -> length mean             157.142857
# -> strand + / - / .         4 / 2 / 1
# -> genome-sorted                 True
```

{: .warning }
> Each `Loci` caches its uid map and its lookup indexes. Writing through a
> `LocusView` (`peaks[0].start = 950`) or assigning `L.starts = ...` directly
> invalidates them on the next access, but a `Loci` is meant to be built, not
> edited: prefer making a new table over mutating one.

---

## Set algebra
{: .sec-green }

The operators work on **overlap**: a row of `a` is kept when it overlaps at
least one row of `b` under the half-open rule (`s1 < e2 and s2 < e1`, same
chromosome). Book-ended intervals (`[100, 200)` and `[200, 300)`) do not
overlap. Rows come back whole, in their original order, with their columns.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-setops.svg %}
</div><figcaption>
<strong>Set operators keep whole rows.</strong> <code>a &amp; b</code> keeps each row of <code>a</code> that touches any row of <code>b</code>; <code>a - b</code> keeps the rest. Neither clips, splits or reorders an interval, so a mask over <code>a</code> (<code>a.overlap_any(b)</code>) gives the same rows. <code>(a + b).merge()</code> is the one call here that changes coordinates: overlapping or book-ended rows fuse into one.
</figcaption></figure>

```python
a = gb.as_loci([("chr1", 100, 200), ("chr1", 300, 400), ("chr1", 500, 600), ("chr2", 0, 100)])
b = gb.as_loci([("chr1", 200, 300), ("chr1", 350, 360), ("chr2", 50, 60)])

(a & b).to_records()        # rows of a that overlap b          (bedtools intersect -u)
# -> [('chr1', 300, 400, '.'), ('chr2', 0, 100, '.')]
(a - b).to_records()        # rows of a that overlap nothing in b (bedtools intersect -v)
# -> [('chr1', 100, 200, '.'), ('chr1', 500, 600, '.')]
(a | b).to_records()        # concatenation: every row of a, then every row of b
# -> [('chr1', 100, 200, '.'), ('chr1', 300, 400, '.'), ('chr1', 500, 600, '.'), ('chr2', 0, 100, '.'), ('chr1', 200, 300, '.'), ('chr1', 350, 360, '.'), ('chr2', 50, 60, '.')]
(a ^ b).to_records()        # (a - b) + (b - a)
# -> [('chr1', 100, 200, '.'), ('chr1', 500, 600, '.'), ('chr1', 200, 300, '.')]
```

| Operator | Method | Meaning |
|---|---|---|
| `a & b` | `a.intersect(b)` | rows of `a` overlapping any row of `b` |
| `a - b`, `a / b` | `a.difference(b)` | rows of `a` overlapping no row of `b` |
| `a + b`, `a \| b` | `a.union(b)` | concatenation (duplicates kept; `merge()` collapses them). Only the columns both sides have are kept |
| `a ^ b` | | rows of either set that do not overlap the other |

The named methods take `backend=`. Two lower-level calls return arrays instead
of a table:

```python
a.overlap_any(b)            # bool per row of a
# -> array([False,  True, False,  True])
a.overlap_pairs(b)          # (rows of a, rows of b) for every overlapping pair, sorted
# -> (array([1, 3]), array([1, 2]))
```

`overlap_pairs` is what you want when one row can hit several (the 1.x `map`);
`overlap_any` is the mask behind `&` and `-`.

The right-hand side can be anything `as_loci` takes — a frame, a path, a region
string, a list of tuples — and it is re-coded onto the left side's `Genome`
automatically:

```python
len(a & b.to_pandas()), (a - [("chr1", 150, 160)]).to_records()
# -> (2, [('chr1', 300, 400, '.'), ('chr1', 500, 600, '.'), ('chr2', 0, 100, '.')])
```

---

## slop, sort, merge
{: .sec-green }

```python
a.slop(50).to_records()                    # ± n bp on each side, clipped at 0
# -> [('chr1', 50, 250, '.'), ('chr1', 250, 450, '.'), ('chr1', 450, 650, '.'), ('chr2', 0, 150, '.')]

u = gb.as_loci([("chr10", 5, 6), ("chr2", 5, 6), ("chr1", 9, 10), ("chr1", 1, 2), ("chrX", 0, 1)])
u.is_sorted
# -> False
u.sort().to_records()                      # natural chromosome order, then start, end
# -> [('chr1', 1, 2, '.'), ('chr1', 9, 10, '.'), ('chr2', 5, 6, '.'), ('chr10', 5, 6, '.'), ('chrX', 0, 1, '.')]

m = gb.as_loci([("chr1", 10, 20, "+"), ("chr1", 20, 30, "-"), ("chr1", 25, 28), ("chr1", 40, 50, "-"), ("chr2", 0, 5)])
m.merge().to_records()                     # overlapping or book-ended rows fuse
# -> [('chr1', 10, 30, '+'), ('chr1', 40, 50, '-'), ('chr2', 0, 5, '.')]
```

`merge()` has bedtools-merge semantics: overlapping **or book-ended** rows fuse
into one; the result is sorted, carries the strand of the first row of each
block, and drops the extra columns. It does not need sorted input. Chaining is
idiomatic:

```python
cre = gb.Loci.make("peaks.bed").slop(100).merge()
```

`sort()` returns a new table in genome order and marks it sorted. Sorted tables
expose their chromosome blocks:

```python
peaks.is_sorted, peaks.chrom_offsets
# -> (True, {'chr1': (0, 5), 'chr2': (5, 7)})
peaks.by_chrom("chr2").to_records()        # a slice on a sorted table, a selection otherwise
# -> [('chr2', 500, 600, '+'), ('chr2', 5000, 5100, '+')]
```

---

## nearest
{: .sec-green }

`q.nearest(r)` returns two arrays aligned to the rows of `q`: the row of `r`
that is nearest, and the distance.

```python
q = gb.as_loci([("chr1", 100, 110), ("chr1", 200, 210), ("chr1", 150, 155), ("chr3", 0, 10)])
r = gb.as_loci([("chr1", 110, 120), ("chr1", 150, 160)])
rows, dist = q.nearest(r)
rows, dist
# -> (array([ 0,  1,  1, -1]), array([ 0, 40,  0, -1]))
```

| Case | row | distance |
|---|---|---|
| the query overlaps a row of `r` | the overlapping row with the lowest start | `0` |
| book-ended (`[100, 110)` next to `[110, 120)`) | that row | `0` |
| a gap | the closer of the left and right neighbour (the left one on an exact tie) | the gap in bases (`[200, 210)` to `[150, 160)` is 40) |
| nothing on the query's chromosome | `-1` | `-1` |

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-nearest.svg %}
</div><figcaption>
<strong><code>nearest()</code> measures the gap between interval ends.</strong> <code>q[0] = [100, 110)</code> is book-ended with <code>r[0] = [110, 120)</code>: row 0, distance 0. <code>q[2] = [150, 155)</code> lies inside <code>r[1]</code>: row 1, distance 0. <code>q[1] = [200, 210)</code> touches nothing; the gap to <code>r[1]</code> is 200 − 160 = 40 and the gap to <code>r[0]</code> would be 80, so the closer row wins (the left one on an exact tie). <code>q[3]</code> is on chr3, where <code>r</code> has no rows: row −1, distance −1. Both output arrays are in <code>q</code>'s row order.
</figcaption></figure>

Distances are gaps between interval ends, not centre-to-centre. Rows with
`-1` must be masked before indexing `r`:

```python
hit = rows >= 0
r_uid = np.where(hit, r.uid[np.maximum(rows, 0)], "")
```

[`Genes.nearest_tss`]({{ '/guide/genes/' | relative_url }}) is this call against
the promoter windows of an annotation.

---

## Point lookups
{: .sec-green }

`overlaps` returns the rows overlapping **one** region as a `Loci`;
`overlap_rows` returns their row numbers. The region can be a string, a
`Locus`, or `(chrom, start, end)` arguments. The lookup index is built on the
first call and cached until the rows change, so repeated lookups are cheap.

```python
peaks.overlaps("chr1:1,000-2,000").to_records()
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '-')]
peaks.overlap_rows("chr1", 1000, 2000), peaks.overlap_rows(gb.Locus("chr1", 1000, 2000))
# -> (array([0, 1]), array([0, 1]))
peaks.overlaps("chr1:4.9-5.1 kb").to_records()
# -> [('chr1', 4900, 5100, '+')]
peaks.overlap_rows("chr9:1-100")           # an unknown chromosome is simply empty
# -> array([], dtype=int64)
```

Region strings take thousands separators and `kb` / `Mb` units on either
number: `'chr1:1,000-2,000'`, `'chr8:127.7-128.1 Mb'`, `'chr2:5kb-12kb'`.
`gb.Locus.parse(text)` gives the `Locus` for one.

---

## Tiles
{: .sec-green }

Uniform windows over a chromosome or a genome, for binning signal or counting
contacts. The last tile of each chromosome is clipped to its length. Sizes come
from a `{chrom: length}` dict, a `.chrom.sizes` path, a pandas Series, a
`Genome`, or a cooler.

```python
gb.Loci.tile_genome({"chr1": 1000, "chr2": 250}, 400).to_records()
# -> [('chr1', 0, 400, '.'), ('chr1', 400, 800, '.'), ('chr1', 800, 1000, '.'), ('chr2', 0, 250, '.')]
gb.Loci.tile_genome("genome.chrom.sizes", 5000)
# -> Loci(n=6, chroms=2)
gb.Loci.tile_genome(gb.Genome.from_fasta("genome.fa"), 5000, chroms=["chr2"]).to_records()
# -> [('chr2', 0, 5000, '.'), ('chr2', 5000, 8000, '.')]
gb.Loci.tile("chr2", 3000, "genome.chrom.sizes").to_records()     # one chromosome; a length works too
# -> [('chr2', 0, 3000, '.'), ('chr2', 3000, 6000, '.'), ('chr2', 6000, 8000, '.')]
```

---

## Sequences
{: .sec-green }

`sequences` fetches the bases of every row — or of `center ± r` — from a FASTA
path (its `.fai` is used when present, else built), a `{chrom: seq}` dict, or
an open pysam / pyfaidx / Biopython handle. Windows are clipped to the
chromosome, so a sequence can come back shorter; `strand=True`
reverse-complements `-` rows; `upper=True` upper-cases.

```python
seqs = peaks.sequences("genome.fa")
[len(s) for s in seqs[:3]]
# -> [200, 200, 200]
site = gb.as_loci([("chr1", 1000, 1012, "-")])
site.sequences("genome.fa"), site.sequences("genome.fa", strand=True)
# -> (['ACGTACGTACGT'], ['ACGTACGTACGT'])
peaks.sequences("genome.fa", r=5)[:2]      # 10 bp around each centre
# -> ['TCCCTACGTA', 'TTATCGACCA']
```

`to_fasta(path, fasta, r=None)` writes them as FASTA with uids as headers;
`to_seqrecords(fasta)` gives Biopython `SeqRecord`s. The engine is the fasta
backend (`backend='pysam'`, `'pyfaidx'`, `'biopython'`, `'memory'`); the
default is the package's own indexed reader. `gb.read_fasta(path)` opens a
handle you can pass around, and `gb.Genome.from_fasta(path)` reads chromosome
names and sizes:

```python
peaks.to_fasta("peaks.fa", "genome.fa", r=10)
# -> 'peaks.fa'
gb.Genome.from_fasta("genome.fa").sizes
# -> {'chr1': 20000, 'chr2': 8000}
```

---

## Liftover
{: .sec-green }

`liftover` maps a table to another assembly through a UCSC chain file
(pyliftover: `pip install pyliftover`). Each end that lands in a chain gap walks
inward by up to `1 - min_match` of the interval length before giving up, the
UCSC base-fraction rule. Rows that do not lift are dropped; the new table's
`source_row` column says where each row came from, so other arrays can be
subset to match.

```python
lifted = peaks.liftover("old_to_new.chain", min_match=0.95)
# -> [INFO] liftover: 7 → 5 (2 dropped, min_match=0.95)
lifted["source_row"]
# -> array([0, 1, 2, 3, 4])
```

---

## Save, load, convert
{: .sec-green }

`save` writes parquet (Arrow); `load` reads it back onto any `Genome`. Extra
columns round-trip, including nulls in text columns. `to_bed` writes BED6 —
the name is `cols[name]` or the uid, the score `cols[score]` or 0:

```python
peaks.save("peaks.parquet")
gb.Loci.load("peaks.parquet").equals(peaks, cols=True)
# -> True
peaks.to_bed("peaks6.bed", name="name", score="score")
print(peaks.to_bed().splitlines()[0])            # no path: the text
# -> chr1	900	1100	chr1:900-1100(+)	0	+
```

Converters, in and out (see [Interoperability]({{ '/interoperability/' | relative_url }})
for the details and what copies):

| Library | In | Out |
|---|---|---|
| pandas | `Loci.from_frame(df)` | `L.to_pandas()` — `chrom` and `strand` categorical, extra columns follow; `uid=True` adds the uid |
| polars | `Loci.from_frame(df)` (lazy frames too) | `L.to_polars()` |
| pyarrow | `Loci.from_frame(table)` | `L.to_arrow()` — dictionary-encoded `chrom` / `strand` |
| bioframe | `Loci.from_frame(df)` | `L.to_bioframe()` — plain string `chrom` |
| pyranges | `Loci.from_pyranges(gr)` | `L.to_pyranges()` — `Chromosome / Start / End / Strand` (0.x groups rows by chromosome: sort first if order matters) |
| pybedtools | `Loci.from_bedtool(bt)` | `L.to_bedtool()` — BED6 |
| cgranges | — | `L.to_cgranges()` — a built index, label = row |
| AnnData | `Loci.from_anndata(adata)` | `L.to_anndata(X, obs=...)` — loci as `var`, `X` is obs x loci |
| numpy | `as_loci((chroms, starts, ends))` | `L.to_numpy()` — a structured array; `np.asarray(L)` |
| records | `Loci.from_records`, `from_uids` | `L.to_records()` — `(chrom, start, end, strand)` tuples |
| BED | `Loci.make(path)` | `L.to_bed(path)` |
| parquet | `Loci.load(path)` | `L.save(path)` |

```python
peaks.to_pandas().dtypes.tolist()[:4]
# -> [CategoricalDtype(categories=['chr1', 'chr2'], ordered=False, categories_dtype=object), dtype('int64'), dtype('int64'), CategoricalDtype(categories=['.', '+', '-'], ordered=False, categories_dtype=object)]
X = np.random.default_rng(0).random((3, len(peaks)))
peaks.to_anndata(X, obs=["s1", "s2", "s3"])
# -> AnnData object with n_obs × n_vars = 3 × 7
# ->     var: 'chrom', 'start', 'end', 'strand', 'name', 'score', 'log_score', 'set'
```

And through the protocols, with no conversion call at all:

```python
import duckdb, polars as pl
pl.DataFrame(peaks).shape
# -> (7, 8)
duckdb.sql("select chrom, count(*) n from peaks group by chrom order by chrom").fetchall()
# -> [('chr1', 5), ('chr2', 2)]
```

---

## Picking the interval engine
{: .sec-green }

Overlap, merge, nearest and point lookups run through the **intervals**
backend. The default is the package's own numpy kernels, which install with
`pip install genomeblocks`. The others are there for comparison or because a
pipeline already standardises on one; every engine returns the same rows
(the test suite checks them against each other, and
[Benchmarks]({{ '/benchmarks/' | relative_url }}) times them).

| Backend | overlap | nearest | merge | point lookup |
|---|---|---|---|---|
| `genomeblocks` (default) | yes | yes | yes | yes |
| `cgranges` | yes | — | — | yes |
| `ncls` | yes | — | — | yes |
| `bioframe` | yes | yes | yes | yes |
| `pyranges` | yes | yes | yes | yes |
| `bedtools` | yes | yes | yes | yes |

Pass `backend=` on a call, or wrap a stretch of code:

```python
a.intersect(b, backend="bioframe").equals(a & b)
# -> True
with gb.use_backend(intervals="pyranges"):
    hits = a & b                                  # every interval call in the block
```

A requested backend that is missing raises `ImportError` with the install
command; one that cannot do the operation raises `NotImplementedError` naming
the ones that can. Nothing falls back silently.

```python
try:
    a.intersect(b, backend="cgranges")
except ImportError as e:
    print(e)
# -> the 'cgranges' intervals backend is not installed: conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)
try:
    a.nearest(b, backend="ncls")
except NotImplementedError as e:
    print(e)
# -> the 'ncls' intervals backend has no nearest; use one of: genomeblocks, bioframe, pyranges, bedtools
```

`gb.backends()` lists every family with what is installed and in use; see
[Backends]({{ '/backends/' | relative_url }}).

---

## Example: a CRE catalogue
{: .sec-green }

```python
import genomeblocks as gb

atac    = gb.Loci.make("atac.narrowPeak", keep=["signalValue"])
h3k27ac = gb.Loci.make("H3K27ac.bed")
h3k4me1 = gb.Loci.make("H3K4me1.bed")

# CRE = ATAC peak widened by 100 bp, merged
cre = atac.slop(100).merge()

# tiers by mark presence
active   = cre & h3k27ac
primed   = (cre & h3k4me1) - h3k27ac
inactive = cre - h3k27ac - h3k4me1

# one label column on the catalogue, by row
cre["tier"] = "inactive"
cre["tier"][cre.overlap_any(h3k4me1)] = "primed"
cre["tier"][cre.overlap_any(h3k27ac)] = "active"

cre.to_bed("cre.bed", name="tier")
cre.save("cre.parquet")
```

The label column lines up with the rows, so a signal cube pulled with
`cre.signal([...])` or a motif matrix from `cre.scan_motifs_matrix(...)` can be
grouped by `cre["tier"]` straight away — see
[Signal]({{ '/guide/signal/' | relative_url }}) and
[Motifs]({{ '/guide/motifs/' | relative_url }}).

---

## See also
{: .sec-green }

- [Loci API]({{ '/api/loci/' | relative_url }}) — every method and argument.
- [Locus API]({{ '/api/locus/' | relative_url }}) — the one-interval class and region parsing.
- [Design: Loci]({{ '/design/loci/' | relative_url }}) — how the kernels and the cached indexes work.
- [Concepts]({{ '/concepts/' | relative_url }}) — tables, the `Genome`, rows as the join key, coordinates.
- [Interoperability]({{ '/interoperability/' | relative_url }}) and [Backends]({{ '/backends/' | relative_url }}).
