---
title: locus
parent: API Reference
layout: default
nav_order: 1
---

# `genomeblocks.locus`
{: .no_toc }

One genomic interval, and the region parsers every entry point uses. Tables
([`Loci`]({{ '/api/loci/' | relative_url }}), [`Genes`]({{ '/api/genes/' | relative_url }}))
hold intervals as columns; a `Locus` is what you get when you look at one row
(`L[i]` is a `LocusView`, a `Locus` whose fields read the columns) and what you
pass to say "this region". Coordinates are 0-based, half-open everywhere, as in
[Concepts]({{ '/concepts/' | relative_url }}).
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
from genomeblocks import Locus
from genomeblocks.locus import parse_region, parse_regions
```

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Locus`
{: .sec-navy }

```python
@dataclass
class Locus:
    chrom: str
    start: int          # 0-based
    end: int            # exclusive
    strand: str = "."   # '.', '+' or '-'
```

A plain dataclass: four fields, no index, no genome. Two loci are equal when
their `uid` is equal, so a `Locus` works as a dict key or set member.

### Properties

| Name | Type | Description |
|---|---|---|
| `uid` | `str` | `"chrom:start-end(strand)"`; `Locus.from_uid` reverses it. |
| `length` | `int` | `end - start`. |
| `center` | `int` | `(start + end) // 2`. |

### Constructors

```python
Locus.from_uid(uid) -> Locus    # 'chr1:100-200(+)' or 'chr1:100-200' (strand '.')
Locus.parse(region) -> Locus    # anything parse_region reads; strand kept from a Locus-like input
```

```python
a = Locus("chr1", 100, 200, "+")
a.uid, a.length, a.center
# -> ('chr1:100-200(+)', 100, 150)
Locus.from_uid("chr1:100-200(+)") == a
# -> True
Locus.from_uid("chr1:100-200")
# -> Locus(chrom='chr1', start=100, end=200, strand='.')
Locus.parse("chr8:127.7-128.1 Mb")
# -> Locus(chrom='chr8', start=127700000, end=128100000, strand='.')
Locus.parse(("chr3", 5, 9))
# -> Locus(chrom='chr3', start=5, end=9, strand='.')
```

### Methods

```python
Locus.overlaps(other) -> bool       # same chromosome and intersecting half-open spans
Locus.distance_to(other) -> int     # |center - other.center|; NotImplemented across chromosomes
Locus.sequence(fasta, r=None) -> str
Locus.copy() -> Locus               # a plain Locus (also from a LocusView)
```

`sequence` reads the bases of the locus, or of `center ± r` when `r` is
given, from anything the [fasta backend]({{ '/api/backends/' | relative_url }})
opens: a FASTA path (indexed on first use), a `{chrom: str}` dict (what
`gb.read_fasta` returns), or an open pyfaidx / pysam / Biopython handle. The
window is clipped at 0.

```python
b = Locus("chr1", 150, 250)
a.overlaps(b), a.overlaps(Locus("chr1", 200, 300))      # [100, 200) and [200, 300) do not touch
# -> (True, False)
a.distance_to(Locus("chr1", 300, 400))
# -> 200
Locus("chr1", 1000, 1012).sequence("genome.fa")
# -> 'ACGTACGTACGT'
Locus("chr1", 1006, 1006).sequence("genome.fa", r=6)   # the same 12 bases around the centre
# -> 'ACGTACGTACGT'
Locus("chr1", 1000, 1012).sequence({"chr1": "A" * 2000})
# -> 'AAAAAAAAAAAA'
```

### Ordering and hashing

`__lt__`, `__le__`, `__gt__`, `__ge__` compare by `(chrom, start)`, so a list
of `Locus` objects sorts into chromosome-then-position order (chromosome
names compare as strings here; `Loci.sort()` uses natural order, chr2 before
chr10). `__eq__` and `__hash__` use the `uid`, so strand matters for equality.

```python
sorted([Locus("chr2", 5, 9), Locus("chr1", 50, 60), Locus("chr1", 5, 9)])
# -> [Locus('chr1', 5, 9), Locus('chr1', 50, 60), Locus('chr2', 5, 9)]
a == Locus("chr1", 100, 200, "+"), a != b, hash(a) == hash(Locus("chr1", 100, 200, "+"))
# -> (True, True, True)
```

{: .note }
A `Locus` is accepted wherever a table is: `gb.as_loci(Locus("chr1", 5, 9, "-"))`
is a one-row `Loci`, and `cre.overlaps(locus)`, `cre.nearest(locus)`,
`genes.labels(locus)` all take one.

---

## `LocusView(Locus)`
{: .sec-navy }

```python
from genomeblocks.loci import LocusView
LocusView(owner: Loci, i: int)       # made by Loci[i], Loci.__iter__, Pairs[i]; not built by hand
```

A `Locus` whose fields live in a `Loci`: nothing is copied. Reading `chrom`,
`start`, `end`, `strand` reads the table's columns at row `i`; writing them
writes through to the table (an unseen chromosome name is added to the table's
`Genome`) and invalidates its cached uids and indexes. Extra columns of the
table are attributes (`L[i].score`, `L[i].name`), returned as Python scalars.

| Name | Description |
|---|---|
| `row` | The row number in the owning table. |
| `chrom`, `start`, `end`, `strand` | Read / write-through properties. |
| `uid`, `length`, `center`, `overlaps`, `distance_to`, `sequence`, ordering | Inherited from `Locus`. |
| `copy()` | A detached, plain `Locus`. |
| `repr` | `Locus[row](uid)`. |

```python
cre = gb.Loci.make("peaks.bed", keep=True)
v = cre[0]
v, type(v).__name__, isinstance(v, Locus)
# -> (Locus[0](chr1:900-1100(+)), 'LocusView', True)
v.row, v.name, v.score                 # extra columns come along
# -> (0, 'p1', 10.0)
v.start = 850                          # writes into cre.starts[0]
cre.starts[0], cre[0]
# -> (850, Locus[0](chr1:850-1100(+)))
v.strand = "-"
cre.strand[:2]
# -> array(['-', '-'], dtype=object)
cre[-1]                                # negative indices count from the end
# -> Locus[6](chr2:5000-5100(+))
cre[2].copy()                          # a plain Locus, detached from the table
# -> Locus(chrom='chr1', start=4900, end=5100, strand='+')
```

{: .warning }
Writing through a `LocusView` is for the odd fix-up. Edits that change many
rows belong in a new table (`Loci.take`, `slop`, `Loci.from_frame`), and the
coordinate columns of a `Loci` cannot be replaced with `L['start'] = ...`.

---

## `parse_region(region) -> (chrom, start, end)`
{: .sec-green }

```python
parse_region(region) -> Tuple[str, int, int]
```

One region from a string, a `(chrom, start, end)` tuple or list, or any
object with `chrom` / `start` / `end` attributes. Strings take thousands
separators (`,` or `_`), an en dash or a hyphen, and `kb` / `Mb` units on
either number (the unit of the second number applies to the first when it
has none). Everything that takes a region (`Loci.overlaps`, `gb.browser`,
`gb.as_loci("chr1:1-2,000")`, `gb.coverage`, `gb.igv_html`) goes through it.

```python
parse_region("chr1:1,000-2,000")
# -> ('chr1', 1000, 2000)
parse_region("chr8:127.7-128.1 Mb")
# -> ('chr8', 127700000, 128100000)
parse_region("chr2:5kb-12kb")
# -> ('chr2', 5000, 12000)
parse_region("chr1:1_000–2_000")
# -> ('chr1', 1000, 2000)
parse_region(("chr3", 5, 9))
# -> ('chr3', 5, 9)
parse_region(Locus("chrX", 1, 2, "-"))
# -> ('chrX', 1, 2)
parse_region("chr1")
# -> ValueError: cannot parse region 'chr1'; use 'chr1:1,000-2,000', 'chr1:1.2-1.5 Mb' or (chrom, start, end)
```

---

## `parse_regions(text) -> list`
{: .sec-green }

```python
parse_regions(text) -> list[Tuple[str, int, int]]
```

Every region in a string — `'chr8:127.7-128.0 Mb chr1:1-2 Mb'` holds two,
which is how the [browser]({{ '/api/browser/' | relative_url }}) gets a split
view — or one region per item of a list of strings, tuples and `Locus`
objects. A single tuple or `Locus` gives a one-item list.

```python
parse_regions("chr8:127.7-128.0 Mb chr1:1,000-2,000")
# -> [('chr8', 127700000, 128000000), ('chr1', 1000, 2000)]
parse_regions(["chr1:1-2", ("chr2", 3, 4), Locus("chr3", 5, 6)])
# -> [('chr1', 1, 2), ('chr2', 3, 4), ('chr3', 5, 6)]
parse_regions(("chr1", 1, 2))
# -> [('chr1', 1, 2)]
```

{: .tip }
To turn regions into a table rather than tuples, use
`gb.Loci.from_records(["chr1:1-2", ("chr2", 3, 4)])` or `gb.as_loci([...])`;
both accept the same spellings, plus uids such as `'chr1:5-15(+)'`.
