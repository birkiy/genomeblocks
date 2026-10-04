---
title: Loci & Locus
parent: Design
layout: default
nav_order: 1
---

# Loci & Locus
{: .no_toc }

`Locus` is one interval; `Loci` is a list of them with two indexes it builds on
demand. Everything else in the package is built on these two classes.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## Locus
{: .sec-green }

A `Locus` is a small dataclass: `chrom`, `start`, `end` (half-open, BED
convention) and `strand` (`.` by default). Three derived values do most of the
work:

| Attribute | Value | Used for |
|---|---|---|
| `uid` | `"chr1:100-200(+)"` | equality, hashing, dictionary keys, graph vertices |
| `center` | `(start + end) // 2` | signal windows, motif windows, loop distances |
| ordering | by `chrom`, then `start` | `sort()`, `merge()` |

`Exon`, `CDS` and `UTR` subclass `Locus`, so gene models are intervals too.

## Loci and its two indexes
{: .sec-green }

`Loci` subclasses `list`, so indexing, slicing and iteration behave as usual.
On top of that it lazily builds two indexes and keeps them until the list is
copied:

- `uids`: `{uid: row}`, so `loci["chr1:100-200(+)"]` is a dictionary lookup.
- `cgr`: an interval index over every row. It is `cgranges` (Heng Li's C
  library) when it is installed, and a pure-Python index otherwise, with
  identical results. `cgranges` is not on PyPI, so pip installs use the
  Python index.

## Set algebra keeps whole rows
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-setops.svg %}
</div><figcaption>
<strong>Operators filter rows; they do not cut intervals.</strong> <code>A &amp; B</code> keeps every interval of A that touches B (like <code>bedtools intersect -u</code>), <code>A − B</code> keeps the rest, and <code>merge()</code> fuses overlapping or book-ended intervals after sorting. A row always comes back with its own uid, so the result can index any table keyed on A.
</figcaption></figure>

| Expression | Result |
|---|---|
| `A & B`, `A.intersect(B)` | rows of A overlapping at least one row of B |
| `A - B`, `A / B`, `A.difference(B)` | rows of A overlapping nothing in B |
| `A + B`, `A | B` | concatenation (no deduplication; follow with `.sort().merge()`) |
| `A ^ B` | `(A - B) + (B - A)` |
| `A.slop(n)` | every row widened by `n` bp (start clamped at 0) |
| `A.sort()`, `A.merge()` | genome order; fused overlapping / book-ended rows |
| `A.nearest(B)` | nearest row of B per row of A (pyranges) |
| `A.map(B)` | `{a.uid: [b.uid, ...]}` for every overlap |

Every overlap test goes through B's `cgr`, so `A & B` costs one index build
for B (once) plus one lookup per row of A.

## How an overlap lookup works
{: .sec-green }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/loci-index.svg %}
</div><figcaption>
<strong>The pure-Python index scans a bounded window.</strong> Intervals are kept sorted by start, per chromosome, along with the length of the longest one. An interval that starts before <code>qs − max_len</code> ends before <code>qs</code>, so a query only tests starts in <code>[qs − max_len, qe)</code>. Two <code>bisect</code> calls find that window, so a lookup costs O(log n + k) for k candidates. Before 1.1 the scan started at the first interval, which made <code>A &amp; B</code> quadratic on pip installs.
</figcaption></figure>

The bound is tight for peak-like data, where the longest interval is a few
kilobases. A set that mixes peaks with megabase domains widens every window to
the longest domain; for such sets install `cgranges`, whose implicit interval
tree does not depend on the longest interval.

## Costs
{: .sec-green }

- **Index build:** one sort per chromosome, on first use; reused by every later
  query on the same `Loci`.
- **`A & B`:** |A| lookups into B's index; the result holds references to A's
  `Locus` objects (nothing is copied).
- **`sort()` / `merge()`:** a key sort on `(chrom, start)`, then one pass.
- **Memory:** one Python object per interval. For millions of intervals the
  [columnar Loci]({{ '/design/columnar/' | relative_url }}) stores four numpy columns instead.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/#loci' | relative_url }}) page.
