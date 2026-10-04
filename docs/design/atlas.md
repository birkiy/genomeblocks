---
title: Atlas
parent: Design
layout: default
nav_order: 6
---

# Atlas
{: .no_toc }

`Atlas` answers "which of these thousands of ChIP-seq tracks overlap my regions
more than expected?" with one sparse matrix and one vectorised test.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## One bit per bin and track
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/atlas-index.svg %}
</div><figcaption>
<strong>A query is a column sum.</strong> The genome is tiled into fixed bins, and every track becomes one column of a sparse bins × tracks matrix with a 1 wherever the track has a peak. A query's peaks are turned into the set of bins they touch, and <code>M[query_bins].sum(axis=0)</code> gives the overlap count for every track in one sparse operation. That fills the 2×2 table of a Fisher test for all tracks at once.
</figcaption></figure>

## Building the index
{: .sec-navy }

`Atlas.make(beds, chromsizes, bin_size=1000)` tiles each chromosome at
`bin_size` and lays the chromosomes end to end with fixed offsets. It reads
each BED file into bin ids (in parallel worker processes) and assembles one
CSC matrix, which it stores as CSR for fast row slicing. Metadata (cell type,
antibody, ...) can be attached from a table and travels with the results.
`save()` and `load()` round-trip the whole index through one `.npz` file.

## Searching
{: .sec-navy }

| Mode | Cells of the 2×2 table |
|---|---|
| `search(query)` | a = query bins in the track, b = rest of the query, c = rest of the track, d = rest of the genome |
| `search(query, ref=other)` | a / b from the query, c / d from the reference set: "more in query than in ref" |
| `bootstrap(query, n)` | an empirical null from `n` position shuffles of the query (each interval stays on its chromosome by default) |

The result is a DataFrame per track with overlaps, log2 odds, p-value and a
GIGGLE-style score (sorted by score), joined with the track metadata.
`Loci.enrich(atlas)` and `Loci.enrich_mc(atlas, n)` call `search` and
`bootstrap` on a `Loci` directly.

Working in bins rather than base pairs is what makes a query one sparse
product. The cost is resolution: two peaks in the same bin count as one bin.
Choose `bin_size` near your peak width.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/#atlas' | relative_url }}) page.
