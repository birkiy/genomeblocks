---
title: BEDPE & pairs
parent: Design
layout: default
nav_order: 5
---

# BEDPE & pairs
{: .no_toc }

`bedpe` reads loop files into `Pair` objects and streams Hi-C read-pair files
into window counts without loading them into memory.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## One BEDPE reader
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/bedpe-reader.svg %}
</div><figcaption>
<strong>Every loop file goes through <code>read_bedpe</code>.</strong> It returns <code>Pair</code> objects with both anchors, their midpoints and the anchor distance (∞ for inter-chromosomal loops). <code>Architecture.make</code>, the browser's arc track and <code>Loci.pair_to_bed</code> all read loops this way and never parse BEDPE themselves.
</figcaption></figure>

`read_bedpe` takes optional `min_score` and `max_distance` filters and counts
malformed lines instead of failing. `pair_to_bed(loci, bedpe, either=True)`
keeps loops with at least one anchor (or `both=True`: both anchors) on the
given loci, using the loci's interval index.

## Counting read pairs in windows
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/bedpe-count.svg %}
</div><figcaption>
<strong>One sorted key, one binary search per read end.</strong> Every window is placed on a single genome-wide key, <code>chromosome code · 2⁴⁰ + start</code>, and the keys are sorted once. The pairs file is read in chunks with categorical chromosome columns, so chromosome names are turned into codes once per category rather than once per row. Each read end is then found with one <code>searchsorted</code>, and the counts are a <code>bincount</code>.
</figcaption></figure>

- **Formats.** HiC-Pro `allValidPairs`, 4DN `.pairs` and Juicer medium are
  recognised from the first data line (`format="auto"`), and `columns=` reads
  any other layout. Gzipped files are read directly when `format=` is given
  (auto-detection peeks at plain text).
- **`count_pairs(windows, pairs)`** returns one row per window and one column per
  partner chromosome (or one `count` column for `target_chrom=`). Each read end
  is counted once, so a cis pair with both ends in windows adds two counts.
- **`count_pairs_2d(windows, pairs, loci_b=None)`** returns a sparse CSR matrix of
  window × window counts (symmetric when `loci_b` is omitted).
  `pair_2d_block` and `pair_2d_to_frame` slice or flatten it.
- **Windows** must not overlap within a chromosome (tiles from
  `Loci.tile_genome` never do). Gaps are fine; a read end in a gap is not
  counted.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/#pairs' | relative_url }}) page.
