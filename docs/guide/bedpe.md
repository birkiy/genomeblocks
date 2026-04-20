---
title: BEDPE
parent: User Guide
layout: default
nav_order: 7
---

# BEDPE
{: .no_toc }

Lightweight paired-interval support — reading, filtering, and intersecting loop files with `Loci`.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The `Pair` dataclass

```python
from genomeblocks.bedpe import Pair, read_bedpe

p = Pair(
    chrom1="chr1", start1=100, end1=200,
    chrom2="chr1", start2=900, end2=1000,
    name="loop1", score=4.2, strand1=".", strand2=".",
)

p.mid1        # 150
p.mid2        # 950
p.distance    # 800  (∞ for inter-chromosomal pairs)
```

Each `Pair` is a dataclass with both standard BEDPE columns and a free-form `extra_fields` tuple for anything beyond column 10.

---

## Reading BEDPE files

```python
pairs = read_bedpe(
    "loops.bedpe",
    min_score=3.0,       # drop weak loops
    max_distance=2e6,    # drop long-range (e.g. keep cis within 2 Mb)
    verbose=True,
)
```

Malformed / short lines are silently skipped and the count reported at the end.

---

## Intersecting with Loci (like `bedtools pairtoBed`)

```python
hits = loci.pair_to_bed(
    bedpe="loops.bedpe",   # or a pre-loaded list of Pair
    r=1000,                # slop anchors by ±1 kb before the overlap test
    either=True,           # keep pairs where ≥1 anchor overlaps
    both=False,            # if True, require both anchors to overlap
    min_score=3.0,
    max_distance=1e6,
)
# → list[Pair]
```

Use `both=True, either=False` to restrict to pairs where both ends fall inside the input loci (e.g. CRE–CRE loops only).

---

## Export back to disk

```python
from genomeblocks.bedpe import pairs_to_bedpe, pairs_to_frame

pairs_to_bedpe(hits, "filtered_loops.bedpe")   # tab-separated output
df = pairs_to_frame(hits)                       # pandas DataFrame
```

---

## Integration points

- **`Architecture.make(loci, bedpe, r=...)`** consumes a BEDPE file path directly — no need to pre-parse.
- **`Architecture.make_spread(source_loci, bedpe, ...)`** same.
- **`browser(..., tracks={"loops": list[Pair]})`** uses `read_bedpe` internally when you pass a string path.

---

## End-to-end: enhancer–promoter loops only

```python
from genomeblocks import Loci, Genes

genes = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)
promoter = Loci(list(genes.annot["prom"]))
enhancer = Loci.make("enhancers.bed")

# Keep only loops where one anchor is a promoter and the other an enhancer
promoter_loops = promoter.pair_to_bed("loops.bedpe", r=2500, both=False, either=True)
ep_loops = [p for p in promoter_loops
            if any(l.chrom == p.chrom1 and l.start <= p.mid1 <= l.end for l in enhancer)
            or any(l.chrom == p.chrom2 and l.start <= p.mid2 <= l.end for l in enhancer)]
```

(For the real pipeline use `Architecture.make_spread(promoter, bedpe, hops=1)` — same result, faster.)
