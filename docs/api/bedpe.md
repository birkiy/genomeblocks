---
title: bedpe
parent: API Reference
layout: default
nav_order: 8
---

# `genomeblocks.bedpe`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Pair`

```python
@dataclass
class Pair:
    chrom1: str; start1: int; end1: int
    chrom2: str; start2: int; end2: int
    name: str = ""
    score: float = 0.0
    strand1: str = "."
    strand2: str = "."
    extra_fields: tuple = ()

    # properties
    mid1 -> int        # (start1 + end1) // 2
    mid2 -> int        # (start2 + end2) // 2
    distance -> int    # |mid2 - mid1| for intra-chrom; float('inf') otherwise
```

---

## `read_bedpe(filename, *, min_score=None, max_distance=None, verbose=True) -> list[Pair]`

Parser. Filters applied in-read if provided.

---

## `pair_to_bed(loci, bedpe, *, r=0, either=True, both=False, min_score=None, max_distance=None, verbose=True) -> list[Pair]`

Equivalent to `bedtools pairtoBed`. Also attached as `Loci.pair_to_bed(...)`.

Arguments:

| Arg | Meaning |
|---|---|
| `loci` | `Loci` (or single `Locus`) to intersect against. |
| `bedpe` | Path to a BEDPE file or a pre-loaded `list[Pair]`. |
| `r` | Slop radius per anchor before the overlap test. |
| `either` | Keep pair if ≥ 1 anchor overlaps (default `True`). |
| `both` | Require both anchors overlap (default `False`). |
| `min_score`, `max_distance` | In-read filters (passed to `read_bedpe`). |

---

## `pairs_to_frame(pairs) -> pandas.DataFrame`

Columns: `chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`.

---

## `pairs_to_bedpe(pairs, filename) -> None`

Writes TSV BEDPE (10 or more columns when `extra_fields` is set).
