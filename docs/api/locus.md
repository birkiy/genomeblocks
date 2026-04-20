---
title: locus
parent: API Reference
layout: default
nav_order: 1
---

# `genomeblocks.locus`
{: .no_toc }

Core interval dataclasses.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Locus`

```python
@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    strand: str = "."
```

### Properties

| Name | Type | Description |
|---|---|---|
| `uid` | `str` | Canonical `"chrom:start-end(strand)"` identifier. |
| `length` | `int` | `end - start`. |
| `center` | `int` | `(start + end) // 2`. |

### Methods

```python
Locus.overlaps(other)       # bool — both on same chrom AND ranges intersect
Locus.distance_to(other)    # int  — |center - other.center|, NotImplemented if cross-chrom
Locus.sequence(genome, r=None)
    # Extract sequence from a {chrom: str} genome dict.
    # If r is None, returns genome[chrom][start:end].
    # Otherwise, returns a ±r bp window around center.
Locus.copy()                # new Locus with same fields
```

### Ordering

Locus defines `__lt__`, `__le__`, `__gt__`, `__ge__` so lists of Locus objects sort by `(chrom, start)`.

---

## `Exon(Locus)`

```python
@dataclass
class Exon(Locus):
    exon_number: int = 0
```

## `CDS(Exon)`

Alias; same fields as `Exon`, used as a type marker.

## `UTR(Exon)`

```python
@dataclass
class UTR(Exon):
    type: str = None      # "5'" or "3'"
```

All three inherit `Locus`'s UID / overlap / distance methods, so they participate in every interval API.
