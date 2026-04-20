---
title: loci
parent: API Reference
layout: default
nav_order: 2
---

# `genomeblocks.loci`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Loci(list)`

List subclass of `Locus` intervals with interval algebra and cached indexes.

### Construction

```python
Loci(iterable=(), *, filename=None)         # bare constructor
Loci.make(filename, filetype=None)          # parse a BED / narrowPeak file
```

File types auto-detected from extension:

- `.narrowPeak` → BED6+N (columns 1..3 required; col 6 used for strand if present).
- anything else → BED-ish (tab-separated; chrom, start, end, optional strand col 6).

### Properties

| Name | Description |
|---|---|
| `cgr` | Lazy `cgranges` interval tree. |
| `uids` | Lazy `{uid → index}` dict. |

### Set-algebra operators

| Operator | Method | Semantics |
|---|---|---|
| `a & b` | `a.intersect(b)` | Loci in `a` overlapping any `b`. |
| `a - b` | `a.difference(b)` | Loci in `a` **not** overlapping any `b`. |
| `a / b` | `a.difference(b)` | Alias for `-`. |
| `a + b` | concat | Like `Loci([*a, *b])`. |
| `a \| b` | concat | Same as `+`. |
| `a ^ b` | `(a-b) + (b-a)` | Symmetric difference. |

### Methods

```python
Loci.overlaps(key, start=None, end=None)     # key is chrom str or Locus
Loci.slop(n)                                  # expand ±n bp
Loci.sort()                                   # by (chrom, start)
Loci.merge()                                  # collapse overlapping / bookended
Loci.nearest(other, s_names=None, o_names=None)  # → pandas.DataFrame (pyranges)
Loci.subloci(uids: list[str])                 # subset by UID list
Loci.map(other)                                # {a_uid: [b_uid, ...]}
Loci.tag(other, tag_name: str)                 # → Tags
Loci.copy()                                    # shallow copy preserving filename
```

### Conversion

```python
Loci.to_frame(names=None)     # pandas.DataFrame (Chr, Start, End, Strand, Name)
Loci.to_pyranges(names=None)  # pyranges.PyRanges
Loci.to_bed(path=None)        # write BED6 ('\t' separated) to path
```

### Attached methods (from other modules)

| Source module | Method |
|---|---|
| `signal` | `Loci.signal(bigwigs, ...)`, `Loci.plot_heatmap(...)`, `Loci.plot_profiles(...)` |
| `motifs` | `Loci.scan_motifs(genome, motif_path, ...)` |
| `bedpe` | `Loci.pair_to_bed(bedpe, ...)` |

See the individual module pages for signatures.

### Indexing

```python
loci[0]                  # int → Locus
loci[1:5]                # slice → Loci
loci["chr1:100-500(.)"]  # str  → Locus (UID lookup, O(1))
```

### Pickling

`__getstate__` / `__setstate__` drop the cached `_cgr` index so pickles don't leak stale trees.
