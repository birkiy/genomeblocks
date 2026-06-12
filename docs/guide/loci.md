---
title: Loci
parent: User Guide
layout: default
nav_order: 1
---

# Loci
{: .no_toc }

`Loci` is the container at the heart of `genomeblocks`: a list of `Locus` intervals with fast overlap queries, set algebra, and a dozen convenience methods.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Construction

```python
from genomeblocks import Loci, Locus

# From a BED / narrowPeak file (tab-separated; chrom, start, end, optional strand in col 6)
loci = Loci.make("peaks.bed")

# From a list of Locus objects
loci = Loci([Locus("chr1", 100, 500),
             Locus("chr1", 800, 1200, "+")])
```

`Loci` is a `list` subclass — iteration, `len()`, slicing, `append` all behave like a plain list. Each element is a `Locus` (dataclass) with `.chrom`, `.start`, `.end`, `.strand`, plus computed `.uid`, `.length`, `.center`.

{: .note }
> `Loci.make()` does not support GTF/GFF directly — those should go through `Genes.make()` which returns a `Genes` container. For BED-like files (`.bed`, `.narrowPeak`, `.broadPeak`, custom) `Loci.make()` is the right entry point.

---

## Indexing

`Loci` supports integer, slice, and UID-based indexing:

```python
loci[0]                          # first Locus
loci[:10]                        # first 10 as a Loci
loci["chr1:100-500(.)"]          # lookup by UID (O(1) via cached dict)
```

UID lookup uses the `Loci.uids` property (lazily populated `{uid → index}` dict). If you mutate the list after building it, invalidate by assigning `loci._uids = None`.

---

## Set algebra

All operators work on *overlap* semantics: a locus in `a` is kept iff it overlaps at least one locus in `b`.

```python
a & b    # intersection — loci from a that overlap any b
a | b    # union (concatenation, duplicates preserved)
a + b    # same as |
a - b    # difference — loci from a that do NOT overlap any b
a / b    # same as -
a ^ b    # symmetric difference: (a-b) + (b-a)
```

The named methods mirror the operators and are the preferred form in chains:

```python
specific = cre.difference(super_enhancers)
shared   = cre.intersect(super_enhancers)
```

---

## Interval transforms

```python
loci.slop(500)    # extend each interval by ±500 bp
loci.sort()       # chrom, start; returns a new Loci
loci.merge()      # collapse overlapping / bookended intervals
```

Chaining is idiomatic:

```python
cre = Loci.make("peaks.bed").slop(100).sort().merge()
```

{: .warning }
> `merge()` operates on sorted input. Call `.sort()` first (or chain `.sort().merge()`). `merge()` itself sorts internally but returns a fresh `Loci`.

---

## Fast overlap queries

Under the hood, `Loci` builds a [cgranges](https://github.com/lh3/cgranges) tree (Heng Li) the first time you run an overlap — see [Credits](../credits#cgranges) for the citation:

```python
# Query by chrom / start / end
hits = loci.overlaps("chr1", 100_000, 200_000)

# Query by another Locus
query = Locus("chr1", 150_000, 160_000)
hits = loci.overlaps(query)
```

Returns a new `Loci` of matching intervals. All of `intersect`, `difference`, `map`, and `Architecture.make` use this same index.

---

## Nearest-neighbor lookup

```python
# Uses pyranges under the hood (Stovner & Sætrom 2020 — see Credits)
df = loci.nearest(other_loci)
# → DataFrame with columns Chr, Start, End, ... and Name_b (neighbor UID)
```

Optionally pass `s_names` / `o_names` to label the output columns with anything you like (e.g. gene names).

---

## Conversion

```python
loci.to_frame()        # pandas DataFrame: Chr, Start, End, Strand, Name (=uid)
loci.to_pyranges()     # pyranges.PyRanges (for joins / overlap reports)
loci.to_bed("out.bed") # BED6: chrom, start, end, name, 0, strand
```

---

## Mapping one set against another

```python
m = a.map(b)
# → {a_uid: [b_uid, b_uid, ...], ...}
```

Useful when you need to preserve *all* overlaps (not just "is there any"). For one-to-one nearest-feature assignment, prefer `nearest()`.

---

## Signal extraction & plotting

`Loci` carries signal methods contributed by the `signal` / `signal_draw`
modules:

```python
# threaded bigWig extraction
cube = loci.signal(["a.bw", "b.bw"],
                   n_bins=200, flank=3_000)   # shape (n_loci, 2, 200)

# matching visualizations — grouping is a plain dict[str, Loci]
loci.plot_heatmap(cube, groups={"up": up_loci, "down": down_loci})
loci.plot_profiles(cube, groups={"up": up_loci, "down": down_loci})
```

See the [Signal](signal) guide and the [AR & FOXA1 walkthrough](../walkthrough/) for the full story.

---

## Example: building a CRE catalogue

```python
from genomeblocks import Loci

atac = Loci.make("atac.narrowPeak")
h3k4me1 = Loci.make("H3K4me1.narrowPeak")
h3k27ac = Loci.make("H3K27ac.narrowPeak")

# CRE = ATAC peak extended by 100bp, merged with any active mark
cre = (atac
       .slop(100)
       .sort()
       .merge())

# Tier by mark presence
typical    = cre & h3k27ac            # active
primed     = (cre & h3k4me1) - h3k27ac  # primed but not active
inaccessible = atac - cre

cre.to_bed("cre.bed")
```
