---
title: Tags
parent: User Guide
layout: default
nav_order: 2
---

# Tags
{: .no_toc }

`Tags` lets you attach arbitrary labels and numeric values to a `Loci` set, then query them with readable boolean / numeric expressions.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Why Tags?

Writing `cre.intersect(h3k27ac).intersect(atac).difference(promoters)` for every combination gets old. `Tags` hoists the membership into named sets and lets you express the combination directly:

```python
from genomeblocks import Tags

tags = (Tags.make(cre)
            .add({"atac":   atac_peaks,
                  "h3k27ac": h3k27ac_peaks,
                  "promoter": promoter_cre}))

active_enhancers = tags.query(lambda l: (l.atac & l.h3k27ac) - l.promoter)
```

`active_enhancers` is a `Loci` containing exactly those CREs whose UIDs satisfy the expression.

---

## Building a Tags object

```python
tags = Tags.make(loci, verbose=True)
```

`tags` is anchored to `loci` — it remembers the UID order so queries always return `Loci` in that canonical order.

Add sets and mappings:

```python
tags.add({
    "atac":    atac_peaks,                     # Loci → set-tag
    "uids":    ["chr1:100-500(.)", ...],       # iterable of UIDs → set-tag
    "h3k27ac_cpm": {l.uid: float(s[i])         # UID → value mapping → numeric tag
                    for i, l in enumerate(loci)},
})
```

Three payload shapes are supported:

| Payload | Stored as | Used in queries via |
|---|---|---|
| `Loci` / UID iterable | `TagSelection` (set of UIDs) | boolean ops (`&`, `|`, `-`, `^`) |
| `{uid: value}` mapping | `TagNumericView` | comparison ops (`<`, `<=`, `==`, `>`) |

---

## Querying

All queries are lambdas that take a tag accessor and return a selection:

```python
# Set-only
sel = tags.query(lambda l: l.atac & l.h3k27ac)

# Numeric
strong = tags.query(lambda l: l.h3k27ac_cpm > 3.0)

# Mixed — intersect a set with a numeric threshold
active_strong = tags.query(lambda l: (l.atac & l.h3k27ac) & (l.h3k27ac_cpm > 3.0))
```

The accessor resolves case-insensitive aliases and auto-split tokens, so all of these work for a tag named `"H3K27ac CPM"`:

```python
l.h3k27ac_cpm
l["H3K27ac CPM"]
l.h3k27accpm
```

Ambiguous aliases raise `AttributeError` — use the full name.

---

## Inspection

```python
print(tags.table())
# Name               Type       Count
# ---------------------------------------
# atac               set        73412
# h3k27ac_cpm        mapping    73412
# promoter           set        12804
```

`tags.keys()`, `tags.items()`, and `tags.uids()` behave like their dict counterparts.

---

## Coercion back to Loci

Every `TagSelection` can become a `Loci` directly:

```python
sel = tags["atac"]
active = sel.to_loci()
```

This is what `query(...)` does internally.

---

## Interop with heatmaps & profiles

`signal.plot_heatmap` / `signal.plot_profiles` accept a `Tags` object for row grouping:

```python
tags = Tags.make(union)
tags.add({"A-specific": a - b, "shared": a & b, "B-specific": b - a})

loci.plot_heatmap(cube, tags=tags, sets=["A-specific", "shared", "B-specific"])
```

This is exactly what `compare_heatmap` does under the hood — see the [Signal guide](signal).

---

## Under the hood

- Set tags are stored as `frozenset[uid]` — queries use native set ops.
- Numeric tags are `dict[uid, value]` — comparisons return a `TagSelection`.
- Mixed expressions like `l.atac & (l.cpm > 2)` go through a `PendingNumericCombination` trampoline that defers the comparison until a comparator is applied.
- All UIDs are validated on insert: registering a tag with an unknown UID raises `KeyError` immediately.
