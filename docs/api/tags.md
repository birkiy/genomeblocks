---
title: tags
parent: API Reference
layout: default
nav_order: 3
---

# `genomeblocks.tags`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Tags`

```python
Tags(loci: Loci, *, verbose=True)                     # constructor
Tags.make(loci: Loci, *, verbose=True) -> Tags        # factory (prefer this)
```

### `.add(mapping) -> Tags`

Register one or more tags. Each value in the mapping is:

| Payload | Storage kind |
|---|---|
| `Loci` | set tag (of UIDs) |
| iterable of UIDs | set tag |
| `{uid: value}` | numeric tag |

Returns `self` so you can chain.

### `.query(predicate) -> Loci`

`predicate` is a callable `l -> TagSelection` where `l` is a tag accessor. Returns a fresh `Loci` of matching members.

### Introspection

```python
Tags.keys()         # iterable of tag names
Tags.items()        # (name, data) pairs
Tags.uids()         # full UID list in insertion order
Tags[name]          # TagSelection or TagNumericView
name in tags        # membership test (alias-resolved)
Tags.table()        # plain-text summary
```

---

## `TagAccessor`

Returned by `tags.query`'s lambda argument. Supports `l.name` and `l["name with space"]` — case-insensitive alias resolution.

---

## `TagSelection`

Result of a set-tag lookup or a fully-resolved query. Behaves like a set of UIDs.

Operators:

| Operator | Meaning |
|---|---|
| `s & t` | intersection |
| `s \| t` | union |
| `s - t` | difference |
| `s ^ t` | symmetric difference |

```python
TagSelection.uids        # set of UIDs
TagSelection.to_loci()   # materialize back to a Loci
len(sel), iter(sel)
```

---

## `TagNumericView`

Result of a numeric-tag lookup. Supports comparison operators against a scalar:

```python
l.cpm > 3.0
l.cpm >= 3.0
l.cpm < 3.0
l.cpm <= 3.0
l.cpm == 3.0
l.cpm != 3.0
```

Each returns a `TagSelection`. Combine with `&` / `|` against other selections.

---

## `PendingNumericCombination`

Internal bookkeeping type returned when you `&` a `TagSelection` and a `TagNumericView` before applying a comparator:

```python
(l.atac & l.cpm) > 3.0   # TagSelection & TagNumericView → PendingNumericCombination
                         # > 3.0 resolves → TagSelection
```

You rarely need to touch it directly — it exists so `l.atac & l.cpm > 3.0` parses cleanly.
