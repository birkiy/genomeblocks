---
title: Release Notes
layout: default
nav_order: 80
permalink: /release-notes/
---

# Release Notes
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Unreleased — ATAC-supported isoforms

`Genes.make()` / `Genes.make_ucsc()` take `cre` (peaks: a path, a `Loci`, or a list)
and/or `bw` (bigwigs), plus `r` (TSS half-window, defaults to `promoter_r`) and `kw`;
the new `Genes.select_isoforms()` does the same job on an already-parsed object. Each gene is pointed at the longest isoform whose TSS is open
in your data instead of its longest annotated isoform, which keeps long silent
isoforms (TGFBR3 and friends) from dragging a gene's body and TSS across the locus.
Transcripts gain `tss_score` / `tss_support` / `tss`, genes gain `canonical` /
`canonical_transcript` / `set_span()`, and no transcript is dropped. See the
[Genes guide](guide/genes#picking-the-isoform-your-cells-actually-use).

## v1.0.1 — pip-install fix

`1.0.0` listed `cgranges` (Heng Li's C interval index) as a dependency, but
**`cgranges` is not published on PyPI**, so `pip install genomeblocks` failed
with `No matching distribution found for cgranges`.

- `cgranges` is no longer a pip dependency. When it isn't importable,
  genomeblocks now uses a **pure-Python interval index** with identical overlap
  results, so `pip install genomeblocks` works out of the box. Install
  `cgranges` (conda, or `pip install git+https://github.com/lh3/cgranges`) for
  the fast C path on large sets.

If you installed `1.0.0` from PyPI and hit the error, upgrade:
`pip install -U genomeblocks`.

---

## v1.0.0 — first stable release

`genomeblocks` v1.0 is the first release with a frozen public surface, a full
synthetic test suite, and complete module documentation. The public API is the
set of names re-exported from the top-level package plus the `Atlas` class:

```python
from genomeblocks import (
    Locus, Exon, CDS, UTR,          # single intervals
    Loci,                            # interval container + set algebra
    Gene, Transcript, Genes,         # gene models / GTF parsing
    Architecture,                    # chromatin-contact graph
    Atlas,                           # GIGGLE-style enrichment
    signal, tmm, compare_heatmap,    # bigWig signal
    scan_motifs, make_genome,        # motif scanning
    browser,                         # region viewer
    coverage,                        # BAM per-base coverage
)
```

### Public surface

The building blocks and their entry points:

| Block | Build from | Key methods |
|---|---|---|
| `Loci` | `.make(bed)`, `.tile`, `.from_frame` | `& \| - ^`, `slop`, `sort`, `merge`, `nearest`, `signal`, `enrich`, `liftover` |
| `Genes` | `.make(gtf)`, `.make_ucsc` | `annotations`, `nearest_genes` |
| `Architecture` | `.make(bedpe)` | `add_mcool`, `normalize`, `annotate`, `strength`, `prime_hubs` |
| `Atlas` | `.make(beds)` | `search`, `bootstrap`, `save`/`load` |
| `signal` | `loci.signal(bigwigs)` | `tmm`, `plot_heatmap`, `compare_heatmap` |
| `browser` | `browser(region, tracks)` | mixed bigWig / BAM / bed / bedpe / gene tracks |

### Correctness fixes since the pre-1.0 code

- **`Locus` ordering across chromosomes.** `>` and `>=` returned the *less-than*
  result across chromosomes, so `a > b` and `b > a` could disagree with `<`.
  Cross-chromosome comparisons are now consistent. Within-chromosome ordering
  (and therefore `Loci.sort` / `merge`) was already correct.
- **`Locus` is now hashable** (by its UID), so a `Locus` can be used in a `set`
  or as a `dict` key — matching its value-like `__eq__`.
- **`Architecture.copy()` / `subgraph()` preserve property types.** Types were
  inferred from a sampled Python value, but graph-tool returns booleans as
  ints, so `bool` vertex/edge properties were silently copied as `double`.
  Property maps now keep their exact graph-tool value type (and empty graphs
  keep their properties too).
- **`Genes.make` no longer crashes on a generic `UTR` feature without a CDS.**
  A non-coding transcript (or a `UTR` line preceding its `CDS`) used to raise
  `IndexError`; the 5′/3′ inference now falls back safely.
- **Signal memory guard checks RAM, not disk.** The `signal()` size check now
  compares the requested cube against available physical memory (it previously
  measured free space on `/`), and raises a clearer `MemoryError`.

### `browser` import is now unambiguous

The region-viewer drawing code moved from `genomeblocks/browser.py` to
`genomeblocks/browserview.py`, so the public name `browser` is **always the
callable** — the old `'module' object is not callable` footgun (where a stray
`import genomeblocks.browser` shadowed the function) is gone:

```python
from genomeblocks import browser              # the function
from genomeblocks.browserview import browser  # equivalent
```

### Testing

The suite is fully synthetic — no external genome / bigWig / ChIP-Atlas data —
and covers both unit helpers and end-to-end pipelines (a real bigWig written
in-test and read back through both signal backends, motif scanning on a planted
PWM, `Atlas` enrichment, a `Loci → Genes → Architecture` pipeline, and a
headless `browser()` render):

```bash
pip install -e ".[test]"
pytest
```

`graph-tool` must be importable for the `Architecture` tests (it is the one
conda-only dependency).

---

## Upgrade notes

There is no prior stable release to migrate from — v1.0 is the baseline. If you
were tracking `main` before the 1.0 tag, note that several exploratory
`Architecture` methods (`make_spread`, `make_clique`, `cluster`, `focus`, …),
the `Tags` layer, and a handful of `Genes` helpers were removed in the run-up to
1.0. Use the public surface in the table above.
