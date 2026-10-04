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

## v1.1.0 — columnar foundation

`1.1.0` adds `genomeblocks.columnar`, a table-first foundation that sits next to
the classic API (which is unchanged). It also adds a benchmark suite and fixes
every slow path that suite found. Upgrade: `pip install -U genomeblocks`
(`pip install -U "genomeblocks[columnar]"` adds polars and pyarrow).

### `genomeblocks.columnar` (new)

CREs, genes and the Architecture are numpy tables that share one `Genome`, and
the row number is the join key: row `i` of the CREs is row `i` of the labels,
the signal cube and every graph column.

```python
import genomeblocks.columnar as gbc
cre   = gbc.Loci.make("atac.narrowPeak")
genes = gbc.Genes.make("gencode.gtf")              # GTF or GFF3
A = (gbc.Architecture.make(cre, "loops.bedpe")
       .add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes).strength())
A.chrom("chr8"); A.cis; A.trans                    # zero-copy views
```

- `Loci`: genome-sorted columns, vectorised overlap / merge / nearest, and
  hand-off to pandas, polars, Arrow and parquet. Classic functions such as
  `signal()` run on it unchanged.
- `Genes`: genes ⇄ transcripts ⇄ features linked by row numbers, GTF and GFF3
  (GENCODE and Ensembl style), and a vectorised `select_isoforms`.
- `Architecture`: one edge table sorted into per-chromosome cis blocks plus a
  trans block. Trans loops are kept and normalised against the mean trans
  weight. Includes `near` / `support`, graph-tool on demand, and parquet
  save / load.
- `hichip` (short-range HiChIP ends → coverage bigWig / MACS3 input), `se`
  (ROSE-style super-enhancers), and `view` / `igv` (one-file shareable browsers).
- Measured against 1.0.1, the whole pipeline in a fresh process takes 43 s →
  2.7 s with identical edges, weights, O/E, labels and hubs, and peak memory
  drops from 1.26 GB to 0.78 GB.
- polars and pyarrow are optional: without them parsing falls back to pandas
  with identical results (parquet save / load needs pyarrow).

See `examples/columnar_prototype` and `examples/case_study_se` in the repository.

### Faster, and one crash fixed (classic API)

The new `benchmarks/` suite times every block against the tools people would
otherwise use (report: `benchmarks/report/index.html`). Every fix below was
checked to give output identical to 1.0.1.

- **`Architecture.normalize` no longer crashes on inter-chromosomal edges.**
  Trans edges get `ep.d = inf` and an expectation equal to the mean trans
  weight. The power law is fitted on cis edges only, so cis-only graphs give
  bit-identical O/E, and `prune()` keeps trans edges.
- The pure-Python overlap index (used whenever `cgranges` is missing, i.e. on
  pip installs) is now O(log n + k) per lookup: `A & B` on 100k × 100k peaks
  takes 24.5 s → 0.24 s.
- `add_mcool` 45×, `strength` 150× and `normalize` 14× faster (array work
  instead of per-edge loops).
- `count_pairs` / `count_pairs_2d` are 3–4× faster, `scan_motifs_matrix`
  3.5×, and `Loci.sort` / `merge` 2–3×. The pure-Python bigWig reader is 1.7×
  faster on sorted peaks.
- `from genomeblocks import Loci` no longer imports matplotlib, pandas or
  scipy: ~1.1 s → ~0.2 s.
- `signal()` uses `fillna=` on pybigtools ≥ 0.3, so there is no deprecation
  warning on every call. `Loci[np.int64(i)]` works.
- **Behaviour change:** `signal(workers=n)` now honours `n` up to the core
  count; it used to be silently capped at half the cores. `workers=None`
  asks for half the cores.

### ATAC-supported isoforms

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
