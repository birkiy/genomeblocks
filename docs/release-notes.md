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

## v2.0.0 — columnar only

`2.0.0` makes the columnar tables introduced in 1.1 **the** API. `genomeblocks.Loci`, `Genes`, `Pairs`, `Architecture` and `Atlas` are tables of numpy columns whose row number is the join key; the 1.x object API (`Loci` as a list of `Locus`, `Gene` / `Transcript` / `Exon` trees, `list[Pair]`, `Architecture` as a `graph_tool.Graph`) is gone, and so is the `genomeblocks.columnar` subpackage — its contents are now the top-level names. There are no compatibility shims. Upgrade: `pip install -U genomeblocks`. The 1.1 documentation stays online: pick **1.1** in the version menu at the top of every page.

### What 2.0 is

- **Tables.** `Loci` (`codes` / `starts` / `ends` / `strands` + `cols`), `Genes` (three linked `Loci`: genes, transcripts, features), `Pairs` (BEDPE: anchors `a`, `b` + columns), `Architecture` (vertices = a `Loci`, an edge table with `ep` columns, `vp` vertex columns), `Atlas` (bins × tracks), plus `Locus` / `LocusView`, `GeneView` and `Genome`. Every container has `shape`, `columns`, `head`, `tail`, `describe`, `_repr_html_`.
- **One boundary, explicit exits.** `gb.as_loci(x)` turns a path (BED / narrowPeak / CSV / TSV / parquet), a region string, a `Locus`, pandas / polars / pyarrow / bioframe / PyRanges / pybedtools frames, dicts of columns, structured arrays, an AnnData or a list of regions into a `Loci`, and every public function calls it on its inputs. Out: `to_pandas` / `to_polars` / `to_arrow` / `to_bioframe` / `to_pyranges` / `to_bedtool` / `to_cgranges` / `to_anndata` / `to_records` / `to_numpy` / `to_bed` / `save` (parquet); `Architecture` adds `to_networkx` / `to_igraph` / `to_graph_tool` / `to_scipy`.
- **Protocols.** Arrow C stream, dataframe interchange and narwhals on every container, so polars, DuckDB, seaborn, plotly and altair take the objects as they are.
- **Backends.** Six families — intervals, bigwig, motifs, fasta, tables, graph — with a pip-installable default each and other engines on request (`backend=` on a call, `with gb.use_backend(...)`). `gb.backends()` lists them. A requested backend that is missing raises with the install command; there is no silent fallback. `benchmarks/bench_backends.py` times every engine; the test suite checks they agree.
- **Default install is pure pip.** `pip install genomeblocks` gives every default backend except a motif engine (`genomeblocks[motifs]` adds MOODS, `genomeblocks[lightmotif]` the alternative engine; the conda package ships MOODS), `Architecture` included (graph algorithms run on scipy; graph-tool is picked up automatically when installed). Extras: `fast`, `bam`, `viz`, `interop`, `all`.
- **New dependencies:** `narwhals`, `pyarrow`. **Dropped from core:** `pyranges` (now an optional interval backend).
- **Also new:** `Loci.from_anndata` / `to_anndata`, `Loci.liftover` (pyliftover), `Loci.sequences` / `to_seqrecords` / `to_fasta` on any FASTA backend, `Genes.from_frame` / `to_gtf` / `to_bed12` / `representative`, `Pairs.overlapping` / `filter`, `Architecture.from_edges` / `from_frame` / `from_scipy` / `pagerank` / `vertices_frame` / `draw` on any graph backend, `gb.load_motifs` for JASPAR / jaspar16 / TRANSFAC / uniprobe / MEME / Biopython / arrays, motif p-value thresholds, `scan_motifs_profile`, `gb.plot_motif_heatmap`, `gb.igv_html`, `gb.View`.

### Breaking changes

- **Coordinates.** Every table is 0-based, half-open. `Genes.make` converts GTF / GFF3 starts (`start - 1`); 1.x kept raw 1-based GTF starts. A gene's TSS is the 1-bp interval `[t, t+1)` with `t = start` on `+` and `t = end - 1` on `-`; 1.x stored `-` strand TSSs end-before-start.
- **No global state.** The session-wide default `Genome` is gone: each table carries its own, and operations between tables re-code the other side. `set_backend` is gone: use `backend=` per call or `with gb.use_backend(...)`.
- **No legacy API.** `Loci` is not a `list`; `Locus` objects are views (`L[i]`), not stored. `Gene`, `Transcript`, `Exon`, `CDS`, `UTR`, `Pair`, `make_genome`, `scan_motifs` (top level), `Architecture.to_legacy`, `Loci.cgr` and `genomeblocks.columnar` no longer exist.
- **`Architecture` methods take no `loci` argument.** `make(loci, bedpe)` fixes the vertices; `add_mcool(mcool, resolution=...)`, `normalize()`, `annotate(genes)`, `strength()`, `prime_hubs()` work on them. `make(..., trans=True)` keeps inter-chromosomal loops (they get `ep.d = inf` and are normalised against the mean trans weight); `dmax` applies to cis loops only.
- **Motif functions take `motifs` and `format=`**, not a JASPAR path: `motifs` is a path, a `Library`, Biopython motifs or matrices; `format` is `'jaspar'` (default), `'jaspar16'`, `'transfac'`, `'uniprobe'` or `'meme'`. Results are aligned to the input rows (index = uid); windows that cannot be read count 0. `scan_motifs_profile` returns `(M, names)`.
- **`signal(backend=...)`** names the bigWig engine: `'pybigtools'`, `'pybigwig'` or `'python'` (the pure-Python reader; 1.x called it `'bigwig'`).
- **The browser draws from the tables**: `browser(region, tracks)` takes bigWig paths or handles, BAM paths, `Pairs` or `.bedpe` paths, `Genes`, and anything `as_loci` takes. `architecture_draw.draw(arch, loci, region)` became `A.draw(region, layout=...)`.
- **Graph-tool is optional.** `Architecture` is not a `graph_tool.Graph`; `A.graph()` / `A.to_graph_tool()` build one on demand when graph-tool is installed.

### Migration: 1.x → 2.0

| 1.x | 2.0 |
|---|---|
| `Loci([Locus(...), ...])`, `Loci(list_of_locus)` | `Loci.from_records([...])` or `gb.as_loci([...])` (tuples, region strings and `Locus` objects all work) |
| `Loci(list(genes.annot['prom']))` | `genes.get_tss().slop(genes.promoter_r)` (already a `Loci`; `genes.annot['prom']` is the merged version) |
| `loci.cgr` (the cgranges index) | gone — the interval engine is a backend: `loci.intersect(o, backend="cgranges")` or `with gb.use_backend(intervals="cgranges")`; `loci.to_cgranges()` builds the index for your own code |
| `loci.uids` → `{uid: Locus}` | `loci.uids` → `{uid: row}`; `loci[uid]` → a `LocusView`; `loci.row(uid)` |
| `for locus in loci:` to build columns | vectorised: `loci.starts`, `loci.lengths`, `loci.centers`, `loci["score"]`; `loci["new"] = values` |
| `Pair`, `read_bedpe(path)` → `list[Pair]` | `Pairs.make(path)` (a table; `read_bedpe` returns it too); `P.a`, `P.b`, `P["score"]`, `P.distance`, `P.is_cis` |
| `pair_to_bed(pairs, loci)` | `Pairs.overlapping(loci, r=..., both=...)`; `Pairs.anchors_overlap(loci)` for the two masks |
| `make_genome(fasta)` | `gb.read_fasta(fasta)` (a `{chrom: sequence}` dict) or `Genome.from_fasta(fasta)` (names and sizes); FASTA paths / dicts / handles are accepted directly by `sequences`, `scan_motifs*`, `browser` |
| `scan_motifs(loci, fasta, jaspar_path)` | `loci.scan_motifs(fasta, motifs, format="jaspar")` — `motifs` is a path, a `Library` (`gb.load_motifs`), Biopython motifs or matrices |
| `scan_motifs_matrix(..., jaspar=...)` | `scan_motifs_matrix(loci, fasta, motifs, format=..., threshold=... or pvalue=...)`; the result index is the uid, aligned to the rows |
| `Architecture.make(cre, bedpe).add_mcool(cre, mcool, resolution=r)` | `Architecture.make(cre, bedpe).add_mcool(mcool, resolution=r)` |
| `arch.normalize(cre, source="w", name="n")` | `A.normalize(source="w", name="n")` |
| `arch.annotate(cre, genes)` | `A.annotate(genes)` |
| `arch.strength(key="n")`, `arch.prime_hubs(key="n")` | unchanged names; `A.vp.strength`, `A.vp.gene`, `A.vp.annot` are arrays aligned to `A.loci` |
| `arch.vp.uid[v]`, `arch.vertex(uid)` | vertex `i` is `A.loci` row `i`: `A.loci.uid[i]`, `A.loci.row(uid)`; `A[uid]` and `A.neighbors(uid)` still answer by uid |
| `graph_tool.all` on `arch` | `A.graph()` / `A.to_graph_tool()` (graph-tool), `A.to_networkx()`, `A.to_igraph()`, `A.to_scipy()`; `A.components()`, `A.pagerank()` on any graph backend |
| `architecture_draw.draw(arch, loci, region)` | `A.draw(region, layout="spring" | "circular" | "genomic", backend=None)` |
| `Gene` / `Transcript` / `Exon` / `CDS` / `UTR` objects, `genes.genes[name].transcripts[...]` | the three tables `genes.genes`, `genes.transcripts`, `genes.features` (0-based) + `GeneView`: `genes["TP53"].tss` / `.transcripts` / `.exons` / `.canonical`; `genes.to_pandas("transcripts")` |
| 1-based GTF starts in `Genes` | 0-based: `start - 1` is applied at parse time; `Genes.from_frame(df, one_based=True)` for a frame that still carries GTF starts |
| `-` strand TSS stored as `end` with `start > end` | `[end - 1, end)`: `genes.get_tss()`, `genes["X"].tss` |
| `genes.annotations(loci)` / `nearest_genes(loci)` | unchanged; they take anything `as_loci` takes and return one row per input row |
| `import genomeblocks.columnar as gbc; gbc.Loci` | `import genomeblocks as gb; gb.Loci` — `genomeblocks.columnar.*` → `genomeblocks.*` |
| `gbc.Architecture.make(cre, "loops.bedpe").add_mcool(...)` | unchanged shape; `Pairs` replaces the BEDPE list, `trans=True` by default |
| `set_backend("pyranges")` | `backend="pyranges"` on the call, or `with gb.use_backend(intervals="pyranges"):` |
| `signal(loci, bws, backend="bigwig")` (the pure-Python reader) | `loci.signal(bws, backend="python")`; `"pybigtools"` (default) and `"pybigwig"` are the others |
| `signal(..., bigwig=...)` with the `genomeblocks.bigwig` module | the reader moved to `genomeblocks.backends._bbi`; use the `backend=` name |
| `compare_heatmap(a, b, bigwigs)` | unchanged; `a` and `b` are anything `as_loci` takes |
| `browser(region, {"loops": read_bedpe(p)})` | `browser(region, {"loops": "loops.bedpe"})` or a `Pairs` |
| `Atlas.make(beds)`, `loci.enrich(atlas)` | unchanged; inputs go through `as_loci`, `atlas.to_pandas()` is the track table |
| `Loci.make(path)` returning unsorted rows in file order | `Loci.make` sorts into genome order (`sort=False` keeps file order); `as_loci(path)` keeps file order |

{: .tip }
> If a 1.x script breaks on `AttributeError: module 'genomeblocks' has no attribute ...`, the name is in the table above. `gb.__all__` lists the whole 2.0 surface.

### Also in 2.0

- `Loci.make` reads narrowPeak / broadPeak columns under their standard names with `keep=True` (or a list), refuses `.gtf` / `.bedpe` with a pointer to the right constructor, and reads `.gz`.
- `Genes.make(..., chr_map={'1': 'chr1'})` renames chromosomes while parsing; GFF3 and Ensembl spellings (`five_prime_utr`, `gene_biotype`) are recognised; duplicate gene ids are handled.
- `Loci.tile` / `tile_genome` take sizes from a dict, a `.chrom.sizes` path, a Series, a `Genome` or a cooler.
- `nearest` reports the gap in bases (0 for overlaps and book-ended intervals) and `(-1, -1)` when the chromosome has nothing, on every backend.
- `Loci.liftover(chain)` walks each end inward by up to `1 - min_match` of the length before dropping a row (UCSC's base-fraction rule); `cols['source_row']` says where each row came from.
- `select_isoforms` and `Genes.make(cre=, bw=)` accept any peak input (a list of peak sets is merged) and bigWig paths or handles.
- `count_pairs` / `count_pairs_2d` detect HiC-Pro, 4DN `.pairs` and Juicer layouts; `read_pairs_chunks` streams any of them.
- `se.call_se` returns the stitched regions as a `Loci` with `rank`, `signal` and `is_se` columns; `hichip.shortrange_track` writes coverage bigWigs / MACS3 input from pairs files.
- The one-file viewers (`igv_html`, `View`) draw genes from the 0-based tables and signal through the bigWig backends.

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

(In 2.0 these tables became the top-level API; `genomeblocks.columnar` no longer exists.)

### Faster, and one crash fixed (classic API)

The new `benchmarks/` suite times every block against the tools people would
otherwise use (results: [Benchmarks]({{ '/benchmarks/' | relative_url }})). Every fix below was
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
[Genes guide]({{ '/guide/genes/' | relative_url }}).

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

## Upgrade notes (1.0)

There is no prior stable release to migrate from — v1.0 is the baseline. If you
were tracking `main` before the 1.0 tag, note that several exploratory
`Architecture` methods (`make_spread`, `make_clique`, `cluster`, `focus`, …),
the `Tags` layer, and a handful of `Genes` helpers were removed in the run-up to
1.0. Use the public surface in the table above.
