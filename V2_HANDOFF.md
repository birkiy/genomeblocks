# genomeblocks 2.0 — handoff notes (branch `v2-columnar`)

Work in progress. This file says what 2.0 is, what is done, what is left, and
how to pick it up. **Delete it before the 2.0 release.** Rollback point: tag
`v1.1.0`.

Last updated: 2026-10-05 (second checkpoint: usage limit reached).

## What 2.0 is (the owner's brief)

1. **Columnar only.** One flat package: `genomeblocks.Loci / Genes / Architecture
   / Pairs / Atlas` *are* the columnar tables (numpy columns, row number = join
   key). No classic object API, no `to_legacy()`, no `genomeblocks.columnar`
   subpackage, no compatibility shims ("there are not really users, no need to
   confuse"). Breaking changes are fine.
2. **Interoperable with everything** — pandas, polars, pyarrow, bioframe,
   pyranges, pybedtools, cgranges, AnnData, Biopython, MOODS, pyBigWig,
   pybigtools, networkx, igraph, ... — in and out.
3. **Backends the user can pick**, defaulting to our columnar code / the
   fastest pip-installable engine, with **benchmarks of every backend inside
   genomeblocks**. Graph default: graph-tool (owner's wish) when installed,
   else scipy (keeps the pip default path working; results are identical).
4. **Our own pure-Python implementations must be the fastest possible.**
5. **Documentation site rewritten from scratch** for 2.0.

### The governing guideline (owner's text, verbatim)

> **Core principle** — Compatibility is explicit and visible. Backend choice is
> automatic and invisible. Build the boundary first; add compute backends only
> when a real dataset forces it.
>
> **1. Inputs** — Every public function/constructor accepts the field's standard
> objects: pandas DataFrame, polars DataFrame/LazyFrame, pyarrow Table, numpy
> arrays, and the common file formats (CSV/TSV, BED-like, parquet). Normalize at
> the boundary once (e.g. via Narwhals or `__dataframe__`), then work on one
> internal representation. Never branch on input type deep inside the code.
> Column-name conventions are lenient on input (`chrom`/`chr`/`seqname`),
> strict internally.
>
> **2. Outputs** — Explicit escape hatches on every container: `.to_pandas()`,
> `.to_polars()`, `.to_arrow()`, `.to_numpy()`, plus domain ones
> (`.to_networkx()`, `.to_bed()`, …). Zero-copy where possible; document when
> they copy. Implement `__dataframe__`, `__arrow_c_stream__`, `__array__`,
> `__len__`, `__iter__`, `__getitem__`. Round-trip test: `from_X(obj).to_X()`
> equals `obj` for every supported X.
>
> **3. Notebook ergonomics** — short `__repr__`; `_repr_html_` head/summary;
> `.head()`, `.describe()`/`.summary()`, `.columns`, `.shape`. No hidden global
> state; no required setup before first use. Errors name the expected input and
> the fix. Import time < 1 s; heavy optional deps imported lazily.
>
> **4. Compute backends (only if needed)** — One default backend, chosen for
> installability and correctness; default must be pure-pip installable. Extra
> backends only for a documented, measured need / a user who asked. Selection is
> automatic from the input where safe; `backend=` is an override, never a
> requirement. No silent fallback: a requested backend that is missing raises
> with the install command. Backends implement a small internal interface;
> public API never references a backend object; features not supported by all
> backends are rejected with a clear error. Test matrix runs every public
> function against every backend; results equal within documented tolerance.
>
> **5. Anti-patterns** — our object as the only accepted input; returning only
> our object with no `.to_pandas()`; conda-only/compiled deps on the default
> path; backend differences leaking into results (nulls, row order, dtypes);
> building the backend layer before the boundary layer.
>
> **Acceptance checklist** — pandas and polars frames accepted everywhere;
> `.to_pandas()`/`.to_polars()` on every container; `__dataframe__` / Arrow
> export so seaborn/plotly/altair accept our objects directly; `_repr_html_`
> renders in Jupyter and VS Code; round-trip tests for every format; default
> install is `pip install genomeblocks`; with >1 backend: parity tests,
> explicit `backend=`, informative missing-backend error.

## Decisions already taken

- **Coordinates:** every table is 0-based half-open. `Genes` converts GTF/GFF3
  starts (`start - 1`); UCSC genePred is already 0-based. A gene's TSS is the
  1-bp interval `[t, t+1)` on its strand, `t = start` (+) or `end - 1` (−).
  (1.x kept raw 1-based GTF starts and stored '−' TSS end-before-start.)
- **No global state:** the session-wide default `Genome` is gone. Each table
  gets its own `Genome` unless `genome=` is passed; operations between tables
  re-code the other side through a small lookup (`Loci._check`). There is no
  global `set_backend`; use `backend=` per call or `with gb.use_backend(...)`.
- **Backends** (`genomeblocks/backends/`): families intervals, bigwig, motifs,
  fasta, tables, graph. Automatic choice only where answers are identical
  (`AUTO` in `backends/__init__.py`). Unknown/missing requested backend raises.
  `gb.backends()` (the module is callable) lists them.
- **Boundary:** `interop.as_loci()` / `interop.frame()` (narwhals) normalise every
  input once; `Loci.from_frame` finds columns by name, falls back to position
  only for a text+int+int header-less frame, otherwise raises naming the
  accepted spellings.
- **Protocols:** `_table.TableMixin` gives `__arrow_c_stream__`,
  `__dataframe__`, `__narwhals_dataframe__` (makes altair work — altair 6.3
  cannot consume interchange-only objects), `to_polars`, `shape`.
- `make_genome` → `read_fasta` (+ `Genome.from_fasta`); motif functions take
  `motifs` (path / `Library` / Biopython motifs / matrices) and `format=`.
- Motif results are aligned to the input rows (index = uid); windows that
  cannot be read count 0. `scan_motifs_profile` returns `(M, names)`.
- `Pairs` (BEDPE table, `bedpe.py`) replaces `list[Pair]`; `bedpe.py` stays the
  one BEDPE reader. `Architecture.make(loci, bedpe, trans=True)`; `dmax` now
  applies to cis loops only.
- New core deps planned: `narwhals`, `pyarrow`. Drop `pyranges` from core
  (optional backend only).

## Status

### Done (written; quick checks passed, not yet covered by the test suite)

| file | notes |
|---|---|
| `backends/__init__.py` | registry, `resolve`, `use_backend`, `backends()` |
| `backends/intervals.py` | genomeblocks / cgranges / ncls / bioframe / pyranges / bedtools; parity verified on random data (pairs, any, nearest, merge, point lookups) |
| `backends/bigwig.py` | pybigtools / pyBigWig / pure-Python handles; accepts open handles |
| `backends/_bbi.py` | the pure-Python bigWig reader (moved from `bigwig.py`) |
| `backends/fasta.py` | our `.fai` reader (default) / pysam / pyfaidx / memory / biopython; parity verified |
| `backends/motifs.py` | `Library`, `load_motifs` (jaspar, jaspar16, transfac, uniprobe, meme, Bio.motifs, arrays), `Block` scanner for lightmotif / MOODS / Biopython — identical hits on probes |
| `backends/graph.py` | native graphs, components (canonical labels), pagerank, layouts |
| `backends/tables.py` | polars / pandas TSV parsing |
| `genome.py`, `locus.py`, `_intervals.py` | `Genome.from_fasta/from_sizes/set_sizes`, `Locus.parse/from_uid`, nearest kernel + `PointIndex` |
| `loci.py` | columnar `Loci` with converters, protocols, `head/tail/describe/columns/shape/to_numpy` |
| `interop.py` | boundary + converters (pandas, polars, arrow, bioframe, pyranges, pybedtools, cgranges, AnnData, Biopython, xarray cubes), liftover |
| `_table.py` | protocol mixin |
| `genes.py` | 0-based tables, `make` (polars/pandas, `chr_map`), `make_ucsc`, `from_frame`, `to_gtf`, `to_bed12`, `representative`, `to_pandas(table)` |
| `bedpe.py` | `Pairs`, `count_pairs`, `count_pairs_2d`, `pair_2d_*` |
| `architecture.py` | Pairs input, interval/graph backends, `from_edges/from_frame/from_scipy`, `to_networkx/igraph/scipy/anndata`, `pagerank`, `vertices_frame`; `to_legacy` removed |
| `signal.py`, `signal_draw.py` | bigwig backends, handles, groups as Loci/mask/rows |
| `motifs.py` | backend scanning, p-value thresholds, masked scan, positional profiles (the owner's WIP), stats and archetypes kept |
| `motifs_draw.py` | `plot_motif_heatmap(loci, M, names, ...)` added (the owner's WIP, 2.0 API) |
| `atlas.py` | `as_loci` inputs, vectorised bin ranges + bootstrap set-up, half-core workers, `to_pandas/to_arrow/columns/head/describe/_repr_html_` (TableMixin) |
| `bam.py` | region parsing from `locus.parse_region` |

### Left to do (in this order)

1. ~~motifs_draw~~, ~~atlas~~, ~~bam~~ — done.
4. **NEXT** `browserview.py` (not started; `backends.bigwig.open_bigwig` now
   returns any object with `stats_array` as is, so test stubs work): tracks may be Loci / anything `as_loci` takes, `Pairs` /
   BEDPE, `Genes` (draw from the three tables: transcripts per gene, exons and
   CDS from `features`), bigWig paths or open handles via `backends.bigwig`.
5. `architecture_draw.py`: `draw(A, region, layout='spring'|'circular'|'genomic',
   backend=None, merge_distance=...)` with `backends.graph.layout` + plain
   matplotlib (no graph-tool needed).
6. `hichip.py`, `se.py`, `view.py`, `igv.py`: imports (no `.columnar`), Genes are
   0-based now (drop the `- 1` GTF fixes in `view.py` / `igv.py`; use
   `Genes.to_bed12` / `representative`), signal via `Loci.signal`.
7. TableMixin + `head/describe/columns/_repr_html_` on `Pairs`, `Genes` (genes
   table), `Architecture` (edge table), `Atlas` (track table), `Library`.
8. `pyproject.toml`: version 2.0.0; deps add `narwhals`, `pyarrow`, drop
   `pyranges`; extras: `fast` (polars), `bam` (pysam), `viz` (logomaker),
   `interop` (bioframe, pyranges, pybedtools, anndata, biopython, networkx,
   igraph, MOODS-python, pyBigWig, ncls); package-data path is now
   `genomeblocks/_view.js`, `_view.css`. Update `environment.yml`,
   `conda-recipe/`, `meta.yaml`.
9. Tests (rewrite `tests/` for 2.0): parity across every backend of every
   family; round trips (pandas, polars, arrow, bioframe, pyranges, bedtool,
   anndata, parquet, BED, CSV); protocols (seaborn / plotly / altair / duckdb);
   the 1.x behaviours worth keeping (set ops, merge, nearest, liftover, tile,
   count_pairs, select_isoforms, annotate, prime_hubs, support, normalize with
   trans edges, atlas, motif profile binning, browser, bam).
10. Benchmarks: `benchmarks/bench_backends.py` timing every backend of every
    family through the genomeblocks API (data from `benchmarks/make_data.py`),
    results JSON → `benchmarks/build_site.py` → the docs Benchmarks page.
    Retire the 1.x-only scripts (`bench_prototype.py`, `prototypes/`, report
    builders) or move them to `benchmarks/archive/`.
11. Examples: re-run `examples/ar_foxa1_lncap` on the 2.0 API (data via its
    `download_data.sh`; hg38 references at the lackgrp paths in the notebook);
    update `examples/case_study_se`; drop `examples/columnar_prototype`.
12. Docs (`docs/`, Jekyll just-the-docs, deployed by `.github/workflows/pages.yml`
    on merge to `main`): rewrite every page for 2.0 — home, install,
    quickstart, concepts (tables, Genome, rows as join key, coordinates),
    interoperability, backends, one guide + API page per module, design pages
    (regenerate diagrams with `docs/_tools/build_diagrams.py`), walkthrough from
    the re-run notebook, benchmarks, release notes 2.0 with a 1.x → 2.0
    migration table. Keep the theme, fonts and colour scheme.
13. README, CONTRIBUTORS; delete this file; open the PR (the owner merges;
    merging through the API is blocked for agents).

## How to work on it

```bash
git fetch origin && git switch v2-columnar
python -m venv .venv && . .venv/bin/activate
pip install -e . polars pyarrow narwhals bioframe pyranges pybedtools anndata xarray \
    biopython MOODS-python pyBigWig ncls pyfaidx pysam networkx igraph duckdb \
    seaborn plotly altair logomaker pytest
# graph-tool and cgranges are conda-only:
#   conda install -c conda-forge graph-tool ; conda install -c bioconda cgranges bedtools
pytest -q
```

On the owner's machine the worktree is `.claude/worktrees/v2` and a ready
environment with every backend is `.claude/tools/v2-venv/bin/python`
(run with `PYTHONPATH=$PWD`). pybedtools writes temp files: point
`pybedtools.set_tempdir()` / `TMPDIR` at a project directory (the owner does
not want files in /tmp).

Commit as the owner (no AI co-author trailer). Push to `v2-columnar`; never
push to `main` directly.
