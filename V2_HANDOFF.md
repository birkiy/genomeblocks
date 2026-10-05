# genomeblocks 2.0 — handoff notes (branch `v2-columnar`)

This file says what 2.0 is, what is done, what is left, and how to pick it
up. **Delete it (or move it to `docs/archive/`) before the 2.0 release.**
Rollback point: tag `v1.1.0`.

Last updated: 2026-10-05 (checkpoint 4, see the log at the end).

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

### Working rules (the owner's)

- Keep the README, the documentation, the benchmarks and the coding hygiene
  intact with every change (short docstrings that say what the function does,
  no comments narrating a change, errors that name the fix, no legacy shims,
  `ruff check --select F` clean).
- Tests must be designed to find real bugs: assert the values (against a
  brute-force oracle or a hand-computed expectation), not shapes or "import
  works".
- Commit as the owner (no AI co-author trailer, no model names anywhere).
  Push to `v2-columnar`; never push to `main` directly. The owner merges the PR.
- Keep the checkpoints in this file for the next session.

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
  (`AUTO` in `backends/__init__.py`). Unknown/missing requested backend raises
  with the install command. `gb.backends()` (the module is callable) lists them.
  Parity rules every engine is normalised to: half-open overlaps (zero-length
  intervals included), bigWig bins on integer edges `floor(n*b/n_bins)`,
  population std, networkx parallel edges summed, canonical component labels.
- **Motif engine: MOODS by default, lightmotif optional** (checkpoint 4, from
  the `bioconda-moods` branch). lightmotif is PyPI-only, so it cannot be a
  Bioconda dependency; MOODS reports the same hits. Motif files (JASPAR raw and
  bracketed, TRANSFAC, uniprobe, MEME, `.gz`) are parsed by genomeblocks
  itself, log-odds reproduce lightmotif's float32 arithmetic bit for bit
  (`backends.motifs.logodds_matrix`), and p-value cutoffs come from our own
  exact score distribution (`threshold_from_pvalue`, within 0.004 bits of
  MOODS). `pip install genomeblocks` therefore has **no motif engine**:
  `genomeblocks[motifs]` (MOODS-python, compiled) or `genomeblocks[lightmotif]`
  (prebuilt wheels); the conda package ships `moods`. `AUTO['motifs'] =
  ['moods', 'lightmotif']`; biopython is the third engine.
- **Boundary:** `interop.as_loci()` / `interop.frame()` (narwhals) normalise every
  input once; `Loci.from_frame` finds columns by name, falls back to position
  only for a text+int+int header-less frame, otherwise raises naming the
  accepted spellings. A protocol speaker that is not a table (a pandas 3
  Series streams one column) is not a frame: `as_loci` names the input type.
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
- Core deps: `numpy`, `pandas>=2`, `pyarrow`, `narwhals>=1`, `scipy`,
  `matplotlib`, `pybigtools`, `cooler`, `tqdm`. `pyranges` and `lightmotif`
  are optional backends. Extras: `motifs`, `lightmotif`, `fast`, `bam`, `viz`,
  `interop`, `test`, `all` (see `pyproject.toml`).
- Signal windows are not clipped at chromosome edges: the bigWig handle bins
  the full grid and bases outside the chromosome have no data, so every
  backend agrees (`signal._windows`).
- `se.call_se` takes `backend=` and `return_all=`; `stitch_peaks` is a public
  function; `nearest_gene_within` re-codes the TSS onto the regions' Genome.

## Status

### Done

Everything on the 1.x surface is ported to the tables and covered by the
test suite (`tests/`, 179 tests, every one asserting values):

| area | what is there |
|---|---|
| core tables | `Genome`, `Loci`, `Locus`/`LocusView`, `Genes` (three linked tables), `Pairs`, `Architecture` (vertex + sorted edge table, cis blocks + trans block), `Atlas` (bins × tracks CSR) — all with `columns / shape / head / tail / describe / _repr_html_` and the Arrow / interchange / narwhals protocols |
| boundary | `interop.as_loci` on every input kind; `to_*` converters; file readers for BED / narrowPeak / CSV / TSV / parquet with `#`, `track`, `browser` lines skipped and a `#chrom` header used as names; `liftover` |
| backends | six families, parity-tested against brute force (`tests/test_loci.py`, `test_signal.py`, `test_motifs.py`, `test_architecture.py`); bedtools / cgranges / graph-tool paths exercised only through their ImportError here (conda-only; see "On the owner's machine") |
| modules | `signal` (+ `tmm`, `signal_draw`), `motifs` (+ `motifs_draw`, `plot_motif_heatmap`), `bedpe` (pairs files, `count_pairs`, `count_pairs_2d`), `atlas`, `se`, `hichip`, `bam`, `browserview`, `architecture_draw` (any graph backend), `view`, `igv` |
| packaging | `pyproject.toml` 2.0.0 with the extras above; `environment.yml` (conda engines); `conda-recipe/` (ready Bioconda recipe with `run_test.py` and `build-local.sh`) |
| tests + CI | `tests/` rewritten; `.github/workflows/tests.yml` runs Python 3.10 and 3.12 with the `test` extra and with every extra |
| benchmarks | `benchmarks/bench_backends.py` (every backend of every family through the API, with an agreement column), the 1.x scripts ported to 2.0, `build_site.py` with a Backends section, 1.x-only studies in `benchmarks/archive/` |
| docs | site rewritten for 2.0 (home, installation, quickstart, concepts, release notes with the 1.x → 2.0 migration table, credits, Interoperability and Backends pages, design pages with regenerated diagrams, guide and API pages, walkthrough snippets on the 2.0 API); every code block executed against synthetic fixtures |
| README | rewritten for 2.0 |

### Left to do

Owner's machine or real data needed:

1. **Confirm CI is green on the current head** (run 2 on `4bcf98a` was in
   progress at this checkpoint; run 1 failed only on the pandas 3 Series case
   fixed by the `bioconda-moods` merge).
2. **Examples** (`examples/ar_foxa1_lncap`, `examples/case_study_se`): re-run on
   the 2.0 API with the data from `download_data.sh` and the hg38 references
   at the lackgrp paths; drop `examples/columnar_prototype`. The walkthrough
   pages keep the 1.x figures and numbers until then (a note says so on
   `docs/walkthrough/index.md`).
3. **Benchmarks**: `benchmarks/make_data.py`, then `benchmarks/run_all.sh`
   (needs pyjaspar and the conda engines for the full table) and
   `python benchmarks/build_site.py` to regenerate `docs/benchmarks/index.md`,
   which is still the 1.1 page.
4. **Run the suite with the conda-only engines** (graph-tool, cgranges,
   bedtools) in `.claude/tools/v2-venv`: their parity is asserted by the tests
   but could not be run in the cloud environment.
5. **Release**: `python -m build`, upload 2.0.0 to PyPI, fill the sdist sha256
   in `conda-recipe/genomeblocks/meta.yaml`, submit to bioconda-recipes
   (steps in `conda-recipe/README.md`).
6. **Open the PR** `v2-columnar → main`; the owner merges (`pages.yml` deploys
   the docs on merge to `main`). Delete or archive this file first.

In progress in the current session (see the checkpoint log; whatever is not
listed as done below was interrupted and should be checked first):

7. Package defects found while the docs examples were run, each with a
   regression test: `igv_html` chrom sizes from the bigWigs when nothing else
   gives them; `plot_profiles` kb labels; `compare_heatmap` length error naming
   the plotted columns; `Architecture.__or__` column order; a warning when
   `normalize()`'s power-law fit fails; `architecture_draw` axis limits;
   `Atlas.search` duplicate column; `Atlas.load` NaN round trip; `se.call_se`
   `progress=` conflict; `interop.frame` on a file-backed BedTool with header
   lines; two docstrings (`_bbi.stats_array`, `compare_motifs`).
8. Docs part 2: `docs/api/{index,loci,genes,bedpe}.md` rewritten,
   `docs/api/{interop,backends}.md` new; the guide and design pages the cut-off
   writers left reviewed; diagrams for the concepts (row as join key, Genome
   re-coding, coordinates and the TSS rule, the boundary), the backend
   dispatch and the interop path, plus the figures the reviewers asked for;
   every "lightmotif is the default" statement corrected to MOODS.
9. Docs verification pass (links, front matter, includes, examples re-run,
   identifiers checked against the source) and the fixes it produces.

## How to work on it

```bash
git fetch origin && git switch v2-columnar
python -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"        # every pip engine, both motif engines, pytest
# graph-tool, cgranges and bedtools are conda-only:
#   conda env create -f environment.yml && conda activate genomeblocks && pip install -e .
export TMPDIR=$PWD/.tmp && mkdir -p .tmp     # pybedtools temp files stay in the project
pytest -q
ruff check --select F genomeblocks tests benchmarks docs/_tools
python docs/_tools/build_diagrams.py         # regenerates docs/_includes/diagrams/*.svg
```

On the owner's machine the worktree is `.claude/worktrees/v2` and a ready
environment with every backend is `.claude/tools/v2-venv/bin/python`
(run with `PYTHONPATH=$PWD`). pybedtools writes temp files: point
`pybedtools.set_tempdir()` / `TMPDIR` at a project directory (the owner does
not want files in /tmp).

Docs: Jekyll just-the-docs under `docs/` (`docs/README.md` has the
conventions); links only as `{{ '/guide/loci/' | relative_url }}`; diagrams
are generated, never hand-edited. Benchmarks: `benchmarks/README.md`.

## Checkpoint log

| checkpoint | commits (oldest → newest) | state |
|---|---|---|
| 1 — owner's machine | `63d5cc3` columnar-only core, swappable backends, interop boundary | WIP: core and backends written, no tests |
| 2 — owner's machine | `3a2c9c9` motif heatmap, columnar Atlas, bam region parsing | WIP: usage limit reached; items 1–3 of the old list done |
| 3 — cloud session | `e0cd91c` audit defects fixed · `4961e7b` remaining modules ported, table protocols on every container · `2ea848a` packaging 2.0.0 · `32db476` the 2.0 test suite + CI · `00c9f17`, `ac6e407` se knee on short curves · `87672a2` `bench_backends.py`, 1.x studies archived · `b9c3070` merge of checkpoint 2 | every module on the tables; suite green locally; CI red on 3.12 (pandas 3 Series) |
| 4 — cloud session | `65358b0` value-level tests, lint · `2c72c80` docs first pass · `cadef94` benchmarks lint · `3d69265` merge of `bioconda-moods` (`71c9f39`: MOODS the motif engine) · `4bcf98a` README extras, Arrow-stream regression test | 179 tests green locally; CI run 2 in progress; docs part 2 and the defect fixes running |

Branches: `v2-columnar` (this work), `bioconda-moods` (merged at checkpoint 4,
can be deleted), `main` (1.1.0).
