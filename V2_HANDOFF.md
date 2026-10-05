# genomeblocks 2.0 — handoff notes (branch `v2-columnar`)

This file says what 2.0 is, what is done, what is left, and how to pick it
up. **Delete it (or move it to `docs/archive/`) before the 2.0 release.**
Rollback point: tag `v1.1.0`.

Last updated: 2026-10-05 (checkpoint 6, see the log at the end).

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
test suite (`tests/`, every test asserting values; 204 with every backend
installed, as the parity tests run once per installed engine), and every item
of the previous "Left to do" that the agent can do is done (checkpoint 6):

| area | what is there |
|---|---|
| core tables | `Genome`, `Loci`, `Locus`/`LocusView`, `Genes` (three linked tables), `Pairs`, `Architecture` (vertex + sorted edge table, cis blocks + trans block), `Atlas` (bins × tracks CSR) — all with `columns / shape / head / tail / describe / _repr_html_` and the Arrow / interchange / narwhals protocols |
| boundary | `interop.as_loci` on every input kind; `to_*` converters; file readers for BED / narrowPeak / CSV / TSV / parquet with `#`, `track`, `browser` lines skipped and a `#chrom` header used as names; peak files named `.bed` (ChIP-Atlas, ENCODE) get the narrowPeak / broadPeak columns; `liftover` |
| backends | six families, parity-tested against brute force; the whole suite also passes with graph-tool, cgranges and bedtools installed (owner's machine, `.claude/tools/v2-venv`) and in a conda environment with pandas 3 / numpy 2.5 |
| modules | `signal` (+ `tmm`, `signal_draw`), `motifs` (+ `motifs_draw`, `plot_motif_heatmap`), `bedpe` (pairs files, `count_pairs`, `count_pairs_2d`), `atlas`, `se`, `hichip`, `bam`, `browserview`, `architecture_draw` (any graph backend), `view`, `igv` |
| packaging | `pyproject.toml` 2.0.0 with the extras; `environment.yml`; `conda-recipe/` (ready Bioconda recipe); `python -m build` gives an sdist and a wheel that pass `twine check` |
| tests + CI | `.github/workflows/tests.yml` (Python 3.10 / 3.12, `test` extra and every extra); green on every push of checkpoint 6 |
| docs | site rewritten for 2.0; every python block outside the walkthrough and benchmarks pages re-run against the test-suite fixtures and its `# ->` outputs fixed where they were stale; every documented signature checked against `inspect.signature`; identifiers in the prose and numbers in the figure captions checked against the code and the SVGs; the site builds with 0 broken internal links. A version menu in the header keeps the 1.1 docs online under `/1.1/` (built by `pages.yml` from `f17d343`, see `docs/README.md` → Versions) |
| examples | `ar_foxa1_lncap` re-run on 2.0 with the ChIP-Atlas data and the hg38 references (walkthrough page, figures and numbers from that run; every walkthrough snippet runs on that data); `case_study_se` ported to 2.0 and re-executed (synthetic benchmark data + MACS3); `columnar_prototype` removed |
| benchmarks | every bench re-run on 2.0 on the owner's machine (one at a time, ≤ 8 workers) and `docs/benchmarks/index.md` regenerated; external tools are found on `PATH`; the motif bench times genomeblocks on MOODS and on lightmotif |
| README | rewritten for 2.0 |

### Left to do — the owner's

1. **Merge the pull request** `v2-columnar → main`. `pages.yml` then deploys
   2.0 at the site root and 1.1 under `/1.1/`. This file goes before the
   release (delete it, or move it to `docs/archive/`).
2. **Release 2.0.0**: tag `v2.0.0` on `main`, `python -m build`,
   `twine upload dist/*`; fill the sdist sha256 in
   `conda-recipe/genomeblocks/meta.yaml`; submit the recipe to
   bioconda-recipes (steps in `conda-recipe/README.md`).
3. Optional: build GIGGLE from source and set `GIGGLE=` for the Atlas
   benchmark's external baseline.

The merged branches (`benchmark-fixes`, `benchmarks-suite`,
`columnar-prototype`, `docs-design-benchmarks`, `release-1.1.0`,
`bioconda-moods`) stay: the owner decided on 2026-10-05 not to delete them.
Do not delete them.

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
| 4 — cloud session | `65358b0` value-level tests, lint · `2c72c80` docs first pass · `cadef94` benchmarks lint · `3d69265` merge of `bioconda-moods` (`71c9f39`: MOODS the motif engine) · `4bcf98a` README extras, Arrow-stream regression test · `832be19` this file | 179 tests green locally; CI run 2 green |
| 5 — cloud session | `acf0843`, `f1e99ca` eleven defects found by the docs examples, each with a failing-first test · `b602de4`, `379b735` CI builds MOODS-python on the runner (a wheel cached from another runner's CPU crashed with SIGILL on the 3.12 jobs) · `d718de1`, `db5c615` docs second pass: API table pages, interop and backends API pages, guide and design review, 14 new and 4 corrected diagrams, MOODS statements · `779eb6c` to_cgranges install hint, Pairs.describe integers, docs follow-ups · `6876bcf` .tmp ignored · docs structure fixes (credits.md BibTeX in raw tags, nav orders, stale motif lines) | 189 tests green locally; CI runs 6 and 7 green on both Python versions and both extras sets; docs verification interrupted (items 7–9 above) |
| 6 — owner's machine | `eab7634` graph-tool PageRank read freed memory · `fb9e722` describe() without a blank header row, `Locus.overlaps(other)` · `b160658` every docs example re-run against the test fixtures, stale outputs, signatures and descriptions fixed · `559e150`, `47aaba9` docs version menu (1.1 under `/1.1/`) · `61f4e1e` peak files named `.bed` load with `keep=True` · `1b6a87b`, `e65598c` benchmarks on the 2.0 names, tools on `PATH` · `f7f4948` AR / FOXA1 re-run on the real data, walkthrough regenerated, `columnar_prototype` removed · `ff771ba` `case_study_se` on 2.0 · `03812f7` every benchmark re-run, Benchmarks page regenerated · this file | 204 tests green with every backend (and in a pandas 3 conda env); CI green on every push; docs, examples and benchmarks verified; the merged branches stay (the owner's decision) |

Branches: `v2-columnar` (this work), `main` (1.1.0). Every other branch is
fully merged and kept on purpose (the owner's decision); do not delete them.
