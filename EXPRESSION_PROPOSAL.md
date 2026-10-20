# Proposal: expression tables and differential expression (genomeblocks 2.1)

*Status: draft for discussion, 2026-10-20. Branch `expression`. Nothing here is built yet.*

genomeblocks 2.0 has no place for RNA-seq. A `Genes` object knows where every
gene and transcript is, `Architecture.annotate` names the gene of every CRE,
`select_isoforms` picks the isoform open chromatin supports — and the one
measurement that says whether any of that *matters*, expression, has to be
joined in by hand from a pandas frame. This proposal adds one table,
`Expression`, that reads the output of the common quantifiers, lines its rows
up with a `Genes` object (or any `Loci`), summarises transcripts to genes the
way tximport does, normalises, and runs DESeq2 / edgeR / limma through the
package's backend mechanism, so a differential-expression result is an array
aligned to gene rows like every other result in the package.

```python
import genomeblocks as gb

genes = gb.Genes.make("gencode.v49.annotation.gtf")
E = gb.Expression.make("salmon/*/quant.sf", genes=genes, level="gene")   # tximport, aligned to genes.genes
E.samples["condition"] = ["WT", "WT", "WT", "KO", "KO", "KO"]

res = E.de("~ condition", contrast=("condition", "KO", "WT"))             # DESeq2 (pydeseq2); one row per gene row
up = (res.padj < 0.05) & (res.log2fc > 1)                                 # a mask over gene rows ...

tss = genes.get_tss()
tss["log2fc"] = res.log2fc.to_numpy()                                      # ... so no join anywhere
cre_up = cre & tss[up].slop(1000)                                          # CREs at promoters of induced genes
A.vp["log2fc"] = E.lookup(A.vp.gene, res.log2fc)                           # after A.annotate(genes)
```

## Contents

1. [Why it fits 2.0](#1-why-it-fits-20)
2. [The `Expression` table](#2-the-expression-table)
3. [Inputs](#3-inputs)
4. [Transcripts to genes](#4-transcripts-to-genes-tximport-on-the-link-column)
5. [Normalisation](#5-normalisation)
6. [Differential expression: the `de` backend family](#6-differential-expression-the-de-backend-family)
7. [Integration with the rest of the package](#7-integration-with-the-rest-of-the-package)
8. [Drawing](#8-drawing)
9. [Module layout, dependencies, packaging](#9-module-layout-dependencies-packaging)
10. [Tests and benchmarks](#10-tests-and-benchmarks)
11. [Documentation and example](#11-documentation-and-example)
12. [Phases](#12-phases)
13. [Non-goals](#13-non-goals)
14. [What I checked locally](#14-what-i-checked-locally)
15. [Decisions needed](#15-decisions-needed)

---

## 1. Why it fits 2.0

Three of the 2.0 principles decide the design almost on their own.

**The row number is the join key.** `Genes` is three `Loci` linked by row:
`transcripts['gene'][k]` is the gene row of transcript `k`. An expression
matrix is naturally a fourth table in that family: *one row per gene row* (or
per transcript row), one column per sample. Once `E.counts[i]` is gene row
`i`, every downstream join is array indexing: `genes.get_tss()[mask]`,
`res.log2fc.to_numpy()[genes.rows(names)]`, `A.vp["x"] = ...`. No
`gene_id` merge, no version-suffix mismatch at analysis time (that is resolved
once, when the table is built, and reported).

**Boundary first.** Quantifier outputs come in a dozen layouts. One
constructor, `Expression.make`, sniffs them (as `Loci.make` sniffs BED /
narrowPeak / bedGraph) and every public function takes the finished table.
Nothing deeper in the package learns what a `quant.sf` is.

**Backends: explicit choice, no silent fallback.** DESeq2, edgeR and limma
are *different methods* and give different answers; pydeseq2 and R's DESeq2
are *different engines* of the same method and give nearly, not exactly, the
same answer. The existing rules handle both: the method is a named argument
(never guessed), the engine is a `backend=` (automatic only where one
candidate exists), and a requested engine that is missing raises with its
install command. Results are normalised to one schema whatever engine ran, as
interval results are.

Two smaller principles carry over unchanged: processing stays matplotlib-free
(`expression.py` + `expression_draw.py`, see
`project_genomeblocks_layout`), and heavy imports stay lazy (pandas,
pydeseq2, subprocess-to-R are imported inside the functions that need them).

## 2. The `Expression` table

```python
class Expression(TableMixin):
    rows      : Genes | Loci | None   # what the rows are aligned to
    level     : "gene" | "transcript" | "loci"
    ids       : object[n]             # feature ids as read (versions kept)
    names     : object[n] | None      # gene_name when known
    layers    : dict[str, float64[n, m]]   # "counts" always; "tpm", "length", ... when the input has them
    samples   : pandas.DataFrame      # m rows, index = sample names, columns = metadata
    filename  : str | None
```

- **Shape is features × samples** (n × m), like the signal cube
  (rows × tracks × bins), `Atlas` (bins × tracks) and Bioconductor's
  SummarizedExperiment. `to_anndata()` transposes to samples × features
  (AnnData's obs × var), `from_anndata()` transposes back.
- **Layers.** `counts` is always present (float64: salmon's `NumReads` are
  fractional; the DE engines receive them rounded, as tximport does). `tpm`
  and `length` (effective length, per feature per sample) are filled from
  salmon / kallisto / RSEM; featureCounts gives a per-feature `Length`
  (broadcast to n × m). `E.counts`, `E.tpm` are attribute shortcuts; `E.X`
  is the default layer (`counts`).
- **Alignment.** With `genes=` (or `loci=`), `len(E) == len(rows table)` and
  row `i` is row `i` of that table. Features of the input that match no row
  are dropped and counted; rows of the table absent from the input are NaN in
  every layer ("not quantified", distinct from 0 "quantified, nothing
  there"). Both counts are printed (`verbose=True`) and kept as `E.unmatched`
  (ids) and `E.missing` (a mask over the rows). `strict=True` raises instead.
  Without `genes=`/`loci=`, rows are the input's features in file order and
  `E.align(genes)` produces the aligned table later.
- **Id matching** (`key="auto"`): exact `gene_id` / `transcript_id`; if
  fewer than half match, both sides are stripped of Ensembl version suffixes
  (`ENSG00000141510.17 → ENSG00000141510`; GENCODE's `_PAR_Y` kept apart);
  GENCODE-style pipe headers from salmon indexes built without `--gencode`
  (`ENST…|ENSG…|…|TP53-201|TP53|1234|protein_coding|`) are cut at the first
  `|`. `key="gene_name"` matches a name-keyed table (duplicated names: first
  gene row wins, reported). The rule that was used is printed and stored in
  `E.meta`.
- **Samples.** `E.samples` is a pandas DataFrame indexed by sample name in
  column order; `samples=` in `make` takes a DataFrame or a TSV/CSV path
  (joined on the index, missing samples raise). Column order of `X` and row
  order of `samples` are one thing.
- **Selection.** `E[:, mask]` / `E.take_samples(idx)` keeps the alignment;
  `E.take(rows)` returns a table whose `rows` link is dropped (ids kept),
  since it no longer lines up with the full gene table. Masks over genes are
  meant to be applied to results and TSS tables, not to `E`.
- **Lookups.** `E["TP53"]` / `E["ENSG…"]` → a pandas Series over samples
  (through `Genes._key_index`); `E[["TP53", "MYC"]]` → a small frame;
  `E.lookup(keys, values=None, layer="counts", agg="mean")` → an array
  aligned to `keys` (gene names / ids, `""` or unknown → NaN), which is the
  one helper the CRE-side joins need.
- **Table protocols** (via `TableMixin.to_arrow`): the wide table of the
  default layer, `gene_id`, `gene_name`, then one column per sample, so
  `pl.DataFrame(E)`, `duckdb.sql("select * from E")`, `sns.*` work.
  `to_pandas(layer=, long=False)`, `to_polars`, `to_arrow`, `to_anndata`,
  `describe()`, `head()`, `_repr_html_` (n × m, layers, matched / missing).
- **Persistence.** `save(dir)` → `samples.parquet`, `features.parquet`
  (ids, names, row link), one parquet per layer, `meta.json` (level, key
  rule, the `Genes` filename). `Expression.load(dir, genes=genes)` re-checks
  the alignment by ids.

Memory: GENCODE v49 is 254k transcripts × 6 samples × 3 layers ≈ 37 MB;
63k genes × 100 samples × 3 layers ≈ 150 MB. float64 throughout.

## 3. Inputs

`Expression.make(source, *, genes=None, loci=None, level=None, samples=None,
key="auto", sample_names=None, strand=None, verbose=True)`.

`source` is one of:

| Source | How samples are found |
|---|---|
| a glob or list of per-sample files (`"salmon/*/quant.sf"`, `["a/quant.sf", "b/quant.sf"]`) | the parent directory name (salmon / kallisto layout), else the file stem |
| a `{sample: path}` dict | the keys |
| one multi-sample file (featureCounts, STAR-merged, nf-core matrices, any matrix) | the header (`sample_names=` to rename; BAM paths are reduced to their stem) |
| a frame (pandas, polars, pyarrow, anything `interop.frame` takes) via `from_frame(df, ids=, names=)` | the numeric columns |
| an AnnData via `from_anndata(adata, layer=None)` | `obs_names`; `obs` becomes `samples` |

Formats are sniffed from the header line (the `tables` backend parses; polars
when installed, else pandas, identical columns either way):

| Format | Recognised by | Layers | Notes |
|---|---|---|---|
| salmon `quant.sf` / `quant.genes.sf` | `Name Length EffectiveLength TPM NumReads` | counts, tpm, length | level from the file name; `quant.genes.sf` is salmon's own summary |
| kallisto `abundance.tsv` | `target_id length eff_length est_counts tpm` | counts, tpm, length | `abundance.h5` not read (bootstraps out of scope) |
| RSEM `.genes.results` / `.isoforms.results` | `gene_id transcript_id(s) length effective_length expected_count TPM FPKM` | counts, tpm, length | level from the header |
| featureCounts | first line `# Program:featureCounts`, header `Geneid Chr Start End Strand Length …` | counts, length | all samples in one file; the `Chr/Start/End` columns are ignored (the `Genes` has the coordinates); an SAF of peaks gives a CRE count matrix (`loci=`) |
| STAR `ReadsPerGene.out.tab` | four `N_*` rows, 4 columns | counts | `strand="none" / "forward" / "reverse"`; `"auto"` picks the stranded column with the larger total and says so |
| HTSeq-count | 2 columns, trailing `__no_feature …` rows | counts | the `__*` rows are dropped and reported |
| nf-core/rnaseq matrices | `gene_id gene_name <samples>` / `tx gene_id <samples>` | counts or tpm (from the file name: `*_counts*`, `*_tpm*`, `*_length_scaled*`) | several files of one run can be passed together to fill several layers |
| a plain matrix (TSV / CSV / parquet) | first column = ids, the rest numeric | counts | `ids=`, `names=` to name columns |

Sniffing failures raise with the recognised headers and the `from_frame`
escape hatch, in the style of `backends.tables._parse_error`.

## 4. Transcripts to genes: tximport on the link column

With `genes=` and a transcript-level input, `level="gene"` (the default
for a `quant.sf`) summarises exactly as `tximport::summarizeToGene` does, and
the three tables make it a few `np.add.at` calls over
`transcripts['gene']`:

```python
E_tx = gb.Expression.make("salmon/*/quant.sf", genes=genes, level="transcript")   # 254k × 6
E = E_tx.to_genes(how="length_scaled_tpm")                                       # 63k × 6
```

| `how` | gene counts | matches tximport |
|---|---|---|
| `"counts"` | plain sum of transcript counts; gene length = abundance-weighted mean of transcript effective lengths (a gene with zero TPM in a sample gets the plain mean of its transcripts' lengths) | `countsFromAbundance="no"` — needs length offsets in the model, see §6 |
| `"scaled_tpm"` | gene TPM scaled to the sample's library size | `"scaledTPM"` |
| `"length_scaled_tpm"` (default) | gene TPM × mean gene length, scaled to the library size | `"lengthScaledTPM"` — the recommended input for DESeq2 / edgeR without offsets |

Gene TPM is the sum of transcript TPM in every mode. Transcripts whose gene
row is unknown (unmatched ids) are dropped and counted. The result is
aligned to `genes.genes`; `E_tx` stays aligned to `genes.transcripts`, and
both can be kept. Correctness is pinned against tximport's numbers on the
test fixture (§10).

## 5. Normalisation

Deterministic, numpy, no backend: every engine would give the same number.

| Call | Gives |
|---|---|
| `E.cpm()`, `E.log_cpm(prior=2)` | counts per million by library size, optionally TMM-scaled (`norm="tmm"`) |
| `E.tpm` | the layer when read; computed from `counts` and `length` otherwise |
| `E.size_factors()` | DESeq2's median-of-ratios factors |
| `E.tmm()` | edgeR's TMM factors — `signal._tmm_norm_factors` already implements the rule (30 % M trim, 5 % A trim, inverse-variance weights, reference = upper quartile closest to the mean); it moves to a shared `_norm.py` and `signal.tmm` keeps calling it |
| `E.filter_by_expr(group=None, min_count=10, min_total=15, ...)` | edgeR's `filterByExpr` rule as a mask over the rows |
| `E.vst()` | *not* reimplemented: dispatches to the `de` engine (pydeseq2 and DESeq2 both have it) |

`cpm` / `log_cpm` return (n × m) arrays; `size_factors` / `tmm` return
length-m arrays. Nothing is stored unless asked (`E.layers["logcpm"] = ...`).

## 6. Differential expression: the `de` backend family

```python
res = E.de("~ batch + condition", contrast=("condition", "KO", "WT"),
           method="deseq2", backend=None, filter=None, **engine_kwargs)
```

**Method vs engine.** `method` names the model and is one of `"deseq2"`
(Wald test), `"edger"` (quasi-likelihood F-test, `glmQLFit` /
`glmQLFTest`), `"limma"` (voom + `eBayes`; `trend=True` for limma-trend).
`backend` names the engine; a new family `de` in `backends._FAMILIES`:

| Engine | Runs | Install | Notes |
|---|---|---|---|
| `pydeseq2` | deseq2 | `pip install 'genomeblocks[de]'` (pydeseq2 ≥ 0.5; bioconda has it) | pure Python, the pip default; apeGLM-style shrinkage through `lfc_shrink=True` |
| `r` | deseq2, edger, limma | R + Bioconductor packages; `Rscript` on `PATH`, or `rscript="/path/to/Rscript"` / `GENOMEBLOCKS_RSCRIPT` | three bundled scripts under `genomeblocks/backends/_r/`; counts, samples and the design cross as TSV in a temporary directory; the result and `sessionInfo` versions come back the same way. No rpy2. |
| `inmoose` (phase 3) | deseq2, edger (LRT only), limma (trend; voom not confirmed in the port) | `pip install inmoose` | pure-Python ports; useful as a second Python engine for parity tests |

`AUTO["de"]` is per method and has exactly one entry each — `deseq2 →
pydeseq2`, `edger → r`, `limma → r` — because no two engines of one method
are numerically identical, so rule 1 ("automatic only where the answers are
identical") forbids an "else" chain. A method the engine cannot run raises
`backends.unsupported("de", "pydeseq2", "edger", ["r", "inmoose"])`, like
`ncls` and `merge` today. `installed("de", "r")` needs a probe that is not an
import: `_FAMILIES` grows a third form, `"rscript:DESeq2,edgeR,limma"`, which
`installed()` resolves by running `Rscript -e 'requireNamespace(...)'` once
and caching per package. `gb.backends()` lists the family; the method each
engine runs is in the docs and in the error message.

**Design and contrasts.** `design` is a formula over `E.samples` columns
(Wilkinson syntax, as pydeseq2 ≥ 0.5 and R take it). `contrast=(factor,
a, b)` gives the log2 fold change of `a` over `b` and sets `b` as the
reference level; `coef="conditionKO"` names a coefficient instead;
`contrasts=[...]` runs several in one fit and returns a dict. Numeric contrast
vectors are phase 3. Continuous covariates pass through unchanged.

**Filtering.** DESeq2 does its own independent filtering, so
`method="deseq2"` sends every row with counts. edgeR and limma expect
`filterByExpr`; for them `filter=None` means "apply `E.filter_by_expr`
grouped by the contrast factor" and says how many rows were kept; `filter=`
a mask overrides, `filter=False` sends everything. Rows that are NaN
(not quantified) are never sent.

**Length offsets.** For `to_genes(how="counts")` tables the per-gene
length matrix is passed as DESeq2 `normalizationFactors` / edgeR `offset`
(what tximport's `DESeqDataSetFromTximport` does); for the default
`length_scaled_tpm` tables nothing extra is needed. pydeseq2 has no offset
matrix, so `how="counts"` + `pydeseq2` raises and names the alternatives.

**One result schema.** A pandas DataFrame with one row per row of `E`
(index = row number, NaN where a gene was not tested), columns:

| Column | deseq2 | edger | limma |
|---|---|---|---|
| `gene_id`, `gene_name` | from the `Genes` | | |
| `base_mean` | `baseMean` | mean CPM | mean log-CPM (`AveExpr`), back-transformed |
| `log2fc` | `log2FoldChange` (shrunk when asked) | `logFC` | `logFC` |
| `lfc_se` | `lfcSE` | NaN | NaN |
| `stat` | Wald statistic | F | moderated t |
| `pvalue`, `padj` | `pvalue`, `padj` (BH) | `PValue`, `FDR` | `P.Value`, `adj.P.Val` |

`res.attrs` records `method`, `backend`, `design`, `contrast`, engine
versions, `n_tested`, and the filter. Because the index is the gene row, a
mask on the frame is a mask on `genes.genes`.

**The R bridge in one paragraph.** `backends/de.py` writes `counts.tsv`
(rounded), `samples.tsv`, and a `job.json` (design, contrasts, method
options) to a `tempfile.TemporaryDirectory()`, runs `Rscript <method>.R
<dir>` with a timeout, reads `result_<contrast>.tsv` and `versions.json`
back, and re-raises R's stderr as a `RuntimeError` naming the script. The
scripts are plain R with no package besides the method's; they are run by
the test suite against the Python engines where both exist.

## 7. Integration with the rest of the package

- **`Genes`.** `E = Expression.make(..., genes=genes)` is the link, one
  direction (as `Architecture` holds its `loci`); `Genes` gets no
  back-reference. `Genes.rows(names, missing=None)` gains `missing=-1`
  so vectorised lookups need no try/except. `select_isoforms(expr=E_tx,
  min_tpm=1.0, min_frac=0.5)` takes transcript TPM as isoform evidence
  with the same support rule as the bigWig path (`min_signal` / `min_frac`
  → `tss_score`, `tss_support`, `canonical`), alone or together with ATAC
  peaks. That is the most direct use of the transcripts table the package
  has, and it is a few lines once `E_tx` is row-aligned.
- **`Loci`.** `genes.get_tss()` carries a `gene` row column already, so
  `tss["tpm"] = E.tpm.mean(axis=1)` and `tss[mask]` are the whole story.
  featureCounts over an SAF of CREs (`Expression.make("peaks.counts",
  loci=cre)`) gives a count matrix aligned to CRE rows and `E.de(...)`
  differential accessibility, with no new code: the row axis is "a `Genes`
  table or a `Loci`".
- **`Architecture`.** After `A.annotate(genes)`, `A.vp.gene` is a name
  per vertex; `A.vp["log2fc"] = E.lookup(A.vp.gene, res.log2fc)` puts the
  assigned gene's fold change on every CRE; `A.prime_hubs()` and
  `A.support(genes)` are unchanged. A `Architecture.add_expression(E, res,
  name=)` convenience can wrap that one line if wanted.
- **`View` / `browser` (phase 3).** `View.genes(color=res.log2fc)` colours
  gene models by a per-gene-row vector; `browser` gets a `genes_color=`.
- **`interop`.** `to_anndata` / `from_anndata` (samples obs × genes var,
  layers kept), `to_pandas(long=True)` for seaborn, Arrow / narwhals
  protocols through `TableMixin`. `as_loci(E)` is *not* defined: an
  `Expression` is not an interval table; `E.rows` is.

## 8. Drawing

`expression_draw.py` (matplotlib, imported only here), attached as methods
where the receiver is obvious:

| Function | Method | Shows |
|---|---|---|
| `volcano(res, padj=0.05, lfc=1, label=[...])` | — | log2fc vs −log10 padj, named genes labelled |
| `ma_plot(res)` | — | log2fc vs base_mean |
| `pca(E, color="condition", n_top=500)` | `E.pca()` | samples on PC1/2 of log-CPM of the most variable rows |
| `plot_genes(E, ["TP53", "MYC"], by="condition")` | `E.plot_genes()` | per-gene strip plots by group |
| `heatmap(E, rows=mask, by="condition", scale="row")` | `E.heatmap()` | z-scored log-CPM of selected rows |

Same look as `signal_draw` (one axis per panel, `ax=` accepted, returns
the Figure).

## 9. Module layout, dependencies, packaging

| File | Contents | Size |
|---|---|---|
| `genomeblocks/expression.py` | `Expression`, readers and sniffers, `to_genes`, normalisation, `de()` dispatch | ~900 lines |
| `genomeblocks/_norm.py` | `tmm_factors`, `size_factors`, `cpm`, `filter_by_expr` (shared with `signal.py`) | ~150 |
| `genomeblocks/backends/de.py` | engine adapters, the result normaliser, the R bridge | ~400 |
| `genomeblocks/backends/_r/{deseq2,edger,limma}.R` | the R side | ~80 each |
| `genomeblocks/expression_draw.py` | §8 | ~250 |
| `genomeblocks/backends/__init__.py` | the `de` family, the `rscript:` probe, per-method `AUTO` | +40 |
| `genomeblocks/genes.py` | `rows(missing=)`, `select_isoforms(expr=)` | +60 |
| `genomeblocks/__init__.py` | `Expression`, `volcano`, `ma_plot` exports | +3 |

Dependencies: nothing new in core (numpy, pandas, pyarrow are there; the
`tables` backend parses). New extra `de = ["pydeseq2>=0.5"]` (pydeseq2 pulls
formulaic for the design formula); `all` and `test` include it. R is never a
Python dependency: `backend="r"` finds `Rscript` at run time. `environment.yml`
gains `pydeseq2`, `r-base`, `bioconductor-deseq2`, `bioconductor-edger`,
`bioconductor-limma`, `bioconductor-tximport` so the dev env can run the
parity tests; the bioconda recipe lists `pydeseq2` under `test` only.
Package data: `backends/_r/*.R` added to `[tool.setuptools.package-data]`.

Version: 2.1.0 (additive, no breaking change).

## 10. Tests and benchmarks

Fixtures (synthetic, written to `tmp_path` like the rest of `conftest.py`):
six salmon directories (3 WT, 3 KO) over the four transcripts of the test
GTF, with version suffixes and one GENCODE pipe header; the same numbers as
kallisto, RSEM, featureCounts (one file), STAR, HTSeq, an nf-core matrix,
and a `samples.tsv`.

| Test | Checks |
|---|---|
| readers | every format yields identical layers for the same numbers; polars and pandas `tables` backends agree (`test_parsers_agree` pattern); sample naming rules; sniff errors |
| alignment | exact / version-stripped / pipe-header matching; unmatched and missing counts; `strict`; `align` after the fact; `loci=` on a peak matrix |
| tximport | `to_genes` in the three modes equals numbers pinned once from `tximport` on the fixture |
| normalisation | TMM factors and size factors equal numbers pinned from edgeR / DESeq2 on a 20 × 6 matrix; `filter_by_expr` equals `filterByExpr` on the same |
| DE, statistical | a seeded negative-binomial simulation (500 genes, 20 induced) — pydeseq2 recovers ≥ 18 at padj < 0.05 with ≤ 2 false positives; the R engines the same, skipped when `Rscript` or the package is absent (as graph-tool / cgranges tests are) |
| DE, parity | where both run: pydeseq2 vs R DESeq2 log2fc Pearson > 0.99, padj < 0.05 sets Jaccard > 0.8 (documented as *close*, not identical) |
| schema | every engine returns the §6 columns, index = row number, NaN for untested rows, `attrs` filled |
| round trips | `save` / `load`, `to_anndata` / `from_anndata`, `to_pandas` / `from_frame`, protocols |
| isoforms | `select_isoforms(expr=)` alone and with peaks reproduces the bigWig path's bookkeeping |

CI: the `test,fast,interop,bam,viz` matrix entry adds `de`; R tests skip
there. A `workflow_dispatch` job with `r-lib/actions/setup-r` and the three
Bioconductor packages runs the parity tests on demand (installing them takes
several minutes, too slow for every push).

Benchmarks: `benchmarks/bench_expression.py` — read 6 × GENCODE-size
`quant.sf` (polars vs pandas), `to_genes`, TMM, `de` on pydeseq2 vs R, on
the synthetic data `make_data.py` already produces; results into
`benchmarks/results/expression.json` and the benchmarks page.

## 11. Documentation and example

- `docs/guide/expression.md` (read → align → summarise → normalise → DE →
  back to CREs), `docs/api/expression.md`, `docs/design/expression.md` with
  three diagrams (rows aligned to gene rows; tximport over the link column;
  method / engine dispatch), built by `docs/_tools/build_diagrams.py`.
- Edits: `concepts.md` (the table list and the boundary table gain a row),
  `backends.md` ("the seven families"), `api/index.md` module table,
  `guide/index.md`, `installation.md` (the `de` extra and R), `README.md`
  (one bullet), `release-notes.md` (2.1.0).
- `examples/rnaseq_de/`: one notebook, one analysis — RNA-seq of the same
  LNCaP ± DHT system as `ar_foxa1_lncap` (quantified with nf-core/rnaseq →
  salmon), `Expression` → DESeq2 → induced genes' promoters ∩ AR / FOXA1
  CREs → a `View`. Timing and `tqdm` as in the other examples
  (`feedback_example_code_style`). Accession to be chosen (§15).

## 12. Phases

| Phase | Delivers | PR |
|---|---|---|
| **1 — the table** | `Expression`, every reader in §3, alignment, `to_genes`, normalisation (`_norm.py` shared with `signal.tmm`), protocols, save / load, `lookup`, `Genes.rows(missing=)`, tests, guide + API pages | PR 1 |
| **2 — differential expression** | the `de` family: pydeseq2 and the R bridge (DESeq2, edgeR-QL, limma-voom), design / contrasts / filtering, the result schema, `volcano` / `ma_plot` / `pca` / `plot_genes` / `heatmap`, tests with the R skips, backends page | PR 2 |
| **3 — integration and release** | `select_isoforms(expr=)`, `View` / `browser` colouring, `inmoose` engine, numeric contrasts and offsets, the example notebook, benchmarks, design page, release notes, 2.1.0 tag | PR 3 |

Phase 1 is useful on its own (an aligned, normalised table that joins to
CREs) and has no new dependency, so it can merge before any DE engine
exists. Each PR carries its docs. I merge nothing: PRs are opened for you to
review and merge (`reference_github_access`).

## 13. Non-goals

- **Single-cell matrices.** AnnData / scanpy own that; `from_anndata` is for
  pseudobulk or bulk stored as AnnData.
- **Batch correction** (ComBat-seq, RUV): pass `batch` in the design instead;
  a corrected matrix from elsewhere can be loaded with `from_frame`.
- **Gene-set enrichment**: the result frame goes straight into decoupler /
  gseapy; not reimplemented.
- **Differential transcript usage** (DRIMSeq, DEXSeq, sleuth) and salmon /
  kallisto inferential replicates. Transcript-level DE with DESeq2 on `E_tx`
  works, DTU does not.
- **rpy2.** An in-process R would pin genomeblocks to one R installation; the
  subprocess bridge takes any `Rscript`, which matters here (§14).

## 14. What I checked locally

- genomeblocks 2.0.1 on `main` (`afd9c41`), clean tree; branch `expression`
  created from it.
- The test env (`~/micromamba/envs/notebook_cn1`) has polars 1.44, pandas
  2.3, anndata 0.13, statsmodels, scipy, no pydeseq2 / rpy2 / inmoose.
- R: the system R 4.1 has none of the packages; `notebook_r` (R 4.4) has
  DESeq2, tximport, apeglm, ashr but not edgeR / limma; `hichip` has edgeR
  and limma but not DESeq2. That split is why `backend="r"` takes an
  explicit `rscript=` and why the method support of one `Rscript` is probed
  per package.
- No salmon / featureCounts outputs under the usual data directories, so the
  fixtures are synthetic and the example needs a download step.
- pydeseq2's latest release is 0.5.4 (2026-01-23); its `DeseqDataSet` takes
  samples × genes counts, a metadata frame and a `design` formula. inmoose
  (0.9.x) ports DESeq2, a partial edgeR (`glmFit` / `glmLRT`; no
  quasi-likelihood fit documented) and limma; I could not confirm `voom`
  in it, hence phase 3 and "edger (LRT only)" above.

## 15. Decisions needed

My recommendation first in each.

1. **Name.** `Expression` (what you asked for; what people search for) —
   or `Counts` / `Assay`, which read better for the CRE-count use in §7.
2. **Rows strictly aligned to the `Genes` / `Loci` table** (NaN for
   unquantified rows; unmatched input features dropped and reported) — or
   an own feature table with a `gene` row-link column like
   `transcripts['gene']`. Strict alignment is what makes every join array
   indexing; the link column keeps every input row but puts a lookup back
   into every use.
3. **Layout features × samples** (signal cube, SummarizedExperiment) — or
   samples × features (AnnData). Conversion is a transpose either way.
4. **`method` defaults to `"deseq2"`** (the pure-pip path) — or is
   required, so a user always names the model. Both obey the backend
   rules; this is about ergonomics.
5. **R through an `Rscript` subprocess with bundled scripts** — or rpy2.
6. **Default summarisation `length_scaled_tpm`** (no offsets needed, every
   engine accepts it) — or `counts` + offsets (closer to the raw model, but
   pydeseq2 cannot take the offset matrix).
7. **All eight input formats in phase 1** (each reader is ~15 lines once
   the sniffer exists) — or salmon + featureCounts + matrix first.
8. **The example dataset**: LNCaP ± DHT RNA-seq to pair with the existing
   example (I would look for an nf-core-processed or GEO-deposited salmon
   run), or something you already have quantified.
9. **Three PRs as in §12** — or one PR per phase 1 + 2 together.
10. **`Genes` stays link-free** (no `genes.expression` back-reference) —
    or `Genes.attach(E)` so `genes["TP53"].expression` works. The former
    keeps the no-hidden-state rule; the latter is the nicer notebook line.
