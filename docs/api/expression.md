---
title: expression
parent: API Reference
layout: default
nav_order: 3.5
---

# `genomeblocks.expression`
{: .no_toc }

A features × samples matrix whose rows are the rows of a
[`Genes`]({{ '/api/genes/' | relative_url }}) table (`genes.genes` or
`genes.transcripts`) or of any [`Loci`]({{ '/api/loci/' | relative_url }}).
Layers are float64 arrays of shape `(n, m)`; `counts` is the default layer
when present. With `genes=` the rows line up exactly: an input feature that
matches no row is dropped and counted, a row the input did not quantify is
NaN. Transcript tables summarise to genes as tximport does; normalisation
follows edgeR and DESeq2. See the [Expression guide]({{ '/guide/expression/' | relative_url }})
for the workflow.
{: .fs-5 .fw-300 }

```python
import numpy as np
import genomeblocks as gb
from genomeblocks import Expression
from genomeblocks.expression import sniff_format, sample_name, read_table, LEVELS, HOW
from genomeblocks._norm import lib_sizes, cpm, tmm_factors, size_factors, filter_by_expr
```

The examples use the three-gene GTF of the [Genes API]({{ '/api/genes/' | relative_url }})
page and six salmon runs of it (`WT_1`–`WT_3`, `KO_1`–`KO_3`), passed as a
list so the samples keep that order.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `Expression`
{: .sec-navy }

```python
Expression(layers, samples, *, ids, names=None, rows=None, level=None, filename=None, meta=None, unmatched=None)
```

`layers` is a dict of `(n, m)` arrays (at least one; normally `counts`),
`samples` a DataFrame indexed by sample name or a list of names, `ids` one
id per row. `rows` / `level` are set by `make` / `align`, not by hand.

| Attribute | Type | One line |
|---|---|---|
| `layers` | `dict[str, ndarray]` | Every layer, `(n, m)` float64. |
| `default_layer` | `str` | `'counts'` when present, else the first layer. |
| `X` | `ndarray` | The default layer. |
| `counts`, `tpm`, `length` | `ndarray` | The layers by name; `tpm` is computed from `counts / length` when the quantifier gave lengths but no TPM; `counts` raises `AttributeError` on a TPM-only table. |
| `ids`, `names` | object arrays | Feature ids (`gene_id`, `transcript_id` or uid once aligned) and gene names (`None` when unknown). |
| `samples` | `DataFrame` | Sample metadata indexed by sample name, in column order; assignable (a permutation of the names is put back in column order). |
| `sample_names`, `n_samples` | `list`, `int` | The index of `samples` and its length. |
| `rows`, `table`, `level` | `Genes` / `Loci` / `None`, `Loci`, `str` | What the rows are aligned to, its row table, and `'gene'` / `'transcript'` / `'loci'` (`None` before alignment when unknown). |
| `missing` | bool array | Rows that are NaN in every layer. |
| `unmatched` | object array | Input ids that matched no row. |
| `meta` | `dict` | `format`, `key_rule`, `rows`, `counts_from`, … |
| `lib_size` | `ndarray` | Column totals of `counts` over the quantified rows. |
| `shape`, `columns`, `len(E)` | | `(n, m)`, the sample names, `n`. |

```python
genes = gb.Genes.make("genes.gtf")
files = [f"salmon/{s}/quant.sf" for s in ("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")]
E = Expression.make(files, genes=genes, samples="samples.tsv", verbose=False)
E
# -> Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes)
E.shape, E.columns, E.default_layer, E.level, E.table is genes.genes
# -> ((3, 6), ['WT_1', 'WT_2', 'WT_3', 'KO_1', 'KO_2', 'KO_3'], 'counts', 'gene', True)
E.ids, E.names, E.missing, E.unmatched
# -> (array(['G1', 'G2', 'G3'], dtype=object), array(['GENE_A', 'GENE_B', 'GENE_C'], dtype=object),
#     array([False, False, False]), array(['ERCC-00002'], dtype=object))
E.meta
# -> {'format': 'salmon', 'sources': 6, 'key_rule': 'transcript_id, pipe header cut + version stripped', 'rows': 'genes.genes', 'counts_from': 'length_scaled_tpm'}
Expression({"counts": np.array([[1.0, 2.0], [3.0, 4.0]])}, ["a", "b"], ids=["x", "y"])
# -> Expression(2 features × 2 samples, layers=[counts])
Expression({}, ["a"], ids=["x"])
# -> ValueError: Expression: at least one layer (normally 'counts') is required
```

### Construction

| Method | Signature | One line |
|---|---|---|
| `Expression.make` | `make(source, *, genes=None, loci=None, level=None, how="length_scaled_tpm", samples=None, key="auto", sample_names=None, strand="auto", strict=False, verbose=True, backend=None)` | Read quantifications — a glob / list of per-sample files (salmon, kallisto, RSEM, STAR, HTSeq), a `{sample: path}` dict, or one or more multi-sample matrices (featureCounts, nf-core, any id + samples table) — and align them. |
| `Expression.from_frame` | `from_frame(df, *, ids=None, names=None, layer="counts", genes=None, loci=None, level=None, how=..., samples=None, key="auto", strict=False, verbose=True)` | From any table: `ids` the id column (default the first), `names` a gene-name column, numeric columns are samples, `layer` names what the numbers are. |
| `Expression.from_anndata` | `from_anndata(adata, *, layer=None, genes=None, loci=None, level=None, how=..., key="auto", strict=False, verbose=True)` | From samples (obs) × features (var): `X` or `layer` becomes `counts`, other layers come along, `obs` becomes `samples`. |
| `Expression.load` | `load(path, *, genes=None, loci=None, key="auto", verbose=False)` | Read a directory written by [`save`](#export-and-persistence); `genes=` / `loci=` re-align by id. |

`source` formats are sniffed from the header: salmon `quant.sf` /
`quant.genes.sf`, kallisto `abundance.tsv`, RSEM `.isoforms.results` /
`.genes.results`, featureCounts, STAR `ReadsPerGene.out.tab` (`strand`:
`"auto"`, `"none"`, `"forward"`, `"reverse"`), HTSeq-count, nf-core matrices
(layer from the file name), plain TSV / CSV / parquet matrices. `level` is
the level of the result (`"gene"` by default with `genes=`; a transcript
input is summarised with `how`). `samples` is a DataFrame or TSV / CSV path
joined on the sample names. `backend` picks the tables parser
(`'polars'`, `'pandas'`; identical numbers).

```python
Expression.make("counts.featureCounts.txt", genes=genes)
# -> [INFO] Expression: 3 features × 6 samples from featurecounts (counts.featureCounts.txt); layers counts, length
# -> [INFO] Aligned to genes.genes by gene_id: 3/3 input ids matched; 0 of 3 rows not quantified
# -> Expression(3 genes × 6 samples, layers=[counts, length], aligned to genes.genes)
Expression.make("salmon/*/quant.sf", verbose=False)                   # no genes=: file order, a glob sorts
# -> Expression(5 transcripts × 6 samples, layers=[counts, tpm, length])
Expression.make(["nfcore/salmon.merged.gene_counts.tsv", "nfcore/salmon.merged.gene_tpm.tsv"], verbose=False).layers.keys()
# -> dict_keys(['counts', 'tpm'])
Expression.make({"a": "salmon/WT_1/quant.sf", "b": "salmon/KO_1/quant.sf"}, verbose=False).sample_names
# -> ['a', 'b']
Expression.from_frame({"symbol": ["GENE_A", "GENE_B"], "a": [1, 2], "b": [3, 4]}, ids="symbol", layer="tpm", verbose=False)
# -> Expression(2 features × 2 samples, layers=[tpm])
Expression.from_anndata(E.to_anndata(), genes=genes, verbose=False)
# -> Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes)
```

### Alignment and summarisation

| Method | Signature | One line |
|---|---|---|
| `align` | `align(genes=None, loci=None, *, level=None, how="length_scaled_tpm", key="auto", strict=False, verbose=True) -> Expression` | The same numbers on the rows of `genes.genes` / `genes.transcripts` or of `loci`: the best of the id rules (exact, pipe header cut, version stripped, `gene_name`; for a Loci `name`, uid, `chrom:start-end`) is used and recorded in `meta['key_rule']`; unmatched ids dropped (`strict=True` raises), unquantified rows NaN, several inputs on one row summed; a transcript table asked for `level="gene"` is then summarised. |
| `to_genes` | `to_genes(how="length_scaled_tpm", *, verbose=True) -> Expression` | Transcripts to genes over `transcripts['gene']` as `tximport::summarizeToGene`: counts and TPM sum; the gene length is the TPM-weighted mean of its transcripts' effective lengths, filled with the geometric mean over the other samples (or the plain mean when no sample has abundance). `how`: `"counts"`, `"scaled_tpm"`, `"length_scaled_tpm"`. |

```python
U = Expression.make(files, verbose=False)
U.align(genes=genes, level="transcript", verbose=False)
# -> Expression(4 transcripts × 6 samples, layers=[counts, tpm, length], aligned to genes.transcripts)
T = Expression.make(files, genes=genes, level="transcript", verbose=False)
T.to_genes(how="counts", verbose=False).counts
# -> array([[150., 160., 170., 300., 310., 275.],
#           [300., 280., 310.,   0.,   0., 160.],
#           [  0.,   0.,   0.,   0.,   0.,   0.]])
T.to_genes(how="counts", verbose=False).length.round(3)
# -> array([[1872.408, 1807.398, 1900.234, 2234.001, 2293.408, 2229.889],
#           [ 620.   ,  625.   ,  630.   ,  629.931,  629.931,  645.   ],
#           [9632.5  , 9632.5  , 9632.5  , 9632.5  , 9632.5  , 9632.5  ]])
T.to_genes(verbose=False).counts.round(3)
# -> array([[157.866, 172.517, 178.762, 300.   , 310.   , 269.14 ],
#           [292.134, 267.483, 301.238,   0.   ,   0.   , 165.86 ],
#           [  0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ]])
Expression.make(files, genes=genes, strict=True, verbose=False)
# -> ValueError: 1 of 5 input ids match no transcript row (first: ['ERCC-00002']); strict=False keeps going
cre = gb.Loci.make("peaks.bed")
Expression.from_frame({"Geneid": ["chr1:900-1100", "chr1_1900_2100", "chr9:1-2"], "a.bam": [12, 3, 1]}, loci=cre, verbose=False)
# -> Expression(7 loci × 1 samples, layers=[counts], aligned to Loci(n=7))
```

### Normalisation

| Method | Signature | Returns |
|---|---|---|
| `cpm` | `cpm(*, log=False, prior=2.0, norm=None) -> ndarray` | Counts per million by library size (edgeR's `cpm`); `norm="tmm"` / `"size"` scales the library sizes by `tmm()` / `size_factors()`; `log=True` adds the library-scaled prior count before log2. |
| `log_cpm` | `log_cpm(*, prior=2.0, norm=None) -> ndarray` | `cpm(log=True, ...)`. |
| `tmm` | `tmm(*, trim_lfc=0.3, trim_mag=0.05, ref=None, lib_size=None) -> ndarray` | edgeR's TMM factors, one per sample, geometric mean 1. |
| `size_factors` | `size_factors() -> ndarray` | DESeq2's median-of-ratios size factors. |
| `filter_by_expr` | `filter_by_expr(group=None, *, lib_size=None, min_count=10, min_total_count=15, large_n=10, min_prop=0.7) -> bool array` | edgeR's `filterByExpr` mask over the rows; `group` is a `samples` column name or an array. |

NaN rows take no part and stay NaN. The test suite pins every rule against
edgeR 4.4 / DESeq2 1.46 on a reference matrix.

```python
E.lib_size, E.tmm().round(3), E.size_factors().round(3)
# -> (array([450., 440., 480., 300., 310., 435.]), array([0.711, 0.697, 0.703, 2.025, 2.025, 0.701]),
#     array([0.709, 0.775, 0.803, 1.347, 1.392, 1.209]))
E.cpm().round(1)[0], E.log_cpm(norm="tmm").round(2)[1]
# -> (array([ 350813.4,  392083. ,  372421.3, 1000000. , 1000000. ,  618711.6]),
#     array([19.8 , 19.73, 19.76, 12.21, 12.21, 19.05]))
E.filter_by_expr("condition"), E.filter_by_expr()
# -> (array([ True,  True, False]), array([ True, False, False]))
```

### Lookups and subsets

| Name | Returns | One line |
|---|---|---|
| `E['TP53']`, `E['ENSG...']`, `E[uid]` | `Series` | One feature over the samples, from the default layer (gene names, ids and Loci uids are keys). |
| `E[['TP53', 'MYC']]` | `DataFrame` | Several features × samples. |
| `E[i]` | `Series` | Row `i`. |
| `E[mask]`, `E[rows]`, `E[a:b]`, `E.take(idx)` | `Expression` | A row subset; the alignment is dropped (ids and names stay). |
| `E[:, mask]`, `E[:, names]`, `E[:, a:b]`, `E.take_samples(idx)` | `Expression` | A sample subset (mask, names, rows or slice); the alignment is kept. |
| `E.get(key, layer=None)` | `Series` | One feature from `layer`. |
| `E.row(key)` | `int` | The row of a key (`KeyError` otherwise). |
| `E.lookup(keys, values=None, *, layer=None, agg="mean", missing=nan)` | `ndarray` | `values` (one per row; default the `agg` of `layer` over the samples) picked for each key; empty / unknown keys give `missing`. |
| `key in E`, `for id in E` | `bool`, `str` | Membership; iterate the ids. |

```python
E["GENE_A"].round(2).tolist(), E[["GENE_A", "G2"]].shape, E.get("GENE_B", "tpm").round(1).tolist()
# -> ([157.87, 172.52, 178.76, 300.0, 310.0, 269.14], (2, 6), [824756.8, 793282.9, 822351.9, 0.0, 0.0, 623464.7])
E.lookup(["GENE_C", "", "nope", "GENE_A"]), E.lookup(["GENE_B", "GENE_A"], E.tpm.mean(axis=1))
# -> (array([  0.        ,          nan,          nan, 231.38071973]), array([510642.7322565 , 411571.04314533]))
E[:, E.samples.condition == "KO"]
# -> Expression(3 genes × 3 samples, layers=[counts, tpm, length], aligned to genes.genes)
E[E.counts.mean(axis=1) > 100]
# -> Expression(2 genes × 6 samples, layers=[counts, tpm, length])
E.row("G3"), "GENE_B" in E, list(E)
# -> (2, True, ['G1', 'G2', 'G3'])
```

### Export and persistence

| Method | Signature | One line |
|---|---|---|
| `to_pandas` | `to_pandas(layer=None, *, long=False) -> DataFrame` | The wide table of `layer` (default layer): id column, `gene_name`, one column per sample; `long=True` one row per (feature, sample) with the sample metadata joined on. |
| `to_polars`, `to_arrow` | `to_polars(layer=None)`, `to_arrow(layer=None)` | The wide table through Arrow; the Arrow / dataframe / narwhals protocols read the default layer. |
| `to_anndata` | `to_anndata(layer=None) -> AnnData` | Samples (obs) × features (var): `X` is `layer`, the other layers kept, `obs` is `samples`, `var` has `gene_name` and `missing`; `uns['genomeblocks']` remembers which layer `X` was, so `from_anndata` is exact. |
| `save` | `save(path)` | `path/` gets `samples.parquet`, `features.parquet` (id, gene_name, missing), one `<layer>.parquet` per layer and `meta.json`. |
| `head`, `tail` | `head(n=5)`, `tail(n=5) -> DataFrame` | The wide table as pandas. |
| `describe` | `describe() -> DataFrame` | Rows, level, samples, layers, alignment, quantified / missing rows, dropped ids, the id rule, library sizes. |

```python
E.to_pandas().round(2)
# ->   gene_id gene_name    WT_1    WT_2    WT_3   KO_1   KO_2    KO_3
#    0      G1    GENE_A  157.87  172.52  178.76  300.0  310.0  269.14
#    1      G2    GENE_B  292.13  267.48  301.24    0.0    0.0  165.86
#    2      G3    GENE_C    0.00    0.00    0.00    0.0    0.0    0.00
E.to_pandas(long=True).columns.tolist()
# -> ['gene_id', 'gene_name', 'sample', 'counts', 'condition', 'batch']
import polars as pl
pl.DataFrame(E).shape, E.to_arrow().schema.names[:3]
# -> ((3, 8), ['gene_id', 'gene_name', 'WT_1'])
E.to_anndata()
# -> AnnData object with n_obs × n_vars = 6 × 3
#        obs: 'condition', 'batch'
#        var: 'gene_name', 'missing'
#        uns: 'genomeblocks'
#        layers: 'tpm', 'length'
import os
E.save("expr")
sorted(os.listdir("expr"))
# -> ['counts.parquet', 'features.parquet', 'length.parquet', 'meta.json', 'samples.parquet', 'tpm.parquet']
Expression.load("expr"), Expression.load("expr", genes=genes)
# -> (Expression(3 genes × 6 samples, layers=[counts, tpm, length]),
#     Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes))
E.describe().loc[["level", "aligned to", "counts from"], "value"].tolist()
# -> ['gene', 'genes.genes', 'length_scaled_tpm']
```

---

## Module helpers
{: .sec-navy }

| Name | Signature | One line |
|---|---|---|
| `sniff_format` | `sniff_format(path) -> dict` | `{'format', 'header', 'skip', 'sep', 'level'}` of a quantification file, from its header line. |
| `sample_name` | `sample_name(path) -> str` | The sample a per-sample file belongs to: the directory for salmon / kallisto layouts, else the file name without its format suffix. |
| `read_table` | `read_table(path, *, strand="auto", sample=None, backend=None)` | One file as an internal table (ids, layers, samples); `Expression.make` is the API. |
| `LEVELS` | `('gene', 'transcript', 'loci')` | The row levels. |
| `HOW` | `('counts', 'scaled_tpm', 'length_scaled_tpm')` | The summarisation modes. |

```python
sniff_format("salmon/WT_1/quant.sf")
# -> {'format': 'salmon', 'header': True, 'skip': 0, 'sep': '\t', 'level': 'transcript'}
sniff_format("counts.featureCounts.txt")
# -> {'format': 'featurecounts', 'header': True, 'skip': 1, 'sep': '\t', 'level': None}
sample_name("salmon/WT_1/quant.sf"), sample_name("star/KO_3.ReadsPerGene.out.tab"), sample_name("rsem/WT_2.isoforms.results")
# -> ('WT_1', 'KO_3', 'WT_2')
```

---

## `genomeblocks._norm`
{: .sec-navy }

The normalisation rules on plain `(features, samples)` arrays — what the
`Expression` methods call, and what `signal.tmm` uses for signal cubes.

| Function | Signature | One line |
|---|---|---|
| `lib_sizes` | `lib_sizes(counts) -> ndarray` | Column totals, ignoring NaN rows. |
| `cpm` | `cpm(counts, *, lib_size=None, factors=None, log=False, prior=2.0) -> ndarray` | edgeR's `cpm`; `factors` scale the library sizes. |
| `tmm_factors` | `tmm_factors(counts, *, trim_lfc=0.3, trim_mag=0.05, ref=None, lib_size=None) -> ndarray` | edgeR's `calcNormFactors(method="TMM")`: reference = the sample whose upper quartile of counts per library size is closest to the mean; rank trimming; asymptotic-variance weights. |
| `size_factors` | `size_factors(counts) -> ndarray` | DESeq2's `estimateSizeFactorsForMatrix`. |
| `filter_by_expr` | `filter_by_expr(counts, *, group=None, lib_size=None, min_count=10, min_total_count=15, large_n=10, min_prop=0.7) -> bool array` | edgeR's `filterByExpr`. |

```python
x = np.array([[10., 20.], [5., 10.], [1., 2.]])
lib_sizes(x), tmm_factors(x), size_factors(x).round(4)
# -> (array([16., 32.]), array([1., 1.]), array([0.7071, 1.4142]))
cpm(x)[0], cpm(x, log=True).round(3)[0], filter_by_expr(x)
# -> (array([625000., 625000.]), array([19.212, 19.212]), array([ True, False, False]))
```

---

## Changes elsewhere
{: .sec-navy }

`Genes.rows(names, *, missing=None)` gained `missing=`: an unknown name
raises `KeyError` as before, or becomes that value when given.

```python
genes.rows(["GENE_A", "nope"], missing=-1)
# -> array([ 0, -1])
```
