---
title: Expression
parent: User Guide
layout: default
nav_order: 3.5
---

# Expression
{: .no_toc }

RNA-seq quantifications as one features × samples table whose rows *are* the
rows of your `Genes` tables — read from salmon, kallisto, RSEM, featureCounts,
STAR, HTSeq, nf-core matrices or any frame, summarised from transcripts to genes
exactly as tximport does, normalised as edgeR and DESeq2 do, and joined back to
CREs by row number.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is
{: .sec-navy }

An `Expression` holds one or more **layers** — float64 arrays of shape
`(features, samples)` — plus a sample table. Built with `genes=`, its rows are
exactly the rows of `genes.genes` (or `genes.transcripts`), so `E.counts[i]`
is gene row `i`, a mask over genes selects the same rows in `genes.get_tss()`,
and nothing is ever merged on an id again.

| Attribute | What it holds |
|---|---|
| `E.counts`, `E.tpm`, `E.length`, `E.layers[name]` | `(n, m)` arrays: `counts` always; `tpm` and `length` (effective length per feature per sample) when the quantifier gives them |
| `E.samples` | a pandas DataFrame indexed by sample name, in column order — your metadata lives here |
| `E.rows`, `E.table`, `E.level` | the `Genes` or `Loci` the rows line up with, its row table, and `'gene'` / `'transcript'` / `'loci'` |
| `E.ids`, `E.names` | the table's ids (`gene_id`, `transcript_id` or uid) and gene names, one per row |
| `E.missing`, `E.unmatched` | rows the input did not quantify (NaN in every layer); input ids that matched no row (dropped, counted) |

```python
import genomeblocks as gb

genes = gb.Genes.make("gencode.v44.annotation.gtf")
E = gb.Expression.make([f"salmon/{s}/quant.sf" for s in ("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")],
                       genes=genes, samples="samples.tsv")
# -> [INFO] Expression: 5 transcripts × 6 samples from 6 salmon files; layers counts, tpm, length
# -> [INFO] Aligned to genes.transcripts by transcript_id, pipe header cut + version stripped: 4/5 input ids matched (1 dropped); 0 of 4 rows not quantified
# -> [INFO] Summarised 4 transcripts to 3 genes (length_scaled_tpm)
E
# -> Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes)
E.shape, list(E.layers), E.table is genes.genes
# -> ((3, 6), ['counts', 'tpm', 'length'], True)
E.ids, E.names
# -> (array(['G1', 'G2', 'G3'], dtype=object), array(['GENE_A', 'GENE_B', 'GENE_C'], dtype=object))
E.counts.round(2)
# -> array([[157.87, 172.52, 178.76, 300.  , 310.  , 269.14],
#           [292.13, 267.48, 301.24,   0.  ,   0.  , 165.86],
#           [  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ]])
E.unmatched, E.missing
# -> (array(['ERCC-00002'], dtype=object), array([False, False, False]))
E.samples
# ->        condition batch
#    sample
#    WT_1          WT    b1
#    WT_2          WT    b2
#    WT_3          WT    b3
#    KO_1          KO    b1
#    KO_2          KO    b2
#    KO_3          KO    b3
```

The examples on this page use the three-gene GTF of the [Genes
guide]({{ '/guide/genes/' | relative_url }}) (`GENE_A` with isoforms `T1` and
`T1b`, `GENE_B`, `GENE_C`) and six salmon runs of it in which `GENE_B` has no
reads in two KO samples and `GENE_C` none anywhere — the cases the gene-length
rule below has to handle. A spike-in, `ERCC-00002`, is in the quantification
but not in the GTF.

**Zero is a measurement, NaN is not.** `GENE_C` is quantified as 0 everywhere.
A gene the input never mentions (see the gene-name CSV below) is NaN in every
layer and `True` in `E.missing`; library sizes and normalisation skip it.

{: .note }
> Differential expression — `E.de(...)` and the `de` backend family (DESeq2,
> edgeR, limma) — is the next release. This one reads, aligns, summarises and
> normalises.

---

## Reading quantifications: `Expression.make`
{: .sec-navy }

```python
E = gb.Expression.make(source, *, genes=None, loci=None, level=None, how="length_scaled_tpm",
                       samples=None, key="auto", sample_names=None, strand="auto", strict=False,
                       verbose=True, backend=None)
```

`source` is a glob or list of per-sample files, a `{sample: path}` dict, or
one or more multi-sample matrices. The format is sniffed from the header
line; the `tables` backend (polars, else pandas) parses, with identical
numbers either way.

| Format | File | Layers | Notes |
|---|---|---|---|
| salmon | `quant.sf`, `quant.genes.sf` | counts, tpm, length | sample = the directory name; `quant.genes.sf` is gene level |
| kallisto | `abundance.tsv` | counts, tpm, length | sample = the directory name |
| RSEM | `*.isoforms.results`, `*.genes.results` | counts, tpm, length | level from the header |
| featureCounts | one file, all samples | counts, length | BAM paths become sample names (`bam/WT_1.sorted.bam` → `WT_1`); the `Chr/Start/End` columns are ignored |
| STAR | `*ReadsPerGene.out.tab` | counts | `strand="auto"` picks the stranded column with far more reads and says so; `"none"`, `"forward"`, `"reverse"` to choose |
| HTSeq-count | two columns | counts | the `__no_feature …` rows are dropped and reported |
| nf-core/rnaseq | `salmon.merged.gene_counts.tsv`, `…gene_tpm.tsv`, `…transcript_counts.tsv` | from the file name | pass several files of one run together to fill several layers |
| any matrix | TSV / CSV / parquet: an id column, then numeric sample columns | counts (or from the name) | a `gene_name` / `symbol` column is kept as the names, other text columns are dropped and reported |

```python
F = gb.Expression.make("counts.featureCounts.txt", genes=genes)
# -> [INFO] Expression: 3 features × 6 samples from featurecounts (counts.featureCounts.txt); layers counts, length
# -> [INFO] Aligned to genes.genes by gene_id: 3/3 input ids matched; 0 of 3 rows not quantified
F.sample_names, F.length[:, 0]
# -> (['WT_1', 'WT_2', 'WT_3', 'KO_1', 'KO_2', 'KO_3'], array([1800.,  800., 9800.]))

gb.Expression.make("star/*.ReadsPerGene.out.tab", genes=genes)
# -> [INFO] Expression: 3 features × 6 samples from 6 star files; layers counts
# -> [INFO]   STAR strand=auto: reverse-stranded counts (column 4; the other strand has far fewer reads) — 6 files
# -> [INFO] Aligned to genes.genes by gene_id: 3/3 input ids matched; 0 of 3 rows not quantified

N = gb.Expression.make(["nfcore/salmon.merged.gene_counts.tsv", "nfcore/salmon.merged.gene_tpm.tsv"], genes=genes)
# -> [INFO] Expression: 3 genes × 6 samples from 2 matrix files; layers counts, tpm
# -> [INFO] Aligned to genes.genes by gene_id: 3/3 input ids matched; 0 of 3 rows not quantified
```

A TPM-only matrix (nf-core's `gene_tpm.tsv` alone, a GTEx table) is a valid
table whose default layer is `tpm`; `E.counts` on it raises and says what
needs counts.

**Samples.** A glob sorts alphabetically, so KO samples come before WT ones;
a list keeps its order and a dict names the samples. `sample_names=` renames
them after the fact. `samples=` is a DataFrame or a TSV / CSV whose first
column is the sample name; it is joined to the table's samples (a missing one
raises) and lives in `E.samples`, where you can add columns at any time.

```python
gb.Expression.make("salmon/*/quant.sf", verbose=False).sample_names
# -> ['KO_1', 'KO_2', 'KO_3', 'WT_1', 'WT_2', 'WT_3']
gb.Expression.make({"a": "salmon/WT_1/quant.sf", "b": "salmon/KO_1/quant.sf"}, verbose=False).sample_names
# -> ['a', 'b']
E.samples["batch"] = ["b1", "b2", "b3", "b1", "b2", "b3"]
```

Without `genes=` / `loci=` the rows are the file's features in file order and
`E.rows` is `None`; [`align`](#lining-up-with-the-gene-tables) lines them up
later.

```python
U = gb.Expression.make("salmon/*/quant.sf", verbose=False)
U, U.rows, U.level
# -> (Expression(5 transcripts × 6 samples, layers=[counts, tpm, length]), None, 'transcript')
```

---

## From a frame or an AnnData
{: .sec-navy }

`Expression.from_frame(df, *, ids=None, names=None, layer="counts", genes=None, loci=None, ...)`
takes any table (pandas, polars, pyarrow, a dict of columns): `ids` names
the feature-id column (default: the first), `names` a gene-name column (found
by name otherwise), the numeric columns are the samples and `layer` says
what the numbers are.

```python
import polars as pl
df = pl.DataFrame({"gene_id": ["G1.5", "G2.1", "G9"], "gene_name": ["GENE_A", "GENE_B", "X"],
                   "s1": [10.0, 20.0, 1.0], "s2": [12.0, 18.0, 2.0]})
P = gb.Expression.from_frame(df, genes=genes)
# -> [INFO] Expression: 3 features × 2 samples from a frame; layers counts
# -> [INFO] Aligned to genes.genes by gene_id, version stripped: 2/3 input ids matched (1 dropped); 1 of 3 rows not quantified
P.counts, P.unmatched
# -> (array([[10., 12.],
#            [20., 18.],
#            [nan, nan]]), array(['G9'], dtype=object))
```

`Expression.from_anndata(adata, *, layer=None, genes=None, ...)` reads an
AnnData of samples (obs) × features (var): `X` (or `layer`) becomes `counts`,
the other layers come along, `obs` becomes `samples`. `to_anndata()` is the
exact inverse (it remembers which layer `X` was).

```python
a = E.to_anndata()
a
# -> AnnData object with n_obs × n_vars = 6 × 3
#        obs: 'condition', 'batch'
#        var: 'gene_name', 'missing'
#        uns: 'genomeblocks'
#        layers: 'tpm', 'length'
gb.Expression.from_anndata(a, genes=genes, verbose=False)
# -> Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes)
```

---

## Lining up with the gene tables
{: .sec-navy }

`genes=` aligns the rows to `genes.genes` (`level="gene"`, the default) or
`genes.transcripts` (`level="transcript"`); a transcript-level input asked
for gene level is aligned to the transcripts and then
[summarised](#transcripts-to-genes-to_genes). Every one of these rules is
tried and the one that matches the most ids wins (the first on a tie), so a
file that mixes versioned and pipe-headed ids still lands on one rule:

| Rule | Matches |
|---|---|
| `gene_id` / `transcript_id` | exactly as written |
| pipe header cut | `ENST…\|ENSG…\|…\|TP53-201\|TP53\|…` (a salmon index built without `--gencode`) cut at the first `\|` |
| version stripped | `ENSG00000141510.17` ↔ `ENSG00000141510`, on both sides; GENCODE's `_PAR_Y` genes stay apart |
| `gene_name` | a symbol-keyed table (gene level only; a duplicated symbol takes its first gene row) |

`key=` forces one (`"gene_id"`, `"transcript_id"`, `"gene_name"`); the rule
used is printed and kept in `E.meta["key_rule"]`. Several input rows that
land on one gene are summed. `strict=True` turns unmatched ids into an error.

```python
B = gb.Expression.make("by_name.csv", genes=genes)         # symbol,biotype,WT_1,...; GENE_C absent
# -> [INFO] Expression: 2 features × 6 samples from matrix (by_name.csv); layers counts
# -> [INFO]   ignored non-numeric column(s) ['biotype']
# -> [INFO] Aligned to genes.genes by gene_name: 2/2 input ids matched; 1 of 3 rows not quantified
B.missing, B.lib_size
# -> (array([False, False,  True]), array([450., 440., 480., 300., 310., 435.]))
gb.Expression.make(files, genes=genes, strict=True, verbose=False)
# -> ValueError: 1 of 5 input ids match no transcript row (first: ['ERCC-00002']); strict=False keeps going
gb.Expression.make("by_name.csv", genes=genes, key="gene_id", verbose=False)
# -> ValueError: none of the 2 input ids matches the gene table (tried 'gene_id'); first ids: ['GENE_A', 'GENE_B'] — pass key= to name the id kind
```

`align(genes=…)` does the same to a table read without `genes=`, and
re-aligns an aligned one to another annotation:

```python
U.align(genes=genes)
# -> [INFO] Aligned to genes.transcripts by transcript_id, pipe header cut + version stripped: 4/5 input ids matched (1 dropped); 0 of 4 rows not quantified
# -> [INFO] Summarised 4 transcripts to 3 genes (length_scaled_tpm)
# -> Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes)
```

**Any `Loci` as the rows.** `loci=` aligns to a CRE table instead — a
featureCounts run over an SAF of peaks gives a count matrix on the CRE rows,
matched by the `name` column, the uid or a `chrom:start-end` region in any
spelling (`key="name"` / `"uid"` / `"region"`):

```python
import pandas as pd
cre = gb.Loci.make("peaks.bed")
pk = pd.DataFrame({"Geneid": ["chr1:900-1100", "chr1_1900_2100", "chr2:500-600", "chr9:1-2"],
                   "a.bam": [12, 3, 7, 1], "b.bam": [30, 2, 9, 1]})
Pk = gb.Expression.from_frame(pk, loci=cre)
# -> [INFO] Expression: 4 features × 2 samples from a frame; layers counts
# -> [INFO] Aligned to Loci(n=7) by region (chrom:start-end): 3/4 input ids matched (1 dropped); 4 of 7 rows not quantified
Pk.counts[:3], Pk.table is cre
# -> (array([[12., 30.],
#            [ 3.,  2.],
#            [nan, nan]]), True)
```

---

## Transcripts to genes: `to_genes`
{: .sec-navy }

A transcript-level table aligned to a `Genes` summarises over
`transcripts['gene']` in a few `np.add.at` calls, and reproduces
`tximport::summarizeToGene` to floating-point precision (the test suite pins
the numbers against tximport, and the [Benchmarks]({{ '/benchmarks/' | relative_url }})
page measures it at GENCODE scale):

| `how` | gene counts | tximport |
|---|---|---|
| `"counts"` | the summed transcript counts | `countsFromAbundance="no"` — needs length offsets in a model |
| `"scaled_tpm"` | gene TPM scaled to the library size | `"scaledTPM"` |
| `"length_scaled_tpm"` (default) | gene TPM × the gene's mean length, scaled to the library size | `"lengthScaledTPM"` — the counts tximport recommends for DESeq2 / edgeR without offsets |

Gene TPM is the sum of transcript TPM in every mode. The gene length is the
TPM-weighted mean of its transcripts' effective lengths; where a gene has no
abundance in a sample it takes the geometric mean of its lengths in the other
samples, and where it has none anywhere the plain mean of its transcripts'
lengths — tximport's `replaceMissingLength`.

```python
T = gb.Expression.make(files, genes=genes, level="transcript", verbose=False)
T.counts
# -> array([[100., 120., 110., 100.,  90.,  95.],
#           [ 50.,  40.,  60., 200., 220., 180.],
#           [300., 280., 310.,   0.,   0., 160.],
#           [  0.,   0.,   0.,   0.,   0.,   0.]])
G = T.to_genes(how="counts")
# -> [INFO] Summarised 4 transcripts to 3 genes (counts)
G.counts
# -> array([[150., 160., 170., 300., 310., 275.],
#           [300., 280., 310.,   0.,   0., 160.],
#           [  0.,   0.,   0.,   0.,   0.,   0.]])
G.length.round(3)                      # GENE_B in KO_1 / KO_2: geometric mean of the rest; GENE_C: plain mean
# -> array([[1872.408, 1807.398, 1900.234, 2234.001, 2293.408, 2229.889],
#           [ 620.   ,  625.   ,  630.   ,  629.931,  629.931,  645.   ],
#           [9632.5  , 9632.5  , 9632.5  , 9632.5  , 9632.5  , 9632.5  ]])
T.to_genes(how="length_scaled_tpm", verbose=False).counts.round(3)
# -> array([[157.866, 172.517, 178.762, 300.   , 310.   , 269.14 ],
#           [292.134, 267.483, 301.238,   0.   ,   0.   , 165.86 ],
#           [  0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ]])
E.meta["counts_from"]                  # what make() did by default
# -> 'length_scaled_tpm'
```

`make(..., how=)` picks the mode when summarising straight from the files;
`level="transcript"` keeps the transcript table (both can be kept).

---

## Normalisation
{: .sec-navy }

Pure numpy, one rule each, and the test suite pins every one against the R
original on a reference matrix (edgeR's `calcNormFactors`, `cpm`,
`filterByExpr`; DESeq2's `estimateSizeFactorsForMatrix`). Rows that are NaN
take no part. Nothing is stored unless you put it in a layer.

| Call | Returns | Rule |
|---|---|---|
| `E.lib_size` | `(m,)` | column totals of `counts` over the quantified rows |
| `E.cpm(log=False, prior=2.0, norm=None, factors=None)` | `(n, m)` | counts per million; `norm="tmm"` / `"size"` scales the library sizes (or `factors=` ones computed once); `log=True` adds edgeR's library-scaled prior count before log2 |
| `E.log_cpm(prior=2.0, norm=None, factors=None)` | `(n, m)` | `cpm(log=True)` |
| `E.tmm(trim_lfc=0.3, trim_mag=0.05, ref=None)` | `(m,)` | edgeR's trimmed mean of M-values, geometric mean 1 |
| `E.size_factors()` | `(m,)` | DESeq2's median of ratios |
| `E.filter_by_expr(group=None, min_count=10, min_total_count=15, large_n=10, min_prop=0.7)` | `(n,)` bool | edgeR's `filterByExpr`; `group` is a column of `E.samples` by name or an array |
| `E.tpm` | `(n, m)` | the layer, or counts / effective length rescaled to a million when the quantifier gave lengths but no TPM (featureCounts) |

```python
E.lib_size
# -> array([450., 440., 480., 300., 310., 435.])
E.tmm().round(3), E.size_factors().round(3)
# -> (array([0.711, 0.697, 0.703, 2.025, 2.025, 0.701]), array([0.709, 0.775, 0.803, 1.347, 1.392, 1.209]))
E.log_cpm(norm="tmm").round(2)
# -> array([[18.91, 19.1 , 19.02, 18.91, 18.91, 19.75],
#           [19.8 , 19.73, 19.76, 12.21, 12.21, 19.05],
#           [12.21, 12.21, 12.21, 12.21, 12.21, 12.21]])
E.filter_by_expr("condition"), E.filter_by_expr()
# -> (array([ True,  True, False]), array([ True, False, False]))
```

(Three genes make for odd-looking factors; the rules are the point.
`F.tpm` on the featureCounts table above is computed from its counts and
lengths.)

---

## Looking things up and subsetting
{: .sec-navy }

A gene name, a gene id or a uid is a key; a layer name picks the values.

```python
E["GENE_A"]                            # a Series over the samples, from the default layer
# -> sample
#    WT_1    157.866049
#    WT_2    172.516506
#    WT_3    178.762212
#    KO_1    300.000000
#    KO_2    310.000000
#    KO_3    269.139551
#    Name: GENE_A, dtype: float64
E[["GENE_A", "G2"]].round(2)           # a frame, features × samples
# -> sample    WT_1    WT_2    WT_3   KO_1   KO_2    KO_3
#    GENE_A  157.87  172.52  178.76  300.0  310.0  269.14
#    G2      292.13  267.48  301.24    0.0    0.0  165.86
E.get("GENE_B", "tpm").round(1).tolist()
# -> [824756.8, 793282.9, 822351.9, 0.0, 0.0, 623464.7]
```

`lookup(keys, values=None, *, layer=None, agg="mean", missing=nan)` picks a
per-row vector (default: `agg` of `layer` over the samples) for an array of
keys, NaN for an empty or unknown one — the one helper the CRE side needs:

```python
E.lookup(["GENE_C", "", "nope", "GENE_A"])
# -> array([  0.        ,          nan,          nan, 231.38071973])
E.lookup(["GENE_B", "GENE_A"], E.tpm.mean(axis=1))
# -> array([510642.7322565 , 411571.04314533])
```

**Subsets.** `E[:, mask]`, `E[:, names]` and `take_samples` keep the row
alignment; `E[mask]`, `E[rows]` and `take` drop it, because the result no
longer lines up with the whole gene table (ids and names stay). Masks over
genes are meant for the TSS table and for results, not for `E`.

```python
K = E[:, E.samples.condition == "KO"]
K, K.rows is genes
# -> (Expression(3 genes × 3 samples, layers=[counts, tpm, length], aligned to genes.genes), True)
S = E[E.counts.mean(axis=1) > 100]
S, S.rows, S.ids
# -> (Expression(2 genes × 6 samples, layers=[counts, tpm, length]), None, array(['G1', 'G2'], dtype=object))
```

---

## Back to the genome
{: .sec-navy }

This is what the alignment buys: every per-gene vector is already in
`genes.genes` order, so it drops onto the TSS table as a column, a mask over
genes selects promoters, and set algebra does the rest.

```python
tss = genes.get_tss()
tss["tpm"] = E.tpm.mean(axis=1)        # one value per gene row — no join
tss.to_pandas()
# ->   chrom  start    end strand gene_name gene_id  gene            tpm
#    0  chr1   1000   1001      +    GENE_A      G1     0  411571.043145
#    1  chr1  10999  11000      -    GENE_B      G2     1  510642.732256
#    2  chr2  20000  20001      +    GENE_C      G3     2       0.000000
on = E.tpm.mean(axis=1) > 1e5
(cre & tss[on].slop(1000)).to_records()               # CREs at promoters of expressed genes
# -> [('chr1', 900, 1100, '.'), ('chr1', 1900, 2100, '.'), ('chr1', 9950, 10050, '.'), ('chr1', 10900, 11100, '.')]
```

After `A.annotate(genes)` every vertex of an
[`Architecture`]({{ '/guide/architecture/' | relative_url }}) carries a gene
name in `A.vp.gene`; `lookup` turns that into a vertex column:

```python
A = (gb.Architecture.make(cre, "loops.bedpe", r=100).add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes))
A.vp.gene
# -> array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_A', ''], dtype=object)
A.vp["tpm"] = E.lookup(A.vp.gene, E.tpm.mean(axis=1))
A.vp["tpm"]
# -> array([411571.04314533, 411571.04314533, 411571.04314533, 510642.7322565 , 510642.7322565 , 411571.04314533, nan])
genes.rows(["GENE_A", "nope"], missing=-1)             # new: unknown names as -1 instead of KeyError
# -> array([ 0, -1])
```

---

## As a table
{: .sec-navy }

`shape` is the matrix shape `(features, samples)` and `columns` the sample
names. As a frame the object is the wide table of its default layer — the
id column, `gene_name`, then one column per sample — through `to_pandas`,
`to_polars`, `to_arrow` and the Arrow / dataframe / narwhals protocols;
`to_pandas(long=True)` is one row per (feature, sample) with the sample
metadata joined on, for seaborn.

```python
E.shape, E.columns, len(E)
# -> ((3, 6), ['WT_1', 'WT_2', 'WT_3', 'KO_1', 'KO_2', 'KO_3'], 3)
E.to_pandas().round(2)
# ->   gene_id gene_name    WT_1    WT_2    WT_3   KO_1   KO_2    KO_3
#    0      G1    GENE_A  157.87  172.52  178.76  300.0  310.0  269.14
#    1      G2    GENE_B  292.13  267.48  301.24    0.0    0.0  165.86
#    2      G3    GENE_C    0.00    0.00    0.00    0.0    0.0    0.00
E.to_pandas(long=True).round(2).head(4)
# ->   gene_id gene_name sample  counts condition batch
#    0      G1    GENE_A   WT_1  157.87        WT    b1
#    1      G1    GENE_A   WT_2  172.52        WT    b2
#    2      G1    GENE_A   WT_3  178.76        WT    b3
#    3      G1    GENE_A   KO_1  300.00        KO    b1
import polars as pl
pl.DataFrame(E).shape, E.to_polars("length").shape
# -> ((3, 8), (3, 8))
E.describe().loc[["aligned to", "id rule", "counts from"], "value"].tolist()
# -> ['genes.genes', 'transcript_id, pipe header cut + version stripped', 'length_scaled_tpm']
```

`save(dir)` writes a directory of parquet files; `load(dir, genes=…)` reads
it back and re-aligns by id, `load(dir)` leaves it unaligned.

```python
import os
E.save("expr")
sorted(os.listdir("expr"))
# -> ['counts.parquet', 'features.parquet', 'length.parquet', 'meta.json', 'samples.parquet', 'tpm.parquet']
gb.Expression.load("expr"), gb.Expression.load("expr", genes=genes)
# -> (Expression(3 genes × 6 samples, layers=[counts, tpm, length]),
#     Expression(3 genes × 6 samples, layers=[counts, tpm, length], aligned to genes.genes))
```

---

## See also
{: .sec-navy }

- [Expression API]({{ '/api/expression/' | relative_url }}) — every method and argument.
- [Genes]({{ '/guide/genes/' | relative_url }}) — the tables the rows line up with; `get_tss`, `annotations`, `select_isoforms`.
- [Architecture]({{ '/guide/architecture/' | relative_url }}) — `annotate` and `vp` columns.
- [Benchmarks]({{ '/benchmarks/' | relative_url }}) — reading, summarising and normalising against tximport, edgeR and DESeq2.
- [Interoperability]({{ '/interoperability/' | relative_url }}) — frames, AnnData and the protocols.
