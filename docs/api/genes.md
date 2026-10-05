---
title: genes
parent: API Reference
layout: default
nav_order: 3
---

# `genomeblocks.genes`
{: .no_toc }

A gene annotation as three linked [`Loci`]({{ '/api/loci/' | relative_url }})
tables — `genes`, `transcripts`, `features` — in the same 0-based half-open
coordinates as every other table (a GTF's 1-based start becomes `start - 1`
once, while parsing). The links are row numbers: `transcripts['gene'][k]`
is the gene row of transcript `k`, `features['transcript'][j]` the
transcript row of exon / CDS / UTR `j`. A gene's TSS is the 1-bp interval
`[t, t + 1)` with `t = start` on `+` and `t = end - 1` on `-`. See the
[Genes guide]({{ '/guide/genes/' | relative_url }}) for the workflow and
[Design: genes]({{ '/design/genes/' | relative_url }}) for the parser.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
from genomeblocks import Genes
from genomeblocks.genes import GeneView, LABELS, KINDS, tss_base
```

The examples use a 20-line GTF with three genes: `GENE_A` (chr1:1000-5000,
`+`, isoforms `T1` with three exons and `T1b` with one), `GENE_B`
(chr1:10000-11000, `-`) and `GENE_C` (chr2:20000-30000, `+`).

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The three tables
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">{% include diagrams/genes-model.svg %}</div><figcaption><strong>Three row-aligned tables, linked by row number.</strong> <code>genes</code> holds one row per gene, <code>transcripts</code> one per isoform with <code>gene</code> pointing at the gene's row, <code>features</code> one per exon / CDS / UTR with <code>transcript</code> pointing at the transcript's row. Every table is a <code>Loci</code> with 0-based starts, so <code>start − 1</code> is applied once while reading the GTF.</figcaption></figure>

| Table | Columns beyond `chrom / start / end / strand` | Notes |
|---|---|---|
| `G.genes` | `gene_id`, `gene_name`, `gene_type` | One row per gene, in file order. `gene_name` falls back to `gene_id`; `gene_type` to `''`. After `select_isoforms`: `canonical` (transcript row, -1 = none). |
| `G.transcripts` | `transcript_id`, `gene` (int32 gene row) | One row per transcript. After `select_isoforms`: `tss_score` (float, NaN = not scored), `tss_support` (bool). |
| `G.features` | `kind` (int8), `transcript` (int32 transcript row), `exon_number` (int32) | One row per exon / CDS / UTR record. `kind`: 0 exon, 1 CDS, 2 5'UTR, 3 3'UTR (an Ensembl-style generic `UTR` is resolved to 5' or 3' against the transcript's CDS). |

```python
genes = Genes.make("genes.gtf")
genes
# -> Genes(3 genes, 4 transcripts, 7 exons)
genes.genes.columns, genes.transcripts.columns, genes.features.columns
# -> (['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type'],
#     ['chrom', 'start', 'end', 'strand', 'transcript_id', 'gene'],
#     ['chrom', 'start', 'end', 'strand', 'kind', 'transcript', 'exon_number'])
genes.genes.to_records()                              # GTF 1001-5000 -> [1000, 5000)
# -> [('chr1', 1000, 5000, '+'), ('chr1', 10000, 11000, '-'), ('chr2', 20000, 30000, '+')]
genes.transcripts["gene"], genes.transcripts["transcript_id"]
# -> (array([0, 0, 1, 2], dtype=int32), array(['T1', 'T1b', 'T2', 'T3'], dtype=object))
genes.features["kind"][:8], genes.features["transcript"][:8], genes.features["exon_number"][:8]
# -> (array([0, 0, 0, 1, 1, 1, 2, 3], dtype=int8), array([0, 0, 0, 0, 0, 0, 0, 0], dtype=int32), array([1, 2, 3, 1, 2, 3, 1, 3], dtype=int32))
```

### Module constants and helpers

| Name | Value | One line |
|---|---|---|
| `KINDS` | `('exon', 'CDS', '5UTR', '3UTR')` | Names of the feature kinds 0–3. |
| `LABELS` | `array(['Intergenic', 'Intronic', 'Exonic', '3UTR', '5UTR', 'Promoter-TSS'])` | Region labels, indexed by what `labels()` returns. |
| `tss_base(starts, ends, strands)` | `ndarray` | 0-based TSS position per row: `start` on `+` and unstranded rows, `end - 1` where `strands == 2`. |

```python
tss_base([10, 10], [20, 20], [1, 2])
# -> array([10, 19])
```

---

## `Genes`
{: .sec-navy }

```python
Genes(genes, transcripts, features, *, promoter_r=1000, filename=None)
```

The constructor takes the three `Loci` (each must carry its link columns —
it raises otherwise, pointing at `make` / `make_ucsc` / `from_frame`).
`promoter_r` is the TSS half-window used by the annotation index,
`nearest_tss` and `select_isoforms`.

The object behaves like its genes table: `len`, `shape`, `columns`,
`head()`, `tail()`, `describe()`, `to_pandas()` and the Arrow / dataframe /
narwhals protocols read the genes; `to_pandas('transcripts')` and
`to_pandas('features')` give the others.

```python
Genes(gb.Loci(), gb.Loci(), gb.Loci())
# -> ValueError: Genes: the genes table lacks the column(s) ['gene_id', 'gene_name', 'gene_type']; build with Genes.make / make_ucsc / from_frame
genes.columns, genes.shape, len(genes)
# -> (['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type'], (3, 7), 3)
genes.describe()
# ->                                  value
#    genes                                3
#    transcripts                          4
#    exons                                7
#    CDS                                  4
#    UTR                                  2
#    chromosomes                          2
#    gene types                           2
#    canonical isoforms        not selected
#    promoter_r                        1000
#    coordinates         0-based, half-open
genes.counts()
# -> {'genes': 3, 'transcripts': 4, 'exons': 7, 'CDS': 4, 'UTR': 2}
```

### Construction

| Method | Signature | One line |
|---|---|---|
| `Genes.make` | `make(filename, *, gene_name_key="gene_name", gene_type_key="gene_type", promoter_r=1000, genome=None, chr_map=None, cre=None, bw=None, r=None, kw=None, backend=None)` | Parse a GTF or GFF3 (GENCODE / Ensembl / RefSeq, `.gz` too). |
| `Genes.make_ucsc` | `make_ucsc(filename, *, promoter_r=1000, genome=None, chr_map=None, keep_alt_contigs=False, cre=None, bw=None, r=None, kw=None)` | Parse a UCSC genePred table (refGene, ncbiRefSeq, knownGene, refFlat; with or without the leading `bin` column). |
| `Genes.from_frame` | `from_frame(df, *, one_based=None, promoter_r=1000, genome=None)` | From a GTF-like table — one row per gene / transcript / exon / CDS / UTR record (`pyranges.read_gtf`, gtfparse, a polars frame). |
| `Genes.load` | `load(path, *, genome=None)` | Read back a directory written by [`save`](#export-and-persistence). |

**`make`** detects GFF3 by extension or by sniffing the attribute column;
GFF3 links follow `ID` / `Parent` and use `gene_id` / `transcript_id` /
`gene_name` / `gene_type` when present (GENCODE), else `ID` / `Name` /
`biotype` (Ensembl, RefSeq). `gene_name_key` / `gene_type_key` name the
attributes to read. `chr_map` renames chromosomes while parsing
(`{'1': 'chr1'}`). `backend` picks the [tables backend]({{ '/api/backends/' | relative_url }})
(`'polars'` by default, else `'pandas'`; identical tables). `cre` / `bw` /
`r` / `kw` run [`select_isoforms`](#select_isoforms) right after parsing
(`kw` is a dict of its keywords).

```python
Genes.make("genes.gtf", chr_map={"chr1": "1"}).genes.chroms
# -> array(['1', '1', 'chr2'], dtype=object)
Genes.make("genes.gtf", backend="pandas").genes.equals(genes.genes, cols=True)
# -> True
g3 = Genes.make("genes.gtf", cre=[("chr1", 1900, 2100)], r=200, kw={"verbose": False})
g3.genes["canonical"]
# -> array([1, 2, 3])
```

**`make_ucsc`** keeps one gene per symbol and chromosome: copies on other
chromosomes (alt haplotypes) are skipped unless `keep_alt_contigs=True`,
which keys them `SYMBOL__chrom`. CDS and UTR features are cut from the
exons with `cdsStart` / `cdsEnd`; genePred is already 0-based.

```python
open("refFlat.txt", "w").write("GENE_A\tNM_1\tchr1\t+\t1000\t5000\t1100\t4300\t3\t1000,2000,4000,\t1200,2500,5000,\n")
u = Genes.make_ucsc("refFlat.txt")
u, u.counts(), u.genes.to_records()
# -> (Genes(1 genes, 1 transcripts, 3 exons), {'genes': 1, 'transcripts': 1, 'exons': 3, 'CDS': 3, 'UTR': 2}, [('chr1', 1000, 5000, '+')])
```

**`from_frame`** finds columns by name (`Chromosome / seqname / chrom`,
`Feature / feature`, `Start / start`, `End / end`, `Strand / strand`,
`gene_id`, `transcript_id`, `gene_name`, `gene_type` or `gene_biotype`,
`exon_number`) and needs `gene_id` and `transcript_id` to link the records.
`one_based=None` treats PyRanges objects and capitalised `Start` as 0-based
and everything else as GTF-style 1-based.

```python
import pandas as pd
df = pd.read_csv("genes.gtf", sep="\t", header=None,
                 names=["chrom", "src", "feature", "start", "end", "score", "strand", "frame", "attr"])
for k in ("gene_id", "transcript_id", "gene_name"):
    df[k] = df["attr"].str.extract(rf'{k} "([^"]+)"')
Genes.from_frame(df).genes.to_records()[0]           # 1-based starts, converted
# -> ('chr1', 1000, 5000, '+')
```

### Lookups

| Name | Returns | One line |
|---|---|---|
| `G['TP53']`, `G['ENSG...']`, `G[row]` | `GeneView` | One gene by name, id or row. |
| `name in G` | `bool` | Known gene name or id. |
| `G.find(pattern)` | `list[GeneView]` | Genes whose name contains `pattern` (case-insensitive). |
| `G.rows(names)` | `int64` array | Gene rows for names / ids (`KeyError` on an unknown one). |
| `for gene in G` | `GeneView` | Iterate the genes. |

```python
genes["G2"].tss, genes[1].gene_name, "GENE_C" in genes, genes.rows(["GENE_A", "G2"])
# -> (Locus(chrom='chr1', start=10999, end=11000, strand='-'), 'GENE_B', True, array([0, 1]))
[g.gene_name for g in genes.find("gene_")]
# -> ['GENE_A', 'GENE_B', 'GENE_C']
```

### TSS and annotation

| Method | Signature | Returns |
|---|---|---|
| `get_tss` | `get_tss(gene_type=None) -> Loci` | The 1-bp TSS of every gene (or only those of `gene_type`, one type or a list), on the gene's strand, with `gene_name`, `gene_id` and `gene` (row) columns. |
| `annot` | property → `dict` | The annotation index, built once: `body` (the genes table), `prom` (TSS ± `promoter_r`, merged), `exon`, `utr5`, `utr3` (merged). |
| `labels` | `labels(L, *, backend=None) -> int8 array` | Region label code per row of `L` (index into `LABELS`); the highest class wins: Promoter-TSS > 5UTR > 3UTR > Exonic > Intronic > Intergenic. |
| `annotations` | `annotations(L, *, backend=None) -> DataFrame` | `uid`, `annotation` per row of `L`. |
| `nearest_tss` | `nearest_tss(L, *, backend=None) -> (names, distances)` | Nearest promoter window (TSS ± `promoter_r`) per row: `('', -1)` where the chromosome has no gene. |
| `nearest_genes` | `nearest_genes(L, *, backend=None) -> DataFrame` | `uid`, `gene_name`, `distance` for the rows of `L` that have a gene on their chromosome. |

`L` is anything [`as_loci`]({{ '/api/interop/' | relative_url }}) takes;
`backend` picks the intervals engine.

```python
genes.get_tss().to_records(), genes.get_tss().columns
# -> ([('chr1', 1000, 1001, '+'), ('chr1', 10999, 11000, '-'), ('chr2', 20000, 20001, '+')],
#     ['chrom', 'start', 'end', 'strand', 'gene_name', 'gene_id', 'gene'])
genes.get_tss("protein_coding").to_records()
# -> [('chr1', 1000, 1001, '+'), ('chr2', 20000, 20001, '+')]
genes.annot.keys(), genes.annot["prom"].to_records()
# -> (dict_keys(['body', 'prom', 'exon', 'utr5', 'utr3']), [('chr1', 0, 2001, '+'), ('chr1', 9999, 12000, '-'), ('chr2', 19000, 21001, '+')])
```

```python
cre = gb.Loci.make("peaks.bed")
genes.labels(cre)
# -> array([5, 5, 3, 5, 5, 0, 0], dtype=int8)
genes.annotations(cre)
# ->                    uid    annotation
#    0     chr1:900-1100(.)  Promoter-TSS
#    1    chr1:1900-2100(.)  Promoter-TSS
#    2    chr1:4900-5100(.)          3UTR
#    3   chr1:9950-10050(.)  Promoter-TSS
#    4  chr1:10900-11100(.)  Promoter-TSS
#    5      chr2:500-600(.)    Intergenic
#    6    chr2:5000-5100(.)    Intergenic
genes.nearest_tss(cre)
# -> (array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_C', 'GENE_C'], dtype=object),
#     array([    0,     0,  2899,     0,     0, 18400, 13900]))
genes.nearest_genes(cre).tail(2)
# ->                  uid gene_name  distance
#    5    chr2:500-600(.)    GENE_C     18400
#    6  chr2:5000-5100(.)    GENE_C     13900
genes.labels(cre, backend="bioframe").tolist()
# -> [5, 5, 3, 5, 5, 0, 0]
```

The example `cre` here is the BED read without `keep`, so its uids carry
strand `.`; `annotations` labels by position only.

### `select_isoforms`

```python
G.select_isoforms(cre=None, bw=None, *, r=None, agg="max", min_signal=0.0, min_frac=0.5,
                  rank="longest", collapse=True, verbose=True) -> Genes
```

Point each gene at the isoform its cells actually use. A transcript is
*supported* when its TSS window (TSS ± `r`, default `promoter_r`) overlaps a
peak in `cre` and, with `bw`, its window score (`agg` over the window, max
over the bigWigs) is `> min_signal` and `>= min_frac` × the best candidate
score of its gene. Per gene the winner is the longest supported isoform
(`rank='longest'`) or the strongest (`'signal'`), the other breaking ties,
then the transcript id; genes with no supported isoform fall back to all of
theirs, so every gene gets a canonical transcript. `cre` is anything
`as_loci` takes (a list of peak sets is merged); `bw` one or several bigWig
paths / open handles or a `{name: path}` dict. At least one of `cre` / `bw`
is required.

Writes `transcripts['tss_score']` (NaN = not scored),
`transcripts['tss_support']` and `genes['canonical']` (transcript row, -1 =
none). With `collapse=True` the gene body — and so its TSS, `annot`,
`get_tss()` and `nearest_*` — moves onto the winner. No transcript is
dropped. Returns `self`.

```python
genes = Genes.make("genes.gtf")
genes.select_isoforms([("chr1", 1900, 2100)], r=200)         # a peak at T1b's TSS
# -> [INFO] Isoform support: 1/3 genes with an open TSS (1/4 isoforms), 2 fell back to the longest isoform.
genes.genes["canonical"], genes.transcripts["tss_support"], genes.transcripts["tss_score"]
# -> (array([1, 2, 3]), array([False,  True, False, False]), array([nan, nan, nan, nan]))
genes["GENE_A"], genes["GENE_A"].canonical, genes["GENE_A"].tss     # collapsed onto T1b
# -> (Gene(GENE_A, chr1:2,000-5,000(+), protein_coding), 'T1b', Locus(chrom='chr1', start=2000, end=2001, strand='+'))
genes.to_pandas("genes")[["gene_name", "start", "end", "canonical"]]
# ->   gene_name  start    end canonical
#    0    GENE_A   2000   5000       T1b
#    1    GENE_B  10000  11000        T2
#    2    GENE_C  20000  30000        T3
g2 = Genes.make("genes.gtf").select_isoforms([("chr1", 1900, 2100)], {"ATAC": "signal.bw"}, r=200,
                                             rank="signal", verbose=False)
g2.transcripts["tss_score"]                                   # only the candidate was scored
# -> array([       nan, 7.03561974,        nan,        nan])
```

{: .warning }
Peaks and bigWigs must use the same chromosome names as the annotation (i.e.
after `chr_map`). When nothing is supported, `verbose` prints a warning to
that effect and every gene keeps its longest isoform.

### `representative`

```python
G.representative(*, canonical=True) -> int64 array
```

One transcript row per gene: the canonical one when `select_isoforms` has
run (and `canonical=True`), else the longest. This is the row set `to_bed12`
draws.

```python
genes.representative(), genes.representative(canonical=False)
# -> (array([1, 2, 3]), array([0, 2, 3]))
```

### Export and persistence

| Method | Signature | One line |
|---|---|---|
| `to_pandas` | `to_pandas(table="genes") -> DataFrame` | One table with its links spelled out: `'genes'` (+ `canonical` transcript id when selected), `'transcripts'` (+ `gene_id`, `gene_name`), `'features'` (+ `feature` name, `transcript_id`, `gene_id`). |
| `to_polars`, `to_arrow` | `to_polars(table="genes")`, `to_arrow(table="genes")` | The same through pandas (a copy). |
| `to_gtf` | `to_gtf(path) -> path` | Write the tables back as a GTF (1-based): genes, then transcripts, then features. |
| `to_bed12` | `to_bed12(path=None, *, canonical=True)` | One BED12 line per gene — its representative isoform — with exons as blocks and the CDS as thick region (`thickStart == thickEnd` when non-coding); text when `path` is `None`. |
| `save` | `save(path)` | `path/` gets `genes.parquet`, `transcripts.parquet`, `features.parquet` and `meta.json` (`promoter_r`, `filename`, coordinate convention). |
| `Genes.load` | `load(path, *, genome=None)` | Read it back. |
| `head`, `tail` | `head(n=5)`, `tail(n=5) -> DataFrame` | The genes table as pandas. |

```python
genes = Genes.make("genes.gtf")
genes.head(2)
# ->   chrom  start    end strand gene_id gene_name       gene_type
#    0  chr1   1000   5000      +      G1    GENE_A  protein_coding
#    1  chr1  10000  11000      -      G2    GENE_B          lncRNA
genes.to_pandas("transcripts")
# ->   chrom  start    end strand transcript_id  gene gene_id gene_name
#    0  chr1   1000   5000      +            T1     0      G1    GENE_A
#    1  chr1   2000   5000      +           T1b     0      G1    GENE_A
#    2  chr1  10000  11000      -            T2     1      G2    GENE_B
#    3  chr2  20000  30000      +            T3     2      G3    GENE_C
genes.to_pandas("features").head(4)
# ->   chrom  start   end strand  kind  transcript  exon_number feature transcript_id gene_id
#    0  chr1   1000  1200      +     0           0            1    exon            T1      G1
#    1  chr1   2000  2500      +     0           0            2    exon            T1      G1
#    2  chr1   4000  5000      +     0           0            3    exon            T1      G1
#    3  chr1   1100  1200      +     1           0            1     CDS            T1      G1
genes.to_polars("genes").shape, genes.to_arrow("features").num_rows
# -> ((3, 7), 13)
import polars as pl
pl.DataFrame(genes).shape                       # the protocols read the genes table
# -> (3, 7)
```

```python
genes.to_bed12().splitlines()
# -> ['chr1\t1000\t5000\tGENE_A\t0\t+\t1100\t4300\t0\t3\t200,500,1000\t0,1000,3000',
#     'chr1\t10000\t11000\tGENE_B\t0\t-\t11000\t11000\t0\t2\t400,400\t0,600',
#     'chr2\t20000\t30000\tGENE_C\t0\t+\t20100\t29900\t0\t1\t10000\t0']
genes.to_gtf("out.gtf")
open("out.gtf").readline()
# -> 'chr1\tgenomeblocks\tgene\t1001\t5000\t.\t+\t.\tgene_id "G1"; gene_name "GENE_A"; gene_type "protein_coding";\n'
import os
genes.save("genes_dir")
sorted(os.listdir("genes_dir"))
# -> ['features.parquet', 'genes.parquet', 'meta.json', 'transcripts.parquet']
Genes.load("genes_dir")
# -> Genes(3 genes, 4 transcripts, 7 exons)
```

---

## `GeneView(Locus)`
{: .sec-navy }

```python
GeneView(genes: Genes, i: int)          # made by Genes[key] and iteration; not built by hand
```

One gene as a [`Locus`]({{ '/api/locus/' | relative_url }}) reading the
genes table at row `i`, with its transcripts and exons one attribute away.
Nothing is copied.

| Name | Type | One line |
|---|---|---|
| `row` | `int` | Row in the genes table. |
| `chrom`, `start`, `end`, `strand` | | Read from the genes table (read-only). |
| `gene_id`, `gene_name`, `gene_type` | `str` | The gene columns. |
| `tss` | `Locus` | The 1-bp TSS on the gene's strand. |
| `transcripts` | `Loci` | The rows of `G.transcripts` whose `gene` is this row. |
| `exons` | `Loci` | The exon rows (`kind == 0`) of those transcripts. |
| `canonical` | `str` or `None` | `transcript_id` chosen by `select_isoforms` (`None` if not run). |
| `uid`, `length`, `center`, `overlaps`, `distance_to`, `sequence` | | Inherited from `Locus`. |

```python
g = genes["GENE_A"]
g, type(g).__name__, g.row
# -> (Gene(GENE_A, chr1:1,000-5,000(+), protein_coding), 'GeneView', 0)
g.gene_id, g.gene_name, g.gene_type, g.tss, g.canonical
# -> ('G1', 'GENE_A', 'protein_coding', Locus(chrom='chr1', start=1000, end=1001, strand='+'), None)
g.transcripts, g.transcripts["transcript_id"]
# -> (Loci(n=2, chroms=1, cols=[transcript_id, gene]), array(['T1', 'T1b'], dtype=object))
g.exons.to_records()
# -> [('chr1', 1000, 1200, '+'), ('chr1', 2000, 2500, '+'), ('chr1', 4000, 5000, '+'), ('chr1', 2000, 5000, '+')]
genes["GENE_B"].exons.to_records()
# -> [('chr1', 10000, 10400, '-'), ('chr1', 10600, 11000, '-')]
for gv in genes:
    print(gv)
# -> Gene(GENE_A, chr1:1,000-5,000(+), protein_coding)
#    Gene(GENE_B, chr1:10,000-11,000(-), lncRNA)
#    Gene(GENE_C, chr2:20,000-30,000(+), protein_coding)
```

{: .tip }
Whole-annotation work belongs on the tables, not on views: `genes.get_tss()`
gives every TSS as one `Loci`, `genes.transcripts.take(genes.representative())`
the representative isoforms, and `genes.features.take(genes.features["kind"] == 1)`
every CDS.
