---
title: Genes
parent: User Guide
layout: default
nav_order: 3
---

# Genes
{: .no_toc }

Gene annotations as three linked tables — genes, transcripts, features — parsed
from GTF / GFF3 or UCSC genePred, in the same 0-based coordinates as every
other table, plus the toolkit for tying CREs to genes: region classes, nearest
TSS, and isoform selection from your own ATAC data.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## What it is
{: .sec-navy }

A `Genes` object holds three [`Loci`]({{ '/guide/loci/' | relative_url }})
tables on one `Genome`. The links between them are row numbers:

| Table | One row per | Columns | Link |
|---|---|---|---|
| `genes.genes` | gene | `gene_id`, `gene_name`, `gene_type` (+ `canonical` after isoform selection) | |
| `genes.transcripts` | transcript | `transcript_id`, `gene` | `gene` = row in `genes.genes` |
| `genes.features` | exon / CDS / UTR | `kind` (0 exon, 1 CDS, 2 5'UTR, 3 3'UTR), `transcript`, `exon_number` | `transcript` = row in `genes.transcripts` |

```python
import genomeblocks as gb

genes = gb.Genes.make("gencode.v44.annotation.gtf")
genes
# -> Genes(3 genes, 4 transcripts, 7 exons)
genes.counts()
# -> {'genes': 3, 'transcripts': 4, 'exons': 7, 'CDS': 4, 'UTR': 2}
genes.transcripts.cols
# -> {'transcript_id': array(['T1', 'T1b', 'T2', 'T3'], dtype=object), 'gene': array([0, 0, 1, 2], dtype=int32)}
```

**Coordinates.** Every table is 0-based, half-open. A GTF / GFF3 record
`start end` (1-based, closed) is stored as `[start - 1, end)`; UCSC genePred
tables are already 0-based and are stored as they are. So gene bodies, exons
and promoters compare directly with peaks, with no off-by-one.

**The TSS rule.** A gene's TSS is the 1-bp interval `[t, t + 1)` on its strand,
with `t = start` on `+` (and unstranded) and `t = end - 1` on `-`:

```python
genes["GENE_A"].tss                 # GTF: chr1 1001-5000 (+)  ->  [1000, 5000)
# -> Locus(chrom='chr1', start=1000, end=1001, strand='+')
genes["GENE_B"].tss                 # GTF: chr1 10001-11000 (-) -> [10000, 11000)
# -> Locus(chrom='chr1', start=10999, end=11000, strand='-')
```

As a table, the object is its genes table: `len`, `shape`, `columns`,
`head()`, `tail()`, `describe()`, iteration, and the Arrow / dataframe
protocols all read `genes.genes`. `head()` / `tail()` return pandas frames;
`to_pandas('transcripts')` and friends give the other tables.

```python
genes.shape, genes.columns
# -> ((3, 7), ['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type'])
genes.head()
# ->   chrom  start    end strand gene_id gene_name       gene_type
# -> 0  chr1   1000   5000      +      G1    GENE_A  protein_coding
# -> 1  chr1  10000  11000      -      G2    GENE_B          lncRNA
# -> 2  chr2  20000  30000      +      G3    GENE_C  protein_coding
genes.describe()
# ->                                  value
# -> genes                                3
# -> transcripts                          4
# -> exons                                7
# -> CDS                                  4
# -> UTR                                  2
# -> chromosomes                          2
# -> gene types                           2
# -> canonical isoforms        not selected
# -> promoter_r                        1000
# -> coordinates         0-based, half-open
```

---

## Parsing GTF / GFF3: `Genes.make`
{: .sec-navy }

```python
genes = gb.Genes.make("gencode.v44.annotation.gtf", promoter_r=1000)
```

| arg | default | meaning |
|---|---|---|
| `promoter_r` | 1000 | half-width of the promoter window (TSS ± `promoter_r`) used by `annot`, `labels`, `nearest_tss`, and the default `r` of `select_isoforms` |
| `gene_name_key`, `gene_type_key` | `'gene_name'`, `'gene_type'` | the attributes that hold the symbol and the biotype in a GTF (Ensembl GTFs use `gene_biotype`) |
| `chr_map` | `None` | rename chromosomes while parsing: `{'1': 'chr1', 'MT': 'chrM'}` |
| `genome` | `None` | a `Genome` to share with other tables |
| `cre`, `bw`, `r`, `kw` | `None` | run `select_isoforms` right after parsing (below) |
| `backend` | `None` | the table parser: `'polars'` (default when installed) or `'pandas'`; both give identical tables |

The format is decided by extension (`.gtf`, `.gff`, `.gff3`, `.gz` of each)
or, failing that, by sniffing the attribute column. Feature rows of type
`gene`, `transcript` / `mRNA` (any non-feature child of a gene in GFF3), `exon`,
`CDS`, `five_prime_UTR`, `three_prime_UTR` and the generic `UTR` are read;
a generic `UTR` is classed 5' or 3' from its position relative to the
transcript's first CDS base. Other feature types (`start_codon`,
`Selenocysteine`, ...) are ignored. Duplicate `gene_id`s keep their first
record. Transcripts whose gene has no `gene` record get a gene made from them.

The spellings of the three sources all parse without options except
`chr_map`:

| Source | File | Gene id / symbol / type | Links |
|---|---|---|---|
| GENCODE GTF | `gencode.v44.annotation.gtf` | `gene_id`, `gene_name`, `gene_type` | `transcript_id` attribute |
| GENCODE GFF3 | `gencode.v44.annotation.gff3` | `gene_id`, `gene_name`, `gene_type` attributes (the `ID` is the fallback) | `ID` / `Parent`; `exon_number` |
| Ensembl GTF | `Homo_sapiens.GRCh38.110.gtf` | `gene_id`, `gene_name`, `gene_biotype` → pass `gene_type_key="gene_biotype"` | `transcript_id`; chromosomes `1`, `2`, `X` → `chr_map` |
| Ensembl GFF3 | `Homo_sapiens.GRCh38.110.gff3` | `ID=gene:ENSG...`, `Name`, `biotype` (the `gene:` / `transcript:` prefixes are stripped) | `ID` / `Parent`; `rank` as the exon number |
| RefSeq GFF3 | `GCF_000001405.40_GRCh38.p14_genomic.gff` | `ID=gene-TP53`, `Name`, `gene_biotype`; transcripts take `transcript_id` (`NM_...`) | `ID` / `Parent`; `NC_000001.11` names → `chr_map` |

```python
ens = gb.Genes.make("Homo_sapiens.GRCh38.110.gtf", chr_map={"1": "chr1", "2": "chr2"},
                    gene_type_key="gene_biotype")
ens.genes.to_pandas()
# ->   chrom  start    end strand gene_id gene_name       gene_type
# -> 0  chr1   1000   5000      +      G1    GENE_A  protein_coding
# -> 1  chr1  10000  11000      -      G2    GENE_B          lncRNA
# -> 2  chr2  20000  30000      +      G3    GENE_C  protein_coding

refseq = gb.Genes.make("GCF_000001405.40_GRCh38.p14_genomic.gff", chr_map={"NC_000001.11": "chr1"})
refseq.to_pandas("transcripts")
# ->   chrom  start   end strand transcript_id  gene      gene_id gene_name
# -> 0  chr1   1000  5000      +   NM_000001.1     0  gene-GENE_A    GENE_A
```

{: .warning }
> Peaks and bigWigs you compare with the annotation must use the **same
> chromosome names** as the parsed tables — that is, the names after
> `chr_map`. An Ensembl annotation read without `chr_map` has chromosomes
> `1`, `2`, ...; a `chr1` peak set then overlaps nothing.

---

## Parsing UCSC tables: `Genes.make_ucsc`
{: .sec-navy }

UCSC genePred dumps — `refGene.txt`, `ncbiRefSeq.txt`, `knownGene.txt`,
`refFlat.txt`, plain genePred and genePredExt — are one transcript per line
with the exon starts and ends as comma lists. The layout is found from the
position of the strand column, with or without the leading `bin` column.

```python
ref = gb.Genes.make_ucsc("ncbiRefSeq.txt", promoter_r=1000)
ref
# -> Genes(2 genes, 2 transcripts, 5 exons)
ref.to_pandas("transcripts")
# ->   chrom  start    end strand transcript_id  gene gene_id gene_name
# -> 0  chr1   1000   5000      +          NM_1     0  GENE_A    GENE_A
# -> 1  chr1  10000  11000      -          NM_2     1  GENE_B    GENE_B
```

- Genes are built from the transcripts: one gene per symbol (`name2` in
  refGene / ncbiRefSeq / genePredExt, `geneName` in refFlat; knownGene and
  plain genePred have no symbol, so the transcript name is used), spanning
  the union of its transcripts, keyed by the symbol. `gene_type` is empty: the
  tables carry none.
- CDS and UTR features are cut from the exons using `cdsStart` / `cdsEnd`
  (a non-coding transcript has `cdsStart == cdsEnd` and gets exons only).
  `exon_number` counts from the 5' end, so it runs backwards on `-`.
- Copies of a symbol on another chromosome (alt haplotypes such as
  `chr6_GL000251v2_alt`) are **skipped**; `keep_alt_contigs=True` keeps them
  as separate genes keyed `SYMBOL__chrom`.
- A transcript id repeated within one gene (MHC paralogs) gets a `__2`,
  `__3` suffix — nothing is dropped silently.
- `chr_map`, `genome`, `cre` / `bw` / `r` / `kw` work as in `make`.

A GTF handed to `make_ucsc` (or a genePred to `make`) raises an error that
names the right constructor.

---

## From a frame: `Genes.from_frame`
{: .sec-navy }

A GTF-like frame — one row per gene / transcript / exon / CDS / UTR record —
from `pyranges.read_gtf`, gtfparse, a polars frame or a DataFrame you parsed
yourself. Columns are found by name: `Chromosome` / `seqname` / `chrom`,
`Feature` / `feature`, `Start`, `End`, `Strand`, `gene_id`, `transcript_id`,
`gene_name`, `gene_type` / `gene_biotype`, `exon_number`. `gene_id` and
`transcript_id` are required (they link the records).

```python
import pyranges as pr
gr = pr.read_gtf("gencode.v44.annotation.gtf")
gb.Genes.from_frame(gr).genes.equals(genes.genes)
# -> True
```

`one_based` says whether starts are GTF-style. By default a PyRanges object
(or a frame with a capitalised `Start`) is taken as 0-based and anything else
as 1-based; pass `one_based=True` / `False` to be explicit.

---

## Looking up genes
{: .sec-navy }

`genes[key]` with a gene name, a gene id or a row number returns a `GeneView`:
a `Locus` reading the genes table, with the gene's attributes and its
transcripts and exons one attribute away.

```python
g = genes["GENE_A"]                 # or genes["G1"], genes[0]
g
# -> Gene(GENE_A, chr1:1,000-5,000(+), protein_coding)
g.row, g.gene_id, g.gene_type, g.uid, g.length
# -> (0, 'G1', 'protein_coding', 'chr1:1000-5000(+)', 4000)
g.tss
# -> Locus(chrom='chr1', start=1000, end=1001, strand='+')
g.transcripts                       # a Loci: the rows of genes.transcripts for this gene
# -> Loci(n=2, chroms=1, cols=[transcript_id, gene])
g.transcripts["transcript_id"]
# -> array(['T1', 'T1b'], dtype=object)
g.exons.to_records()                # every exon of every isoform
# -> [('chr1', 1000, 1200, '+'), ('chr1', 2000, 2500, '+'), ('chr1', 4000, 5000, '+'), ('chr1', 2000, 5000, '+')]
g.canonical                         # the transcript_id chosen by select_isoforms, else None
# -> None
```

```python
"GENE_C" in genes, "nope" in genes
# -> (True, False)
genes.rows(["GENE_A", "G2"])        # gene rows for names / ids (KeyError on an unknown one)
# -> array([0, 1])
genes.find("gene_")                 # substring match on the name, case-insensitive
# -> [Gene(GENE_A, chr1:1,000-5,000(+), protein_coding), Gene(GENE_B, chr1:10,000-11,000(-), lncRNA), Gene(GENE_C, chr2:20,000-30,000(+), protein_coding)]
[x.gene_name for x in genes]
# -> ['GENE_A', 'GENE_B', 'GENE_C']
```

For whole-table work use the tables directly: `genes.genes["gene_type"]` is a
numpy array, `genes.genes[mask]` a `Loci`, and
`genes.transcripts["gene"]` tells you which gene row each transcript belongs
to.

---

## TSS
{: .sec-navy }

`get_tss()` returns the 1-bp TSS of every gene as a `Loci` on the gene's
strand, with `gene_name`, `gene_id` and `gene` (the gene row) columns. A
`gene_type` (or a list of them) filters:

```python
tss = genes.get_tss()
tss.to_records()
# -> [('chr1', 1000, 1001, '+'), ('chr1', 10999, 11000, '-'), ('chr2', 20000, 20001, '+')]
genes.get_tss("protein_coding").to_records()
# -> [('chr1', 1000, 1001, '+'), ('chr2', 20000, 20001, '+')]
promoters = tss.slop(1000)          # TSS ± 1 kb windows, same rows as tss
```

Because the TSS table has one row per gene in gene order, `tss.slop(r)` is the
promoter table, and its `gene` column maps any result back to `genes.genes`.

---

## Annotating loci
{: .sec-navy }

### The annotation windows

`genes.annot` is a dict of `Loci` built on first use and cached: the raw
material of the region classes. They are ordinary tables you can intersect
with, use as promoter sources for an
[Architecture]({{ '/guide/architecture/' | relative_url }}), or pass to
`Pairs.overlapping`.

| Key | Content |
|---|---|
| `annot['body']` | the genes table itself |
| `annot['prom']` | TSS ± `promoter_r`, merged |
| `annot['exon']` | all exons of all isoforms, merged |
| `annot['utr5']`, `annot['utr3']` | all 5' / 3' UTRs, merged |

```python
genes.annot["prom"].to_records()
# -> [('chr1', 0, 2001, '+'), ('chr1', 9999, 12000, '-'), ('chr2', 19000, 21001, '+')]
enhancers = gb.Loci.make("peaks.bed") - genes.annot["prom"]
```

### Region class per locus: `labels` / `annotations`

```python
cre = gb.Loci.make("peaks.bed")
genes.annotations(cre)
# ->                    uid    annotation
# -> 0     chr1:900-1100(+)  Promoter-TSS
# -> 1    chr1:1900-2100(-)  Promoter-TSS
# -> 2    chr1:4900-5100(+)          3UTR
# -> 3   chr1:9950-10050(.)  Promoter-TSS
# -> 4  chr1:10900-11100(-)  Promoter-TSS
# -> 5      chr2:500-600(+)    Intergenic
# -> 6    chr2:5000-5100(+)    Intergenic
```

`annotations(L)` returns a DataFrame (`uid`, `annotation`); `labels(L)` the
same thing as an int8 code per row, indexing `genomeblocks.genes.LABELS`
(`Intergenic`, `Intronic`, `Exonic`, `3UTR`, `5UTR`, `Promoter-TSS`), which is
what you want to store as a column. The highest-priority class a row overlaps
wins: **Promoter-TSS > 5UTR > 3UTR > Exonic > Intronic > Intergenic**. A row
inside a gene body that touches no exon is Intronic.

```python
from genomeblocks.genes import LABELS
cre["annot"] = LABELS[genes.labels(cre)]
narrow = gb.Genes.make("gencode.v44.annotation.gtf", promoter_r=100)
LABELS[narrow.labels(cre)]
# -> array(['Promoter-TSS', 'Exonic', '3UTR', 'Exonic', 'Promoter-TSS', 'Intergenic', 'Intergenic'], dtype=object)
```

`promoter_r` is the one knob: with ±1 kb promoters two of the rows above are
promoter; with ±100 bp they fall into exons.

### Nearest gene: `nearest_tss` / `nearest_genes`

Both look for the nearest **promoter window** (TSS ± `promoter_r`), so a row
anywhere inside a promoter is 0 bp from that gene; otherwise the distance is
the gap in bases to the window's edge (a row 1 kb past the window is 1000, not
2000).

```python
names, dist = genes.nearest_tss(cre)        # two arrays aligned to the rows of cre
names, dist
# -> (array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_C', 'GENE_C'], dtype=object), array([    0,     0,  2899,     0,     0, 18400, 13900]))
genes.nearest_genes(cre)                    # a DataFrame; rows with no gene on their chromosome are left out
# ->                    uid gene_name  distance
# -> 0     chr1:900-1100(+)    GENE_A         0
# -> 1    chr1:1900-2100(-)    GENE_A         0
# -> 2    chr1:4900-5100(+)    GENE_A      2899
# -> 3   chr1:9950-10050(.)    GENE_B         0
# -> 4  chr1:10900-11100(-)    GENE_B         0
# -> 5      chr2:500-600(+)    GENE_C     18400
# -> 6    chr2:5000-5100(+)    GENE_C     13900
```

A row on a chromosome with no gene gets `('', -1)` from `nearest_tss`. Ties
between two equally near windows go to the one with the lower start.

Every annotation call accepts anything `as_loci` takes (a frame, a path, a
list of regions) and takes `backend=` to pick the interval engine; the labels
are identical on every engine, and `nearest_*` on the engines that implement
nearest (`genomeblocks`, `bioframe`, `pyranges`, `bedtools`).

{: .note }
> `annotations()` and `nearest_tss()` are exactly what
> [`Architecture.annotate`]({{ '/guide/architecture/' | relative_url }}) writes
> into its `vp.annot` and `vp.gene` vertex columns.

---

## Picking the isoform your cells use: `select_isoforms`
{: .sec-navy }

A gene in a GTF spans the **union** of its isoforms, so `gene.start`,
`gene.end` and the TSS follow the longest *annotated* transcript. For a gene
with a long, rarely used isoform that TSS can sit tens of kilobases from the
promoter that is open in your cell type — which skews promoter labels,
nearest-gene calls and the browser view.

`select_isoforms` keeps the isoforms whose TSS window (TSS ± `r`) overlaps a
peak in `cre` and, with `bw`, carries signal; it then points each gene at one
of them.

```python
peaks = gb.as_loci([("chr1", 1900, 2100)])          # open over T1b's TSS, not T1's
genes = gb.Genes.make("gencode.v44.annotation.gtf").select_isoforms(peaks, r=200)
# -> [INFO] Isoform support: 1/3 genes with an open TSS (1/4 isoforms), 2 fell back to the longest isoform.
genes["GENE_A"], genes["GENE_A"].canonical, genes["GENE_A"].tss
# -> (Gene(GENE_A, chr1:2,000-5,000(+), protein_coding), 'T1b', Locus(chrom='chr1', start=2000, end=2001, strand='+'))
genes.to_pandas("genes")
# ->   chrom  start    end strand gene_id gene_name       gene_type canonical
# -> 0  chr1   2000   5000      +      G1    GENE_A  protein_coding       T1b
# -> 1  chr1  10000  11000      -      G2    GENE_B          lncRNA        T2
# -> 2  chr2  20000  30000      +      G3    GENE_C  protein_coding        T3
```

| arg | default | meaning |
|---|---|---|
| `cre` | `None` | peaks: anything `as_loci` takes; a list of several sets is merged |
| `bw` | `None` | a bigWig path or open handle, a list of them, or a `{name: path}` dict; a TSS keeps its **highest** score across them |
| `r` | `promoter_r` | TSS half-window |
| `agg` | `'max'` | the bigWig statistic over the window |
| `min_signal` | 0.0 | a supported TSS must score above this |
| `min_frac` | 0.5 | ... and at least this fraction of its gene's best TSS score |
| `rank` | `'longest'` | the winner among supported isoforms: the longest, or `'signal'` for the strongest (the other breaks ties, then the transcript id) |
| `collapse` | `True` | move the gene body (and so its TSS) onto the winner |
| `verbose` | `True` | the `[INFO]` line; a `[WARN]` when nothing at all is supported |

Either evidence argument works alone: peaks gate first, then the signal cut.
It runs in place and returns the same object, so it chains, and the
constructors run it for you:

```python
genes = gb.Genes.make("gencode.v44.annotation.gtf", cre="atac.narrowPeak", bw="atac.bw", r=200)
# -> [INFO] Isoform support: 2/3 genes with an open TSS (3/4 isoforms), 1 fell back to the longest isoform.
```

What it writes:

- `genes.transcripts['tss_support']` — bool per isoform; `['tss_score']` — the
  window score (NaN when no bigWig was given or the isoform was not a candidate).
- `genes.genes['canonical']` — the chosen transcript **row** per gene (`-1`
  for none); `genes['TP53'].canonical` gives its id.
- With `collapse=True`, the gene's `start` / `end` move to the winner, so
  `annot`, `get_tss()`, `labels()` and `nearest_tss()` all follow the
  supported isoform. `annot['exon' / 'utr5' / 'utr3']` still pool every isoform.
- **Nothing is dropped.** Every isoform stays in `transcripts`; genes with no
  supported isoform fall back to all of theirs, so every gene keeps a canonical
  transcript.

`representative()` returns one transcript row per gene — the canonical one
when selected, else the longest — and is what `to_bed12` and the browser draw:

```python
genes.representative(), genes.representative(canonical=False)
# -> (array([1, 2, 3]), array([0, 2, 3]))
```

---

## Export
{: .sec-navy }

`to_pandas(table)`, `to_polars(table)` and `to_arrow(table)` give one table
with its links spelled out: `'genes'` (+ `canonical` id when selected),
`'transcripts'` (+ `gene_id`, `gene_name`), `'features'` (+ `feature` name,
`transcript_id`, `gene_id`).

```python
genes.to_pandas("features").head(4)
# ->   chrom  start   end strand  kind  transcript  exon_number feature transcript_id gene_id
# -> 0  chr1   1000  1200      +     0           0            1    exon            T1      G1
# -> 1  chr1   2000  2500      +     0           0            2    exon            T1      G1
# -> 2  chr1   4000  5000      +     0           0            3    exon            T1      G1
# -> 3  chr1   1100  1200      +     1           0            1     CDS            T1      G1
```

The genes table also flows through the protocols: `pl.DataFrame(genes)`,
`duckdb.sql("select * from genes")`.

```python
genes.to_gtf("out.gtf")              # 1-based GTF: genes, then transcripts, then features
print(genes.to_bed12().splitlines()[0])    # one BED12 line per gene: its representative isoform, exons as blocks, CDS as thick
# -> chr1	2000	5000	GENE_A	0	+	5000	5000	0	1	3000	0
genes.save("genes_v44")              # a directory: genes / transcripts / features .parquet + meta.json
gb.Genes.load("genes_v44")
# -> Genes(3 genes, 4 transcripts, 7 exons)
```

`to_bed12(path=None, canonical=True)` writes to `path` when given, else
returns the text; `canonical=False` always draws the longest isoform. Saved
tables load onto any `Genome` (`load(path, genome=g)`) with `promoter_r` and
the source filename restored.

---

## End-to-end: annotate CREs and write a table
{: .sec-navy }

```python
import genomeblocks as gb
from genomeblocks.genes import LABELS

cre   = gb.Loci.make("atac.narrowPeak")
genes = gb.Genes.make("gencode.v44.annotation.gtf", promoter_r=1000,
                      cre=cre, bw="atac.bw")           # isoforms with an open TSS

cre["annot"] = LABELS[genes.labels(cre)]
cre["gene"], cre["distance"] = genes.nearest_tss(cre)

cre.to_pandas(uid=True).to_csv("cre_annotated.csv", index=False)
cre.save("cre_annotated.parquet")                     # columns included
```

Three columns, one per row of the catalogue, ready to group a signal cube or
a motif matrix by region class or by gene.

---

## See also
{: .sec-navy }

- [Genes API]({{ '/api/genes/' | relative_url }}) — every method and argument.
- [Design: Genes]({{ '/design/genes/' | relative_url }}) — the three-table model and the parsers.
- [Loci]({{ '/guide/loci/' | relative_url }}) — the tables underneath.
- [Architecture]({{ '/guide/architecture/' | relative_url }}) — `annotate` and hub genes.
- [Browser]({{ '/guide/browser/' | relative_url }}) — drawing genes and transcripts.
