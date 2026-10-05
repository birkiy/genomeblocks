---
title: Genes
parent: Design
layout: default
nav_order: 2
---

# Genes
{: .no_toc }

`Genes` parses a GTF, GFF3 or UCSC genePred table into three linked `Loci`
tables (genes, transcripts, features) in the same 0-based half-open
coordinates as every other table, then builds a cached index of promoter,
UTR, exon and gene-body intervals that labels any set of regions.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## From GTF lines to three tables
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/genes-model.svg %}
</div><figcaption>
<strong>Three Loci, linked by row number.</strong> <code>genes</code> has one row per gene (<code>gene_id</code>, <code>gene_name</code>, <code>gene_type</code>), <code>transcripts</code> one per transcript with <code>gene</code> = the row of its gene, and <code>features</code> one per exon / CDS / UTR record with <code>kind</code>, <code>transcript</code> = the row of its transcript, and <code>exon_number</code>. The GTF's 1-based starts become <code>start − 1</code> once, while parsing, so genes and CREs compare directly.
</figcaption></figure>

| Table | Columns beyond the coordinates | Row count |
|---|---|---|
| `genes` | `gene_id`, `gene_name`, `gene_type` (+ `canonical` after isoform selection) | one per `gene_id`, first record wins |
| `transcripts` | `transcript_id`, `gene` (int32 row into `genes`) (+ `tss_score`, `tss_support`) | one per transcript record |
| `features` | `kind` (int8: 0 exon, 1 CDS, 2 5′UTR, 3 3′UTR), `transcript` (int32 row), `exon_number` | one per exon / CDS / UTR record with a known transcript |

```python
import genomeblocks as gb

genes = gb.Genes.make("genes.gtf")
genes
# -> Genes(3 genes, 4 transcripts, 7 exons)
G, T, F = genes.genes, genes.transcripts, genes.features
G.to_records()                                  # GTF 1001-5000 became [1000, 5000)
# -> [('chr1', 1000, 5000, '+'), ('chr1', 10000, 11000, '-'), ('chr2', 20000, 30000, '+')]
T.cols["transcript_id"], T.cols["gene"]
# -> (array(['T1', 'T1b', 'T2', 'T3'], dtype=object), array([0, 0, 1, 2], dtype=int32))
F.cols["kind"], F.cols["transcript"]
# -> (array([0, 0, 0, 1, 1, 1, 2, 3, 0, 0, 0, 0, 1], dtype=int8), array([0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 3, 3], dtype=int32))
G.genome is T.genome is F.genome
# -> True
```

**The parser is the tables backend.** With polars (the default when
installed) the file is read once into a frame of `chrom, feature, start, end,
strand` plus one column per attribute extracted by regex, and the three tables
come out of polars joins on `gene_id` / `transcript_id`; only numbers cross
into numpy. Without polars, pandas reads the same columns and the joins are
`pandas.Index.get_indexer` calls. Both produce identical tables
(`Genes.make(path, backend="pandas")` is the parity check).

- **GFF3** is told apart by extension or by sniffing the attribute column
  (`key=value` vs `key "value"`); links follow `ID` / `Parent`, and GENCODE's
  `gene_id` / `transcript_id` attributes are used when present, else `ID` /
  `Name` / `biotype` (Ensembl, RefSeq). Ensembl's `rank=` is read as the exon
  number, and a generic `UTR` feature is resolved to 5′ or 3′ by its position
  relative to the transcript's first CDS base.
- **UCSC genePred** (`Genes.make_ucsc`: refGene, ncbiRefSeq, knownGene,
  refFlat, with or without the leading `bin` column) is already 0-based; the
  layout is found from the strand column, exons come from the
  `exonStarts` / `exonEnds` lists, and CDS / UTR pieces are cut from them.
- **A GTF-like frame** (`pyranges.read_gtf`, gtfparse, a polars frame) goes
  in through `Genes.from_frame`, which decides `one_based` from the input
  (PyRanges and a capitalised `Start` are 0-based) unless told.
- `chr_map` renames chromosomes while parsing; `gene_name_key` /
  `gene_type_key` pick the attribute names (GENCODE defaults); transcripts
  whose gene has no `gene` record get a gene row of their own.

**TSS.** A gene's or transcript's TSS is the 1-bp interval `[t, t + 1)` with
`t = start` on `+` (and unstranded) and `t = end − 1` on `−`
(`genes.tss_base`). `get_tss()` returns every TSS as a `Loci` with
`gene_name`, `gene_id` and `gene` (row) columns, on the gene's strand.

```python
genes["GENE_A"].tss, genes["GENE_B"].tss
# -> (Locus(chrom='chr1', start=1000, end=1001, strand='+'), Locus(chrom='chr1', start=10999, end=11000, strand='-'))
genes.get_tss().to_records()
# -> [('chr1', 1000, 1001, '+'), ('chr1', 10999, 11000, '-'), ('chr2', 20000, 20001, '+')]
```

**One gene as an object.** `genes["TP53"]`, `genes["ENSG..."]` or
`genes[row]` return a `GeneView`: a `Locus` reading the genes table, with
`transcripts`, `exons`, `tss` and `canonical` one attribute away. Nothing is
copied until you ask for a sub-table.

```python
g = genes["GENE_A"]
g, g.row, len(g.transcripts), len(g.exons)
# -> (Gene(GENE_A, chr1:1,000-5,000(+), protein_coding), 0, 2, 4)
```

## Labelling regions: the annotation index
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/genes-annot.svg %}
</div><figcaption>
<strong>The highest-priority class wins.</strong> <code>genes.annot</code> holds five merged <code>Loci</code>: promoters (every TSS ± <code>promoter_r</code>), 5′ UTRs, 3′ UTRs, exons and gene bodies. <code>labels(cres)</code> tests every CRE against each set with one <code>overlap_any</code> and labels it with the first set it touches in that order; a CRE that touches none is <em>Intergenic</em>.
</figcaption></figure>

The index is a dict of five `Loci` built on first access and cached on the
object (`select_isoforms` resets it): `body` is the genes table itself,
`prom` is `get_tss().slop(promoter_r).merge()`, and `exon`, `utr5`, `utr3`
are merged slices of the features table by `kind`.

```python
cre = gb.Loci.make("peaks.bed")
{k: len(v) for k, v in genes.annot.items()}
# -> {'body': 3, 'prom': 3, 'exon': 5, 'utr5': 1, 'utr3': 1}
genes.labels(cre)                               # int8 codes into genes.LABELS
# -> array([5, 5, 3, 5, 5, 0, 0], dtype=int8)
genes.annotations(cre.to_pandas()).annotation.tolist()
# -> ['Promoter-TSS', 'Promoter-TSS', '3UTR', 'Promoter-TSS', 'Promoter-TSS', 'Intergenic', 'Intergenic']
```

`labels` writes the codes in increasing priority (body → exon → utr3 → utr5
→ prom), so the last write is the highest class: five vectorised
`overlap_any` calls against merged sets, whatever the number of CREs.
`nearest_tss(cres)` is a separate path: one `nearest` against every promoter
window (TSS ± `promoter_r`), returning the gene name and the gap in bases
(0 inside the window, -1 when the chromosome has no gene); `nearest_genes`
is the same as a DataFrame of the rows that have one.

```python
genes.nearest_tss(cre)
# -> (array(['GENE_A', 'GENE_A', 'GENE_A', 'GENE_B', 'GENE_B', 'GENE_C', 'GENE_C'], dtype=object),
#     array([    0,     0,  2899,     0,     0, 18400, 13900]))
```

Both go through the `intervals` backend: `labels(cre, backend="pyranges")`
gives the same codes, and `nearest_tss` on an engine without `nearest`
raises `NotImplementedError` naming the ones that have it.

## Picking the isoform your cells use
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/genes-isoforms.svg %}
</div><figcaption>
<strong>Open chromatin decides the TSS.</strong> An isoform is <em>supported</em> when its TSS window (± <code>r</code>) overlaps a peak. When bigWigs are also given, its TSS signal must exceed <code>min_signal</code> and reach <code>min_frac</code> of the gene's best TSS. The longest supported isoform becomes canonical (<code>rank="signal"</code> picks the strongest TSS instead). With <code>collapse=True</code> the gene's span and TSS move to it, which every downstream label and nearest-gene call then reads.
</figcaption></figure>

`select_isoforms(cre, bw, r=, ...)` is one pass over the transcripts table:

1. the TSS windows of all transcripts are one `Loci`, and `overlap_any`
   against the peaks (anything `as_loci` takes; a list of peak sets is merged)
   gives the candidates;
2. with `bw`, the distinct windows are scored by one `signal(..., n_bins=1,
   span=True)` call (max over the bigWigs), and the support rule is applied
   with a `maximum.at` per gene;
3. the winner per gene is the first row of a `lexsort` on (gene, −length,
   −score, transcript id), or (gene, −score, −length, ...) for
   `rank="signal"`; genes with no supported isoform fall back to all of
   theirs, so every gene keeps a canonical transcript and no transcript is
   removed.

The evidence stays in the tables: `transcripts['tss_score']` (NaN = not
scored), `transcripts['tss_support']`, `genes['canonical']` (transcript row,
-1 = none); `collapse=True` rewrites the gene's `starts` / `ends` and drops
the annotation index.

```python
genes.select_isoforms(gb.as_loci([("chr1", 1900, 2100)]), r=200, verbose=False)
G.cols["canonical"], T.cols["tss_support"]
# -> (array([1, 2, 3]), array([False,  True, False, False]))
G.to_records()[0]                               # GENE_A's body moved onto T1b
# -> ('chr1', 2000, 5000, '+')
genes.representative(), genes.representative(canonical=False)
# -> (array([1, 2, 3]), array([0, 2, 3]))
```

`representative()` (canonical when selected, else longest) is what
`to_bed12`, the browser's collapsed gene track, `View` and `igv_html` draw.

## Exports and persistence
{: .sec-navy }

| Call | Gives |
|---|---|
| `to_pandas(table)`, `to_polars(table)`, `to_arrow(table)` | one of the three tables with its links spelled out (`transcripts` gains `gene_id` / `gene_name`, `features` gains `feature`, `transcript_id`, `gene_id`) |
| `to_gtf(path)` | a GTF again, 1-based (`start + 1`), genes then transcripts then features |
| `to_bed12(path=None, canonical=True)` | one BED12 line per gene from its representative transcript; exons as blocks, CDS as the thick part |
| `save(dir)` / `Genes.load(dir)` | three parquet files plus `meta.json` (`promoter_r`, coordinates) |
| `shape`, `columns`, `head()`, `describe()`, protocols | the genes table, through `TableMixin` |

```python
genes.to_bed12().splitlines()[0]
# -> 'chr1\t2000\t5000\tGENE_A\t0\t+\t5000\t5000\t0\t1\t3000\t0'
genes.save("genes_tables"); gb.Genes.load("genes_tables").genes.equals(G, cols=True)
# -> True
import polars as pl
pl.DataFrame(genes).shape
# -> (3, 8)
```

## Costs
{: .sec-navy }

- **Parsing** is one polars (or pandas) read with regex attribute extraction
  and three joins; the Python side only builds the `Loci` from finished
  columns, so a GENCODE GTF takes seconds, not minutes.
- **The annotation index** is five `merge` calls, built once per `Genes`.
- **Labelling** is five `overlap_any` calls over all CREs at once; nearest
  genes is one `nearest`.
- **Isoform selection** is one `overlap_any` over all transcripts plus one
  bigWig summary per distinct TSS window when `bw` is given.
- **Memory:** 21 bytes per row per table for the coordinates plus the id
  columns (object arrays), shown by `_repr_html_`.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
