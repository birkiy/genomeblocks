---
title: interop
parent: API Reference
layout: default
nav_order: 10
---

# `genomeblocks.interop`
{: .no_toc }

The boundary between genomeblocks and the rest of the Python ecosystem.
Anything interval-like comes in through [`as_loci`](#as_loci) — which every
public function calls on its inputs — and every table goes back out through
the `to_*` converters on the containers, all of which are one-line wrappers
of the functions on this page. Frames are read as 0-based half-open
intervals (BED / bioframe / pyranges convention) and row order is always
kept, so a frame's rows, an AnnData's `var` and the resulting `Loci` line up
one to one. See [Interoperability]({{ '/interoperability/' | relative_url }})
for the overview and [Design: tables]({{ '/design/tables/' | relative_url }})
for the protocols.
{: .fs-5 .fw-300 }

```python
import numpy as np, pandas as pd
import genomeblocks as gb
from genomeblocks import as_loci                 # re-exported at the top level
from genomeblocks import interop
```

| Library / format | In | Out |
|---|---|---|
| BED / narrowPeak / broadPeak | `Loci.make(path)`, `as_loci(path)` | `L.to_bed(path)` |
| CSV / TSV with a header, parquet | `as_loci(path)`, `Loci.load(path)` | `L.save(path)` |
| pandas, polars, pyarrow | `Loci.from_frame(df)`, `as_loci(df)` | `L.to_pandas()`, `L.to_polars()`, `L.to_arrow()` |
| bioframe | `Loci.from_frame(df)` | `L.to_bioframe()` |
| pyranges | `Loci.from_pyranges(gr)` | `L.to_pyranges()` |
| pybedtools | `Loci.from_bedtool(bt)` | `L.to_bedtool()` |
| cgranges | — | `L.to_cgranges()` |
| AnnData (scATAC) | `Loci.from_anndata(adata)` | `L.to_anndata(X)`, `cube_to_anndata(S, L)` |
| xarray | — | `cube_to_xarray(S, L)` |
| Biopython | FASTA handles everywhere | `L.to_seqrecords(fasta)`, `write_fasta` |
| UCSC chain files | — | `liftover(L, chain)` |

The examples use the seven-row `peaks.bed` (`cre`, columns `name` and
`score`), a synthetic `genome.fa` and two bigWigs.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The boundary
{: .sec-navy }

### `as_loci`

```python
as_loci(x, *, genome=None, keep=True) -> Loci
```

Anything interval-like as a `Loci`:

| Input | Read as |
|---|---|
| a `Loci` | Returned as is (re-coded onto `genome` when one is given and differs). |
| a `Locus` | A one-row `Loci`. |
| a path | BED / narrowPeak / broadPeak / bedGraph (`Loci.make`, rows in **file order**); `.csv` / `.tsv` / `.txt` with a header row (a UCSC `#chrom chromStart ...` header counts) or header-less with BED column names; `.parquet` / `.pq` (`Loci.load`); `.gz` everywhere. |
| a region string | `'chr1:1-2,000'`, `'chr8:127.7-128.1 Mb'` — anything [`parse_region`]({{ '/api/locus/' | relative_url }}) reads. |
| a table | pandas, polars (eager or lazy), pyarrow, bioframe or PyRanges frames, a `BedTool`, a dict of equal-length columns, a numpy structured array, any object with `__arrow_c_stream__` or `__dataframe__` — through [`frame`](#frame) and [`loci_from_frame`](#loci_from_frame). |
| an AnnData | Its `var` (`loci_from_anndata`). |
| a `(chroms, starts, ends[, strands])` tuple of arrays | Columns. |
| a list | Of `Locus` objects, `(chrom, start, end[, strand])` tuples, region strings or uids → `Loci.from_records`; or of any of the above, concatenated in order. |

`keep=True` keeps a table's other columns (a list keeps those). A path that
does not exist and is not a region raises `FileNotFoundError`; anything
else unreadable raises `TypeError` / `ValueError` naming what was expected.

```python
cre = gb.Loci.make("peaks.bed", keep=True)
df = cre.to_pandas()
import polars as pl
for x in (cre, df, pl.from_pandas(df), pl.from_pandas(df).lazy(), cre.to_arrow(), cre.to_bioframe(),
          cre.to_pyranges(), cre.to_bedtool(), cre.to_numpy(), cre.to_records(), cre.uid.tolist(),
          cre.to_anndata(np.zeros((2, 7)))):
    assert as_loci(x).equals(cre)                       # twelve spellings of the same seven rows
as_loci({"chrom": cre.chroms, "start": cre.starts, "end": cre.ends})
# -> Loci(n=7, chroms=2)
as_loci((cre.chroms, cre.starts, cre.ends)).equals(cre)       # no strand column: strands differ
# -> False
```

```python
as_loci("chr1:1,000-2,000")
# -> Loci(n=1, chroms=1)
as_loci(gb.Locus("chr1", 5, 9, "-")).to_records()
# -> [('chr1', 5, 9, '-')]
as_loci("peaks.bed").is_sorted, as_loci("cre.parquet")
# -> (True, Loci(n=7, chroms=2, cols=[name, score]))
as_loci([cre.head(2), "chr3:1-2", ("chr4", 1, 5)]).to_records()
# -> [('chr1', 900, 1100, '+'), ('chr1', 1900, 2100, '-'), ('chr3', 1, 2, '.'), ('chr4', 1, 5, '.')]
as_loci([]), as_loci(["chr1:1-2", "chr1:5-9"]).to_records()
# -> (Loci(n=0, chroms=0), [('chr1', 1, 2, '.'), ('chr1', 5, 9, '.')])
as_loci(df, keep=["name"]).columns
# -> ['chrom', 'start', 'end', 'strand', 'name']
g = cre.genome
as_loci(df, genome=g).genome is g, as_loci(cre, genome=gb.Genome()).genome is g
# -> (True, False)
```

Text files with a header, and the column names `from_frame` recognises:

```python
open("u.tsv", "w").write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\nchr1\t5\t15\tb\t0\t-\n")
as_loci("u.tsv").columns, as_loci("u.tsv").strand
# -> (['chrom', 'start', 'end', 'strand', 'name', 'score'], array(['-'], dtype=object))
open("c.csv", "w").write("chrom,start,end\nchr2,100,200\nchr1,1,2\n")
as_loci("c.csv").chroms                                   # file order, not sorted
# -> array(['chr2', 'chr1'], dtype=object)
for cols in (("Chromosome", "Start", "End"), ("seqnames", "start", "end"), ("#chrom", "chromStart", "chromEnd"), ("chr", "begin", "stop")):
    as_loci(pd.DataFrame({cols[0]: ["chr1"], cols[1]: [5], cols[2]: [9]})).to_records()
# -> [('chr1', 5, 9, '.')]      (each of the four)
as_loci(pd.DataFrame([["chr1", 5, 9, "x"]])).to_records()  # header-less: text + int + int by position
# -> [('chr1', 5, 9, '.')]
```

The errors name the expected input and the fix:

```python
as_loci({"a": 1})
# -> TypeError: a dict of columns needs an array per key, but ['a'] hold scalars; expected e.g. {'chrom': [...], 'start': [...], 'end': [...]}
as_loci({"chrom": ["chr1", "chr1"], "start": [1], "end": [2, 3]})
# -> ValueError: the columns have different lengths: {'chrom': 2, 'start': 1, 'end': 2}
as_loci(pd.Series([1, 2]))
# -> TypeError: cannot read intervals from a Series: pass a DataFrame / Table with chrom, start and end columns
as_loci("no_such_file.bed")
# -> FileNotFoundError: no such file: 'no_such_file.bed' (and not a region like 'chr1:1,000-2,000')
as_loci(pd.DataFrame({"a": ["chr1"], "b": ["x"], "c": [1]}))
# -> ValueError: expected columns for the chromosome, start and end, e.g. chrom/start/end (bioframe), Chromosome/Start/End (pyranges), seqnames/start/end or chr/chromStart/chromEnd — or pass chrom=, start=, end= to name them; got columns ['a', 'b', 'c']
as_loci(object())
# -> TypeError: cannot read intervals from object: pass a Loci, a BED / CSV / parquet path, a region string, a pandas / polars / pyarrow frame, or a list of regions
```

{: .note }
`as_loci(path)` keeps BED rows in file order, like every other input, while
`Loci.make(path)` sorts into genome order by default. Call `.sort()` when
you need `chrom_offsets` or a genome-ordered table.

### `frame`

```python
frame(obj) -> narwhals.DataFrame | None
```

An eager [narwhals](https://narwhals-dev.github.io/narwhals/) DataFrame for
anything table-like, else `None`: pandas, polars (lazy frames are collected),
pyarrow, modin, cuDF, DuckDB and the other frames narwhals knows; PyRanges
and pybedtools objects (through their pandas view); dicts of equal-length
columns; numpy structured arrays; and any object that speaks the Arrow C
stream (`__arrow_c_stream__`) or dataframe interchange (`__dataframe__`)
protocol. This is the one place genomeblocks branches on input type; from
here on everything is a narwhals frame.

```python
f = interop.frame(df)
type(f).__name__, f.shape
# -> ('DataFrame', (7, 6))
interop.frame(pl.from_pandas(df).lazy()).shape, interop.frame(cre.to_pyranges()).columns
# -> ((7, 6), ['Chromosome', 'Start', 'End', 'Strand', 'name', 'score'])
interop.frame(cre.to_numpy()).columns, interop.frame({"a": [1, 2]}).shape, interop.frame(object())
# -> (['chrom', 'start', 'end', 'strand', 'name', 'score'], (2, 1), None)
```

---

## Into Loci
{: .sec-green }

| Function | Signature | One line |
|---|---|---|
| `loci_from_frame` | `loci_from_frame(df, chrom=None, start=None, end=None, strand=None, *, keep=True, genome=None)` | `Loci.from_frame`: any table `frame` accepts, columns found by name or given, row order kept. |
| `loci_from_pyranges` | `loci_from_pyranges(gr, **kw)` | `Loci.from_pyranges` (pyranges 0.x or 1.x). |
| `loci_from_bedtool` | `loci_from_bedtool(bt, **kw)` | `Loci.from_bedtool`. |
| `loci_from_anndata` | `loci_from_anndata(adata, axis="var", *, genome=None, keep=True)` | `Loci.from_anndata`: the regions of one axis, from chrom/start/end columns when present, else parsed from names like `chr1:100-200`, `chr1-100-200`, `chr1_100_200`. |

### `loci_from_frame`

Column names are matched case-insensitively: chromosome `chrom /
chromosome / chr / seqnames / seqname / seqid / contig / #chrom / ...`,
start `start / chromStart / begin / start_position / txStart`, end `end /
chromEnd / stop / end_position / txEnd`, strand `strand` (`+`, `-`, `1`,
`-1`; anything else is `.`). `chrom=` / `start=` / `end=` / `strand=` name
them explicitly. A frame with none of these but a text column followed by
two integer columns (`read_csv(header=None)` on a BED) is read by position.
Missing values in the three coordinate columns raise.

```python
L = interop.loci_from_frame(pd.DataFrame({"c": ["chr1"], "s": [1], "e": [9], "x": [2.0]}), chrom="c", start="s", end="e")
L, L.cols
# -> (Loci(n=1, chroms=1, cols=[x]), {'x': array([2.])})
interop.loci_from_frame(pd.DataFrame({"chrom": ["chr1", None], "start": [1, 2], "end": [5, 6]}))
# -> ValueError: column 'chrom' has 1 missing value(s) in 2 rows; chrom, start and end must be complete — drop or fill those rows first
interop.loci_from_frame(pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [9]}), chrom="nope")
# -> KeyError: "column 'nope' not found; columns are ['chrom', 'start', 'end']"
interop.loci_from_pyranges(cre.to_pyranges()).equals(cre), interop.loci_from_bedtool(cre.to_bedtool()).columns
# -> (True, ['chrom', 'start', 'end', 'strand', 'name', 'score'])
```

### `loci_from_anndata`

```python
import anndata as ad
A = ad.AnnData(np.zeros((2, 3)), var=pd.DataFrame(index=["chr1:1-2", "chr1_3_4", "chr2-5-6"]))
interop.loci_from_anndata(A).to_records()
# -> [('chr1', 1, 2, '.'), ('chr1', 3, 4, '.'), ('chr2', 5, 6, '.')]
A2 = ad.AnnData(np.zeros((3, 2)), obs=pd.DataFrame({"chrom": ["chr1"] * 3, "start": [1, 2, 3], "end": [5, 6, 7]}, index=list("abc")))
interop.loci_from_anndata(A2, axis="obs").to_records()
# -> [('chr1', 1, 5, '.'), ('chr1', 2, 6, '.'), ('chr1', 3, 7, '.')]
interop.loci_from_anndata(ad.AnnData(np.zeros((1, 1)), var=pd.DataFrame(index=["peak1"])))
# -> ValueError: cannot read a region from the var name 'peak1': expected chr1:100-200, chr1-100-200 or chr1_100_200 names, or chrom/start/end columns in .var
```

---

## Out of Loci
{: .sec-green }

Each of these is what the `Loci.to_*` method of the same name calls.
`chrom` goes out as a categorical / dictionary array ordered by the genome's
natural order and `strand` as a categorical, so pandas, polars and Arrow
keep the small-integer representation; bioframe and pyranges get plain
strings, which is what they expect.

| Function | Returns | One line |
|---|---|---|
| `loci_to_pandas(L, uid=False)` | `DataFrame` | `chrom` (categorical), `start`, `end`, `strand` (categorical), the extra columns, `uid` when asked. |
| `loci_to_arrow(L)` | `pyarrow.Table` | Dictionary-encoded `chrom` / `strand`; NaN / None in object columns become Arrow nulls. The source of `__arrow_c_stream__`, `__dataframe__`, `to_polars` and `save`. |
| `loci_to_polars(L)` | `polars.DataFrame` | `pl.from_arrow(loci_to_arrow(L))`. |
| `loci_to_bioframe(L)` | `DataFrame` | Plain `chrom / start / end / strand` strings + extra columns. |
| `loci_to_pyranges(L)` | `PyRanges` | `Chromosome / Start / End / Strand` + extra columns; `.` strands are kept (pyranges then calls the object unstranded). |
| `loci_to_bedtool(L)` | `BedTool` | BED6: `name` = `cols['name']` or the uid, `score` = `cols['score']` or 0. |
| `loci_to_cgranges(L)` | `cgranges` index | Built, label = row. |
| `loci_to_anndata(L, X=None, *, obs=None, layers=None, **kw)` | `AnnData` | Loci as `var` (index `chrom:start-end`); `X` is `(n_obs, n_loci)`; `obs` a DataFrame or a list of names; `**kw` go to `AnnData`. |

```python
interop.loci_to_pandas(cre).dtypes.to_dict()
# -> {'chrom': CategoricalDtype(categories=['chr1', 'chr2'], ordered=False, categories_dtype=object), 'start': dtype('int64'), 'end': dtype('int64'),
#     'strand': CategoricalDtype(categories=['.', '+', '-'], ordered=False, categories_dtype=object), 'name': dtype('O'), 'score': dtype('float64')}
interop.loci_to_arrow(cre).schema
# -> chrom: dictionary<values=string, indices=int32, ordered=0>
#    start: int64
#    end: int64
#    strand: dictionary<values=string, indices=int8, ordered=0>
#    name: string
#    score: double
interop.loci_to_bioframe(cre).dtypes.to_dict()
# -> {'chrom': dtype('O'), 'start': dtype('int64'), 'end': dtype('int64'), 'strand': dtype('O'), 'name': dtype('O'), 'score': dtype('float64')}
import bioframe as bf
bf.cluster(cre.to_bioframe(), min_dist=1000).shape            # straight into bioframe
# -> (7, 9)
str(interop.loci_to_bedtool(cre)).splitlines()[0]
# -> 'chr1\t900\t1100\tp1\t10.0\t+'
a = interop.loci_to_anndata(cre, np.ones((2, 7)), obs=["s1", "s2"])
a, a.var.index[:2].tolist()
# -> (AnnData object with n_obs × n_vars = 2 × 7
#         var: 'chrom', 'start', 'end', 'strand', 'name', 'score', ['chr1:900-1100', 'chr1:1900-2100'])
interop.loci_to_anndata(cre, np.ones((2, 5)))
# -> ValueError: X has 5 columns for 7 loci (X is obs x loci)
interop.loci_to_cgranges(cre)                               # without cgranges
# -> ImportError: No module named 'cgranges'
```

---

## Signal cubes
{: .sec-purple }

[`Loci.signal`]({{ '/api/signal/' | relative_url }}) returns a `(rows,
tracks, bins)` cube aligned to the loci. These three functions attach the
loci (and the track names) to it for the tools that want labelled arrays.

| Function | Signature | Returns |
|---|---|---|
| `cube_to_xarray` | `cube_to_xarray(S, L, tracks=None, *, flank=None, name="signal")` | `xarray.DataArray` with dims `(region, track, bin)`; coordinates `region` (names), `chrom` / `start` / `end` along `region`, `track`, and `bin` as bin centres in bp relative to the locus centre when `flank` is given, else bin indices. |
| `cube_to_anndata` | `cube_to_anndata(S, L, tracks=None, *, agg="mean")` | `AnnData` with tracks as `obs`, loci as `var`, `X` the per-locus `agg` (`mean`, `sum`, `max`) over bins; with more than one bin each track's profile is kept in `varm['bins:<track>']` (loci x bins). |
| `cube_to_pandas` | `cube_to_pandas(S, L, tracks=None)` | `DataFrame` indexed by region name: one column per track for a single bin, else `(track, bin)` MultiIndex columns. |

`tracks` defaults to `track_0, track_1, ...`.

```python
S = cre.signal(["signal.bw", "signal2.bw"], n_bins=5, flank=250, verbose=False, progress=False)
da = interop.cube_to_xarray(S, cre, ["ATAC", "H3K27ac"], flank=250)
da.dims, da.shape, da.coords["bin"].values, da.name
# -> (('region', 'track', 'bin'), (7, 2, 5), array([-200., -100.,    0.,  100.,  200.]), 'signal')
da.sel(track="ATAC").mean("region").values
# -> array([3.102884 , 5.5290074, 4.218671 , 1.5205952, 5.7615027], dtype=float32)
interop.cube_to_xarray(S, cre).coords["bin"].values, interop.cube_to_xarray(S, cre).coords["track"].values
# -> (array([0, 1, 2, 3, 4]), array(['track_0', 'track_1'], dtype='<U7'))
```

```python
A3 = interop.cube_to_anndata(S, cre, ["ATAC", "H3K27ac"])
A3
# -> AnnData object with n_obs × n_vars = 2 × 7
#        var: 'chrom', 'start', 'end', 'strand', 'name', 'score'
#        varm: 'bins:ATAC', 'bins:H3K27ac'
A3.obs_names.tolist(), A3.X.shape
# -> (['ATAC', 'H3K27ac'], (2, 7))
interop.cube_to_anndata(S, cre, agg="max").X[:, :2]
# -> array([[9.415652 , 7.0356197],
#           [7.6968775, 9.237916 ]], dtype=float32)
d = interop.cube_to_pandas(S, cre, ["ATAC", "H3K27ac"])
d.shape, d.columns[:3].tolist(), d.index.name
# -> ((7, 10), [('ATAC', 0), ('ATAC', 1), ('ATAC', 2)], 'region')
interop.cube_to_pandas(S[:, :, :1], cre, ["ATAC", "H3K27ac"]).head(2)
# ->                     ATAC   H3K27ac
#    region
#    chr1:900-1100   5.892625  1.650387
#    chr1:1900-2100  3.697407  5.250703
```

---

## Sequences
{: .sec-purple }

| Function | Signature | One line |
|---|---|---|
| `sequences` | `sequences(L, fasta, r=None, *, strand=False, upper=False, backend=None) -> list[str]` | Sequence of every row, or of `center ± r`, through the [fasta backend]({{ '/api/backends/' | relative_url }}). `fasta` is a path (indexed on first use), a `{chrom: str}` dict or an open pysam / pyfaidx / Biopython handle. Windows are clipped to the chromosome; rows on a chromosome the FASTA lacks give `''`; `strand=True` reverse-complements `-` rows; `upper=True` upper-cases. |
| `windows` | `windows(L, r=None) -> (starts, ends)` | The windows `sequences` reads: the rows themselves, or `center - r`, `center + r`. |
| `revcomp` | `revcomp(s) -> str` | Reverse complement; IUPAC codes and case are kept. |
| `to_seqrecords` | `to_seqrecords(L, fasta, r=None, *, strand=False, **kw) -> list[SeqRecord]` | Biopython records, `id` = uid, empty description. |
| `write_fasta` | `write_fasta(L, path, fasta, r=None, *, strand=False, width=60, **kw) -> path` | Write them as FASTA (`>uid` headers, `width` bases per line). |

These are `Loci.sequences`, `Loci.to_seqrecords` and `Loci.to_fasta`.

```python
interop.revcomp("ACGTN"), interop.revcomp("acgtRY")
# -> ('NACGT', 'RYacgt')
interop.windows(cre.head(2)), interop.windows(cre.head(2), r=10)
# -> ((array([ 900, 1900]), array([1100, 2100])), (array([ 990, 1990]), array([1010, 2010])))
interop.sequences(cre, "genome.fa", r=5)
# -> ['TCAGCACGTA', 'TGATGTTGTG', 'CGTTAGGAAG', 'TCTCAACCGG', 'GCGCCGGGTC', 'CGTAATCTCC', 'TTGGCTGATT']
interop.sequences(as_loci([("chr1", 1000, 1012, "-")]), "genome.fa", strand=True)
# -> ['ACGTACGTACGT']
interop.sequences(as_loci([("chr1", 19_990, 20_050), ("chr9", 0, 10)]), "genome.fa")   # clipped / missing chromosome
# -> ['AGACCATCTG', '']
interop.sequences(cre.head(2), {"chr1": "N" * 20000}, r=2)
# -> ['NNNN', 'NNNN']
interop.sequences(cre.head(1), "genome.fa", backend="pysam") == interop.sequences(cre.head(1), "genome.fa")
# -> True
recs = interop.to_seqrecords(cre, "genome.fa", r=10)
recs[0].id, str(recs[0].seq)
# -> ('chr1:900-1100(+)', 'GCTGATCAGCACGTACGTAC')
interop.write_fasta(cre, "out.fa", "genome.fa", r=10, width=30)
# -> 'out.fa'
```

---

## Liftover
{: .sec-purple }

```python
liftover(L, chain_file, *, min_match=0.95, verbose=True) -> Loci
```

UCSC-style liftover through `pyliftover` (`pip install pyliftover`; the
`LiftOver` object is cached per chain file). Each end that falls in a chain
gap walks inward by up to `1 - min_match` of the interval length before
giving up (UCSC's base-fraction rule); a row is kept when both ends land on
one chromosome and strand and the walked bases stay within budget. Rows
that do not lift are dropped; `cols['source_row']` says where each output
row came from, so columns can be carried over with `L.cols[k][out['source_row']]`.
The result lives on `L`'s `Genome` (new chromosome names are added). This is
`Loci.liftover`.

```python
# a chain that maps chr1 onto chrA shifted by +10,000 and says nothing about chr2
open("t.chain", "w").write("chain 1000 chr1 20000 + 0 20000 chrA 30000 + 10000 30000 1\n20000\n\n")
out = interop.liftover(cre, "t.chain")
# -> [INFO] liftover: 7 → 5 (2 dropped, min_match=0.95)
out, out.to_records()[:2], out["source_row"]
# -> (Loci(n=5, chroms=1, cols=[source_row]), [('chrA', 10900, 11100, '+'), ('chrA', 11900, 12100, '-')], array([0, 1, 2, 3, 4]))
out.genome is cre.genome
# -> True
interop.liftover(cre, "t.chain")                          # without pyliftover
# -> ImportError: Loci.liftover needs pyliftover: pip install pyliftover
```
