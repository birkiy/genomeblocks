---
title: Interoperability
layout: default
nav_order: 4.5
permalink: /interoperability/
---

# Interoperability
{: .no_toc }

genomeblocks tables go in and out of the Python ecosystem without glue code.
Anything interval-like goes in through one function, `gb.as_loci`, which every
public function calls on its inputs. Every container converts back out with an
explicit `to_*` method, and speaks the Arrow and dataframe protocols so polars,
DuckDB, seaborn, plotly and altair take it as it is.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Two rules
{: .sec-navy }

**Compatibility is explicit and visible.** Inputs are lenient: a pandas or
polars frame, a pyarrow Table, a bioframe or PyRanges object, a dict of
columns, a file path or a list of region strings all work, and column names
may be `chrom` / `chr` / `Chromosome` / `seqnames`. Internally there is one
representation: numpy columns on one `Genome`, 0-based half-open, strict names.
Outputs are explicit: `to_pandas()`, `to_polars()`, `to_arrow()` and the domain
converters on every container.

**Row order is kept.** A frame's rows, an AnnData's `var`, a BED file's lines
and the resulting `Loci` line up one to one. Only `Loci.make` sorts (into
genome order, because the index and the Architecture need it); `as_loci(path)`
keeps file order like every other input.

```python
import numpy as np
import pandas as pd
import genomeblocks as gb

df = pd.DataFrame({"chrom": ["chr1", "chr1", "chr2"], "start": [900, 1900, 500],
                   "end": [1100, 2100, 600], "name": ["p1", "p2", "p6"],
                   "score": [10.0, 20.0, 60.0]})
L = gb.as_loci(df)
L                    # -> Loci(n=3, chroms=2, cols=[name, score])
L.columns            # -> ['chrom', 'start', 'end', 'strand', 'name', 'score']
L.to_polars()        # the same rows, chrom and strand as categoricals
```

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/interop-path.svg %}
</div><figcaption>
<strong>In through narwhals, out through Arrow.</strong> <code>interop.frame</code> turns any frame — pandas, polars, pyarrow, DuckDB, bioframe, PyRanges, a BedTool, a dict or a structured array — into one eager narwhals DataFrame; <code>loci_from_frame</code> finds the chromosome, start, end and strand columns by name and casts them into the numpy columns of a <code>Loci</code>. On the way out, <code>to_arrow</code> wraps the numeric columns without copying, and the Arrow C stream, interchange and narwhals protocols give polars, DuckDB, seaborn, plotly and altair that table as it is.
</figcaption></figure>

---

## In: `as_loci` takes everything
{: .sec-green }

`gb.as_loci(x)` is the boundary. It accepts:

| Input | How it is read |
|---|---|
| a `Loci` | returned as is (re-coded onto `genome=` if one is given) |
| a `Locus`, a region string `"chr1:1,000-2,000"`, a uid `"chr1:900-1100(+)"` | one row |
| a list of regions, tuples `(chrom, start, end[, strand])` or `Locus` objects | one row each, in order |
| a list of any of the above (frames, paths, Loci) | concatenated in order |
| pandas DataFrame, polars DataFrame or LazyFrame, pyarrow Table | columns found by name (see below) |
| bioframe frame, PyRanges (0.x or 1.x), pybedtools `BedTool` | the same, through their own frame |
| a dict of equal-length columns, a numpy structured array | the same |
| a `(chroms, starts, ends[, strands])` tuple of arrays | by position |
| any object with `__arrow_c_stream__` or `__dataframe__` (DuckDB relations, modin, cuDF, ...) | through narwhals |
| an AnnData | its `var` (scATAC peaks): chrom/start/end columns, else parsed from the names |
| a path to BED / narrowPeak / broadPeak / bedGraph (`.gz` too) | file order kept |
| a path to CSV / TSV / `.txt` | header row or BED-like header-less file, see [Files](#files) |
| a path to parquet | what `Loci.save` wrote |

```python
import polars as pl, pyarrow as pa
for x in (pl.from_pandas(df), pl.from_pandas(df).lazy(), pa.Table.from_pandas(df),
          df.to_dict("list"), df.to_records(index=False)):
    assert gb.as_loci(x).equals(L, cols=True)

gb.as_loci([("chr1", 900, 1100), "chr1:1,900-2,100", L[2]]).to_records()
# -> [('chr1', 900, 1100, '.'), ('chr1', 1900, 2100, '.'), ('chr2', 500, 600, '.')]

import duckdb
gb.as_loci(duckdb.sql("select * from df where score > 15"))     # a DuckDB relation
gb.as_loci(pl.scan_parquet("peaks.parquet").filter(pl.col("score") > 25))   # lazy polars
```

`Loci.from_frame` / `from_pandas` / `from_polars` / `from_arrow` /
`from_bioframe` / `from_pyranges` / `from_bedtool` / `from_anndata` /
`from_records` / `from_uids` are the same readers with a name, for code that
wants to say what it expects. `keep=True` (the default) keeps every other
column; a list keeps those; `keep=False` keeps only the coordinates.

### Column names: lenient on input, strict internally

Names are matched case-insensitively against these spellings:

| Field | Accepted names |
|---|---|
| chromosome | `chrom`, `chromosome`, `chr`, `seqnames`, `seqname`, `seqid`, `contig`, `#chrom`, `chrom_name`, `chromosome_name`, `ref`, `reference` |
| start | `start`, `chromStart`, `begin`, `start_position`, `txStart` |
| end | `end`, `chromEnd`, `stop`, `end_position`, `txEnd` |
| strand | `strand` (values `+`, `-`, `.`, `*`, `1`, `-1`; anything else is `.`) |

So bioframe (`chrom/start/end`), PyRanges (`Chromosome/Start/End/Strand`), UCSC
(`#chrom/chromStart/chromEnd`), R (`seqnames`) and GTF-style (`seqname`)
frames all read without renaming. A frame with none of these names but a text
column and two integer columns first (`pd.read_csv(path, header=None)` on a
BED) is read by position. Anything else raises and names the accepted
spellings; `chrom=`, `start=`, `end=`, `strand=` name the columns explicitly.

```python
for cols in (("Chromosome", "Start", "End"), ("seqnames", "start", "end"),
             ("#chrom", "chromStart", "chromEnd"), ("chr", "begin", "stop")):
    d = pd.DataFrame({cols[0]: ["chr1"], cols[1]: [5], cols[2]: [9]})
    assert gb.as_loci(d).to_records() == [("chr1", 5, 9, ".")]

gb.Loci.from_frame(pd.DataFrame({"contig": ["chr1"], "lo": [5], "hi": [9]}), start="lo", end="hi")

gb.as_loci(pd.DataFrame({"a": ["chr1"], "b": ["x"], "c": [1]}))
# ValueError: expected columns for the chromosome, start and end, e.g. chrom/start/end
#   (bioframe), Chromosome/Start/End (pyranges), seqnames/start/end or
#   chr/chromStart/chromEnd — or pass chrom=, start=, end= to name them; got columns ['a', 'b', 'c']
```

Inside, the columns are always `chrom`, `start`, `end`, `strand` (plus the
kept columns under their own names), and every output uses those names.
Coordinates are 0-based half-open, the BED / bioframe / pyranges convention; a
frame with GTF-style 1-based starts goes to `Genes.from_frame`, which converts.

{: .warning }
> Missing chromosomes or coordinates are an error, not a silent drop:
> `column 'chrom' has 1 missing value(s) in 2 rows; chrom, start and end must be
> complete — drop or fill those rows first`.

### Files

| Format | Read with | Notes |
|---|---|---|
| BED, narrowPeak, broadPeak, bedGraph (`.gz`) | `Loci.make(path)` (sorted), `as_loci(path)` (file order) | `#`, `track`, `browser` lines skipped; `keep=True` names the standard columns (`name`, `score`, `signalValue`, `pValue`, `qValue`, `peak`) |
| CSV / TSV / `.txt` with a header row | `as_loci(path)` | a UCSC `#chrom chromStart ...` header counts as a header |
| header-less `.tsv` / `.txt` with BED columns | `as_loci(path)` | columns get the BED names |
| parquet (`Loci.save`) | `Loci.load(path)`, `as_loci(path)` | chrom / strand are dictionary columns |
| BEDPE (`.gz`) | `Pairs.make(path)` | six anchor columns, then name, score, strand1, strand2 |
| GTF / GFF3 (`.gz`), GENCODE / Ensembl / RefSeq | `Genes.make(path)` | 1-based starts become 0-based |
| UCSC genePred (refGene, refFlat, knownGene) | `Genes.make_ucsc(path)` | already 0-based |
| 4DN `.pairs`, HiC-Pro allValidPairs, Juicer medium | `genomeblocks.bedpe.count_pairs`, `count_pairs_2d`, `read_pairs_chunks` | format detected from the first data line |
| FASTA (`.fai` built when missing; `.gz` read into memory) | `Loci.sequences`, `gb.read_fasta`, `Genome.from_fasta` | see [Backends]({{ '/backends/' | relative_url }}) for the readers |
| bigWig | `Loci.signal`, `genomeblocks.browser` | paths or open handles |
| `.chrom.sizes` | `Genome.from_sizes`, `Loci.tile_genome` | also a dict, a Series or a cooler |
| `.mcool` / `.cool` | `Architecture.add_mcool` | through cooler |
| BAM | `gb.coverage`, `genomeblocks.browser` | needs pysam (`pip install "genomeblocks[bam]"`) |
| UCSC chain | `Loci.liftover` | needs pyliftover |
| JASPAR, jaspar16, TRANSFAC, uniprobe, MEME | `gb.load_motifs(path, format=...)` | a `Library` |

```python
cre = gb.as_loci("peaks.bed")              # file order
cre = gb.Loci.make("peaks.bed", keep=True) # genome order, name and score kept
gb.as_loci("regions.csv")                  # chrom,start,end header
gb.as_loci("peaks.parquet")                # what cre.save() wrote
```

---

## Out: every container converts
{: .sec-green }

| Container | Table it exposes | Converters |
|---|---|---|
| `Loci` | chrom, start, end, strand + extra columns | `to_pandas(uid=False)`, `to_polars`, `to_arrow`, `to_bioframe`, `to_pyranges`, `to_bedtool`, `to_cgranges`, `to_anndata(X, obs=)`, `to_records`, `to_numpy`, `to_bed(path, name=, score=)`, `save(path)` (parquet), `to_seqrecords(fasta)`, `to_fasta(path, fasta)` |
| `Genes` | the genes table (`to_pandas("transcripts")` / `"features"` for the others, links spelled out) | `to_pandas(table)`, `to_polars(table)`, `to_arrow(table)`, `to_gtf(path)`, `to_bed12(path)`, `save(dir)`; `genes` / `transcripts` / `features` are plain `Loci` |
| `Pairs` | the BEDPE table | `to_pandas`, `to_polars`, `to_arrow`, `to_numpy`, `to_bedpe(path)`, `save(path)`; `P.a` / `P.b` are `Loci` |
| `Architecture` | the edge table (src, tgt, uid1, uid2, chrom1, chrom2, cis + edge columns) | `to_pandas` / `edges_frame`, `to_polars`, `to_arrow`, `vertices_frame`, `to_networkx`, `to_igraph`, `to_scipy`, `to_graph_tool`, `to_anndata`, `save(dir)` |
| `Atlas` | one row per track (name, n_peaks, n_bins + metadata) | `to_pandas`, `to_polars`, `to_arrow`, `save(path)` |
| `Library` (motifs) | name, description, width, consensus | `to_pandas`, `to_polars`, `to_arrow`, `to_biopython`, `to_moods`, `to_meme(path)`, `to_jaspar(path)` |

Every container also has `shape`, `columns`, `head()`, `tail()`,
`describe()`, a short `repr` and an HTML table in notebooks.

```python
L.to_pandas().dtypes
# chrom     category   (categories in genome order)
# start     int64
# end       int64
# strand    category   ('.', '+', '-')
# name      object
# score     float64
L.to_arrow().schema
# chrom: dictionary<values=string, indices=int32>, start: int64, end: int64,
# strand: dictionary<values=string, indices=int8>, name: string, score: double
L.to_bioframe().dtypes            # chrom and strand as plain object columns
L.to_bed(name="name", score="score")
# chr1	900	1100	p1	10.0	.
# chr1	1900	2100	p2	20.0	.
# chr2	500	600	p6	60.0	.
```

---

## Protocols: hand the table over without converting
{: .sec-green }

Every container implements `__arrow_c_stream__` (the Arrow PyCapsule
interface), `__dataframe__` (the dataframe interchange protocol) and
`__narwhals_dataframe__`; `Loci` also implements `__array__`, `__len__`,
`__iter__` and `__getitem__`. Libraries that look for these take the object
directly:

```python
import polars as pl, pyarrow as pa, duckdb, narwhals as nw
pl.DataFrame(L).shape                              # -> (3, 6)   Arrow C stream
pa.table(L).num_rows                               # -> 3
pd.api.interchange.from_dataframe(L).shape         # -> (3, 6)   interchange protocol
duckdb.sql("select chrom, count(*) n, avg(score) s from L group by chrom order by chrom").fetchall()
# -> [('chr1', 2, 15.0), ('chr2', 1, 60.0)]
nw.from_native(L, eager_only=True).shape           # -> (3, 6)
np.asarray(L)["start"]                             # -> array([ 900, 1900,  500])

import seaborn as sns, plotly.express as px, altair as alt
sns.scatterplot(data=L, x="start", y="score")
px.scatter(L, x="start", y="score", color="chrom")
alt.Chart(L).mark_point().encode(x="start:Q", y="score:Q")
```

The same holds for `Genes` (its genes table), `Pairs`, `Architecture` (its
edge table), `Atlas` (its track table) and `Library`:
`duckdb.sql("select count(*) from A where cis")` counts the cis edges of an
Architecture.

{: .note }
> altair 6 reads frames through narwhals and cannot consume
> interchange-only objects, which is why the containers also expose
> `__narwhals_dataframe__`.

---

## Library by library
{: .sec-navy }

| Library / object | In | Out |
|---|---|---|
| **pandas** DataFrame | `as_loci(df)`, `Loci.from_frame(df)`, `Genes.from_frame(gtf_df)`, `Pairs.from_frame(bedpe_df)`, `Architecture.from_frame(loci, edges_df)` | `to_pandas()` on every container; `Loci.to_pandas(uid=True)` adds the uid |
| **polars** DataFrame / LazyFrame | the same readers (a LazyFrame is collected once) | `to_polars()` (through Arrow; numeric columns are not copied) |
| **pyarrow** Table | the same readers, `Loci.from_arrow` | `to_arrow()`, `save()` (parquet) |
| **bioframe** | `Loci.from_bioframe(df)` (plain chrom/start/end frame) | `to_bioframe()`: chrom / strand as object columns, bioframe's names |
| **pyranges** 0.x / 1.x | `Loci.from_pyranges(gr)`, `Genes.from_frame(pr.read_gtf(path))` | `to_pyranges()`: Chromosome / Start / End / Strand + extra columns |
| **pybedtools** `BedTool` | `Loci.from_bedtool(bt)`: a BedTool on a BED-like file is read as `as_loci(bt.fn)` reads it (header lines skipped, file order kept); others through their frame | `to_bedtool()`: a BED6 file (name = `cols['name']` or the uid, score = `cols['score']` or 0) |
| **cgranges** | as an interval backend (`backend='cgranges'`) | `to_cgranges()`: a built index whose label is the row number |
| **AnnData** (scATAC) | `Loci.from_anndata(adata, axis='var')`: chrom/start/end columns, else names `chr1:100-200`, `chr1-100-200`, `chr1_100_200` | `to_anndata(X, obs=)` (loci as `var`, names `chrom:start-end`); `Architecture.to_anndata()` (vertices as `obs`, one adjacency per edge column in `obsp`); `interop.cube_to_anndata(S, L, tracks)` |
| **Biopython** | `SeqIO.index` / `SeqIO.to_dict` mappings as FASTA sources; `Bio.motifs` objects in `load_motifs` and every motif function; `backend='biopython'` for scanning and FASTA | `to_seqrecords(fasta)`, `to_fasta(path, fasta)`, `Library.to_biopython()` |
| **xarray** | — | `interop.cube_to_xarray(S, L, tracks, flank=)`: `(region, track, bin)` with chrom / start / end coordinates and bin centres in bp |
| **MOODS** | the default scanning engine (`backend='moods'`; `genomeblocks[motifs]` or conda `moods`) | `Library.to_moods()`: 4 x W log-odds lists |
| **lightmotif** | the alternative scanning engine (`backend='lightmotif'`; `genomeblocks[lightmotif]`); lightmotif motif objects in `load_motifs` | — |
| **pyBigWig / pybigtools** | open handles anywhere a bigWig path goes (`Loci.signal`, `browser`, `call_se`) | — (`backends.bigwig.open_bigwig` wraps a path or a handle) |
| **pysam / pyfaidx** | `pysam.FastaFile` / `pyfaidx.Fasta` handles as FASTA sources; pysam for BAM | — |
| **networkx / igraph / scipy / graph-tool** | `Architecture.from_scipy(loci, M)`; `backend=` for algorithms | `to_networkx()`, `to_igraph()`, `to_scipy(weight)`, `to_graph_tool()`; `A.graph(backend=)` |
| **duckdb, seaborn, plotly, altair, marimo** | — | take any container as is through the protocols |
| **R / anything that reads parquet** | — | `save()` on `Loci`, `Genes`, `Pairs`, `Architecture`, `Atlas` |

### pandas, polars, pyarrow

```python
G = gb.Genes.make("genes.gtf")
G.to_pandas("transcripts").columns.tolist()
# -> ['chrom', 'start', 'end', 'strand', 'transcript_id', 'gene', 'gene_id', 'gene_name']
P = gb.Pairs.make("loops.bedpe")
pl.DataFrame(P).shape                 # -> (4, 10)
A = gb.Architecture.make(cre, P, r=100, verbose=False)
A.to_pandas().columns.tolist()
# -> ['src', 'tgt', 'uid1', 'uid2', 'chrom1', 'chrom2', 'cis', 'w']
B = gb.Architecture.from_frame(cre, A.edges_frame())      # src/tgt rows, or uid1/uid2
list(B) == list(A)                    # -> True
```

### bioframe and pyranges

```python
import bioframe as bf, pyranges as pr
bf.merge(L.to_bioframe())             # any bioframe op on the frame
gb.Loci.from_bioframe(L.to_bioframe()).equals(L, cols=True)     # -> True

gr = L.to_pyranges()                  # PyRanges, extra columns included
gb.Loci.from_pyranges(gr).equals(L, cols=True)                  # -> True (L is sorted)
G2 = gb.Genes.from_frame(pr.read_gtf("genes.gtf"))             # 0-based Starts, as PyRanges stores them
G2.genes.starts.tolist() == G.genes.starts.tolist()             # -> True
```

### pybedtools

```python
import pybedtools
bt = L.to_bedtool()                   # a BED6 on disk (set pybedtools.set_tempdir first)
open(bt.fn).read().splitlines()[0]    # -> 'chr1\t900\t1100\tp1\t10.0\t.'
gb.Loci.from_bedtool(bt).equals(L)    # -> True
```

{: .note }
> `pybedtools` writes temporary files; point `pybedtools.set_tempdir()` (or
> `TMPDIR`) at a project directory before the first call.

### AnnData

```python
import numpy as np
ad = L.to_anndata(np.arange(6).reshape(2, 3), obs=["s1", "s2"])
# AnnData object with n_obs × n_vars = 2 × 3
#     var: 'chrom', 'start', 'end', 'strand', 'name', 'score'
ad.var.index.tolist()                 # -> ['chr1:900-1100', 'chr1:1900-2100', 'chr2:500-600']
gb.Loci.from_anndata(ad).equals(L)    # -> True

import anndata
peaks = anndata.AnnData(np.zeros((1, 2)), var=pd.DataFrame(index=["chr1-100-200", "chr3_5_9"]))
gb.Loci.from_anndata(peaks).to_records()
# -> [('chr1', 100, 200, '.'), ('chr3', 5, 9, '.')]
```

### Biopython

```python
from Bio import SeqIO, motifs as bm
recs = cre.to_seqrecords("genome.fa")                # SeqRecord per row, id = uid
idx = SeqIO.index("genome.fa", "fasta")              # a Biopython index is a FASTA source
cre.sequences(idx) == cre.sequences("genome.fa")     # -> True

lib = gb.load_motifs("motifs.jaspar")                # Library(2 motifs: M1, M2)
bio = lib.to_biopython()                             # list of Bio.motifs.Motif (counts)
gb.load_motifs(bio).names                            # -> ['M1', 'M2']
with open("motifs.jaspar") as fh:
    gb.load_motifs(list(bm.parse(fh, "jaspar")))     # Biopython's own parser, same Library

from genomeblocks import motifs as gm
M = gm.scan_motifs_matrix(cre, idx, bio, r=50, threshold=7.0, norm=False, verbose=False)
M.shape, M.index.name                                # -> ((7, 2), 'uid')
```

`load_motifs` also takes lightmotif motifs, a `{name: matrix}` dict (4 x W or
W x 4, counts or probabilities) or bare arrays (`motif1`, `motif2`, ...).

### xarray and signal cubes

`Loci.signal` returns a `(rows, tracks, bins)` numpy cube aligned to the
rows. Three converters give it coordinates:

```python
from genomeblocks import interop
S = cre.signal("atac.bw", n_bins=10, flank=250, verbose=False, progress=False)   # (7, 1, 10) float32
da = interop.cube_to_xarray(S, cre, ["ATAC"], flank=250)
da.dims                                # -> ('region', 'track', 'bin')
da.coords["bin"].values[:3]            # -> array([-225., -175., -125.])   bin centres in bp
da.mean("bin").to_pandas()             # per-region means, indexed by region name

interop.cube_to_pandas(S, cre, ["ATAC"])    # (track, bin) columns; one column per track if n_bins == 1
interop.cube_to_anndata(S, cre, ["ATAC"])   # tracks as obs, loci as var, profiles in varm['bins:<track>']
```

### bigWig and FASTA handles

```python
import pyBigWig, pybigtools, pysam, pyfaidx
from genomeblocks.signal import signal
h = pyBigWig.open("atac.bw")
S1 = signal(cre, [h], n_bins=10, flank=250, verbose=False, progress=False)
S2 = signal(cre, {"ATAC": pybigtools.open("atac.bw", "r")}, n_bins=10, flank=250,
            verbose=False, progress=False)
np.allclose(S, S1) and np.allclose(S, S2)       # -> True

cre.sequences(pysam.FastaFile("genome.fa"))     # pysam handle
cre.sequences(pyfaidx.Fasta("genome.fa"))       # pyfaidx handle
cre.sequences(gb.read_fasta("genome.fa"))       # a plain {chrom: str} dict
```

An open handle is used as it is; its engine is whatever opened it, so
`backend=` on such a call is an error rather than a silent switch (see
[Backends]({{ '/backends/' | relative_url }})).

### Graph libraries

```python
A.ep["w"][:] = [5.0, 3.0, 2.0, 1.0]   # weights (add_mcool fills them from a .cool)
A.normalize(verbose=False)            # adds ep.d and ep.n
g = A.to_networkx()                   # node i = row i; chrom/start/end on nodes, edge columns on edges
dict(g.nodes[0])                      # -> {'chrom': 'chr1', 'start': 900, 'end': 1100}
dict(g.edges[0, 2])                   # -> {'w': 5.0, 'd': 4000.0, 'n': 1.34...}
ig = A.to_igraph()                    # ig.es.attributes() -> ['w', 'd', 'n']
M = A.to_scipy("w")                   # symmetric CSR over the rows; to_scipy(None) is 1 per edge
gb.Architecture.from_scipy(cre, M).ep["w"].tolist()     # -> [5.0, 3.0, 2.0, 1.0]
adA = A.to_anndata()                  # obs = vertices, obsp['w'], obsp['d'], obsp['n']
```

---

## What copies and what does not
{: .sec-navy }

A `Loci` owns its columns. Every reader copies out of the input, so a table you
pass is never aliased: changing the frame afterwards does not change the Loci,
and `L["score"] = ...` does not change your frame. Views within genomeblocks
are a different matter and are listed in [Design → Tables]({{ '/design/tables/' | relative_url }}).

| Converter | Copies? | Why |
|---|---|---|
| `from_frame` and every `from_*` / `as_loci` | yes | starts and ends are cast to int64, chromosomes are encoded to codes |
| `to_arrow()` | **no** for start, end and numeric extra columns | pyarrow wraps the numpy buffers; chrom / strand are small dictionary arrays; object columns become Arrow strings (copied, `None` / `NaN` → null) |
| `to_polars()` | **no** for the numeric columns | built from `to_arrow()` |
| `__arrow_c_stream__`, `__dataframe__`, `__narwhals_dataframe__` | as `to_arrow()` | |
| `to_pandas()` | yes | pandas copies dict input into its blocks; chrom and strand are categoricals |
| `to_bioframe()` | yes | chromosome and strand names are decoded to strings |
| `to_pyranges()` | yes | pyranges splits the frame per chromosome |
| `to_bedtool()` | yes, to disk | pybedtools works on files |
| `to_cgranges()` | yes | an index is built, one `add` per row |
| `to_anndata()` | `var` yes; `X` as given (a numpy array or sparse matrix is not copied) | |
| `to_numpy()`, `to_records()` | yes | a structured array / a list of tuples |
| `cube_to_xarray` | **no** | the DataArray wraps the cube |
| `cube_to_pandas`, `cube_to_anndata` | yes / reduced | reshaped or aggregated over bins |
| `save()` / `load()` | to and from disk | parquet |

```python
np.shares_memory(L.starts, L.to_arrow().column("start").to_numpy())   # -> True
np.shares_memory(L.starts, L.to_polars()["start"].to_numpy())         # -> True
np.shares_memory(L.starts, L.to_pandas()["start"].to_numpy())         # -> False
np.shares_memory(df["start"].to_numpy(), gb.Loci.from_frame(df).starts)  # -> False
```

---

## Round-trip rules
{: .sec-navy }

`from_X(L.to_X())` gives back the same rows, in the same order, with the same
extra columns, for pandas, polars (eager and lazy), pyarrow, bioframe,
pyranges (sorted input), pybedtools, AnnData, parquet, BED and records. The
test suite checks each one. What each format keeps:

| Format | Coordinates and strand | Extra columns | Order |
|---|---|---|---|
| pandas, polars, pyarrow, parquet | all; `.` strands kept | all, dtypes kept (`None` / `NaN` in object columns survive as nulls) | kept |
| bioframe | all | all | kept |
| pyranges | all; `.` rows kept in `Strand` (pyranges calls the object unstranded) | all | grouped by chromosome on pyranges 0.x, see below |
| pybedtools | all | `name` and `score` only (BED6) | kept |
| AnnData | all; `var` carries chrom / start / end / strand | all, as `var` columns | kept |
| BED (`to_bed` / `Loci.make`) | all | `name` and `score` by choice | `Loci.make` sorts; `as_loci(path)` keeps |
| records / uids | all | none | kept |
| `Genes.to_gtf` → `Genes.make` | 0-based ↔ 1-based handled | gene_id, transcript_id, gene_name, gene_type, exon_number | — |
| `Pairs.to_bedpe` → `Pairs.make` | both anchors and strands | name, score | kept |
| `Architecture.edges_frame` → `from_frame` | vertex rows (or uids) | numeric edge columns | canonical (sorted into cis blocks) |

Chromosome names round-trip as given (`chr1` stays `chr1`, `1` stays `1`); a
`Genome` only maps names to codes.

### The pyranges 0.x row-order caveat

pyranges 0.x stores one DataFrame per chromosome (and strand), so a PyRanges
hands rows back grouped by chromosome, in the order the chromosomes were first
seen. Rows within a chromosome keep their order. Sort the Loci first when the
row order matters, or go through pandas / polars, which keep it:

```python
U = gb.as_loci([("chr2", 500, 600), ("chr1", 900, 1100), ("chr1", 100, 200)])
gb.Loci.from_pyranges(U.to_pyranges()).to_records()
# -> [('chr1', 900, 1100, '.'), ('chr1', 100, 200, '.'), ('chr2', 500, 600, '.')]
gb.Loci.from_pyranges(U.sort().to_pyranges()).equals(U.sort())      # -> True
gb.Loci.from_frame(U.to_pandas()).equals(U)                          # -> True
```

{: .tip }
> Everything computed from a `Loci` (a signal cube, a motif matrix, labels, an
> Architecture's vertex columns) is aligned to its rows. Convert the `Loci`
> last, or carry the uid (`to_pandas(uid=True)`) when you must reorder.

---

## Worked examples
{: .sec-purple }

### scATAC peaks: AnnData in, signal back into AnnData

```python
import anndata, numpy as np
adata = anndata.read_h5ad("atac.h5ad")             # var names like chr1:100-200
peaks = gb.as_loci(adata)                          # var order kept
S = peaks.signal({"H3K27ac": "k27ac.bw"}, n_bins=1, span=True, verbose=False, progress=False)
adata.var["H3K27ac"] = S[:, 0, 0]                  # row i of S is var i
adata.varm["H3K27ac_profile"] = peaks.signal("k27ac.bw", n_bins=50, flank=1000,
                                             verbose=False, progress=False)[:, 0, :]
```

### polars and bioframe in one pipeline

```python
import polars as pl, bioframe as bf
cre = gb.as_loci(pl.scan_csv("peaks.tsv", separator="\t").filter(pl.col("score") > 25))
near = cre.nearest(gb.Genes.make("genes.gtf").get_tss())       # (row in tss, distance)
cre["tss_distance"] = near[1]
report = bf.cluster(cre.to_bioframe(), min_dist=1000)           # bioframe's clustering
tidy = cre.to_polars().with_columns(pl.col("tss_distance").log1p().alias("log_d"))
```

### An Architecture as a networkx graph, and back

```python
A = (gb.Architecture.make(cre, "loops.bedpe", r=2500, verbose=False)
       .normalize(verbose=False))
g = A.to_networkx()                                # node i = cre row i, edges carry w, d, n
import networkx as nx
bc = nx.betweenness_centrality(g, weight="n")
A.vp["betweenness"] = np.array([bc[i] for i in range(len(cre))])   # back as a vertex column
A.vertices_frame()                                 # chrom, start, end, ..., betweenness per vertex
```

Related pages: [Loci guide]({{ '/guide/loci/' | relative_url }}),
[Backends]({{ '/backends/' | relative_url }}),
[API → interop]({{ '/api/interop/' | relative_url }}),
[Design → Tables]({{ '/design/tables/' | relative_url }}).
