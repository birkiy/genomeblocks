---
title: backends
parent: API Reference
layout: default
nav_order: 11
---

# `genomeblocks.backends`
{: .no_toc }

Which library does the work. genomeblocks keeps its data as columnar
tables; the heavy lifting on those tables runs through a *backend*, one per
family of operations. The registry (this module) resolves names, checks
installs and scopes choices; one submodule per family wraps the engines
behind a small interface that returns plain numpy, so a caller never sees
which engine ran. See [Backends]({{ '/backends/' | relative_url }}) for the
user-level overview and [Benchmarks]({{ '/benchmarks/' | relative_url }})
for the timings.
{: .fs-5 .fw-300 }

```python
import os, numpy as np
import genomeblocks as gb
from genomeblocks import backends                        # the module; gb.backends() lists the engines
from genomeblocks.backends import families, installed, resolve, use_backend, unsupported, AUTO
from genomeblocks.backends import intervals, bigwig, motifs, fasta, tables, graph
```

| Family | What it does | Default (automatic) | Others on request |
|---|---|---|---|
| `intervals` | overlap, nearest, merge, point lookups | `genomeblocks` (numpy) | `cgranges`, `ncls`, `bioframe`, `pyranges`, `bedtools` |
| `bigwig` | reading bigWig signal | `pybigtools` | `pybigwig`, `python` |
| `motifs` | scoring motif matrices along sequences | `moods`, else `lightmotif` | `lightmotif`, `biopython` |
| `fasta` | fetching sequence | `genomeblocks` (indexed reader) | `pysam`, `pyfaidx`, `memory`, `biopython` |
| `tables` | parsing text tables (BED, GTF, pairs) | `polars`, else `pandas` | `pandas` |
| `graph` | graph algorithms on an `Architecture` | `graph-tool`, else `scipy` | `igraph`, `networkx`, `scipy` |

Every default installs with `pip install genomeblocks` except a motif
engine: MOODS ships with the conda package (bioconda `moods`); on pip it is
the `genomeblocks[motifs]` extra (MOODS-python, compiled at install) or
`genomeblocks[lightmotif]` (prebuilt wheels), and a bare pip install has no
motif engine — `scan_*` then raise `ImportError` naming both extras. An
engine is picked automatically only where every candidate gives the same
answer (`AUTO` below). Anything else is a request: `backend=` on the call,
or `with gb.use_backend(...)` for a stretch of code. A requested backend
that is not installed raises `ImportError` with the install command — there
is no silent switch.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The registry
{: .sec-navy }

| Name | Signature | One line |
|---|---|---|
| `families` | `families() -> dict` | `{family: [backend, ...]}`, default first. |
| `installed` | `installed(family, backend) -> bool` | Whether that engine can be imported here (and, for `bedtools`, whether the binary is on `PATH`). Unknown names raise `ValueError`. |
| `resolve` | `resolve(family, backend=None) -> str` | The backend a call uses: `backend` if given, else the one chosen by an enclosing `use_backend`, else the automatic default. Unknown → `ValueError`; known but missing → `ImportError` with the install command. |
| `use_backend` | `use_backend(**choices)` (context manager) | `with use_backend(intervals="pyranges", graph="networkx"):` — names are checked on entry, so a missing backend fails before any work; `None` or `"auto"` restores the default. |
| `backends` | `backends() -> DataFrame` | Every family and backend: `installed`, `default`, `in use`, `install` hint. Also `gb.backends()` — the subpackage is callable. |
| `unsupported` | `unsupported(family, backend, op, supported) -> NotImplementedError` | The error a family module raises when an engine lacks an operation. |
| `AUTO` | `dict` | What `"auto"` may pick, in order: only engines with identical answers. |

Aliases are accepted anywhere a backend name is: `numpy` → `genomeblocks`,
`graph_tool` / `gt` → `graph-tool`, `nx` → `networkx`, `pybedtools` →
`bedtools`, `pyBigWig` → `pybigwig`, `Bio` → `biopython`, `MOODS` → `moods`,
`dict` → `memory`.

```python
gb.backends()[["family", "backend", "installed", "default", "in use"]].head(8)
# ->       family       backend  installed  default  in use
#    0  intervals  genomeblocks       True     True    True
#    1  intervals      cgranges      False    False   False
#    2  intervals          ncls       True    False   False
#    3  intervals      bioframe       True    False   False
#    4  intervals      pyranges       True    False   False
#    5  intervals      bedtools      False    False   False
#    6     bigwig    pybigtools       True     True    True
#    7     bigwig      pybigwig       True    False   False
families()
# -> {'intervals': ['genomeblocks', 'cgranges', 'ncls', 'bioframe', 'pyranges', 'bedtools'], 'bigwig': ['pybigtools', 'pybigwig', 'python'],
#     'motifs': ['moods', 'lightmotif', 'biopython'], 'fasta': ['genomeblocks', 'pysam', 'pyfaidx', 'memory', 'biopython'],
#     'tables': ['polars', 'pandas'], 'graph': ['graph-tool', 'scipy', 'igraph', 'networkx']}
AUTO
# -> {'intervals': ['genomeblocks'], 'bigwig': ['pybigtools', 'python'], 'motifs': ['moods', 'lightmotif'],
#     'fasta': ['genomeblocks'], 'tables': ['polars', 'pandas'], 'graph': ['graph-tool', 'scipy']}
installed("intervals", "bioframe"), installed("intervals", "cgranges"), installed("graph", "graph-tool"), installed("intervals", "numpy")
# -> (True, False, False, True)
```

```python
resolve("intervals"), resolve("graph"), resolve("tables"), resolve("bigwig"), resolve("graph", "nx")
# -> ('genomeblocks', 'scipy', 'polars', 'pybigtools', 'networkx')
resolve("intervals", "nope")
# -> ValueError: unknown intervals backend 'nope'; choose from: genomeblocks, cgranges, ncls, bioframe, pyranges, bedtools
resolve("nope")
# -> ValueError: unknown backend family 'nope'; families: intervals, bigwig, motifs, fasta, tables, graph
resolve("intervals", "cgranges")
# -> ImportError: the 'cgranges' intervals backend is not installed: conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)
resolve("graph", "graph-tool")
# -> ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool
with use_backend(intervals="bioframe", graph="networkx"):
    print(resolve("intervals"), resolve("graph"))
    with use_backend(intervals="auto"):
        print(resolve("intervals"))
# -> bioframe networkx
#    genomeblocks
resolve("intervals")
# -> 'genomeblocks'
with use_backend(intervals="cgranges"):                   # checked on entry
    pass
# -> ImportError: the 'cgranges' intervals backend is not installed: ...
unsupported("intervals", "ncls", "nearest", ("genomeblocks", "bioframe"))
# -> NotImplementedError("the 'ncls' intervals backend has no nearest; use one of: genomeblocks, bioframe")
```

The examples below use `cre = gb.Loci.make("peaks.bed")` (seven rows on
chr1 / chr2), `genome.fa`, `signal.bw`, a two-motif JASPAR file (`M1` =
ACGT, `M2` = GGAA) and `loops.bedpe`.

---

## `backends.intervals`
{: .sec-green }

Overlap, nearest, merge and point lookups. Every function takes columnar
`Loci` that share one `Genome` and returns row numbers, normalised to
genomeblocks' rules whatever the engine: intervals are half-open (`s1 < e2
and s2 < e1`); `overlap_pairs` is sorted by (query row, reference row);
`nearest` gives 0 to an overlap, else the gap in bases (book-ended = 0), -1
when the chromosome has nothing; `merge` fuses overlapping and book-ended
intervals. Zero-length intervals are points between two bases — `[p, p)`
overlaps `[s, e)` exactly when `s < p < e` — and every engine's pairs are
filtered with that rule. The chromosome *code* is handed to the engines as
the chromosome label on both sides, so no decoding happens.

| Function | Signature | Returns |
|---|---|---|
| `overlap_pairs` | `overlap_pairs(q, r, *, backend=None)` | `(query rows, reference rows)` of every overlapping pair. |
| `overlap_any` | `overlap_any(q, r, *, backend=None)` | Boolean mask over the query rows. |
| `nearest` | `nearest(q, r, *, backend=None)` | `(reference row, distance)` per query row. Ties among overlapping intervals follow the genomeblocks rule (lowest start, then lowest row) on every engine. |
| `merge` | `merge(L, *, backend=None)` | `(codes, starts, ends, first)` of the merged blocks; `first` is the row of the first interval (lowest start) of each block. |
| `point_rows` | `point_rows(L, code, start, end, *, backend=None)` | Rows of `L` overlapping one window (sorted); the lookup index is cached on `L` per backend. |

| | genomeblocks | cgranges | ncls | bioframe | pyranges | bedtools |
|---|---|---|---|---|---|---|
| overlap | yes | yes | yes | yes | yes | yes |
| nearest | yes | – | – | yes | yes | yes |
| merge | yes | – | – | yes | yes | yes |
| point lookup | yes | yes | yes | yes | yes | yes |

These are what `Loci.overlap_pairs`, `overlap_any`, `nearest`, `merge`,
`overlap_rows` and the set operators call; at the `Loci` level the other
side is first brought onto the same `Genome` (`Loci._check`).

```python
cre = gb.Loci.make("peaks.bed")
other = cre._check(gb.as_loci([("chr1", 1000, 1050), ("chr2", 0, 10_000)]))     # same Genome
intervals.overlap_pairs(cre, other), intervals.overlap_pairs(cre, other, backend="ncls")
# -> ((array([0, 5, 6]), array([0, 1, 1])), (array([0, 5, 6]), array([0, 1, 1])))
intervals.overlap_any(cre, other)
# -> array([ True, False, False, False, False,  True,  True])
intervals.nearest(cre, other)
# -> (array([0, 0, 0, 0, 0, 1, 1]), array([   0,  850, 3850, 8900, 9850,    0,    0]))
intervals.nearest(cre, other, backend="pyranges")
# -> (array([0, 0, 0, 0, 0, 1, 1]), array([   0,  850, 3850, 8900, 9850,    0,    0]))
intervals.nearest(cre, other, backend="ncls")
# -> NotImplementedError: the 'ncls' intervals backend has no nearest; use one of: genomeblocks, bioframe, pyranges, bedtools
m = gb.as_loci([("chr1", 10, 20), ("chr1", 20, 30), ("chr1", 40, 50), ("chr2", 0, 5)])
intervals.merge(m)
# -> (array([0, 0, 1], dtype=int32), array([10, 40,  0]), array([30, 50,  5]), array([0, 2, 3]))
intervals.merge(m, backend="bioframe")
# -> (array([0, 0, 1]), array([10, 40,  0]), array([30, 50,  5]), array([0, 2, 3]))
intervals.point_rows(cre, cre.genome["chr1"], 1000, 2000), intervals.point_rows(cre, cre.genome["chr1"], 1000, 2000, backend="ncls")
# -> (array([0, 1]), array([0, 1]))
```

{: .note }
`bedtools` widens empty intervals on its side, so it may miss pairs that
need the widened base; every other engine agrees exactly after the
half-open filter. The test suite checks each engine against the default.

---

## `backends.bigwig`
{: .sec-green }

One handle interface over pybigtools (Rust, the default), pyBigWig (C,
libBigWig) and a pure-Python reader (numpy + mmap, no compiled dependency).
Bin `b` of a window of `n` bases is `[floor(n*b/n_bins), floor(n*(b+1)/n_bins))`
for every engine, bins outside the chromosome are `missing`, and windows
that leave the chromosome (or name one the file lacks) are clipped and
padded the same way, so the three give the same numbers.

| Name | Signature | One line |
|---|---|---|
| `open_bigwig` | `open_bigwig(src, *, backend=None) -> handle` | A uniform handle for a bigWig path or an already-open pyBigWig / pybigtools / `BigWigReader` handle (used as is, wrapped in its own backend). Any object with a `stats_array` method is accepted unchanged. |
| `handle_backend` | `handle_backend(src) -> str | None` | The backend an open handle belongs to (`None` for a path). |
| `PyBigToolsHandle`, `PyBigWigHandle`, `PythonHandle` | classes | The three handles; `.backend` names the engine. |

The handle interface:

| Method | Returns |
|---|---|
| `h.chroms()` | `{chrom: size}` (cached). |
| `h.stats_array(chrom, start, end, *, n_bins=1, stat="mean", exact=True, missing=0.0)` | `float64` array of `n_bins` values; `stat` is `mean`, `min`, `max`, `sum`, `std` (population) or `coverage` (fraction of bases with data). |
| `h.values(chrom, start, end)` | Per-base `float64` values, NaN where there is no data (also outside the chromosome). |
| `h.close()` | Closes the file when the handle opened it. |

`Loci.signal` / `signal.signal` open every track through `open_bigwig`
and pass `backend=`; passing `backend=` together with someone else's open
handle of another engine raises `ValueError`.

```python
h = bigwig.open_bigwig("signal.bw")
type(h).__name__, h.backend, h.chroms()
# -> ('PyBigToolsHandle', 'pybigtools', {'chr1': 20000})
h.stats_array("chr1", 0, 1000, n_bins=5).round(3)
# -> array([5.08 , 1.259, 8.767, 7.013, 7.654])
h.stats_array("chr1", 0, 1000, n_bins=5, stat="coverage")        # 50 bp of data per 100 bp
# -> array([0.5, 0.5, 0.5, 0.5, 0.5])
h.values("chr1", 45, 55)                                         # data ends at base 50
# -> array([6.73265505, 6.73265505, 6.73265505, 6.73265505, 6.73265505, nan, nan, nan, nan, nan])
h.stats_array("chr9", 0, 100, n_bins=2, missing=np.nan)          # unknown chromosome
# -> array([nan, nan])
h.close()
for b in ("pybigwig", "python"):
    hb = bigwig.open_bigwig("signal.bw", backend=b)
    assert np.allclose(hb.stats_array("chr1", 0, 1000, n_bins=5), bigwig.open_bigwig("signal.bw").stats_array("chr1", 0, 1000, n_bins=5))
    hb.close()
import pyBigWig
raw = pyBigWig.open("signal.bw")
bigwig.handle_backend(raw), bigwig.handle_backend("signal.bw"), type(bigwig.open_bigwig(raw)).__name__
# -> ('pybigwig', None, 'PyBigWigHandle')
bigwig.open_bigwig(raw, backend="pybigtools")
# -> ValueError: the track is an open pybigwig handle, so backend='pybigtools' cannot apply: pass the file path instead, or drop backend=
```

---

## `backends.motifs`
{: .sec-green }

One motif library, three scanning engines. A [`Library`](#library) holds
motifs as `(W x 4)` count matrices (columns A C G T) whatever they were read
from — JASPAR (raw or bracketed), TRANSFAC, UniPROBE or MEME files, parsed
here in pure Python (`.gz` too); Biopython `Bio.motifs` objects; lightmotif
motifs; plain arrays. Each motif is scored with one log-odds matrix,
`log2((count + 0.1) / (column total + 0.4) / 0.25)`, computed by
[`logodds_matrix`](#logodds_matrix-and-threshold_from_pvalue) bit for bit
the way lightmotif computes it (float32, its summation order, the platform
`log2f`), and every engine scans that same matrix, so **MOODS** (C++, the
default), **lightmotif** (SIMD, used automatically only when MOODS is
absent) and **Biopython** (numpy) report the same hits. A hit is a window
position whose score is at least the threshold; windows containing N never
hit. `threshold_from_pvalue` turns a p-value into a per-motif score cutoff,
the same whichever engine scans.

| Name | Signature | One line |
|---|---|---|
| `Library` | class | Motifs as count matrices + names; a table of `name, description, width, consensus`. |
| `load_motifs` | `load_motifs(src, format="jaspar") -> Library` | A `Library` from a motif file or motif objects (also `gb.load_motifs`). |
| `logodds_matrix` | `logodds_matrix(counts, pseudocount=0.1) -> ndarray` | `(W x 4)` log2-odds against a uniform background, float64 holding float32-exact values. |
| `threshold_from_pvalue` | `threshold_from_pvalue(lo_matrix, pvalue, *, resolution=1e-3) -> float` | The smallest score with `P(score >= s) <= pvalue` under a uniform background. |
| `write_meme` | `write_meme(pfms, path, *, alphabet="ACGT", bg=(0.25, 0.25, 0.25, 0.25))` | `{name: (4 x W) probability matrix}` as MEME minimal format. |
| `Block` | class | Windows laid end to end once, scanned per motif by one engine. |
| `FORMATS` | `('jaspar', 'jaspar16', 'transfac', 'uniprobe', 'meme')` | Accepted `format=` values. |
| `PSEUDOCOUNT` | `0.1` | The default pseudocount. |

### `Library`

```python
Library(names, counts, descriptions=None, *, pseudocount=0.1)
```

| Name | Returns | One line |
|---|---|---|
| `len(lib)`, `lib.shape`, `lib.columns` | | `['name', 'description', 'width', 'consensus']`. |
| `lib.names`, `lib.descriptions`, `lib.counts` | lists | Names, header descriptions (TF symbols in JASPAR), `(W x 4)` count matrices. |
| `lib.widths` | `int64` array | Motif widths. |
| `lib[i]`, `lib['CTCF']` | `ndarray` | The count matrix by position or name (a unique substring match works too). |
| `for name in lib` | | Iterates the names. |
| `consensus(i)` | `str` | Most frequent base per position. |
| `logodds(i)` | `(W x 4)` | Log2-odds matrix (computed once, cached). |
| `pfm(i, pseudo=0.01)` | `(4 x W)` | Probability matrix, rows A C G T — the archetype functions' format. |
| `select(motifs=None, match="substring")` | `Library` | The motifs named; `substring` is case-insensitive and searches descriptions too, `exact` compares names. |
| `indices(motifs=None, match="substring")` | `list[int]` | Their positions. |
| `take(idx)` | `Library` | A subset by position. |
| `to_pandas()`, `to_arrow()`, `head()`, `tail()`, `describe()`, `_repr_html_()` | | Like every table; the Arrow / dataframe / narwhals protocols read `to_arrow()`. |
| `to_biopython()` | `list[Bio.motifs.Motif]` | Counts, `matrix_id` = name, `name` = description. |
| `to_moods()` | `list` | One `4 x W` log-odds list of lists per motif (A C G T rows), what `MOODS.scan.Scanner.set_motifs` takes. |
| `to_meme(path, pseudo=0.01)`, `to_jaspar(path)` | | Write the library (MEME probabilities / bracketed JASPAR counts). |

```python
lib = gb.load_motifs("motifs.jaspar")
lib, lib.columns, lib.shape
# -> (Library(2 motifs: M1, M2), ['name', 'description', 'width', 'consensus'], (2, 4))
lib.to_pandas()
# ->   name description  width consensus
#    0   M1         TFA      4      ACGT
#    1   M2         TFB      4      GGAA
lib["M1"]
# -> array([[10.,  0.,  0.,  0.],
#           [ 0., 10.,  0.,  0.],
#           [ 0.,  0., 10.,  0.],
#           [ 0.,  0.,  0., 10.]])
lib[0].shape, lib.widths, lib.consensus(1), lib.pfm(0).shape
# -> ((4, 4), array([4, 4]), 'GGAA', (4, 4))
lib.logodds(0).round(2)
# -> array([[ 1.96, -4.7 , -4.7 , -4.7 ],
#           [-4.7 ,  1.96, -4.7 , -4.7 ],
#           [-4.7 , -4.7 ,  1.96, -4.7 ],
#           [-4.7 , -4.7 , -4.7 ,  1.96]])
lib.select("TF"), lib.select("M1", match="exact"), lib.indices("tfb"), lib.take([1])
# -> (Library(2 motifs: M1, M2), Library(1 motifs: M1), [1], Library(1 motifs: M2))
list(lib), lib.to_moods()[0][0][:2], lib.to_biopython()[0].matrix_id
# -> (['M1', 'M2'], [1.9577716588974, -4.700439929962158], 'M1')
lib.to_meme("lib.meme"); lib.to_jaspar("lib.jaspar")
gb.load_motifs("lib.meme", format="meme"), gb.load_motifs("lib.jaspar", format="jaspar16")
# -> (Library(2 motifs: M1, M2), Library(2 motifs: M1, M2))
```

### `load_motifs`

```python
load_motifs(src, format="jaspar") -> Library
```

`src` may be a path (`format` = `'jaspar'` for raw four-row counts,
`'jaspar16'` for the bracketed layout, `'transfac'`, `'uniprobe'` or
`'meme'`), a `Library` (returned as is), one or many Biopython
`Bio.motifs.Motif`, lightmotif motifs, or a `{name: matrix}` dict (or a bare
list) of counts / probabilities, `4 x W` or `W x 4`, rows or columns in A C
G T order. Probability matrices (rows summing to 1) are scaled to 100 sites;
MEME motifs by their `nsites` (100 when absent). Names follow lightmotif's
conventions: the JASPAR id and the rest of the header line as description;
TRANSFAC `NA`, else `ID`, else `AC`.

```python
gb.load_motifs({"X": np.array([[10, 0, 0, 0], [0, 10, 0, 0]])}), gb.load_motifs([np.eye(4)])
# -> (Library(1 motifs: X), Library(1 motifs: motif1))
gb.load_motifs(lib.to_biopython()).names, gb.load_motifs(lib) is lib
# -> (['M1', 'M2'], True)
gb.load_motifs("motifs.jaspar", format="nope")
# -> ValueError: unknown motif format 'nope'; use one of jaspar, jaspar16, transfac, uniprobe, meme
motifs.write_meme({"CTCF": lib.pfm(0)}, "ctcf.meme")
open("ctcf.meme").read().count("MOTIF")
# -> 1
```

### `logodds_matrix` and `threshold_from_pvalue`

```python
logodds_matrix(counts, pseudocount=0.1) -> ndarray            # (W x 4), A C G T
threshold_from_pvalue(lo_matrix, pvalue, *, resolution=1e-3) -> float
```

`threshold_from_pvalue` computes the exact score distribution of one
log-odds matrix under a uniform background on a `resolution`-bit grid (a
coarse pass bounds the answer from below so only the upper tail is built)
and returns the smallest score whose tail probability is at most `pvalue`.
It agrees with MOODS' `threshold_from_p` to a few thousandths of a bit, and
is what every `scan_*(pvalue=...)` call uses, so the cutoffs do not depend
on the engine. A motif too short to reach `pvalue` gets its best possible
score (only perfect matches count).

```python
motifs.logodds_matrix(lib["M1"]).round(3)
# -> array([[ 1.958, -4.7  , -4.7  , -4.7  ],
#           [-4.7  ,  1.958, -4.7  , -4.7  ],
#           [-4.7  , -4.7  ,  1.958, -4.7  ],
#           [-4.7  , -4.7  , -4.7  ,  1.958]])
for p in (0.1, 0.01, 1e-3):
    print(p, round(motifs.threshold_from_pvalue(lib.logodds(0), p), 3))
# -> 0.1 -5.483
#    0.01 1.175
#    0.001 7.831                  # 4 bp: only the perfect match (score 7.83) is that rare
import MOODS.tools
round(MOODS.tools.threshold_from_p(lib.logodds(0).T.tolist(), [0.25] * 4, 0.01), 3)
# -> 1.174
```

### `Block`

```python
Block(seqs, backend=None)
```

The scanning unit: the windows are upper-cased, those that are empty or
hold letters other than A C G T N are dropped (`rows` maps block windows
back to input rows), and the rest are laid end to end in `text` with
`offsets`. Each motif is scanned once over the concatenation and the hits
split back per window — a hit straddling two windows is dropped — so
per-window results equal scanning each window alone. `backend` is resolved
on construction.

| Method | Returns |
|---|---|
| `hits(mat, threshold, both=False)` | `(input row, position in the window, strand)` of every hit of a `(W x 4)` log-odds matrix; strand 0 `+`, 1 `-` (the reverse complement of the matrix is scanned when `both`). |
| `counts(mat, threshold, both=False)` | Hits per input window (`int64`, length = number of input windows). |
| `moods_counts(mats, thresholds, both=False)` | MOODS only: many matrices in one pass, `(windows x motifs)` counts. |

```python
seqs = cre.sequences("genome.fa", r=50)                   # seven 100-bp windows
blk = motifs.Block(seqs)
blk.backend, blk.n, len(blk.text), blk.offsets
# -> ('moods', 7, 700, array([  0, 100, 200, 300, 400, 500, 600, 700]))
r, p, s = blk.hits(lib.logodds(0), 7.0, both=True)        # ACGT x 3 at the centre of window 0 ...
r, p, s
# -> (array([0, 0, 0, 5, 0, 0, 0, 5]), array([50, 54, 58, 44, 50, 54, 58, 44]), array([0, 0, 0, 0, 1, 1, 1, 1], dtype=int8))
blk.counts(lib.logodds(0), 7.0, both=True)                # ... ACGT is its own reverse complement
# -> array([6, 0, 0, 0, 0, 2, 0])
motifs.Block(seqs, backend="lightmotif").counts(lib.logodds(0), 7.0, both=True)
# -> array([6, 0, 0, 0, 0, 2, 0])
motifs.Block(seqs, backend="biopython").counts(lib.logodds(0), 7.0, both=True)
# -> array([6, 0, 0, 0, 0, 2, 0])
blk.moods_counts([lib.logodds(0), lib.logodds(1)], [7.0, 7.0], both=True).T
# -> array([[6, 0, 0, 0, 0, 2, 0],
#           [1, 0, 3, 0, 1, 0, 2]])
motifs.Block(["ACGTNACGT", ""]).rows                       # N is fine, an empty window is dropped
# -> array([0])
```

---

## `backends.fasta`
{: .sec-green }

Where sequence comes from. `genomeblocks` (the default) is a pure-Python
reader on the samtools `.fai` index — built next to the file when missing —
through a memory map, so only the requested bases are read; `pysam` (htslib),
`pyfaidx`, `memory` (the whole file read into a `{chrom: str}` dict) and
`biopython` (`Bio.SeqIO.index`) return the same bases. A gzip (not bgzip)
FASTA cannot be indexed and is read into memory.

| Name | Signature | One line |
|---|---|---|
| `open_fasta` | `open_fasta(src, *, backend=None) -> source` | A sequence source for a FASTA path, a `{chrom: str}` dict, a Biopython `SeqIO.index` / `to_dict` mapping, a `pyfaidx.Fasta` or a `pysam.FastaFile` (open objects and dicts are used as they are). Also `gb.read_fasta`'s dict. |
| `read_fasta` | `read_fasta(path) -> dict[str, str]` | Every record of a FASTA (plain or `.gz`) as `{name: sequence}`; `gb.read_fasta`. |
| `IndexedSource`, `PysamSource`, `PyfaidxSource`, `MemorySource`, `BiopythonSource` | classes | The five sources; `.backend` names the engine. |

The source interface:

| Method | Returns |
|---|---|
| `src.sizes()` | `{chrom: length}`. |
| `src.fetch(chrom, start, end)` | `str`, as stored (case kept); clipped to the record. |
| `src.fetch_many(chroms, starts, ends)` | `list[str]`. |
| `src.close()` | Closes what the source opened. |

`Loci.sequences`, `Locus.sequence`, the motif scanners (`fasta_backend=`)
and `Genome.from_fasta` all read through `open_fasta`.

```python
src = fasta.open_fasta("genome.fa")
type(src).__name__, src.backend, src.sizes()
# -> ('IndexedSource', 'genomeblocks', {'chr1': 20000, 'chr2': 8000})
src.fetch("chr1", 1000, 1012), src.fetch_many(["chr1", "chr2"], [1000, 0], [1004, 4])
# -> ('ACGTACGTACGT', ['ACGT', 'ATGT'])
src.close()
os.path.exists("genome.fa.fai")                           # the index was written next to the file
# -> True
for b in ("pysam", "pyfaidx", "memory", "biopython"):
    s = fasta.open_fasta("genome.fa", backend=b)
    print(b, type(s).__name__, s.fetch("chr1", 1000, 1012))
    s.close()
# -> pysam PysamSource ACGTACGTACGT
#    pyfaidx PyfaidxSource ACGTACGTACGT
#    memory MemorySource ACGTACGTACGT
#    biopython BiopythonSource ACGTACGTACGT
type(fasta.open_fasta({"chr1": "ACGT"})).__name__, fasta.open_fasta({"chr1": "ACGT"}).fetch("chr1", 1, 3)
# -> ('MemorySource', 'CG')
import pysam
type(fasta.open_fasta(pysam.FastaFile("genome.fa"))).__name__
# -> 'PysamSource'
d = fasta.read_fasta("genome.fa")
list(d), len(d["chr1"])
# -> (['chr1', 'chr2'], 20000)
type(fasta.open_fasta("genome.fa.gz")).__name__          # gzip: read into memory
# -> 'MemorySource'
fasta.open_fasta("genome.fa.gz", backend="genomeblocks")
# -> ValueError: a gzip FASTA cannot be indexed
gb.Genome.from_fasta("genome.fa")
# -> Genome(2 chroms: chr1, chr2)
```

{: .warning }
The indexed reader needs fixed-width records (every line of a record the
same length except the last); a FASTA with uneven lines raises a
`ValueError` suggesting `seqkit seq -w 60` or `backend='memory'` /
`'pysam'` / `'pyfaidx'`.

---

## `backends.tables`
{: .sec-green }

Parsing tab-separated text (BED, BEDPE, narrowPeak, pairs) with polars (the
default when installed) or pandas. Both return the same thing — a dict of
numpy columns, positions as `int64`, floats as `float64`, everything else
as object strings — so callers never touch the frame library. Leading `#` /
`track` / `browser` lines are skipped and `.gz` files are read directly.
`Loci.make`, `Pairs.make` and `Genes.make` (GTF / GFF3 parsing, in `genes.py`)
take `backend=` for this family.

| Function | Signature | Returns |
|---|---|---|
| `sniff` | `sniff(path, max_lines=50)` | `(header lines to skip, number of columns of the first data line)`. |
| `read_columns` | `read_columns(path, columns, *, ints=(), floats=(), backend=None)` | `{column index: array}` for the 0-based `columns` of a header-less TSV; `ints` become int64, `floats` float64, the rest strings. A `#` inside a data field is data. |

```python
tables.sniff("peaks.bed"), tables.sniff("loops.bedpe")    # two header lines, six columns
# -> ((2, 6), (0, 10))
d = tables.read_columns("peaks.bed", [0, 1, 2, 4], ints=[1, 2], floats=[4])
{k: v[:2] for k, v in d.items()}
# -> {0: array(['chr1', 'chr1'], dtype=object), 1: array([ 900, 1900]), 2: array([1100, 2100]), 4: array([10., 20.])}
tables.read_columns("peaks.bed", [0, 1], ints=[1], backend="pandas")[1]
# -> array([  900,  1900,  4900,  9950, 10900,   500,  5000])
open("h.bed", "w").write("chrom\tstart\tend\nchr1\t1\t5\n")
tables.read_columns("h.bed", [0, 1, 2], ints=[1, 2])
# -> ValueError: h.bed: cannot read as a tab-separated BED-like table (column 2): ... Header lines must start with '#', 'track' or 'browser';
#    a file with a column-name header is read with Loci.from_frame(pandas.read_csv(path, sep='\t')) or as_loci(path) for .tsv / .csv
```

---

## `backends.graph`
{: .sec-green }

Graph algorithms on an [`Architecture`]({{ '/api/architecture/' | relative_url }}):
graph-tool (the default when installed), scipy (the pip default), igraph
and networkx. The Architecture is two tables (vertices = `Loci` rows, edges =
`src` / `tgt` + edge columns); a backend turns them into its own graph
object when an algorithm needs one. Every function returns plain numpy
aligned to the `Loci` rows, so results do not depend on the engine.

| Function | Signature | Returns |
|---|---|---|
| `to_scipy` | `to_scipy(A, weight="w", *, n=None)` | Symmetric CSR adjacency `(n x n)` over the `Loci` rows; `weight` an edge column name, an array, or `None` for 1 per edge; parallel edges are summed. |
| `native` | `native(A, backend=None, *, vprops=(), loci_cols=False)` | The engine's own graph (vertex `i` = row `i`): a `graph_tool.Graph`, `igraph.Graph`, `networkx.Graph` or the scipy matrix. Edge columns become edge attributes; `vprops` names vertex columns to copy; `loci_cols=True` adds `chrom` / `start` / `end` per vertex. networkx holds one edge per pair, so parallel edges are merged with their numeric columns summed. |
| `components` | `components(A, backend=None) -> int64 array` | Connected component per row, numbered by first row (so every engine agrees); -1 for rows without links. |
| `pagerank` | `pagerank(A, weight=None, *, damping=0.85, backend=None, tol=1e-10) -> float64 array` | PageRank per row (unlinked rows included, as every engine does). |
| `layout` | `layout(A, rows, *, kind="spring", backend=None, seed=0) -> (len(rows) x 2)` | 2-D positions for the subgraph induced by `rows`: `'spring'` (the engine's force-directed layout; spectral for scipy), `'circular'`, or `'genomic'` (x = genomic centre, y = 0). |
| `layout_edges` | `layout_edges(m, s, t, *, kind="spring", backend=None, seed=0, x=None)` | The same for `m` vertices joined by edges `(s[k], t[k])`; `x` feeds the genomic layout. |

`Architecture.graph(backend=)`, `A.components()`, `A.pagerank()` and
`architecture_draw.draw` call these.

```python
A = gb.Architecture.make(cre, gb.Pairs.make("loops.bedpe"), r=100, verbose=False)
A, A.src.tolist(), A.tgt.tolist()
# -> (Architecture(name='Skeleton', loci=6, links=4 [3 cis · 1 trans], edge_props=[w], vertex_props=[]), [0, 1, 5, 0], [2, 4, 6, 5])
M = graph.to_scipy(A)
type(M).__name__, M.shape, M.nnz
# -> ('csr_matrix', (7, 7), 8)
graph.to_scipy(A, weight=None).toarray().astype(int)[:3]
# -> array([[0, 0, 1, 0, 0, 1, 0],
#           [0, 0, 0, 0, 1, 0, 0],
#           [1, 0, 0, 0, 0, 0, 0]])
G = graph.native(A, "networkx")
print(G); list(G.edges(data=True))[:2]
# -> Graph with 7 nodes and 4 edges
#    [(0, 2, {'w': 0.0}), (0, 5, {'w': 0.0})]
graph.native(A, "igraph").summary()
# -> 'IGRAPH U--- 7 4 -- \n+ attr: w (e)'
graph.native(A, "networkx", loci_cols=True).nodes[0]
# -> {'chrom': 'chr1', 'start': 900, 'end': 1100}
graph.components(A).tolist(), graph.components(A, backend="networkx").tolist()
# -> ([0, 1, 0, -1, 1, 0, 0], [0, 1, 0, -1, 1, 0, 0])
graph.pagerank(A).round(3).tolist()
# -> [0.211, 0.163, 0.114, 0.024, 0.163, 0.211, 0.114]
np.allclose(graph.pagerank(A), graph.pagerank(A, backend="igraph"), atol=1e-6)
# -> True
graph.layout(A, [0, 2, 5], kind="genomic").tolist()
# -> [[1000.0, 0.0], [5000.0, 0.0], [550.0, 0.0]]
graph.layout(A, [0, 2, 5], kind="circular").round(2).tolist()
# -> [[1.0, 0.0], [-0.5, 0.87], [-0.5, -0.87]]
graph.layout(A, np.arange(len(cre))).shape, graph.layout_edges(3, [0, 1], [1, 2], backend="networkx").shape
# -> ((7, 2), (3, 2))
graph.to_scipy(A, weight="nope")
# -> ValueError: edge column 'nope' not found; have: w
graph.pagerank(A, backend="graph-tool")                    # on a pip install
# -> ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool
```

{: .note }
`make` leaves the edge weight `w` at 0 until `add_mcool` fills it, which is
why the weighted PageRank above would be uniform; `pagerank(A, weight="w")`
is the call to make after normalisation.
