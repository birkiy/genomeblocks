---
title: Backends
layout: default
nav_order: 4.6
permalink: /backends/
---

# Backends
{: .no_toc }

The heavy work on genomeblocks tables runs through a *backend*: one library
per family of operations. The defaults install with `pip install genomeblocks`
and are chosen for correctness and installability. Any other engine is a
request, per call or per block of code, and a request that cannot be met
raises with the install command. Whatever engine runs, the answer is the same.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The six families
{: .sec-navy }

| Family | What it does | Default | On request |
|---|---|---|---|
| `intervals` | overlap, nearest, merge, point lookups | `genomeblocks` (numpy kernels) | `cgranges`, `ncls`, `bioframe`, `pyranges`, `bedtools` |
| `bigwig` | reading bigWig signal | `pybigtools` | `pybigwig`, `python` |
| `motifs` | scoring motif matrices along sequences | `lightmotif` | `moods`, `biopython` |
| `fasta` | fetching sequence | `genomeblocks` (indexed `.fai` reader) | `pysam`, `pyfaidx`, `memory`, `biopython` |
| `tables` | parsing text tables (BED, BEDPE, GTF, pairs) | `polars`, else `pandas` | `pandas` |
| `graph` | graph algorithms on an `Architecture` | `graph-tool`, else `scipy` | `igraph`, `networkx`, `scipy` |

```python
import genomeblocks as gb
gb.backends()                     # the module is callable: a table of every family and backend
```

```
   family      backend  installed  default  in use
intervals genomeblocks       True     True    True
intervals     cgranges      False    False   False
intervals         ncls       True    False   False
intervals     bioframe       True    False   False
intervals     pyranges       True    False   False
intervals     bedtools      False    False   False
   bigwig   pybigtools       True     True    True
   bigwig     pybigwig       True    False   False
   bigwig       python       True    False   False
   motifs   lightmotif       True     True    True
   motifs        moods       True    False   False
   motifs    biopython       True    False   False
    fasta genomeblocks       True     True    True
    fasta        pysam       True    False   False
    fasta      pyfaidx       True    False   False
    fasta       memory       True    False   False
    fasta    biopython       True    False   False
   tables       polars       True     True    True
   tables       pandas       True    False   False
    graph   graph-tool      False    False   False
    graph        scipy       True     True    True
    graph       igraph       True    False   False
    graph     networkx       True    False   False
```

The table has a sixth column, `install`, with the command for each missing
backend. `gb.backends.families()` returns the same names as a dict,
`gb.backends.installed(family, name)` says whether one can run here, and
`gb.backends.resolve(family, name=None)` returns the backend a call would use.

### Why these defaults

Every default is pure pip, so `pip install genomeblocks` gives a working
package with no conda step and no compiler:

- **intervals → genomeblocks.** Overlap, merge and nearest are a handful of
  numpy calls on one genome-wide axis (`code · 2⁴⁰ + position`), so they need
  nothing beyond numpy and are the fastest option for whole-set work. Point
  lookups use a cached index built on first use. cgranges (C) and bedtools are
  conda-only; ncls, bioframe and pyranges are pip-installable alternatives for
  people who already have them.
- **bigwig → pybigtools.** A Rust reader that does I/O, decompression and
  binning in one call per region; pyBigWig (libBigWig, C) gives identical
  numbers, and the pure-Python reader (`python`, numpy + mmap) exists so the
  package works where neither wheel does.
- **motifs → lightmotif.** SIMD scanning of a striped sequence; MOODS and
  Biopython score the same log-odds matrices.
- **fasta → genomeblocks.** A reader on the samtools `.fai` index (built next
  to the file when missing), so only the requested bases are read. pysam,
  pyfaidx and Biopython give the same bases; `memory` reads the whole file
  into a dict (and is what a gzip FASTA, which cannot be indexed, gets).
- **tables → polars, else pandas.** Both parse to the same numpy columns.
  polars is an optional extra (`pip install "genomeblocks[fast]"`); the
  fallback to pandas is automatic and loses nothing but speed.
- **graph → graph-tool, else scipy.** graph-tool is the owner's choice when it
  is installed (conda only). scipy gives identical components and PageRank, so
  the pip install never lacks the `Architecture`.

---

## The automatic choice: only identical answers
{: .sec-navy }

`backend=None` (the default on every call) means *auto*. Auto never guesses
between engines that could disagree: it walks a short list of candidates known
to give the same result and takes the first one installed.

| Family | Auto tries, in order |
|---|---|
| `intervals` | `genomeblocks` |
| `bigwig` | `pybigtools`, `python` |
| `motifs` | `lightmotif` |
| `fasta` | `genomeblocks` |
| `tables` | `polars`, `pandas` |
| `graph` | `graph-tool`, `scipy` |

So auto is `genomeblocks` for intervals even when pyranges is installed, and
`scipy` for graphs unless graph-tool is present. The list is
`genomeblocks.backends.AUTO`.

```python
gb.backends.resolve("intervals")      # -> 'genomeblocks'
gb.backends.resolve("graph")          # -> 'scipy'        (graph-tool not installed here)
gb.backends.resolve("tables")         # -> 'polars'
```

---

## Choosing a backend
{: .sec-green }

There is no global switch and no hidden state. Two scopes exist:

**One call.** Every function that runs through a backend takes `backend=`:

```python
a = gb.as_loci("peaks.bed")
b = gb.as_loci([("chr1", 1000, 1050), ("chr2", 0, 10_000), ("chr1", 2000, 2000)])
a.overlap_pairs(b)                        # -> (array([0, 1, 5, 6]), array([0, 2, 1, 1]))
a.overlap_pairs(b, backend="bioframe")    # the same pairs
cre = gb.Loci.make("peaks.bed", backend="pandas")               # table parser
S = cre.signal("atac.bw", backend="pybigwig", verbose=False)    # bigWig reader
A.components(backend="networkx")                                # graph algorithm
```

The motif functions take two: `backend=` for the scanner and
`fasta_backend=` for the sequence source.

**A block of code.** `gb.use_backend(**families)` sets backends inside a
`with` statement and restores the previous choice on exit, even on an error.
Blocks nest; `None` or `"auto"` restores the automatic choice inside a block.

```python
with gb.use_backend(intervals="pyranges", bigwig="pybigwig"):
    gb.backends.resolve("intervals")      # -> 'pyranges'
    a.intersect(b)                        # runs on pyranges
    with gb.use_backend(intervals="auto"):
        gb.backends.resolve("intervals")  # -> 'genomeblocks'
gb.backends.resolve("intervals")          # -> 'genomeblocks'   (restored)
```

Names are checked when the block is entered, so a missing backend fails
before any work starts. An explicit `backend=` on a call inside a block wins
over the block.

Aliases are accepted wherever a name is: `numpy` → `genomeblocks`,
`graph_tool` / `gt` → `graph-tool`, `nx` → `networkx`, `pybedtools` →
`bedtools`, `pyBigWig` → `pybigwig`, `Bio` → `biopython`, `MOODS` → `moods`,
`dict` → `memory`.

### Open handles pick their own engine

A bigWig or FASTA argument may be an already-open object: a `pyBigWig` or
`pybigtools` handle, a `pysam.FastaFile`, a `pyfaidx.Fasta`, a Biopython
`SeqIO.index` or a `{chrom: str}` dict. It is used as it is, with the engine
that opened it. Asking for a different backend on such a call is an error, not
a silent switch:

```python
import pyBigWig
h = pyBigWig.open("atac.bw")
cre.signal([h], verbose=False)                      # fine: reads through pyBigWig
cre.signal([h], backend="python", verbose=False)
# ValueError: the track is an open pybigwig handle, so backend='python' cannot apply:
#   pass the file path instead, or drop backend=
```

---

## Errors name the fix
{: .sec-green }

A requested backend that is not installed raises `ImportError` with the
install command. There is no fallback to another engine:

```python
gb.backends.resolve("intervals", "cgranges")
# ImportError: the 'cgranges' intervals backend is not installed: conda install -c bioconda cgranges
#   (or pip install git+https://github.com/lh3/cgranges)
with gb.use_backend(graph="graph-tool"):
    ...
# ImportError: the 'graph-tool' graph backend is not installed: conda install -c conda-forge graph-tool
gb.backends.resolve("intervals", "bedtools")
# ImportError: the 'bedtools' intervals backend is not installed: pip install pybedtools;
#   conda install -c bioconda bedtools  (the bedtools binary must be on PATH)
```

An unknown name lists the choices; an operation a backend lacks lists the
backends that have it:

```python
gb.backends.resolve("intervals", "nope")
# ValueError: unknown intervals backend 'nope'; choose from: genomeblocks, cgranges, ncls, bioframe, pyranges, bedtools
gb.backends.resolve("nope")
# ValueError: unknown backend family 'nope'; families: intervals, bigwig, motifs, fasta, tables, graph
a.nearest(b, backend="ncls")
# NotImplementedError: the 'ncls' intervals backend has no nearest; use one of: genomeblocks, bioframe, pyranges, bedtools
```

{: .tip }
> `pip install "genomeblocks[interop]"` installs every pip-installable
> engine (ncls, bioframe, pyranges, pybedtools, pyBigWig, MOODS, Biopython,
> pysam, pyfaidx, networkx, igraph). cgranges, bedtools and graph-tool are
> conda extras; see [Installation]({{ '/installation/' | relative_url }}).

---

## Parity: the rules every engine is held to
{: .sec-purple }

Backend differences must not leak into results: no engine-specific nulls, row
orders or dtypes. Each family normalises its engines to one set of rules, and
the test suite runs every public function against every installed backend and
checks they agree (`tests/`, `conftest.installed_backends`).

### Intervals

- **Half-open.** `[s1, e1)` and `[s2, e2)` overlap when `s1 < e2 and s2 < e1`.
  Every engine's pairs are filtered with this rule, so engines that treat ends
  differently still agree.
- **Zero-length intervals** are points between two bases: `[p, p)` overlaps
  `[s, e)` exactly when `s < p < e`. A point at the start of an interval does
  not overlap it. bedtools widens empty intervals on its side and may miss
  pairs that need the widened base; the filter removes anything it adds.
- **Pairs are sorted** by (query row, reference row).
- **Nearest** gives distance 0 to an overlap and otherwise the gap in bases
  (book-ended intervals are 0 apart); `-1` when the chromosome has nothing.
  Among overlapping candidates every engine returns the one with the lowest
  start, then the lowest row. On an exact tie between two equally near
  non-overlapping rows the engines may pick either.
- **Merge** fuses overlapping and book-ended intervals (bedtools `merge`) and
  reports the first row of each block, whatever the engine.

```python
r = gb.as_loci([("chr1", 5, 10)])
inside, at_start = gb.as_loci([("chr1", 7, 7)]), gb.as_loci([("chr1", 5, 5)])
for be in ("genomeblocks", "ncls", "bioframe", "pyranges"):
    inside.overlap_pairs(r, backend=be)[0].tolist()      # -> [0]  on every engine
    at_start.overlap_pairs(r, backend=be)[0].tolist()    # -> []
    r.overlap_rows("chr1", 10, 10, backend=be).tolist()  # -> []   the end is exclusive

q = gb.as_loci([("chr1", 0, 100), ("chr1", 100, 200), ("chr1", 250, 300), ("chr2", 0, 10)])
ref = gb.as_loci([("chr1", 200, 210), ("chr1", 400, 500)])
for be in ("genomeblocks", "bioframe", "pyranges"):
    q.merge(backend=be).to_records()
    # -> [('chr1', 0, 200, '.'), ('chr1', 250, 300, '.'), ('chr2', 0, 10, '.')]
    q.nearest(ref, backend=be)
    # -> (array([ 0,  0,  0, -1]), array([100,   0,  40,  -1]))
```

### bigWig

- **Integer-edge binning.** Bin `b` of a window of `n` bases covers
  `[floor(n·b/n_bins), floor(n·(b+1)/n_bins))`, the edges libBigWig and
  pybigtools use, so fractional bins split the same way everywhere.
- **Stats** are `mean`, `min`, `max`, `sum`, `std` and `coverage`. `std` is the
  population standard deviation; `coverage` is the fraction of bases with
  data. Stats an engine lacks natively are computed from the per-base values
  with the same edges.
- **Missing data.** A bin without data is `missing` (0 by default in
  `signal`); per-base values are NaN where the file has none. Windows that
  leave the chromosome, or name one the file lacks, are clipped and padded the
  same way by every handle.

```python
from genomeblocks.backends.bigwig import open_bigwig
for be in ("pybigtools", "pybigwig", "python"):
    h = open_bigwig("atac.bw", backend=be)
    h.stats_array("chr1", 0, 1003, n_bins=7, stat="mean").round(3)
    # -> [5.205 1.622 3.781 8.95  6.728 6.57  9.355]   identical on all three
    h.stats_array("chr1", 19_950, 20_050, n_bins=2, stat="std", missing=-1)
    # -> [-1. -1.]   past the chromosome end
```

### FASTA

Every source returns the same bases, case kept, clipped to the chromosome.
The indexed reader needs fixed-width records (the `.fai` contract) and says
so, naming `memory` / `pysam` / `pyfaidx` as the way out, when a file is not.

### Motifs

Every engine scans the same `(W × 4)` log-odds matrix,
`log2((count + 0.1) / (total + 0.4) / 0.25)`, lightmotif's own
`normalize(0.1).log_odds()`. A hit is a window position whose score reaches
the threshold; windows containing `N` never hit; a hit that would straddle two
windows is dropped. lightmotif, MOODS and Biopython report the same counts:

```python
from genomeblocks import motifs as gm
{be: gm.scan_motifs_matrix(cre, "genome.fa", "motifs.jaspar", r=50, threshold=7.0,
                           norm=False, backend=be, verbose=False).sum().to_dict()
 for be in ("lightmotif", "moods", "biopython")}
# -> {'lightmotif': {'M1': 8, 'M2': 7}, 'moods': {'M1': 8, 'M2': 7}, 'biopython': {'M1': 8, 'M2': 7}}
```

### Tables

polars and pandas return the same dict of numpy columns: positions as int64,
everything else as strings. Header lines (`#`, `track`, `browser`) are skipped
and `.gz` files read directly by both.
`Loci.make(path, backend="pandas").equals(Loci.make(path, backend="polars"), cols=True)`
holds for every file the suite reads.

### Graph

- **Canonical component labels.** Components are numbered by their first
  row; rows without links get `-1`. graph-tool, scipy, igraph and networkx
  produce the same vector.
- **Parallel edges are summed.** `to_scipy` sums the weights of parallel
  edges into one adjacency cell; the networkx export does the same (one edge
  per pair, numeric columns summed), so PageRank agrees across engines.
  igraph and graph-tool keep the parallel edges and reach the same PageRank.
- **PageRank** is computed to `tol=1e-10` on every engine (scipy runs its own
  power iteration with dangling mass spread uniformly).
- **Layouts are not parity-checked.** `spring` is each engine's own
  force-directed layout (spectral for scipy); only `circular` and `genomic`
  are identical everywhere.

```python
for be in ("scipy", "igraph", "networkx"):
    A.components(backend=be).tolist()         # -> [0, 1, 0, -1, 1, 0, 0]
    A.pagerank("w", backend=be).round(4)      # -> the same vector on every engine
B = gb.Architecture.from_edges(cre, [0, 0, 1], [2, 2, 4], w=[1.0, 2.0, 5.0])   # 0–2 twice
B.to_scipy("w")[0, 2]                         # -> 3.0
B.to_networkx().edges[0, 2]["w"]              # -> 3.0, and number_of_edges() == 2
B.to_igraph().ecount()                        # -> 3
```

---

## Support matrix
{: .sec-purple }

Not every engine does everything. A request for an operation an engine lacks
raises `NotImplementedError` naming the engines that have it.

**intervals**

| Operation | genomeblocks | cgranges | ncls | bioframe | pyranges | bedtools |
|---|---|---|---|---|---|---|
| overlap (`overlap_pairs`, `overlap_any`, `intersect`, `difference`, `Architecture.make`, `Genes.labels`) | yes | yes | yes | yes | yes | yes |
| nearest (`nearest`, `nearest_tss`, `nearest_genes`) | yes | — | — | yes | yes | yes |
| merge (`merge`, `stitch_peaks`) | yes | — | — | yes | yes | yes |
| point lookup (`overlap_rows`, `overlaps`, `A.region`) | yes, cached index | yes, cached index | yes, cached index | via overlap | via overlap | via overlap |

**bigwig**

| | pybigtools | pybigwig | python |
|---|---|---|---|
| `mean`, `min`, `max` | native | native | native |
| `sum`, `coverage` | from per-base values | native | native |
| `std` (population) | from per-base values | from per-base values | native |
| fractional bins, chromosome edges | same reduce path | same reduce path | same reduce path |
| `exact=False` (zoom levels) | yes | yes | yes, when the file has them |
| open handles | yes | yes | `BigWigReader` |

**fasta**: all five sources support `sizes`, `fetch` and `fetch_many`; only
`genomeblocks` and `memory` have no dependency; `memory` is the only one for
gzip files.

**motifs**: all three engines support single-matrix scans and both strands;
MOODS additionally scans many matrices in one pass, which
`scan_motifs_matrix` uses when `backend="moods"`. A `pvalue=` threshold is
turned into a score with lightmotif's exact score distribution for every
engine, so it needs lightmotif installed.

**tables**: polars and pandas parse the BED-like files of `Loci.make`, the
BEDPE of `Pairs.make` and the GTF / GFF3 of `Genes.make`; polars is faster on
large files. Pairs files (`count_pairs`, `read_pairs_chunks`) stream through
pandas in chunks, and the HiChIP short-range reader uses polars when installed
and pandas otherwise, outside the registry.

**graph**: all four engines support `graph()` (a native object), `components`,
`pagerank` and `layout`. `degree` and `strength` need no engine.

---

## Adding a backend
{: .sec-navy }

Each family keeps a small internal interface in its `genomeblocks/backends/*.py`
module. The public API never references a backend object; a new engine is a
registry entry plus the functions that family expects.

1. **Register it** in `_FAMILIES` in `backends/__init__.py`:
   `"name": ("module_to_probe", "pip install ...")`. The module name is what
   `installed()` probes with `importlib.util.find_spec`; `"module+binary"`
   also requires that executable on `PATH` (bedtools); `None` means always
   available. Add the name to `AUTO[family]` only if it gives results identical
   to the default on the whole test matrix.
2. **Implement the family's hooks:**

| Family | What to add |
|---|---|
| `intervals` | a function `(q, r) -> (query rows, reference rows)` in `_PAIRS`; optionally `(q, r) -> (query rows, nearest reference rows)` in `_NEAREST_IMPL` and `(L) -> (codes, starts, ends)` in `_MERGE_IMPL`, with the name added to `_NEAREST` / `_MERGE`. Engines receive the integer chromosome codes as labels (both sides share one Genome), and `overlap_pairs` applies the half-open filter and sorting for you. A cached point index is optional (`L._index_cache(name, build)`). |
| `bigwig` | a `_Handle` subclass with `backend`, `__init__(src)`, `_stats(chrom, start, end, n_bins, stat, exact, missing)` and `_values(chrom, start, end)` (per base, NaN where no data); register it in `_HANDLES` and, if it wraps a third-party handle, teach `handle_backend` to recognise it. `_Handle` adds the clipping and the shared `_reduce` for stats the engine lacks. |
| `fasta` | a `_Source` subclass with `backend`, `sizes()`, `fetch(chrom, start, end)` and optionally `fetch_many`; register it in `_SOURCES` and in `_kind` for open handles. |
| `motifs` | a branch in `Block._prepare` (how the concatenated text is held) and `Block._positions(mat, threshold)` (start positions of hits of one `(W × 4)` log-odds matrix). The split back into windows and strand handling are shared. |
| `tables` | a branch in `read_columns` returning `{column index: numpy array}` with the requested int / float / string dtypes, after `sniff` has counted the header lines. |
| `graph` | branches in `native` (build the engine's graph: vertex `i` = row `i`, edge columns as attributes), `components` (return raw labels; `_canonical` renumbers them), `pagerank` and `layout_edges`. |

3. **Run the parity tests** (`pytest tests/`): every test parametrised over
   `installed_backends(family)` picks the new engine up automatically and
   compares it with the default on random data, zero-length intervals,
   chromosome edges and the planted motif sites.
4. **Time it** with `benchmarks/bench_backends.py`, which runs every installed
   engine of every family through the public API on the same data and records
   the agreement alongside the timings.

The [Benchmarks]({{ '/benchmarks/' | relative_url }}) page shows what each
engine costs or buys; [Design → Tables]({{ '/design/tables/' | relative_url }})
explains the columns the engines receive, and
[API → backends]({{ '/api/backends/' | relative_url }}) lists the registry
functions.
