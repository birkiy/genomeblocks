---
title: Motifs
parent: User Guide
layout: default
nav_order: 8
---

# Motifs
{: .no_toc }

Scan transcription-factor motifs (JASPAR / TRANSFAC / MEME / your own matrices)
over the windows of a `Loci` table, read straight from a FASTA. genomeblocks
parses the motif files, builds the log-odds matrices and the p-value cutoffs
itself, and hands them to a scanning engine: [MOODS](https://github.com/jhkorhonen/MOODS)
by default, [lightmotif](https://github.com/althonos/lightmotif) when only it
is installed, Biopython on request (see [Credits]({{ '/credits/' | relative_url }})).
Every engine reports the same hits.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The pieces
{: .sec-navy }

| piece | what it is |
|---|---|
| `gb.load_motifs(src, format=)` | any motif source → a `Library` (count matrices + names) |
| a FASTA path, dict or open handle | where the sequence comes from (the `fasta` backend family) |
| `scan_motifs_matrix(loci, fasta, motifs)` | hits per locus and motif, a DataFrame aligned to the rows |
| `scan_motifs_profile(...)` | *where* in the window the hits fall: a `(rows, motifs, bins)` cube |
| `bootstrap_enrichment`, `compare_motifs`, `compare_motifs_to_ref` | statistics on those matrices |
| `pwm_distance_matrix`, `cluster_motifs`, `build_archetypes` | collapse a redundant library into archetypes |

Every function takes its loci as anything [`as_loci`]({{ '/interoperability/' | relative_url }})
accepts — a `Loci`, a BED path, a pandas / polars frame, a list of region
strings — and its motifs as anything `load_motifs` accepts. Results are aligned
to the input rows: row *i* of a matrix or cube is locus *i*, and the index is
the locus `uid`.

```python
import genomeblocks as gb
from genomeblocks import motifs as gm
```

{: .warning }
> **A bare `pip install genomeblocks` has no motif engine.** MOODS is C++ and
> lightmotif is Rust, so neither is a core dependency. Install one:
>
> - `conda install -c bioconda moods` — the conda package of genomeblocks already depends on it;
> - `pip install 'genomeblocks[motifs]'` — MOODS-python, compiled at install time;
> - `pip install 'genomeblocks[lightmotif]'` — prebuilt wheels, used automatically when MOODS is absent.
>
> Without one, every `scan_*` call raises
> `ImportError: no motifs backend is installed: conda install -c bioconda moods  (or pip install 'genomeblocks[motifs]', which compiles MOODS-python) or pip install 'genomeblocks[lightmotif]'  (prebuilt wheels)`.
> Reading libraries, the statistics and the archetypes need no engine.

---

## Motif libraries — `load_motifs` and `Library`
{: .sec-navy }

`gb.load_motifs` reads one motif collection into a `Library`: a list of
`(W x 4)` count matrices (columns A C G T) with names and descriptions.

```python
lib = gb.load_motifs("JASPAR2024_CORE.jaspar")                  # raw 4-row counts (JASPAR's .pfm layout)
lib = gb.load_motifs("JASPAR2024_CORE.txt", format="jaspar16")  # bracketed  A [ 1 2 3 ]  layout
lib = gb.load_motifs("matrix.dat", format="transfac")
lib = gb.load_motifs("uniprobe.txt", format="uniprobe")         # probabilities (scaled to 100 sites)
lib = gb.load_motifs("archetypes.meme", format="meme")          # MEME minimal or full
```

| `format` | layout |
|---|---|
| `jaspar` (default) | `>name description` then four rows of counts, A C G T order |
| `jaspar16` | the bracketed JASPAR-2016 layout, `A  [ 10 0 0 ]` |
| `transfac` | TRANSFAC `ID / P0 / 01 ...` blocks |
| `uniprobe` | UniPROBE `A:` `C:` `G:` `T:` probability rows |
| `meme` | MEME `MOTIF` blocks with a `letter-probability matrix` |

All five are parsed by genomeblocks itself (`genomeblocks.backends.motifs`),
not by the scanning engine, so a library reads the same way whichever engine
is installed; `.gz` files are fine.

{: .warning }
> `jaspar` is the **raw four-line count** layout, not the bracketed JASPAR-2016
> one. A bracketed file read as `jaspar` fails with a parse error that names the
> fix — pass `format='jaspar16'`. Counts must be integers in either layout.

Objects work too, so a library can come from anywhere:

```python
import numpy as np
from Bio import motifs as bm

lib = gb.load_motifs(bm.parse(open("jaspar.txt"), "jaspar"))     # Biopython motifs (one or a list)
lib = gb.load_motifs({"CTCF": pfm, "GATA1": pfm2})                # {name: matrix}, 4 x W or W x 4
lib = gb.load_motifs([np.ones((5, 4)), np.ones((6, 4))])           # bare matrices → motif1, motif2
lib = gb.load_motifs(lib)                                          # a Library is returned as is
```

A `{name: matrix}` dict accepts counts or probabilities (a matrix whose rows sum
to 1 is scaled to 100 sites) and either orientation.

### A `Library` is a table

The rest of this page uses a two-motif library — `motifs.jaspar`, the same file
the test suite uses:

```
>M1 TFA
10 0 0 0
0 10 0 0
0 0 10 0
0 0 0 10
>M2 TFB
0 0 10 10
0 0 0 0
10 10 0 0
0 0 0 0
```

```python
lib = gb.load_motifs("motifs.jaspar")
lib
# -> Library(2 motifs: M1, M2)
lib.shape, lib.columns
# -> ((2, 4), ['name', 'description', 'width', 'consensus'])
lib.head()
# ->   name description  width consensus
# -> 0   M1         TFA      4      ACGT
# -> 1   M2         TFB      4      GGAA
lib.describe()
# ->                                                value
# -> motifs                                              2
# -> width min                                           4
# -> width median                                      4.0
# -> width max                                           4
# -> pseudocount                                       0.1
# -> log-odds      log2((count + p) / (total + 4p) / 0.25)
```

`len(lib)`, `list(lib)` (the names), `lib.head()` / `lib.tail()`,
`lib.describe()`, `lib.to_pandas()` / `to_polars()` / `to_arrow()` and the
Arrow / dataframe-interchange protocols are all there, so polars, DuckDB and
seaborn take a `Library` as they are:

```python
import polars as pl, duckdb, seaborn as sns

pl.DataFrame(lib)                                       # shape: (2, 4)
duckdb.sql("select name, width from lib").fetchall()    # [('M1', 4), ('M2', 4)]
sns.barplot(lib, x="name", y="width")
```

Indexing gives matrices:

```python
lib["M1"]                 # the (W x 4) count matrix, by name ...
lib[0]                    # ... or by position
lib.names, lib.descriptions, lib.widths
# -> (['M1', 'M2'], ['TFA', 'TFB'], array([4, 4]))
lib.consensus(1)          # -> 'GGAA'
lib.logodds(0)            # (W x 4) log2-odds matrix — what the scanners score
lib.pfm(0)                # (4 x W) probabilities (rows A C G T), for logos and archetypes
```

Every motif is scored with one log-odds matrix,
`log2((count + 0.1) / (column total + 0.4) / 0.25)`, which genomeblocks
computes (`genomeblocks.backends.motifs.logodds_matrix`: float32, summed in
the same order lightmotif uses, so it matches that library bit for bit) and
every engine scans that same matrix — MOODS, lightmotif and Biopython never
see the counts, only this matrix and a cutoff.

### Selecting motifs

```python
lib.select("tfa")                          # case-insensitive substring over names and descriptions
# -> Library(1 motifs: M1)
lib.select(["M2"], match="exact")          # verbatim names
lib.indices("tf")                          # -> [0, 1]   (positions, same matching rules)
```

`select` searches the description too, so a TF symbol finds its JASPAR
matrix (`lib.select("CTCF")` matches `MA0139.1 CTCF`). `lib["ctcf"]` works the
same way when exactly one motif matches.

### Writing a library out

```python
lib.to_jaspar("lib.jaspar16")              # bracketed JASPAR counts (read back with format='jaspar16')
lib.to_meme("lib.meme")                    # MEME minimal format, probabilities
lib.to_biopython()                         # list[Bio.motifs.Motif] (counts; name = description, matrix_id = name)
lib.to_moods()                             # MOODS matrices: one 4 x W log-odds list of lists per motif
```

`gm.write_meme({name: pfm}, path)` writes any `{name: (4 x W) probabilities}`
dict — the archetypes below, for example — as MEME.

---

## Sequence — the FASTA side
{: .sec-navy }

Scanning reads each window from a FASTA. Pass the path and the fasta backend
does the rest:

```python
M = gm.scan_motifs_matrix(loci, "hg38.fa", "JASPAR2024_CORE.jaspar")
```

| `fasta` argument | how it is read |
|---|---|
| a plain FASTA path | genomeblocks' indexed reader: a samtools `.fai` (built next to the file when missing) and a memory map, so only the requested bases are read |
| a `.gz` FASTA path | gzip cannot be indexed: the file is read into memory once |
| `{chrom: str}` dict | used as is (`gb.read_fasta(path)` makes one) |
| an open `pyfaidx.Fasta`, `pysam.FastaFile` or `Bio.SeqIO.index` | used as is |

The other readers are backends — `pysam`, `pyfaidx`, `memory`, `biopython` —
chosen per call with `fasta_backend=` (or `backend=` on `Loci.sequences`), or
for a block with `gb.use_backend(fasta=...)`. All return the same bases.

```python
seqs = gb.read_fasta("hg38.fa")             # {chrom: str}, the whole file in memory
genome = gb.Genome.from_fasta("hg38.fa")    # chromosome names and sizes from the .fai
genome.sizes                                # {'chr1': <length>, 'chr2': <length>, ...}

L = gb.as_loci(["chr1:950-1050", "chr1:4950-5050", "chr1:7000-7100"])
L.sequences("genome.fa")                    # one string per row, as stored
L.sequences("genome.fa", r=6)               # centre ± 6 bp
# -> ['ATCAGCACGTAC', 'TCGTTAGGAAGG', 'CGTATATGAGGA']
L.sequences("genome.fa", r=6, backend="pysam")       # same bases, another reader
L.sequences("genome.fa", strand=True, upper=True)    # reverse-complement '-' rows, upper-case
L.to_seqrecords("genome.fa", r=6)           # Biopython SeqRecords, id = uid
L.to_fasta("windows.fa", "genome.fa", r=6)  # write the windows as FASTA
```

Windows are clipped to the chromosome, so a window at a chromosome end comes
back shorter; the scanners treat such windows (and chromosomes missing from the
FASTA) as zero hits and say how many rows that affected.

{: .note }
> The indexed reader needs fixed-width records (every line of a chromosome the
> same length, as `samtools faidx` does). A FASTA with uneven lines raises an
> error that names the fix: reformat it, or pass `backend='memory'`, `'pysam'`
> or `'pyfaidx'`.

---

## Per-locus x per-motif matrix — `scan_motifs_matrix`
{: .sec-navy }

The workhorse. One row per locus, one column per motif, each cell the number
of positions in the locus's **centre ± `r`** window whose log-odds score
reaches the threshold.

```python
M = gm.scan_motifs_matrix(
    loci, "hg38.fa", "JASPAR2024_CORE.jaspar",
    r=250,                 # window = centre ± r (500 bp)
    threshold=13.0,        # log2-odds cutoff; or pvalue=1e-4 for a per-motif cutoff
    norm=True,             # divide counts by motif width
    both_strands=False,    # also count reverse-strand matches
    workers=None,          # processes over motifs; None = half the cores
    backend=None,          # 'moods' (default; lightmotif when MOODS is absent), 'lightmotif', 'biopython'
    fasta_backend=None,    # 'genomeblocks' (default), 'pysam', 'pyfaidx', 'memory', 'biopython'
)
# -> pandas.DataFrame, index = uid (in the loci's order), columns = motif names
```

The same call is a method on `Loci`: `loci.scan_motifs_matrix(fasta, motifs, **kw)`.

On the two-motif library and three windows around the planted sites
(`ACGT` x 3 at chr1:1000, `GGAA` x 3 at chr1:5000):

```python
L = gb.as_loci(["chr1:950-1050", "chr1:4950-5050", "chr1:7000-7100"])
M = gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
M
# ->                    M1  M2
# -> uid
# -> chr1:950-1050(.)    3   1
# -> chr1:4950-5050(.)   0   3
# -> chr1:7000-7100(.)   2   0
```

### Threshold or p-value

`threshold` is a log2-odds score: a scalar for every motif, or one value per
motif. `pvalue` replaces it with a per-motif cutoff — the score each motif
reaches with that probability under a uniform background — so motifs of
different widths and information content are cut at a comparable stringency.
The cutoff comes from genomeblocks' own `threshold_from_pvalue`: the exact
score distribution of the log-odds matrix on a 0.001-bit grid (it agrees with
MOODS' `threshold_from_p` to a few thousandths of a bit), computed once per
motif and handed to whichever engine scans, so the hits do not depend on the
engine.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/motifs-threshold.svg %}
</div><figcaption>
<strong>One matrix, one cutoff, every engine.</strong> The count matrix becomes a log-odds matrix through <code>logodds_matrix()</code> — <code>log2((c + 0.1) / (total + 0.4) / 0.25)</code> in float32, summed in lightmotif's order — and that is the only matrix any engine sees. <code>threshold_from_pvalue()</code> builds the matrix's exact score distribution on a 0.001-bit grid and returns the smallest score whose upper tail is at most <code>p</code>: 7.49 bits for <code>p = 0.001</code> on this matrix (MOODS' <code>threshold_from_p</code> gives 7.4875). MOODS, lightmotif and Biopython each receive the same matrix and the same cutoff, so their hits are identical; <code>threshold=</code> skips the distribution and uses the number as given.
</figcaption></figure>

```python
gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=[7.0, 5.0], norm=False)   # per motif
gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, pvalue=0.01, both_strands=True)       # per-motif cutoff
```

A 4-mer cannot reach `pvalue=1e-3` (one exact match has p = 1/256), so a
very strict p-value on a short motif gives zero hits — pick the p-value for
the widths in your library.

### `norm`, `both_strands`, `workers`

- `norm=True` (default) divides each count by the motif width, so long motifs
  do not dominate a sum over the library. Set `norm=False` for raw counts.
- `both_strands=True` also scans the reverse complement of each matrix and
  adds those hits; a hit that straddles two windows is dropped.
- `workers` spreads the motifs over a process pool. `None` means half the
  cores; small jobs (fewer than 16 motifs or 200 windows) stay serial because
  process start-up would dominate.

Windows that run off a chromosome end or sit on a chromosome missing from the
FASTA count 0 and are reported:

```python
gm.scan_motifs_matrix(gb.as_loci(["chr1:10-30", "chrZ:100-200"]), "genome.fa", "motifs.jaspar", r=50)
# -> [motifs] 2 of 2 windows run off a chromosome or are missing from the FASTA; their rows count 0.
```

---

## Totals per motif — `scan_motifs`
{: .sec-navy }

When you only need the per-motif total over all windows:

```python
counts = gm.scan_motifs(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
counts
# -> {'M1': 5.0, 'M2': 4.0}
import pandas as pd
pd.Series(counts).sort_values(ascending=False).head(20)
```

It is `scan_motifs_matrix(...).sum(axis=0)` with the same arguments (minus
`workers`), also available as `loci.scan_motifs(fasta, motifs)`.

---

## What the scan does
{: .sec-navy }

1. `load_motifs` turns the source into a `Library`; each motif gets one
   log-odds matrix and one cutoff (`threshold`, or `pvalue` through
   `threshold_from_pvalue`), both computed by genomeblocks.
2. Every window (centre ± `r`) is read from the FASTA once, upper-cased, and
   the windows are laid end to end into one block.
3. Each motif is scanned once over that block by the engine — MOODS scans the
   whole library in one pass, lightmotif and Biopython one matrix at a time;
   hits are split back per window (a hit spanning two windows is dropped), so
   the result equals scanning each window on its own.
4. Counts are divided by the motif width when `norm=True` and returned as a
   DataFrame indexed by `uid`.

A window containing a letter other than A C G T N never hits; `N` positions
score `-inf`.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/motifs-block.svg %}
</div><figcaption>
<strong>One engine call per motif over all windows.</strong> The windows are joined into one block once; each motif is scanned over the whole block (MOODS takes the whole library in one call), and every hit is mapped back to its window by position, with hits that straddle a boundary dropped. The result equals scanning each window alone, at a fraction of the calls; with <code>workers</code> the motifs are split across processes that each hold the block.
</figcaption></figure>

---

## Masking an anchor motif — `scan_motifs_matrix_masked`
{: .sec-navy }

"Which co-factors set these regions apart *besides* CTCF?" Mask every match of
the anchor motifs before scanning the library:

```python
masked = gm.scan_motifs_matrix_masked(
    loci, "hg38.fa", "JASPAR2024_CORE.jaspar",
    anchors=["CTCF"],        # case-insensitive substring against motif names and descriptions
    r=250, threshold=13.0,
    window=10,               # ±bp around each anchor hit replaced by random bases
    anchor_threshold=None,   # a separate cutoff for finding the anchor hits (default: threshold / pvalue)
    skip_anchors=True,       # leave the anchor motifs out of the result
    seed=0,
)
```

On the small library, masking `M1` removes its `ACGT` sites and the column:

```python
gm.scan_motifs_matrix_masked(L, "genome.fa", "motifs.jaspar", ["M1"], r=50, threshold=7.0, norm=False,
                             seed=0, verbose=False)
# ->                    M2
# -> uid
# -> chr1:950-1050(.)    1
# -> chr1:4950-5050(.)   3
# -> chr1:7000-7100(.)   1
```

(The third window had two `M1` sites; the random bases that replaced them
happen to contain one `GGAA`.) `threshold` may be a list with one value per library motif, and `pvalue` works
as in `scan_motifs_matrix`. With `skip_anchors=False` the anchor columns stay
(and show the hits left after masking).

{: .warning }
> **Mask the reference the same way** before a differential test. A masked
> query against an unmasked pool looks depleted for every motif that
> co-occurs with the anchor, for reasons unrelated to biology.

---

## Positional profiles — `scan_motifs_profile` and `plot_motif_heatmap`
{: .sec-navy }

Where in the window do the hits fall? `scan_motifs_profile` returns a
`(rows, motifs, bins)` cube — the layout of a [signal cube]({{ '/guide/signal/' | relative_url }}),
so the same heatmap code draws it — plus the names of the motifs on the middle
axis.

```python
M, names = gm.scan_motifs_profile(
    loci, "hg38.fa", "JASPAR2024_CORE.jaspar",
    select=["CTCF", "GATA", "SOX"],   # motifs by name (substring; match='exact' for verbatim); None = all
    r=500, n_bins=100,                # window centre ± r, split into n_bins (10 bp/bin here)
    threshold=13.0,                   # or pvalue=
    both_strands=True,                # default here: reverse matches in forward coordinates
    smooth=1.0,                       # Gaussian sigma in bins along each row; 0 = raw counts
    norm=False,                       # divide by motif width
)
fig = gb.plot_motif_heatmap(loci, M, names, r=500, groups={"A": a, "B": b})
```

Each hit is credited to the bin of its centre. On the planted sites:

```python
M, names = gm.scan_motifs_profile(L, "genome.fa", "motifs.jaspar", ["M1", "M2"], r=60, n_bins=12,
                                  threshold=7.0, smooth=0, verbose=False)
M.shape, names
# -> ((3, 2, 12), ['M1', 'M2'])
M[0, 0]                      # M1 hits along the first window: three overlapping ACGT at the centre
# -> array([0., 0., 0., 0., 0., 0., 4., 2., 0., 0., 0., 0.], dtype=float32)
```

`gb.plot_motif_heatmap(loci, M, names, *, r=500, groups=None, vmax=None, ymax=None, cmap="Purples")`
draws one column per motif with `signal_draw.plot_heatmap`: `groups` splits the
rows like the signal heatmap (`{name: Loci | mask | row numbers}`), `r` labels
the x axis, and `vmax` / `ymax` default to a per-motif scale (98th percentile
of the non-zero cells; the top of the group mean profiles). The method form is
`loci.scan_motifs_profile(fasta, motifs, select, **kw)`.

---

## Enrichment and differential tests
{: .sec-navy }

All three take the matrices `scan_motifs_matrix` (or the masked variant)
returns. Columns are aligned by name, so the matrices need not share a motif
set — missing columns count 0.

### Bootstrapped log fold change — `bootstrap_enrichment`

Point estimates of each group's mean hit rate against a reference pool, with
bootstrap resampling for stability:

```python
mat_a    = gm.scan_motifs_matrix(a_cre,   "hg38.fa", motif_db, r=250)
mat_b    = gm.scan_motifs_matrix(b_cre,   "hg38.fa", motif_db, r=250)
mat_pool = gm.scan_motifs_matrix(all_cre, "hg38.fa", motif_db, r=250)

enr = gm.bootstrap_enrichment({"A": mat_a, "B": mat_b}, ref=mat_pool,
                              boot=100, sample=500, pseudo=0.1, seed=0)
```

One row per motif: `Factor`, `mean_ref`, then `mean_<g>` and
`LFC_<g> = log2((mean_g + pseudo) / (mean_ref + pseudo))` per group. With
exactly two groups there is also `LFC = LFC_<first> - LFC_<second>`, the
motif's preference for one set over the other. `boot` iterations each draw
`sample` rows without replacement (capped at the group size).

```python
enr.sort_values("LFC", ascending=False)[["Factor", "LFC", "LFC_A", "LFC_B"]]
```

### Two sets against each other — `compare_motifs`

A per-motif Mann-Whitney U test on the per-locus counts (rank-based, right for
zero-inflated counts) plus a log2 fold change of the means, with
Benjamini-Hochberg FDR:

```python
diff = gm.compare_motifs(mat_a, mat_b, pseudo=0.1, alternative="two-sided")
diff.head()
# columns: Factor, mean_A, mean_B, LFC, U, p, p_adj  — sorted by LFC
```

### Sets against a background — `compare_motifs_to_ref`

The same test of one or several query sets against a reference pool, when you
want per-motif significance rather than the bootstrap's point estimates:

```python
res = gm.compare_motifs_to_ref({"A": mat_a, "B": mat_b}, ref=mat_pool)
# columns: Factor, mean_ref, then per group: mean_<g>, LFC_<g>, U_<g>, p_<g>, p_adj_<g>;
# with two groups also LFC = LFC_A - LFC_B
res = gm.compare_motifs_to_ref(mat_a, ref=mat_pool)        # one query: the group is called 'query'
```

`pseudo` is added to both means before the log ratio; set it near the median
non-zero column mean so it does not dominate small rates.

---

## Motif archetypes — clustering to consensus matrices
{: .sec-navy }

A motif library is redundant (dozens of forkhead matrices). Collapse it:
pairwise Sandelin-Wasserman similarity → hierarchical clustering → one
consensus matrix per cluster.

```python
res = gm.build_archetypes("JASPAR2024_CORE.jaspar", cutoff=0.3, workers=None)
res.keys()
# -> dict_keys(['D', 'Z', 'labels', 'names', 'pwms', 'archetypes', 'members'])
res["archetypes"]           # {'ARCH_01': (4 x W) PFM, ...}
res["members"]              # {'ARCH_01': ['MA0139.1 CTCF', ...], ...}
gm.write_meme(res["archetypes"], "archetypes.meme")       # scan with the archetypes instead
```

`cutoff` is a distance (`1 - similarity`): motifs with similarity ≥
`1 - cutoff` merge. 0.25–0.35 is a typical range; tune it on the dendrogram.
The steps are available on their own:

```python
D, names, pwms = gm.pwm_distance_matrix("motifs.jaspar", min_overlap=3)   # (n x n) distances in [0, 1]
labels, Z = gm.cluster_motifs(D, cutoff=0.5, linkage_method="average")    # 1-based cluster ids, linkage matrix
pfm = gm.archetype([pwms[0], pwms[1]], min_overlap=3)                     # one consensus (4 x W) from members
pfm.shape
# -> (4, 5)
```

`archetype` picks the medoid, aligns every member to it (best offset and
strand) and averages the aligned columns, weighting each member by its mean
information content (`weight_by_ic=True`).

A consensus from a named subset — say the top hits of a differential test —
without clustering the whole library:

```python
fox = gm.archetype_from_names("JASPAR2024_CORE.jaspar", ["FOXA1", "FOXA2", "FOXA3"])
fox["archetype"], fox["members"], fox["pwms"]      # the (4 x W) consensus, the matched names, their PFMs
```

### Logos

`genomeblocks.motifs_draw` renders logos with `logomaker`
(`pip install genomeblocks[viz]`):

```python
from genomeblocks.motifs_draw import plot_archetype, plot_archetypes, plot_cluster_members, plot_dendrogram

plot_archetype(lib.pfm(0), title=lib.names[0])                      # one (4 x W) PFM as a bits logo
plot_archetypes(res["archetypes"], members=res["members"])          # a grid, sized by cluster
plot_cluster_members(res["archetypes"]["ARCH_01"], member_pwms, member_names)   # archetype + members for QC
plot_dendrogram(Z, names=names, cutoff=0.5)                          # the clustering with the cutoff line
```

---

## The three engines agree
{: .sec-navy }

The `motifs` backend family has three members. genomeblocks parses the
library, builds the log-odds matrices and turns a p-value into a cutoff; the
engine only slides those matrices along the block. So all three report the
same hits, and they differ in speed and in what they need installed.

| backend | engine | install |
|---|---|---|
| `moods` (default) | MOODS, C++ — scans the whole library in one pass | `conda install -c bioconda moods`, or `pip install 'genomeblocks[motifs]'` (MOODS-python, compiled at install) |
| `lightmotif` | SIMD scanner, Rust — one matrix per call | `pip install 'genomeblocks[lightmotif]'` (prebuilt wheels); picked automatically only when MOODS is absent |
| `biopython` | `Bio.motifs` PSSM, numpy | `pip install biopython`; never picked automatically |

With `backend=None` the choice is automatic between the two engines that are
interchangeable (`moods`, then `lightmotif`); Biopython is always a request.

```python
ref = gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, verbose=False)
for b in ["moods", "lightmotif", "biopython"]:
    assert gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, backend=b,
                                 verbose=False).equals(ref)

with gb.use_backend(motifs="lightmotif", fasta="pyfaidx"):     # for a block of code
    M = gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, verbose=False)
M.equals(ref)
# -> True
gb.backends().query("family == 'motifs'")
# ->     family     backend  installed  default  in use install
# -> 9   motifs       moods       True     True    True
# -> 10  motifs  lightmotif       True    False   False
# -> 11  motifs   biopython       True    False   False
```

A backend that is requested but not installed raises `ImportError` with its
install command (`the 'lightmotif' motifs backend is not installed: pip install
'genomeblocks[lightmotif]'  (prebuilt wheels)`); an unknown name raises
`ValueError` listing the choices; with no engine at all the error names both
extras (see the top of this page). There is no silent switch. `gb.backends()`
shows what is installed and which engine each family uses; see
[Backends]({{ '/backends/' | relative_url }}).

---

## Tips
{: .sec-navy }

- `threshold=13.0` is a reasonable default for JASPAR counts (log2-odds in
  bits). Pick an empirical value by scanning a shuffled control, or use
  `pvalue=` so every motif is cut at the same stringency.
- Motif counts scale roughly with `r`; use the same `r` in every downstream
  analysis so sets stay comparable.
- Scanning the background pool is usually the slow step (it is the largest
  set). `workers` parallelises over motifs; subsample a pool of hundreds of
  thousands of regions to a few thousand first.
- `loci.sequences(fasta, r=...)` is the windowed slice the scanners use. For
  repeated scans, load the motifs once (`lib = gb.load_motifs(path)`) and pass
  the `Library` instead of the path; the FASTA index is built once and read
  back from the `.fai` on every later call.
- See [API → motifs]({{ '/api/motifs/' | relative_url }}) for every signature
  and [Design → Motifs]({{ '/design/motifs/' | relative_url }}) for how the
  block scanner works.
