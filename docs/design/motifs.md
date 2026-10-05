---
title: Motifs
parent: Design
layout: default
nav_order: 4
---

# Motifs
{: .no_toc }

`scan_motifs_matrix` counts position-weight-matrix hits in a fixed window
around every interval, for every motif in a `Library`, and returns a rows ×
motifs table aligned to the input. Sequence comes from the `fasta` backend,
scoring from the `motifs` backend, and every engine reports the same hits.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## The Library: one matrix format for every source
{: .sec-navy }

`gb.load_motifs(src, format=)` turns any motif source into a `Library`: motif
names, descriptions and one `(W × 4)` count matrix per motif in A C G T
order.

| Source | How it is read |
|---|---|
| JASPAR (raw counts), `jaspar16` (bracketed), TRANSFAC, uniprobe | lightmotif's parsers; uniprobe probabilities become counts at 100 sites |
| MEME (minimal or full) | the package's own reader; probabilities × `nsites` (100 when absent) |
| Biopython `Bio.motifs.Motif` objects, lightmotif motifs | their count dictionaries / arrays |
| `{name: matrix}` dicts, lists of arrays, one array | `4 × W` or `W × 4`; a matrix whose rows sum to 1 is scaled to 100 sites |

Each motif is scored with one log-odds matrix,
`log2((count + 0.1) / (column total + 0.4) / 0.25)`, which is lightmotif's
own `normalize(0.1).log_odds()`. Every engine scans that same matrix, so
their hits agree; `lib.logodds(i)` is cached per motif.

```python
import genomeblocks as gb

lib = gb.load_motifs("motifs.jaspar")
lib
# -> Library(2 motifs: M1, M2)
lib.to_pandas()
# ->   name description  width consensus
# -> 0   M1         TFA      4      ACGT
# -> 1   M2         TFB      4      GGAA
lib.logodds(0).round(2)[0]
# -> array([ 1.96, -4.7 , -4.7 , -4.7 ])
gb.load_motifs(lib.to_biopython()).names
# -> ['M1', 'M2']
```

As a table the `Library` is `(name, description, width, consensus)`, with
`select(names)` (case-insensitive substring over names and descriptions, or
`match="exact"`), `take`, `lib["CTCF"]`, and `to_jaspar` / `to_meme` /
`to_biopython` / `to_moods` out.

## Windows come from the fasta backend
{: .sec-navy }

Each interval contributes `centre ± r` (`2r` bases). `Loci.sequences(fasta,
r=)` reads them through `backends.fasta.open_fasta`, which wraps a path, a
`{chrom: str}` dict, or an open pyfaidx / pysam / Biopython handle in one
`fetch_many` interface:

| Backend | How |
|---|---|
| `genomeblocks` (default) | the samtools `.fai` index (built next to the file when missing) plus a memory map: only the requested bases are read |
| `pysam`, `pyfaidx` | the libraries' indexed readers |
| `memory` | the whole file in a dict (what a gzip FASTA gets, since it cannot be indexed) |
| `biopython` | `Bio.SeqIO.index` |

All return the same bases. A window that runs off a chromosome, or sits on a
chromosome the FASTA lacks, comes back shorter than `2r`; such rows are kept
in the result and count 0 for every motif (as a signal cube leaves such rows
at 0), and a message says how many there were. Windows containing letters
other than A C G T N are skipped by the scanner and count 0 too.

```python
from genomeblocks.backends.fasta import open_fasta
for b in ("genomeblocks", "pysam", "pyfaidx", "memory", "biopython"):
    print(b, open_fasta("genome.fa", backend=b).fetch("chr1", 1000, 1012))
# -> genomeblocks ACGTACGTACGT
# -> pysam ACGTACGTACGT
# -> pyfaidx ACGTACGTACGT
# -> memory ACGTACGTACGT
# -> biopython ACGTACGTACGT
```

## Scanning many short windows as one block
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/motifs-block.svg %}
</div><figcaption>
<strong>One scan per motif, not per window.</strong> A 500-bp window is so short that calling the scanner once per window costs more in call overhead than in scoring. The windows are joined once into a <code>Block</code>, and each motif is scanned over the whole block in a single engine call. Each hit is assigned to its window by binary search on the window offsets; a hit that would span two windows is dropped, so the counts equal scanning each window on its own.
</figcaption></figure>

`backends.motifs.Block(seqs, backend)` keeps the valid windows' text, their
offsets and the map back to input rows. `block.hits(mat, threshold,
both)` returns `(row, position in window, strand)` for every hit;
`block.counts` is the `bincount`.

| Backend | Engine call | Notes |
|---|---|---|
| `lightmotif` (default) | `lightmotif.scan(ScoringMatrix, striped block, threshold=)` | SIMD; the striped sequence is prepared once per block |
| `moods` | `MOODS.scan.Scanner.scan` | scans all matrices of a batch in one pass (`moods_counts`) |
| `biopython` | `PositionSpecificScoringMatrix.calculate` | numpy; no compiled dependency beyond Biopython |

```python
from genomeblocks.backends.motifs import Block
L = gb.as_loci([("chr1", 950, 1050), ("chr1", 4950, 5050), ("chr1", 7000, 7100)])
seqs = L.sequences("genome.fa", r=50, upper=True)
blk = Block(seqs)
blk.n, blk.offsets, len(blk.text)
# -> (3, array([  0, 100, 200, 300]), 300)
blk.hits(lib.logodds(0), 7.0)              # rows, positions, strands
# -> (array([0, 2, 0, 2, 0]), array([58, 35, 54, 43, 50]), array([0, 0, 0, 0, 0], dtype=int8))
for b in ("lightmotif", "moods", "biopython"):
    print(b, Block(seqs, b).counts(lib.logodds(0), 7.0).tolist())
# -> lightmotif [3, 0, 2]
# -> moods [3, 0, 2]
# -> biopython [3, 0, 2]
```

- **Thresholds.** A hit is any position scoring at least `threshold` (13
  log2-odds by default, scalar or one per motif); `pvalue=` instead sets a
  per-motif cutoff from lightmotif's exact score distribution under a uniform
  background. `norm=True` divides counts by motif width; `both_strands=True`
  also scans the reverse complement.
- **Parallelism.** Motifs are independent, so they are split into batches
  across a process pool whose initializer builds the block once per worker.
  Jobs with fewer than 16 motifs or 200 windows run serially (`_workers`),
  because process start-up would dominate.
- **Output.** A DataFrame with one row per input interval (index = uid, in
  the input's order) and one column per motif.

```python
from genomeblocks import motifs as gm
M = gm.scan_motifs_matrix(L, "genome.fa", lib, r=50, threshold=7.0, norm=False, verbose=False)
M
# ->                    M1  M2
# -> uid
# -> chr1:950-1050(.)    3   1
# -> chr1:4950-5050(.)   0   3
# -> chr1:7000-7100(.)   2   0
M.equals(gm.scan_motifs_matrix(L, "genome.fa", lib, r=50, threshold=7.0, norm=False, backend="moods", verbose=False))
# -> True
edge = gb.as_loci([("chr1", 10, 60), ("chr1", 950, 1050), ("chrZ", 0, 100)])
gm.scan_motifs_matrix(edge, "genome.fa", lib, r=50, threshold=7.0, norm=False).M1.tolist()
# -> [motifs] 2 of 3 windows run off a chromosome or are missing from the FASTA; their rows count 0.
# -> [0, 3, 0]
```

## Masking an anchor motif
{: .sec-navy }

`scan_motifs_matrix_masked(..., anchors=["CTCF"])` finds every hit of the
anchor motifs on the block first (`block.hits`), overwrites `centre ±
window` bp around each with random A/C/G/T (seeded by `seed=`), and then scans
the library on the masked windows. Co-factors are counted without the
anchor's own sites, and the anchor motifs are left out of the result unless
`skip_anchors=False`. `anchor_threshold=` sets a separate cutoff for finding
the anchors.

```python
gm.scan_motifs_matrix_masked(L, "genome.fa", lib, ["M1"], r=50, threshold=7.0, norm=False, seed=0, verbose=False).columns.tolist()
# -> ['M2']
```

## Positional profiles
{: .sec-navy }

`scan_motifs_profile(loci, fasta, motifs, select=[...], r=500, n_bins=100)`
keeps the position of every hit: each window is split into `n_bins` bins and
a hit is credited to the bin of its centre (both strands by default, reverse
matches in forward coordinates). The result is a `(rows, motifs, bins)`
float32 cube with the layout of a signal cube, smoothed along the bins by a
Gaussian of `smooth` bins, plus the motif names of its middle axis, so
`gb.plot_motif_heatmap(loci, M, names, r=)` draws it with
`signal_draw.plot_heatmap`.

```python
prof, names = gm.scan_motifs_profile(L, "genome.fa", lib, ["M1"], r=60, n_bins=12, threshold=7.0, smooth=0, verbose=False)
prof.shape, names, prof[0, 0]
# -> ((3, 1, 12), ['M1'], array([0., 0., 0., 0., 0., 0., 4., 2., 0., 0., 0., 0.], dtype=float32))
fig = gb.plot_motif_heatmap(L, prof, names, r=60)
```

## From matrices to enrichment
{: .sec-navy }

- `bootstrap_enrichment(groups, ref, boot=, sample=)` resamples rows of each
  matrix and reports per-motif log2 fold change of the bootstrap means
  against a reference set (and `LFC = LFC_a − LFC_b` when exactly two groups
  are given).
- `compare_motifs(mat_a, mat_b)` runs a per-motif Mann-Whitney U test on the
  per-sequence counts, with the log2 fold change of the means and
  Benjamini-Hochberg `p_adj`.
- `compare_motifs_to_ref(query, ref)` does the same for one or several query
  sets against a reference pool.

All three take the DataFrames `scan_motifs_matrix` returns, align the motif
columns, and work on arrays.

## Motif archetypes
{: .sec-navy }

Large libraries repeat themselves (dozens of near-identical FOX or ETS
matrices). `pwm_distance_matrix(motifs)` scores every pair of PWMs with the
Sandelin-Wasserman similarity over all offsets and both strands (parallel
over pairs when there are enough), `cluster_motifs(D, cutoff=)` cuts a
hierarchical clustering, `archetype(pwms)` aligns each cluster to its medoid
into one information-weighted consensus PFM, and `build_archetypes` runs the
whole chain. `archetype_from_names` builds one archetype from a named subset.
`Library.to_meme` / `write_meme` export for MEME tools, and `motifs_draw`
draws logos, cluster members and the dendrogram (`logomaker` is optional).

```python
D, names, pwms = gm.pwm_distance_matrix(lib, min_overlap=3, verbose=False)
D.round(3), names
# -> (array([[0.   , 0.661], [0.661, 0.   ]]), ['M1', 'M2'])
arch = gm.build_archetypes(lib, min_overlap=3, cutoff=0.5, verbose=False)
arch["members"]
# -> {'ARCH_01': ['M1'], 'ARCH_02': ['M2']}
```

## Costs
{: .sec-navy }

- **Sequence:** one `fetch_many` over the valid windows, reading only those
  bases with the indexed reader; the `.fai` is built once.
- **Scanning:** one engine call per motif over the whole block, then one
  `searchsorted` and one `bincount` per motif; the block itself is built once
  per worker.
- **Thresholds from p-values:** one exact score distribution per motif.
- **Masking:** one extra pass of `hits` per anchor motif before the scan.
- **Profiles:** one `add.at` per motif into the cube, one Gaussian filter
  along the bins.

Measured numbers are on the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page.
