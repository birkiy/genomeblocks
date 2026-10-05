---
title: motifs
parent: API Reference
layout: default
nav_order: 9
---

# `genomeblocks.motifs`
{: .no_toc }

Motif scanning over loci windows, enrichment statistics on the resulting
count matrices, and motif clustering into archetypes. Scanning reads each
window (centre ± `r`) through the fasta backend and scores a motif library
with the motifs backend (MOODS by default, lightmotif when MOODS is absent,
Biopython on request); every engine scores the same log-odds matrices and
reports the same hits. See the [Motifs guide]({{ '/guide/motifs/' | relative_url }}) for the
workflow.
{: .fs-5 .fw-300 }

```python
import genomeblocks as gb
from genomeblocks import motifs as gm
```

Motif libraries (`Library`, `gb.load_motifs`, the JASPAR / jaspar16 / TRANSFAC
/ uniprobe / MEME readers, `to_meme` / `to_jaspar` / `to_biopython` /
`to_moods`) are documented on the
[backends page]({{ '/api/backends/' | relative_url }}); `genomeblocks.motifs`
re-exports `Library`, `load_motifs` and `write_meme`.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Common arguments
{: .sec-navy }

Every scanning function takes the same first three arguments and shares
these keywords.

| Argument | Meaning |
|---|---|
| `loci` | Anything [`as_loci`]({{ '/api/interop/' | relative_url }}) takes. Results line up with its rows. |
| `fasta` | A FASTA path (indexed on first use), a `{chrom: sequence}` dict (what `gb.read_fasta` returns), or an open pyfaidx / pysam / Biopython handle. |
| `motifs` | A motif file, a `Library`, Biopython motifs, or `(4, W)` / `(W, 4)` matrices (a list or a `{name: matrix}` dict). |
| `format` | The file format when `motifs` is a path: `'jaspar'` (default, raw 4-line counts), `'jaspar16'` (bracketed), `'transfac'`, `'uniprobe'`, `'meme'`. |
| `r` | Half-window around each locus centre; windows are `2r` bp. |
| `threshold` | Log2-odds cutoff, a scalar or one value per motif. |
| `pvalue` | Instead of `threshold`: the score each motif reaches with that probability under a uniform background — the exact score distribution computed by [`threshold_from_pvalue`]({{ '/api/backends/' | relative_url }}), the same cutoff whichever engine scans — one value per motif. |
| `norm` | Divide counts by motif width. |
| `both_strands` | Also count reverse-strand matches. |
| `workers` | Processes over motifs; `None` = half the cores. Small jobs (fewer than 16 motifs or 200 windows) stay serial. |
| `backend` | The [motifs backend]({{ '/backends/' | relative_url }}): `'lightmotif'`, `'moods'` or `'biopython'`. |
| `verbose` | Progress bars and a summary line. |

Windows that run off a chromosome end, sit on a chromosome missing from the
FASTA, or contain letters other than A C G T N count 0; `verbose` reports how
many.

The examples below use a FASTA with `ACGT` planted three times at chr1:1000
and `GGAA` three times at chr1:5000, and a two-motif JASPAR file (`M1` =
ACGT, `M2` = GGAA).

```python
L = gb.as_loci([("chr1", 950, 1050), ("chr1", 4950, 5050), ("chr1", 7000, 7100)])
```

---

## Scanning
{: .sec-navy }

| Function | Returns | One line |
|---|---|---|
| `scan_motifs` | `dict[str, float]` | Total hits per motif over all windows. |
| `scan_motifs_matrix` | `DataFrame` (loci x motifs) | Hits of every motif in every window. |
| `scan_motifs_matrix_masked` | `DataFrame` | The same after hiding the matches of anchor motifs. |
| `scan_motifs_profile` | `(ndarray, names)` | Where the hits fall in each window: a `(rows x motifs x bins)` cube. |

All four are also `Loci` methods: `L.scan_motifs(fasta, motifs, **kw)`,
`L.scan_motifs_matrix(...)`, `L.scan_motifs_matrix_masked(fasta, motifs, anchors, **kw)`,
`L.scan_motifs_profile(fasta, motifs, select=None, **kw)`.

### `scan_motifs`

```python
gm.scan_motifs(loci, fasta, motifs, *, format="jaspar", r=250, threshold=13.0, pvalue=None,
               norm=True, both_strands=False, backend=None, verbose=True) -> dict[str, float]
```

```python
gm.scan_motifs(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
# -> {'M1': 5.0, 'M2': 4.0}
gm.scan_motifs(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, verbose=False)   # norm: / width 4
# -> {'M1': 1.25, 'M2': 1.0}
```

### `scan_motifs_matrix`

```python
gm.scan_motifs_matrix(loci, fasta, motifs, *, format="jaspar", r=250, threshold=13.0, pvalue=None,
                      norm=True, both_strands=False, workers=None, backend=None,
                      fasta_backend=None, verbose=True) -> pandas.DataFrame
```

One row per locus (index `uid`, in the loci's order), one column per motif.
Windows are read once; motifs are spread over a process pool.
`fasta_backend` picks the [fasta backend]({{ '/backends/' | relative_url }})
(`'genomeblocks'`, `'pysam'`, `'pyfaidx'`, `'memory'`, `'biopython'`).

```python
M = gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
M
# ->                    M1  M2
#    uid
#    chr1:950-1050(.)    3   1
#    chr1:4950-5050(.)   0   3
#    chr1:7000-7100(.)   2   0

gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, pvalue=0.01, norm=False, both_strands=True, verbose=False)
# ->                    M1  M2
#    uid
#    chr1:950-1050(.)    6   1
#    chr1:4950-5050(.)   0   3
#    chr1:7000-7100(.)   4   0

for b in ("lightmotif", "moods", "biopython"):           # every engine, the same hits
    assert gm.scan_motifs_matrix(L, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False,
                                 backend=b, verbose=False).equals(M)

lib  = gb.load_motifs("motifs.jaspar")                   # a Library, a dict of sequences, a frame: same result
seqs = gb.read_fasta("genome.fa")
gm.scan_motifs_matrix(L.to_pandas(), seqs, lib, r=50, threshold=7.0, norm=False, verbose=False).equals(M)
# -> True

gm.scan_motifs_matrix(gb.as_loci([("chr1", 10, 60), ("chrZ", 100, 200)]), "genome.fa", "motifs.jaspar", r=50)
# -> [motifs] 2 of 2 windows run off a chromosome or are missing from the FASTA; their rows count 0.
#                      M1   M2
#    uid
#    chr1:10-60(.)    0.0  0.0
#    chrZ:100-200(.)  0.0  0.0
```

### `scan_motifs_matrix_masked`

```python
gm.scan_motifs_matrix_masked(loci, fasta, motifs, anchors, *, format="jaspar", r=250, window=10,
                             threshold=13.0, pvalue=None, anchor_threshold=None, norm=True,
                             both_strands=False, skip_anchors=True, seed=None, workers=None,
                             backend=None, verbose=True) -> pandas.DataFrame
```

Every match of a motif whose name contains one of `anchors`
(case-insensitive; names and descriptions are searched) gets centre ±
`window` bp replaced by random bases, then the library is scanned on the
masked windows — "which motifs set these regions apart *besides* the
anchor?" (e.g. mask CTCF). `anchor_threshold` overrides the cutoff used to
find the anchor matches; `seed` fixes the random bases. Anchor motifs are
left out of the result unless `skip_anchors=False`. No matching anchor raises
`ValueError`.

```python
gm.scan_motifs_matrix_masked(L, "genome.fa", "motifs.jaspar", ["M1"], r=50, threshold=7.0, norm=False, seed=0)
# -> [mask] 1 anchor PSSM(s): M1
#                       M2
#    uid
#    chr1:950-1050(.)    2
#    chr1:4950-5050(.)   3
#    chr1:7000-7100(.)   0
gm.scan_motifs_matrix_masked(L, "genome.fa", "motifs.jaspar", ["TFA"], r=50, threshold=7.0, norm=False,
                             skip_anchors=False, seed=0, verbose=False)        # matched by description
# ->                    M1  M2
#    uid
#    chr1:950-1050(.)    0   2
#    chr1:4950-5050(.)   0   3
#    chr1:7000-7100(.)   0   0
gm.scan_motifs_matrix_masked(L, "genome.fa", "motifs.jaspar", ["CTCF"], r=50, verbose=False)
# -> ValueError: No motifs matched anchors ['CTCF']
```

{: .note }
> When the queries were masked, mask the reference with the same anchors
> before `compare_motifs_to_ref` — otherwise every query looks depleted for
> reasons unrelated to biology.

### `scan_motifs_profile`

```python
gm.scan_motifs_profile(loci, fasta, motifs, select=None, *, format="jaspar", match="substring",
                       r=500, n_bins=100, threshold=13.0, pvalue=None, both_strands=True,
                       smooth=1.0, norm=False, workers=None, backend=None, verbose=True)
    -> (np.ndarray, list[str])
```

Each window (centre ± `r`) is split into `n_bins` bins and every hit is
credited to the bin of its centre; both strands by default (reverse matches
in forward coordinates, so they share the frame). `select` picks motifs by
name (case-insensitive substring unless `match='exact'`) — usually a handful;
`None` takes the whole library. `smooth` is a Gaussian sigma in bins along
each row (`0` keeps raw counts); `norm` divides by motif width. Returns the
`float32` cube `(rows, motifs, bins)` and the motif names of its middle axis —
the layout of a signal cube, so
[`gb.plot_motif_heatmap`]({{ '/api/signal/' | relative_url }}#plot_motif_heatmap)
(or `signal_draw.plot_heatmap`) draws it like one.

```python
prof, names = gm.scan_motifs_profile(L, "genome.fa", "motifs.jaspar", ["M1"], r=60, n_bins=12, threshold=7.0, smooth=0)
# -> [profile] 3 loci x 1 motifs x 12 bins | 120 bp window, 10 bp/bin: M1
#    [profile] mean hits/locus: M1=3.33
prof.shape, prof.dtype, names
# -> ((3, 1, 12), dtype('float32'), ['M1'])
prof[0, 0].tolist()                       # the three ACGT copies at the window centre
# -> [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0]
gm.scan_motifs_profile(L, "genome.fa", "motifs.jaspar", ["CTCF"], verbose=False)
# -> ValueError: No motifs matched ['CTCF'] (match='substring'); check the names exist in the library.
```

---

## Enrichment and differential statistics
{: .sec-navy }

These take the matrices `scan_motifs_matrix` / `scan_motifs_matrix_masked`
return (rows = sequences, columns = motifs). Missing columns are filled with
0, so the matrices need not share a motif set.

| Function | Returns | One line |
|---|---|---|
| `bootstrap_enrichment` | `DataFrame` | Resampled mean counts per group against a reference pool; LFC point estimates. |
| `compare_motifs` | `DataFrame` | Per-motif Mann-Whitney U between two matrices + LFC + BH FDR. |
| `compare_motifs_to_ref` | `DataFrame` | One or more query groups against a reference pool, with per-motif significance. |

### `bootstrap_enrichment`

```python
gm.bootstrap_enrichment(groups, ref, *, boot=100, sample=500, pseudo=0.1, seed=None, verbose=True)
    -> pandas.DataFrame
```

`groups` is `{name: matrix}`; `ref` the pool (e.g. all CREs). Each of `boot`
iterations draws `sample` rows without replacement (capped at the matrix
size) and averages them; the LFC is `log2((pseudo + mean_g) / (pseudo +
mean_ref))`. Columns: `Factor`, `mean_ref`, then `mean_<g>` and `LFC_<g>`
per group; with exactly two groups also `LFC = LFC_<first> - LFC_<second>`.

```python
Lb = gb.as_loci([("chr1", s, s + 100) for s in range(500, 6500, 200)])
Mb = gm.scan_motifs_matrix(Lb, "genome.fa", "motifs.jaspar", r=50, threshold=7.0, norm=False, verbose=False)
gm.bootstrap_enrichment({"a": Mb.iloc[:10], "b": Mb.iloc[10:]}, Mb, boot=20, sample=5, seed=0, verbose=False)
# ->   Factor  mean_ref  mean_a     LFC_a  mean_b     LFC_b       LFC
#    0     M1      0.41    0.48  0.185556    0.29 -0.387023  0.572579
#    1     M2      0.36    0.30 -0.201634    0.32 -0.131245 -0.070389
```

### `compare_motifs`

```python
gm.compare_motifs(mat_a, mat_b, *, pseudo=0.1, alternative="two-sided", verbose=True) -> pandas.DataFrame
```

Per motif: a Mann-Whitney U test on the per-sequence counts (rank-based, no
normality assumption — right for zero-inflated counts) and the log2 fold
change of the column means. `alternative` goes to
`scipy.stats.mannwhitneyu`. Columns `Factor`, `mean_A`, `mean_B`, `LFC`, `U`,
`p`, `p_adj` (Benjamini-Hochberg), sorted by `LFC` descending.

```python
gm.compare_motifs(Mb.iloc[:10], Mb.iloc[10:], pseudo=0.1)
# -> [compare_motifs] 2 motifs | sig (p_adj<0.01): 0
#      Factor  mean_A  mean_B       LFC      U         p     p_adj
#    0     M1     0.5    0.35  0.415037  113.5  0.477653  0.748428
#    1     M2     0.3    0.40 -0.321928   93.5  0.748428  0.748428
```

### `compare_motifs_to_ref`

```python
gm.compare_motifs_to_ref(query, ref, *, pseudo=0.1, alternative="two-sided", verbose=True) -> pandas.DataFrame
```

`query` is one matrix or a `{name: matrix}` dict. Columns: `Factor`,
`mean_ref`, and per group `g`: `mean_<g>`, `LFC_<g>`, `U_<g>`, `p_<g>`,
`p_adj_<g>` (BH within the group); with exactly two groups also
`LFC = LFC_<first> - LFC_<second>`. A single matrix is the group `query`.
Sorted by `LFC` (or the first group's LFC) descending.

```python
gm.compare_motifs_to_ref({"a": Mb.iloc[:10], "b": Mb.iloc[10:]}, Mb, verbose=False).columns.tolist()
# -> ['Factor', 'mean_ref', 'mean_a', 'LFC_a', 'U_a', 'p_a', 'p_adj_a', 'mean_b', 'LFC_b', 'U_b', 'p_b', 'p_adj_b', 'LFC']
gm.compare_motifs_to_ref(Mb.iloc[:10], Mb, verbose=False).columns.tolist()
# -> ['Factor', 'mean_ref', 'mean_query', 'LFC_query', 'U_query', 'p_query', 'p_adj_query']
```

---

## Archetypes
{: .sec-navy }

Pairwise Sandelin-Wasserman PWM similarity, hierarchical clustering, and a
consensus PFM per cluster — the recipe of STAMP, RSAT matrix-clustering and
Vierstra et al. 2020, in numpy / scipy. PFMs are `(4, W)` probability
matrices, rows A C G T. Logos are drawn by
[`genomeblocks.motifs_draw`]({{ '/api/signal/' | relative_url }}#genomeblocksmotifs_draw).

| Function | Returns | One line |
|---|---|---|
| `pwm_distance_matrix` | `(D, names, pwms)` | Pairwise distances `1 - SW similarity` of a library. |
| `cluster_motifs` | `(labels, Z)` | Hierarchical clustering of `D`. |
| `archetype` | `(4, W)` PFM | Consensus of a list of PFMs aligned to their medoid. |
| `archetype_from_names` | `dict` | One archetype from a named subset of a library, no clustering. |
| `build_archetypes` | `dict` | End to end: load, distances, cluster, one archetype per cluster. |
| `write_meme` | | `write_meme(pfms, path, *, alphabet="ACGT", bg=(0.25, 0.25, 0.25, 0.25))`: a `{name: pfm}` dict as a MEME file, so archetypes can be re-scanned. |

### `pwm_distance_matrix`

```python
gm.pwm_distance_matrix(motifs, *, format="jaspar", pseudo=0.01, min_overlap=5, workers=None, verbose=True)
    -> (np.ndarray, list[str], list[np.ndarray])
```

`motifs` is anything `load_motifs` takes; `pseudo` is added before counts
are turned into probabilities; `min_overlap` is the smallest alignment
considered. Pairs are spread over `workers` processes when there are at least
1000 of them.

```python
D, names, pwms = gm.pwm_distance_matrix("motifs.jaspar", min_overlap=3, verbose=False)
D.round(3), names, [p.shape for p in pwms]
# -> (array([[0.   , 0.661],
#            [0.661, 0.   ]]), ['M1', 'M2'], [(4, 4), (4, 4)])
```

### `cluster_motifs`

```python
gm.cluster_motifs(D, *, cutoff=0.3, linkage_method="average") -> (np.ndarray, np.ndarray)
```

`labels` are 1-based cluster ids (`scipy.cluster.hierarchy.fcluster` with
the `'distance'` criterion at `cutoff`); `Z` is the linkage matrix for
`plot_dendrogram`. A cutoff of `c` merges motifs with SW similarity of at
least `1 - c`; 0.25–0.35 is a typical range.

```python
labels, Z = gm.cluster_motifs(D, cutoff=0.5)
labels, Z.shape
# -> (array([1, 2], dtype=int32), (1, 4))
```

### `archetype`

```python
gm.archetype(pwms, *, min_overlap=5, weight_by_ic=True) -> np.ndarray
```

Picks the medoid (lowest mean distance to the rest), aligns every other
member to it (best offset and strand), and averages the aligned PFMs
column-wise in a frame that grows to fit members past either edge, weighting
each motif by its mean per-column information content when `weight_by_ic`.
Columns no member covers fall back to uniform.

```python
gm.archetype(pwms, min_overlap=3).round(2)
# -> array([[1. , 0. , 0. , 0.5, 1. ],
#           [0. , 0.5, 0. , 0. , 0. ],
#           [0. , 0.5, 1. , 0. , 0. ],
#           [0. , 0. , 0. , 0.5, 0. ]])
```

### `archetype_from_names`

```python
gm.archetype_from_names(motifs, names, *, format="jaspar", match="substring", pseudo=0.01,
                        min_overlap=5, weight_by_ic=True, verbose=True) -> dict
```

`names` are matched case-insensitively as substrings of the motif names and
descriptions (`match='substring'`) or verbatim (`'exact'`). Returns
`{'archetype': (4, W) PFM, 'members': [motif names], 'pwms': [their PFMs]}` —
the input of `plot_cluster_members`.

```python
one = gm.archetype_from_names("motifs.jaspar", ["TFA", "TFB"], min_overlap=3)
# -> [archetype_from_names] 2 motifs → 1 archetype: M1, M2
one["members"], one["archetype"].shape
# -> (['M1', 'M2'], (4, 5))
```

### `build_archetypes`

```python
gm.build_archetypes(motifs, *, format="jaspar", cutoff=0.3, min_overlap=5, linkage_method="average",
                    pseudo=0.01, weight_by_ic=True, workers=None, name_prefix="ARCH", verbose=True) -> dict
```

Returns a dict with `D`, `Z`, `labels`, `names`, `pwms`, `archetypes`
(`{ARCH_01: pfm, ...}`) and `members` (`{ARCH_01: [motif names], ...}`).

```python
res = gm.build_archetypes("motifs.jaspar", cutoff=0.5, min_overlap=3)
# -> [archetypes] 2 motifs → 2 archetypes (cutoff=0.5, linkage=average). Top sizes: [1, 1]
res["members"]
# -> {'ARCH_01': ['M1'], 'ARCH_02': ['M2']}
gm.write_meme(res["archetypes"], "archetypes.meme")
gm.scan_motifs_matrix(L, "genome.fa", "archetypes.meme", format="meme", r=50, threshold=5.0, norm=False, verbose=False)
# ->                    ARCH_01  ARCH_02
#    uid
#    chr1:950-1050(.)         3        1
#    chr1:4950-5050(.)        0        3
#    chr1:7000-7100(.)        2        0
```
