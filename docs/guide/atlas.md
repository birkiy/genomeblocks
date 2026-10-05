---
title: Atlas
parent: User Guide
layout: default
nav_order: 8
---

# Atlas
{: .no_toc }

GIGGLE-style enrichment of a region set against a large collection of BED
tracks (all of ChIP-Atlas, say). One sparse `bin x track` index turns each
query into a single sparse matrix-vector product, so thousands of tracks are
scored in one vectorised pass.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## The idea
{: .sec-navy }

Tile the genome at a fixed resolution (default 1 kb) and store one bit per
`(bin, track)` cell in a single CSR matrix. A query set is reduced to the bins
it covers; `M[query_bins].sum(0)` then returns the shared-bin count for
**every track at once**. From those counts the `Atlas` computes a Fisher 2x2
enrichment (the GIGGLE score) or an empirical Monte-Carlo null.

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/atlas-index.svg %}
</div><figcaption>
<strong>One sparse product scores every track.</strong> The index holds one bit per (bin, track). A query becomes the bins it covers; <code>M[query_bins].sum(axis=0)</code> gives, for every track at once, the shared-bin count <em>a</em>, and the track's size and the query's size fill the other three cells of the Fisher table. <code>ref=</code> swaps the genome for a second set in those cells; <code>bootstrap()</code> reshuffles the query positions instead.
</figcaption></figure>

At ChIP-Atlas scale (~25k tracks, 1 kb bins) the index is a few GB resident
and each query is sub-second.

> Layer et al., *GIGGLE: a search engine for large-scale integrated genome
> analysis.* Nat. Methods 15, 123–126 (2018).

Every query, reference or pool below is anything
[`as_loci`]({{ '/interoperability/' | relative_url }}) takes: a `Loci`, a BED
path, a pandas / polars frame, a list of region strings.

```python
import genomeblocks as gb
from genomeblocks import Atlas
```

---

## Building an index — `Atlas.make`
{: .sec-navy }

```python
atlas = Atlas.make(
    "chipatlas/*.bed.gz",               # glob, directory, one path, or a list of paths
    chromsizes="hg38.chrom.sizes",     # path, {chrom: length}, Genome, pandas Series, or a cooler
    bin_size=1000,
    names=None,                        # track names; default: the file name without its BED extension
    name_pattern=r"^[^.]+",            # or a regex over the file name: keep just the SRX accession
    workers=None,                      # None = half the cores, one BED per task
    meta="chipatlas_meta.tsv",         # optional per-track table (see below)
    meta_id_col="srx",
)
```

`paths` is a glob string, a directory (every `.bed` / `.narrowPeak` /
`.broadPeak`, gzipped or not), a single file, or an explicit list.
`chromsizes` is a UCSC `.chrom.sizes` path, a `{chrom: length}` dict, a
`gb.Genome` with sizes (`Genome.from_fasta("hg38.fa")`), a pandas Series, or
anything with a `.chromsizes` mapping such as a `cooler.Cooler`. Intervals on
chromosomes not in `chromsizes` are dropped.

Track names default to the file name with its extension stripped
(`SRX23002840.05.bed.gz` → `SRX23002840.05`). `name_pattern` is applied with
`re.search` to the file name and keeps the first group (or the whole match):
`r"^[^.]+"` or `r"([^.]+)"` gives `SRX23002840`, which is what joins against a
metadata table. `names=` overrides them one by one.

On three small tracks (the fixture of the test suite):

```python
atlas = Atlas.make(["t0.bed", "t1.bed", "t2.bed"], chromsizes={"chr1": 20_000, "chr2": 8_000},
                   bin_size=500, verbose=False)
atlas
# -> Atlas(tracks=3, bins=56, bin_size=500, nnz=21)
atlas.track_names, atlas.n_bins, atlas.chrom_offsets
# -> (['t0', 't1', 't2'], 56, {'chr1': 0, 'chr2': 40})
```

A glob with no match raises `ValueError: No BED files matched ...`; a list with
a missing file names the first missing path.

### Save and load

Persist the index as a compressed `.npz` — plain arrays, never a pickle, so
it is safe to share and loads without the code that built it:

```python
atlas.save("chipatlas_hg38_1kb.npz")
atlas = Atlas.load("chipatlas_hg38_1kb.npz")       # tracks, bins, matrix and metadata come back
```

---

## The track table
{: .sec-navy }

An `Atlas` is a table of its tracks — one row per track with `name`,
`n_peaks` (intervals read), `n_bins` (bins covered) and any metadata columns:

```python
atlas.shape, atlas.columns
# -> ((3, 3), ['name', 'n_peaks', 'n_bins'])
atlas.head(2)
# ->   name  n_peaks  n_bins
# -> 0   t0        7       7
# -> 1   t1        7       7
atlas.describe()
# ->                        value
# -> tracks                     3
# -> chromosomes                2
# -> bin size                 500
# -> bins                      56
# -> nonzero cells             21
# -> density                0.125
# -> bins per track min         7
# -> bins per track median    7.0
# -> bins per track max         7
# -> metadata columns           —
atlas["t1"]                       # one track's row as a dict (by name or position)
# -> {'name': 't1', 'n_peaks': 7, 'n_bins': 7}
list(atlas)                       # -> ['t0', 't1', 't2']
```

`to_pandas()`, `to_polars()`, `to_arrow()` and the Arrow / dataframe-interchange
protocols are there too, so polars, DuckDB and seaborn take the atlas as it is:

```python
import polars as pl, duckdb, seaborn as sns

pl.DataFrame(atlas)                                                      # shape: (3, 3)
duckdb.sql("select name, n_peaks from atlas order by n_peaks desc").fetchall()
sns.barplot(atlas, x="name", y="n_bins")
```

The matrix itself is `atlas.M` (scipy CSR, `bins x tracks`, uint8).

---

## Enrichment search — `atlas.search(query)`
{: .sec-navy }

Fisher 2x2 of each track against the genome null, in bin units:

```python
cre = gb.Loci.make("my_peaks.narrowPeak")
res = atlas.search(cre)                        # DataFrame sorted by giggle_score
res.head()[["name", "overlaps", "log2_odds", "p", "giggle_score"]]
```

| column | meaning |
|---|---|
| `name` | track |
| `n_query_bins` | bins the query covers |
| `overlaps` | query bins the track also covers (cell *a* of the 2x2) |
| `track_n_bins`, `track_n_peaks` | the track's own size |
| `log2_odds` | log2 odds ratio (0.5 added to every cell) |
| `p` | exact hypergeometric p (`alternative='two-sided'`, `'greater'` or `'less'`) |
| `giggle_score` | `-log10(p) · log2(OR)` — negative for depletion; the sort key |

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/atlas-search-table.svg %}
</div><figcaption>
<strong>Three numbers fill the four cells, for every track at once.</strong> <code>a</code> is the number of query bins the track also covers, read off <code>M[query_bins].sum(0)</code>; <code>b = n_query_bins − a</code>, <code>c = track_n_bins − a</code> and <code>d = n_bins − a − b − c</code> follow from the query's bins, the track's bins and the genome's bins, as vectors over all tracks. From them come <code>log2_odds</code> (0.5 added to every cell), the exact hypergeometric <code>p</code> and <code>giggle_score = −log10(p) · log2_odds</code>, the sort key. With <code>ref=</code> the reference set takes the genome's place in <code>c</code> and <code>d</code>: <code>c</code> becomes the reference bins the track covers and <code>d = n_ref_bins − c</code>, so the test asks whether the track is more enriched in the query than in <code>ref</code>.
</figcaption></figure>

Metadata columns are merged in when attached. The query is anything
`as_loci` takes, so these all work and give the same table:

```python
atlas.search(cre)
atlas.search(cre.to_pandas())
atlas.search(cre.to_polars())
atlas.search("my_peaks.narrowPeak")
atlas.search(["chr1:900-1100", "chr1:4900-5100"])
```

Pass a **reference set** to ask "more enriched in the query than in the
reference?" — the *c* / *d* cells of the 2x2 come from `ref` instead of the
genome, and `n_ref_bins` replaces `n_track_bins`:

```python
res = atlas.search(query=up_cre, ref=all_cre)
```

A query with no bins on the atlas's chromosomes raises
`ValueError: Query has no bins on the atlas's chromosomes.` — usually a
chromosome naming mismatch (`1` vs `chr1`).

The same call is a method on every `Loci`:

```python
res = cre.enrich(atlas)                # == atlas.search(cre)
res = up.enrich(atlas, ref=all_cre)
```

---

## Monte-Carlo null — `atlas.bootstrap(query, n=...)`
{: .sec-navy }

When you want an empirical null instead of the analytic Fisher test:

```python
res = atlas.bootstrap(cre, n=1000, seed=0)    # position-shuffle null, chromosome kept
res[["name", "observed", "expected", "z", "p_emp"]].head()
```

One row per track: `observed` and `expected` overlap counts (means over the
iterations), their `obs_std` / `null_std`, `log2fc = log2((obs + 1) / (exp + 1))`,
the paired `z` of `observed − null` per iteration (the sort key) and
`p_emp = (1 + #iterations with obs ≤ null) / (n + 1)`.

Three interacting knobs control the null:

| argument | null model |
|---|---|
| `pool=None` (default) | per-interval position shuffle; `keep_chrom=True` keeps each interval on its chromosome, `False` draws a chromosome in proportion to its length |
| `pool=<intervals>` | draw the null from a curated universe (LOLA / regioneR style) — controls for the universe's own bias; each group draws its own size from the pool |
| `sample=<int>` | per iteration, subsample every query group **and** the pool to `sample` regions; the pool draw is shared across groups so comparisons are paired |

`replace=True` allows drawing with replacement when a group or the pool is
smaller than `sample` (otherwise that raises with the sizes).

### One query or a dict of groups

`query` may be a single interval set **or** a `{group: intervals}` dict. With
a dict the result is long-format with a `group` column, and all groups share
one null:

```python
res = atlas.bootstrap({"up": up_cre, "down": down_cre}, n=1000, seed=0)
res[res.group == "up"].head()

res = atlas.bootstrap({"up": up_cre, "down": down_cre}, pool=all_cre, n=1000, sample=500, seed=0)
```

{: .note }
> A dict whose keys are column names (`chrom` / `chr` / `start` / `end` …) is
> one query made of columns, not a dict of groups:
> `atlas.bootstrap({"chrom": ["chr1"], "start": [100], "end": [900]}, n=100)`
> returns a single-set result without a `group` column.

Fluent form: `cre.enrich_mc(atlas, n=1000, seed=0)`.

---

## Attaching metadata — `attach_meta`
{: .sec-navy }

Track names alone are rarely enough — you want antigen / cell-type labels
alongside the scores. Attach a table keyed on the track name:

```python
atlas.attach_meta(
    "chipatlas_meta.tsv",
    id_col="srx",                          # the column that joins against track_names (default: the first)
    columns=["srx", "antigen", "cell"],    # supply headers for a header-less file
    sep="\t",
)
atlas.attach_meta(meta_df, id_col="id")    # or an in-memory DataFrame
```

Tracks with no metadata row get `NaN`; metadata rows for unknown tracks are
dropped. The columns join the track table and every `search` / `bootstrap`
result:

```python
import pandas as pd
meta = pd.DataFrame({"id": ["t0", "t1", "t9"], "antigen": ["AR", "FOXA1", "X"],
                     "cell_line": ["LNCAP", "LNCAP", "VCaP"]})
atlas.attach_meta(meta, id_col="id")
atlas.columns
# -> ['name', 'n_peaks', 'n_bins', 'antigen', 'cell_line']
atlas.meta
# ->    antigen cell_line
# -> id
# -> t0      AR     LNCAP
# -> t1   FOXA1     LNCAP
# -> t2     NaN       NaN
```

{: .warning }
> Two arguments decide whether the metadata lines up. **`name_pattern`** (or
> `names`) must turn `SRX23002840.05.bed.gz` into the id used in the table
> (`SRX23002840`); otherwise every metadata column comes back `NaN`. And a
> **header-less** table needs `columns=` (`meta_columns=` on `Atlas.make`), or
> its first row is taken as the header.

`Atlas.make(..., meta=, meta_id_col=, meta_columns=, meta_sep=)` does the same
at build time, and `save` / `load` carry the metadata along.

---

## See also
{: .sec-navy }

- [API → Atlas]({{ '/api/atlas/' | relative_url }}) for full signatures.
- [Design → Atlas]({{ '/design/atlas/' | relative_url }}) for the sparse index.
- [Enrichment walkthrough]({{ '/walkthrough/enrichment/' | relative_url }}) for a worked ChIP-Atlas example.
