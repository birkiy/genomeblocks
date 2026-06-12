---
title: Quickstart
layout: default
nav_order: 3
---

# Quickstart
{: .no_toc }

A 10-minute tour of `genomeblocks`. We'll load ATAC peaks, annotate them with a GTF, build a chromatin-contact graph from HiChIP loops, compute O/E weights from a Hi-C matrix, and render a multi-track browser view.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Load peaks as `Loci`

```python
from genomeblocks import Loci

peaks = Loci.make("atac_peaks.narrowPeak")
print(peaks)                     # Loci(n=73412)
print(peaks[0])                  # Locus(chrom='chr1', start=9800, end=10100, strand='.')
print(peaks[0].uid)              # 'chr1:9800-10100(.)'
```

`Loci` is a `list` subclass, so indexing, iteration, and `len()` work as expected. It also exposes an on-demand `cgranges` index for O(log n) overlap queries.

---

## 2. Interval algebra

```python
# Extend each peak by ±200 bp, then collapse overlapping peaks.
cre = peaks.slop(200).sort().merge()

# Boolean set operations
superenhancers = Loci.make("H3K27ac_SE.bed")
se_cre = cre & superenhancers    # intersection
only_cre = cre - superenhancers  # difference
combined = cre | superenhancers  # union
xor = cre ^ superenhancers       # symmetric difference
```

See the [Loci guide](guide/loci) for the full API.

---

## 3. Annotate with a GTF

```python
from genomeblocks import Genes

genes = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)

# Per-CRE genomic context (Promoter-TSS / 5UTR / 3UTR / Exonic / Intronic / Intergenic)
annot_df = genes.annotations(cre)
print(annot_df.groupby("annotation").size())

# Nearest gene via pyranges
near_df = genes.nearest_genes(cre)
```

`Genes` also parses UCSC RefSeq tables (`Genes.make_ucsc(...)`) and supports alt-promoter-aware transcript-level operations. See [Genes guide](guide/genes).

---

## 4. Annotate and group CREs

Label CREs by gene context, and split them with set algebra. Grouping for
plots is just a plain `dict[str, Loci]` — no separate annotation object:

```python
prom = Loci(list(genes.annot['prom']))

annot  = genes.annotations(cre)                  # region class per CRE uid
counts = annot["annotation"].value_counts()

groups = {"promoters": cre & prom,               # CREs overlapping a promoter
          "enhancers": cre - prom}               # the rest
```

`groups` feeds straight into `plot_heatmap(..., groups=groups)` (see the
[AR & FOXA1 walkthrough](walkthrough/signal-heatmaps)).

---

## 5. Build a chromatin-contact graph

```python
from genomeblocks import Architecture

arch = (Architecture.make(cre, "HiChIP_loops.bedpe", r=2500)
                    .add_mcool(cre, "cohesin.mcool", resolution=5000)
                    .normalize(cre, source="w", name="n"))

print(arch)
# Architecture(name='Skeleton', loci=42311, links=185093, vertex_props=[uid], edge_props=[w, n, d])
```

Annotate the graph and extract hub genes in one line:

```python
arch.annotate(cre, genes)
result = arch.prime_hubs(key="n")
print(sorted(result["prime_genes"]))
```

Full API: [Architecture guide](guide/architecture).

---

## 6. Compare two Loci sets as a heatmap

```python
from genomeblocks import compare_heatmap

bigwigs = ["ATAC.bw", "H3K4me1.bw", "H3K27ac.bw"]

fig, union, S, tags = compare_heatmap(
    a=cre,
    b=superenhancers,
    bigwigs=bigwigs,
    a_name="CRE-only",
    b_name="SE-only",
    common_name="shared",
    n_bins=200,
    flank=3000,
    normalize=True,
    cmap=["Blues", "Reds", "Greens"],
    vmax=[10, 5, 8],
)
fig.savefig("compare.pdf")
```

See [Signal guide](guide/signal) for the full pipeline.

---

## 7. Render a browser view

```python
from genomeblocks import browser
from genomeblocks.bedpe import read_bedpe

fig, _ = browser(
    region=("chr6", 122_600_000, 122_800_000),
    tracks={
        "HiChIP loops": read_bedpe("loops.bedpe"),
        "ATAC signal":  "ATAC.bw",
        "ATAC peaks":   peaks,
        "Genes":        genes,
    },
    figsize=(12, None),      # auto-height
    bw_n_bins=120,
)
fig.savefig("browser.svg")
```

A full Nanog-locus example is in [`examples/browser_example.py`](https://github.com/birkiy/genomeblocks/blob/main/examples/browser_example.py) and renders into an IGV-like, fully-vectorial SVG.

---

## Next steps

- [Concepts →](concepts) — the mental model behind the API.
- [Example: AR & FOXA1 →](walkthrough/) — a complete real-data walkthrough, concept by concept.
- [Architecture guide →](guide/architecture) — build & mine chromatin networks.
