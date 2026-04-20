---
title: Examples
layout: default
nav_order: 7
---

# Examples
{: .no_toc }

End-to-end scripts you can copy and adapt. The files live under [`examples/`](https://github.com/birkiy/genomeblocks/tree/main/examples) in the repository.
{: .fs-5 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Browser view of the Nanog locus

**[`examples/browser_example.py`](https://github.com/birkiy/genomeblocks/blob/main/examples/browser_example.py)**

Renders a four-track IGV-like SVG around the mouse Nanog locus:

- HiChIP loops (half-sine arcs)
- ATAC bigWig coverage (binned)
- ATAC narrowPeak rectangles
- GENCODE vM25 protein-coding gene models

```python
from genomeblocks import browser, Loci, Genes
from genomeblocks.bedpe import read_bedpe

peaks = Loci.make("atac_peaks.narrowPeak")
loops = read_bedpe("nanog_hichip.bedpe")
genes = Genes.make("gencode.vM25.gtf")

fig, _ = browser(
    ("chr6", 122_286_666, 122_902_344),
    tracks={
        "HiChIP loops": loops,
        "ATAC signal":  "atac.bw",
        "ATAC peaks":   peaks,
        "Genes":        genes,
    },
    figsize=(12, None),
)
fig.savefig("nanog_locus.svg")
```

Output: [`examples/browser_Nanog.svg`](https://github.com/birkiy/genomeblocks/blob/main/examples/browser_Nanog.svg).

---

## 2. Cancer-gene regulatory network

**[`examples/real_world_example.py`](https://github.com/birkiy/genomeblocks/blob/main/examples/real_world_example.py)**

Two-hop spreading network starting from a list of cancer-gene promoters:

```python
from genomeblocks import Loci, Genes, Architecture

all_cres = Loci.make("cres.bed")
genes    = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)

cancer = {"MYC", "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA"}
prom   = Loci([g.tss for g in genes.values() if g.gene_name in cancer]).slop(1000)

arch = (Architecture.make_spread(prom, "hichip.bedpe",
                                 hops=2, r=2500, directed=True)
                     .add_mcool(all_cres, "hic.mcool", resolution=5000)
                     .normalize(all_cres)
                     .annotate(all_cres, genes))

# Find hubs and focus genes
result = arch.prime_hubs(key="n")
print(result["prime_genes"])
```

The script also exports edge and vertex tables, and draws the MYC sub-network.

---

## 3. ROSE-style enhancer-to-gene mapping

**[`examples/rose_region_to_gene.py`](https://github.com/birkiy/genomeblocks/blob/main/examples/rose_region_to_gene.py)**

ROSE produces super-enhancer calls but assigns genes by enhancer-center → TSS distance. `Genes.enhancer_to_genes` reproduces the overlap / proximal / closest categories with **real interval-to-interval distances** and alt-promoter resolution.

```python
df = genes.enhancer_to_genes(super_enhancers, prox=50_000, level="transcript")
df.to_csv("SE_region_to_gene.tsv", sep="\t", index=False)
```

Columns: `uid, overlap, proximal, closest`.

---

## 4. Compare CRE sets across conditions

```python
from genomeblocks import Loci, compare_heatmap

cre_a = Loci.make("cre_mESC.bed")
cre_b = Loci.make("cre_hESC.bed")

fig, union, S, tags = compare_heatmap(
    a=cre_a, b=cre_b,
    bigwigs=["ATAC_mESC.bw", "ATAC_hESC.bw",
             "H3K27ac_mESC.bw", "H3K27ac_hESC.bw"],
    a_name="mESC", b_name="hESC",
    samples={"ATAC": [0, 1], "H3K27ac": [2, 3]},
    normalize=True,
    cmap=["Blues", "Reds"],
    vmax=[10, 6],
)
fig.savefig("mESC_vs_hESC.pdf")
```

---

## 5. One-line pipelines

### CRE catalogue from ATAC

```python
cre = Loci.make("atac.narrowPeak").slop(100).sort().merge()
```

### Active-enhancer set via Tags

```python
tags = Tags.make(cre).add({"atac": atac, "h3k27ac": h3k27ac, "prom": promoters})
active_enhancers = tags.query(lambda l: (l.atac & l.h3k27ac) - l.prom)
```

### Nearest-gene annotation

```python
df = Genes.make("gencode.gtf").nearest_genes(cre)
```

### Motif scan across a CRE set

```python
counts = cre.scan_motifs(make_genome("hg38.fa"), "jaspar.jaspar", r=250)
```

---

## Longer walk-throughs

The repository also ships a few markdown guides alongside the scripts:

- [`MAKE_SPREAD_GUIDE.md`](https://github.com/birkiy/genomeblocks/blob/main/examples/MAKE_SPREAD_GUIDE.md) — deep dive into the `make_spread` algorithm.
- [`MAKE_SPREAD_NEW_API.md`](https://github.com/birkiy/genomeblocks/blob/main/examples/MAKE_SPREAD_NEW_API.md) — API reference for `make_spread`.
- [`ROSE_dff.md`](https://github.com/birkiy/genomeblocks/blob/main/examples/ROSE_dff.md) — notes on how `enhancer_to_genes` differs from classical ROSE.
