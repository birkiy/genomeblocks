---
title: "Example: AR & FOXA1"
layout: default
nav_order: 7
has_children: true
permalink: /walkthrough/
---

# Example: AR & FOXA1 in LNCaP
{: .no_toc }

A complete, real-data walkthrough of `genomeblocks` — from raw peak calls to a
publication-style figure set — built around one biological question.
{: .fs-5 .fw-300 }

The runnable notebook lives at
[`examples/ar_foxa1_lncap/`](https://github.com/birkiy/genomeblocks/tree/main/examples/ar_foxa1_lncap).
These pages explain the *concepts* behind each step: the biology, the method,
and the `genomeblocks` API that implements it.

## The question

In prostate cells, the **androgen receptor (AR)** is a hormone-activated
transcription factor: androgen (here **DHT**, dihydrotestosterone) triggers AR
to enter the nucleus and bind chromatin. But AR rarely works alone — it leans on
**FOXA1**, a *pioneer factor* that can engage closed chromatin and pre-license
sites for other factors.

We use **LNCaP** cells at **0 h** and **4 h** of DHT and ask:

> Of the AR sites that appear after androgen stimulation, which are
> **FOXA1-dependent** and which are **FOXA1-independent** — and how do those two
> classes differ in accessibility, signal, genomic location, sequence motifs,
> and co-factor binding?

We split the 4 h AR cistrome into two sets and characterise each:

| set | definition | interpretation |
|---|---|---|
| **AR+F** | AR 4 h peaks that overlap a FOXA1 peak (0 h or 4 h) | FOXA1-dependent AR |
| **AR−F** | AR 4 h peaks that overlap **no** FOXA1 peak | FOXA1-independent AR |

## The pipeline

Each step maps to one page in this section:

1. **[Peaks, set algebra & categories](peaks-and-sets)** — load peaks as `Loci`,
   build accessible chromatin from ATAC, and derive **AR+F** / **AR−F** with set
   operations and a Venn diagram.
2. **[Signal heatmaps](signal-heatmaps)** — extract bigWig signal into a cube,
   average ATAC replicates, and draw grouped heatmaps.
3. **[Genomic annotation](annotation)** — label each set by gene context
   (promoter / intron / intergenic …) with a `Genes` model.
4. **[Motif & ChIP-Atlas enrichment](enrichment)** — what sequence motifs and
   what published TF datasets are enriched in each set.
5. **[Genome browser](browser)** — an IGV-like view of all conditions at a
   locus, with replicate averaging and shared y-axes.

## Data

ChIP-Atlas (hg38) peak + bigWig files, fetched by the example's
`download_data.sh`:

| factor | condition | accession |
|---|---|---|
| FOXA1 | 0 h / 4 h DHT | `SRX23002839` / `SRX23002841` |
| AR | 0 h / 4 h DHT | `SRX23002834` / `SRX23002836` |
| ATAC | 0 h DHT (rep1/rep2) | `SRX23002894` / `SRX23002895` |
| ATAC | 4 h DHT (rep1/rep2) | `SRX23002898` / `SRX23002899` |

Each accession has a **peak** file (`.05.bed`, q < 1e‑5) and a **signal** file
(`.bw`). Peaks tell you *where* a factor binds; bigWigs give the *continuous
coverage* you average for heatmaps and browser tracks.

Annotation, motif scanning, and the ChIP-Atlas index also need genome-scale
references (hg38 FASTA, a GTF, a JASPAR motif database, and a GIGGLE index);
their paths are set once at the top of the notebook.

{: .note }
> Every long step in the notebook is wrapped in a small `timer()` context manager
> and loops use `tqdm`, so you can see where the time goes on real data.

Start with **[Peaks, set algebra & categories →](peaks-and-sets)**
