---
title: "Full run & figures"
parent: "Example: AR & FOXA1"
layout: default
nav_order: 6
---

{: .note }
> This is the AR & FOXA1 notebook executed end to end, with its figures and
> result tables. The previous pages explain the *concepts*; this one shows the
> actual run. Source:
> [`examples/ar_foxa1_lncap/`](https://github.com/birkiy/genomeblocks/tree/main/examples/ar_foxa1_lncap).
> The outputs shown are from the 1.x run; the code cells are the 2.0 API and
> the notebook has not yet been re-run on 2.0.

# AR & FOXA1 in LNCaP (±DHT)

How the **androgen receptor (AR)** cistrome depends on the pioneer factor
**FOXA1** after androgen (DHT) stimulation in LNCaP cells.

Data are ChIP-Atlas (hg38) peaks + bigwigs — run `./download_data.sh` first.
LNCaP, 0h vs 4h DHT, for FOXA1, AR, and ATAC-seq (two ATAC replicates / condition).

> **Note** — FOXA1 4h peaks use `SRX23002841.05.bed` (matching the FOXA1 4h
> bigwig). The reference paths in the next cell point at the lackgrp cluster;
> edit them for your environment.

## 0. Setup — paths & helpers


```python
import os, time
from contextlib import contextmanager

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import genomeblocks as gb
from genomeblocks import Loci, Genes, Atlas, tmm, browser
from genomeblocks.signal_draw import plot_heatmap
from genomeblocks.motifs import scan_motifs_matrix, bootstrap_enrichment


@contextmanager
def timer(label):
    """Print a step's wall-clock time."""
    t0 = time.perf_counter()
    print(f"\u25b6 {label} ...")
    try:
        yield
    finally:
        print(f"\u2713 {label} \u2014 {time.perf_counter() - t0:.1f}s")


# --- downloaded data (see download_data.sh) ---
DATA = "data"
DATA
SRX = {
    "FOXA1_0h": "SRX23002839", "FOXA1_4h": "SRX23002841",
    "AR_0h":    "SRX23002834", "AR_4h":    "SRX23002836",
    "ATAC_0h_r1": "SRX23002894", "ATAC_0h_r2": "SRX23002895",
    "ATAC_4h_r1": "SRX23002898", "ATAC_4h_r2": "SRX23002899",
}
PEAK = {k: f"{DATA}/{s}.05.bed" for k, s in SRX.items()}
BW   = {k: f"{DATA}/{s}.bw"     for k, s in SRX.items()}

# --- hg38 references (edit for your environment) ---
GENOME_FA   = "/groups/lackgrp/genome_annotations/hg38/hg38.fa"
GTF         = "/groups/lackgrp/genome_annotations/hg38/gencode.v49.annotation_protein_coding.gtf"
MOTIF_DB    = "/groups/lackgrp/databases/motifs/motif-db/H14CORE_jaspar_format_lightmotif.txt"
CHROMSIZES  = "/groups/lackgrp/genome_annotations/hg38/hg38_chr.chrom.sizes"
GIGGLE_META = "/groups/lackgrp/databases/giggle/giggle_hg38/hg38_tfs_meta.tsv"
GIGGLE_DIR  = "/groups/lackgrp/databases/giggle/giggle_hg38/ChIP-Atlas-ALL"
```

## 1. Load peaks


```python
with timer("load peak files"):
    peaks = {k: Loci.make(v) for k, v in PEAK.items()}

for k, v in peaks.items():
    print(f"{k:12} {len(v):>8,} peaks")
```

    ▶ load peak files ...
    ✓ load peak files — 0.6s
    FOXA1_0h        8,657 peaks
    FOXA1_4h       10,870 peaks
    AR_0h             910 peaks
    AR_4h           3,645 peaks
    ATAC_0h_r1     53,671 peaks
    ATAC_0h_r2     62,741 peaks
    ATAC_4h_r1     47,780 peaks
    ATAC_4h_r2     51,196 peaks


## 2. Accessible chromatin — union of ATAC peaks

Concatenate the four ATAC peak sets (both conditions, both reps) and `merge()`
overlapping intervals into one accessible-region set (the result comes back in
genome order).


```python
with timer("ATAC union (accessible regions)"):
    accessible = (peaks["ATAC_0h_r1"] + peaks["ATAC_0h_r2"]
                  + peaks["ATAC_4h_r1"] + peaks["ATAC_4h_r2"]).merge()

print(f"accessible regions (merged): {len(accessible):,}")
```

    ▶ ATAC union (accessible regions) ...
    ✓ ATAC union (accessible regions) — 0.7s
    accessible regions (merged): 72,470


## 3. AR / FOXA1 binding inside vs outside accessible chromatin


```python
rows = []
for name in ["AR_0h", "AR_4h", "FOXA1_0h", "FOXA1_4h"]:
    s = peaks[name]
    inside = len(s & accessible)              # peaks overlapping accessible
    total = len(s)
    rows.append((name, inside, total - inside))
    print(f"{name:10} {100 * inside / total:5.1f}% inside accessible ({inside:,}/{total:,})")

names = [r[0] for r in rows]
ins   = np.array([r[1] for r in rows])
outs  = np.array([r[2] for r in rows])
fig, ax = plt.subplots(figsize=(5, 3))
ax.bar(names, ins,  label="inside accessible", color="#4c78a8")
ax.bar(names, outs, bottom=ins, label="outside", color="#d6d6d6")
ax.set_ylabel("peaks"); ax.legend(frameon=False, fontsize=8)
ax.set_title("TF binding vs accessibility")
plt.xticks(rotation=30, ha="right"); plt.tight_layout()
```

    AR_0h       31.2% inside accessible (284/910)
    AR_4h       83.0% inside accessible (3,027/3,645)
    FOXA1_0h    85.0% inside accessible (7,362/8,657)
    FOXA1_4h    88.0% inside accessible (9,561/10,870)



    
![figure]({{ '/assets/images/ar_foxa1/ar_foxa1_lncap_8_1.png' | relative_url }})
    


## 4. Keep only accessible TF peaks


```python
with timer("filter TF peaks to accessible"):
    F0 = peaks["FOXA1_0h"] & accessible
    F4 = peaks["FOXA1_4h"] & accessible
    A4 = peaks["AR_4h"]    & accessible

print(f"accessible  FOXA1 0h={len(F0):,}  FOXA1 4h={len(F4):,}  AR 4h={len(A4):,}")
```

    ▶ filter TF peaks to accessible ...
    ✓ filter TF peaks to accessible — 0.0s
    accessible  FOXA1 0h=7,362  FOXA1 4h=9,561  AR 4h=3,027


## 5. Venn — FOXA1 (0h, 4h) vs AR (4h)

To venn genomic intervals we merge all peaks into a shared region *universe*,
then label each region by which input set overlaps it (`overlap_any` gives one
boolean per universe row; the row numbers are the region ids). From this:

- **AR+F** = AR 4h peaks that overlap FOXA1 (0h or 4h) — FOXA1-dependent AR
- **AR−F** = AR 4h peaks that overlap no FOXA1 peak — FOXA1-independent AR


```python
from matplotlib_venn import venn3

def venn_id_sets(sets):
    """Merge all peaks into a region universe; return one set of region-ids per
    input set (the ids it overlaps) so matplotlib_venn can count the regions."""
    universe = sum(sets[1:], sets[0]).merge()
    return [set(np.flatnonzero(universe.overlap_any(s)).tolist()) for s in sets]

with timer("venn3 membership"):
    ids = venn_id_sets([F0, F4, A4])

fig, ax = plt.subplots(figsize=(5, 5))
venn3(ids, set_labels=["FOXA1 0h", "FOXA1 4h", "AR 4h"], ax=ax)
ax.set_title("Accessible peaks: FOXA1 (0h, 4h) vs AR (4h)")
```

    ▶ venn3 membership ...


    venn membership: 100%|██████████| 11624/11624 [00:00<00:00, 238333.00it/s]

    ✓ venn3 membership — 0.1s


    





    Text(0.5, 1.0, 'Accessible peaks: FOXA1 (0h, 4h) vs AR (4h)')




    
![figure]({{ '/assets/images/ar_foxa1/ar_foxa1_lncap_12_5.png' | relative_url }})
    



```python
with timer("define AR+F / AR-F"):
    F_any = (F0 + F4).merge()            # FOXA1-bound at either timepoint
    ARpF  = A4 & F_any                   # AR 4h that IS FOXA1-bound
    ARmF  = A4 - F_any                   # AR 4h that is NOT FOXA1-bound

print(f"AR+F (FOXA1-dependent): {len(ARpF):,}   AR-F (FOXA1-independent): {len(ARmF):,}")
```

    ▶ define AR+F / AR-F ...
    ✓ define AR+F / AR-F — 0.0s
    AR+F (FOXA1-dependent): 2,515   AR-F (FOXA1-independent): 512


## 6. Signal heatmaps across conditions

Extract signal for AR+F and AR−F over all 8 bigwigs, average the two ATAC
replicates per condition (the same averaging the browser does), and draw a
grouped heatmap: ATAC 0h, ATAC 4h, AR 0h, AR 4h, FOXA1 0h, FOXA1 4h.


```python
regions = ARpF + ARmF
groups  = {"AR+F": ARpF, "AR-F": ARmF}

bw_order = ["ATAC_0h_r1", "ATAC_0h_r2", "ATAC_4h_r1", "ATAC_4h_r2",
            "AR_0h", "AR_4h", "FOXA1_0h", "FOXA1_4h"]
bigwigs = [BW[k] for k in bw_order]

with timer("extract signal cube"):
    S = np.nan_to_num(regions.signal(bigwigs, n_bins=200, flank=2000, workers=4))
#S = tmm(np.nan_to_num(S))                       # per-track normalization

# group ATAC replicates by averaging their columns -> 6 conditions
S6 = np.stack([
    S[:, [0, 1], :].mean(1),   # ATAC 0h (rep1+rep2)
    S[:, [2, 3], :].mean(1),   # ATAC 4h (rep1+rep2)
    S[:, 4, :], S[:, 5, :],    # AR 0h, AR 4h
    S[:, 6, :], S[:, 7, :],    # FOXA1 0h, FOXA1 4h
], axis=1)
samples = ["ATAC 0h", "ATAC 4h", "AR 0h", "AR 4h", "FOXA1 0h", "FOXA1 4h"]
vmax = float(np.percentile(S6, 99))

with timer("draw heatmap"):
    fig = plot_heatmap(regions, S6, groups=groups, sets=["AR+F", "AR-F"],
                       samples=samples, vmax=vmax, ymax=vmax)
```

    ▶ extract signal cube ...
    [INFO] Extracting 8 bigwigs for 3027 loci into 200 bins (span=False, agg='mean', backend='pybigtools', exact=True, workers=4).


    chunks: 100%|██████████| 4/4 [00:03<00:00,  1.16it/s]


    ✓ extract signal cube — 3.6s
    ▶ draw heatmap ...
    ✓ draw heatmap — 0.1s



    
![figure]({{ '/assets/images/ar_foxa1/ar_foxa1_lncap_15_3.png' | relative_url }})
    


## 7. Genomic annotation of AR+F vs AR−F


```python
with timer("load genes (GTF)"):
    genes = Genes.make(GTF)

with timer("annotate AR+F / AR-F"):
    annot_pf = genes.annotations(ARpF)
    annot_mf = genes.annotations(ARmF)

fig, axes = plt.subplots(1, 2, figsize=(9, 4))
for ax, df, title in [(axes[0], annot_pf, "AR+F"), (axes[1], annot_mf, "AR-F")]:
    vc = df["annotation"].value_counts()
    ax.pie(vc.values, labels=vc.index, autopct="%1.0f%%", textprops={"fontsize": 7})
    ax.set_title(f"{title}  (n={len(df):,})")
plt.tight_layout()
```

    ▶ load genes (GTF) ...


    [INFO] Parsing GTF/GFF file 🧩: 6731674it [01:05, 103025.37it/s]


    [INFO] Unmapped feature types: start_codon, Selenocysteine, stop_codon
    ✓ load genes (GTF) — 65.5s
    ▶ annotate AR+F / AR-F ...
    ✓ annotate AR+F / AR-F — 16.5s



    
![figure]({{ '/assets/images/ar_foxa1/ar_foxa1_lncap_17_3.png' | relative_url }})
    


## 8. Motif enrichment — AR+F vs AR−F

Scan JASPAR motifs over AR+F, AR−F, and a subsample of the accessible (ATAC)
pool as background, then bootstrap the log-fold-change. Sorting by `LFC`
(= LFC_AR+F − LFC_AR−F) gives motifs preferential to each set. The scanners
read the windows from the FASTA through its index; to hold the genome in
memory instead, pass `gb.read_fasta(GENOME_FA)` in place of the path.


```python
with timer("index genome FASTA"):
    genome = gb.Genome.from_fasta(GENOME_FA)     # builds hg38.fa.fai once; names and sizes

# background pool = subsample of accessible chromatin (keeps scanning cheap)
pool = accessible

with timer("scan motifs (AR+F, AR-F, pool)"):
    mat_pf   = scan_motifs_matrix(ARpF, GENOME_FA, MOTIF_DB, r=250, workers=8)
    mat_mf   = scan_motifs_matrix(ARmF, GENOME_FA, MOTIF_DB, r=250, workers=8)
    mat_pool = scan_motifs_matrix(pool, GENOME_FA, MOTIF_DB, r=250, workers=8)

with timer("bootstrap enrichment"):
    enr = bootstrap_enrichment({"AR+F": mat_pf, "AR-F": mat_mf},
                               ref=mat_pool, boot=100, sample=500, seed=0)

show = ["Factor", "LFC", "LFC_AR+F", "LFC_AR-F"]
ranked = enr.sort_values("LFC", ascending=False)
print("Top motifs in AR+F (FOXA1-dependent):")
display(ranked.head(10)[show])
print("Top motifs in AR-F (FOXA1-independent):")
display(ranked.tail(10)[show].iloc[::-1])
```

    ▶ index genome FASTA ...
    ✓ index genome FASTA — 20.4s
    ▶ scan motifs (AR+F, AR-F, pool) ...


    [motifs]: 100%|██████████| 33/33 [00:01<00:00, 29.99it/s]
    [motifs]: 100%|██████████| 33/33 [00:00<00:00, 118.08it/s]
    [motifs]: 100%|██████████| 33/33 [00:30<00:00,  1.08it/s]


    ✓ scan motifs (AR+F, AR-F, pool) — 55.5s
    ▶ bootstrap enrichment ...


    100%|██████████| 100/100 [00:00<00:00, 121.79it/s]
    100%|██████████| 100/100 [00:00<00:00, 140.44it/s]
    100%|██████████| 100/100 [00:00<00:00, 261.65it/s]


    ✓ bootstrap enrichment — 2.4s
    Top motifs in AR+F (FOXA1-dependent):



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Factor</th>
      <th>LFC</th>
      <th>LFC_AR+F</th>
      <th>LFC_AR-F</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1439</th>
      <td>ZN671.H14CORE.0.P.C</td>
      <td>0.200354</td>
      <td>0.189548</td>
      <td>-0.010806</td>
    </tr>
    <tr>
      <th>1098</th>
      <td>TSH2.H14CORE.0.SG.A</td>
      <td>0.149624</td>
      <td>0.145360</td>
      <td>-0.004264</td>
    </tr>
    <tr>
      <th>240</th>
      <td>FOXA2.H14CORE.0.PSM.A</td>
      <td>0.123507</td>
      <td>0.073603</td>
      <td>-0.049904</td>
    </tr>
    <tr>
      <th>242</th>
      <td>FOXA3.H14CORE.0.PS.A</td>
      <td>0.106075</td>
      <td>0.077160</td>
      <td>-0.028915</td>
    </tr>
    <tr>
      <th>265</th>
      <td>FOXL2.H14CORE.0.PSM.A</td>
      <td>0.094850</td>
      <td>0.072495</td>
      <td>-0.022355</td>
    </tr>
    <tr>
      <th>266</th>
      <td>FOXM1.H14CORE.0.P.B</td>
      <td>0.091278</td>
      <td>0.073294</td>
      <td>-0.017984</td>
    </tr>
    <tr>
      <th>271</th>
      <td>FOXP1.H14CORE.0.PS.A</td>
      <td>0.088453</td>
      <td>0.072825</td>
      <td>-0.015627</td>
    </tr>
    <tr>
      <th>262</th>
      <td>FOXK1.H14CORE.0.PS.A</td>
      <td>0.087973</td>
      <td>0.059407</td>
      <td>-0.028566</td>
    </tr>
    <tr>
      <th>238</th>
      <td>FOXA1.H14CORE.0.P.B</td>
      <td>0.079758</td>
      <td>0.051366</td>
      <td>-0.028392</td>
    </tr>
    <tr>
      <th>268</th>
      <td>FOXO3.H14CORE.0.PS.A</td>
      <td>0.068372</td>
      <td>0.043914</td>
      <td>-0.024458</td>
    </tr>
  </tbody>
</table>
</div>


    Top motifs in AR-F (FOXA1-independent):



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Factor</th>
      <th>LFC</th>
      <th>LFC_AR+F</th>
      <th>LFC_AR-F</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7</th>
      <td>ANDR.H14CORE.0.P.B</td>
      <td>-0.153016</td>
      <td>0.078426</td>
      <td>0.231442</td>
    </tr>
    <tr>
      <th>836</th>
      <td>PRGR.H14CORE.0.P.B</td>
      <td>-0.090640</td>
      <td>0.068276</td>
      <td>0.158916</td>
    </tr>
    <tr>
      <th>293</th>
      <td>GCR.H14CORE.0.PS.A</td>
      <td>-0.085943</td>
      <td>0.034013</td>
      <td>0.119956</td>
    </tr>
    <tr>
      <th>491</th>
      <td>KLF16.H14CORE.1.P.B</td>
      <td>-0.082473</td>
      <td>-0.244051</td>
      <td>-0.161577</td>
    </tr>
    <tr>
      <th>553</th>
      <td>MAZ.H14CORE.1.P.B</td>
      <td>-0.072905</td>
      <td>-0.176444</td>
      <td>-0.103539</td>
    </tr>
    <tr>
      <th>1559</th>
      <td>ZNF48.H14CORE.0.PSG.A</td>
      <td>-0.071781</td>
      <td>-0.042426</td>
      <td>0.029356</td>
    </tr>
    <tr>
      <th>1122</th>
      <td>VEZF1.H14CORE.1.P.B</td>
      <td>-0.070181</td>
      <td>-0.129904</td>
      <td>-0.059723</td>
    </tr>
    <tr>
      <th>506</th>
      <td>KMT2A.H14CORE.0.P.B</td>
      <td>-0.057221</td>
      <td>-0.300223</td>
      <td>-0.243002</td>
    </tr>
    <tr>
      <th>557</th>
      <td>MCR.H14CORE.0.S.B</td>
      <td>-0.057109</td>
      <td>0.027812</td>
      <td>0.084921</td>
    </tr>
    <tr>
      <th>1342</th>
      <td>ZN467.H14CORE.0.P.C</td>
      <td>-0.056477</td>
      <td>-0.088766</td>
      <td>-0.032288</td>
    </tr>
  </tbody>
</table>
</div>


## 9. ChIP-Atlas (atlas) differential enrichment

Build a 1 kb bin × track index over the ChIP-Atlas hg38 collection, then use
each region set as the other's reference to find TF tracks differentially
enriched between AR+F and AR−F.

> Building the index over the full ChIP-Atlas is heavy (minutes, several GB
> RAM). `atlas.save("chipatlas_hg38_1kb.npz")` / `Atlas.load(...)` reuse it
> across sessions.


```python
with timer("build ChIP-Atlas index (1kb)"):
    atlas = Atlas.make(
        GIGGLE_DIR, chromsizes=CHROMSIZES,
        meta=GIGGLE_META,
        name_pattern=r"([^.]+)",                       # SRX23002840.20.bed.gz -> SRX23002840
        meta_columns=["id", "antigen", "class", "cell_line"],  # meta TSV is header-less
        meta_id_col="id",
    )

with timer("atlas search AR+F vs AR-F"):
    enr_pf = atlas.search(ARpF, ref=ARmF)    # enriched in AR+F over AR-F
    enr_mf = atlas.search(ARmF, ref=ARpF)    # enriched in AR-F over AR+F

cols = ["name", "antigen", "class", "cell_line", "overlaps", "log2_odds", "giggle_score"]
print("TF tracks enriched in AR+F (FOXA1-dependent):")
display(enr_pf.head(15)[cols])
print("TF tracks enriched in AR-F (FOXA1-independent):")
display(enr_mf.head(15)[cols])
```

    ▶ build ChIP-Atlas index (1kb) ...


    [atlas index]: 100%|██████████| 33368/33368 [00:19<00:00, 1726.86it/s]


    ✓ build ChIP-Atlas index (1kb) — 60.5s
    ▶ atlas search AR+F vs AR-F ...
    ✓ atlas search AR+F vs AR-F — 4.6s
    TF tracks enriched in AR+F (FOXA1-dependent):



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>antigen</th>
      <th>class</th>
      <th>cell_line</th>
      <th>overlaps</th>
      <th>log2_odds</th>
      <th>giggle_score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SRX23002840</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1329</td>
      <td>7.519583</td>
      <td>964.923631</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SRX23002841</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1089</td>
      <td>7.056433</td>
      <td>700.212091</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SRX18285237</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1847</td>
      <td>4.287255</td>
      <td>618.551089</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SRX14353424</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1990</td>
      <td>4.058525</td>
      <td>604.652460</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SRX18285235</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1722</td>
      <td>4.333244</td>
      <td>580.711304</td>
    </tr>
    <tr>
      <th>5</th>
      <td>SRX062360</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1564</td>
      <td>4.546911</td>
      <td>564.964852</td>
    </tr>
    <tr>
      <th>6</th>
      <td>SRX14353423</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1756</td>
      <td>4.155838</td>
      <td>548.736062</td>
    </tr>
    <tr>
      <th>7</th>
      <td>SRX5577144</td>
      <td>Epitope tags</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1779</td>
      <td>4.090175</td>
      <td>540.034633</td>
    </tr>
    <tr>
      <th>8</th>
      <td>SRX1885188</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>2044</td>
      <td>3.775184</td>
      <td>533.691700</td>
    </tr>
    <tr>
      <th>9</th>
      <td>SRX18285236</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1470</td>
      <td>4.498439</td>
      <td>515.002750</td>
    </tr>
    <tr>
      <th>10</th>
      <td>SRX1885186</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1936</td>
      <td>3.774294</td>
      <td>503.551086</td>
    </tr>
    <tr>
      <th>11</th>
      <td>SRX1885187</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1878</td>
      <td>3.818223</td>
      <td>499.236335</td>
    </tr>
    <tr>
      <th>12</th>
      <td>SRX1885185</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1972</td>
      <td>3.724534</td>
      <td>499.063986</td>
    </tr>
    <tr>
      <th>13</th>
      <td>SRX1212235</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1936</td>
      <td>3.749926</td>
      <td>496.615360</td>
    </tr>
    <tr>
      <th>14</th>
      <td>SRX6878606</td>
      <td>FOXA1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>1877</td>
      <td>3.763680</td>
      <td>484.060928</td>
    </tr>
  </tbody>
</table>
</div>


    TF tracks enriched in AR-F (FOXA1-independent):



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>antigen</th>
      <th>class</th>
      <th>cell_line</th>
      <th>overlaps</th>
      <th>log2_odds</th>
      <th>giggle_score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SRX3070471</td>
      <td>NR3C1</td>
      <td>Breast</td>
      <td>MCF 10A</td>
      <td>151</td>
      <td>1.985641</td>
      <td>57.274594</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SRX3070475</td>
      <td>NR3C1</td>
      <td>Breast</td>
      <td>MCF 10A</td>
      <td>120</td>
      <td>2.021899</td>
      <td>49.016228</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SRX21439743</td>
      <td>ESR1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>452</td>
      <td>1.518987</td>
      <td>46.612336</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SRX3070473</td>
      <td>NR3C1</td>
      <td>Breast</td>
      <td>MCF 10A</td>
      <td>120</td>
      <td>1.969408</td>
      <td>45.866518</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SRX306516</td>
      <td>AR</td>
      <td>Prostate</td>
      <td>DU 145</td>
      <td>161</td>
      <td>1.690385</td>
      <td>39.808143</td>
    </tr>
    <tr>
      <th>5</th>
      <td>SRX19970679</td>
      <td>AR</td>
      <td>Prostate</td>
      <td>PC-346C</td>
      <td>167</td>
      <td>1.537133</td>
      <td>31.894465</td>
    </tr>
    <tr>
      <th>6</th>
      <td>SRX8520802</td>
      <td>NR3C1</td>
      <td>Breast</td>
      <td>HCC1187</td>
      <td>173</td>
      <td>1.506006</td>
      <td>31.057426</td>
    </tr>
    <tr>
      <th>7</th>
      <td>SRX21439744</td>
      <td>ESR1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>336</td>
      <td>1.292770</td>
      <td>30.240255</td>
    </tr>
    <tr>
      <th>8</th>
      <td>SRX3070469</td>
      <td>NR3C1</td>
      <td>Breast</td>
      <td>MCF 10A</td>
      <td>53</td>
      <td>2.255766</td>
      <td>30.173945</td>
    </tr>
    <tr>
      <th>9</th>
      <td>SRX21439742</td>
      <td>ESR1</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>295</td>
      <td>1.284194</td>
      <td>28.358706</td>
    </tr>
    <tr>
      <th>10</th>
      <td>SRX083218</td>
      <td>AR</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>236</td>
      <td>1.337458</td>
      <td>28.128496</td>
    </tr>
    <tr>
      <th>11</th>
      <td>SRX1067071</td>
      <td>AR</td>
      <td>Prostate</td>
      <td>LHSAR</td>
      <td>319</td>
      <td>1.253582</td>
      <td>27.290928</td>
    </tr>
    <tr>
      <th>12</th>
      <td>SRX1067070</td>
      <td>AR</td>
      <td>Prostate</td>
      <td>LHSAR</td>
      <td>340</td>
      <td>1.245091</td>
      <td>27.235381</td>
    </tr>
    <tr>
      <th>13</th>
      <td>SRX083219</td>
      <td>AR</td>
      <td>Prostate</td>
      <td>LNCAP</td>
      <td>233</td>
      <td>1.311355</td>
      <td>26.426950</td>
    </tr>
    <tr>
      <th>14</th>
      <td>SRX3630817</td>
      <td>NR3C1</td>
      <td>Uterus</td>
      <td>Ishikawa</td>
      <td>53</td>
      <td>2.092061</td>
      <td>25.168441</td>
    </tr>
  </tbody>
</table>
</div>


## 10. Browser view — `chr19:50,792,009-50,923,669`

All six conditions as signal tracks (ATAC replicates passed as a **list** of
bigwigs are averaged into one track), their peak calls, and gene models. The
region string is 0-based, half-open like every table.


```python
region = "chr19:50,792,009-50,923,669"
tracks = {
    "ATAC 0h":   [BW["ATAC_0h_r1"], BW["ATAC_0h_r2"]],   # list -> averaged
    "ATAC 4h":   [BW["ATAC_4h_r1"], BW["ATAC_4h_r2"]],
    "AR 0h":     BW["AR_0h"],
    "AR 4h":     BW["AR_4h"],
    "FOXA1 0h":  BW["FOXA1_0h"],
    "FOXA1 4h":  BW["FOXA1_4h"],
    "AR 4h peaks":    PEAK["AR_4h"],
    "FOXA1 4h peaks": PEAK["FOXA1_4h"],
    "genes":     genes,
}
# share the y-axis within each assay so the 0h -> 4h gain is honest
with timer("browser render"):
    fig, axes = browser(region, tracks, bw_n_bins=2000, figsize=(11, None),
                        bw_share=[["ATAC 0h", "ATAC 4h"],
                                  ["AR 0h", "AR 4h"],
                                  ["FOXA1 0h", "FOXA1 4h"]])
```

    ▶ browser render ...
    ✓ browser render — 1.3s



    
![figure]({{ '/assets/images/ar_foxa1/ar_foxa1_lncap_23_1.png' | relative_url }})
    



```python

```
