"""Example: IGV-like browser view around the Nanog locus (mm10).

Renders a multi-track figure at ``chr6:122,600,000-122,800,000`` with:

    - mESC ATAC-seq peaks          (narrowPeak, sampled)
    - mESC ATAC-seq coverage       (bigWig, read directly by random access)
    - mESC H3K27ac HiChIP loops    (bedpe, sampled)
    - GENCODE vM25 protein-coding  (gtf, sampled)

The small text tracks under ``examples/data/`` were sampled from the
original files on the lab share; the bigWig is referenced by absolute
path because the reader supports fast random access, so no sampling is
needed.

Run:
    micromamba run -n notebook_cn1 python examples/browser_example.py
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")

from genomeblocks import browser, Loci, Genes
from genomeblocks.bedpe import read_bedpe


HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# ── sampled text tracks (mm10, chr6:122,500,000-122,900,000) ────────────
NARROWPEAK = '/groups/lackgrp/projects/inv-berkay-nfatacseq/mm10-atac/results/bwa/merged_replicate/macs2/narrow_peak/C57BL6_ESC_GSE113431.mRp.clN_peaks.narrowPeak'
BEDPE      = os.path.join(DATA, "Nanog-promoter-HiChIP.bedpe")
GTF        = '/groups/lackgrp/genomeAnnotations/mm10/gencode.vM25.protein_coding.annotation.gtf'

# ── full bigWig on the lab share (random-access read, no sampling) ──────
ATAC_BW = (
    "/groups/lackgrp/projects/inv-berkay-nfatacseq/mm10-atac/"
    "results/bwa/merged_replicate/bigwig/"
    "C57BL6_ESC_GSE113431.mRp.clN.bigWig"
)

# Region: 200 kb window around Nanog (chr6:122,707,489-122,714,633)
REGION = ('chr6',122286666, 122902344)


def main(out_svg: str = os.path.join(HERE, "browser_Nanog.svg")) -> None:
    # Load objects once so they can be reused across regions.
    peaks = Loci.make(NARROWPEAK)
    loops = read_bedpe(BEDPE, verbose=False)
    genes = Genes.make(GTF)

    fig, _ = browser(
        REGION,
        tracks={
            "HiChIP loops": loops,
            
            "ATAC signal":  ATAC_BW,
            "ATAC peaks":   peaks,

            "Genes":        genes,
        },
        figsize=(12, None),                 # auto-height
        colors={
            "HiChIP loops": "#DA0000",
            "ATAC peaks":   "#DA0000",
            "ATAC signal":  "#DA0000",
            
        },
        bw_n_bins=120,
    )
    fig.savefig(out_svg, bbox_inches="tight")
    print(f"wrote {out_svg}")


if __name__ == "__main__":
    main()
