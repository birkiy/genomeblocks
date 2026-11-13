#!/usr/bin/env python3
"""Add the confirmed jobs that had no checklist key to landscape/tools.json.

    python merge_inventory.py tools.json part1.json ...   # the first survey (2026-10-22)
    python resurvey_jobs.py                                # this step, idempotent

The 160 jobs in jobs.json (10 core + 10 auxiliary per theme, confirmed by the
author on 2026-10-23) were chosen from the field, not from genomeblocks; 62 of
them had no item in the first survey's checklists. This script appends those 62
keys to each theme's checklist, fills them for every tool listed in the theme,
and writes the key names back into jobs.json.

One code per key, in NEW_KEYS order:
  y  true   the tool documents a function, command or option that does the job
  n  false  it does not (composing the job from the tool's primitives does not count)
  ?  null   not established: the docs were not checked closely enough to say, or the
            tool is offline; unknown, never counted as "no"

The values were filled on 2026-10-23 from model knowledge (Claude) of each tool's
documented interface, not from a fresh crawl of every repository; "?" wherever
that knowledge was not firm. genomeblocks' values come from its source code; its
partial jobs (jobs.json) are "n" here, by the same no-composition rule.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

NEW_KEYS = {  # theme -> [(key, jobs.json id)]
    "intervals": [("shuffle_regions", "A5"), ("multi_set_intersect", "A6"), ("window_search", "A7")],
    "annotation": [("derived_gene_regions", "C2"), ("genes_in_window", "C5"), ("transcript_level_annotation", "C8"),
                   ("gene_set_function_enrichment", "C10"), ("gene_id_mapping", "A1"), ("gene_filters", "A5"),
                   ("custom_annotation_sets", "A6"), ("chrom_naming_conversion", "A8"), ("sequence_context_gc_cpg", "A10")],
    "signal": [("multi_bam_region_counts", "C4"), ("sample_correlation_qc", "C6"), ("track_comparison", "C7"),
               ("paired_end_fragments", "A5"), ("blacklist_masking", "A6"), ("peak_set_qc", "A7"), ("signal_smoothing", "A8"),
               ("bigwig_merge", "A9")],
    "motifs": [("library_scale_scan", "C2"), ("background_construction", "C4"), ("per_sample_motif_activity", "C9"),
               ("motif_format_export", "A3"), ("motif_model_conversion", "A4"), ("motif_hits_genomic", "A5"),
               ("motif_cooccurrence_spacing", "A6"), ("variant_motif_effect", "A7"), ("masked_scanning", "A9"),
               ("attribution_motifs", "A10")],
    "cooccupancy": [("matched_background", "C6"), ("multiple_testing_correction", "C10"), ("local_zscore_shift", "A2"),
                    ("block_bootstrap", "A3"), ("effect_sizes", "A4"), ("set_similarity_search", "A5"),
                    ("overlap_set_visualisation", "A6"), ("enrichment_plots", "A7"), ("aggregate_by_metadata", "A8"),
                    ("public_atlas_bundles", "A9"), ("per_track_overlap_counts", "A10")],
    "graph3d": [("pairs_parse_dedup", "C1"), ("matrix_balancing", "C3"), ("tad_calling", "C7"), ("compartment_calling", "C8"),
                ("differential_contacts", "C9"), ("aggregate_pileups", "A1"), ("virtual_4c", "A8"), ("contact_qc", "A9")],
    "expression": [("count_reads_from_bam", "C1"), ("dispersion_estimation", "C6"), ("variance_stabilising_transform", "C9"),
                   ("sample_qc_pca", "A3"), ("ma_volcano_plots", "A4"), ("independent_filtering_fdr", "A5"),
                   ("gene_set_tests", "A6"), ("pseudobulk", "A10")],
    "viewers": [("interactive_browsing", "C7"), ("variant_tracks", "A7"), ("region_highlights", "A8"),
                ("linked_views", "A9"), ("publication_export", "A10")],
}

FILL = {
    # shuffle_regions · multi_set_intersect · window_search
    "intervals": {
        "genomeblocks": "nnn", "AIList": "nnn", "ailist (Python)": "nnn",
        "BEDOPS": "nyy",                       # bedops -i/-u/-p over many files; bedmap --range
        "bedr": "?y?",                         # bedr.join.multiple; other bedtools methods via the generic wrapper
        "bedtk": "nnn", "bedtools": "yyy", "biobear": "nnn",
        "BiocPy genomicranges / iranges": "??y",
        "bioframe": "?n?", "bx-python": "??y", "cgranges": "nnn", "coitrees": "nnn",
        "data.table": "nny",                   # foverlaps(maxgap=)
        "GenomeKit": "nn?",
        "GenomicRanges + IRanges": "nny",      # findOverlaps(maxgap=)
        "GIGGLE": "nnn", "granges": "?n?", "gtars": "?n?", "HTSeq GenomicArray": "nnn", "IGD": "nnn",
        "intervaltree": "nnn", "ncls": "nnn",
        "plyranges": "nny",                    # join_overlap_*(maxgap=)
        "polars-bio": "nn?", "pybedlite": "nnn", "pybedtools": "yyy", "pyfaidx": "nnn",
        "pyranges (0.x)": "?yy",               # count_overlaps over a dict of ranges; join(slack=)
        "pyranges (1.x)": "??y",
        "rtracklayer": "nnn", "superintervals": "nnn", "tabix / htslib": "nnn",
        "valr": "y?y",                         # bed_shuffle, bed_window
    },
    # derived_gene_regions · genes_in_window · transcript_level_annotation · gene_set_function_enrichment ·
    # gene_id_mapping · gene_filters · custom_annotation_sets · chrom_naming_conversion · sequence_context_gc_cpg
    "annotation": {
        "genomeblocks": "ynnnnnnnn",
        "AGAT": "ynnnnyn?n",
        "annotatr": "ynyny?ynn",
        "ArchR": "y?nn?n?n?",
        "BCBio.GFF": "nnnnnnnnn",
        "bedtools": "nynnnnyny",               # window -w; annotate -files; nuc
        "biobear": "nnnnn?nnn",
        "BiocPy genomicranges / iranges": "nynnnnnnn",
        "chipenrich": "n?ny?n?nn",
        "ChIPpeakAnno": "yyyyy?y??",           # assignChromosomeRegion, bindingRegion, getEnrichedGO, addGeneIDs
        "ChIPseeker": "yyynynnnn",             # addFlankGeneInfo, level=, annoDb; GO is delegated to clusterProfiler
        "ensembldb": "ynnnyynyn",              # GeneBiotypeFilter / TxSupportLevel; seqlevelsStyle on EnsDb
        "genomation": "ynynnnynn",
        "GenomeKit": "y?yn??n?n",
        "GenomicDistributions": "ynnnnnyny",   # calcPartitions with custom partitions; calcGCContent
        "GenomicFeatures + txdbmaker": "ynynnnnyn",
        "gffpandas": "nnnnnynnn",
        "gffread": "?nnnnynyn",                # -C / --ids / region filters; -m chromosome mapping
        "gffutils": "ynnnnynnn",               # create_introns, interfeatures
        "goldmine": "y??n?nynn",
        "GREAT (web)": "nynynnnnn",
        "greatpy": "nynynnnnn",
        "gtfparse": "nnnnn?nnn",
        "GTFtools": "yn?nn?nnn",
        "HOMER": "y?yyynyny",                  # -go, gene ID columns, -ann, -CpG / -gc
        "HTSeq GenomicArray": "nnnnnnnnn",
        "PAVIS": "y??n?nnnn",
        "polars-bio": "n?nnnnnnn",
        "pyensembl": "??yny?n?n",
        "pygtftk": "y?yn?y?yn",
        "pyranges (0.x)": "yy?nnnnnn",
        "pyranges (1.x)": "?y?nnnnnn",
        "rGREAT": "nyny?nnnn",
        "rtracklayer": "nnnnnnnnn",
        "Signac": "nnnnnnnny",                 # RegionStats GC per peak
        "UROPA": "nyynny?nn",
        "valr": "nynnnnnnn",
    },
    # multi_bam_region_counts · sample_correlation_qc · track_comparison · paired_end_fragments ·
    # blacklist_masking · peak_set_qc · signal_smoothing · bigwig_merge
    "signal": {
        "genomeblocks": "nnnnnnnn",
        "bamsignals": "ynnynnnn",
        "bedtools": "ynnynnn?",                # multicov; bamtobed -bedpe / genomecov -pc
        "bigwig-loader": "nnnnnnnn",
        "bwtool": "nnnnyn?n",                  # bwtool remove mask
        "ChIPseeker": "nnnnnnnn",
        "CoverageView": "nn?nnnnn",
        "CREAM": "nnnnnnnn",
        "csaw": "yn?yy?nn",                    # windowCounts / regionCounts; readParam(pe=, discard=)
        "deepTools": "yyyyyyyy",
        "DiffBind": "yynyyynn",                # dba.count, dba.plotPCA, dba.blacklist, FRiP
        "EnrichedHeatmap": "nnnnnnyn",
        "FitHiChIP": "nnnnnnnn",
        "genomation": "?nn?nnnn",
        "HiC-Pro": "nnnnnnnn", "HiChIP-Peaks": "nnnnnnnn", "hichipper": "nnnnnnnn",
        "HOMER": "yn????nn",
        "imPROSE": "nnnnnnnn", "LILY": "nnnnnnnn",
        "MAnorm2": "n?nnnnnn",
        "megadepth": "ynn?nnnn",
        "metagene2": "?n??nn?n",
        "mosdepth": "?nnynnnn",                # --fragment-mode
        "ngs.plot": "nny?nnyn",
        "profileplyr": "nnnnnnnn",
        "pybigtools": "nnnnnnnn", "pyBigWig": "nnnnnnnn",
        "pysam": "nnn?nnnn",
        "ROSE": "nnnnnnnn", "ROSE2": "nnnnnnnn", "rtracklayer": "nnnnnnnn",
        "samtools": "ynnynnnn",                # bedcov over many BAMs; fixmate / stats insert sizes
        "SEanalysis 2.0": "nnnnnnnn",
        "seqplots": "nnnnnn?n",
        "soGGi": "nnnynnnn",
        "TOBIAS": "nnn?yn??",
        "trackplot": "nnnnnnnn",
        "UCSC kent utilities": "nynnnnny",     # bigWigCorrelate, bigWigMerge
        "WiggleTools": "nyyn?nyy",             # pearson, diff/ratio, smooth, mean/sum
    },
    # library_scale_scan · background_construction · per_sample_motif_activity · motif_format_export ·
    # motif_model_conversion · motif_hits_genomic · motif_cooccurrence_spacing · variant_motif_effect ·
    # masked_scanning · attribution_motifs
    "motifs": {
        "genomeblocks": "ynnnynnnyn",
        "ArchR": "yyynnynnnn",                 # addBgdPeaks, chromVAR deviations, getPositions
        "Biopython": "nnnyynnnnn",
        "Biostrings": "nnnnyynn?n",
        "chromVAR": "yyynnnynnn",              # getBackgroundPeaks; getAnnotationSynergy
        "ggseqlogo": "nnnnnnnnnn",
        "GimmeMotifs": "yyy?yy?nnn",           # gimme background, maelstrom, scan -b
        "gReLU": "?nnnn?n?ny",
        "HOMER": "yynnny?nyn",                 # GC-matched background, scanMotifGenomeWide, -mask
        "JASPAR web tools": "nnny??nnnn",
        "lightmotif": "?nnnynnnnn",
        "Logomaker": "nnnnynnnnn",
        "MACRO-APE / PERFECTOS-APE": "nnnn?nnynn",  # PERFECTOS-APE: SNP effect on motif affinity
        "MEME suite": "yynyyyyn?n",            # fasta-shuffle-letters, *2meme converters, FIMO genomic coords, SpaMo
        "memes": "yynnnynnnn",
        "memesuite-lite": "?nnnn?nnnn",
        "monaLisa": "yynnnynnnn",
        "MOODS": "ynnnynnnnn",
        "MotifDb": "nnnynnnnnn",
        "motifmatchr": "ynnnnynnnn",
        "motifStack": "nnnn?nnnnn",
        "motifTestR": "?ynnnynnnn",
        "oPOSSUM 3": "n?nnnnynnn",
        "Pscan": "nnnnnnnnnn",
        "PWMEnrich": "yynny?nnnn",
        "PWMScan": "?nnn?ynnnn",
        "pycistarget": "yynnn?nnnn",
        "pyJASPAR": "nnnnnnnnnn",
        "pymemesuite": "??nn?nnnnn",
        "RcisTarget": "ynnnn?nnnn",
        "RSAT": "yynyy??y?n",                  # convert-matrix, variation-scan
        "seqLogo": "nnnnnnnnnn",
        "Signac": "yyynn?nnnn",                # MatchRegionStats, RunChromVAR
        "SnapATAC2": "??nnn?nnnn",
        "STAMP": "nnnnnnnnnn",
        "tangermeme": "nynnnn?yny",
        "TF-MoDISco-lite": "nnnnnnnnny",
        "TFBSTools": "?nn?yynnnn",
        "TOBIAS": "y??y?y?nnn",                # TFBScan, FormatMotifs
        "universalmotif": "yynyy??nyn",        # shuffle_sequences, write_*, convert_type, mask_seqs
        "Vierstra motif archetypes": "nnn?nynnnn",
    },
    # matched_background · multiple_testing_correction · local_zscore_shift · block_bootstrap · effect_sizes ·
    # set_similarity_search · overlap_set_visualisation · enrichment_plots · aggregate_by_metadata ·
    # public_atlas_bundles · per_track_overlap_counts
    "cooccupancy": {
        "genomeblocks": "nnnnyynnnny",
        "BART2": "??nn?nn??yn",
        "Bedshift": "nnnnnnnnnnn",
        "bedtk": "nnnnnnnnnnn",
        "bedtools": "nnnny?nnnny",             # jaccard / fisher; intersect -c
        "ChIP-Atlas Enrichment Analysis": "nynnynn?yyy",
        "Cistrome DB Toolkit": "n?nnyyn??y?",
        "Cistrome-GO": "nnnnnnnnnnn",
        "Coloc-stats": "????yn??n??",
        "ColoWeb": "???????????",
        "GAT": "yynnynnnnny",                  # isochores; q-values; fold / l2fold
        "geniml / BEDbase": "nnnnnynnnyn",
        "GenometriCorr": "nnnnyn??nnn",
        "Genomic HyperBrowser": "????y??????",
        "GenomicDistributions": "nnnnnnnnnnn",
        "GIGGLE": "nnnnyynnn?y",
        "IGD": "nnnn?ynnnny",
        "IntervalStats": "nnnn?nnnnnn",
        "Intervene": "nnnnynynnny",            # UpSet / Venn / pairwise heatmaps
        "LOLA": "nynnynnyyyy",
        "LOLAweb": "nynnynnyyyy",
        "nullranges": "ynnynnnnnnn",           # matchRanges, bootRanges
        "OLOGRAM (pygtftk)": "??n?ynyynny",
        "poverlap": "?nnn?nnnnn?",
        "pyranges": "nnnnynnnnny",
        "regioneR": "nnynynnynny",             # localZScore
        "regioneReloaded": "n?ynyn?ynny",
        "RegRegSEA": "???????????",
        "RELI": "yynnynnn?yy",
        "ReMapEnrich": "nynnynny?yy",
        "rGREAT": "nynnynnynnn",
        "TF-COMB": "??nnyny?nny",
        "TFEA.ChIP": "nynnynnyyyn",
        "UniBind Enrichment": "nynnynnyyyy",
    },
    # pairs_parse_dedup · matrix_balancing · tad_calling · compartment_calling · differential_contacts ·
    # aggregate_pileups · virtual_4c · contact_qc
    "graph3d": {
        "genomeblocks": "nnnnnnnn",
        "3D-GNOME": "nnnnnnnn",
        "Activity-by-Contact (ABC)": "nnnnnnnn",
        "bedtools": "nnnnnnnn",
        "C3D": "nnnnnnnn",
        "CHiCAGO": "nnnnnny?",
        "chromosight": "nnynnynn",
        "Cicero": "nnnnnnnn",
        "cLoops2": "?nynyyyy",                 # callDomains, callDiffLoops, agg, estDis
        "cooler": "nynnnnnn",
        "coolpup.py": "nnnnnynn",
        "cooltools": "nnyynyyy",               # insulation, eigs / saddle, pileup, virtual4c, expected
        "dcHiC": "nnnyynnn",
        "diffHic": "yn?nynny",
        "FAN-C": "yyyyyy?y",
        "FitHiChIP": "nnnnynn?",
        "GenomicInteractions": "nnnnnny?",     # viewPoint / plotViewpoint
        "GENOVA": "nnyynyyy",                  # insulation, saddle, APA / ATA, virtual4C, RCP
        "graph-tool": "nnnnnnnn",
        "Hi-Cociety": "n???nnnn",
        "HiC-DC+": "nn?nynnn",
        "HiC-Pro": "yynnnnny",
        "HiC1Dmetrics": "nnyn?nnn",
        "HiCcompare": "nnnnynn?",
        "HiCExperiment / HiContacts": "??yy?yyy",
        "HiCExplorer": "yyyyyyyy",
        "hichipper": "nnnnnnn?",
        "HiCHub": "nnnnynnn",
        "HiCseg": "nnynnnnn",
        "hicstuff": "yy??nnny",
        "hictk": "nynnnnnn",
        "hictkpy": "n?nnnnnn",
        "HiGlass": "nnnnnnnn",
        "Hub-enhancer / enhancer-community analyses (method only)": "nnnnnnnn",
        "iced": "nynnnnnn",
        "igraph": "nnnnnnnn",
        "Juicebox": "nnnnnnnn",
        "Juicer / Juicer Tools": "yyyyyyny",   # KR / SCALE, Arrowhead, eigenvector, HiCCUPSDiff, APA
        "MAPS": "ynnnnnn?",
        "MINE": "????????",
        "multiHiCcompare": "nnnnynn?",
        "Mustache": "nnnnynnn",
        "networkx": "nnnnnnnn",
        "pairtools": "ynnnnnny",
        "Peakachu": "nnnnnnnn",
        "pgltools": "nnnnnnnn",
        "plotgardener": "nnnnnnnn",
        "pyGenomeTracks": "nnnnnnnn",
        "QuIN": "nnnnnnnn",
        "ROSE": "nnnnnnnn",
        "SciPy (csgraph)": "nnnnnnnn",
        "straw": "nnnnnnnn",
    },
    # count_reads_from_bam · dispersion_estimation · variance_stabilising_transform · sample_qc_pca ·
    # ma_volcano_plots · independent_filtering_fdr · gene_set_tests · pseudobulk
    "expression": {
        "genomeblocks": "nnnnnnnn",
        "apeglm": "nnnnnnnn", "ashr": "nnnnnnnn",
        "baySeq": "nynnyynn",
        "BiocPy summarizedexperiment": "nnnnnnnn",
        "bioinfokit": "nnnyynnn",
        "conorm": "nnnnnnnn",
        "dearseq": "nnnnnyyn",
        "DEBrowser": "ny?yyyyn",
        "decoupler": "nnnnynyy",
        "DEGseq": "nnnnyynn",
        "Degust": "ny?yyynn",
        "DESeq (legacy)": "nyyyyynn",
        "DESeq2": "nyyyyynn",
        "diffxpy": "nynnyynn",
        "EBSeq": "nynnnynn",
        "edgePython": "nyn??y??",
        "edgeR": "nynyyyyy",                   # plotMDS, plotMD, camera / roast / goana, Seurat2PB
        "fishpond (swish)": "nnnnyynn",
        "GEO2R": "ny??yynn",
        "glmGamPoi": "ny?nnyny",
        "HTSeq": "ynnnnnnn",
        "iDEP": "nyyyyyyn",
        "InMoose": "ny???ynn",
        "limma": "nyyyyyyn",                   # eBayes variance moderation counted as dispersion estimation
        "nf-core/rnaseq": "ynyynnnn",
        "NOISeq": "nnnyy?nn",
        "OmicVerse": "ny?yyyy?",
        "pertpy": "nynnyyyy",
        "PyDESeq2": "nyynyynn",
        "pytximport": "nnnnnnnn",
        "RNAlysis": "yy?yyyyn",
        "RNAnorm": "nnnnnnnn",
        "Rsubread featureCounts": "ynnnnnnn",
        "RUVSeq": "nnnynnnn",
        "scanpy": "nn?ynyny",
        "sleuth": "n?nyyynn",
        "SummarizedExperiment": "nnnnnnnn",
        "sva / ComBat-seq": "nnnnnnnn",
        "tximeta": "nnnnnnnn", "tximport": "nnnnnnnn",
        "Wilcoxon rank-sum DE (Li et al. 2022)": "nnnnn?nn",
    },
    # interactive_browsing · variant_tracks · region_highlights · linked_views · publication_export
    "viewers": {
        "genomeblocks": "ynyny",
        "Biopython GenomeDiagram": "nn?ny",
        "ChIPseeker": "nnnny",
        "Circos": "n?yny",
        "CoolBox": "y?yny",
        "FAN-C (fancplot)": "?n?ny",
        "GenomeSpy": "yy?y?",
        "genomeview": "n??ny",
        "GENOVA": "nn?ny",
        "ggbio": "ny?ny",
        "ggcoverage": "nyyny",
        "gghic": "nn?ny",
        "Gosling / gos": "yyyyy",
        "Gviz": "n?yny",                       # HighlightTrack
        "HiCExplorer (hicPlotMatrix / hicPlotTADs)": "nn?ny",
        "HiCPlotter": "nnyny",
        "HiGlass": "y?yyy",
        "IGV desktop": "yyyny",
        "igv-notebook": "yyyny",
        "igv-reports": "yy?nn",
        "igv.js": "yyyny",
        "JBrowse 2": "yyyyy",
        "jbrowse-jupyter": "yy???",
        "Juicebox": "ynyyy",
        "karyoploteR": "n?yny",
        "LocusZoom": "yy?yy",
        "plotgardener": "n?yny",               # annoHighlight
        "pyGenomeTracks": "n?yny",
        "pyGenomeViz": "yn?ny",
        "rGenomeTracks": "n??ny",
        "Sushi": "n??ny",
        "svist4get": "nnyny",
        "trackc": "n??ny",
        "trackplot": "nn?ny",
        "tracktables": "?nnnn",
        "trackViewer": "yy?ny",
        "UCSC Genome Browser": "yyyny",
        "WashU Epigenome Browser": "y?y?y",
    },
}

CODE = {"y": True, "n": False, "?": None}


def dump_jobs(jobs):
    """jobs.json in its hand-written layout: one job per line."""
    out = ['{', ' "meta": ' + json.dumps(jobs["meta"], indent=1).replace("\n", "\n ") + ",", ' "themes": {']
    for i, (theme, d) in enumerate(jobs["themes"].items()):
        out.append(f'  "{theme}": {{')
        for tier in ("core", "auxiliary"):
            out.append(f'   "{tier}": [')
            out.append(",\n".join("    " + json.dumps(j, ensure_ascii=False) for j in d[tier]))
            out.append("   ]" + ("," if tier == "core" else ""))
        out.append("  }" + ("," if i < len(jobs["themes"]) - 1 else ""))
    return "\n".join(out + [" }", "}", ""])


def main():
    tp, jp = HERE / "tools.json", HERE / "jobs.json"
    tools, jobs = json.loads(tp.read_text()), json.loads(jp.read_text())
    clash = {k for theme, keys in NEW_KEYS.items() for k, _ in keys
             for other, th in tools["themes"].items() if other != theme and k in th["checklist"]}
    assert not clash, f"new keys already used by another theme: {clash}"
    for theme, keys in NEW_KEYS.items():
        names = [k for k, _ in keys]
        members = [t for t in tools["tools"] if theme in (t.get("themes") or [])]
        fill = FILL[theme]
        missing = [t["name"] for t in members if t["name"] not in fill]
        extra = sorted(set(fill) - {t["name"] for t in members})
        bad = {n: c for n, c in fill.items() if len(c) != len(names) or set(c) - set(CODE)}
        assert not missing and not extra and not bad, (theme, missing, extra, bad)
        cl = tools["themes"][theme]["checklist"]
        cl.extend(k for k in names if k not in cl)
        for t in members:
            caps = t.setdefault("capabilities", {})
            for k, c in zip(names, fill[t["name"]]):
                caps[k] = CODE[c]
        ids = dict(keys)
        for j in jobs["themes"][theme]["core"] + jobs["themes"][theme]["auxiliary"]:
            for k, jid in ids.items():
                if j["id"] == jid and k not in j["checklist"]:
                    j["checklist"].append(k)
    tools["meta"]["resurvey"] = {
        "date": "2026-10-23", "script": "landscape/resurvey_jobs.py", "keys": {th: [k for k, _ in v] for th, v in NEW_KEYS.items()},
        "how": "the jobs.json jobs that had no checklist key, filled for every tool in the theme from model knowledge of each "
               "tool's documented interface; null = not established (unknown), never counted as no"}
    tp.write_text(json.dumps(tools, indent=1) + "\n")
    jp.write_text(dump_jobs(jobs))
    n = sum(len(FILL[th]) * len(v) for th, v in NEW_KEYS.items())
    u = sum(c.count("?") for th in FILL for c in FILL[th].values())
    y = sum(c.count("y") for th in FILL for c in FILL[th].values())
    print(f"[resurvey] {sum(len(v) for v in NEW_KEYS.values())} keys, {n} cells: {y} yes, {n - y - u} no, {u} unknown")


if __name__ == "__main__":
    main()
