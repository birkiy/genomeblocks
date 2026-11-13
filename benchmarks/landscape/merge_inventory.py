#!/usr/bin/env python3
"""Merge the per-theme inventories into landscape/tools.json.

    python merge_inventory.py out.json part1.json part2.json ...

Tools that appear in several parts (bedtools, rtracklayer, HOMER, ...) are
merged by name: themes and capability flags are unioned, the other fields
keep the first non-null value. The genomeblocks entry is added last from
GENOMEBLOCKS below (the only hand-written record: its flags come from the API).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ALIASES = {  # spelling differences between the parts -> one name
    "genomicranges + iranges": "GenomicRanges + IRanges",
    "genomicranges": "GenomicRanges + IRanges",
    "biocpy genomicranges": "BiocPy genomicranges / iranges",
    "biocpy genomicranges / iranges": "BiocPy genomicranges / iranges",
    "bigtools / pybigtools": "pybigtools",
    "bigtools/pybigtools": "pybigtools",
    "pybigtools / bigtools": "pybigtools",
    "ucsc kent utilities": "UCSC kent utilities",
    "ucsc bigwigaverageoverbed": "UCSC kent utilities",
    "meme suite": "MEME suite",
    "meme-suite": "MEME suite",
    "biopython bio.motifs": "Biopython",
    "biopython": "Biopython",
    "homer annotatepeaks.pl": "HOMER",
    "homer": "HOMER",
    "genomicinteractions": "GenomicInteractions",
    "interactionset": "GenomicInteractions",
    "pydeseq2": "PyDESeq2",
    "rnanorm": "RNAnorm",
    "ailist (python)": "ailist (Python)",
    "chipseeker (covplot)": "ChIPseeker",
    "scipy": "SciPy (csgraph)",
    "scipy (csgraph)": "SciPy (csgraph)",
    "igraph (r)": "igraph",
    "data.table::foverlaps": "data.table",
    "data.table": "data.table",
    "python-igraph": "igraph",
    "networkx": "networkx",
    "graph-tool": "graph-tool",
}

GENOMEBLOCKS = {
    "name": "genomeblocks", "aka": "", "themes": ["intervals", "annotation", "signal", "motifs", "cooccupancy", "graph3d", "expression", "viewers"],
    "platform": "Python", "interface": "library", "install": ["pip", "bioconda"],
    "data_model": "columnar tables (numpy columns; the row is the join key) with swappable backends",
    "license": "MIT", "first_published": {"year": 2025, "venue": "", "doi": ""},
    "repo": "https://github.com/birkiy/genomeblocks", "last_release": {"version": "2.0.1", "date": "2026-10"},
    "active_12m": True, "docs": {"site": True, "tutorial": True, "api_reference": True}, "tests_ci": {"tests": True, "ci": True},
    "peer_reviewed": False, "parallel": True,
    "capabilities": {
        # intervals
        "read_bed": True, "overlap": True, "set_algebra": True, "merge": True, "slop": True, "sort": True, "nearest": True,
        "point_lookup_index": True, "tile": True, "liftover": True, "getfasta": True, "dataframe_interop": True, "coverage": True,
        "jaccard_fisher": False, "out_of_core": False, "parallel": True, "columnar_storage": True,
        # annotation
        "parse_gtf": True, "parse_genepred": True, "feature_labels": True, "nearest_tss": True, "nearest_gene": True,
        "custom_promoter_window": True, "isoform_selection": True, "tss_table": True, "export_bed12": True, "export_gtf": True,
        "region_gene_enrichment": False, "distance_plots": False, "gtf_parse_speed_claim": True, "uses_txdb": False, "dataframe_output": True,
        # signal
        "bigwig_read": True, "bigwig_write": True, "region_matrix": True, "region_mean": True, "profile_plot": True, "heatmap_plot": True,
        "multi_track": True, "track_normalisation": True, "bam_coverage": True, "bam_pileup": True, "bedgraph_to_bigwig": True,
        "super_enhancer": True, "hichip_shortrange_track": True, "strand_aware": False, "in_memory_api": True, "gpu": False,
        # motifs
        "parse_jaspar": True, "parse_meme": True, "parse_transfac": True, "pwm_scan": True, "pvalue_threshold": True,
        "per_position_profile": True, "enrichment_vs_background": True, "differential_enrichment": True, "motif_comparison": True,
        "motif_clustering": True, "logo_plot": True, "de_novo_discovery": False, "multi_engine": True, "web_only": False,
        # cooccupancy
        "index_collection": True, "enrichment_search": True, "fisher_test": True, "permutation_null": True, "reference_background": True,
        "metadata_search": True, "pairwise_cooccupancy": False, "giggle_score": True, "scale_10k_tracks": True, "local_install": True,
        "bootstrap_resampling_of_query": True,
        # graph3d
        "read_bedpe": True, "loop_region_overlap": True, "pairs_streaming": True, "bin_pairs_to_matrix": True, "read_cool_mcool": True,
        "read_hic_juicer": False, "oe_normalisation": True, "loop_calling": True, "contact_graph": True, "graph_algorithms": True,
        "hub_mining": True, "anchor_gene_annotation": True, "hichip_support": True, "graph_export": True, "arc_visualisation": True,
        # expression
        "read_salmon_kallisto": True, "read_rsem": True, "read_featurecounts": True, "read_star_htseq": True, "tx_to_gene": True,
        "tmm": True, "size_factors": True, "cpm_logcpm": True, "filter_by_expr": True, "de_test": False, "lfc_shrinkage": False,
        "aligned_to_gene_models": True, "anndata_or_sce": True, "batch_correction": False, "python_native": True, "r_native": False,
        "multi_format_sniffing": True,
        # viewers
        "static_region_plot": True, "bigwig_tracks": True, "bam_tracks": True, "gene_tracks": True, "loop_arcs": True,
        "hic_matrix_track": False, "interactive_html_export": True, "igv_session_or_igvjs": True, "notebook_inline": True,
        "config_file_driven": False, "cli": False, "programmatic_api": True, "multi_region_layout": True, "graph_network_view": True,
        "web_app": False,
    },
    "notes": "Columnar tables for peaks, genes, signal, motifs, enrichment, loops / contact graphs and RNA-seq in one API; de testing is planned for 2.1.",
    "confidence": "high",
}

EXTRA = [  # engines the benchmarks run that the per-theme research did not list as tools of their own
    {"name": "Biostrings", "themes": ["motifs"], "platform": "R", "interface": "library", "install": ["Bioconductor", "bioconda"],
     "data_model": "DNAStringSet / XStringViews", "license": "Artistic-2.0", "first_published": {"year": None, "venue": "", "doi": ""},
     "repo": "https://github.com/Bioconductor/Biostrings", "last_release": {"version": "2.78", "date": "2026-04"}, "active_12m": True,
     "docs": {"site": True, "tutorial": True, "api_reference": True}, "tests_ci": {"tests": True, "ci": True}, "peer_reviewed": False, "parallel": False,
     "capabilities": {"parse_jaspar": False, "parse_meme": False, "parse_transfac": False, "pwm_scan": True, "pvalue_threshold": False,
                      "per_position_profile": False, "enrichment_vs_background": False, "differential_enrichment": False, "motif_comparison": False,
                      "motif_clustering": False, "logo_plot": False, "de_novo_discovery": False, "multi_engine": False, "parallel": False, "gpu": False,
                      "in_memory_api": True, "web_only": False},
     "notes": "matchPWM / countPWM: one PWM against one sequence at an absolute score; the primitive under motifmatchr and TFBSTools.", "confidence": "high"},
    {"name": "SciPy (csgraph)", "themes": ["graph3d"], "platform": "Python", "interface": "library", "install": ["pip", "conda-forge"],
     "data_model": "sparse adjacency matrix", "license": "BSD-3-Clause", "first_published": {"year": 2020, "venue": "Nature Methods", "doi": "10.1038/s41592-019-0686-2"},
     "repo": "https://github.com/scipy/scipy", "last_release": {"version": "1.17", "date": "2026"}, "active_12m": True,
     "docs": {"site": True, "tutorial": True, "api_reference": True}, "tests_ci": {"tests": True, "ci": True}, "peer_reviewed": True, "parallel": False,
     "capabilities": {"read_bedpe": False, "loop_region_overlap": False, "pairs_streaming": False, "bin_pairs_to_matrix": False, "read_cool_mcool": False,
                      "read_hic_juicer": False, "oe_normalisation": False, "loop_calling": False, "contact_graph": False, "graph_algorithms": True,
                      "hub_mining": False, "anchor_gene_annotation": False, "hichip_support": None, "graph_export": False, "arc_visualisation": False,
                      "parallel": False, "in_memory_api": True},
     "notes": "connected_components and sparse linear algebra on an adjacency matrix: genomeblocks' pip-default graph backend.", "confidence": "high"},
    {"name": "pyfaidx", "themes": ["intervals"], "platform": "Python", "interface": "library", "install": ["pip", "bioconda"],
     "data_model": "indexed FASTA (dict-like)", "license": "BSD-3-Clause", "first_published": {"year": 2017, "venue": "PeerJ PrePrints", "doi": "10.7287/peerj.preprints.970v1"},
     "repo": "https://github.com/mdshw5/pyfaidx", "last_release": {"version": "0.9", "date": "2025"}, "active_12m": True,
     "docs": {"site": True, "tutorial": True, "api_reference": True}, "tests_ci": {"tests": True, "ci": True}, "peer_reviewed": False, "parallel": False,
     "capabilities": {"read_bed": True, "overlap": False, "set_algebra": False, "merge": False, "slop": False, "sort": False, "nearest": False,
                      "point_lookup_index": False, "tile": False, "liftover": False, "getfasta": True, "dataframe_interop": False, "coverage": False,
                      "jaccard_fisher": False, "out_of_core": True, "parallel": False, "columnar_storage": False},
     "notes": "Random access to FASTA by .fai index; its `faidx` CLI extracts BED regions (bedtools getfasta's library-side counterpart).", "confidence": "high"},
    {"name": "GenomicAlignments", "themes": ["signal", "expression"], "platform": "R", "interface": "library", "install": ["Bioconductor", "bioconda"], "data_model": "GAlignments (S4) from BAM", "license": "Artistic-2.0", "first_published": {"year": 2013, "venue": "PLoS Computational Biology", "doi": "10.1371/journal.pcbi.1003118"}, "repo": "https://github.com/Bioconductor/GenomicAlignments", "last_release": {"version": None, "date": ""}, "active_12m": True, "docs": {"site": True, "tutorial": True, "api_reference": True}, "tests_ci": {"tests": True, "ci": True}, "peer_reviewed": True, "parallel": True, "capabilities": {"bam_coverage": True, "multi_bam_region_counts": True, "count_reads_from_bam": True, "bigwig_read": False, "bigwig_write": False, "region_matrix": False, "in_memory_api": True}, "notes": "Reads BAM alignments into R; summarizeOverlaps counts reads per feature (Union / IntersectStrict / IntersectNotEmpty modes). Added 2026-10-23 as an engine the landscape benchmark runs.", "confidence": "high"},
    {"name": "pychromVAR", "themes": ["motifs"], "platform": "Python", "interface": "library", "install": ["pip"], "data_model": "AnnData / MuData (peaks x cells)", "license": "MIT", "first_published": {"year": None, "venue": "", "doi": ""}, "repo": "https://github.com/pinellolab/pychromVAR", "last_release": {"version": "0.0.4", "date": ""}, "active_12m": None, "docs": {"site": True, "tutorial": True, "api_reference": True}, "tests_ci": {"tests": None, "ci": None}, "peer_reviewed": None, "parallel": True, "capabilities": {"per_sample_motif_activity": True, "pwm_scan": True, "pvalue_threshold": True, "in_memory_api": True, "de_novo_discovery": False}, "notes": "Python port of chromVAR (GC-matched background peaks, motif deviation z-scores). Added 2026-10-23 as an engine the landscape benchmark runs.", "confidence": "medium"},
]


def canon(name: str) -> str:
    key = re.sub(r"\s+", " ", name.strip()).lower()
    return ALIASES.get(key, name.strip())


def merge(parts):
    themes, tools = {}, {}
    for p in parts:
        d = json.load(open(p))
        for k, v in d["themes"].items():
            themes.setdefault(k, v)
        for t in d["tools"]:
            n = canon(t["name"])
            if n not in tools:
                t = dict(t); t["name"] = n
                t["themes"] = list(dict.fromkeys(t.get("themes") or []))
                tools[n] = t
                continue
            cur = tools[n]
            cur["themes"] = list(dict.fromkeys((cur.get("themes") or []) + (t.get("themes") or [])))
            caps = dict(t.get("capabilities") or {})
            caps.update({k: v for k, v in (cur.get("capabilities") or {}).items() if v is not None})
            cur["capabilities"] = caps
            for k, v in t.items():
                if k in ("name", "themes", "capabilities"):
                    continue
                if cur.get(k) in (None, "", [], {}) and v not in (None, "", [], {}):
                    cur[k] = v
                elif isinstance(v, dict) and isinstance(cur.get(k), dict):
                    cur[k] = {kk: (cur[k].get(kk) if cur[k].get(kk) is not None else vv) for kk, vv in {**v, **cur[k]}.items()}
            if t.get("notes") and t["notes"] not in (cur.get("notes") or ""):
                cur["notes"] = (cur.get("notes") or "").rstrip(". ") + ("; " if cur.get("notes") else "") + t["notes"]
    for t in EXTRA:
        tools.setdefault(t["name"], t)
    tools["genomeblocks"] = GENOMEBLOCKS
    return {"themes": themes, "tools": [tools["genomeblocks"]] + [tools[k] for k in sorted(tools, key=str.lower) if k != "genomeblocks"],
            "meta": {"compiled": "2026-10-22", "parts": [Path(p).name for p in parts]}}


if __name__ == "__main__":
    out, *parts = sys.argv[1:]
    d = merge(parts)
    json.dump(d, open(out, "w"), indent=1)
    print(f"{len(d['tools'])} tools, themes {list(d['themes'])}")
