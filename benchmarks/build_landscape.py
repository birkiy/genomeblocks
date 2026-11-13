#!/usr/bin/env python3
"""Build the Landscape page and figure: genomeblocks against the field, theme by
theme, scored the way Saelens et al. (Nat. Biotechnol. 2019) scored
trajectory-inference methods — accuracy, scalability, stability, usability,
and one overall score.

    python benchmarks/build_landscape.py

Inputs
  landscape/tools.json     the tool inventory (properties, capability checklists)
  landscape/jobs.json      ten core + ten auxiliary jobs per theme, the tasks
                           that time each, and genomeblocks' status read from
                           its source code (the grade's coverage)
  results/*.json           the main suite (loci, signal, atlas, motifs, pairs,
                           genes, architecture, expression, backends, shortrange,
                           and the loop callers: loops, maps, mustache, hiccups)
                           and results/landscape.json (bench_landscape.py)
Outputs
  docs/benchmarks/landscape.md                the page
  docs/_includes/diagrams/landscape-figure.svg  the figure, inline on the page
  figures/landscape.{svg,png,pdf}              the figure as files
  results/landscape_scores.json                every score, for the record

Scores (all in [0, 1]):
  accuracy     mean agreement of the tool's answers with the reference answer
               of each task it ran (1 = identical; else Pearson / Spearman /
               fraction of rows agreeing — see each task's detail column)
  scalability  per task, 1 - log10(t / t_fastest) / 3 at the task's largest
               size, floored at 0: a decade slower than the fastest tool costs
               a third of the score; averaged over the tasks the tool ran
  stability    fraction of the perturbation tests (rows shuffled, chromosomes
               in another order, samples / genes reordered) that returned the
               identical answer
  usability    a ten-item checklist from the inventory: packaged on an index,
               documentation site, tutorial, API reference, tests, CI, a release
               or commit in the last 12 months, peer-reviewed paper, open
               license, a programmatic (library) interface
  overall      the geometric mean of the aspects a tool has

The figure's per-theme grade (theme_grades): S = C x Q, C the share of the
theme's core jobs a tool does (over the jobs whose status is known), Q its
accuracy-weighted speed on contested tasks (two or more tools timed); tiers
relative to the theme's best S.
"""
from __future__ import annotations

import json
import math
import re
from collections import defaultdict
from html import escape
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
DOCS = HERE.parent / "docs"
TOOLS = json.loads((HERE / "landscape" / "tools.json").read_text())

THEMES = {   # key: (label, colour, short)
    "intervals": ("Interval operations", "#2a9d8f", "Intervals"),
    "annotation": ("Gene annotation", "#264653", "Annotation"),
    "signal": ("Signal (bigWig / BAM)", "#1d7a46", "Signal"),
    "motifs": ("Motifs", "#e76f51", "Motifs"),
    "cooccupancy": ("Co-occupancy / enrichment", "#8a5cf6", "Co-occupancy"),
    "graph3d": ("3D genome & contact graphs", "#5b3fa0", "3D / graph"),
    "expression": ("Expression / DE", "#d4a017", "Expression"),
    "viewers": ("Viewers", "#6c757d", "Viewers"),
}
ASPECTS = [("overall", "Overall", "#7f7f7f"), ("accuracy", "Accuracy", "#3b7dd8"), ("scalability", "Scalability", "#e4572e"),
           ("stability", "Stability", "#e8b923"), ("usability", "Usability", "#3f9b4a")]

# engine strings in the results -> inventory tool names (first match wins)
ENGINE_TOOL = [
    (r"^genomeblocks|^Atlas\.|^Expression\.|^Genes\.|^Pairs\.|^Architecture\.|^per-track overlap_any", "genomeblocks"),
    (r"bedtools|^shell$|awk", "bedtools"),
    (r"pyranges 1\.|pyranges \(1", "pyranges (1.x)"),
    (r"^pyranges", "pyranges (0.x)"),
    (r"bioframe", "bioframe"),
    (r"intervaltree", "intervaltree"),
    (r"cgranges", "cgranges"),
    (r"^ncls", "ncls"),
    (r"superintervals", "superintervals"),
    (r"AIList|ailist", "ailist (Python)"),
    (r"polars-bio", "polars-bio"),
    (r"BiocPy genomicranges", "BiocPy genomicranges / iranges"),
    (r"^GenomicRanges", "GenomicRanges + IRanges"),
    (r"^plyranges", "plyranges"),
    (r"^valr", "valr"),
    (r"data\.table", "data.table"),
    (r"BEDOPS", "BEDOPS"),
    (r"txdbmaker|GenomicFeatures", "GenomicFeatures + txdbmaker"),
    (r"rtracklayer \+ EnrichedHeatmap|EnrichedHeatmap", "EnrichedHeatmap"),
    (r"rtracklayer", "rtracklayer"),
    (r"ChIPseeker", "ChIPseeker"),
    (r"ChIPpeakAnno", "ChIPpeakAnno"),
    (r"HOMER", "HOMER"),
    (r"UROPA", "UROPA"),
    (r"genomation", "genomation"),
    (r"bigWigAverageOverBed|UCSC", "UCSC kent utilities"),
    (r"pybigtools", "pybigtools"),
    (r"pyBigWig|pybigwig", "pyBigWig"),
    (r"deepTools", "deepTools"),
    (r"FIMO|MEME", "MEME suite"),
    (r"^MOODS|^moods", "MOODS"),
    (r"lightmotif", "lightmotif"),
    (r"Biopython|biopython", "Biopython"),
    (r"universalmotif", "universalmotif"),
    (r"Biostrings", "Biostrings"),
    (r"motifmatchr", "motifmatchr"),
    (r"GIGGLE", "GIGGLE"),
    (r"LOLA", "LOLA"),
    (r"regioneR", "regioneR"),
    (r"GenomicInteractions|InteractionSet", "GenomicInteractions"),
    (r"^igraph", "igraph"),
    (r"graph-tool", "graph-tool"),
    (r"networkx", "networkx"),
    (r"^scipy", "SciPy (csgraph)"),
    (r"cooler", "cooler"),
    (r"(?i)pytximport", "pytximport"),             # before tximport, which it contains
    (r"tximport", "tximport"),
    (r"edgeR", "edgeR"),
    (r"(?i)pydeseq2", "PyDESeq2"),                 # before DESeq2, which it contains
    (r"InMoose", "InMoose"),
    (r"DESeq2", "DESeq2"),
    (r"^limma", "limma"),
    (r"^sva ", "sva / ComBat-seq"),
    (r"(?i)subread|featureCounts", "Rsubread featureCounts"),
    (r"HTSeq", "HTSeq"),
    (r"GenomicAlignments", "GenomicAlignments"),
    (r"csaw", "csaw"),
    (r"megadepth", "megadepth"),
    (r"WiggleTools", "WiggleTools"),
    (r"(?i)pychromvar", "pychromVAR"),             # before chromVAR, which it contains
    (r"chromVAR", "chromVAR"),
    (r"(?i)intervene", "Intervene"),
    (r"hictk", "hictk"),
    (r"(?i)^mustache", "Mustache"),
    (r"(?i)juicer", "Juicer / Juicer Tools"),
    (r"^FitHiChIP", "FitHiChIP"),
    (r"^MAPS\b", "MAPS"),
    (r"HiCExplorer", "HiCExplorer"),
    (r"cooltools", "cooltools"),
    (r"rnanorm", "RNAnorm"),
    (r"pysam", "pysam"),
    (r"pyfaidx", "pyfaidx"),
]
PART_THEME = {"graph": "graph3d", "hic": "graph3d", "hic_loops": "graph3d", "signal_bigwig": "signal", "cooccupancy_bootstrap": "cooccupancy"}       # landscape parts whose name is not their theme
SKIP_ENGINES = re.compile(r"naive|numpy sliding|pandas\.read_csv|pandas chunked|R read\.delim|sparse row-sum only|warm index|"
                          r"zoom|pure-Python|threads \(same|index build|query \(warm\)|^memory$|^python$|^polars$|^pandas$|"
                          r"re-striped|exact interval Fisher|stats\(nBins\)|\(reference\)$")


def tool_of(engine: str):
    if engine is None or SKIP_ENGINES.search(engine):
        return None
    for pat, name in ENGINE_TOOL:
        if re.search(pat, engine):
            return name
    return None


def load(name):
    p = RES / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


# ── unified benchmark rows ───────────────────────────────────────────────────

def _acc_from(r):
    if r.get("accuracy") is not None:
        return float(r["accuracy"])
    if r.get("agrees") is not None:
        return 1.0 if r["agrees"] else 0.0
    if r.get("max_abs_diff") is not None:
        return 1.0 if r["max_abs_diff"] < 1e-3 else max(0.0, float(r.get("pearson_r", 0)))
    if r.get("ref_hits") is not None and r.get("total_hits") is not None:
        return 1.0 if r["ref_hits"] == r["total_hits"] else max(0.0, float(r.get("pearson_r", 0)))
    return None


def collect():
    """Every (theme, task, tool, engine, n, seconds, accuracy, detail, source) across the results."""
    rows, stab, skipped = [], [], []

    def add(theme, task, engine, n, seconds, acc, detail="", source="", **kw):
        # genomeblocks run through one of its backends is genomeblocks' answer (its time is the time of its fastest
        # backend), not the backend library's: the library is credited only where it was run on its own
        t = "genomeblocks" if "through genomeblocks' backend" in str(engine) else tool_of(engine)
        if t is None or seconds is None:
            return
        rows.append(dict(theme=theme, task=task, tool=t, engine=engine, n=n, seconds=float(seconds), accuracy=acc,
                         detail=detail, source=source, **kw))

    d = load("landscape")
    if d:
        for r in d["rows"]:
            part = r["part"]
            if r.get("error") and r.get("seconds") is None:
                if part != "stability" and r.get("engine") not in (None, "—"):
                    skipped.append(dict(theme=PART_THEME.get(part, part), task=r.get("task"), engine=r["engine"],
                                        tool=tool_of(r["engine"]) or r["engine"], error=r["error"]))
                continue
            if part == "stability":
                t = tool_of(r["engine"])
                if t:
                    stab.append(dict(theme=r["theme"], task=r["task"], tool=t, engine=r["engine"], variant=r["variant"], stable=bool(r["stable"]),
                                     error=r.get("error")))
                continue
            theme = PART_THEME.get(part, part)
            add(theme, r["task"], r["engine"], r.get("n"), r.get("seconds"), r.get("accuracy"), r.get("detail", ""), "landscape")
    d = load("loci")
    if d:
        for r in d["rows"]:
            task = {"intersect": "intersect", "merge": "merge", "latency": "point_lookup", "make": "read_bed"}.get(r["op"])
            if task is None or r["engine"] == "genomeblocks (cgranges, warm index)":
                continue
            add("intervals", task, r["engine"], r["n"], r["seconds"], _acc_from(r), "", "loci")
    d = load("backends")
    if d:
        for r in d["rows"]:
            if r.get("seconds") is None:
                continue
            p, op = r["part"], r["op"]
            eng = r["engine"] if r["engine"] != "genomeblocks" else "genomeblocks"
            via = f"{eng} (through genomeblocks' backend)" if eng != "genomeblocks" else "genomeblocks"
            if p == "intervals" and op == "A & B":
                add("intervals", "overlap_pairs", via, r["n"], r["seconds"], _acc_from(r), "", "backends")
            elif p == "intervals" and op.startswith("nearest"):
                add("intervals", "nearest", via, r["n"], r["seconds"], _acc_from(r), "", "backends")
            elif p == "bigwig" and eng != "genomeblocks":
                add("signal", "region_matrix", via, r["n"], r["seconds"], _acc_from(r), "", "backends")
            elif p == "motifs":
                add("motifs", "scan_logodds", via, r["n"], r["seconds"], _acc_from(r), "", "backends")
            elif p == "fasta":
                add("intervals", "getfasta", "genomeblocks (indexed FASTA reader)" if eng == "genomeblocks" else via, r["n"], r["seconds"], _acc_from(r), "", "backends")
            elif p == "graph" and not op.startswith("build"):
                add("graph3d", op.split(" ")[0], via, r["n"], r["seconds"], _acc_from(r), "", "backends")
    d = load("signal")
    if d:
        for r in d["rows"]:
            if r["part"] == "engines":
                add("signal", "region_matrix", r["engine"], r["n_loci"], r["seconds"], _acc_from(r),
                    f"max |Δ| {r.get('max_abs_diff', float('nan')):.2g}, Pearson {r.get('pearson_r', float('nan')):.4f}", "signal")
            elif r["part"] == "parallel" and r["workers"] == 8:
                add("signal", "region_matrix_16_tracks_8_workers", r["engine"], r["n_loci"], r["seconds"], None, "", "signal")
    d = load("atlas")
    if d:
        T = d["meta"].get("n_tracks", 500)
        for r in d["rows"]:
            if r["part"] == "query" and r.get("n_tracks") == T and "Atlas.search" not in (r.get("engine") or ""):
                add("cooccupancy", "overlap_counts_500_tracks", r["engine"], T, r["seconds"], None, "", "atlas")
    d = load("motifs")
    if d:
        for r in d["rows"]:
            if r["part"] == "engines":
                add("motifs", "scan_logodds", r["engine"], r["n_seqs"], r["seconds"], _acc_from(r),
                    f"{r['total_hits']:,} hits vs {r['ref_hits']:,}, Pearson {r['pearson_r']:.4f}", "motifs")
            elif r["part"] == "library" and r["n_seqs"] == 1000:
                add("motifs", "scan_library_1019", r["engine"], r["n_seqs"], r["seconds"], None, "", "motifs")
    d = load("pairs")
    if d:
        for r in d["rows"]:
            e = r["engine"]
            if "count_pairs_2d" in e or "cload" in e:
                add("graph3d", "bin_pairs_to_matrix", e, r["n_pairs"], r["seconds"], _acc_from(r), "", "pairs")
            elif "count_pairs (50 kb" in e:
                add("graph3d", "pairs_streaming", e, r["n_pairs"], r["seconds"], None, "", "pairs")
    d = load("genes")
    if d:
        for r in d["rows"]:
            s = r["step"]
            if s.startswith("Genes.make (GTF"):
                add("annotation", "parse_gtf", "genomeblocks Genes.make", r["n"], r["seconds"], None, "", "genes")
            elif s.startswith("annotations ("):
                eng = s[len("annotations ("):-1]
                # genomeblocks' labels on each of its interval backends: genomeblocks' answer, not the backend library's
                add("annotation", "labels", "genomeblocks labels" if eng == "genomeblocks" else f"genomeblocks labels ({eng} backend)", r["n"], r["seconds"], None, "", "genes")
            elif s == "nearest_genes":
                add("annotation", "nearest_tss", "genomeblocks nearest_genes", r["n"], r["seconds"], None, "", "genes")
            elif s.startswith("select_isoforms"):
                add("annotation", "isoform_selection", "genomeblocks " + s, r["n"], r["seconds"], None, "", "genes")
    d = load("architecture")
    if d:
        for r in d["rows"]:
            s = r["step"]
            if r.get("engine"):
                continue                                   # the graph engines come from backends.json
            task = {"make": "contact_graph", "add_mcool": "read_cool_mcool", "normalize": "oe_normalisation",
                    "annotate": "anchor_gene_annotation", "prime_hubs": "hub_mining"}.get(s.split(" ")[0])
            if task:
                add("graph3d", task, "genomeblocks Architecture." + s.split(" ")[0], r.get("n_edges") or r.get("n_vertices"), r["seconds"], None, "", "architecture")
    d = load("expression")
    if d:
        big = max(d["rows"], key=lambda r: (r.get("n_genes") or 0) * (r.get("n_samples") or 0))
        key = (big.get("n_genes"), big.get("n_samples"))
        for r in d["rows"]:
            if (r.get("n_genes"), r.get("n_samples")) != key or r.get("seconds") is None:
                continue
            task = {"read": "read_salmon", "to_genes": "tx_to_gene", "normalise": r["step"].split(" ")[0], "featurecounts": "read_featurecounts"}.get(r["part"])
            if task is None or (r["part"] == "to_genes" and "length_scaled_tpm" not in r["step"]) or r["step"].startswith("read + align"):
                continue
            add("expression", task, r["engine"], r.get("n_genes"), r["seconds"], _acc_from(r), "", "expression", n_samples=r.get("n_samples"))
    # the loop callers genomeblocks has ports of, each against the original on the same input (the benches' synthetic
    # runs, both timed there): reference = the original's output, accuracy = the share of its loops the port calls
    for name, task, ref_impl, gb_impl, n_key, shared_key, extra in (
            ("loops", "loops_fithichip", "FitHiChIP", "genomeblocks", "loops_ref", "loops_shared", "FitHiChIP 11.0 peak-to-all"),
            ("maps", "loops_maps", "MAPS", "genomeblocks", "loops_ref", "loops_shared", "MAPS (MAPS.py + regression; feather's file split not counted)"),
            ("mustache", "loops_mustache", "Mustache", "genomeblocks", "loops_ref", "loops_shared", "Mustache 1.3.3, 4 processes each"),
            ("hiccups", "loops_hiccups", "juicer_tools 2.13.06 hiccups --cpu", "genomeblocks (Juicer's SCALE + expected)", "merged_loops_juicer",
             "merged_loops_shared", "Juicer HiCCUPS, 8 threads each, the port on Juicer's SCALE vectors")):
        d = load(name)
        if not d:
            continue
        syn = [r for r in d["rows"] if r.get("part") == "synthetic"]
        agr = next((r for r in syn if r.get("impl") == "agreement"), None)
        ref = next((r for r in syn if r.get("impl") == ref_impl), None)
        port = next((r for r in syn if r.get("impl") == gb_impl), None)
        if not (agr and ref and port):
            continue
        n = int(agr[n_key])
        acc = agr[shared_key] / max(n, 1)
        orig = ref_impl if ref_impl != "Mustache" else "Mustache 1.3.3"
        add("graph3d", task, orig, n, ref["seconds"], 1.0, f"reference: the original's {n:,} loops ({extra})", name)
        add("graph3d", task, f"genomeblocks {name}()" if name != "loops" else "genomeblocks fithichip()", n, port["seconds"], acc,
            f"{agr[shared_key]:,} of the original's {n:,} loops, same input ({extra})", name)
    d = load("shortrange")
    if d:
        chk = next((r for r in d["rows"] if r["part"] == "check"), {})
        ok = 1.0 if chk.get("ends") and chk.get("bedgraph") else None
        for r in d["rows"]:
            if r["part"] == "total":
                add("signal", "hichip_shortrange_track", "genomeblocks hichip" if r["impl"] == "genomeblocks" else "awk | sort | bedtools genomecov (shell)",
                    chk.get("n_ends"), r["seconds"], ok, "", "shortrange")
    return rows, stab, skipped


# ── scoring ──────────────────────────────────────────────────────────────────

def usability(tool):
    inst = set(tool.get("install") or [])
    items = {
        "packaged": bool(inst & {"pip", "bioconda", "conda-forge", "Bioconductor", "CRAN", "npm", "cargo"}),
        "docs site": (tool.get("docs") or {}).get("site"),
        "tutorial": (tool.get("docs") or {}).get("tutorial"),
        "API reference": (tool.get("docs") or {}).get("api_reference"),
        "tests": (tool.get("tests_ci") or {}).get("tests"),
        "CI": (tool.get("tests_ci") or {}).get("ci"),
        "active (12 mo)": tool.get("active_12m"),
        "peer-reviewed": tool.get("peer_reviewed"),
        "open license": (bool(tool.get("license")) and not re.search(r"propriet|commercial|unknown", tool["license"], re.I)) if tool.get("license") is not None else None,
        "library API": tool.get("interface") in ("library", "mixed"),
    }
    known = {k: v for k, v in items.items() if v is not None}
    score = sum(1 for v in known.values() if v) / len(known) if known else None
    return score, items


def theme_coverage(tool, theme):
    keys = TOOLS["themes"][theme]["checklist"]
    caps = tool.get("capabilities") or {}
    return sum(1 for k in keys if caps.get(k) is True) / len(keys)


def score_tools(rows, stab):
    tools = {t["name"]: t for t in TOOLS["tools"]}
    by = defaultdict(list)
    for r in rows:
        by[(r["tool"], r["theme"])].append(r)
    # scalability: per (theme, task), the time of each tool at the largest size at which two or more tools ran
    tasks = defaultdict(list)
    for r in rows:
        tasks[(r["theme"], r["task"])].append(r)
    task_scores = {}
    for key, rs in tasks.items():
        sizes = sorted({r["n"] for r in rs if r["n"] is not None}, reverse=True)
        n_ref = next((n for n in sizes if len({r["tool"] for r in rs if r["n"] == n}) >= 2), sizes[0] if sizes else None)
        at = defaultdict(list)
        for r in rs:
            if r["n"] == n_ref:
                at[r["tool"]].append(r["seconds"])
        if not at:
            continue
        best = {t: min(v) for t, v in at.items()}
        t_min = min(best.values())
        task_scores[key] = {"n": n_ref, "t_min": t_min, "times": best,
                            "scores": {t: max(0.0, 1 - math.log10(v / t_min) / 3) for t, v in best.items()}}
    scores = {}
    for (tool, theme), rs in by.items():
        acc = [r["accuracy"] for r in rs if r["accuracy"] is not None]
        sc = [task_scores[(theme, r["task"])]["scores"][tool] for r in rs if (theme, r["task"]) in task_scores and tool in task_scores[(theme, r["task"])]["scores"]]
        sc = sorted(set(sc))
        st = [s["stable"] for s in stab if s["tool"] == tool and s["theme"] == theme]
        u = usability(tools[tool])[0] if tool in tools else None
        parts = {"accuracy": float(np.mean(acc)) if acc else None, "scalability": float(np.mean(sc)) if sc else None,
                 "stability": float(np.mean(st)) if st else None, "usability": u}
        have = [v for v in parts.values() if v is not None]
        # no overall without a measured accuracy: a tool whose tasks have no reference answer is not ranked
        parts["overall"] = float(np.exp(np.mean(np.log(np.maximum(have, 1e-3))))) if (have and parts["accuracy"] is not None) else None
        parts["n_tasks"] = len({r["task"] for r in rs})
        parts["n_stability_tests"] = len(st)
        scores[(tool, theme)] = parts
    return scores, task_scores


def tool_summary(scores):
    """One row per tool: the mean of its theme scores (usability is tool-level)."""
    out = defaultdict(lambda: defaultdict(list))
    for (tool, theme), p in scores.items():
        for k in ("overall", "accuracy", "scalability", "stability", "usability"):
            if p[k] is not None:
                out[tool][k].append(p[k])
        out[tool]["themes"].append(theme)
    return {t: {k: (float(np.mean(v)) if k != "themes" else v) for k, v in d.items()} for t, d in out.items()}



# ── the per-theme grade (2026-10-23) ─────────────────────────────────────────
JOBS = json.loads((HERE / "landscape" / "jobs.json").read_text())["themes"]
TIERS = [("best", 0.9), ("good", 0.6), ("moderate", 0.3), ("low", 0.0)]    # S relative to the theme's top S
MIN_KNOWN = 7                                    # core jobs whose status must be known before a tool is graded


def job_capable(tool_name, tool_rec, job, measured):
    """True / False / None (unknown) for one tool and one core job. genomeblocks: its status read from the source
    (jobs.json; 'partial' counts as no, like a recipe composed from other functions). Others: measured here, else the
    inventory's checklist keys for the job (any true -> yes; all false -> no; otherwise unknown)."""
    if tool_name == "genomeblocks":
        return job["genomeblocks"] == "yes"
    if measured & set(job["tasks"]):
        return True
    caps = (tool_rec or {}).get("capabilities") or {}
    vals = [caps.get(k) for k in job["checklist"]]
    if any(v is True for v in vals):
        return True
    if vals and all(v is False for v in vals):
        return False
    return None


def theme_grades(rows, task_scores):
    """{(tool, theme): grade} with coverage C (core jobs the tool does / core jobs whose status is known), quality Q
    (mean over its measured jobs of speed x accuracy on contested tasks: two or more tools timed at the scored size),
    S = C x Q and a tier; 'untested' when it does some of the jobs but none was measured on a contested task,
    'unable' when it does none."""
    tools = {t["name"]: t for t in TOOLS["tools"]}
    acc = defaultdict(list)
    measured = defaultdict(set)
    for r in rows:
        measured[(r["tool"], r["theme"])].add(r["task"])
        if r["accuracy"] is not None:
            acc[(r["tool"], r["theme"], r["task"])].append(r["accuracy"])
    out = {}
    names = set(tools) | {r["tool"] for r in rows}
    for theme, d in JOBS.items():
        core = d["core"]
        for name in names:
            m = measured.get((name, theme), set())
            cap = [job_capable(name, tools.get(name), j, m) for j in core]
            known = [c for c in cap if c is not None]
            n_yes = sum(1 for c in known if c)
            if n_yes == 0 and not m:
                if known:
                    out[(name, theme)] = {"state": "unable", "C": 0.0, "Q": None, "S": None, "unknown": len(cap) - len(known), "jobs": []}
                continue
            C = n_yes / len(known) if known else None
            qs, jobs_q = [], []
            for j, c in zip(core, cap):
                if not c:
                    continue
                q = []
                for task in j["tasks"]:
                    ts = task_scores.get((theme, task))
                    if ts and len(ts["scores"]) >= 2 and name in ts["scores"]:
                        a = acc.get((name, theme, task))
                        q.append(ts["scores"][name] * (float(np.mean(a)) if a else 1.0))
                if q:
                    qs.append(float(np.mean(q)))
                    jobs_q.append((j["id"], float(np.mean(q))))
            Q = float(np.mean(qs)) if qs else None
            S = C * Q if (C is not None and Q is not None and len(known) >= MIN_KNOWN) else None
            out[(name, theme)] = {"state": "untested", "C": C, "Q": Q, "S": S, "n_yes": n_yes, "n_known": len(known),
                                  "unknown": len(cap) - len(known), "jobs": jobs_q}
        top = max((v["S"] for (n, t), v in out.items() if t == theme and v["S"] is not None), default=None)
        for (n, t), v in out.items():
            if t == theme and v["S"] is not None:
                v["rel"] = v["S"] / top if top else 0.0
                v["state"] = next(tier for tier, lo in TIERS if v["rel"] >= lo)
    return out

# ── figure ───────────────────────────────────────────────────────────────────

FEATURED = [  # widely used tools shown off-the-shelf even without a benchmark row here
    "pybedtools", "GenomicRanges + IRanges", "plyranges", "bioframe", "pyranges (1.x)", "polars-bio", "BEDOPS",
    "ChIPseeker", "HOMER", "annotatr", "rGREAT", "GenomicFeatures + txdbmaker", "ChIPpeakAnno", "UROPA",
    "deepTools", "pyBigWig", "pybigtools", "EnrichedHeatmap", "genomation", "ROSE", "bwtool",
    "MEME suite", "MOODS", "lightmotif", "motifmatchr", "universalmotif", "GimmeMotifs", "chromVAR", "monaLisa", "tangermeme",
    "GIGGLE", "LOLA", "regioneR", "ChIP-Atlas Enrichment Analysis", "Cistrome DB Toolkit", "ReMapEnrich", "nullranges",
    "cooler", "cooltools", "pairtools", "HiCExplorer", "GenomicInteractions", "FitHiChIP", "hictk", "FAN-C", "cLoops2", "Juicer / Juicer Tools",
    "tximport", "DESeq2", "edgeR", "limma", "PyDESeq2", "pytximport", "RNAnorm", "InMoose",
    "pyGenomeTracks", "IGV desktop", "igv.js", "Gviz", "plotgardener", "CoolBox", "HiGlass",
]


def platform_label(t):
    p = t.get("platform") or "?"
    return {"C/C++": "C", "JavaScript": "JS", "mixed": "mixed"}.get(p, p)


def install_glyphs(t):
    inst = set(t.get("install") or [])
    g = []
    if inst & {"pip"}: g.append("pip")
    if inst & {"bioconda", "conda-forge"}: g.append("conda")
    if inst & {"Bioconductor"}: g.append("Bioc")
    if inst & {"CRAN"}: g.append("CRAN")
    if inst & {"npm"}: g.append("npm")
    if inst & {"web"} and not g: g.append("web")
    if not g: g.append("src")
    return " · ".join(g)


def figure_rows(scores, summary, grades=None):
    tools = {t["name"]: t for t in TOOLS["tools"]}
    bench = set(summary)
    shown = [n for n in tools if n in bench or n in FEATURED]
    order = list(THEMES)
    grades = grades or {}

    def primary(n):
        t = tools[n]
        th = [x for x in (t.get("themes") or []) if x in THEMES]
        graded = sorted((x for x in order if (grades.get((n, x)) or {}).get("S") is not None), key=lambda x: -grades[(n, x)]["S"])
        th = graded + th
        cov = {x: theme_coverage(t, x) for x in order}
        return th[0] if th else max(cov, key=cov.get)

    def n_themes(n):
        """themes where the tool does at least 3 of the 10 core jobs"""
        return sum(1 for x in order if (grades.get((n, x)) or {}).get("n_yes", 0) >= 3)

    groups = defaultdict(list)
    for n in shown:
        g = "frameworks" if n_themes(n) >= 3 else primary(n)
        groups[g].append(n)
    seq = ["frameworks"] + order
    out = []
    for g in seq:
        names = groups.get(g, [])
        th = g if g in THEMES else None
        names.sort(key=lambda n: (-((grades.get((n, th)) or {}).get("S") or -1) if th else -max([(grades.get((n, x)) or {}).get("S") or -1 for x in order]),
                                  -(summary.get(n, {}).get("overall") or -1), n.lower()))
        if names:
            out.append(("group", "Multi-theme frameworks" if g == "frameworks" else THEMES[g][0]))
            out.extend(("tool", n) for n in names)
    return out


GRADE_R = {"best": 0.16, "good": 0.125, "moderate": 0.09, "low": 0.055}       # dot radius per tier


def draw_figure(scores, summary, out_svg, out_files, grades=None):
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Circle
    plt.rcParams.update({"svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 8})
    tools = {t["name"]: t for t in TOOLS["tools"]}
    grades = grades or {}
    rows = figure_rows(scores, summary, grades)
    rh = 0.36
    cols = [("name", 2.9), ("platform", 0.75), ("interface", 0.7), ("install", 1.05), ("since", 0.5), ("active", 0.5)]
    cols += [(f"th:{k}", 0.42) for k in THEMES]
    cols += [("gap", 0.25)]
    cols += [(f"sc:{k}", 1.25) for k, _, _ in ASPECTS]
    xs, x = {}, 0.0
    for k, w in cols:
        xs[k] = (x, w); x += w
    W = x
    n_tool_rows = sum(1 for k, _ in rows if k == "tool")
    n_groups = sum(1 for k, _ in rows if k == "group")
    head = 2.3                                                  # header bars + rotated labels
    foot = 2.6                                                  # footnotes
    H = head + rh * n_tool_rows + rh * 1.3 * n_groups + foot
    fig = plt.figure(figsize=(W + 0.4, H + 0.4))
    ax = fig.add_axes([0.2 / (W + 0.4), 0.2 / (H + 0.4), W / (W + 0.4), H / (H + 0.4)])
    ax.set_xlim(0, W); ax.set_ylim(-H, 0); ax.axis("off")
    a0, a1 = xs["name"][0], xs["gap"][0]
    b0, b1 = xs["sc:overall"][0], W
    th0, th1 = xs["th:intervals"][0], xs["th:viewers"][0] + xs["th:viewers"][1]
    # panel header bars at the top, the theme bar under the first, the rotated labels below both
    bar_y = -0.36
    for (x0, x1, lab, letter) in ((a0, a1, "Tool", "a"), (b0, b1, "Summary · aggregated scores per theme experiment", "b")):
        ax.add_patch(Rectangle((x0, bar_y - 0.3), x1 - x0 - 0.05, 0.3, color="#555", lw=0))
        ax.text(x0 + 0.08, bar_y - 0.15, letter, color="white", fontweight="bold", va="center", fontsize=11)
        ax.text((x0 + x1) / 2, bar_y - 0.15, lab, color="white", va="center", ha="center", fontsize=9)
    ax.add_patch(Rectangle((th0, bar_y - 0.72), th1 - th0 - 0.05, 0.34, color="#d9d9d9", lw=0))
    ax.text((th0 + th1) / 2, bar_y - 0.55, "Per-theme grade (core jobs)", va="center", ha="center", fontsize=7.5)
    lab_y = -head + 0.15                                        # the rotated labels' baseline
    labels = {"platform": "Platform", "interface": "Interface", "install": "Install", "since": "Since", "active": "Active 12 mo"}
    for k, _ in cols:
        if k in labels:
            x0, w = xs[k]
            ax.text(x0 + w / 2, lab_y, labels[k], rotation=45, ha="left", va="bottom", fontsize=8, rotation_mode="anchor")
    for k, (lab, colr, short) in THEMES.items():
        x0, w = xs[f"th:{k}"]
        ax.text(x0 + w / 2, lab_y, short, rotation=45, ha="left", va="bottom", fontsize=8, color=colr, fontweight="bold", rotation_mode="anchor")
    for k, lab, colr in ASPECTS:
        x0, w = xs[f"sc:{k}"]
        ax.text(x0 + 0.1, lab_y, lab, rotation=45, ha="left", va="bottom", fontsize=8, rotation_mode="anchor")
    y = -head
    i_tool = 0
    for kind, name in rows:
        if kind == "group":
            y -= rh * 0.55
            ax.text(a0, y - rh * 0.15, name, fontsize=9, fontweight="bold", va="center")
            y -= rh * 0.75
            continue
        t = tools[name]
        s = summary.get(name)
        if i_tool % 2 == 0:
            ax.add_patch(Rectangle((a0, y - rh), W - a0, rh, color="#efefef", lw=0, zorder=0))
        i_tool += 1
        cy = y - rh / 2
        bold = name == "genomeblocks"
        ax.text(a0 + 0.05, cy, name, va="center", fontsize=8.5, fontweight="bold" if bold else "normal")
        ax.text(xs["platform"][0] + xs["platform"][1] / 2, cy, platform_label(t), va="center", ha="center")
        ax.text(xs["interface"][0] + xs["interface"][1] / 2, cy, (t.get("interface") or "?")[:7], va="center", ha="center")
        ax.text(xs["install"][0] + xs["install"][1] / 2, cy, install_glyphs(t), va="center", ha="center", fontsize=7)
        fp = (t.get("first_published") or {}).get("year")
        ax.text(xs["since"][0] + xs["since"][1] / 2, cy, str(fp) if fp else "–", va="center", ha="center")
        act = t.get("active_12m")
        ax.text(xs["active"][0] + xs["active"][1] / 2, cy, "✓" if act else ("–" if act is None else "✗"), va="center", ha="center",
                color="#3f9b4a" if act else "#999")
        for k, (lab, colr, short) in THEMES.items():
            x0, w = xs[f"th:{k}"]
            gr = grades.get((name, k))
            state = gr["state"] if gr else "unable"
            if state in GRADE_R:
                ax.add_patch(Circle((x0 + w / 2, cy), GRADE_R[state], color=colr, lw=0))
            elif state == "untested":
                ax.add_patch(Circle((x0 + w / 2, cy), 0.075, fill=False, ec=colr, lw=1.1))
            else:
                ax.add_patch(Circle((x0 + w / 2, cy), 0.035, color="#cfcfcf", lw=0))
        if s is None:
            u = usability(t)[0]
            x0, w = xs["sc:overall"]
            ax.text(x0 + 0.05, cy, "not benchmarked here", va="center", fontsize=7.5, color="#666", style="italic")
            x0, w = xs["sc:usability"]
            ax.add_patch(Rectangle((x0 + 0.05, y - rh + 0.06), w - 0.15, rh - 0.12, fill=False, ec="#d0d0d0", lw=0.5))
            if u is not None:
                ax.add_patch(Rectangle((x0 + 0.05, y - rh + 0.06), (w - 0.15) * u, rh - 0.12, color="#3f9b4a", alpha=0.35 + 0.65 * u, lw=0))
                ax.text(x0 + w - 0.12, cy, f"{u:.2f}", va="center", ha="right", fontsize=6.5, color="#222")
        else:
            for k, lab, colr in ASPECTS:
                x0, w = xs[f"sc:{k}"]
                v = s.get(k)
                ax.add_patch(Rectangle((x0 + 0.05, y - rh + 0.06), w - 0.15, rh - 0.12, fill=False, ec="#d0d0d0", lw=0.5))
                if v is None:
                    ax.text(x0 + w / 2, cy, "n/a", va="center", ha="center", fontsize=7, color="#999")
                    continue
                ax.add_patch(Rectangle((x0 + 0.05, y - rh + 0.06), (w - 0.15) * v, rh - 0.12, color=colr, alpha=0.3 + 0.7 * v, lw=0))
                ax.text(x0 + w - 0.12, cy, f"{v:.2f}", va="center", ha="right", fontsize=6.5, color="#222")
        y -= rh
    # footnotes
    shown = {n for k, n in rows if k == "tool"}
    rest = sorted(n for n in tools if n not in shown and n != "genomeblocks")
    note = ("Themes: the per-theme grade S = C x Q, where C is the share of the theme's 10 core jobs the tool does (of those whose status is known; at least 7 must be) and Q its mean "
            "speed x accuracy on the core jobs it was timed on against at least one other tool; dot size = S relative to the theme's best (best >= 0.9, good >= 0.6, moderate >= 0.3, low); "
            "hollow ring = does some core jobs but was not graded (not timed on a contested task, or too many unknowns); grey dot = does none of them. "
            "Scores: accuracy = agreement with each task's reference answer; scalability = 1 − log10(t / t_fastest)/3 at the largest size; "
            "stability = share of the perturbation tests (rows shuffled, chromosomes reordered, samples reordered) answered identically; "
            "usability = a ten-item checklist (packaged, docs, tutorial, API reference, tests, CI, active, peer-reviewed, open license, library API); "
            "overall = geometric mean, given only when accuracy was measured. Genomeblocks' usability lacks only the peer-reviewed item.")
    ax.text(a0, y - 0.3, "\n".join(textwrap.wrap(note, 230)), va="top", fontsize=7, color="#333", linespacing=1.4)
    ax.text(a0, y - 1.25, "\n".join(textwrap.wrap("Not shown (in the inventory, not benchmarked here): " + ", ".join(rest), 250)),
            va="top", fontsize=6.3, color="#555", linespacing=1.35)
    for p in out_files:
        fig.savefig(p, bbox_inches="tight", dpi=170)
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


# ── page ─────────────────────────────────────────────────────────────────────

def ftime(s):
    if s is None or not math.isfinite(s):
        return "–"
    if s < 1e-3:
        return f"{s * 1e6:.0f} µs"
    if s < 1:
        return f"{s * 1e3:.1f} ms" if s < 0.01 else f"{s * 1e3:.0f} ms"
    if s < 120:
        return f"{s:.1f} s" if s < 10 else f"{s:.0f} s"
    return f"{s / 60:.1f} min"


def fnum(x):
    if x is None:
        return "–"
    for d, u in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if x >= d:
            v = x / d
            t = f"{v:.0f}" if v >= 10 else f"{v:.1f}"
            return (t[:-2] if t.endswith(".0") else t) + u
    return f"{x:.0f}"


def fsc(v):
    return "–" if v is None else f"{v:.2f}"


def page(rows, stab, scores, task_scores, summary, skipped=(), grades=None):
    grades = grades or {}
    tools = {t["name"]: t for t in TOOLS["tools"]}
    md = []
    md.append("""---
title: Landscape
parent: Benchmarks
layout: default
nav_order: 1
permalink: /benchmarks/landscape/
---

# Landscape: genomeblocks against the field
{: .no_toc }

Every building block of genomeblocks does a job other tools also do. This page
lists those tools theme by theme, runs the ones that install here on the same
data and the same task, and scores each of them the way Saelens *et al.*
(*Nat. Biotechnol.* 2019) scored trajectory-inference methods: **accuracy**,
**scalability**, **stability** and **usability**, aggregated into one overall
score. The figure follows that paper's Fig. 2: tool properties and capability
icons on the left, the aggregated scores on the right.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## The figure
{: .sec-navy }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/landscape-figure.svg %}
</div><figcaption>
<strong>Tools by theme, graded per theme.</strong> Rows are tools, grouped
by the theme they serve (tools doing three or more core jobs in three or more
themes come first). Panel <strong>a</strong> gives each tool's platform,
interface, package index, first publication year, whether it saw a release or
commit in the last twelve months, and one dot per theme: its grade on that
theme's ten core jobs (see <a href="#the-grade">The grade</a>) — a filled dot
sized best / good / moderate / low, a hollow ring when the tool does some core
jobs but was not graded, a grey dot when it does none. Panel <strong>b</strong>
gives the four aspect scores and the overall score, averaged over the themes a
tool was benchmarked in; "not benchmarked here" rows show only usability.
</figcaption></figure>

Download: [SVG]({{ '/assets/images/landscape.svg' | relative_url }}) · [PNG]({{ '/assets/images/landscape.png' | relative_url }}) · [PDF]({{ '/assets/images/landscape.pdf' | relative_url }})

## The grade
{: .sec-navy #the-grade }

Each theme has ten **core jobs** (and ten auxiliary ones), chosen on
2026-10-23 from what the field's main tools are built around, before looking at
what genomeblocks does (`benchmarks/landscape/jobs.json`). A tool's grade in a
theme combines how much of the theme it covers with how well it does what it
was timed on:

- **C, coverage**: the share of the ten core jobs the tool does, among those
  whose status is known. A job counts when the tool has its own function for
  it (measured here, or documented in the inventory); a recipe composed from
  other functions does not count. genomeblocks' status comes from its source
  code, and its "partial" jobs count as no. A tool is graded only when at least
  seven of the ten statuses are known; unknown is never read as no.
- **Q, quality**: over the core jobs the tool was timed on *against at least one
  other tool* (a task only one tool ran gives no comparison and is left out),
  the mean of speed x accuracy, where speed is `1 − log10(t / t_fastest) / 3`
  and accuracy the task's agreement with its reference or simulated truth.
- **S = C x Q**, and the tier is S relative to the best S in the theme: best
  (≥ 0.9 of it), good (≥ 0.6), moderate (≥ 0.3), low.

The grade does not average across themes, and usability (a property of the
whole tool) is not part of it.

## The metrics
{: .sec-navy }

| Aspect | What is measured | How it is scored |
|---|---|---|
| **Accuracy** | Does the tool return the reference answer for the task? The reference is an independent tool or an exact computation: bedtools for interval sets and nearest distances, GenomicFeatures/TxDb for region classes, per-base bigWig means for signal matrices, an exact numpy sliding-window scan for log-odds motif hits and MEME FIMO for p-value hits (genomeblocks' own scan used to be the reference until 2026-10-23), the exact interval-level Fisher test for enrichment rankings, bedtools `pairtobed` for loop overlaps, R tximport / edgeR / DESeq2 for expression; and, for the jobs added on 2026-10-23, a simulated truth wherever one exists — simulated reads for BAM counting and coverage, simulated fold changes, dispersions and batch effects for the DE family, planted motifs, a Hi-C map with planted TADs, loops, compartments and bias. | 1 when the normalised answer is identical; otherwise the agreement — Pearson *r* for matrices and vectors, Spearman *ρ* for rankings, the fraction of rows agreeing for labels and distances — averaged over the tasks the tool ran in the theme. |
| **Scalability** | Wall time of the task at its largest size (1 M peaks, 100 k regions, 10 k windows × 200 bins, 1,000 windows × 100 motifs, 500 tracks, 50 k loops, 60 k genes × 24 samples), medians of repeated runs after a warm-up, one tool at a time. | `1 − log10(t / t_fastest) / 3`, floored at 0: each decade behind the fastest tool costs a third of the score; averaged over the tool's tasks. The slopes across sizes are in the tables. |
| **Stability** | The same data presented differently — rows shuffled, chromosomes in another order, samples or genes reordered — must give the identical answer (after undoing the permutation). | The fraction of perturbation tests passed. Tools that were not perturbed show "n/a". |
| **Usability** | A transparent ten-item checklist from the inventory, in the spirit of the paper's quality-control sheet: packaged on an index (pip / conda / Bioconductor / CRAN / npm), a documentation site, a tutorial, a function-level API reference, automated tests, continuous integration, a release or commit in the last twelve months, a peer-reviewed paper, an open license, a programmatic (library) interface. | The share of the known items a tool fulfils, equal weights. |
| **Overall** | | The geometric mean of the aspects a tool has, so a tool cannot buy back a failed aspect with the others. |

Timings of the main suite (bedtools, pyranges 0.x, bioframe, intervaltree,
cgranges, ncls, deepTools, pyBigWig, pybigtools, FIMO, MOODS, lightmotif,
Biopython, cooler, graph-tool, igraph, networkx, scipy, tximport, edgeR,
DESeq2) come from the [Benchmarks]({{ '/benchmarks/' | relative_url }}) page's
results; the tools added for this page (the Bioconductor stack, pyranges 1.x,
polars-bio, BiocPy, BEDOPS, HOMER, UROPA, EnrichedHeatmap, genomation, UCSC
tools, motifmatchr, universalmotif, Biostrings, LOLA, regioneR, GIGGLE,
GenomicInteractions, igraph, pytximport, PyDESeq2, RNAnorm) were run by
`benchmarks/bench_landscape.py` on the same data, as one SLURM job per theme on
cluster nodes with the same CPU model as the main suite's machine (Intel Xeon
Gold 6130, 6 CPUs per job; the parts are single-threaded except where noted).
""")
    # per-theme sections
    tasks_by_theme = defaultdict(set)
    for r in rows:
        tasks_by_theme[r["theme"]].add(r["task"])
    for theme, (label, colr, short) in THEMES.items():
        md.append(f"\n## {label}\n{{: .sec-navy #{theme} }}\n")
        inv = [t for t in TOOLS["tools"] if theme in (t.get("themes") or [])]
        gbt = tools.get("genomeblocks")
        if gbt:
            inv = [gbt] + [t for t in inv if t["name"] != "genomeblocks"]
        # the grade: the core jobs, then every graded / untested tool
        core = JOBS[theme]["core"]
        md.append("### Core jobs and grades\n{: .no_toc }\n")
        md.append("Core jobs: " + "; ".join(f"**{j['id']}** {escape(j['job'])}" + (" (timed)" if j["tasks"] else "") for j in core) + ".\n")
        gr = sorted(((n, v) for (n, th), v in grades.items() if th == theme and v["state"] != "unable"),
                    key=lambda x: (x[1]["S"] is None, -(x[1]["S"] or 0), x[0].lower()))
        if gr:
            md.append("| Tool | grade | S | C (core jobs done / known) | Q | timed core jobs (speed x accuracy) |\n|---|---|---|---|---|---|")
            for n, v in gr:
                if v["S"] is None and not v["jobs"]:
                    continue
                b = "**" if n == "genomeblocks" else ""
                cc = f"{v['n_yes']}/{v['n_known']}" + (f" ({v['unknown']} unknown)" if v["unknown"] else "")
                jq = ", ".join(f"{jid} {q:.2f}" for jid, q in v["jobs"]) or "–"
                md.append(f"| {b}{escape(n)}{b} | {v['state']} | {fsc(v['S'])} | {cc} | {fsc(v['Q'])} | {jq} |")
            md.append("")
        # scores table
        sc = [(t, scores[(t, theme)]) for (t, th) in scores if th == theme]
        sc.sort(key=lambda x: -(x[1]["overall"] or 0))
        if sc:
            md.append("### Scores\n{: .no_toc }\n")
            md.append("| Tool | Overall | Accuracy | Scalability | Stability | Usability | tasks | perturbation tests |\n|---|---|---|---|---|---|---|---|")
            for t, p in sc:
                b = "**" if t == "genomeblocks" else ""
                md.append(f"| {b}{escape(t)}{b} | {fsc(p['overall'])} | {fsc(p['accuracy'])} | {fsc(p['scalability'])} | {fsc(p['stability'])} | {fsc(p['usability'])} | {p['n_tasks']} | {p['n_stability_tests']} |")
            md.append("")
        # task tables
        for task in sorted(tasks_by_theme.get(theme, [])):
            rs = [r for r in rows if r["theme"] == theme and r["task"] == task]
            ts = task_scores.get((theme, task))
            n_ref = ts["n"] if ts else None
            md.append(f"### `{task}`\n{{: .no_toc }}\n")
            if n_ref is not None:
                md.append(f"Scored at n = {fnum(n_ref)}; fastest {ftime(ts['t_min'])}.\n")
            md.append("| Engine (as run) | Tool | n | time | accuracy | detail | source |\n|---|---|---|---|---|---|---|")
            for r in sorted(rs, key=lambda r: (-(r['n'] or 0), r["seconds"])):
                b = "**" if r["tool"] == "genomeblocks" else ""
                md.append(f"| {b}{escape(r['engine'])}{b} | {escape(r['tool'])} | {fnum(r['n'])} | {ftime(r['seconds'])} | {fsc(r['accuracy'])} | {escape(r.get('detail') or '')} | {r['source']} |")
            md.append("")
        sk = [x for x in skipped if x["theme"] == theme]
        if sk:
            md.append("### Not run or did not finish\n{: .no_toc }\n")
            md.append("| Engine | task | why |\n|---|---|---|")
            for x in sk:
                md.append(f"| {escape(x['engine'])} | {x['task']} | {escape(str(x['error'])[:300])} |")
            md.append("")
        st = [s for s in stab if s["theme"] == theme]
        if st:
            md.append("### Perturbation tests\n{: .no_toc }\n")
            md.append("| Tool | Engine | task | variant | identical answer |\n|---|---|---|---|---|")
            for s in sorted(st, key=lambda s: (s["tool"] != "genomeblocks", s["tool"], s["task"], s["variant"])):
                why = f" ({escape(str(s['error'])[:120])})" if (not s["stable"] and s.get("error")) else ""
                md.append(f"| {escape(s['tool'])} | {escape(s['engine'])} | {s['task']} | {s['variant']} | {'yes' if s['stable'] else '**no**' + why} |")
            md.append("")
        # inventory
        keys = TOOLS["themes"][theme]["checklist"]
        md.append(f"### Every tool in this theme\n{{: .no_toc }}\n\nCapability checklist: {', '.join(f'`{k}`' for k in keys)}.\n")
        md.append("| Tool | Platform | Interface | Install | Since | Last release | Active | Checklist | Usability | Notes |\n|---|---|---|---|---|---|---|---|---|---|")
        for t in sorted(inv, key=lambda t: (t["name"] != "genomeblocks", -theme_coverage(t, theme), t["name"].lower())):
            u = usability(t)[0]
            fp = (t.get("first_published") or {}).get("year") or "–"
            lr = t.get("last_release") or {}
            lrs = " ".join(x for x in (lr.get("version"), lr.get("date")) if x) or "–"
            act = t.get("active_12m")
            b = "**" if t["name"] == "genomeblocks" else ""
            md.append(f"| {b}{escape(t['name'])}{b} | {escape(t.get('platform') or '')} | {escape(t.get('interface') or '')} | {escape(install_glyphs(t))} | {fp} | {escape(lrs)} | "
                      f"{'yes' if act else ('?' if act is None else 'no')} | {theme_coverage(t, theme) * 100:.0f}% | {fsc(u)} | {escape(t.get('notes') or '')} |")
        md.append("")
    # usability checklist detail
    md.append("\n## Usability checklist\n{: .sec-navy #usability }\n")
    md.append("| Tool | " + " | ".join(usability(TOOLS['tools'][0])[1].keys()) + " | score |\n|---|" + "---|" * 11)
    shown = sorted({t for t, th in scores}, key=lambda n: (n != "genomeblocks", n.lower()))
    for n in shown:
        u, items = usability(tools[n])
        md.append(f"| {escape(n)} | " + " | ".join("yes" if v else ("?" if v is None else "no") for v in items.values()) + f" | {fsc(u)} |")
    md.append("""
## Reproduce
{: .sec-navy #reproduce }

```bash
cd benchmarks
python make_data.py                       # the synthetic hg38-shaped data
GB_LANDSCAPE_R=/path/to/Rscript \\
GB_LANDSCAPE_PY=/path/to/python \\
GB_LANDSCAPE_CLI=/path/to/cli-env/bin GB_LANDSCAPE_HIC=/path/to/hic-env \\
GIGGLE=/path/to/giggle HOMER_BIN=/path/to/homer/bin \\
MUSTACHE=/path/to/mustache JUICER_TOOLS=/path/to/juicer_tools_2.13.06.jar \\
python bench_landscape.py                 # results/landscape.json
python build_landscape.py                 # this page and the figure
```

One part per call (`GB_BENCH_PARTS=intervals`) keeps the other parts' rows;
with `GB_BENCH_OUT=landscape_<part>` each part writes its own file, so the
parts can run as separate cluster jobs, and `landscape/merge_results.py`
folds them back into `results/landscape.json`.

| Variable | What it points at |
|---|---|
| `GB_LANDSCAPE_R` | an Rscript with GenomicRanges, plyranges, valr, data.table, txdbmaker, ChIPseeker, ChIPpeakAnno, rtracklayer, EnrichedHeatmap, genomation, universalmotif, Biostrings, motifmatchr, TFBSTools, chromVAR, LOLA, regioneR, GenomicInteractions, InteractionSet, igraph, Rsubread, GenomicAlignments, csaw, edgeR, limma, DESeq2, apeglm, ashr, sva and tximport |
| `GB_LANDSCAPE_PY` | a Python with pyranges1, polars-bio, genomicranges, superintervals, ailist, pytximport, pydeseq2, inmoose, rnanorm, pychromvar and uropa, plus BEDOPS and `bigWigAverageOverBed` in its `bin/` |
| `GB_LANDSCAPE_CLI` | a `bin/` with featureCounts, htseq-count, megadepth, wiggletools, intervene and bedGraphToBigWig |
| `GB_LANDSCAPE_HIC` | an environment with cooler, cooltools, hictk and HiCExplorer (the `hic` part) |
| `GIGGLE`, `HOMER_BIN` | the GIGGLE binary (built from source) and HOMER's `bin/` |
| `MUSTACHE`, `JUICER_TOOLS`, `JAVA` | Mustache 1.3.3's command, juicer_tools 2.13.06 and the Java to run it (the `hic_loops` part; genomeblocks' ports run without them) |

The four loop callers genomeblocks has ports of are also compared with the
originals on real HiChIP and Hi-C samples by `bench_loops.py`, `bench_maps.py`,
`bench_mustache.py` and `bench_hiccups.py` (see [Benchmarks]({{ '/benchmarks/' | relative_url }})); their
synthetic runs feed the `loops_*` tasks here. The inventory in
`landscape/tools.json` was compiled on 2026-10-22 from each tool's repository
and package index; corrections are welcome as pull requests.
""")
    return "\n".join(md)


def build():
    rows, stab, skipped = collect()
    scores, task_scores = score_tools(rows, stab)
    summary = tool_summary(scores)
    inc = DOCS / "_includes" / "diagrams" / "landscape-figure.svg"
    img = DOCS / "assets" / "images"
    figs = HERE / "figures"
    figs.mkdir(exist_ok=True)
    grades = theme_grades(rows, task_scores)
    draw_figure(scores, summary, inc, [figs / "landscape.png", figs / "landscape.pdf", figs / "landscape.svg", img / "landscape.png", img / "landscape.pdf", img / "landscape.svg"],
                grades)
    # the inline include must not carry an XML prologue
    txt = inc.read_text()
    txt = txt[txt.index("<svg"):]
    inc.write_text(txt)
    (DOCS / "benchmarks" / "landscape.md").write_text(page(rows, stab, scores, task_scores, summary, skipped, grades))
    (RES / "landscape_scores.json").write_text(json.dumps({
        "grades": {f"{t} | {th}": g for (t, th), g in grades.items()},
        "scores": {f"{t} | {th}": p for (t, th), p in scores.items()},
        "tasks": {f"{th} | {task}": v for (th, task), v in task_scores.items()},
        "summary": summary}, indent=1, default=float))
    unmatched = sorted({r.get("engine") for d in (load(n) for n in ("landscape",)) if d for r in d["rows"]
                        if r.get("part") != "stability" and tool_of(r.get("engine")) is None and not (r.get("engine") and SKIP_ENGINES.search(r["engine"]))})
    if unmatched:
        print("[landscape] engines without an inventory tool:", unmatched)
    print(f"[landscape] {len(rows)} rows, {len(scores)} (tool, theme) scores, figure rows {len(figure_rows(scores, summary, grades))}")


if __name__ == "__main__":
    build()
