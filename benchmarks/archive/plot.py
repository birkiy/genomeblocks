#!/usr/bin/env python3
"""Render benchmark figures (PNG) from ``results/*.json``.

Style: thin marks, hairline grid, log axes where the spread is orders of
magnitude. Engine comparisons use emphasis (genomeblocks in blue, the
alternatives in gray); scaling plots give every engine a fixed hue that is
the same in every figure.
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

from common import FIGURES, RESULTS

# reference palette (validated categorical order) + chrome
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
GRAY = "#b5b4ad"

# colour follows the engine, in every figure
ENGINE_COLOR = {
    "genomeblocks (cgranges)": BLUE,
    "pyranges": ORANGE,
    "bioframe": AQUA,
    "bedtools (CLI)": YELLOW,
    "intervaltree": MAGENTA,
    "genomeblocks (pure-Python fallback)": RED,
    "genomeblocks (fallback, bounded scan)": VIOLET,
}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 9.5, "text.color": INK,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.linewidth": 0.8,
    "axes.titlesize": 11, "axes.titleweight": "semibold", "axes.titlelocation": "left",
    "axes.titlepad": 10, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
    "legend.frameon": False, "legend.fontsize": 9, "lines.linewidth": 2,
    "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
})


def load(name):
    return json.loads((RESULTS / f"{name}.json").read_text())


def fmt_time(s):
    s = float(f"{s:.6g}")
    if s < 1e-3:
        return f"{s * 1e6:.0f} µs" if s >= 1e-5 else f"{s * 1e6:.1f} µs"
    if s < 1:
        return f"{s * 1e3:.0f} ms" if s >= 0.01 else f"{s * 1e3:.1f} ms"
    if s < 120:
        return f"{s:.1f} s" if s < 10 else f"{s:.0f} s"
    return f"{s / 60:.1f} min"


def fmt_num(x):
    if x == 0:
        return "0"
    for div, suf in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if x >= div:
            v = x / div
            return f"{v:.0f}{suf}" if v >= 10 else f"{v:.1f}{suf}".replace(".0" + suf, suf)
    return f"{x:.0f}" if x >= 10 else f"{x:.1f}"


def save(fig, name):
    p = FIGURES / f"{name}.png"
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {p}")


def emphasis_bars(ax, labels, values, highlight, fmt, *, log=True, xlabel=""):
    """Horizontal bars, best at top; genomeblocks rows blue, others gray."""
    n = len(labels)
    y = list(range(n))[::-1]
    colors = [BLUE if h else GRAY for h in highlight]
    ax.barh(y, values, height=0.56, color=colors, zorder=3)
    ax.set_yticks(y, labels)
    ax.tick_params(axis="y", length=0)
    if log:
        ax.set_xscale("log")
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(xlabel)
    for yi, v in zip(y, values):
        ax.text(v * 1.06 if log else v, yi, " " + fmt(v), va="center",
                ha="left", fontsize=8.8, color=INK2)
    lo, hi = min(values), max(values)
    if log:
        ax.set_xlim(lo / 2.2, hi * 9)
    else:
        ax.set_xlim(0, hi * 1.25)
    for t in ax.get_yticklabels():
        if t.get_text().startswith("genomeblocks"):
            t.set_color(INK); t.set_fontweight("semibold")


def loglog_lines(ax, series, *, xlabel, ylabel, yfmt=fmt_time, label_ends=True):
    """series: [(label, xs, ys, color)]"""
    for label, xs, ys, color in series:
        ax.plot(xs, ys, color=color, label=label, zorder=3, marker="o", markersize=4.5,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_num(v)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: yfmt(v)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())


# ── Loci ─────────────────────────────────────────────────────────────────────

def fig_loci():
    d = load("loci")
    rows = d["rows"]

    def ser(op, eng):
        pts = sorted((r["n"], r["seconds"]) for r in rows if r["op"] == op and r["engine"] == eng)
        return [p[0] for p in pts], [p[1] for p in pts]

    # 1) intersect scaling: genomeblocks index backends | vs other tools
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for eng in ("genomeblocks (cgranges)", "genomeblocks (pure-Python fallback)",
                "genomeblocks (fallback, bounded scan)"):
        xs, ys = ser("intersect", eng)
        lab = {"genomeblocks (cgranges)": "cgranges (C index)",
               "genomeblocks (pure-Python fallback)": "pure-Python fallback (pip default)",
               "genomeblocks (fallback, bounded scan)": "fallback + bounded scan (proposed fix)"}[eng]
        loglog_lines(a1, [(lab, xs, ys, ENGINE_COLOR[eng])], xlabel="peaks per set (n)",
                     ylabel="time for A & B")
    a1.set_title("genomeblocks: which interval index")
    a1.legend(loc="lower right")
    for eng in ("genomeblocks (cgranges)", "pyranges", "bioframe", "bedtools (CLI)",
                "intervaltree"):
        xs, ys = ser("intersect", eng)
        loglog_lines(a2, [(eng, xs, ys, ENGINE_COLOR[eng])], xlabel="peaks per set (n)",
                     ylabel="")
    a2.set_title("vs. other tools")
    a2.legend(loc="upper left")
    fig.suptitle("Loci intersect (A & B), lower is better", x=0.06, ha="left",
                 fontsize=12.5, fontweight="bold", y=1.02)
    save(fig, "loci_intersect")

    # 2) single-query latency (emphasis bars)
    lat = [(r["engine"], r["seconds"]) for r in rows if r["op"] == "latency"]
    lat.sort(key=lambda t: t[1])
    fig, ax = plt.subplots(figsize=(8.5, 3.2))
    emphasis_bars(ax, [l for l, _ in lat], [v for _, v in lat],
                  [l.startswith("genomeblocks") for l, _ in lat], fmt_time,
                  xlabel="time per overlap query against 100k indexed peaks (log)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_time(v)))
    ax.set_title("One overlap lookup, lower is better")
    save(fig, "loci_latency")

    # 3) merge + make at 1M as grouped emphasis panels
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 2.9))
    for ax, op, title, n in ((a1, "merge", "sort + merge (2M intervals)", 1_000_000),
                             (a2, "make", "BED parse (1M lines)", 1_000_000)):
        pts = sorted([(r["engine"], r["seconds"]) for r in rows
                      if r["op"] == op and r["n"] == n], key=lambda t: t[1])
        emphasis_bars(ax, [p[0] for p in pts], [p[1] for p in pts],
                      [p[0].startswith("genomeblocks") for p in pts], fmt_time, log=False)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_time(v) if v else "0"))
        ax.set_title(title)
    save(fig, "loci_merge_make")


# ── signal ───────────────────────────────────────────────────────────────────

SIG_SHORT = {
    "genomeblocks · pybigtools (exact)": "genomeblocks signal()",
    "genomeblocks · pybigtools (zoom, exact=False)": "genomeblocks signal(exact=False)",
    "genomeblocks · pure-Python reader": "genomeblocks, pure-Python reader",
    "pybigtools values() + numpy binning": "pybigtools per-base + numpy",
    "pyBigWig values() + numpy binning": "pyBigWig per-base + numpy",
    "pyBigWig stats(nBins) exact": "pyBigWig stats(nBins=200)",
    "pyBigWig stats(nBins) zoom": "pyBigWig stats(nBins=200), zoom",
    "deepTools computeMatrix (-p 1)": "deepTools computeMatrix",
}


def fig_signal():
    rows = load("signal")["rows"]
    eng = sorted([r for r in rows if r["part"] == "engines"], key=lambda r: -r["rate"])
    fig, ax = plt.subplots(figsize=(9, 3.9))
    emphasis_bars(ax, [SIG_SHORT[r["engine"]] for r in eng], [r["rate"] for r in eng],
                  [r["engine"].startswith("genomeblocks") for r in eng], fmt_num,
                  xlabel="regions per second, 1 track, 1 core (log)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_num(v)))
    ax.set_title("bigWig → heatmap matrix (±3 kb, 200 bins), higher is better")
    save(fig, "signal_engines")

    sc = [r for r in rows if r["part"] == "scaling"]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for eng_, col in (("genomeblocks · pybigtools (exact)", BLUE),
                      ("pyBigWig values() + numpy binning", ORANGE),
                      ("genomeblocks · pure-Python reader", AQUA)):
        pts = sorted((r["n_loci"], r["rate"]) for r in sc if r["engine"] == eng_)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, marker="o",
                markersize=4.5, markeredgecolor=SURFACE, markeredgewidth=1.5, label=SIG_SHORT[eng_])
    ax.set_xscale("log"); ax.set_ylim(0)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_num(v)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_num(v)))
    ax.set_xlabel("loci (sampled from 100k peaks, genome order)"); ax.set_ylabel("regions / s")
    ax.set_title("Denser peak sets reuse inflated blocks")
    ax.legend(loc="upper left")
    save(fig, "signal_scaling")

    par = [r for r in rows if r["part"] == "parallel"]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    proc = sorted((r["workers"], r["rate"]) for r in par if r["engine"].startswith("genomeblocks"))
    thr = [(1, proc[0][1])] + sorted((r["workers"], r["rate"]) for r in par if r["engine"].startswith("threads"))
    dt = sorted((r["workers"], r["rate"]) for r in par if r["engine"].startswith("deepTools"))
    for pts, col, lab in ((proc, BLUE, "genomeblocks, processes + shared memory"),
                          (thr, ORANGE, "same work on threads"),
                          (dt, AQUA, "deepTools computeMatrix -p")):
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, marker="o", markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=1.5, label=lab)
    ax.set_xticks([1, 2, 3, 4]); ax.set_ylim(0)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_num(v)))
    ax.set_xlabel("workers (4-core machine)"); ax.set_ylabel("region-tracks / s")
    ax.set_title("16 tracks × 5,000 loci")
    ax.legend(loc="upper left")
    save(fig, "signal_parallel")


# ── Atlas ────────────────────────────────────────────────────────────────────

ATLAS_COLOR = {"Atlas.search (incl. Fisher + DataFrame)": BLUE,
               "GIGGLE search -s (CLI)": ORANGE,
               "per-track cgranges loop (prebuilt)": AQUA,
               "pyranges overlap(), per track": YELLOW,
               "bedtools intersect -C (CLI)": MAGENTA}


def fig_atlas():
    rows = load("atlas")["rows"]
    q = [r for r in rows if r["part"] == "query"]
    q += [{**r, "engine": "GIGGLE search -s (CLI)"} for r in rows
          if r["part"] == "giggle" and r.get("kind") == "query"]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for eng, col in ATLAS_COLOR.items():
        pts = sorted((r["n_tracks"], r["seconds"]) for r in q if r["engine"] == eng)
        if pts:
            loglog_lines(ax, [(eng, [p[0] for p in pts], [p[1] for p in pts], col)],
                         xlabel="tracks in the collection", ylabel="time per 20k-peak query")
    ax.set_title("Enrichment query vs. collection size, lower is better")
    ax.legend(loc="upper left")
    save(fig, "atlas_query")


# ── motifs ───────────────────────────────────────────────────────────────────

def fig_motifs():
    rows = load("motifs")["rows"]
    eng = sorted([r for r in rows if r["part"] == "engines"], key=lambda r: r["seconds"])
    fig, ax = plt.subplots(figsize=(9, 3.2))
    emphasis_bars(ax, [r["engine"] for r in eng], [r["seconds"] for r in eng],
                  [r["engine"].startswith("genomeblocks") for r in eng], fmt_time,
                  xlabel="time: 1,000 × 500 bp windows × 100 JASPAR motifs (log)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_time(v)))
    ax.set_title("Motif scanning, lower is better")
    save(fig, "motifs_engines")


# ── pairs / genes / architecture / import ────────────────────────────────────

def fig_pairs():
    rows = sorted(load("pairs")["rows"], key=lambda r: -r["rate"])
    fig, ax = plt.subplots(figsize=(9, 3.0))
    emphasis_bars(ax, [r["engine"] for r in rows], [r["rate"] for r in rows],
                  [r["engine"].startswith("genomeblocks") for r in rows], fmt_num,
                  xlabel="read pairs per second (log)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_num(v)))
    ax.set_title("Hi-C pair counting (5M pairs), higher is better")
    save(fig, "pairs")


def fig_architecture():
    rows = load("architecture")["rows"]
    fig, ax = plt.subplots(figsize=(9, 3.4))
    emphasis_bars(ax, [r["step"] for r in rows], [r["seconds"] for r in rows],
                  ["loop" not in r["step"] for r in rows], fmt_time,
                  xlabel="seconds (log)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_time(v)))
    ax.set_title("Architecture pipeline steps")
    save(fig, "architecture")


def fig_genes():
    rows = load("genes")["rows"]
    fig, ax = plt.subplots(figsize=(9, 3.6))
    lab = [f'{r["step"]}' + (f' · n={fmt_num(r["n"])}' if "annotations" in r["step"] or "nearest" in r["step"] else "")
           for r in rows]
    emphasis_bars(ax, lab, [r["seconds"] for r in rows],
                  ["fallback" not in r["step"] for r in rows], fmt_time, xlabel="seconds (log)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_time(v)))
    ax.set_title("Gene models")
    save(fig, "genes")


def fig_import():
    rows = [r for r in load("import")["rows"][1:] if not r["stmt"].startswith("breakdown")]
    fig, ax = plt.subplots(figsize=(8, 2.9))
    emphasis_bars(ax, [r["stmt"] for r in rows], [max(r["seconds"], 1e-4) for r in rows],
                  [r["stmt"] == "import genomeblocks" for r in rows], fmt_time, log=False,
                  xlabel="seconds above bare interpreter start-up")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_time(v) if v else "0"))
    ax.set_title("Import cost")
    save(fig, "import")


FIGS = {"loci": fig_loci, "signal": fig_signal, "atlas": fig_atlas, "motifs": fig_motifs,
        "pairs": fig_pairs, "genes": fig_genes, "architecture": fig_architecture,
        "import": fig_import}

if __name__ == "__main__":
    for name in (sys.argv[1:] or list(FIGS)):
        FIGS[name]()
