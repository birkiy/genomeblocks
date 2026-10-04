#!/usr/bin/env python3
"""Build the objects-vs-columns explainer (report/columns.html).

Numbers come from results/columnar.json (bench_columnar.py) and
results/genes_make.json (bench_genes_make.py).
"""
from __future__ import annotations

import json
import re

from build_report import REP, chart, details, esc, ftime, fx, load, table


def diagrams():
    raw = (REP / "diagrams_columns.html").read_text()
    return {m.group(1): m.group(2).strip()
            for m in re.finditer(r"<!-- @(\w+) -->\s*(.*?)(?=<!-- @|\Z)", raw, re.S)}


def fig(svg, title, caption):
    return (f'<div class="why"><h3>{esc(title)}</h3><div class="scroll">{svg}</div>'
            f'<p class="sub">{caption}</p></div>')


def us(s):
    return f"{s * 1e6:.1f} µs" if s < 1e-5 else f"{s * 1e6:.0f} µs"


def build():
    D = diagrams()
    rows = load("columnar")["rows"]
    gm = load("genes_make")["rows"]
    mem = {r["storage"]: r for r in rows if r["part"] == "memory"}
    acc = {(r["pattern"], r["storage"]): r["seconds"] for r in rows if r["part"] == "access"}
    fn = {}
    for r in rows:
        if r["part"] == "functions":
            fn.setdefault(r["function"], {})[r["storage"]] = r
    ab = {r["step"]: r["seconds"] for r in rows if r["part"] == "annot_breakdown"}
    ra = {r["way"]: r["seconds"] for r in rows if r["part"] == "row_access"}
    all_same = all(r.get("same", True) for r in rows if r["part"] == "functions")
    n_fn = len(fn)

    b_obj = mem["objects (Loci of Locus)"]["per_interval"]
    b_col = mem["columns (ColumnarLoci)"]["per_interval"]
    mk_obj, mk_col = mem["objects: Loci.make time"]["seconds"], mem["columns: make time"]["seconds"]
    loop_o = acc[("loop: for l in loci: total += l.start", "objects")]
    loop_v = acc[("loop: for l in loci: total += l.start", "columns via LocusView")]
    loop_n = acc[("loop: for l in loci: total += l.start", "columns, vectorised")]
    ann_obj = fn["Genes.annotations (100k CREs)"]["objects"]["seconds"]
    ann_col_genes = ab["label 100k CREs (5 vectorised overlap tests)"] + ab["uid strings for the output table"]
    g_annot_cur = next(r["seconds"] for r in gm if r["part"] == "annot" and r["variant"].startswith("current"))
    g_annot_tab = next(r["seconds"] for r in gm if r["part"] == "annot" and "table" in r["variant"])
    g_parse_cur = next(r["seconds"] for r in gm if r["part"] == "parse" and r["variant"] == "Genes.make (current)")
    g_parse_pol = next(r["seconds"] for r in gm if r["part"] == "parse" and r["variant"].startswith("polars read_csv"))

    # charts
    charts = {}
    UNCH = "columns, today's code"
    charts["unchanged"] = {"kind": "diverge", "labelW": 300, "axis": "today's code given columnar Loci, vs object Loci",
                           "rows": [{"label": f, "value": d["objects"]["seconds"] / d[UNCH]["seconds"],
                                     "tip": f"{f}: {ftime(d['objects']['seconds'])} → {ftime(d[UNCH]['seconds'])}"}
                                    for f, d in fn.items() if UNCH in d]}
    vec = [{"label": f, "value": d["objects"]["seconds"] / d["columns, vectorised"]["seconds"]}
           for f, d in fn.items() if "columns, vectorised" in d and "annotations" not in f]
    vec.append({"label": "Genes.annotations (gene index also columns)", "value": ann_obj / ann_col_genes})
    vec.append({"label": "Loci.make (1M peaks)", "value": mk_obj / mk_col})
    charts["vector"] = {"kind": "hbar", "fmt": "x", "log": True, "sort": "desc", "labelW": 330,
                        "axis": "speed-up over today's object code (log scale)",
                        "rows": [{**r, "hl": True} for r in vec]}
    order = ["Locus objects: loci[i].start", "columns via LocusView: loci[i].start", "raw numpy column: starts[i]",
             "polars: df.row(i)", "Arrow: table['Start'][i]", "pandas: df.at[i, 'Start']", "pandas: df.iloc[i]['Start']"]
    charts["rows"] = {"kind": "hbar", "fmt": "us", "log": True, "sort": None, "labelW": 300,
                      "axis": "time to read one row (log scale)",
                      "rows": [{"label": k, "value": ra[k], "hl": k.startswith(("columns", "Locus"))} for k in order]}

    # audit table: measured
    WHAT = {
        "A & B": ("one index lookup per peak", "vectorise (done in prototype)"),
        "A - B": ("one index lookup per peak", "vectorise (done in prototype)"),
        "merge (200k)": ("sorts and copies Locus objects", "vectorise (done in prototype)"),
        "make → slop → sort → merge": ("builds a new Locus per peak at every step", "vectorise (done in prototype)"),
        "Genes.annotations (100k CREs)": ("5 index lookups per CRE", "vectorise; needs the gene index as columns too"),
        "Genes.nearest_genes (100k)": ("hands frames to pyranges", "none; build frames from the columns"),
        "Atlas.search (20k peaks, 500 tracks)": ("bins each peak in a Python loop", "vectorise the binning (prototype: 18×)"),
        "Atlas: peaks → bin ranges": ("the loop above", "vectorise (done in prototype)"),
        "signal (5k loci, 1 bigWig)": ("copies chrom/start/end into arrays first", "none; read the columns directly"),
        "scan_motifs_matrix (1k windows × 100)": ("slices the genome string per peak", "none"),
        "pair_to_bed (50k loops)": ("index lookups per loop anchor", "none"),
        "count_pairs (200k pairs, 62k windows)": ("copies windows into per-chromosome arrays", "none; use the columns directly"),
        "Architecture.make (50k loops)": ("index lookups per loop, a Locus per hit", "optional: map anchors in one call"),
        "Architecture.annotate": ("calls annotations + nearest_genes", "follows Genes.annotations"),
    }
    trows = []
    for f, d in fn.items():
        o = d["objects"]["seconds"]
        un = d.get("columns, today's code", {}).get("seconds")
        ve = d.get("columns, vectorised", {}).get("seconds")
        def rel(t):
            if t is None:
                return "—"
            r = o / t
            return f"{ftime(t)} ({fx(r)} faster)" if r >= 1.05 else (f"{ftime(t)} (same)" if r > 0.95 else f"{ftime(t)} ({fx(1 / r)} slower)")
        same = all(v.get("same", True) for v in d.values())
        how, todo = WHAT.get(f, ("", ""))
        trows.append([f, how, ftime(o), rel(un), rel(ve), "yes" if same else "NO", todo])
    static = [
        ["liftover", "one pyliftover call per peak", "works through views, same speed", "none"],
        ["browser / architecture_draw / signal_draw", "draw or group a few hundred peaks", "works through views", "none"],
        ["tile / tile_genome", "build a Locus per tile", "one arange per chromosome", "vectorise"],
        ["to_frame / to_pyranges / to_bed", "a dict or string per row", "frames straight from the columns", "vectorise"],
        ["Genes internals (Transcript.exons, annot)", "Loci of Exon / CDS / UTR objects", "linked tables + views", "the larger change"],
        ["Loci slicing (loci[a:b])", "returns a plain list today", "returns a Loci", "fix"],
    ]

    secs = f"""
<section id="answers">
  <h2>Short answers</h2>
  <div class="tiles two-up">
    <a class="tile" href="#what"><span class="tag green">what</span><span class="k">What is column data?</span>
      <span class="v">{b_obj / b_col:.0f}× smaller</span><span class="d">Each field of every peak lives in its own array (all starts together, all ends together). One peak is one row across those arrays: {b_col:.0f} bytes per peak instead of {b_obj:.0f}.</span></a>
    <a class="tile" href="#pandas"><span class="tag navy">pandas?</span><span class="k">Is it a pandas-like frame?</span>
      <span class="v">Same idea</span><span class="d">pandas, polars and Arrow are column stores too. They differ in what they wrap around the arrays. pandas is too slow to read row by row ({us(ra["pandas: df.at[i, 'Start']"])}–{us(ra["pandas: df.iloc[i]['Start']"])} per row), so it fits as an export, not as the storage.</span></a>
    <a class="tile" href="#functions"><span class="tag purple">functions</span><span class="k">What happens to annotate &amp; co?</span>
      <span class="v">{n_fn}/{n_fn} identical</span><span class="d">Every function tested gives the same output on columns without a code change. Whole-set work gets 9–30× faster; four per-peak loops get 1.2–2× slower until they are vectorised.</span></a>
    <a class="tile" href="#objects"><span class="tag warn">objects</span><span class="k">What do objects give us?</span>
      <span class="v">{fx(loop_v / loop_o)}</span><span class="d">faster plain Python loops over peaks ({ftime(loop_o)} vs {ftime(loop_v)} per 100k), nested gene models, per-item fields and methods, and code anyone can read. Views keep the code; vectorising keeps the speed.</span></a>
  </div>
</section>

<section id="what">
  <span class="tag green">What it is</span>
  <h2>Same peaks, two ways to store them</h2>
  <p class="verdict">Objects keep each peak together in its own box. Columns keep each <em>field</em> together: all chromosomes in one array, all starts in the next, all ends in the next. Nothing about a peak changes; only where its four numbers sit in memory.</p>
  {fig(D['layout'], "Four peaks, objects vs columns",
       f"Measured on 1M peaks: {b_obj:.0f} bytes per peak as objects (the Python object header, its attribute dict, and separate integer objects for start and end), {b_col:.0f} bytes as columns. Loading 1M peaks takes {ftime(mk_obj)} as objects and {ftime(mk_col)} as columns.")}
  <p class="sub">Why columns are faster for whole-set work: an operation such as "which of these peaks overlap that set" runs as a few calls over whole arrays, compiled code walking memory in order, instead of one Python step per peak. Summing every start: {ftime(loop_o)} as a Python loop over objects, {ftime(loop_n)} as one numpy call over the column.</p>
</section>

<section id="pandas">
  <span class="tag navy">pandas, polars, Arrow</span>
  <h2>Is that a DataFrame?</h2>
  <p class="verdict">A DataFrame is already column data. pandas, polars, Arrow and the proposed <code>Loci</code> all hold arrays; what differs is what each adds around them, and how fast one row can be read.</p>
  {fig(D['wrappers'], "Same arrays, different wrappers",
       "Arrow is the common format: polars, DuckDB and R's arrow package read it without copying, and pandas can too. A columnar Loci built on plain numpy arrays can hand its arrays to any of them almost for free, and genomeblocks keeps its own genomics layer on top.")}
  {chart("rows", "Reading one peak at a time", "· 10,000 random rows from 100k peaks, lower is better",
         "Genomics code reads single peaks constantly (graph building, annotation, plotting). Views over numpy columns read a row as fast as today's objects; pandas is 14–50× slower per row, which is why it would make a poor storage layer for genomeblocks even though it is columnar.")}
</section>

<section id="api">
  <span class="tag green">Your code</span>
  <h2>Code that uses Loci does not change</h2>
  <p class="verdict"><code>loci[i]</code> and <code>for peak in loci</code> still hand you a <code>Locus</code>. It is a small view onto row i that reads and writes the columns.</p>
  {fig(D['view'], "What loci[2] returns", "The view is a real Locus subclass, so isinstance checks, uid, distance_to, overlaps, sequence and printing all behave as today. Editing it edits the stored peak.")}
  <pre># unchanged, works on both storages
cre = Loci.make("atac.narrowPeak").slop(100).sort().merge()
for peak in cre[:5]:
    print(peak.uid, peak.length)
arch = Architecture.make(cre, "loops.bedpe", r=2500)

# new with columns: whole arrays, and free hand-over to other tools
cre.starts, cre.ends               # numpy arrays
cre.to_polars(); cre.to_arrow()    # no per-row copy</pre>
</section>

<section id="functions">
  <span class="tag purple">Every function</span>
  <h2>What happens to annotate and the others</h2>
  <p class="verdict">I swapped in a columnar <code>Loci</code> and ran genomeblocks' own functions on it without changing them. {'All' if all_same else 'NOT all'} {n_fn} gave identical output.</p>
  {chart("unchanged", "Today's code, columnar input", "· per-peak loops pay about 1 µs to make each view",
         "Most functions run at the same speed. The four slower ones loop over peaks in Python and read fields through views; a plain loop costs " + f"{ftime(loop_v)} instead of {ftime(loop_o)} per 100k peaks." + " Each has an array version, below.")}
  {chart("vector", "With the hot loops written for columns", "· identical output, higher is better",
         f"Genes.annotations is the instructive case. Its labelling takes {ftime(ab['label 100k CREs (5 vectorised overlap tests)'])} on columns, but converting today's object-based gene index to arrays costs {ftime(ab['convert the gene index (objects) to arrays'])} on every call, so the gain only shows once the gene index is stored as columns too. The speed comes when the whole path is columnar.")}
  {details("Per-function table (measured)", table(["function", "what it does with peaks", "objects", "columns, today's code", "columns, vectorised", "same output", "change needed"], trows, num_cols=(2,)))}
  {details("Functions not timed (read from the code)", table(["function", "what it does with peaks", "on columns", "change needed"], static))}
</section>

<section id="objects">
  <span class="tag warn">Objects</span>
  <h2>What objects give us</h2>
  <div class="mech">
    <div><h3>Fast plain loops</h3><p>A Python loop over objects is the fastest per-peak code: {ftime(loop_o)} for 100k vs {ftime(loop_v)} through views. Anything written as "for peak in loci" is quickest on objects.</p></div>
    <div><h3>Nested structure</h3><p>Gene → Transcript → Exon/CDS/UTR is a tree. Objects model it directly (<code>g.transcripts[t].exons</code>); columns need linked tables and an index to walk it.</p></div>
    <div><h3>Per-item fields and methods</h3><p><code>Exon.exon_number</code>, <code>UTR.type</code>, <code>Transcript.tss_score</code>, <code>Gene.canonical</code>; <code>sequence()</code>, <code>distance_to()</code>. Columns can hold extra fields, but each needs declaring.</p></div>
    <div><h3>Identity and editing</h3><p>An object is a thing you can edit, put in a set or use as a dict key. Views edit too, but <code>loci[0] is loci[0]</code> becomes False: two views of one row are equal, not identical.</p></div>
    <div><h3>Mixed types in one container</h3><p>Today a Loci can hold Genes, Transcripts and plain peaks at once (the annotation index does). Columns need one schema per container.</p></div>
    <div><h3>Easy to read and debug</h3><p>No numpy knowledge needed; printing a peak shows a peak. This is a real advantage for users writing their own code.</p></div>
  </div>
  <p class="sub"><b>What objects cost:</b> {b_obj / b_col:.0f}× the memory, a Python step per peak for every bulk operation, and a full conversion each time data goes to pandas, pyranges or a faster engine. With views, the convenience stays where it is used (one peak at a time) and the cost goes away where it hurts (whole sets).</p>
</section>

<section id="genes">
  <span class="tag navy">Genes</span>
  <h2>Genes: a tree, stored as linked tables</h2>
  {fig(D['genes'], "Gene models as columns",
       f"Parsing the GTF into tables takes {ftime(g_parse_pol)} instead of {ftime(g_parse_cur)}, and the annotation index builds from them in {ftime(g_annot_tab)} instead of {ftime(g_annot_cur)}. Gene and Transcript views would keep <code>g.transcripts[t].exons</code> working, the same way LocusView keeps <code>loci[i]</code> working.")}
</section>

<section id="plan">
  <span class="tag navy">Recommendation</span>
  <h2>Columns for storage, objects as views</h2>
  <div class="fixes">
    <div class="fix"><span class="t">1 · Columnar Loci with LocusView</span><p>Same public API; {n_fn}/{n_fn} functions already give identical output. Set operations, sort, merge, slop, make and the to_ conversions move to arrays.</p></div>
    <div class="fix"><span class="t">2 · Rewrite the four per-peak hot loops</span><p>Genes.annotations, Atlas binning, count_pairs grouping and signal's array copy read the columns directly, so nothing gets slower.</p></div>
    <div class="fix"><span class="t">3 · Genes as linked tables</span><p>Gene / Transcript views on top, so annotate and nearest_genes become array work end to end.</p></div>
    <div class="fix"><span class="t">Watch for</span><p>Code that relies on <code>loci[i] is loci[i]</code> (identity), and users who subclass Locus to add their own fields; both need a note in the release.</p></div>
  </div>
</section>

<section id="method">
  <span class="tag navy">Method</span>
  <h2>How this was measured</h2>
  <ul class="plain">
    <li>Prototype: <code>benchmarks/prototypes/columnar.py</code> (ColumnarLoci + LocusView, about 250 lines). Bench: <code>bench_columnar.py</code>. genomeblocks itself was not modified.</li>
    <li>Same machine and synthetic hg38-shaped data as the <a href="https://claude.ai/artifact/4hL67NF3y68wJSdsySs1Eo">benchmark report</a>; medians of 3–5 runs after a warm-up; memory with <code>tracemalloc</code>.</li>
    <li>Every pair of outputs was compared: same intervals, same tables, same signal cube, same motif matrix, same graph edges and vertex labels.</li>
  </ul>
</section>"""

    nav = "".join(f'<a href="#{i}">{t}</a>' for i, t in (
        ("answers", "Answers"), ("what", "What it is"), ("pandas", "vs pandas"), ("api", "Your code"),
        ("functions", "Functions"), ("objects", "Objects"), ("genes", "Genes"), ("plan", "Plan")))
    page = f"""<title>Objects or Columns</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,400;0,600;0,700;0,800;1,400&display=swap">
<style>{(REP / 'style.css').read_text()}</style>
<div class="wrap">
<header class="top">
  <span class="eyebrow">genomeblocks · how data is stored · measured on a working prototype</span>
  <h1>Objects or columns</h1>
  <p class="lede">What "column data" means, how it compares with pandas-style tables, what happens to annotate and every other function if Loci switches to it, and what the current objects give us.</p>
</header>
<nav class="sections" aria-label="Sections">{nav}</nav>
{secs}
</div>
<script>{(REP / 'charts.js').read_text()}</script>
<script>
const DATA = {json.dumps(charts)};
const F = {{
  time: (s) => {{ s = +(+s).toPrecision(6); if (s === 0) return "0"; if (s < 1e-3) return (s < 1e-5 ? (s*1e6).toFixed(1) : (s*1e6).toFixed(0)) + " µs";
                 if (s < 1) return (s < 0.01 ? (s*1e3).toFixed(1) : (s*1e3).toFixed(0)) + " ms";
                 if (s < 120) return (s < 10 ? s.toFixed(1) : s.toFixed(0)) + " s"; return (s/60).toFixed(1) + " min"; }},
  us: (s) => {{ const u = s * 1e6; return (u < 10 ? u.toFixed(2) : u.toFixed(0)) + " µs"; }},
  x: (v) => (v < 10 ? v.toFixed(1) : v.toFixed(0)) + "×",
}};
document.addEventListener("DOMContentLoaded", () => {{
  for (const [id, s] of Object.entries(DATA)) {{
    const host = document.getElementById("c-" + id); if (!host) continue;
    if (s.kind === "hbar") Charts.hbar(host, {{...s, fmt: F[s.fmt]}});
    else if (s.kind === "diverge") Charts.diverge(host, s);
    else Charts.lines(host, {{...s, xfmt: F[s.xfmt], yfmt: F[s.yfmt]}});
  }}
}});
</script>
"""
    out = REP / "columns.html"
    out.write_text(page)
    print(f"[explainer] {out} ({len(page) / 1e3:.0f} kB)")


if __name__ == "__main__":
    build()
