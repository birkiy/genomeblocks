#!/usr/bin/env python3
"""Build report/prototype.html: the columnar prototype measured against main.

Numbers come from results/prototype.json (bench_prototype.py).
"""
from __future__ import annotations

import json
import math
import re

from build_report import REP, chart, details, esc, ftime, fx, load, table


def diagrams():
    raw = (REP / "diagrams_prototype.html").read_text()
    return {m.group(1): m.group(2).strip()
            for m in re.finditer(r"<!-- @(\w+) -->\s*(.*?)(?=<!-- @|\Z)", raw, re.S)}


def fig(svg, title, caption):
    return (f'<div class="why"><h3>{esc(title)}</h3><div class="scroll">{svg}</div>'
            f'<p class="sub">{caption}</p></div>')


def blocks_svg(meta):
    """Edge table as a strip of blocks + the chrom x chrom matrix (real counts)."""
    chroms, M = meta["block_counts"]["chroms"], meta["block_counts"]["matrix"]
    n = len(chroms)
    cis = [M[i][i] for i in range(n)]
    trans = meta["n_trans"]
    total = sum(cis) + trans
    W, x0, x1 = 760, 20, 740
    out = [f'<svg viewBox="0 0 {W} 466" role="img" aria-label="The edge table sorted into one block per '
           f'chromosome ({sum(cis):,} cis edges) followed by one trans block ({trans:,} edges), and the '
           f'same edges as a chromosome-by-chromosome matrix: large counts on the diagonal, small counts '
           f'off it.">']
    out.append('<text x="20" y="22" class="head">The edge table, top to bottom (row numbers)</text>')
    x = x0
    for i, c in enumerate(cis + [trans]):
        w = (x1 - x0) * c / total
        cls = "cell on" if i < n else "cell q"
        op = 0.35 + 0.65 * (i % 2) if i < n else 1
        out.append(f'<rect x="{x:.1f}" y="34" width="{max(w, 1.5):.1f}" height="26" class="{cls} strip" '
                   f'style="fill-opacity:{op:.2f}"><title>{"trans" if i == n else chroms[i]}: {c:,} edges</title></rect>')
        if i < 3 or i == n:
            lab = "trans" if i == n else f"cis {chroms[i]}"
            anchor = "end" if i == n else "start"
            lx = x + w if i == n else x + 3
            out.append(f'<text x="{lx:.1f}" y="76" class="small" text-anchor="{anchor}">{lab}</text>')
        x += w
    out.append(f'<text x="{x0 + (x1 - x0) * 0.42:.0f}" y="76" class="small">… one block per chromosome …</text>')
    out.append(f'<text x="20" y="96" class="small">{sum(cis):,} cis edges in {n} blocks, then '
               f'{trans:,} trans edges in one block. Each block is a slice: A.chrom(c), A.cis, A.trans.</text>')
    # matrix: one hue, darker = more edges (log)
    top, size, mx = 150, 300, 240
    cell = size / n
    out.append('<text x="20" y="122" class="head">The same edges as a matrix</text>')
    out.append('<text x="20" y="140" class="small">chromosome × chromosome edge counts, log colour scale</text>')
    vmax = math.log10(max(max(r) for r in M) + 1)
    for i in range(n):
        for j in range(n):
            v = M[i][j]
            op = 0.04 if v == 0 else 0.08 + 0.92 * (math.log10(v + 1) / vmax) ** 2
            out.append(f'<rect x="{mx + j * cell:.1f}" y="{top + i * cell:.1f}" width="{cell:.1f}" '
                       f'height="{cell:.1f}" class="cell on heat" style="fill-opacity:{op:.2f}">'
                       f'<title>{chroms[i]} × {chroms[j]}: {v:,} edges</title></rect>')
    for i, c in enumerate(chroms):
        if i % 2 == 0 or i == n - 1:
            out.append(f'<text x="{mx - 6}" y="{top + (i + 0.7) * cell:.1f}" class="small" '
                       f'text-anchor="end">{c}</text>')
    offd = [M[i][j] for i in range(n) for j in range(n) if i != j]
    lx = mx + size + 24
    out.append(f'<rect x="{lx}" y="{top + 10}" width="14" height="14" class="cell on"/>'
               f'<text x="{lx + 22}" y="{top + 22}" class="small">diagonal: one cis block</text>'
               f'<text x="{lx + 22}" y="{top + 38}" class="small">{min(cis):,}–{max(cis):,} edges each</text>')
    out.append(f'<rect x="{lx}" y="{top + 56}" width="14" height="14" class="cell on" style="fill-opacity:.2"/>'
               f'<text x="{lx + 22}" y="{top + 68}" class="small">off-diagonal: a trans pair</text>'
               f'<text x="{lx + 22}" y="{top + 84}" class="small">{min(offd):,}–{max(offd):,} edges each</text>')
    for k, line in enumerate(("One graph. The diagonal is most of it;", "the trans pairs stay in the same table,",
                              "so neighbours and graph algorithms", "always see them.")):
        out.append(f'<text x="{lx}" y="{top + 120 + 16 * k}" class="small">{line}</text>')
    out.append("</svg>")
    return "".join(out)


def build():
    D = diagrams()
    R = load("prototype")
    rows = R["rows"]
    st = {}
    for r in rows:
        if r["part"] == "steps":
            st.setdefault(r["step"], {})[r["impl"]] = r["seconds"]
    meta = next(r for r in rows if r["part"] == "steps_meta")
    vw = {}
    for r in rows:
        if r["part"] == "views":
            vw.setdefault(r["op"], {})[r["impl"]] = r["seconds"]
    tr = {}
    for r in rows:
        if r["part"] == "trans":
            tr.setdefault(r["step"], {})[r["impl"]] = r
    tmeta = next(r for r in rows if r["part"] == "trans_meta")
    sc = {}
    for r in rows:
        if r["part"] == "scale":
            sc.setdefault(r["impl"], []).append((r["loops"], r["seconds"], r["n_edges"]))
    e2e = {}
    for r in rows:
        if r["part"] == "e2e":
            e2e.setdefault(r["impl"], {})[r["step"]] = r["seconds"]
    mem = {r["impl"]: r for r in rows if r["part"] == "e2e_mem"}
    reload_ = {r["impl"]: r for r in rows if r["part"] == "e2e_load"}
    main_commit = mem["main"]["main_commit"]

    tot_m, tot_c = e2e["main"]["total"], e2e["columnar"]["total"]
    arch_steps = ["Architecture.make (50k loops)", "add_mcool (5 kb)", "normalize (power-law O/E)",
                  "annotate (labels + genes)", "strength", "prime_hubs"]
    arch_m = sum(st[s]["main"] for s in arch_steps)
    arch_c = sum(st[s]["columnar"] for s in arch_steps)
    fail = next((s for s, d in tr.items() if d.get("main", {}).get("status", "ok") not in ("ok", "not reached")), None)
    fail_msg = tr[fail]["main"]["status"] if fail else ""

    # ── charts ───────────────────────────────────────────────────────────
    charts = {}
    order = ["Loci.make (100k CREs)", "Genes.make (GTF, 877k lines)", "gene annotation index",
             "Genes.annotations (100k CREs)"] + arch_steps + ["save", "load"]
    LABEL = {"Architecture.make (50k loops)": "Architecture.make (50k loops)",
             "annotate (labels + genes)": "Architecture.annotate", "save": "save the graph",
             "load": "load it back"}
    charts["steps"] = {"kind": "dumbbell", "names": ["main", "columnar prototype"], "fmt": "time",
                       "labelW": 260, "axis": "seconds per call (log scale) · lower is better",
                       "rows": [{"label": LABEL.get(s, s), "a": st[s]["main"], "b": st[s]["columnar"]}
                                for s in order]}
    groups = [("k1", "add_mcool", ["add_mcool"]), ("k2", "Genes.make", ["Genes.make"]),
              ("k3", "make + annotate", ["make", "annotate"]), ("k4", "save", ["save"]),
              ("k5", "everything else", ["import", "Loci.make", "normalize", "strength", "prime_hubs"])]
    charts["e2e"] = {"kind": "stack", "fmt": "time", "labelW": 150,
                     "axis": "wall time, fresh Python process (seconds)",
                     "keys": [{"key": k, "name": n, "cls": k} for k, n, _ in groups],
                     "rows": [{"label": lab, "parts": [{"key": k, "value": sum(e2e[impl].get(s, 0) for s in ss)}
                                                       for k, _, ss in groups]}
                              for impl, lab in (("main", "main"), ("columnar", "columnar"))]}
    vorder = ["one chromosome (chr8)", "region chr8:20-40 Mb", "neighbours of 1,000 CREs", "copy",
              "connected components (graph-tool)", "pagerank (graph-tool)"]
    charts["views"] = {"kind": "dumbbell", "names": ["main", "columnar prototype"], "fmt": "time",
                       "labelW": 260, "axis": "seconds (log scale) · lower is better",
                       "rows": [{"label": o, "a": vw[o]["main"], "b": vw[o]["columnar"]} for o in vorder]}
    charts["scale"] = {"kind": "lines", "xlog": True, "ylog": True, "xfmt": "count", "yfmt": "time",
                       "xLabel": "loops in the BEDPE (100k CREs)", "yLabel": "pipeline time",
                       "endLabels": True,
                       "series": [{"name": "main", "color": "var(--s2)",
                                   "pts": [[n, s] for n, s, _ in sorted(sc["main"])]},
                                  {"name": "columnar", "color": "var(--s1)",
                                   "pts": [[n, s] for n, s, _ in sorted(sc["columnar"])]}]}
    cmap = {n: t for n, t, _ in sc["columnar"]}
    gap_txt = ", ".join(f"{fx(t / cmap[n])} at {n / 1000:g}k loops" for n, t, _ in sorted(sc["main"]))
    big_c = max(sc["columnar"])
    big_m = max(sc["main"])

    # ── tables ───────────────────────────────────────────────────────────
    step_rows = [[LABEL.get(s, s), ftime(st[s]["main"]), ftime(st[s]["columnar"]),
                  fx(st[s]["main"] / st[s]["columnar"])] for s in order]
    view_rows = [[o, ftime(vw[o]["main"]) if "main" in vw[o] else "—", ftime(vw[o]["columnar"]),
                  fx(vw[o]["main"] / vw[o]["columnar"]) if "main" in vw[o] else "—"]
                 for o in vw]
    trans_rows = []
    for s in ("make", "add_mcool", "normalize", "annotate", "strength"):
        m = tr[s].get("main", {})
        ms = ftime(m["seconds"]) if m.get("seconds") is not None else (
            f"fails ({m['status'].split(':')[0]})" if m.get("status") not in (None, "ok", "not reached")
            else "not reached")
        trans_rows.append([s, ms, ftime(tr[s]["columnar"]["seconds"])])
    e2e_rows = [[k, ftime(e2e["main"].get(k, 0)), ftime(e2e["columnar"].get(k, 0))]
                for k in e2e["columnar"] if k != "total"]
    e2e_rows.append(["total", ftime(tot_m), ftime(tot_c)])

    secs = f"""
<section id="answers">
  <h2>What the prototype shows</h2>
  <div class="tiles two-up">
    <a class="tile" href="#speed"><span class="tag green">speed</span><span class="k">Whole pipeline, fresh process</span>
      <span class="v">{ftime(tot_m)} → {ftime(tot_c)}</span><span class="d">Load CREs and genes, build the graph, add Hi-C, normalise, annotate, find hubs, save. {fx(tot_m / tot_c)} faster end to end; the Architecture steps alone go from {ftime(arch_m)} to {ftime(arch_c)}.</span></a>
    <a class="tile" href="#same"><span class="tag navy">same answers</span><span class="k">Checked against main</span>
      <span class="v">identical</span><span class="d">Same edges, Hi-C weights, O/E, distances, region labels, strengths and hubs. Gene picks match too, except 6 of 73k CREs where two promoters tie on weight and either answer is valid.</span></a>
    <a class="tile" href="#trans"><span class="tag warn">trans</span><span class="k">Inter-chromosomal loops</span>
      <span class="v">{tmeta['n_trans']:,} kept</span><span class="d">On main, any trans edge makes <code>normalize</code> stop with a TypeError. The prototype normalises them against the mean trans weight and keeps them in neighbours and graph algorithms. (Fixed in 1.1.0: the classic <code>normalize</code> now does the same.)</span></a>
    <a class="tile" href="#notebook"><span class="tag purple">notebook</span><span class="k">Cut, look up, reload</span>
      <span class="v">µs–ms</span><span class="d">One chromosome: {ftime(vw['one chromosome (chr8)']['columnar'])} instead of {ftime(vw['one chromosome (chr8)']['main'])}. Reloading a saved graph: {ftime(reload_['columnar']['seconds'])} instead of {ftime(reload_['main']['seconds'])}, from a file {meta['file_mb_main'] / meta['file_mb_columnar']:.0f}× smaller.</span></a>
  </div>
</section>

<section id="design">
  <span class="tag navy">Design</span>
  <h2>Tables that line up by row</h2>
  <p class="verdict">Every table shares one <b>Genome</b>, and the <b>row number</b> links the tables. Row i of the CREs is row i of the labels, of the signal cube and of every graph column. Graph edges store row numbers, so nothing needs a lookup dictionary.</p>
  {fig(D['structure'], "The notebook structure",
       "Loci, Genes and Architecture are plain numpy columns. Looking at one row still gives a normal <code>Locus</code> (a view onto the columns), so code that reads <code>cre[i].start</code> or loops over peaks keeps working. Classic functions such as <code>signal()</code> run unchanged on the new Loci and return the same cube.")}
</section>

<section id="blocks">
  <span class="tag navy">Architecture</span>
  <h2>One graph, sorted into blocks</h2>
  <p class="verdict">Most loops are within one chromosome, and a few cross chromosomes. Keep them all in <b>one</b> edge table, sorted so each chromosome's edges form one block and all trans edges form a final block. Nothing is split, but every part can be sliced out.</p>
  <div class="why"><h3>The edge table of the benchmark graph ({tmeta['n_edges']:,} edges)</h3><div class="scroll">{blocks_svg(tmeta)}</div>
    <p class="sub">Real counts from the trans benchmark: {tmeta['n_cis']:,} cis edges and {tmeta['n_trans']:,} trans edges (5% of loops). A per-chromosome analysis takes one block. A whole-genome analysis takes the whole table. Neither needs a copy.</p></div>
  {fig(D['views'], "Taking one chromosome",
       f"On main this is a Python loop that rebuilds a graph ({ftime(vw['one chromosome (chr8)']['main'])} for chr8). In the prototype it is a slice of the same arrays ({ftime(vw['one chromosome (chr8)']['columnar'])}), and edits to the slice write through to the full graph.")}
  <div class="mech">
    <div><h3>Per chromosome</h3><p><code>A.chrom("chr8")</code>, or <code>for c, view in A.chroms()</code>. Each block is independent, so blocks can go to separate processes.</p></div>
    <div><h3>Cis or trans only</h3><p><code>A.cis</code>, <code>A.trans</code>: slices again. Distance-decay fits use cis; trans gets its own flat expectation.</p></div>
    <div><h3>Whole graph</h3><p><code>A.neighbors(row)</code> uses one adjacency index over all edges, so a CRE's trans partners are always listed. <code>A.graph()</code> builds a graph-tool Graph from the arrays when you need real graph algorithms.</p></div>
  </div>
</section>

<section id="speed">
  <span class="tag green">Speed</span>
  <h2>Every step, main vs prototype</h2>
  <p class="verdict">Same inputs, same outputs, every step faster. The biggest gain is <b>add_mcool</b>: {ftime(st['add_mcool (5 kb)']['main'])} → {ftime(st['add_mcool (5 kb)']['columnar'])}. The prototype maps every edge to its Hi-C pixel with array lookups, reading pixels one chromosome block at a time, instead of one pandas lookup per edge.</p>
  {chart("steps", "Per step", "· 100k CREs, 50k loops, 5 kb Hi-C, medians", "Each row is one call on the same data. The orange dot is main and the blue dot is the prototype; the label on the right is the speed-up. Steps on main that are already fast (strength, prime_hubs) are still loops over graph-tool property maps, and in the prototype they become one or two numpy calls.")}
  {chart("e2e", "Where the time goes", "· one fresh Python process per implementation", f"main imports genomeblocks from a clean checkout of origin/main ({esc(main_commit)}). Peak memory: {mem['main']['peak_mb']:,.0f} MB on main, {mem['columnar']['peak_mb']:,.0f} MB for the prototype.")}
  {details("Per-step table", table(["step", "main", "columnar", "speed-up"], step_rows, num_cols=(1, 2, 3)))}
  {details("End-to-end table (fresh process)", table(["step", "main", "columnar"], e2e_rows, num_cols=(1, 2)))}
</section>

<section id="notebook">
  <span class="tag purple">Notebook</span>
  <h2>The interactive part</h2>
  <p class="verdict">In a notebook you cut, look up and reload all the time. Views and lookups take <b>microseconds to milliseconds</b>. Graph algorithms still run in graph-tool, at the same speed, built from the tables on demand.</p>
  {chart("views", "Everyday operations", "· same graph, medians", f"Building the graph-tool Graph from the tables takes {ftime(vw['graph-tool Graph from tables']['columnar'])} and is cached. Components and PageRank are graph-tool's own code on both sides, so they match.")}
  <div class="tiles">
    <div class="tile"><span class="k">Reload a saved graph (fresh process)</span><span class="v">{ftime(reload_['columnar']['seconds'])}</span><span class="d">main: {ftime(reload_['main']['seconds'])} to unpickle and rebuild the graph.</span></div>
    <div class="tile"><span class="k">File on disk</span><span class="v">{meta['file_mb_columnar']:.1f} MB</span><span class="d">Parquet tables that pandas, polars and R can read; main's pickle: {meta['file_mb_main']:.1f} MB, Python only.</span></div>
    <div class="tile"><span class="k">Peak memory, whole pipeline</span><span class="v">{mem['columnar']['peak_mb']:,.0f} MB</span><span class="d">main: {mem['main']['peak_mb']:,.0f} MB in the same pipeline.</span></div>
  </div>
  {details("All timings", table(["operation", "main", "columnar", "speed-up"], view_rows, num_cols=(1, 2, 3)))}
</section>

<section id="scale">
  <span class="tag green">Scale</span>
  <h2>Bigger loop sets</h2>
  <p class="verdict">The prototype runs <b>{big_c[0]:,} loops</b> ({big_c[2]:,} edges) in {ftime(big_c[1])}. main needs {ftime(big_m[1])} for {big_m[0]:,} loops; it was not run at 1M.</p>
  {chart("scale", "Pipeline time vs loop count", "· log–log, make → add_mcool → normalize → annotate → strength", "The gap widens with size: " + gap_txt + ". main does Python work per edge, so its time grows with the edge count. The prototype's fixed costs (reading the cooler and GTF index, fitting the power law) dominate at small sizes.")}
</section>

<section id="trans">
  <span class="tag warn">Trans loops</span>
  <h2>Inter-chromosomal edges</h2>
  <p class="verdict">With 5% trans loops, main builds the graph and adds Hi-C, then <b>stops at normalize</b> ({esc(fail_msg.split(':')[0]) if fail_msg else 'error'}). <code>Locus.distance_to</code> has no distance between chromosomes, and normalize writes that missing value into a float column. The prototype gives trans edges an infinite distance and normalises them against the mean trans contact. (Fixed in 1.1.0: the classic <code>normalize</code> now does the same.)</p>
  <div class="two">
    <figure><div class="ttl">Same pipeline, 52.5k loops (2.5k trans)</div>{table(["step", "main", "columnar"], trans_rows, num_cols=(1, 2))}</figure>
    <figure><div class="ttl">What trans edges add</div>{table(["", "value"], [
        ["trans edges", f"{tmeta['n_trans']:,}"],
        ["with Hi-C contacts", f"{tmeta['trans_with_contacts']:,}"],
        ["median O/E, trans (with contacts)", f"{tmeta['trans_n_median']:.2f}"],
        ["median O/E, cis (with contacts)", f"{tmeta['cis_n_median']:.2f}"],
        ["connected components, cis only", f"{tmeta['components_cis_only']:,}"],
        ["connected components, with trans", f"{tmeta['components_all']:,}"],
        ["largest component, cis only → with trans", f"{tmeta['largest_component_cis_only']:,} → {tmeta['largest_component_all']:,} CREs"],
    ], num_cols=(1,))}</figure>
  </div>
  <p class="sub">The synthetic Hi-C adds about 40 read pairs at each trans loop, so those edges carry real weight. The flat trans expectation is a first choice. A per-chromosome-pair expectation would be the obvious next refinement.</p>
</section>

<section id="try">
  <span class="tag green">Try it</span>
  <h2>Using it in a notebook</h2>
  <pre>import genomeblocks.columnar as gbc

cre   = gbc.Loci.make("atac.narrowPeak")          # sorted, one block per chromosome
genes = gbc.Genes.make("gencode.gtf")              # polars parse when installed
A = (gbc.Architecture.make(cre, "loops.bedpe")
       .add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes).strength())

A.chrom("chr8"); A.cis; A.trans                    # views, no copy
A.neighbors("chr8:127735000-127736000(.)")         # cis + trans partners
A.vp.strength[cre.row(uid)]                        # every column is aligned with cre
cube = cre.take(rows).signal(bigwigs)              # classic functions still work
g = A.graph()                                      # graph-tool for algorithms
A.save("session/arch"); gbc.Architecture.load("session/arch")</pre>
  <p class="sub">The notebook <code>examples/columnar_prototype/columnar_prototype.ipynb</code> on branch <code>columnar-prototype</code> walks through all of this with outputs: block matrix, views, trans partners, graph-tool components and PageRank, a hub signal profile, a 1,000-permutation genome-wide test, save/reload and hand-off to pandas, polars and Arrow.</p>
</section>

<section id="limits">
  <span class="tag warn">Limits</span>
  <h2>What the prototype does not do yet</h2>
  <div class="fixes">
    <div class="fix"><span class="t">Classic-only features</span><p><code>Genes.select_isoforms</code>, <code>make_ucsc</code> and the Architecture draw helpers are not ported. Use <code>A.to_legacy()</code> to get a classic graph for drawing.</p></div>
    <div class="fix"><span class="t">Coolers</span><p>Fixed-size bins only. As on main, each CRE maps to the bin holding its start.</p></div>
    <div class="fix"><span class="t">Interval kernel</span><p>The all-pairs overlap used by <code>make</code> scans a window as wide as the longest reference interval. That suits CREs; it would need cgranges or a tree for megabase-long intervals.</p></div>
    <div class="fix"><span class="t">Graph-tool vertices</span><p><code>A.graph()</code> has one vertex per CRE row, isolated CREs included, so vertex i is row i. Algorithms that normalise over all vertices (PageRank) are scaled differently from main's graph, which holds linked CREs only.</p></div>
    <div class="fix"><span class="t">Edge cases kept simple</span><p>Duplicate CRE rows stay separate vertices (main merges them by uid). Ensembl-style generic "UTR" lines are classed 5' or 3' by the transcript's first CDS base. Without polars, GTF parsing falls back to pandas, which is slow (about 11 s here).</p></div>
  </div>
</section>

<section id="method">
  <span class="tag navy">Method</span>
  <h2>How this was measured</h2>
  <ul class="plain">
    <li>Code: <code>genomeblocks/columnar/</code> on branch <code>columnar-prototype</code>. The classic modules on that branch are byte-identical to origin/main ({'checked' if meta['identical_to_main'] else 'NOT identical'}), so the in-process "main" numbers run main's code. The fresh-process runs import main from a separate checkout ({esc(main_commit)}).</li>
    <li>Tests: <code>tests/test_columnar.py</code> compares every step with the classic implementation, and the full suite passes. The benchmark repeats the comparison at full size.</li>
    <li>Data: the synthetic hg38-shaped set from the <a href="https://claude.ai/artifact/4hL67NF3y68wJSdsySs1Eo">benchmark report</a>: 100k CREs, 20k genes, 50k loops (52.5k with trans), a 5 kb cooler from 5M read pairs. 4 vCPU Xeon at 2.1 GHz; medians of 3–5 runs after a warm-up (main's add_mcool: one run).</li>
    <li>Bench: <code>benchmarks/bench_prototype.py</code>; page: <code>benchmarks/build_prototype_report.py</code>.</li>
  </ul>
</section>"""

    nav = "".join(f'<a href="#{i}">{t}</a>' for i, t in (
        ("answers", "Summary"), ("design", "Design"), ("blocks", "Blocks"), ("speed", "Speed"),
        ("notebook", "Notebook"), ("scale", "Scale"), ("trans", "Trans"), ("try", "Try it"),
        ("limits", "Limits")))
    page = f"""<title>Columnar Architecture Prototype</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,400;0,600;0,700;0,800;1,400&display=swap">
<style>{(REP / 'style.css').read_text()}</style>
<div class="wrap">
<header class="top">
  <span class="eyebrow">genomeblocks · branch columnar-prototype · measured against main</span>
  <h1>One graph, two tables</h1>
  <p class="lede">A working prototype of the notebook design. CREs, genes and the Architecture are stored as tables that share one genome and line up by row. Intra- and inter-chromosomal edges live in one sorted edge table. Every result is checked against main.</p>
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
  count: (v) => v >= 1e6 ? (v/1e6) + "M" : v >= 1e3 ? (+(v/1e3).toPrecision(3)) + "k" : String(v),
}};
document.addEventListener("DOMContentLoaded", () => {{
  for (const [id, s] of Object.entries(DATA)) {{
    const host = document.getElementById("c-" + id); if (!host) continue;
    if (s.kind === "dumbbell") Charts.dumbbell(host, {{...s, fmt: F[s.fmt]}});
    else if (s.kind === "stack") Charts.stack(host, {{...s, fmt: F[s.fmt]}});
    else Charts.lines(host, {{...s, xfmt: F[s.xfmt], yfmt: F[s.yfmt]}});
  }}
}});
</script>
"""
    out = REP / "prototype.html"
    out.write_text(page)
    print(f"[prototype] {out} ({len(page) / 1e3:.0f} kB)")


if __name__ == "__main__":
    build()
