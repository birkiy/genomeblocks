#!/usr/bin/env python3
"""Build the pluggable-backends design proposal (report/backends.html).

Numbers come from results/backends.json, genes_make.json and
motif_backends.json; see bench_backends.py, bench_genes_make.py and
bench_motif_backends.py for how each was measured.
"""
from __future__ import annotations

import json
import re

from build_report import REP, chart, details, esc, fnum, ftime, fx, load, table, why


def diagrams():
    raw = (REP / "diagrams_backends.html").read_text()
    return {m.group(1): m.group(2).strip()
            for m in re.finditer(r"<!-- @(\w+) -->\s*(.*?)(?=<!-- @|\Z)", raw, re.S)}


def build():
    D = diagrams()
    charts = {}
    bk = load("backends")["rows"]
    gm = load("genes_make")["rows"]
    mb = load("motif_backends")["rows"]

    # ── intervals ───────────────────────────────────────────────────────────
    ops = [r for r in bk if r["part"] == "ops"]
    COL = {"cgranges (today)": "var(--s8)", "numpy, converting from objects": "var(--s4)",
           "numpy, columnar (no conversion)": "var(--s1)", "pyranges, native frames": "var(--s2)",
           "pyranges, converting from Loci": "var(--s5)"}
    LAB = {"cgranges (today)": "today", "numpy, converting from objects": "numpy, converting",
           "numpy, columnar (no conversion)": "numpy on columns", "pyranges, native frames": "pyranges",
           "pyranges, converting from Loci": "pyranges via Loci"}

    def series(op):
        out = []
        for e, c in COL.items():
            pts = sorted((r["n"], r["seconds"]) for r in ops if r["op"] == op and r["engine"] == e)
            if pts:
                out.append({"name": LAB[e], "color": c, "pts": pts})
        return out
    charts["iv_and"] = {"kind": "lines", "xlog": True, "ylog": True, "xfmt": "num", "yfmt": "time",
                        "xLabel": "peaks per set (n)", "yLabel": "time for A & B", "endLabels": True,
                        "series": series("A & B")}
    charts["iv_merge"] = {"kind": "lines", "xlog": True, "ylog": True, "xfmt": "num", "yfmt": "time",
                          "xLabel": "intervals merged", "yLabel": "time for sort + merge", "endLabels": True,
                          "series": series("merge")}
    get = lambda op, e, n: next(r["seconds"] for r in ops if r["op"] == op and r["engine"] == e and r["n"] == n)
    and_today, and_col = get("A & B", "cgranges (today)", 1_000_000), get("A & B", "numpy, columnar (no conversion)", 1_000_000)
    and_conv, and_pr = get("A & B", "numpy, converting from objects", 1_000_000), get("A & B", "pyranges, native frames", 1_000_000)
    m_today, m_col = get("merge", "cgranges (today)", 2_000_000), get("merge", "numpy, columnar (no conversion)", 2_000_000)
    m_pr = get("merge", "pyranges, native frames", 2_000_000)
    all_same = all(r.get("same", True) for r in ops)
    conv = {(r["step"], r["n"]): r["seconds"] for r in bk if r["part"] == "convert"}
    disp = {r["step"]: r["seconds"] for r in bk if r["part"] == "dispatch"}
    look = {r["engine"]: r["seconds"] for r in bk if r["part"] == "lookup"}
    conv_rows = [[s, f"{n:,}", ftime(t)] for (s, n), t in conv.items()]

    # ── genes ───────────────────────────────────────────────────────────────
    parse = {r["variant"]: r for r in gm if r["part"] == "parse"}
    ann = {r["variant"]: r for r in gm if r["part"] == "annot"}
    comp = {r["step"]: r["seconds"] for r in gm if r["part"] == "components"}
    cur = parse["Genes.make (current)"]["seconds"]
    pol = parse["polars read_csv + str.extract (table)"]["seconds"]
    pol_obj = parse["polars table, then the same objects"]["seconds"]
    a_cur = ann["current (Loci.sort/merge)"]["seconds"]
    a_key = ann["key-sort fix (same Loci objects)"]["seconds"]
    a_tab = ann["numpy merge (from columnar table)"]["seconds"]
    charts["genes_parse"] = {"kind": "hbar", "fmt": "time", "log": False, "sort": "asc", "labelW": 330,
                             "axis": "seconds to parse an 877k-line GTF",
                             "rows": [{"label": k.replace(" (table)", " → table"), "value": r["seconds"],
                                       "hl": k.startswith(("polars", "Genes.make (current)"))} for k, r in parse.items()]}
    charts["genes_annot"] = {"kind": "hbar", "fmt": "time", "log": True, "sort": "asc", "labelW": 330,
                             "axis": "seconds to build the promoter / exon / UTR index (log scale)",
                             "rows": [{"label": k, "value": r["seconds"], "hl": "columnar" in k or "current" in k}
                                      for k, r in ann.items()]}

    # ── motifs ──────────────────────────────────────────────────────────────
    sc = [r for r in mb if r["part"] == "scan" and r["n_seqs"] == 5000]
    prune = next(r["columns_scored_fraction"] for r in mb if r["part"] == "pruning")
    mm = {(r["engine"], r["workers"]): r for r in sc}
    lm1, lm4 = mm[("lightmotif (scan_motifs_matrix, today)", 1)]["seconds"], mm[("lightmotif (scan_motifs_matrix, today)", 4)]["seconds"]
    mo1, mo4 = mm[("MOODS backend", 1)]["seconds"], mm[("MOODS backend", 4)]["seconds"]
    la1, la4 = mm[("lookahead scanner (numba prototype)", 1)]["seconds"], mm[("lookahead scanner (numba prototype)", 4)]["seconds"]
    diff_cells = max(r.get("cells_differing", 0) for r in sc)
    cells = sc[0].get("cells", 0) or max(r.get("cells", 0) for r in sc)
    SHORT = {"lightmotif (scan_motifs_matrix, today)": "lightmotif (today)", "MOODS backend": "MOODS backend",
             "lookahead scanner (numba prototype)": "lookahead prototype (numba)"}
    charts["motif_backends"] = {"kind": "hbar", "fmt": "time", "log": False, "sort": None, "labelW": 300,
                                "axis": "seconds for 5,000 windows × 1,019 motifs",
                                "rows": [{"label": f"{SHORT[r['engine']]} · {r['workers']} {'core' if r['workers'] == 1 else 'workers'}",
                                          "value": r["seconds"], "hl": r["engine"] == "MOODS backend"} for r in sc]}

    lo = load("loci")["rows"]
    fb = next(r["seconds"] for r in lo if r["op"] == "intersect" and r["n"] == 100_000 and r["engine"] == "genomeblocks (pure-Python fallback)")
    bd = next(r["seconds"] for r in lo if r["op"] == "intersect" and r["n"] == 100_000 and r["engine"] == "genomeblocks (fallback, bounded scan)")

    # ── interop table ───────────────────────────────────────────────────────
    interop = [
        ["pandas DataFrame", "to_frame / from_frame", f"to_frame builds a dict per row ({ftime(conv[('Loci.to_frame() (pandas)', 1_000_000)])} per 1M); from_frame drops strand", f"{ftime(conv[('numpy arrays -> pandas frame', 1_000_000)])} per 1M from columns"],
        ["pyranges", "to_pyranges", "goes through to_frame; no from_pyranges", "zero-copy both ways"],
        ["polars / Arrow", "none", "", "zero-copy columns"],
        ["bioframe", "none (works via to_frame)", "column names differ (Chr/Start/End)", "chrom/start/end frames"],
        ["pybedtools BedTool", "none (via to_bed file)", "", "to_bedtool / from_bedtool"],
        ["BED / narrowPeak", "Loci.make, to_bed", "no .gz; name, score, signal, summit dropped; to_bed() without a path returns None", "keep extra columns; gzip; return text"],
        ["GTF / GFF3", "Genes.make, make_ucsc", "GFF3 fails (KeyError 'gene_id'); no .gz; GTF starts kept 1-based (1 bp off BED and pyranges)", "GFF3 attributes; 0-based starts"],
        ["AnnData / MuData", "none", "", "peaks as var names ('chr:start-end') ↔ Loci"],
        ["xarray", "none (signal cube is a bare array)", "", "cube with loci / track / bin coordinates"],
        ["graphs", "Architecture is a graph-tool Graph", "GraphML/GML already work through graph-tool's save()", "edge-list export for networkx / igraph"],
        ["cooler, bigWig, BAM, MEME", "add_mcool, signal, browser, write_meme", "", "keep"],
    ]

    secs = f"""
<section id="answers">
  <h2>Short answers</h2>
  <div class="tiles two-up">
    <a class="tile" href="#intervals"><span class="tag green">intervals</span><span class="k">cgranges + pyranges as backends?</span>
      <span class="v">{fx(and_today / and_col)}</span><span class="d">faster A &amp; B at 1M with a numpy engine written from scratch, beating pyranges ({ftime(and_pr)}). The decision costs {ftime(disp['auto policy (size + op check)'])}; converting objects costs more than the engine saves ({ftime(and_conv)}). Columns first, then engines.</span></a>
    <a class="tile" href="#genes"><span class="tag navy">Genes.make</span><span class="k">Does fixing sort/merge help?</span>
      <span class="v">{fx(cur / pol)}</span><span class="d">faster parsing with polars into a table. The sort/merge fix speeds the first annotation call {fx(a_cur / a_key)} but leaves make untouched; a table-based index is {fx(a_cur / a_tab)} faster.</span></a>
    <a class="tile" href="#motifs"><span class="tag green">motifs</span><span class="k">MOODS as the backend?</span>
      <span class="v">{fx(lm1 / mo1)}</span><span class="d">faster than today on one core, identical matrix ({diff_cells} of {cells:,} cells differ). A MOODS-style scanner written in numba is correct but {fx(la1 / lm1)} slower than today.</span></a>
    <a class="tile" href="#interop"><span class="tag warn">interop</span><span class="k">Compatible with other tools?</span>
      <span class="v">5 gaps</span><span class="d">GFF3 fails, .gz is not read, GTF starts are 1 bp off other tools, extra BED columns are dropped, and no zero-copy path to pandas, polars, Arrow, AnnData.</span></a>
  </div>
</section>

<section id="idea">
  <span class="tag navy">The idea</span>
  <h2>Store columns, plug in engines</h2>
  <p class="verdict">Any backend can be swapped in cheaply once genomeblocks keeps its data as columns instead of Python objects. The switch is a {ftime(disp['auto policy (size + op check)'])} decision; the conversion around it is what costs.</p>
  {why("Where the time goes when an engine is swapped in", D['tax'],
       f"Measured on 1M peaks: <code>Loci</code> → arrays takes {ftime(conv[('Loci -> numpy arrays', 1_000_000)])}, arrays → <code>Locus</code> objects {ftime(conv[('numpy arrays -> Loci objects', 1_000_000)])}, while the numpy engine needs {ftime(and_col)} for the whole A &amp; B. Keeping the columns as the storage removes both conversions and makes pandas, polars, Arrow and pyranges views nearly free ({ftime(conv[('numpy arrays -> pandas frame', 1_000_000)])} for a 1M-row pandas frame).")}
  <div class="why"><h3>Proposed shape</h3><div class="scroll">{D['arch']}</div>
    <p class="sub">Each kind of work has a small protocol (for intervals: <code>overlaps_any</code>, <code>merge</code>, <code>index</code>). Users pin a backend globally, for a block of code, or per call; <code>auto</code> picks by operation and size. Third-party packages can register engines through an entry point, so genomeblocks never has to import them up front.</p></div>
  <pre>import genomeblocks as gb

gb.backends.available()        # {{'intervals': ['numpy', 'cgranges', 'pyranges'], 'motifs': ['lightmotif', 'moods', 'fimo'], ...}}
gb.set_backend(motifs="moods")                        # global default
with gb.backend(intervals="pyranges"):                # scoped; contextvars, so thread-safe
    active = atac &amp; chip
m = cre.scan_motifs_matrix(genome, "jaspar.txt", backend="fimo")   # one call

# out to other tools without copying
cre.to_polars();  cre.to_arrow();  cre.to_pyranges();  Loci.from_frame(df, chrom="Chromosome")</pre>
</section>

<section id="intervals">
  <span class="tag green">Intervals</span>
  <h2>cgranges, pyranges, or our own?</h2>
  <p class="verdict">Both, by job: a numpy engine of about a hundred lines for whole-set operations (A &amp; B in <b>{ftime(and_col)}</b> at 1M vs {ftime(and_today)} today and {ftime(and_pr)} in pyranges), and cgranges kept for single lookups that return hits ({ftime(look['cgranges'])} each). pyranges becomes an adapter, not a dependency.</p>
  <p class="sub">The numpy engine places every chromosome on one genome-wide axis, sorts B once, keeps a running maximum of its ends, and answers each A interval with one <code>searchsorted</code>: an interval overlaps something in B exactly when the largest end among B intervals starting before its end lies past its start. Merge is one sort plus the same running maximum. Results are identical to today's{'' if all_same else ' (CHECK: some runs differ)'}.</p>
  {chart("iv_and", "A & B", "· log–log, lower is better", "Converting objects on every call is slower than today at 100k and 1M; the same engine on columns is the fastest at every size.")}
  {chart("iv_merge", "Sort + merge", "· log–log, lower is better", f"Today's merge takes {ftime(m_today)} for 2M intervals; numpy on columns takes {ftime(m_col)}, pyranges {ftime(m_pr)}.")}
  <div class="two">
    <figure><div class="ttl">What a backend choice costs <span>· per call</span></div>
      {table(["decision", "time"], [["auto policy (operation + size)", ftime(disp['auto policy (size + op check)'])], ["user-pinned backend", ftime(disp['user-pinned backend'])]], num_cols=(1,))}
      <figcaption>Negligible next to any operation; even A &amp; B on 1k peaks takes about 1 ms.</figcaption></figure>
    <figure><div class="ttl">One lookup <span>· against 100k peaks</span></div>
      {table(["engine", "per query"], [[k, ftime(v)] for k, v in look.items()], num_cols=(1,))}
      <figcaption>cgranges stays the lookup engine: fastest, and it returns the overlapping intervals, not just yes/no.</figcaption></figure>
  </div>
  {details("Conversion costs", table(["conversion", "n", "time"], conv_rows, num_cols=(1, 2)))}
</section>

<section id="genes">
  <span class="tag navy">Genes</span>
  <h2>Would fixing sort and merge make Genes.make faster?</h2>
  <p class="verdict">No. <code>Genes.make</code> never sorts or merges; its {ftime(cur)} go to splitting lines and building about a million objects. The sort/merge fix speeds up the <em>first annotation call</em>, which builds the promoter/exon/UTR index lazily: {ftime(a_cur)} → {ftime(a_key)}.</p>
  <p class="sub">Breakdown: reading the file takes {ftime(comp['read lines'])}, splitting every line {ftime(comp['read + split'])}, attribute parsing brings it to {ftime(comp['read + split + _parse_attributes'])}, and object creation (a Gene, Transcript and Exon/CDS/UTR dataclass per line, plus three Loci per transcript) the rest. A faster pure-Python parser does not help ({ftime(parse['fast Python parser, same objects']['seconds'])}), and pausing the garbage collector saves only {1 - parse['Genes.make, garbage collector paused']['seconds'] / cur:.0%}. What does: parse into a table with polars ({ftime(pol)}, Rust CSV reader and regex), build the annotation index straight from the table ({ftime(a_tab)}), and create Gene objects only when one is accessed. Building every object from the table still costs {ftime(pol_obj)}, which is the price of the object model.</p>
  {chart("genes_parse", "Parsing an 877k-line GTF", "· same file", "pyranges.read_gtf and pandas are not faster than today. The polars table holds every attribute genomeblocks uses; 'same objects' rows are checked to produce identical Gene/Transcript/Exon trees.")}
  {chart("genes_annot", "Building the annotation index", "· promoters, exons, 5′/3′ UTRs; identical output", "This runs once, on the first annotations() / nearest_genes() call.")}
</section>

<section id="motifs">
  <span class="tag green">Motifs</span>
  <h2>MOODS as a backend, or a MOODS-like scanner?</h2>
  <p class="verdict">MOODS as a backend: yes. It produced the identical 1,019-motif count matrix in <b>{ftime(mo1)}</b> on one core vs {ftime(lm1)} today ({fx(lm1 / mo1)}); with 4 workers {ftime(mo4)} vs {ftime(lm4)}.</p>
  <p class="sub">Writing our own: I built a MOODS-style lookahead scanner in numba (columns scored most-informative first, positions abandoned once the threshold is out of reach). It returns the identical matrix and skips {1 - prune:.0%} of the column work, yet takes {ftime(la1)} on one core ({ftime(la4)} on 4). A full-width SIMD scan (lightmotif) beats branchy pruning in a simple compiled loop; MOODS wins because its filter checks all motifs at once from one window hash per position, which would mean reimplementing MOODS. Wrapping it is the better trade. MOODS is pip-installable (<code>MOODS-python</code>) and dual-licensed GPLv3 / Biopython License, so it can be an optional extra of an MIT package. FIMO stays the choice when you need p-values and MEME-suite compatibility.</p>
  {chart("motif_backends", "Backends, same task", "· 5,000 × 500 bp windows, 1,019 JASPAR motifs, identical output", "The 4-worker MOODS run spends most of its time starting processes and shipping sequences at this size; one core is already the fast path.")}
</section>

<section id="interop">
  <span class="tag warn">Interop</span>
  <h2>Working with other tools today</h2>
  <p class="verdict">Inputs and outputs exist for the core formats, but five gaps stop genomeblocks from slotting in next to other tools; one of them silently shifts gene coordinates by 1 bp.</p>
  <div class="tbl" style="padding:0">{table(["tool / format", "today", "problem found", "proposed"], interop).replace('<div class="tbl">', '').replace('</div>', '')}</div>
  <p class="sub">The 1 bp GTF offset: genomeblocks stores GTF starts as written (1-based), pyranges converts them to 0-based like BED. The same gene reads 43,925,015 in genomeblocks and 43,925,014 in pyranges, so gene and exon starts are off by one against BED peaks. A related quirk: minus-strand TSS windows come out 1 bp shifted because the TSS is stored as <code>Locus(end, end - 1)</code>.</p>
</section>

<section id="plan">
  <span class="tag navy">Plan</span>
  <h2>A path that keeps the public API</h2>
  <div class="fixes">
    <div class="fix"><span class="t">1 · Fixes that need no redesign</span><p>Bounded scan in the pip fallback index ({fx(fb / bd)}), GFF3 attributes, gzip input, 0-based GTF starts, <code>to_bed()</code> returning text, numpy integer indexing, <code>signal()</code> honouring <code>workers</code>, pybigtools <code>fillna</code>, lazy plotting imports, key-based sort ({fx(a_cur / a_key)} on the annotation index).</p></div>
    <div class="fix"><span class="t">2 · Columnar Loci</span><p>Back <code>Loci</code> with columns (chrom codes, start, end, strand, extra BED columns) and create <code>Locus</code> objects on access. Indexing, iteration and set operators keep working; whole-set operations move to the numpy engine; pandas, polars, Arrow and pyranges become views. The main risk is code that mutates a <code>Locus</code> in place, which the 113-test suite can catch.</p></div>
    <div class="fix"><span class="t">3 · Backend registry</span><p><code>set_backend</code> / <code>backend()</code> / per-call <code>backend=</code> with an <code>auto</code> policy: intervals (numpy, cgranges), bigWig (pybigtools, pyBigWig, pure Python), motifs (lightmotif, MOODS, FIMO), annotation files (Python, polars), Hi-C (numpy, cooler). Optional extras such as <code>pip install genomeblocks[moods,polars]</code>.</p></div>
    <div class="fix"><span class="t">4 · Columnar Genes</span><p>Feature table as storage, annotation index built from it ({ftime(a_tab)}), Gene / Transcript objects on demand. Parse time drops from {ftime(cur)} toward {ftime(pol)}.</p></div>
    <div class="fix"><span class="t">5 · Adapters</span><p><code>to_/from_</code> for pandas, polars, Arrow, pyranges, bioframe and BedTool; AnnData peak names; an xarray view of the signal cube; edge lists for networkx / igraph.</p></div>
  </div>
</section>

<section id="method">
  <span class="tag navy">Method</span>
  <h2>How this was measured</h2>
  <ul class="plain">
    <li>Same machine and synthetic hg38-shaped data as the <a href="https://claude.ai/artifact/4hL67NF3y68wJSdsySs1Eo">benchmark report</a>; medians of 3–5 runs after a warm-up.</li>
    <li>Every prototype is checked for identical output before its time counts: same intervals, same merged blocks, same Gene/Transcript/Exon trees, same annotation index, same motif count matrix.</li>
    <li>Prototypes live in <code>benchmarks/prototypes/</code> (numpy interval engine, backend registry, numba lookahead scanner); benches in <code>bench_backends.py</code>, <code>bench_genes_make.py</code>, <code>bench_motif_backends.py</code>.</li>
  </ul>
</section>"""

    nav = "".join(f'<a href="#{i}">{t}</a>' for i, t in (
        ("answers", "Answers"), ("idea", "The idea"), ("intervals", "Intervals"), ("genes", "Genes"),
        ("motifs", "Motifs"), ("interop", "Interop"), ("plan", "Plan"), ("method", "Method")))
    page = f"""<title>genomeblocks Backends</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,400;0,600;0,700;0,800;1,400&display=swap">
<style>{(REP / 'style.css').read_text()}</style>
<div class="wrap">
<header class="top">
  <span class="eyebrow">genomeblocks · design proposal · measured prototypes</span>
  <h1>Interchangeable backends for genomeblocks</h1>
  <p class="lede">Can genomeblocks run on cgranges or pyranges, MOODS or lightmotif, and talk to the rest of the Python genomics stack? Each answer below comes from a working prototype checked against today's output.</p>
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
  num: (x) => {{ for (const [d, u] of [[1e9, "G"], [1e6, "M"], [1e3, "k"]]) if (x >= d) {{ const v = x / d;
                 return (v >= 10 ? v.toFixed(0) : v.toFixed(1)).replace(/\\.0$/, "") + u; }}
                 return x >= 10 ? x.toFixed(0) : (+x.toFixed(1)).toString(); }},
  int: (x) => String(Math.round(x)),
}};
document.addEventListener("DOMContentLoaded", () => {{
  for (const [id, s] of Object.entries(DATA)) {{
    const host = document.getElementById("c-" + id); if (!host) continue;
    if (s.kind === "hbar") Charts.hbar(host, {{...s, fmt: F[s.fmt]}});
    else Charts.lines(host, {{...s, xfmt: F[s.xfmt], yfmt: F[s.yfmt]}});
  }}
}});
</script>
"""
    out = REP / "backends.html"
    out.write_text(page)
    print(f"[proposal] {out} ({len(page) / 1e3:.0f} kB)")


if __name__ == "__main__":
    build()
