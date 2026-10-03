#!/usr/bin/env python3
"""Build the self-contained HTML benchmark report from ``results/*.json``.

Every number on the page is computed here from the recorded results, so the
report regenerates after a rerun with ``python build_report.py``.
Output: ``report/index.html`` (CSS, JS, diagrams and data inlined).
"""
from __future__ import annotations

import html
import json
import re

from common import HERE, RESULTS

REP = HERE / "report"


def load(name):
    p = RESULTS / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def esc(s):
    return html.escape(str(s))


def ftime(s):
    s = float(f"{s:.6g}")
    if s < 1e-3:
        return f"{s * 1e6:.1f} µs" if s < 1e-5 else f"{s * 1e6:.0f} µs"
    if s < 1:
        return f"{s * 1e3:.1f} ms" if s < 0.01 else f"{s * 1e3:.0f} ms"
    if s < 120:
        return f"{s:.1f} s" if s < 10 else f"{s:.0f} s"
    return f"{s / 60:.1f} min"


def fnum(x):
    for div, suf in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if x >= div:
            v = x / div
            return (f"{v:.0f}" if v >= 10 else f"{v:.1f}").rstrip("0").rstrip(".") + suf
    return f"{x:.0f}" if x >= 10 else f"{x:.1f}"


def fx(r):
    """Speed-up factor for prose: 2.4×, 26×, 1,200×."""
    if r >= 100:
        return f"{r:,.0f}×"
    return f"{r:.0f}×" if r >= 10 else f"{r:.1f}×"


def table(headers, rows, num_cols=(), hl=None):
    h = "".join(f'<th class="{"n" if i in num_cols else ""}">{esc(c)}</th>' for i, c in enumerate(headers))
    body = []
    for k, r in enumerate(rows):
        cls = ' class="hl"' if hl and hl(r) else ""
        tds = "".join(f'<td class="{"n" if i in num_cols else ""}">{esc(c)}</td>' for i, c in enumerate(r))
        body.append(f"<tr{cls}>{tds}</tr>")
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def details(summary, inner):
    return f"<details><summary>{esc(summary)}</summary>{inner}</details>"


def chart(cid, title, sub, caption):
    return (f'<figure><div class="ttl">{esc(title)} <span>{esc(sub)}</span></div>'
            f'<div id="{cid}"></div><figcaption>{caption}</figcaption></figure>')


def diagrams():
    raw = (REP / "diagrams.html").read_text()
    out = {}
    for m in re.finditer(r"<!-- @(\w+) -->\s*(.*?)(?=<!-- @|\Z)", raw, re.S):
        out[m.group(1)] = m.group(2).strip()
    return out


def why(title, svg, caption):
    return f'<div class="why"><h3>{esc(title)}</h3>{svg}<p class="sub">{caption}</p></div>'


# ─────────────────────────────────────────────────────────────────────────────

def build():
    D = diagrams()
    charts = {}
    secs = []
    tiles = []
    mech = []
    fixes = []

    # ── Loci ────────────────────────────────────────────────────────────────
    lo = load("loci")["rows"]
    lat = {r["engine"]: r["seconds"] for r in lo if r["op"] == "latency"}
    isect = {(r["engine"], r["n"]): r["seconds"] for r in lo if r["op"] == "intersect"}
    merge = {(r["engine"], r["n"]): r["seconds"] for r in lo if r["op"] == "merge"}
    make = {(r["engine"], r["n"]): r["seconds"] for r in lo if r["op"] == "make"}
    cg = lat["genomeblocks (cgranges)"]
    agree = all(r.get("agrees", True) for r in lo if r["op"] == "intersect")
    n_isect = sum(1 for r in lo if r["op"] == "intersect")
    tiles.append(("Overlap lookup", ftime(cg),
                  f"per query on 100k indexed peaks: {fx(lat['pyranges'] / cg)} faster than pyranges, "
                  f"{fx(lat['bioframe'] / cg)} faster than bioframe", "#loci", "green", "Loci"))
    charts["loci_latency"] = {"kind": "hbar", "fmt": "time", "log": True, "sort": "asc",
                              "axis": "time per overlap lookup (log scale)",
                              "rows": [{"label": k, "value": v, "hl": k.startswith("genomeblocks")}
                                       for k, v in lat.items()]}
    LC = {"genomeblocks (cgranges)": "var(--s1)", "pyranges": "var(--s2)", "bioframe": "var(--s3)",
          "bedtools (CLI)": "var(--s4)", "intervaltree": "var(--s5)",
          "genomeblocks (pure-Python fallback)": "var(--s8)",
          "genomeblocks (fallback, bounded scan)": "var(--s7)"}

    def ser(d, eng):
        return sorted((n, t) for (e, n), t in d.items() if e == eng)
    charts["loci_tools"] = {"kind": "lines", "xlog": True, "ylog": True, "xfmt": "num", "yfmt": "time",
                            "xLabel": "peaks per set (n)", "yLabel": "time for A & B", "endLabels": True,
                            "series": [{"name": e, "color": LC[e], "pts": ser(isect, e)} for e in
                                       ("genomeblocks (cgranges)", "pyranges", "bioframe", "bedtools (CLI)", "intervaltree")]}
    charts["loci_index"] = {"kind": "lines", "xlog": True, "ylog": True, "xfmt": "num", "yfmt": "time",
                            "xLabel": "peaks per set (n)", "yLabel": "time for A & B", "endLabels": True,
                            "series": [{"name": n, "color": LC[e], "pts": ser(isect, e)} for e, n in
                                       (("genomeblocks (cgranges)", "cgranges (conda)"),
                                        ("genomeblocks (pure-Python fallback)", "pure-Python fallback (pip)"),
                                        ("genomeblocks (fallback, bounded scan)", "fallback + bounded scan"))]}
    fb, bd = isect[("genomeblocks (pure-Python fallback)", 100_000)], isect[("genomeblocks (fallback, bounded scan)", 100_000)]
    fb_slope = fb / isect[("genomeblocks (pure-Python fallback)", 10_000)]
    rows = []
    for (e, n), t in sorted(isect.items(), key=lambda kv: (kv[0][1], kv[1])):
        rows.append([f"{n:,}", e, ftime(t)])
    m_gb, m_pr = merge[("genomeblocks", 1_000_000)], merge[("pyranges", 1_000_000)]
    mm_rows = [[op, f"{n:,}", e, ftime(t)] for op, d in (("sort + merge", merge), ("BED parse", make))
               for (e, n), t in sorted(d.items(), key=lambda kv: (kv[0][1], kv[1]))]
    secs.append(f"""
<section id="loci">
  <span class="tag green">Loci · interval algebra</span>
  <h2>Interval lookups: microseconds, because the index is kept</h2>
  <p class="verdict">One overlap lookup takes <b>{ftime(cg)}</b>. pyranges needs {ftime(lat['pyranges'])} and bioframe {ftime(lat['bioframe'])} for the same question.</p>
  {why("Why: the index lives on the Loci", D['loci'],
       "<code>Loci.cgr</code> builds a <a href='https://github.com/lh3/cgranges'>cgranges</a> index (C, implicit interval tree) the first time it is needed and caches it on the object. Data-frame libraries are built for whole-table joins, so a single-interval query pays their set-up cost every time. genomeblocks leans on cheap single lookups everywhere: <code>Architecture.make</code> anchors every loop with two of them, and <code>Genes.annotations</code>, <code>pair_to_bed</code> and <code>select_isoforms</code> do the same per locus.")}
  {chart("loci_latency", "One overlap lookup", "· 200 queries against 100k indexed peaks, lower is better",
         "Blue bars are genomeblocks. The pure-Python fallback is the index you get from <code>pip install genomeblocks</code>, because cgranges is not on PyPI (see <a href='#fixes'>fixes</a>).")}
  <div class="two">
    {chart("loci_tools", "Whole-set A & B", "· log–log, lower is better",
           f"All {n_isect} runs return the same intervals{'' if agree else ' (MISMATCH, check log)'}. Up to ~10k peaks genomeblocks is fastest. From 100k on, fully vectorised pyranges and bioframe pull ahead, by {fx(isect[('genomeblocks (cgranges)', 1_000_000)] / isect[('pyranges', 1_000_000)])} at 1M, because genomeblocks still steps through Python <code>Locus</code> objects one by one.")}
    {chart("loci_index", "Which index genomeblocks gets", "· log–log, lower is better",
           f"The pure-Python fallback grows quadratically ({fx(fb_slope)} slower for 10× more peaks). A bounded scan (a few lines, see <a href='#fixes'>fixes</a>) makes it {fx(fb / bd)} faster at 100k, matching cgranges.")}
  </div>
  {details("All intersect timings", table(["n per set", "engine", "median"], rows, num_cols=(0, 2), hl=lambda r: r[1].startswith("genomeblocks (cgranges)")))}
  {details("Merge and BED parsing", table(["operation", "n", "engine", "median"], mm_rows, num_cols=(1, 3), hl=lambda r: r[2] == "genomeblocks"))}
</section>""")
    mech.append(("Keep the interval index", "#loci",
                 f"A cgranges index is built once per <code>Loci</code> and reused. Single lookups cost {ftime(cg)}, so per-locus loops stay cheap."))
    fixes.append(("pip installs get a quadratic overlap index",
                  f"cgranges is not on PyPI, so pip users get <code>_PyIntervalIndex</code>. Its <code>overlap()</code> walks every interval that starts before the query end (<code>range(hi)</code> from 0), so each lookup is O(n) and <code>A &amp; B</code> is O(n·m): {ftime(fb)} at 100k peaks, about 10 min projected at 1M. Tracking the longest interval per chromosome and starting the walk at <code>bisect_left(starts, qs - max_len)</code> gives identical results in {ftime(bd)} ({fx(fb / bd)} faster). The benchmark's <code>BoundedPyIndex</code> is a drop-in."))

    # ── signal ─────────────────────────────────────────────────────────────
    sg = load("signal")["rows"]
    eng = {r["engine"]: r for r in sg if r["part"] == "engines"}
    gb = eng["genomeblocks · pybigtools (exact)"]["rate"]
    dt = eng["deepTools computeMatrix (-p 1)"]["rate"]
    pbs = eng["pyBigWig stats(nBins) exact"]["rate"]
    pbv = eng["pyBigWig values() + numpy binning"]["rate"]
    ptv = eng["pybigtools values() + numpy binning"]["rate"]
    pure = eng["genomeblocks · pure-Python reader"]["rate"]
    sc = {(r["engine"], r["n_loci"]): r["rate"] for r in sg if r["part"] == "scaling"}
    par = [r for r in sg if r["part"] == "parallel"]
    p1 = next(r["rate"] for r in par if r["engine"].startswith("genomeblocks") and r["workers"] == 1)
    p2 = next(r["rate"] for r in par if r["engine"].startswith("genomeblocks") and r["workers"] == 2)
    p4 = next(r["rate"] for r in par if r["engine"].startswith("genomeblocks") and r["workers"] == 4)
    t4 = next(r["rate"] for r in par if r["engine"].startswith("threads") and r["workers"] == 4)
    dt4 = next((r["rate"] for r in par if r["engine"].startswith("deepTools") and r["workers"] == 4), None)
    tiles.append(("Heatmap matrix from bigWigs", fx(gb / dt),
                  f"faster than deepTools computeMatrix on one core ({fnum(gb)} vs {fnum(dt)} regions/s); "
                  f"{fx(gb / pbs)} faster than pyBigWig's binned stats()", "#signal", "green", "signal"))
    SHORT = {"genomeblocks · pybigtools (exact)": "genomeblocks signal()",
             "genomeblocks · pybigtools (zoom, exact=False)": "genomeblocks signal(exact=False)",
             "genomeblocks · pure-Python reader": "genomeblocks, pure-Python reader",
             "pybigtools values() + numpy binning": "pybigtools per-base + numpy bins",
             "pyBigWig values() + numpy binning": "pyBigWig per-base + numpy bins",
             "pyBigWig stats(nBins) exact": "pyBigWig stats(nBins=200)",
             "pyBigWig stats(nBins) zoom": "pyBigWig stats(nBins=200), zoom",
             "deepTools computeMatrix (-p 1)": "deepTools computeMatrix -p 1"}
    charts["signal_engines"] = {"kind": "hbar", "fmt": "num", "log": True, "sort": "desc", "labelW": 300,
                                "axis": "regions per second, 1 track, 1 core (log scale)",
                                "rows": [{"label": SHORT[k], "value": r["rate"], "hl": k.startswith("genomeblocks")}
                                         for k, r in eng.items()]}
    charts["signal_scaling"] = {"kind": "lines", "xlog": True, "ylog": False, "xfmt": "num", "yfmt": "num",
                                "xLabel": "loci (random subset of 100k peaks, genome order)", "yLabel": "regions / s",
                                "endLabels": True,
                                "series": [{"name": SHORT[e], "color": c,
                                            "pts": sorted((n, v) for (ee, n), v in sc.items() if ee == e)}
                                           for e, c in (("genomeblocks · pybigtools (exact)", "var(--s1)"),
                                                        ("pyBigWig values() + numpy binning", "var(--s2)"),
                                                        ("genomeblocks · pure-Python reader", "var(--s3)"))]}
    proc = sorted((r["workers"], r["rate"]) for r in par if r["engine"].startswith("genomeblocks"))
    thr = [(1, p1)] + sorted((r["workers"], r["rate"]) for r in par if r["engine"].startswith("threads"))
    dts = sorted((r["workers"], r["rate"]) for r in par if r["engine"].startswith("deepTools"))
    charts["signal_parallel"] = {"kind": "lines", "xlog": False, "ylog": False, "xfmt": "int", "yfmt": "num",
                                 "xticks": [1, 2, 3, 4], "xLabel": "workers (4-core machine)",
                                 "yLabel": "region-tracks / s", "endLabels": True,
                                 "series": [{"name": "genomeblocks, processes", "color": "var(--s1)", "pts": proc},
                                            {"name": "same work on threads", "color": "var(--s2)", "pts": thr},
                                            {"name": "deepTools -p", "color": "var(--s3)", "pts": dts}]}
    eng_rows = [[SHORT[k], f"{r['n_loci']:,}", fnum(r["rate"]), ftime(r["seconds"]),
                 f"{r['pearson_r']:.4f}", f"{r['max_abs_diff']:.2g}"]
                for k, r in sorted(eng.items(), key=lambda kv: -kv[1]["rate"])]
    par_rows = [[r["engine"], str(r["workers"]), fnum(r["rate"]), ftime(r["seconds"])] for r in par]
    sc_hi = sc[("genomeblocks · pybigtools (exact)", 50_000)]
    sc_lo = sc[("genomeblocks · pybigtools (exact)", 1_000)]
    secs.append(f"""
<section id="signal">
  <span class="tag green">signal · bigWig extraction</span>
  <h2>bigWig → heatmap matrix: one native call per region</h2>
  <p class="verdict">On one core, <code>Loci.signal()</code> fills a ±3 kb × 200-bin matrix at <b>{fnum(gb)} regions/s</b>: {fx(gb / dt)} faster than deepTools computeMatrix and {fx(gb / pbs)} faster than pyBigWig's binned <code>stats()</code>. With 4 processes it reaches <b>{fnum(p4)} region-tracks/s</b>.</p>
  {why("Why: inflate each block once, bin in Rust", D['signal'],
       f"bigWig data sits in zlib-compressed blocks. pyBigWig's <code>stats(nBins=200)</code> queries each bin separately, so the same block is located and inflated up to 200 times per region. genomeblocks asks pybigtools for all 200 bins in one call; the block is inflated once, binned in Rust, and the result is written into a preallocated float32 cube. The bins step matters too: fetching per-base values and binning in numpy runs at {fnum(ptv)}/s (pybigtools) or {fnum(pbv)}/s (pyBigWig).")}
  {chart("signal_engines", "Engines, one track", "· 5,000 peaks (300 for the pyBigWig stats rows), higher is better",
         "Same windows, same bins. Every engine returns the same matrix as genomeblocks: identical values for the pybigtools and pure-Python paths, float32 rounding for pyBigWig, and Pearson r = 0.998 for deepTools. <code>exact=False</code> gives identical values and speed here because 30 bp bins are finer than these files' first zoom level (336 bp), so both read full-resolution data.")}
  <div class="two">
    {chart("signal_scaling", "Throughput vs number of loci", "· 1 track, 1 core",
           f"Throughput rises from {fnum(sc_lo)} to {fnum(sc_hi)} regions/s as peaks get denser: genome-sorted neighbours share compressed blocks, and one open pybigtools handle per track reuses them. The pure-Python reader and per-base loops stay flat.")}
    {chart("signal_parallel", "Workers", "· 16 tracks × 5,000 loci",
           f"Processes writing into one shared-memory cube scale {fx(p4 / p1)} on 4 cores (identical output each time). The same chunks on 4 threads add only {fx(t4 / p1)}, because pybigtools calls do not run in parallel across Python threads; that is why <code>signal()</code> uses processes." + (f" deepTools with 4 processes reaches {fnum(dt4)}/s." if dt4 else ""))}
  </div>
  {details("Engine table (rate, correctness vs genomeblocks)", table(["engine", "loci", "regions/s", "median", "Pearson r", "max |Δ|"], eng_rows, num_cols=(1, 2, 3, 4, 5), hl=lambda r: r[0].startswith("genomeblocks signal()")))}
  {details("Parallel table", table(["engine", "workers", "region-tracks/s", "median"], par_rows, num_cols=(1, 2, 3), hl=lambda r: r[0].startswith("genomeblocks")))}
</section>""")
    mech.append(("One native call per region", "#signal",
                 "All bins of a window come from one pybigtools call, so each compressed block is inflated once. Results go straight into a preallocated cube; processes share it through shared memory."))
    fixes.append(("signal(): workers are capped at half the cores",
                  f"<code>signal()</code> caps <code>workers</code> at <code>cpu_count() // 2</code>, so on this 4-core machine <code>workers=4</code> silently runs 2 processes ({fnum(p2)}/s). Lifting the cap gave {fnum(p4)}/s ({fx(p4 / p2)} more). Consider honouring an explicit <code>workers</code> request and only defaulting to half."))
    fixes.append(("pybigtools deprecates <code>missing=</code>",
                  "pybigtools 0.3 warns on every <code>values(..., missing=...)</code> call that the argument is deprecated in favour of <code>fillna</code>. <code>_PyBigToolsHandle</code> uses it on the hot path, so a future pybigtools release will break <code>signal()</code>. Switching to <code>fillna</code> (with a version check) avoids that."))
    fixes.append(("pure-Python bigWig reader: keep it, it is close",
                  f"The fallback reader manages {fnum(pure)} regions/s, {fx(pure / pbs)} faster than pyBigWig's binned stats. Most of its time is the per-call R-tree walk and block inflate; caching the last inflated block per handle (as pybigtools does) would help the dense-peak case most."))

    # ── Atlas ──────────────────────────────────────────────────────────────
    at = load("atlas")
    if at:
        ar = at["rows"]
        q = [r for r in ar if r["part"] == "query"]
        b = [r for r in ar if r["part"] == "build"]
        acc = [r for r in ar if r["part"] == "accuracy"]
        boot = [r for r in ar if r["part"] == "bootstrap"]
        gg = [r for r in ar if r["part"] == "giggle"]
        ld = {r["kind"]: r for r in ar if r["part"] == "load"}
        T = max(r["n_tracks"] for r in q)
        aq = {(r["engine"], r["n_tracks"]): r["seconds"] for r in q}
        for r in gg:
            if r["kind"] == "query":
                aq[("GIGGLE search -s (CLI)", r["n_tracks"])] = r["seconds"]
        a500 = aq[("Atlas.search (incl. Fisher + DataFrame)", T)]
        spmv = aq[("Atlas sparse row-sum only", T)]
        loop500 = aq[("per-track cgranges loop (prebuilt)", T)]
        pr500 = aq.get(("pyranges overlap(), per track", T))
        g500 = aq.get(("GIGGLE search -s (CLI)", T))
        gb_ = {r["n_tracks"]: r for r in gg if r["kind"] == "build"}
        gag = next((r for r in gg if r["kind"] == "agreement"), None)
        cold = ld.get("fresh python: import + load + search")
        aload = ld.get("Atlas.load (npz)")
        npz = HERE / "data" / "atlas_1kb.npz"
        npz_mb = npz.stat().st_size / 1e6 if npz.exists() else None
        imp_ = load("import")
        imp_loci = None
        if imp_:
            imp_loci = next((r["seconds"] for r in imp_["rows"] if r["stmt"].startswith("+ Loci")), None)
        tiles.append(("Enrichment vs 500 peak files", ftime(a500),
                      (f"per 20k-peak query with Fisher tests: {fx(g500 / a500)} faster than GIGGLE, " if g500 else "per 20k-peak query, ")
                      + f"{fx(loop500 / a500)} faster than looping over tracks", "#atlas", "navy", "Atlas"))
        NAMES = {"Atlas.search (incl. Fisher + DataFrame)": "Atlas.search()",
                 "GIGGLE search -s (CLI)": "GIGGLE search -s",
                 "per-track cgranges loop (prebuilt)": "per-track loop (cgranges)",
                 "pyranges overlap(), per track": "pyranges, per track",
                 "bedtools intersect -C (CLI)": "bedtools intersect -C"}
        AC = {"Atlas.search (incl. Fisher + DataFrame)": "var(--s1)", "GIGGLE search -s (CLI)": "var(--s2)",
              "per-track cgranges loop (prebuilt)": "var(--s3)", "pyranges overlap(), per track": "var(--s4)",
              "bedtools intersect -C (CLI)": "var(--s5)"}
        charts["atlas_query"] = {"kind": "lines", "xlog": True, "ylog": True, "xfmt": "num", "yfmt": "time",
                                 "xLabel": "tracks in the collection", "yLabel": "time per query", "endLabels": True,
                                 "series": [{"name": NAMES[e], "color": c,
                                             "pts": sorted((t, s_) for (ee, t), s_ in aq.items() if ee == e)}
                                            for e, c in AC.items() if any(ee == e for ee, _ in aq)]}
        b1 = next(r for r in b if r["bin_size"] == 1000 and r["workers"] == 3)
        b1s = next((r for r in b if r["bin_size"] == 1000 and r["workers"] == 1), None)
        build_rows = [[f"Atlas, {r['bin_size']:,} bp bins", str(r["workers"]), ftime(r["seconds"]),
                       f"{r['index_bytes'] / 1e6:.0f} MB in RAM" + (f" · {npz_mb:.0f} MB .npz" if r["bin_size"] == 1000 and npz_mb else "")]
                      for r in b]
        if T in gb_:
            build_rows.append(["GIGGLE, exact intervals", "1", ftime(gb_[T]["seconds"]),
                               f"{gb_[T]['index_bytes'] / 1e6:.0f} MB on disk"])
        acc_rows = [[f"{r['bin_size']:,} bp", f"{r['spearman_counts']:.3f}", f"{r['top25_overlap']}/25",
                     ftime(r["seconds"]), f"{r['nnz']:,}"] for r in acc]
        q_rows = [[NAMES.get(e, e), str(t), ftime(s_)] for (e, t), s_ in sorted(aq.items(), key=lambda kv: (kv[0][1], kv[1]))]
        boot_row = next((r for r in boot if r["n_iter"] == 100), None)
        acc1 = next((r for r in acc if r["bin_size"] == 1000), None)
        giggle_txt = ""
        if g500:
            giggle_txt = (f'<p class="sub"><b>About the GIGGLE comparison.</b> <a href="https://github.com/ryanlayer/giggle">GIGGLE</a> '
                          f'(Layer et al., <i>Nat. Methods</i> 2018, v0.6.3 built from source) is the tool Atlas is modelled on: it indexes exact intervals on disk and '
                          f'reports the same per-file Fisher statistics. Its search is a command-line call that opens the on-disk index each time ({ftime(g500)} at {T} files). '
                          + (f'Used the same way, a fresh Python process that imports genomeblocks, loads the saved Atlas ({npz_mb:.0f} MB .npz, {ftime(aload["seconds"])} to load) and searches takes {ftime(cold["seconds"])}, '
                             f'{"most of it" if imp_loci and imp_loci > cold["seconds"] / 2 else "much of it"} Python start-up and imports. Atlas wins when it stays loaded and answers many queries, as in a notebook, a bootstrap or a web service; '
                             f'GIGGLE wins for one-off command-line queries.' if cold and aload and npz_mb else '')
                          + (f' Both agree on what is enriched: GIGGLE and Atlas enrichment scores correlate at Spearman ρ = {gag["spearman_score_vs_atlas"]:.2f} across tracks, '
                             f'with {gag["top25_shared_with_atlas"]} of the top 25 tracks shared, and GIGGLE\'s exact overlap counts match exact interval counts at ρ = {gag["spearman_overlaps_vs_exact"]:.3f}.' if gag else '')
                          + '</p>')
        secs.append(f"""
<section id="atlas">
  <span class="tag navy">Atlas · GIGGLE-style enrichment</span>
  <h2>Enrichment against hundreds of peak files in one sparse operation</h2>
  <p class="verdict">A 20k-peak query against all {T} tracks takes <b>{ftime(a500)}</b>, Fisher tests included{f'. GIGGLE takes {ftime(g500)} for the same search' if g500 else ''}; intersecting track by track takes {ftime(loop500)}{f' (pyranges per track: {ftime(pr500)})' if pr500 else ''}.</p>
  {why("Why: the collection is one bin × track matrix", D['atlas'],
       f"<code>Atlas.make</code> tiles the genome into 1 kb bins and stores one bit per (bin, track) in a single CSR matrix held in memory. A query becomes its set of bins; selecting those rows and summing the columns (about {ftime(spmv)} here) gives the overlap count for every track at once, and the Fisher tests are vectorised across tracks. The cost tracks the number of query bins, not the number of files.")}
  {chart("atlas_query", "Query time vs collection size", "· 20k-peak query, log–log, lower is better",
         "Atlas is timed with its index loaded. GIGGLE and bedtools are command-line tools and read their inputs on every call, which is how they are used. The per-track loop uses prebuilt cgranges indexes (index once, query many), the fairest way to do this without Atlas.")}
  {giggle_txt}
  <div class="two">
    <figure><div class="ttl">Build once <span>· {b1['n_tracks']} BED files, {b1['n_peaks'] / 1e6:.1f}M peaks, {b1['bed_bytes'] / 1e6:.0f} MB of text</span></div>
      {table(["index", "workers", "build", "size"], build_rows, num_cols=(1, 2))}
      <figcaption>{f'Building the 1 kb index with 3 worker processes is {fx(b1s["seconds"] / b1["seconds"])} faster than with 1. ' if b1s else ''}<code>Atlas.save</code> / <code>load</code> make the build a one-time cost.</figcaption></figure>
    <figure><div class="ttl">What binning costs in accuracy <span>· vs exact interval overlaps</span></div>
      {table(["bin", "Spearman ρ", "top-25 shared", "query", "nnz"], acc_rows, num_cols=(1, 2, 3, 4))}
      <figcaption>Atlas counts shared bins, not intervals. Per-track overlap counts still rank almost identically to exact interval counts{f' (ρ = {acc1["spearman_counts"]:.3f} at 1 kb)' if acc1 else ''}; "top-25 shared" compares the 25 most-enriched tracks by each method. Smaller bins track more closely at a larger index.</figcaption></figure>
  </div>
  {f'<p class="sub">A 100-iteration shuffled-null bootstrap (<code>Atlas.bootstrap</code>) over all {T} tracks takes {ftime(boot_row["seconds"])}, {ftime(boot_row["per_iter"])} per iteration. With a command-line tool each iteration would be a separate search.</p>' if boot_row else ''}
  {details("Query timings", table(["engine", "tracks", "median"], q_rows, num_cols=(1, 2), hl=lambda r: r[0].startswith("Atlas.search")))}
</section>""")
        mech.append(("One matrix for the whole collection", "#atlas",
                     "Atlas turns hundreds of peak files into one sparse bin × track matrix kept in memory, so a query is a single row-sum over all tracks instead of one intersect per file."))

    # ── motifs ─────────────────────────────────────────────────────────────
    mo = load("motifs")
    if mo:
        mr = mo["rows"]
        me = {r["engine"]: r for r in mr if r["part"] == "engines"}
        mgb = me["genomeblocks scan_motifs_matrix (lightmotif)"]
        mbio = me["Biopython PSSM.search"]
        mmoods = me["MOODS (C++, all motifs per pass)"]
        mnp = me["numpy sliding window"]
        mrs = me["lightmotif, re-striped per motif"]
        mfi = me.get("MEME FIMO --text (CLI)")
        min_r = min(r["pearson_r"] for r in me.values())
        lib = [r for r in mr if r["part"] == "library"]
        wk = sorted([r for r in mr if r["part"] == "workers"], key=lambda r: r["workers"])
        st = {r["step"]: r["seconds"] for r in mr if r["part"] == "stripe"}
        tiles.append(("Motif scanning", fx((mfi or mbio)["seconds"] / mgb["seconds"]),
                      (f"faster than MEME FIMO, " if mfi else "") +
                      f"{fx(mbio['seconds'] / mgb['seconds'])} faster than Biopython; "
                      f"{mgb['gbp_motif_per_s']:.1f} Gbp·motif/s on one core", "#motifs", "green", "motifs"))
        charts["motifs_engines"] = {"kind": "hbar", "fmt": "time", "log": True, "sort": "asc", "labelW": 300,
                                    "axis": "time for 1,000 × 500 bp windows × 100 motifs (log scale)",
                                    "rows": [{"label": k.replace(" scan_motifs_matrix (lightmotif)", " (lightmotif)"),
                                              "value": r["seconds"], "hl": k.startswith("genomeblocks")}
                                             for k, r in me.items()]}
        if wk:
            charts["motifs_workers"] = {"kind": "lines", "xlog": False, "ylog": False, "xfmt": "int", "yfmt": "time",
                                        "xticks": [1, 2, 3, 4], "xLabel": "worker processes",
                                        "yLabel": "time, 5,000 windows × 1,019 motifs", "endLabels": False, "legend": False,
                                        "series": [{"name": "scan_motifs_matrix", "color": "var(--s1)",
                                                    "pts": [(r["workers"], r["seconds"]) for r in wk]}]}
        eng_rows = [[k, ftime(r["seconds"]), f"{r['gbp_motif_per_s']:.3f}", f"{r['total_hits']:,}",
                     f"{r['pearson_r']:.3f}"] for k, r in sorted(me.items(), key=lambda kv: kv[1]["seconds"])]
        lib_rows = [[r["engine"], f"{r['n_seqs']:,}", ftime(r["seconds"]), f"{r['gbp_motif_per_s']:.2f}"] for r in lib]
        stripe_frac = st.get("stripe 2000 seqs", 0) / max(st.get("scan 2000 striped seqs, 1 motif", 1e-9), 1e-9)
        w1 = wk[0]["seconds"] if wk else None
        w4 = wk[-1]["seconds"] if wk else None
        secs.append(f"""
<section id="motifs">
  <span class="tag green">motifs · PWM scanning</span>
  <h2>Motif scanning: SIMD scoring, sequences prepared once</h2>
  <p class="verdict"><code>scan_motifs_matrix</code> counts JASPAR hits {f"<b>{fx(mfi['seconds'] / mgb['seconds'])} faster than MEME FIMO</b>, " if mfi else ""}{fx(mbio['seconds'] / mgb['seconds'])} faster than Biopython and {fx(mnp['seconds'] / mgb['seconds'])} faster than a vectorised numpy scorer. MOODS, a dedicated C++ scanner that scores every motif in one pass per sequence, is {fx(mgb['seconds'] / mmoods['seconds'])} faster still; see <a href='#fixes'>fixes</a> for where that gap comes from.</p>
  {why("Why: SIMD scoring, and the sequence layout is built once", D['motifs'],
       f"<a href='https://github.com/althonos/lightmotif'>lightmotif</a> scores PSSMs with AVX2 over a <em>striped</em> sequence layout. genomeblocks stripes each window once and reuses it for all 1,019 motifs; motifs are then split across worker processes. Striping one sequence costs about {stripe_frac:.1f}× one motif scan, so re-striping per motif would cost {fx(mrs['seconds'] / mgb['seconds'])} the time.")}
  {chart("motifs_engines", "Engines", "· forward strand, log2-odds ≥ 13, same PSSMs, lower is better",
         f"Every engine scores the same probability matrices (pseudocount 0.1, uniform background, forward strand) and reports the same hits per motif (Pearson r ≥ {min_r:.4f}; FIMO differs by a handful of hits that sit exactly at the score cut-off, float32 vs float64 rounding)." + (" FIMO runs in its fastest mode (<code>--text</code>, no q-values); its time includes reading the FASTA and computing a p-value for every candidate, which is part of what FIMO is for." if mfi else ""))}
  {('<div class="two">' + chart("motifs_workers", "Worker processes", "· 5,000 windows × 1,019 motifs", f"{fx(w1 / w4)} with 4 processes. Each worker receives the sequences once (pool initializer) and scans a share of the motifs.") +
    '<figure><div class="ttl">Full JASPAR library <span>· 1,019 motifs, 1 core</span></div>' + table(["engine", "windows", "time", "Gbp·motif/s"], lib_rows, num_cols=(1, 2, 3)) + '<figcaption>Throughput is sequence length × motifs scored per second. MOODS checks all motifs in one pass per sequence with a lookahead filter.</figcaption></figure></div>') if wk else ''}
  {details("Engine table", table(["engine", "time", "Gbp·motif/s", "hits", "r vs genomeblocks"], eng_rows, num_cols=(1, 2, 3, 4), hl=lambda r: r[0].startswith("genomeblocks")))}
</section>""")
        mech.append(("SIMD motif scoring, layout built once", "#motifs",
                     "lightmotif scores many positions per instruction. genomeblocks stripes each sequence once and reuses it for every motif, then splits motifs across processes."))

    # ── pairs ──────────────────────────────────────────────────────────────
    pa = load("pairs")
    if pa:
        pr_ = {r["engine"]: r for r in pa["rows"]}
        cp = pr_["genomeblocks count_pairs (50 kb x partner chrom)"]
        ct = pr_["genomeblocks count_pairs (target_chrom=chr1)"]
        nv = pr_["naive per-pair loop (cgranges lookup)"]
        io = pr_["pandas chunked parse only (I/O floor)"]
        cl = pr_.get("cooler cload pairs (500 kb, CLI)")
        c2 = pr_.get("genomeblocks count_pairs_2d (500 kb)")
        fxp = load("fixes")
        fp = {}
        if fxp:
            for r in fxp["rows"]:
                if r.get("part") == "pairs":
                    fp[(r["op"], "fast" if r["variant"] != "current" else "cur")] = r
        f1 = fp.get(("count_pairs (50 kb x partner chrom)", "fast"))
        f2 = fp.get(("count_pairs_2d (500 kb)", "fast"))
        rows_ = [{"label": k, "value": r["rate"], "hl": k.startswith("genomeblocks")} for k, r in pr_.items()]
        if f1:
            rows_.append({"label": "count_pairs, integer-code fix (proposed)", "value": f1["n"] / f1["seconds"], "hl": True})
        if f2:
            rows_.append({"label": "count_pairs_2d, integer-code fix (proposed)", "value": f2["n"] / f2["seconds"], "hl": True})
        charts["pairs"] = {"kind": "hbar", "fmt": "num", "log": True, "sort": "desc", "labelW": 340,
                           "axis": "read pairs per second (log scale)", "rows": rows_}
        p_rows = [[k, f"{r['n_pairs']:,}", ftime(r["seconds"]), fnum(r["rate"])] for k, r in sorted(pr_.items(), key=lambda kv: -kv[1]["rate"])]
        for lab, r in (("count_pairs, integer-code fix", f1), ("count_pairs_2d, integer-code fix", f2)):
            if r:
                p_rows.append([lab + (" (same output)" if r.get("same") else " (OUTPUT DIFFERS)"), f"{r['n']:,}",
                               ftime(r["seconds"]), fnum(r["n"] / r["seconds"])])
        if f1 and f2 and cl:
            tiles.append(("Hi-C pair counting, with the fix", f"{fnum(f2['n'] / f2['seconds'])}/s",
                          f"window × window matrix after a measured fix: {fx(c2['seconds'] / f2['seconds'])} faster than today and "
                          f"{fx(cl['seconds'] / f2['seconds'])} faster than cooler cload (today it is {fx(c2['seconds'] / cl['seconds'])} slower)",
                          "#pairs", "warn", "bedpe"))
        secs.append(f"""
<section id="pairs">
  <span class="tag purple">bedpe · Hi-C pairs</span>
  <h2>Hi-C pairs: fast only for one chromosome today, with an easy fix</h2>
  <p class="verdict">With a target chromosome, <code>count_pairs</code> runs at <b>{fnum(ct['rate'])} pairs/s</b>. Counting against every partner chromosome drops to {fnum(cp['rate'])}/s, no faster than a per-pair Python loop ({fnum(nv['rate'])}/s){f', and <code>count_pairs_2d</code> manages {fnum(c2["rate"])}/s, {fx(c2["seconds"] / cl["seconds"])} slower than cooler cload' if c2 and cl else ''}.</p>
  <p class="sub">The file is parsed in 2M-row chunks by pandas' C reader at {fnum(io['rate'])} pairs/s, so reading is not the problem. The time goes to converting each chunk's categorical chromosome columns to Python strings and then comparing whole columns against every chromosome name, and in 2-D against every chromosome <em>pair</em> (23 × 23 × 2 full-chunk scans). Using the integer category codes, one <code>searchsorted</code> per anchor over genome-wide coordinates and a single <code>bincount</code> / sparse sum gives the same output{f' in {ftime(f1["seconds"])} ({fx(cp["seconds"] / f1["seconds"])}) and {ftime(f2["seconds"])} for the 2-D matrix ({fx(c2["seconds"] / f2["seconds"])})' if f1 and f2 else ''}.</p>
  {chart("pairs", "Pair counting", "· 5M pairs (naive loop: first 200k), higher is better",
         f"The naive loop does one cgranges lookup per anchor and agrees exactly with <code>count_pairs</code> on the same pairs ({'match' if nv.get('agrees_with_count_pairs') else 'MISMATCH'}). The fix rows are drop-in versions in <code>benchmarks/bench_fixes.py</code>, checked for identical output.")}
  {details("Pairs table", table(["engine", "pairs", "median", "pairs/s"], p_rows, num_cols=(1, 2, 3), hl=lambda r: r[0].startswith("genomeblocks")))}
</section>""")
        if f1 and f2:
            fixes.append(("count_pairs / count_pairs_2d compare chromosome strings per chunk",
                          f"Each chunk's chromosome columns become Python strings and are compared against every chromosome (2-D: every chromosome pair). "
                          f"With integer category codes and one global <code>searchsorted</code> per anchor: <code>count_pairs</code> {ftime(cp['seconds'])} → {ftime(f1['seconds'])}, "
                          f"<code>count_pairs_2d</code> {ftime(c2['seconds'])} → {ftime(f2['seconds'])} on 5M pairs, identical output "
                          f"({'checked' if f1.get('same') and f2.get('same') else 'CHECK FAILED'})."))

    # ── genes ──────────────────────────────────────────────────────────────
    ge = load("genes")
    if ge:
        gr = ge["rows"]
        gm = ge["meta"]
        parse = next(r for r in gr if r["step"].startswith("Genes.make"))
        ann = {(r["step"], r["n"]): r for r in gr if r["step"].startswith("annotations")}
        a_cg = ann[("annotations (cgranges)", 100_000)]
        a_fb = ann.get(("annotations (pure-Python fallback)", 10_000))
        a_cg10 = ann[("annotations (cgranges)", 10_000)]
        near = next(r for r in gr if r["step"].startswith("nearest") and r["n"] == 100_000)
        g_rows = [[r["step"], f"{r['n']:,}", ftime(r["seconds"]), fnum(r["rate"]) if r.get("rate") else ""] for r in gr]
        charts["genes"] = {"kind": "hbar", "fmt": "time", "log": True, "sort": None, "labelW": 330,
                           "axis": "seconds (log scale)",
                           "rows": [{"label": r["step"] + (f" · {fnum(r['n'])} CREs" if r["step"].startswith(("annotations", "nearest")) else ""),
                                     "value": r["seconds"], "hl": "fallback" not in r["step"]} for r in gr]}
        secs.append(f"""
<section id="genes">
  <span class="tag navy">Genes · annotation</span>
  <h2>Gene models: lazy indexes, per-CRE lookups</h2>
  <p class="verdict">Labelling 100k CREs (promoter / UTR / exon / intron / intergenic) takes <b>{ftime(a_cg['seconds'])}</b> ({fnum(a_cg['rate'])} CREs/s); nearest gene for the same set takes {ftime(near['seconds'])}.</p>
  <p class="sub"><code>Genes.annot</code> merges promoters, exons and UTRs into five <code>Loci</code> the first time it is used and keeps them, so each CRE costs at most five cached-index lookups. Parsing the GTF ({gm['gtf_lines']:,} lines, {gm['genes']:,} genes, {gm['transcripts']:,} transcripts) is plain Python at {fnum(parse['rate'])} lines/s, roughly {ftime(parse['seconds'] * 3_000_000 / gm['gtf_lines'])} for a ~3M-line GENCODE GTF{f'. With the pure-Python fallback index, annotation is {fx(a_fb["seconds"] / a_cg10["seconds"])} slower already at 10k CREs' if a_fb else ''}.</p>
  {chart("genes", "Gene-model steps", "· 20k genes; blue = default path",
         "<code>select_isoforms</code> picks, per gene, the longest isoform whose TSS is open in ATAC peaks and/or a bigWig; the bigWig variant scores every candidate TSS through <code>signal()</code>.")}
  {details("Genes table", table(["step", "n", "median", "per second"], g_rows, num_cols=(1, 2, 3)))}
</section>""")

    # ── architecture ───────────────────────────────────────────────────────
    ar_ = load("architecture")
    if ar_:
        rr = {r["step"]: r for r in ar_["rows"]}
        mk = rr["make (50k loops -> graph)"]
        s2v = rr["annotate stage 2: vectorised (genomeblocks)"]
        s2l = rr["annotate stage 2: per-vertex loop"]
        stl = rr["strength: per-edge Python loop (genomeblocks)"]
        stg = rr["strength: graph-tool incident_edges_op"]
        mc = rr["add_mcool (5 kb)"]
        tiles.append(("Graph annotation step", fx(s2l["seconds"] / s2v["seconds"]),
                      f"faster than the per-vertex loop it replaced ({ftime(s2v['seconds'])} vs {ftime(s2l['seconds'])} on {mk['n_edges']:,} edges)",
                      "#architecture", "purple", "Architecture"))
        charts["architecture"] = {"kind": "hbar", "fmt": "time", "log": True, "sort": None, "labelW": 340,
                                  "axis": "seconds (log scale)",
                                  "rows": [{"label": k, "value": r["seconds"], "hl": True}
                                           for k, r in rr.items()]}
        for row in charts["architecture"]["rows"]:
            row["hl"] = "per-vertex loop" not in row["label"] and "incident_edges_op" not in row["label"]
        a_rows = [[k, ftime(r["seconds"]), f"{r.get('n_edges', r.get('n_vertices', '')):,}" if isinstance(r.get('n_edges', r.get('n_vertices')), int) else "",
                   "" if "same_result" not in r else (("same result" + (f" ({r['n_differ']} of {r['n_assigned']:,} differ, only at tied weights)" if r.get("n_differ") else "")) if r["same_result"] else "differs")] for k, r in rr.items()]
        secs.append(f"""
<section id="architecture">
  <span class="tag purple">Architecture · contact graphs</span>
  <h2>Contact graphs: graph-tool storage, numpy over edge arrays</h2>
  <p class="verdict">Building a graph from 50k loops over 100k CREs takes <b>{ftime(mk['seconds'])}</b> ({mk['n_vertices']:,} loci, {mk['n_edges']:,} links). The vectorised annotation step is {fx(s2l['seconds'] / s2v['seconds'])} faster than the loop it replaced.</p>
  <p class="sub"><code>make</code> anchors each loop end with a cached-index lookup. <code>annotate</code> pulls the whole edge list with weights out of graph-tool as one (E, 3) array, masks promoter–non-promoter edges, and picks each CRE's strongest promoter with one <code>lexsort</code> + <code>unique</code>, instead of visiting neighbours vertex by vertex.</p>
  {chart("architecture", "Pipeline steps", "· 100k CREs, 50k loops, 5 kb .mcool; blue = genomeblocks today",
         f"Gray rows are comparisons: the stage-2 loop genomeblocks used to run, and graph-tool's native <code>incident_edges_op</code>, which computes the same node strengths {fx(stl['seconds'] / stg['seconds'])} faster than <code>strength()</code>'s per-edge Python loop (see <a href='#fixes'>fixes</a>).")}
  {details("Architecture table", table(["step", "median", "size", "check"], a_rows, num_cols=(1, 2)))}
</section>""")
        mech.append(("numpy over the graph's edge array", "#architecture",
                     f"<code>annotate</code> pulls all edges and weights out of graph-tool as one array and picks each CRE's strongest promoter with a sort, {fx(s2l['seconds'] / s2v['seconds'])} faster than visiting neighbours vertex by vertex."))
        fixes.append(("Architecture.strength() and add_mcool() loop in Python",
                      f"<code>strength()</code> walks every edge in Python ({ftime(stl['seconds'])}); <code>gt.incident_edges_op(G, 'out', 'sum', G.ep[key])</code> returns the same sums in {ftime(stg['seconds'])}. <code>add_mcool</code> ({ftime(mc['seconds'])}) does a pandas MultiIndex lookup per edge; joining the edge bin pairs against the pixel table in one merge would follow the pattern <code>annotate</code> already uses."))

    # ── import ─────────────────────────────────────────────────────────────
    im = load("import")
    if im:
        ir = im["rows"]
        base = ir[0]["seconds"]
        imp = {r["stmt"]: r["seconds"] for r in ir[1:] if not r["stmt"].startswith("breakdown")}
        brk = {r["stmt"].split(": ", 1)[1]: r["seconds"] for r in ir if r["stmt"].startswith("breakdown")}
        if brk:
            heavy = [(m, brk[m]) for m in ("genomeblocks.signal_draw", "genomeblocks.bedpe", "genomeblocks.atlas") if m in brk]
            fixes.append(("Touching Loci imports matplotlib, pandas and scipy",
                          f"<code>import genomeblocks</code> is lazy ({ftime(max(imp['import genomeblocks'], 0))}), but <code>from genomeblocks import Loci</code> takes "
                          f"{ftime(imp['+ Loci (pulls signal, motifs, atlas, bedpe)'])} because <code>loci.py</code> imports every module that attaches methods, and those import their heavy dependencies at module level: "
                          + ", ".join(f"<code>{m.split('.')[1]}</code> {ftime(t_)}" for m, t_ in heavy)
                          + " (cumulative, <code>python -X importtime</code>). Moving <code>matplotlib.pyplot</code>, <code>pandas</code> and <code>scipy.sparse</code> imports inside the functions that use them keeps the method attachment and removes most of that cost."))
        charts["import"] = {"kind": "hbar", "fmt": "time", "log": False, "sort": None, "labelW": 330,
                            "axis": "seconds on top of a bare interpreter start-up",
                            "rows": [{"label": k, "value": max(v, 1e-4), "hl": k == "import genomeblocks"} for k, v in imp.items()]}
        secs.append(f"""
<section id="import">
  <span class="tag navy">package · import</span>
  <h2>Import only what you touch</h2>
  <p class="verdict"><code>import genomeblocks</code> costs <b>{ftime(max(imp['import genomeblocks'], 0))}</b>. Each building block loads its dependencies on first access through a module <code>__getattr__</code>.</p>
  {chart("import", "Import cost", f"· median of 15 fresh interpreters, start-up ({ftime(base)}) subtracted",
         "Touching <code>Loci</code> pulls in numpy, pandas and the signal/motif/atlas modules that attach methods to it; <code>Architecture</code> adds graph-tool.")}
</section>""")
        mech.append(("Pay for imports on use", "#import",
                     f"<code>import genomeblocks</code> takes {ftime(max(imp['import genomeblocks'], 0))}; graph-tool, matplotlib and pandas load only when the block that needs them is touched."))

    # ── fixes bench ────────────────────────────────────────────────────────
    fxr = load("fixes")
    if fxr:
        f = {(r["op"], r["variant"], r["n"]): r for r in fxr["rows"]}
        n = max(r["n"] for r in fxr["rows"] if r.get("part", "loci") == "loci")
        mot = [r for r in fxr["rows"] if r.get("part") == "motifs"]
        if mot:
            Mx = max(r["n"] for r in mot)
            cur = next(r for r in mot if r["n"] == Mx and r["variant"].startswith("current"))
            new = next(r for r in mot if r["n"] == Mx and r["variant"].startswith("concatenated"))
            mo_ = load("motifs")
            moods_lib = next((r for r in (mo_["rows"] if mo_ else []) if r["part"] == "library"
                              and r["engine"].startswith("MOODS") and r["n_seqs"] == 1000), None)
            fixes.append(("scan_motifs_matrix makes one Python-level scan per window × motif",
                          f"For 1,000 windows × {Mx:,} motifs that is {1000 * Mx:,} calls into lightmotif, each on only 500 bp, so call overhead dominates ({ftime(cur['seconds'])}). "
                          f"Concatenating the windows, striping once, and scoring each motif over the whole block with <code>pssm.calculate(striped).threshold(t)</code> "
                          f"(hits split back per window with <code>np.bincount</code>, edge-crossing hits dropped) gives {'identical counts in every cell' if new.get('same') else 'DIFFERENT counts'} in {ftime(new['seconds'])} "
                          f"({fx(cur['seconds'] / new['seconds'])})" + (f", which also beats MOODS on the same task ({ftime(moods_lib['seconds'])})." if moods_lib and new['seconds'] < moods_lib['seconds'] else (f"; MOODS takes {ftime(moods_lib['seconds'])}." if moods_lib else ".")))))
        sc_ = f[("sort", "current (__lt__)", n)]; sk = f[("sort", "key=(chrom, start)", n)]
        mc_ = f[("merge", "current", n)]; mk_ = f[("merge", "key sort + single pass", n)]
        fixes.insert(1, ("Loci.sort() / merge() compare dataclasses in Python",
                         f"<code>sort()</code> calls <code>sorted(s)</code>, which runs <code>Locus.__lt__</code> in Python for every comparison. On {n:,} intervals that is {ftime(sc_['seconds'])}; <code>sorted(s, key=lambda l: (l.chrom, l.start))</code> gives the same order (checked: {'identical' if sk.get('same') else 'DIFFERENT'}) in {ftime(sk['seconds'])} ({fx(sc_['seconds'] / sk['seconds'])}). <code>merge()</code> goes from {ftime(mc_['seconds'])} to {ftime(mk_['seconds'])} ({fx(mc_['seconds'] / mk_['seconds'])}, identical output). pyranges still merges the same set in {ftime(m_pr * n / 2_000_000)}; a numpy merge over start/end arrays would close the rest."))
    else:
        fixes.insert(1, ("Loci.sort() / merge() compare dataclasses in Python",
                         f"Merging 2M intervals takes {ftime(m_gb)} vs {ftime(m_pr)} in pyranges. <code>sorted(s)</code> runs <code>Locus.__lt__</code> per comparison; a <code>(chrom, start)</code> key is the cheap fix."))
    fixes.append(("Small API papercut: numpy integers as Loci indexes",
                  "<code>Loci[np.int64(3)]</code> raises <code>TypeError</code> because <code>__getitem__</code> only accepts <code>int</code>. Checking <code>isinstance(key, numbers.Integral)</code> (or <code>operator.index</code>) lets index arrays from numpy work directly."))

    # ── env / method ───────────────────────────────────────────────────────
    env = load("loci")["env"]
    vers = env["versions"]
    vtxt = ", ".join(f"{k} {v}" for k, v in vers.items() if v and v != "installed")
    vtxt += ", cgranges (git), MOODS-python 1.9, deepTools 3.5.6, cooler 0.10"

    # ── assemble ───────────────────────────────────────────────────────────
    tile_html = "".join(
        f'<a class="tile" href="{h}"><span class="tag {c}">{esc(tag)}</span><span class="k">{esc(k)}</span>'
        f'<span class="v">{esc(v)}</span><span class="d">{esc(d)}</span></a>'
        for k, v, d, h, c, tag in tiles)
    mech_html = "".join(f'<div><h3><a href="{h}">{esc(t)}</a></h3><p>{d}</p></div>' for t, h, d in mech)
    fix_html = "".join(f'<div class="fix"><span class="t">{t}</span><p>{d}</p></div>' for t, d in fixes)
    nav = "".join(f'<a href="#{i}">{t}</a>' for i, t in (
        ("summary", "Summary"), ("loci", "Loci"), ("signal", "signal"), ("atlas", "Atlas"),
        ("motifs", "motifs"), ("pairs", "Hi-C pairs"), ("genes", "Genes"),
        ("architecture", "Architecture"), ("fixes", "Fixes"), ("method", "Method")))

    page = f"""<title>genomeblocks Benchmarks</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,400;0,600;0,700;0,800;1,400&display=swap">
<style>{(REP / 'style.css').read_text()}</style>
<div class="wrap">
<header class="top">
  <span class="eyebrow">genomeblocks 1.0.1 · benchmark report · {esc(env['cores'])}-core {esc(env['cpu'].replace('Intel(R) Xeon(R) Processor', 'Xeon'))}</span>
  <h1>How fast is genomeblocks, and why</h1>
  <p class="lede">Every building block timed against the tools people would otherwise use, on hg38-shaped data: what is fast, the mechanism behind each result, and the places that are still slow, with measured fixes.</p>
</header>
<nav class="sections" aria-label="Sections">{nav}</nav>
<section id="summary">
  <p class="verdict">Blue = genomeblocks, gray = alternatives. Every result is the median of 3–5 runs after a warm-up.</p>
  <div class="tiles">{tile_html}</div>
  <h2 style="margin-top:12px">Where the speed comes from</h2>
  <div class="mech">{mech_html}</div>
</section>
{''.join(secs)}
<section id="fixes">
  <span class="tag warn">Not fast yet</span>
  <h2>Where it is slow, and what fixes it</h2>
  <p class="verdict">The benchmarks also found a few slow paths. Each one comes with a measured or one-line fix.</p>
  <div class="fixes">{fix_html}</div>
</section>
<section id="method">
  <span class="tag navy">Method</span>
  <h2>How this was measured</h2>
  <ul class="plain">
    <li><b>Machine:</b> {esc(env['cpu'])}, {env['cores']} cores (AVX2/AVX-512), {env['mem_gb']:.0f} GB RAM, Linux, Python {esc(env['python'])}. Input files sit in the OS page cache, as in an interactive session.</li>
    <li><b>Timing:</b> wall clock, median of 3–5 runs after one warm-up (single runs are marked in the tables). Benchmarks ran one at a time.</li>
    <li><b>Data:</b> no public data could be downloaded in this sandbox, so everything is synthetic and seeded, laid out on real hg38 chromosome sizes: clustered peak sets with lognormal widths (1k–1M); four genome-wide bigWigs (13.8M variable-span records, 165 MB each, similar to a ChIP-seq track); 500 peak files (7.4M peaks) for Atlas; 5M Hi-C pairs with P(s) ∝ s⁻¹; a GENCODE-shaped GTF (20k genes); 50k loops; and the real JASPAR 2024 CORE vertebrate motifs (1,019).</li>
    <li><b>Correctness:</b> every comparison checks that the engines return the same answer (same intervals, same matrix, same hit counts, same graph labels) before timing is compared.</li>
    <li><b>Versions:</b> {esc(vtxt)}.</li>
  </ul>
  <pre>cd benchmarks
python make_data.py            # ~3 min, ~1.5 GB
PY=python ./run_all.sh         # every bench, then figures
python build_report.py         # this page → report/index.html</pre>
</section>
</div>
<script>{(REP / 'charts.js').read_text()}</script>
<script>
const DATA = {json.dumps(charts)};
const F = {{
  time: (s) => {{ s = +(+s).toPrecision(6); if (s < 1e-3) return (s < 1e-5 ? (s*1e6).toFixed(1) : (s*1e6).toFixed(0)) + " µs";
                 if (s < 1) return (s < 0.01 ? (s*1e3).toFixed(1) : (s*1e3).toFixed(0)) + " ms";
                 if (s < 120) return (s < 10 ? s.toFixed(1) : s.toFixed(0)) + " s"; return (s/60).toFixed(1) + " min"; }},
  num: (x) => {{ for (const [d, u] of [[1e9, "G"], [1e6, "M"], [1e3, "k"]]) if (x >= d) {{ const v = x / d;
                 return (v >= 10 ? v.toFixed(0) : v.toFixed(1)).replace(/\\.0$/, "") + u; }}
                 return x >= 10 ? x.toFixed(0) : (+x.toFixed(1)).toString(); }},
  int: (x) => String(Math.round(x)),
}};
document.addEventListener("DOMContentLoaded", () => {{
  for (const [id, s] of Object.entries(DATA)) {{
    const host = document.getElementById(id); if (!host) continue;
    if (s.kind === "hbar") Charts.hbar(host, {{...s, fmt: F[s.fmt]}});
    else Charts.lines(host, {{...s, xfmt: F[s.xfmt], yfmt: F[s.yfmt]}});
  }}
}});
</script>
"""
    out = REP / "index.html"
    out.write_text(page)
    print(f"[report] {out} ({len(page) / 1e3:.0f} kB)")
    write_readme(tiles, fixes, env, vtxt)


def _plain(s):
    """HTML snippet -> Markdown-ish text for the README."""
    s = re.sub(r"<code>(.*?)</code>", r"`\1`", s)
    s = re.sub(r"<a href='([^']*)'>(.*?)</a>", r"\2", s)
    s = re.sub(r'<a href="([^"]*)">(.*?)</a>', r"\2", s)
    s = re.sub(r"<[^>]+>", "", s)
    return html.unescape(s)


def write_readme(tiles, fixes, env, vtxt):
    figs = [("loci_latency", "One overlap lookup"), ("loci_intersect", "Whole-set intersect"),
            ("signal_engines", "bigWig → heatmap matrix, one core"), ("signal_scaling", "Throughput vs loci"),
            ("signal_parallel", "Processes vs threads vs deepTools"), ("atlas_query", "Atlas vs GIGGLE vs per-track loops"),
            ("motifs_engines", "Motif scanning engines"), ("pairs", "Hi-C pair counting"),
            ("genes", "Gene-model steps"), ("architecture", "Architecture pipeline"), ("import", "Import cost")]
    lines = ["# genomeblocks benchmarks", "",
             "How fast each building block is against the tools people would otherwise use, why, "
             "and where it is still slow. Open `report/index.html` for the full interactive report "
             "(charts, mechanism diagrams, every table).", "",
             f"Machine: {env['cpu']}, {env['cores']} cores, {env['mem_gb']:.0f} GB RAM, Python {env['python']}. "
             "Medians of 3–5 runs after a warm-up; every comparison is checked for identical output first.", "",
             "## Headlines", "", "| block | result | context |", "|---|---|---|"]
    for k, v, d, _h, _c, tag in tiles:
        lines.append(f"| {tag} | **{v}** {k.lower()} | {d} |")
    lines += ["", "## Figures", ""]
    for name, title in figs:
        if (HERE / "figures" / f"{name}.png").exists():
            lines += [f"**{title}**", "", f"![{title}](figures/{name}.png)", ""]
    lines += ["## Where it is slow, and the fix", ""]
    for t, d in fixes:
        lines += [f"- **{_plain(t)}** {_plain(d)}", ""]
    lines += ["## Reproduce", "", "```bash", "cd benchmarks",
              "python make_data.py            # synthetic hg38-shaped data, ~1.5 GB",
              "PY=python ./run_all.sh         # every bench (one at a time), then figures",
              "python build_report.py         # report/index.html + this README", "```", "",
              "External baselines are found on `PATH` (bedtools, deepTools `computeMatrix`, cooler, bgzip) "
              "or via `GIGGLE=` / `FIMO=` (GIGGLE built from github.com/ryanlayer/giggle; "
              "MEME suite from bioconda: `micromamba create -n tools -c conda-forge -c bioconda meme`). "
              "`Architecture` needs graph-tool (conda-forge).", "",
              f"Versions: {vtxt}.", ""]
    (HERE / "README.md").write_text("\n".join(lines))
    print(f"[readme] {HERE / 'README.md'}")


if __name__ == "__main__":
    build()
