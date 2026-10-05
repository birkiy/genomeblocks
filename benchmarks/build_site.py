#!/usr/bin/env python3
"""Build the Benchmarks page of the docs site from results/*.json.

    python benchmarks/build_site.py      # writes docs/benchmarks/index.md

Charts are inline SVG styled by the site's tokens (``.gbc`` rules in
docs/_sass/custom/custom.scss). Form: timings span orders of magnitude, so
engine comparisons are dot plots on a log axis, with genomeblocks in the
module's colour and every other engine in grey (emphasis, not a categorical
palette). Every chart has a native hover title per mark and a table twin.
"""
from __future__ import annotations

import json
import math
from html import escape
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = HERE.parent / "docs" / "benchmarks" / "index.md"


# ── formatting ───────────────────────────────────────────────────────────────

def ftime(s: float) -> str:
    if s is None or not math.isfinite(s):
        return "–"
    if s < 1e-3:
        return f"{s * 1e6:.1f} µs" if s < 1e-5 else f"{s * 1e6:.0f} µs"
    if s < 1:
        return f"{s * 1e3:.1f} ms" if s < 0.01 else f"{s * 1e3:.0f} ms"
    if s < 120:
        return f"{s:.1f} s" if s < 10 else f"{s:.0f} s"
    return f"{s / 60:.1f} min"


def ftick(s: float) -> str:
    """Axis ticks at decades: 1 µs, 10 ms, 1 s, 1,000 s."""
    if s < 1e-3:
        return f"{s * 1e6:g} µs"
    if s < 1:
        return f"{s * 1e3:g} ms"
    return f"{s:,g} s"


def fnum(x: float) -> str:
    for d, u in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if x >= d:
            v = x / d
            t = f"{v:.0f}" if v >= 10 else f"{v:.1f}"
            return (t[:-2] if t.endswith(".0") else t) + u
    t = f"{x:.0f}" if x >= 10 else f"{x:.1f}"
    return t[:-2] if t.endswith(".0") else t


def fx(r: float) -> str:
    return f"{r:,.0f}×" if r >= 10 else f"{r:.1f}×"


def load(name):
    p = RES / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


# ── charts ───────────────────────────────────────────────────────────────────

def _log_ticks(lo, hi):
    a, b = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
    return [10.0 ** k for k in range(a, b + 1)]


def dotplot(rows, *, kind="green", unit="s", label="", width=880, left=300):
    """rows: (label, value, is_genomeblocks, hover_note). Log x axis."""
    rows = [r for r in rows if r[1] and r[1] > 0]
    rows.sort(key=lambda r: -r[1] if unit == "rate" else r[1])     # best first
    rh, top = 28, 8
    h = top + rh * len(rows) + 34
    vals = [r[1] for r in rows]
    ticks = _log_ticks(min(vals) / 1.6, max(vals) * 1.6)
    lo, hi = math.log10(ticks[0]), math.log10(ticks[-1])
    x0, x1 = left, width - 96
    X = lambda v: x0 + (math.log10(v) - lo) / (hi - lo) * (x1 - x0)
    fmt = ftime if unit == "s" else fnum
    tfmt = ftick if unit == "s" else fnum
    out = [f'<svg class="gbc" viewBox="0 0 {width} {h}" role="img" aria-label="{escape(label)}" '
           f'xmlns="http://www.w3.org/2000/svg">']
    for t in ticks:
        out.append(f'<line class="grid" x1="{X(t):.1f}" y1="{top}" x2="{X(t):.1f}" y2="{h - 26}"/>')
        out.append(f'<text class="tick" x="{X(t):.1f}" y="{h - 10}" text-anchor="middle">{tfmt(t)}</text>')
    for i, (lab, v, gb, note) in enumerate(rows):
        y = top + rh * i + rh / 2
        out.append(f'<g class="row{" gb " + kind if gb else ""}"><title>{escape(lab)}: {fmt(v)}'
                   f'{" · " + escape(note) if note else ""}</title>')
        out.append(f'<rect class="hit" x="0" y="{y - rh / 2:.1f}" width="{width}" height="{rh}"/>')
        out.append(f'<line class="rowline" x1="{x0}" y1="{y:.1f}" x2="{x1}" y2="{y:.1f}"/>')
        out.append(f'<text class="lab" x="{left - 14}" y="{y + 4:.1f}" text-anchor="end">{escape(lab)}</text>')
        out.append(f'<circle class="dot" cx="{X(v):.1f}" cy="{y:.1f}" r="5.5"/>')
        out.append(f'<text class="val" x="{X(v) + 11:.1f}" y="{y + 4:.1f}">{fmt(v)}</text>')
        out.append("</g>")
    out.append("</svg>")
    return "\n".join(out)


def linechart(series, *, kind="green", xlab="", ylab="", label="", width=880, h=300, xfmt=fnum, yfmt=fnum,
              ylog=True, ytick=None):
    """series: (name, [(x, y), ...], is_genomeblocks). Log x (and y) axes, end labels."""
    left, right, top, bottom = 64, 230, 14, 44
    xs = [x for _, pts, _ in series for x, _ in pts]
    ys = [y for _, pts, _ in series for _, y in pts if y > 0]
    xt = _log_ticks(min(xs), max(xs))
    yt = _log_ticks(min(ys) / 1.3, max(ys) * 1.3) if ylog else None
    lx0, lx1 = math.log10(xt[0]), math.log10(xt[-1])
    if ylog:
        ly0, ly1 = math.log10(yt[0]), math.log10(yt[-1])
    X = lambda v: left + (math.log10(v) - lx0) / (lx1 - lx0) * (width - left - right)
    Y = lambda v: (h - bottom) - (math.log10(v) - ly0) / (ly1 - ly0) * (h - bottom - top)
    out = [f'<svg class="gbc" viewBox="0 0 {width} {h}" role="img" aria-label="{escape(label)}" '
           f'xmlns="http://www.w3.org/2000/svg">']
    for t in yt:
        out.append(f'<line class="grid" x1="{left}" y1="{Y(t):.1f}" x2="{width - right}" y2="{Y(t):.1f}"/>')
        out.append(f'<text class="tick" x="{left - 8}" y="{Y(t) + 4:.1f}" text-anchor="end">{(ytick or yfmt)(t)}</text>')
    for t in xt:
        out.append(f'<text class="tick" x="{X(t):.1f}" y="{h - bottom + 18}" text-anchor="middle">{xfmt(t)}</text>')
    out.append(f'<line class="axis" x1="{left}" y1="{h - bottom}" x2="{width - right}" y2="{h - bottom}"/>')
    out.append(f'<text class="axl" x="{(left + width - right) / 2}" y="{h - 6}" text-anchor="middle">{escape(xlab)}</text>')
    ends = []
    order = sorted(series, key=lambda s: s[2])          # grey first, emphasis on top
    for name, pts, gb in order:
        pts = sorted(pts)
        d = " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in pts)
        cls = f"line gb {kind}" if gb else "line"
        out.append(f'<g class="{cls}"><title>{escape(name)}</title><polyline points="{d}"/>')
        for x, y in pts:
            out.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="4"><title>{escape(name)} · '
                       f'{xfmt(x)}: {yfmt(y)}</title></circle>')
        out.append("</g>")
        ends.append([Y(pts[-1][1]), X(pts[-1][0]), name, gb])
    ends.sort()
    for i in range(1, len(ends)):                     # keep end labels 14 px apart
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + 14)
    for y, x, name, gb in ends:
        out.append(f'<text class="lab{" gbl" if gb else ""}" x="{width - right + 10}" y="{y + 4:.1f}">{escape(name)}</text>')
        out.append(f'<line class="leader" x1="{x + 6:.1f}" y1="{y:.1f}" x2="{width - right + 6}" y2="{y:.1f}"/>')
    out.append("</svg>")
    return "\n".join(out)


def figure(svg, caption, table=None):
    t = ""
    if table:
        head, rows = table
        th = "".join(f"<th>{escape(h)}</th>" for h in head)
        tr = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
        t = f'<details class="gb-table"><summary>Table</summary><table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></details>'
    return (f'<figure class="gb-fig gb-chart"><div class="gb-fig-body">\n{svg}\n</div>'
            f'<figcaption>{caption}</figcaption>{t}</figure>\n')


def tiles(items):
    out = ['<div class="gb-tiles">']
    for label, value, context, kind in items:
        out.append(f'<div class="gb-tile {kind}"><span class="tile-label">{escape(label)}</span>'
                   f'<span class="tile-value">{escape(value)}</span><span class="tile-ctx">{context}</span></div>')
    out.append("</div>")
    return "\n".join(out)


def rows_where(d, **kw):
    return [r for r in d["rows"] if all(r.get(k) == v for k, v in kw.items())]


def pick(d, **kw):
    r = rows_where(d, **kw)
    return r[0] if r else None


def dumbbell(rows, *, kind="green", label="", width=880, left=250, a_name="classic", b_name="columnar"):
    """rows: (label, a_seconds, b_seconds). Grey dot = a, coloured dot = b, log x axis."""
    rows = [r for r in rows if r[1] and r[2]]
    rh, top = 28, 30
    h = top + rh * len(rows) + 34
    vals = [v for _, a, b in rows for v in (a, b)]
    ticks = _log_ticks(min(vals) / 1.6, max(vals) * 1.6)
    lo, hi = math.log10(ticks[0]), math.log10(ticks[-1])
    x0, x1 = left, width - 70
    X = lambda v: x0 + (math.log10(v) - lo) / (hi - lo) * (x1 - x0)
    out = [f'<svg class="gbc" viewBox="0 0 {width} {h}" role="img" aria-label="{escape(label)}" '
           f'xmlns="http://www.w3.org/2000/svg">']
    # legend
    out.append(f'<circle class="dot a" cx="{x0 + 6}" cy="12" r="5"/><text class="lab" x="{x0 + 16}" y="16">{a_name}</text>')
    out.append(f'<g class="gb {kind}"><circle class="dot" cx="{x0 + 110}" cy="12" r="5"/></g>'
               f'<text class="lab" x="{x0 + 120}" y="16">{b_name}</text>')
    for t in ticks:
        out.append(f'<line class="grid" x1="{X(t):.1f}" y1="{top}" x2="{X(t):.1f}" y2="{h - 26}"/>')
        out.append(f'<text class="tick" x="{X(t):.1f}" y="{h - 10}" text-anchor="middle">{ftick(t)}</text>')
    for i, (lab, a, b) in enumerate(rows):
        y = top + rh * i + rh / 2
        out.append(f'<g class="row gb {kind}"><title>{escape(lab)}: {a_name} {ftime(a)} → {b_name} {ftime(b)} '
                   f'({fx(a / b)} {"faster" if a > b else "slower"})</title>')
        out.append(f'<rect class="hit" x="0" y="{y - rh / 2:.1f}" width="{width}" height="{rh}"/>')
        out.append(f'<text class="lab" x="{left - 14}" y="{y + 4:.1f}" text-anchor="end">{escape(lab)}</text>')
        out.append(f'<line class="span" x1="{X(min(a, b)):.1f}" y1="{y:.1f}" x2="{X(max(a, b)):.1f}" y2="{y:.1f}"/>')
        out.append(f'<circle class="dot a" cx="{X(a):.1f}" cy="{y:.1f}" r="5"/>')
        out.append(f'<circle class="dot" cx="{X(b):.1f}" cy="{y:.1f}" r="5.5"/>')
        out.append(f'<text class="val" x="{x1 + 64}" y="{y + 4:.1f}" text-anchor="end">{fx(a / b)}</text>')
        out.append("</g>")
    out.append("</svg>")
    return "\n".join(out)


# ── engine labels ────────────────────────────────────────────────────────────

SHORT = {
    "pandas chunked parse only (I/O floor)": "parse only (I/O floor)",
    "genomeblocks count_pairs (50 kb x partner chrom)": "count_pairs · 50 kb × partner",
    "genomeblocks count_pairs (target_chrom=chr1)": "count_pairs · one partner",
    "genomeblocks count_pairs_2d (500 kb)": "count_pairs_2d · 500 kb",
    "cooler cload pairs (500 kb, CLI)": "cooler cload · 500 kb (CLI)",
    "naive per-pair loop (cgranges lookup)": "per-pair Python loop",
    "genomeblocks scan_motifs_matrix (lightmotif)": "genomeblocks scan_motifs_matrix",
    "MOODS (C++, all motifs per pass)": "MOODS (C++)",
    "MEME FIMO --text (CLI)": "MEME FIMO (CLI)",
    "Atlas.search (incl. Fisher + DataFrame)": "Atlas.search (with Fisher + table)",
    "make (50k loops -> graph)": "make (50k loops → graph)",
}


def short(name: str) -> str:
    return SHORT.get(name, name)


def is_gb(name: str) -> bool:
    return name.lower().startswith("genomeblocks") or name.startswith("Atlas")


# ── sections ─────────────────────────────────────────────────────────────────

def sec_loci(md, tl):
    d = load("loci")
    if not d:
        return
    dc = load("loci_columnar")
    rows = d["rows"] + (dc["rows"] if dc else [])
    lat = {r["engine"]: r["seconds"] for r in rows if r["op"] == "latency"}
    gbl = lat.get("genomeblocks (cgranges)")
    md.append('## Loci\n{: .sec-green #loci }\n')
    if gbl:
        tl.append(("One overlap lookup", ftime(gbl),
                   f"on 100k indexed peaks: {fx(lat['pyranges'] / gbl)} faster than pyranges, "
                   f"{fx(lat['bioframe'] / gbl)} than bioframe", "green"))
        md.append(f"A single overlap query against an indexed set is where the object model pays off: "
                  f"**{ftime(gbl)}** per lookup with cgranges and "
                  f"{ftime(lat['genomeblocks (pure-Python fallback)'])} with the pure-Python index, against "
                  f"{ftime(lat['pyranges'])} for pyranges, which has to slice a dataframe per call. "
                  f"That is the cost of every interactive question (\"what overlaps this peak?\").\n")
        md.append(figure(dotplot([(k, v, is_gb(k), "per query, 100k-peak index") for k, v in lat.items()],
                                 kind="green", label="Single overlap query latency"),
                         "<strong>One overlap lookup</strong> against 100,000 indexed peaks (median per query; log scale)."))
    big = max(r["n"] for r in rows if r["op"] == "intersect")
    for op, title, cap in (("intersect", "A & B", "rows of A overlapping B"), ("merge", "merge", "sort + merge of A ∪ B")):
        sel = [r for r in rows if r["op"] == op and r["n"] == big and "warm" not in r["engine"]]
        if not sel:
            continue
        t = {r["engine"]: r["seconds"] for r in sel}
        md.append(figure(dotplot([(k, v, is_gb(k), f"n = {big:,}") for k, v in t.items()], kind="green",
                                 label=f"{title} at {big:,} intervals"),
                         f"<strong>{escape(title)}</strong> on {fnum(big)} intervals per set ({cap}); every engine "
                         f"returns the same rows."))
    # verdict on bulk ops
    t_i = {r["engine"]: r["seconds"] for r in rows if r["op"] == "intersect" and r["n"] == big}
    t_m = {r["engine"]: r["seconds"] for r in rows if r["op"] == "merge" and r["n"] == big}
    if "genomeblocks.columnar" in t_i:
        md.append(f"On whole-set operations at {fnum(big)} intervals the object-per-interval `Loci` is not the "
                  f"fastest: `A & B` takes {ftime(t_i['genomeblocks (cgranges)'])} and `merge` "
                  f"{ftime(t_m['genomeblocks'])}, against {ftime(t_i['pyranges'])} and {ftime(t_m['pyranges'])} "
                  f"for pyranges. The [columnar `Loci`]({{{{ '/design/columnar/' | relative_url }}}}) does the same work on numpy columns: "
                  f"{ftime(t_i['genomeblocks.columnar'])} and {ftime(t_m['genomeblocks.columnar'])}.\n")
    sizes = sorted({r["n"] for r in rows if r["op"] == "intersect"})
    engines = []
    for r in rows:
        if r["op"] in ("intersect", "merge", "make") and r["engine"] not in engines:
            engines.append(r["engine"])
    trows = []
    for op in ("intersect", "merge", "make"):
        for e in engines:
            vals = {r["n"]: r["seconds"] for r in rows if r["op"] == op and r["engine"] == e}
            if vals:
                trows.append([op, escape(e)] + [ftime(vals.get(n)) if n in vals else "–" for n in sizes])
    md.append('<details class="gb-table"><summary>All Loci timings</summary><table><thead><tr><th>op</th><th>engine</th>'
              + "".join(f"<th>n = {fnum(n)}</th>" for n in sizes) + "</tr></thead><tbody>"
              + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in trows)
              + "</tbody></table></details>\n")


def sec_signal(md, tl):
    d = load("signal")
    if not d:
        return
    md.append('## Signal\n{: .sec-green #signal }\n')
    eng = rows_where(d, part="engines")
    if eng:
        gb = pick(d, part="engines", engine="genomeblocks · pybigtools (exact)")
        dt = next((r for r in eng if r["engine"].startswith("deepTools")), None)
        pw = next((r for r in eng if r["engine"].startswith("pyBigWig stats") and "exact" in r["engine"]), None)
        if gb and dt:
            tl.append(("bigWig → heatmap matrix", fx(gb["rate"] / dt["rate"]),
                       f"faster than deepTools computeMatrix on one core ({fnum(gb['rate'])} vs "
                       f"{fnum(dt['rate'])} regions/s)", "green"))
            md.append(f"Filling a heatmap matrix (200 bins, ±3 kb) is one native call per region and track: "
                      f"**{fnum(gb['rate'])} regions/s** on one core, {fx(gb['rate'] / dt['rate'])} deepTools "
                      f"`computeMatrix`" + (f" and {fx(gb['rate'] / pw['rate'])} pyBigWig's binned `stats()`" if pw else "")
                      + ". Every engine returns the same matrix (deepTools differs only at bin edges, r ≥ 0.99).\n")
        md.append(figure(dotplot([(r["engine"], r["rate"], is_gb(r["engine"]), f"{r['n_loci']:,} regions")
                                  for r in eng], kind="green", unit="rate", label="Regions per second, one core"),
                         "<strong>Regions per second</strong> into a 200-bin matrix from one bigWig, one core "
                         "(higher is better; log scale)."))
    sc = rows_where(d, part="scaling")
    if sc:
        names = []
        for r in sc:
            if r["engine"] not in names:
                names.append(r["engine"])
        nice = {"genomeblocks · pybigtools (exact)": "genomeblocks (default)",
                "genomeblocks · pybigtools (zoom, exact=False)": "genomeblocks, exact=False",
                "genomeblocks · pure-Python reader": "genomeblocks, pure-Python",
                "pyBigWig values() + numpy binning": "pyBigWig values()"}
        series = [(nice.get(n, n), [(r["n_loci"], r["rate"]) for r in sc if r["engine"] == n],
                   n == "genomeblocks · pybigtools (exact)") for n in names]
        md.append(figure(linechart(series, kind="green", xlab="regions", label="Throughput vs number of regions"),
                         "<strong>Throughput vs number of regions</strong> (regions/s, log–log). The native path speeds "
                         "up with batch size as fixed costs amortise; per-base readers stay flat."))
    par = rows_where(d, part="parallel")
    if par:
        pr = [(f"{r['engine']} · {r['workers']} worker{'s' if r['workers'] > 1 else ''}", r["seconds"],
               is_gb(r["engine"]), f"{r['n_tracks']} tracks × {r['n_loci']:,} regions") for r in par]
        p1 = next((r for r in par if is_gb(r["engine"]) and r["workers"] == 1), None)
        pm = max((r for r in par if is_gb(r["engine"])), key=lambda r: r["workers"], default=None)
        if p1 and pm:
            md.append(f"Many tracks scale with processes: 16 bigWigs × 5,000 regions take {ftime(p1['seconds'])} "
                      f"sequentially and {ftime(pm['seconds'])} with `workers={pm['workers']}`. Threads do not help "
                      f"(pybigtools serialises Python threads), which is why `signal()` uses processes.\n")
        md.append(figure(dotplot(pr, kind="green", label="16 tracks, parallel"),
                         "<strong>16 bigWigs × 5,000 regions</strong>: worker processes vs threads vs deepTools "
                         "<code>-p</code> (wall time; log scale)."))


def sec_atlas(md, tl):
    d = load("atlas")
    if not d:
        return
    md.append('## Atlas\n{: .sec-navy #atlas }\n')
    q = rows_where(d, part="query")
    if q:
        nt = max(r["n_tracks"] for r in q)
        sel = [r for r in q if r["n_tracks"] == nt]
        s = next((r for r in sel if r["engine"].startswith("Atlas.search")), None)
        loop = next((r for r in sel if "cgranges loop" in r["engine"]), None)
        if s:
            tl.append(("Enrichment vs " + f"{nt} peak files", ftime(s["seconds"]),
                       "per 20k-peak query including Fisher tests" +
                       (f": {fx(loop['seconds'] / s['seconds'])} faster than looping over the tracks" if loop else ""),
                       "navy"))
            md.append(f"Testing a 20,000-peak query against {nt} peak files is one sparse row-sum plus a "
                      f"vectorised Fisher test: **{ftime(s['seconds'])}** per query, including building the result "
                      f"table" + (f", {fx(loop['seconds'] / s['seconds'])} faster than a per-track loop over "
                                  f"prebuilt cgranges indexes" if loop else "") + ".\n")
        md.append(figure(dotplot([(short(r["engine"]), r["seconds"], is_gb(r["engine"]), f"{nt} tracks") for r in sel],
                                 kind="navy", label=f"Query against {nt} tracks"),
                         f"<strong>One query against {nt} peak files</strong> (wall time; log scale). GIGGLE was not "
                         f"re-run for 1.1 (it has to be built from source)."))
    b = rows_where(d, part="build")
    if b:
        md.append('<details class="gb-table"><summary>Index build, load, accuracy</summary><table><thead><tr>'
                  '<th>step</th><th>bin size</th><th>workers</th><th>time</th><th>index size</th></tr></thead><tbody>'
                  + "".join(f"<tr><td>build</td><td>{r['bin_size']:,} bp</td><td>{r['workers']}</td>"
                            f"<td>{ftime(r['seconds'])}</td><td>{r.get('index_bytes', 0) / 1e6:.0f} MB</td></tr>"
                            for r in b)
                  + "".join(f"<tr><td>load ({escape(str(r.get('kind', '')))})</td><td></td><td></td>"
                            f"<td>{ftime(r['seconds'])}</td><td></td></tr>" for r in rows_where(d, part="load"))
                  + "".join(f"<tr><td>accuracy vs exact overlaps</td><td>{r['bin_size']:,} bp</td><td></td>"
                            f"<td>Spearman {r['spearman_counts']:.3f}</td><td>top-25 shared {r['top25_overlap']}</td></tr>"
                            for r in rows_where(d, part="accuracy"))
                  + "".join(f"<tr><td>bootstrap ({r['n_iter']} shuffles)</td><td></td><td></td>"
                            f"<td>{ftime(r['per_iter'])} per shuffle</td><td></td></tr>"
                            for r in rows_where(d, part="bootstrap"))
                  + "</tbody></table></details>\n")


def agree_sentence(eng):
    """Which engines return exactly the reference hit count, and by how much the others differ."""
    ref = next((r["ref_hits"] for r in eng if r.get("ref_hits")), None)
    if ref is None:
        return ""
    def short(e):
        for k, v in (("FIMO", "FIMO"), ("MOODS", "MOODS"), ("Biopython", "Biopython"), ("numpy", "numpy"),
                     ("re-striped", "lightmotif per window")):
            if k in e:
                return v
        return e
    same = [short(r["engine"]) for r in eng if r.get("total_hits") == ref and not is_gb(r["engine"])]
    diff = [(short(r["engine"]), r["total_hits"]) for r in eng
            if r.get("total_hits") is not None and r["total_hits"] != ref]
    out = f"The {ref:,} hits are identical to those of {', '.join(same[:-1])} and {same[-1]}"
    if diff:
        out += "; " + "; ".join(f"{n} reports {h:,}" for n, h in diff)
    return out + "."


def moods_vs_gb(lib):
    """Sentence comparing genomeblocks and MOODS at each window count."""
    parts = []
    for n in sorted({r["n_seqs"] for r in lib}):
        g = next((r["seconds"] for r in lib if r["n_seqs"] == n and is_gb(r["engine"])), None)
        m = next((r["seconds"] for r in lib if r["n_seqs"] == n and "MOODS" in r["engine"]), None)
        if g and m:
            faster = "genomeblocks" if g < m else "MOODS"
            parts.append(f"{faster} is {fx(max(g, m) / min(g, m))} faster at {n:,} windows")
    return ("Against MOODS: " + ", ".join(parts) + ".") if parts else ""


def sec_motifs(md, tl):
    d = load("motifs")
    if not d:
        return
    md.append('## Motifs\n{: .sec-navy #motifs }\n')
    eng = rows_where(d, part="engines")
    if eng:
        gb = next((r for r in eng if is_gb(r["engine"])), None)
        fimo = next((r for r in eng if "FIMO" in r["engine"]), None)
        bio = next((r for r in eng if "Biopython" in r["engine"]), None)
        if gb and fimo:
            tl.append(("Motif scanning", fx(fimo["seconds"] / gb["seconds"]),
                       "faster than MEME FIMO" + (f", {fx(bio['seconds'] / gb['seconds'])} faster than Biopython"
                                                  if bio else ""), "navy"))
            md.append(f"Counting hits of {gb['n_motifs']} JASPAR motifs in {gb['n_seqs']:,} windows of 500 bp: "
                      f"**{ftime(gb['seconds'])}** with genomeblocks (one block scan per motif), "
                      f"{fx(fimo['seconds'] / gb['seconds'])} faster than FIMO. "
                      + agree_sentence(eng) + "\n")
        md.append(figure(dotplot([(short(r["engine"]), r["seconds"], is_gb(r["engine"]), f"{r['n_seqs']:,} windows × "
                                  f"{r['n_motifs']} motifs") for r in eng], kind="navy", label="Motif scanning engines"),
                         "<strong>Motif scanning engines</strong> on identical PSSMs and windows (wall time; log scale)."))
    lib = rows_where(d, part="library")
    if lib:
        md.append(figure(dotplot([(f"{short(r['engine'])} · {fnum(r['n_seqs'])} windows", r["seconds"],
                                   is_gb(r["engine"]), f"{r['n_seqs']:,} windows × {r['n_motifs']} motifs")
                                  for r in lib], kind="navy", label="Whole library"),
                         f"<strong>The whole JASPAR library</strong> ({lib[0]['n_motifs']:,} motifs) on one core "
                         f"(wall time; log scale). " + moods_vs_gb(lib) + " FIMO was timed at 1,000 windows only."))
    w = rows_where(d, part="workers")
    if w:
        md.append('<details class="gb-table"><summary>Worker processes</summary><table><thead><tr><th>workers</th>'
                  '<th>time</th><th>speed-up</th></tr></thead><tbody>'
                  + "".join(f"<tr><td>{r['workers']}</td><td>{ftime(r['seconds'])}</td>"
                            f"<td>{fx(w[0]['seconds'] / r['seconds'])}</td></tr>" for r in w)
                  + f"</tbody></table></details>\n")


def sec_pairs(md, tl):
    d = load("pairs")
    if not d:
        return
    md.append('## Hi-C pairs\n{: .sec-purple #pairs }\n')
    rows = d["rows"]
    cp = next((r for r in rows if r["engine"].startswith("genomeblocks count_pairs (50")), None)
    cl = next((r for r in rows if r["engine"].startswith("cooler")), None)
    if cp and cl:
        tl.append(("Hi-C pair counting", f"{fnum(cp['rate'])}/s",
                   f"read pairs into 50 kb windows: {fx(cp['rate'] / cl['rate'])} the rate of cooler cload", "purple"))
        md.append(f"Streaming {fnum(cp['n_pairs'])} read pairs into 50 kb windows × partner chromosome runs at "
                  f"**{fnum(cp['rate'])} pairs/s**, close to the speed of just parsing the file, and "
                  f"{fx(cp['rate'] / cl['rate'])} the rate of `cooler cload`.\n")
    md.append(figure(dotplot([(short(r["engine"]), r["rate"], is_gb(r["engine"]), f"{r['n_pairs']:,} pairs") for r in rows],
                             kind="purple", unit="rate", label="Pairs per second"),
                     "<strong>Read pairs per second</strong> (higher is better; log scale). The naive loop was timed on "
                     "the first 200k pairs."))


def sec_genes(md, tl):
    d = load("genes")
    if not d:
        return
    md.append('## Genes\n{: .sec-navy #genes }\n')
    m = d.get("meta", {})
    md.append(f"A GENCODE-shaped GTF with {m.get('genes', 0):,} genes, {m.get('transcripts', 0):,} transcripts and "
              f"{m.get('gtf_lines', 0):,} lines.\n")
    rows = [(f"{r['step']}" + (f" · {r['n']:,}" if "annotations" in r["step"] or "nearest" in r["step"] else ""),
             r["seconds"], True, "") for r in d["rows"]]
    md.append(figure(dotplot(rows, kind="navy", label="Gene-model steps"),
                     "<strong>Gene-model steps</strong> (wall time; log scale). Labelling uses the cached annotation "
                     "index; the counts are CREs labelled."))


def sec_architecture(md, tl):
    d = load("architecture")
    p = load("prototype")
    if not d and not p:
        return
    md.append('## Architecture\n{: .sec-purple #architecture }\n')
    if d:
        rows = [(short(r["step"]), r["seconds"], not r["step"].startswith("annotate stage 2: per")
                 and "graph-tool" not in r["step"], "") for r in d["rows"]]
        mk = next((r for r in d["rows"] if r["step"].startswith("make")), None)
        if mk:
            md.append(f"The classic pipeline on 100k CREs and 50k loops ({mk.get('n_edges', 0):,} edges), "
                      f"step by step:\n")
        md.append(figure(dotplot(rows, kind="purple", label="Architecture steps"),
                         "<strong>Architecture pipeline steps</strong>, classic module (wall time; log scale). The grey "
                         "rows are the per-vertex loop that annotate's vectorised stage replaced and graph-tool's own "
                         "reduction for comparison."))


def sec_columnar(md, tl):
    p = load("prototype")
    if not p:
        return
    md.append('## Columnar vs classic\n{: .sec-purple #columnar }\n')
    steps = rows_where(p, part="steps")
    names = []
    for r in steps:
        if r["step"] not in names:
            names.append(r["step"])
    pairs = []
    for n in names:
        a = next((r["seconds"] for r in steps if r["step"] == n and r["impl"] == "main"), None)
        b = next((r["seconds"] for r in steps if r["step"] == n and r["impl"] == "columnar"), None)
        if a and b:
            pairs.append((n, a, b))
    e2e = rows_where(p, part="e2e")
    tot = {r["impl"]: r["seconds"] for r in e2e if r["step"] == "total"}
    mem = {r["impl"]: r["peak_mb"] for r in rows_where(p, part="e2e_mem")}
    if tot.get("main") and tot.get("columnar"):
        tl.append(("Whole pipeline, fresh process", f"{ftime(tot['main'])} → {ftime(tot['columnar'])}",
                   "classic vs columnar: load CREs and genes, build the graph, add Hi-C, normalise, annotate, hubs, "
                   "save", "purple"))
        md.append(f"The same pipeline (load CREs and genes, build, add Hi-C, normalise, annotate, strength, hubs, "
                  f"save) in a fresh process: **{ftime(tot['main'])}** with the classic modules, "
                  f"**{ftime(tot['columnar'])}** with `genomeblocks.columnar`"
                  + (f", peak memory {mem['main']:,.0f} MB → {mem['columnar']:,.0f} MB" if mem.get("main") else "")
                  + ". Both give the same edges, weights, O/E, labels and hubs.\n")
        ld = {r["impl"]: r for r in rows_where(p, part="e2e_load")}
        meta = pick(p, part="steps_meta") or {}
        if ld.get("main") and ld.get("columnar"):
            md.append(f"Reloading the saved graph in a new process takes {ftime(ld['main']['seconds'])} from the classic "
                      f"pickle and {ftime(ld['columnar']['seconds'])} from the columnar parquet tables"
                      + (f" ({meta['file_mb_main']:.0f} MB vs {meta['file_mb_columnar']:.1f} MB on disk)"
                         if meta.get("file_mb_main") else "") + ".\n")
    if pairs:
        md.append(figure(dumbbell(pairs, kind="purple", label="Classic vs columnar per step"),
                         "<strong>Per step, classic (grey) vs columnar</strong> on 100k CREs and 50k loops "
                         "(in-process medians; log scale; right column = speed-up)."))
    views = rows_where(p, part="views")
    if views:
        ops = []
        for r in views:
            if r["op"] not in ops:
                ops.append(r["op"])
        trs = []
        for o in ops:
            a = next((r["seconds"] for r in views if r["op"] == o and r["impl"] == "main"), None)
            b = next((r["seconds"] for r in views if r["op"] == o and r["impl"] == "columnar"), None)
            trs.append([escape(o), ftime(a) if a else "–", ftime(b) if b else "–", fx(a / b) if a and b else "–"])
        md.append('<details class="gb-table" open><summary>Interactive operations</summary><table><thead><tr><th>operation</th>'
                  '<th>classic</th><th>columnar</th><th>speed-up</th></tr></thead><tbody>'
                  + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in trs)
                  + "</tbody></table></details>\n")
    sc = rows_where(p, part="scale")
    if sc:
        series = []
        for impl, nm in (("main", "classic"), ("columnar", "columnar")):
            pts = [(r["loops"], r["seconds"]) for r in sc if r["impl"] == impl]
            if pts:
                series.append((nm, pts, impl == "columnar"))
        md.append(figure(linechart(series, kind="purple", xlab="loops", label="Pipeline time vs loop count",
                                   yfmt=ftime, ytick=ftick),
                         "<strong>Pipeline time vs loop count</strong> (make → add_mcool → normalize → annotate → "
                         "strength; log–log)."))


def sec_shortrange(md, tl):
    d = load("shortrange")
    if not d:
        return
    md.append('## HiChIP short-range track\n{: .sec-purple #shortrange }\n')
    tot = {r["impl"]: r for r in rows_where(d, part="total")}
    if len(tot) == 2:
        sh = next(v for k, v in tot.items() if k != "genomeblocks")
        gb = tot.get("genomeblocks")
        if gb:
            md.append(f"From allValidPairs to a coverage bigWig: **{ftime(gb['seconds'])}** and "
                      f"{gb['peak_mb'] / 1e3:.1f} GB peak memory with `columnar.hichip`, against "
                      f"{ftime(sh['seconds'])} and {sh['peak_mb'] / 1e3:.1f} GB for the `awk | sort | bedtools` recipe, "
                      f"with the same ends and a byte-identical bedGraph.\n")
    st = rows_where(d, part="steps")
    if st:
        md.append('<details class="gb-table"><summary>Steps</summary><table><thead><tr><th>implementation</th>'
                  '<th>step</th><th>time</th><th>peak memory</th></tr></thead><tbody>'
                  + "".join(f"<tr><td>{escape(r['impl'])}</td><td>{escape(r['step'])}</td><td>{ftime(r['seconds'])}</td>"
                            + (f"<td>{r['peak_mb']:,.0f} MB</td></tr>" if r.get('peak_mb') else '<td>–</td></tr>')
                            for r in st)
                  + "</tbody></table></details>\n")


def sec_import(md, tl):
    d = load("import")
    if not d:
        return
    md.append('## Import cost\n{: .sec-navy #import }\n')
    base = next((r["seconds"] for r in d["rows"] if r["stmt"] == "python -c pass"), 0)
    rows = [(r["stmt"], r["seconds"], True, "") for r in d["rows"]
            if r["stmt"] != "python -c pass" and not r["stmt"].startswith("breakdown")]
    loci = next((r["seconds"] for r in d["rows"] if r["stmt"].startswith("+ Loci")), None)
    if loci:
        md.append(f"`import genomeblocks` is lazy, and `from genomeblocks import Loci` costs **{ftime(loci)}** on top "
                  f"of starting Python ({ftime(base)}): matplotlib, pandas, scipy and graph-tool load only when a "
                  f"function needs them.\n")
    md.append(figure(dotplot(rows, kind="navy", label="Import cost"),
                     "<strong>Import cost</strong> above <code>python -c pass</code>, fresh process each "
                     "(cumulative statements; log scale)."))


def build():
    md, tl = [], []
    secs = [sec_loci, sec_signal, sec_atlas, sec_motifs, sec_pairs, sec_genes, sec_architecture, sec_columnar,
            sec_shortrange, sec_import]
    for s in secs:
        s(md, tl)
    env = next((load(n)["env"] for n in ("loci", "signal", "motifs") if load(n)), {})
    v = env.get("versions", {})
    vtxt = ", ".join(f"{k} {val}" for k, val in v.items())
    head = f"""---
title: Benchmarks
layout: default
nav_order: 9
permalink: /benchmarks/
---

# Benchmarks
{{: .no_toc }}

How fast each building block is against the tools people would otherwise use,
measured on genomeblocks 1.1. Every comparison first checks that the engines
return the same answer; timings are medians of 3–5 runs after a warm-up.
{{: .fs-5 .fw-300 }}

<p class="gb-env">{escape(env.get('cpu', ''))} · {env.get('cores', '?')} cores · {env.get('mem_gb', 0):.0f} GB RAM ·
Python {escape(env.get('python', ''))} · synthetic hg38-shaped data (seeded) · benchmarks ran one at a time and use at
most 8 worker processes</p>

{tiles(tl)}

1. TOC
{{:toc}}

"""
    tail = f"""
## Reproduce
{{: .sec-navy #reproduce }}

```bash
cd benchmarks
python make_data.py            # synthetic hg38-shaped data, ~1.5 GB
PY=python ./run_all.sh         # every bench, one at a time
python build_site.py           # this page: docs/benchmarks/index.md
```

External baselines are found on `PATH` (bedtools, deepTools `computeMatrix`,
cooler) or via `FIMO=`. The data are synthetic and seeded, laid out on real hg38
chromosome sizes; motifs are the JASPAR CORE vertebrate collection.

<p class="gb-env">Versions: {escape(vtxt)}</p>
"""
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(head + "\n".join(md) + tail)
    print(f"[site] {OUT}")


if __name__ == "__main__":
    build()
