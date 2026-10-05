"""Shareable genome browser: one self-contained HTML file (igv.js + embedded data).

The reader only needs a web browser: no install, no server, no data files.
Every track is embedded in the page as a gzipped data URI, and igv.js
(the JavaScript IGV) draws it. The page opens on a list of regions to look
at, each a one-click button; users can still type any locus or gene name.

    >>> from genomeblocks.columnar.igv import igv_html
    >>> igv_html("share.html", regions=["chr8:127.7-128.0 Mb", ...],
    ...          loci={"CREs": cre}, genes=genes, signal={"ATAC": "atac.bw"},
    ...          architecture=A, title="MYC enhancer hubs")

Size is the only limit: signal and loops are embedded for the listed
regions only (± ``flank``), CREs and genes for the whole genome. The
function prints the size of every track.
"""
from __future__ import annotations

import base64
import gzip
import html
import json
from typing import Dict, Optional, Sequence

import numpy as np

IGV_VERSION = "3.8.9"
IGV_CDN = f"https://cdn.jsdelivr.net/npm/igv@{IGV_VERSION}/dist/igv.min.js"
_PALETTE = ["#1f6f8b", "#d1495b", "#2a9d8f", "#e9a03b", "#6a4c93", "#577590"]


def _uri(text: str) -> str:
    return "data:application/gzip;base64," + base64.b64encode(gzip.compress(text.encode(), 9)).decode()


def _region(r):
    """'chr1:1,000-2,000', 'chr1:1.2-1.5 Mb', (chrom, start, end) or a Locus -> (chrom, start, end)."""
    if hasattr(r, "chrom"):
        return r.chrom, int(r.start), int(r.end)
    if isinstance(r, (tuple, list)):
        return str(r[0]), int(r[1]), int(r[2])
    chrom, rest = r.split(":")
    rest = rest.replace(",", "").strip()
    scale = 1
    if rest.lower().endswith("mb"):
        scale, rest = 1_000_000, rest[:-2]
    elif rest.lower().endswith("kb"):
        scale, rest = 1_000, rest[:-2]
    a, b = rest.split("-")
    return chrom.strip(), int(float(a) * scale), int(float(b) * scale)


def _windows(regions, flank):
    out = []
    for r in regions:
        c, a, b = _region(r)
        out.append((c, max(0, a - flank), b + flank))
    return out


def _bed_loci(L, names=None) -> str:
    ch = L.chroms
    nm = L.uid if names is None else names
    strand = np.array([".", "+", "-"], dtype=object)[L.strands]
    return "\n".join(f"{c}\t{s}\t{e}\t{n}\t0\t{d}" for c, s, e, n, d in
                     zip(ch, L.starts.tolist(), L.ends.tolist(), nm, strand))


def _bed12_genes(genes) -> str:
    """One BED12 line per gene: its longest transcript, exons as blocks."""
    T, F, G = genes.transcripts, genes.features, genes.genes
    tlen = T.ends - T.starts
    order = np.lexsort((-tlen, T.cols["gene"]))
    g_sorted = T.cols["gene"][order]
    first = order[np.r_[True, g_sorted[1:] != g_sorted[:-1]]]          # longest transcript per gene
    ex = F.cols["kind"] == 0
    f_tx, f_s, f_e = F.cols["transcript"][ex], F.starts[ex] - 1, F.ends[ex]   # GTF is 1-based
    o = np.lexsort((f_s, f_tx))
    f_tx, f_s, f_e = f_tx[o], f_s[o], f_e[o]
    lo = np.searchsorted(f_tx, first)
    hi = np.searchsorted(f_tx, first, side="right")
    names = G.cols["gene_name"]
    strand = np.array([".", "+", "-"], dtype=object)[T.strands]
    lines = []
    for t, a, b in zip(first.tolist(), lo.tolist(), hi.tolist()):
        s, e = int(T.starts[t]) - 1, int(T.ends[t])
        if b > a:
            bs, be = f_s[a:b], f_e[a:b]
        else:
            bs, be = np.array([s]), np.array([e])
        sizes = ",".join(str(int(x)) for x in be - bs)
        starts = ",".join(str(int(x)) for x in bs - s)
        lines.append(f"{T.chroms[t]}\t{s}\t{e}\t{names[T.cols['gene'][t]]}\t0\t{strand[t]}\t{s}\t{e}\t0\t"
                     f"{len(bs)}\t{sizes}\t{starts}")
    return "\n".join(lines)


def _bedgraph(bw_path, windows, bin_size) -> str:
    import pybigtools
    out = []
    with pybigtools.open(bw_path) as bw:
        sizes = bw.chroms()
        for c, a, b in windows:
            if c not in sizes:
                continue
            b = min(b, sizes[c])
            n = max(1, (b - a) // bin_size)
            v = bw.values(c, a, a + n * bin_size, bins=n, summary="mean", missing=0.0)
            edges = a + np.arange(n + 1) * bin_size
            out += [f"{c}\t{s}\t{e}\t{x:.3g}" for s, e, x in zip(edges[:-1].tolist(), edges[1:].tolist(), v.tolist())
                    if x > 0]
    return "\n".join(out)


def _interact(A, windows, score) -> str:
    """Edges with at least one end in a window, as BEDPE (score = ep[score])."""
    L = A.loci
    keep = np.zeros(A.n_links, bool)
    for c, a, b in windows:
        code = L.genome.code.get(c)
        if code is None:
            continue
        inw = (L.codes == code) & (L.ends > a) & (L.starts < b)
        keep |= inw[A.src] | inw[A.tgt]
    s, t = A.src[keep], A.tgt[keep]
    w = A.ep[score][keep] if score in A.ep else np.ones(keep.sum())
    ch = L.chroms
    return "\n".join(
        f"{ch[i]}\t{L.starts[i]}\t{L.ends[i]}\t{ch[j]}\t{L.starts[j]}\t{L.ends[j]}\t.\t{x:.3g}"
        for i, j, x in zip(s.tolist(), t.tolist(), w.tolist()))


def igv_html(path: str, *, regions: Sequence, loci: Optional[Dict] = None, genes=None,
             signal: Optional[Dict[str, str]] = None, architecture=None, score: str = "n",
             chrom_sizes: Optional[Dict[str, int]] = None, genome_id: Optional[str] = None,
             flank: int = 250_000, bin_size: int = 50, title: str = "Genome browser",
             notes: Optional[Dict[str, str]] = None, igv_js: str = "cdn",
             standalone: bool = True, subtitle: Optional[str] = None) -> dict:
    """Write a single HTML file that opens an IGV browser on ``regions``.

    Args:
        regions: loci to offer as one-click buttons (strings, tuples or Locus).
            Two regions joined by a space ('chr1:… chr7:…') open side by side,
            which is how a trans loop is shown.
        loci: {track name: columnar Loci} — embedded genome-wide.
        genes: columnar Genes — one gene model per gene, genome-wide.
        signal: {track name: bigWig path} — embedded around ``regions`` only.
        architecture: columnar Architecture — loops touching ``regions`` as arcs.
        score: edge column shown as arc height (default O/E ``n``).
        chrom_sizes: {chrom: length}; default from the Loci's genome / bigWig.
        genome_id: e.g. 'hg38' to use igv.js's hosted genome (sequence,
            ideogram) instead of embedded chromosome sizes. Needs internet.
        notes: optional {region: one-line note} shown next to each button.
        igv_js: 'cdn' (load igv.js from jsDelivr; the page stays small), a path
            to a local igv.min.js to embed it (~1.5 MB; the file then works
            fully offline), or any other URL.
        standalone: write a full HTML document (False: page content only, for
            hosts that add their own <html>/<head>, e.g. a Claude artifact).
    Returns:
        sizes of the embedded tracks in bytes.
    """
    loci = loci or {}
    signal = signal or {}
    flat = [x for r in regions for x in str(r).split()] if regions else []
    windows = _windows(flat, flank)
    tracks, sizes = [], {}

    def add(name, text, cfg):
        uri = _uri(text)
        sizes[name] = len(uri)
        tracks.append({"name": name, "url": uri, **cfg})

    k = 0
    for name, sig in signal.items():
        add(name, _bedgraph(sig, windows, bin_size),
            {"type": "wig", "format": "bedgraph", "color": _PALETTE[k % len(_PALETTE)], "height": 60,
             "autoscale": True})
        k += 1
    for name, L in loci.items():
        add(name, _bed_loci(L), {"type": "annotation", "format": "bed", "displayMode": "SQUISHED",
                                 "color": _PALETTE[k % len(_PALETTE)], "height": 30})
        k += 1
    if architecture is not None:
        add(f"loops ({score})", _interact(architecture, windows, score),
            {"type": "interact", "format": "bedpe", "color": "#d1495b", "arcType": "nested",
             "height": 120, "showBlocks": True})
    if genes is not None:
        add("genes", _bed12_genes(genes), {"type": "annotation", "format": "bed", "displayMode": "EXPANDED",
                                           "searchable": True, "color": "#333333", "height": 120})

    if genome_id:
        genome = genome_id
    else:
        if chrom_sizes is None:
            chrom_sizes = {}
            for L in list(loci.values()) + ([architecture.loci] if architecture is not None else []):
                for c in np.unique(L.chroms):
                    chrom_sizes[c] = max(chrom_sizes.get(c, 0), int(L.ends[L.chroms == c].max()) + flank)
        cs = "\n".join(f"{c}\t{n}" for c, n in chrom_sizes.items())
        # chromosome sizes only (no sequence). igv.js 3.8.9 mishandles a data URI
        # here, so the page hands them over as an in-memory File instead.
        genome = {"id": "custom", "name": title, "format": "chromsizes"}

    first = str(regions[0]) if regions else None
    cfg = {"locus": " ".join(f"{c}:{a + 1}-{b}" for c, a, b in
                             (_region(x) for x in first.split())) if first else None,
           "tracks": tracks, "showSampleNames": False, "showChromosomeWidget": True}
    if isinstance(genome, str):
        cfg["genome"] = genome
    else:
        cfg["reference"] = genome

    def loc(r):
        return " ".join(f"{c}:{a + 1}-{b}" for c, a, b in (_region(x) for x in str(r).split()))
    notes = notes or {}
    buttons = "".join(
        f'<button type="button" class="go" data-locus="{html.escape(loc(r))}">'
        f'<span class="lbl">{html.escape(str(r))}</span>'
        f'{"<span class=note>" + html.escape(notes[str(r)]) + "</span>" if str(r) in notes else ""}</button>'
        for r in regions)
    import os
    if igv_js != "cdn" and os.path.isfile(igv_js):
        with open(igv_js) as f:
            script = f"<script>{f.read()}</script>"
    else:
        script = f'<script src="{IGV_CDN if igv_js == "cdn" else igv_js}"></script>'

    page = _PAGE.format(title=html.escape(title), buttons=buttons, script=script,
                        config=json.dumps(cfg), n_regions=len(regions),
                        sizes=json.dumps(None if genome_id else cs),
                        sub=html.escape(subtitle) if subtitle else "",
                        tracks=", ".join(t["name"] for t in tracks))
    if not standalone:
        import re
        page = re.sub(r"<!doctype html>\s*<html[^>]*>\s*<head>|</head>\s*<body>|</body>\s*</html>\s*$", "", page)
        page = re.sub(r'<meta [^>]*>\s*', "", page)
    with open(path, "w") as f:
        f.write(page)
    sizes["total"] = len(page)
    return sizes


_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{title}</title>
<style>
:root {{ --bg:#f6f4e8; --surface:#fbfaf4; --ink:#1e1e1e; --ink2:#52514b; --rule:#e1dcc5; --brand:#003552; --hl:#e8eff4; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#1b1b1a; --surface:#232220; --ink:#f2f0e6; --ink2:#c7c3b1; --rule:#3a3935; --brand:#8cc3e8; --hl:#1d2a33; color-scheme:dark; }} }}
:root[data-theme="dark"] {{ --bg:#1b1b1a; --surface:#232220; --ink:#f2f0e6; --ink2:#c7c3b1; --rule:#3a3935; --brand:#8cc3e8; --hl:#1d2a33; color-scheme:dark; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
.wrap {{ max-width:1280px; margin:0 auto; padding-inline:16px; padding-block:18px 40px; display:grid; gap:14px; }}
h1 {{ font-size:22px; margin:0; color:var(--brand); text-wrap:balance; }}
.sub {{ color:var(--ink2); font-size:14px; margin:0; }}
.layout {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:14px; align-items:start; }}
@media (max-width:820px) {{ .layout {{ grid-template-columns:1fr; }} }}
.regions {{ display:grid; gap:6px; align-content:start; }}
.regions h2 {{ font-size:12px; letter-spacing:.06em; text-transform:uppercase; color:var(--ink2); margin:0 0 2px; }}
button.go {{ text-align:left; display:grid; gap:2px; background:var(--surface); color:var(--ink); border:1px solid var(--rule);
  border-radius:8px; padding:8px 10px; cursor:pointer; font:inherit; font-size:13.5px; }}
button.go .note {{ font-size:12px; color:var(--ink2); }}
button.go:hover, button.go:focus-visible {{ border-color:var(--brand); outline:none; }}
button.go[aria-pressed="true"] {{ background:var(--hl); border-color:var(--brand); }}
#igv {{ background:#fff; border:1px solid var(--rule); border-radius:8px; padding:6px; min-height:420px; min-width:0; overflow-x:auto; }}
.help {{ font-size:13px; color:var(--ink2); }}
</style></head>
<body><div class="wrap">
<h1>{title}</h1>
<p class="sub">{sub} Click a region on the left, or type a locus or gene name in the search box. Drag to pan, double-click to zoom. Tracks: {tracks}.</p>
<div class="layout">
  <nav class="regions" aria-label="Regions"><h2>Regions to look at ({n_regions})</h2>{buttons}
    <p class="help">All data is inside this one file: nothing to install, nothing to upload.</p></nav>
  <div id="igv"></div>
</div></div>
{script}
<script>
const CONFIG = {config};
const SIZES = {sizes};
if (SIZES) CONFIG.reference.fastaURL = new File([SIZES], "genome.chrom.sizes");
igv.createBrowser(document.getElementById("igv"), CONFIG).then(b => {{
  window.browser = b;
  const btns = document.querySelectorAll("button.go");
  btns.forEach((x, i) => {{
    if (i === 0) x.setAttribute("aria-pressed", "true");
    x.addEventListener("click", () => {{
      btns.forEach(y => y.setAttribute("aria-pressed", "false"));
      x.setAttribute("aria-pressed", "true");
      b.search(x.dataset.locus);
    }});
  }});
  document.body.dataset.ready = "1";
}}).catch(e => {{ document.getElementById("igv").textContent = "Could not start the viewer: " + e; document.body.dataset.ready = "error"; }});
</script>
</body></html>
"""
