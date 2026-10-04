"""genomeblocks view: our own interactive browser for genomeblocks tables.

One HTML file, no install, no server, works offline. Data travels as the
same columns genomeblocks keeps in memory (gzipped typed arrays); a small
canvas viewer (``_view.js`` + ``_view.css``) draws them. Because the viewer
knows the tables it can do what a generic browser cannot: search a gene and
its anchor CREs (TSS ± r) light up, the anchor's contact profile is drawn
per sample, partners are listed with their O/E, and a trans partner opens
side by side.

Build it like a figure, one track per call (tracks appear in call order)::

    v = View(A, genes=genes, samples={"LNCaP": "#46a8e4", "LuCaP35CR": "#ffa600"})
    v.anchor_profile("graph O/E (anchor)")
    v.cre_values("node strength", {s: A.vp[f"strength_{s}"] for s in CELLS})
    v.signal("ATAC", ATBW)                       # {sample: bigWig}
    v.intervals("SE", SE); v.intervals("prime", PRIME)   # {sample: Loci or CRE mask}
    v.points("copy number", CNR, segments=SV)    # {sample: DataFrame}
    v.cres(); v.loops({s: f"n_{s}" for s in CELLS}); v.genes()
    v.region("MYC", gene="MYC", locus="chr8:124.8-129.2 Mb"); v.mark("E-MYC", "chr8", 125_234_800)
    v.save("myc.html")        # or just display `v` in a notebook

For BAM/VCF review or huge remote files, export to IGV instead
(:func:`genomeblocks.columnar.igv.igv_html`).
"""
from __future__ import annotations

import base64
import gzip
import html as _html
import json
import re
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .igv import _region

_HERE = Path(__file__).resolve().parent
_AUTO = ["#46a8e4", "#ffa600", "#2a9d8f", "#d1495b", "#6a4c93", "#577590"]


def _pack(a, dtype) -> dict:
    a = np.ascontiguousarray(a, dtype=dtype)
    return {"t": np.dtype(dtype).name, "n": int(a.size),
            "b": base64.b64encode(gzip.compress(a.tobytes(), 6)).decode()}


def _delta(x) -> np.ndarray:
    x = np.asarray(x, np.int64)
    return np.diff(x, prepend=0) if len(x) else x


def _natural(name):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


class View:
    """An interactive genome view of genomeblocks tables, saved as one HTML file."""

    def __init__(self, architecture=None, *, cre=None, genes=None, samples=None, title="genomeblocks view",
                 subtitle="", anchor_r: int = 5000, anchor_mode: str = "center", gene_half: int = 500_000,
                 chrom_sizes=None):
        self.A = architecture
        self.L = architecture.loci if architecture is not None else cre
        if self.L is None:
            raise ValueError("pass an Architecture or cre=Loci")
        self.G = genes
        if samples is None:
            samples = {}
        if not isinstance(samples, dict):
            samples = {s: _AUTO[i % len(_AUTO)] for i, s in enumerate(samples)}
        self.samples = list(samples.items())
        self.title, self.subtitle = title, subtitle
        self.anchor_r, self.anchor_mode, self.gene_half = anchor_r, anchor_mode, gene_half
        if isinstance(chrom_sizes, str):
            chrom_sizes = {l.split()[0]: int(l.split()[1]) for l in open(chrom_sizes) if l.strip()}
        self.chrom_sizes = dict(chrom_sizes or {})
        self._tracks, self._regions, self._marks, self._hubs = [], [], [], None
        self._scores = None

    # ── helpers ───────────────────────────────────────────────────────────
    def _sample(self, key):
        names = [s for s, _ in self.samples]
        if key in names:
            return names.index(key)
        if key is None:
            return None
        self.samples.append((key, _AUTO[len(self.samples) % len(_AUTO)]))
        return len(self.samples) - 1

    def _per_sample(self, x):
        return x.items() if isinstance(x, dict) else [(None, x)]

    # ── tracks (drawn in the order they are added) ────────────────────────
    def anchor_profile(self, name="contacts of the anchor (O/E)", smooth=1.5, height=None):
        """Per-sample sum of edge weight from the selected gene / CRE to every partner."""
        self._tracks.append({"type": "anchor", "name": name, "smooth": smooth, "height": height})
        return self

    def cre_values(self, name, values, smooth=1.0, height=None):
        """Per-CRE numbers (aligned to the CRE rows), drawn as a profile; {sample: array} overlays."""
        self._tracks.append({"type": "creval", "name": name, "smooth": smooth, "height": height,
                             "raw": [(self._sample(k), np.asarray(v, float)) for k, v in self._per_sample(values)]})
        return self

    def signal(self, name, bigwigs, *, bin_size=50, flank=250_000, genome_bin=20_000, height=None):
        """bigWig signal, fine bins around the listed regions and coarse bins genome-wide."""
        self._tracks.append({"type": "signal", "name": name, "height": height, "bin": bin_size,
                             "flank": flank, "gbin": genome_bin,
                             "raw": [(self._sample(k), v) for k, v in self._per_sample(bigwigs)]})
        return self

    def intervals(self, name, sets, color=None, height=None):
        """Interval rows: {sample or label: Loci / list of Locus / boolean mask over the CRE rows}."""
        raw = []
        for k, v in self._per_sample(sets):
            si = self._sample(k) if (k is None or k in dict(self.samples) or color is None) else None
            raw.append((si, None if si is not None or k is None else str(k), v))
        self._tracks.append({"type": "intervals", "name": name, "color": color, "height": height, "raw": raw})
        return self

    def points(self, name, tables, *, chrom="chromosome", start="start", end="end", value="log2",
               segments=None, seg_cols=("chrom", "start", "end", "label"), ylim=None, height=None):
        """Scatter of values along the genome (e.g. CNVkit .cnr: chromosome/start/end/log2),
        with optional segments drawn as labelled bars (e.g. SV calls)."""
        raw = []
        for k, df in self._per_sample(tables):
            if isinstance(df, str):
                import pandas as pd
                df = pd.read_csv(df, sep="\t")
            pos = ((df[start].to_numpy() + df[end].to_numpy()) // 2) if end in df else df[start].to_numpy()
            raw.append((self._sample(k), df[chrom].astype(str).to_numpy(), pos, df[value].to_numpy(float)))
        segs = []
        for k, df in self._per_sample(segments or {}):
            c, a, b, lab = seg_cols
            for row in df.itertuples(index=False):
                d = row._asdict()
                segs.append((self._sample(k), str(d[c]), int(d[a]), int(d[b]), str(d.get(lab, "")) if lab else ""))
        self._tracks.append({"type": "points", "name": name, "height": height, "ylim": ylim, "raw": raw, "segs": segs})
        return self

    def cres(self, name="CREs", height=None):
        """The CRE rows (promoters darker); the anchor and its partners are highlighted."""
        self._tracks.append({"type": "cre", "name": name, "height": height})
        return self

    def loops(self, score="n", name="loops", height=None):
        """Edges as arcs (trans edges as labelled stubs); score: an edge column or {sample: column}."""
        if self.A is None:
            raise ValueError("loops need an Architecture")
        self._scores = [(self._sample(k), np.asarray(self.A.ep[col], float)) for k, col in self._per_sample(score)]
        self._score_name = "O/E" if all(str(c).startswith("n") for _, c in self._per_sample(score)) else "score"
        self._tracks.append({"type": "loops", "name": name, "height": height})
        return self

    def genes(self, name="genes", height=None):
        self._tracks.append({"type": "genes", "name": name, "height": height})
        return self

    # ── navigation ────────────────────────────────────────────────────────
    def region(self, label, locus=None, *, gene=None, note="", anchor=None):
        """A one-click button. ``locus`` may hold two loci ('chrA:… chrB:…') for a split view;
        ``gene`` anchors that gene, ``anchor`` a CRE row."""
        if locus is None and gene is not None and self.G is not None:
            g = self.G[gene]
            tss = g.start if g.strand == "+" else g.end
            locus = f"{g.chrom}:{max(0, tss - self.gene_half)}-{tss + self.gene_half}"
        self._regions.append({"label": label, "note": note, "locus": locus or label, "gene": gene,
                              "anchorRow": None if anchor is None else int(anchor)})
        return self

    def mark(self, label, chrom, pos):
        self._marks.append({"label": label, "chrom": chrom, "pos": int(pos)})
        return self

    def hubs(self, values, n=12, label=None):
        """'Top hubs' list: the n CRE rows with the largest ``values``."""
        v = np.asarray(values, float)
        top = np.argsort(-v)[:n]
        self._hubs = [{"row": int(i), "label": self._uid(i),
                       "note": f"{label or 'value'} {v[i]:.1f}" + self._cre_note(i)} for i in top]
        return self

    def _uid(self, i):
        L = self.L
        return f"{L.chroms[i]}:{int(L.starts[i]):,}-{int(L.ends[i]):,}"

    def _cre_note(self, i):
        A = self.A
        if A is None:
            return ""
        out = ""
        if "annot" in A.vp:
            out += f" · {A.vp.annot[i]}"
        for k in A.vp:
            if k.startswith("gene") and A.vp[k][i]:
                out += f" · {A.vp[k][i]}"
                break
        return out

    # ── encoding ──────────────────────────────────────────────────────────
    def _payload(self) -> dict:
        from .genes import LABELS
        L, A, G = self.L, self.A, self.G
        g = L.genome
        names = set(L.genome.names[c] for c in np.unique(L.codes))
        if G is not None:
            names |= set(G.genes.chroms)
        for t in self._tracks:
            if t["type"] == "points":
                for _, ch, _, _ in t["raw"]:
                    names |= set(np.unique(ch))
        chroms = sorted(names, key=_natural)
        cidx = {c: i for i, c in enumerate(chroms)}
        remap = np.array([cidx.get(n, -1) for n in g.names], np.int64)
        sizes = []
        for c in chroms:
            n = self.chrom_sizes.get(c) or g.size(c)
            if not n:
                n = 0
                if c in g.code:
                    m = L.ends[L.codes == g.code[c]]
                    n = int(m.max()) if len(m) else 0
                    if G is not None:
                        gm = G.genes.ends[G.genes.codes == g.code[c]]
                        n = max(n, int(gm.max()) if len(gm) else 0)
                n += 250_000
            sizes.append(int(n))
        D = {"title": self.title, "chroms": chroms, "sizes": sizes, "labels": list(LABELS),
             "samples": [{"name": s, "color": c} for s, c in self.samples],
             "anchorR": self.anchor_r, "anchorMode": self.anchor_mode, "geneHalf": self.gene_half,
             "marks": self._marks}
        parts = {}

        def iv(codes, starts, ends):                       # sorted interval set, remapped chroms
            codes = remap[np.asarray(codes)] if not isinstance(codes, np.ndarray) or codes.dtype.kind != "O" else \
                np.array([cidx.get(c, -1) for c in codes])
            ok = codes >= 0
            codes, starts, ends = codes[ok], np.asarray(starts)[ok], np.asarray(ends)[ok]
            o = np.lexsort((starts, codes))
            return {"chrom": _pack(codes[o], np.uint8), "start": _pack(_delta(starts[o]), np.int32),
                    "len": _pack(np.maximum(ends[o] - starts[o], 1), np.uint32)}, o, ok

        # CRE rows (already genome-sorted; natural order = our chrom order)
        cre = {"chrom": _pack(remap[L.codes], np.uint8), "start": _pack(_delta(L.starts), np.int32),
               "len": _pack(L.ends - L.starts, np.uint32)}
        gene_names = []
        if A is not None:
            if "annot" in A.vp:
                lab = {v: i for i, v in enumerate(LABELS)}
                cre["label"] = _pack([lab.get(x, 0) for x in A.vp.annot], np.uint8)
            gk = next((k for k in A.vp if k.startswith("gene")), None)
            if gk:
                gene_names, inv = np.unique(np.asarray(A.vp[gk], dtype=object).astype(str), return_inverse=True)
                gene_names = gene_names.tolist()
                gidx = inv.astype(np.int32)
                if "" in gene_names:
                    gidx[gidx == gene_names.index("")] = -1
                cre["gene"] = _pack(gidx, np.int32)
        D["cre"], D["creGeneNames"] = cre, gene_names
        parts["CREs"] = sum(len(v["b"]) for v in cre.values())

        # edges with per-sample scores
        if A is not None:
            scores = self._scores or [(None, np.asarray(A.ep.get("n", A.ep.get("w", np.zeros(A.n_links))), float))]
            D["edges"] = {"src": _pack(_delta(A.src), np.int32),
                          "dt": _pack(A.tgt.astype(np.int64) - A.src, np.int32),
                          "scores": [{"sample": s, "v": _pack(np.nan_to_num(v), np.float32)} for s, v in scores],
                          "scoreName": getattr(self, "_score_name", "O/E")}
            parts["edges"] = len(D["edges"]["src"]["b"]) + len(D["edges"]["dt"]["b"]) + \
                sum(len(x["v"]["b"]) for x in D["edges"]["scores"])

        # genes: longest transcript per gene
        if G is not None:
            T, F, GG = G.transcripts, G.features, G.genes
            tlen = T.ends - T.starts
            order = np.lexsort((-tlen, T.cols["gene"]))
            gs = T.cols["gene"][order]
            tx = order[np.r_[True, gs[1:] != gs[:-1]]] if len(order) else order
            tx = tx[np.lexsort((T.starts[tx], remap[T.codes[tx]]))]
            tx = tx[remap[T.codes[tx]] >= 0]
            ex = F.cols["kind"] == 0
            f_tx, f_s, f_e = F.cols["transcript"][ex], F.starts[ex] - 1, F.ends[ex]
            o = np.lexsort((f_s, f_tx))
            f_tx, f_s, f_e = f_tx[o], f_s[o], f_e[o]
            lo, hi = np.searchsorted(f_tx, tx), np.searchsorted(f_tx, tx, side="right")
            ex_idx = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)]) if len(tx) else np.zeros(0, int)
            gstart = T.starts[tx] - 1
            D["genes"] = {"chrom": _pack(remap[T.codes[tx]], np.uint8), "start": _pack(_delta(gstart), np.int32),
                          "len": _pack(T.ends[tx] - gstart, np.uint32), "strand": _pack(T.strands[tx], np.uint8),
                          "nex": _pack(hi - lo, np.uint16),
                          "exStart": _pack(f_s[ex_idx] - np.repeat(gstart, hi - lo), np.int32),
                          "exLen": _pack(f_e[ex_idx] - f_s[ex_idx], np.uint32),
                          # anchor base, as Architecture.support: gene start on '+', gene end on '-'
                          "tss": _pack(np.where(GG.strands[T.cols["gene"][tx]] == 2, GG.ends[T.cols["gene"][tx]],
                                                GG.starts[T.cols["gene"][tx]]), np.int32),
                          "names": GG.cols["gene_name"][T.cols["gene"][tx]].astype(str).tolist()}
            parts["genes"] = sum(len(v["b"]) for v in D["genes"].values() if isinstance(v, dict)) + \
                len(json.dumps(D["genes"]["names"]))

        # regions → loci
        regions = []
        for r in self._regions:
            loci = [list(_region(x)) for x in str(r["locus"]).split()]
            regions.append({"label": r["label"], "note": r["note"], "loci": loci, "gene": r["gene"],
                            "anchorRow": r["anchorRow"]})
        D["regions"] = regions
        if self._hubs is not None:
            D["hubs"] = self._hubs
        elif A is not None and "strength" in A.vp:
            self.hubs(A.vp.strength, label="strength")
            D["hubs"] = self._hubs

        # tracks
        creCols, tracks = [], []
        windows = sorted({(c, max(0, a), b) for r in regions for c, a, b in r["loci"]})
        for t in self._tracks:
            o = {"type": t["type"], "name": t["name"]}
            if t.get("height"):
                o["height"] = t["height"]
            if t["type"] == "anchor":
                o["smooth"] = t["smooth"]
            elif t["type"] == "creval":
                o["smooth"], o["series"] = t["smooth"], []
                for s, v in t["raw"]:
                    creCols.append({"name": t["name"], "sample": s, "v": _pack(v, np.float32)})
                    o["series"].append({"sample": s, "col": len(creCols) - 1})
                parts[f"track: {t['name']}"] = sum(len(c["v"]["b"]) for c in creCols[-len(t["raw"]):])
            elif t["type"] == "signal":
                import pybigtools
                o["series"], size = [], 0
                for s, path in t["raw"]:
                    coarse, fine = [], []
                    with pybigtools.open(path) as bw:
                        csz = bw.chroms()
                        for c in chroms:
                            if c not in csz:
                                continue
                            n = max(1, csz[c] // t["gbin"])
                            v = bw.values(c, 0, n * t["gbin"], bins=n, summary="max", missing=0.0)
                            coarse.append({"chrom": cidx[c], "bin": t["gbin"], "v": _pack(np.nan_to_num(v), np.float32)})
                        for c, a, b in windows:
                            if c not in csz or c not in cidx:
                                continue
                            a, b = max(0, a - t["flank"]), min(csz[c], b + t["flank"])
                            n = max(1, (b - a) // t["bin"])
                            v = bw.values(c, a, a + n * t["bin"], bins=n, summary="mean", missing=0.0)
                            fine.append({"chrom": cidx[c], "start": a, "bin": t["bin"],
                                         "v": _pack(np.nan_to_num(v), np.float32)})
                    o["series"].append({"sample": s, "coarse": coarse, "fine": fine})
                    size += sum(len(w["v"]["b"]) for w in coarse + fine)
                parts[f"track: {t['name']}"] = size
            elif t["type"] == "intervals":
                o["series"], size = [], 0
                if t["color"]:
                    o["color"] = t["color"]
                for s, label, v in t["raw"]:
                    if isinstance(v, np.ndarray) and v.dtype == bool:
                        enc, _, _ = iv(L.codes[v], L.starts[v], L.ends[v])
                    elif hasattr(v, "codes"):
                        lut = np.array([g._add(n) for n in v.genome.names], np.int64)
                        if len(remap) < len(g.names):
                            remap = np.concatenate([remap, np.full(len(g.names) - len(remap), -1)])
                        enc, _, _ = iv(lut[v.codes], v.starts, v.ends)
                    else:
                        items = list(v)
                        enc, _, _ = iv(np.array([x.chrom for x in items], dtype=object),
                                       np.array([x.start for x in items]), np.array([x.end for x in items]))
                    o["series"].append({"sample": s, "label": label, **enc})
                    size += sum(len(x["b"]) for x in enc.values())
                parts[f"track: {t['name']}"] = size
            elif t["type"] == "points":
                o["series"], size = [], 0
                for s, ch, pos, val in t["raw"]:
                    codes = np.array([cidx.get(c, -1) for c in ch])
                    ok = (codes >= 0) & np.isfinite(val)
                    oo = np.lexsort((pos[ok], codes[ok]))
                    enc = {"chrom": _pack(codes[ok][oo], np.uint8), "pos": _pack(_delta(pos[ok][oo]), np.int32),
                           "v": _pack(val[ok][oo], np.float32)}
                    o["series"].append({"sample": s, **enc})
                    size += sum(len(x["b"]) for x in enc.values())
                o["segments"] = [{"sample": s, "chrom": cidx[c], "start": a, "end": b, "label": lab}
                                 for s, c, a, b, lab in t["segs"] if c in cidx]
                if t["ylim"]:
                    o["ylim"] = list(t["ylim"])
                parts[f"track: {t['name']}"] = size
            tracks.append(o)
        D["tracks"], D["creCols"] = tracks, creCols
        return D, parts

    def to_html(self, standalone: bool = True):
        if not self._tracks:                                  # a sensible default figure
            if self.A is not None:
                self.anchor_profile()
                if "strength" in self.A.vp:
                    self.cre_values("node strength", self.A.vp.strength)
            self.cres()
            if self.A is not None:
                self.loops()
            if self.G is not None:
                self.genes()
        D, parts = self._payload()
        js = (_HERE / "_view.js").read_text()
        css = (_HERE / "_view.css").read_text()
        body = _BODY.format(title=_html.escape(self.title), sub=_html.escape(self.subtitle))
        blob = json.dumps(D, separators=(",", ":")).replace("</", "<\\/")
        page = (f"<title>{_html.escape(self.title)}</title>\n<style>{css}</style>\n{body}\n"
                f'<script type="application/json" id="gb-data">{blob}</script>\n<script>{js}</script>\n')
        if standalone:
            page = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">\n'
                    '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                    + page.replace("</style>\n", "</style>\n</head><body>\n", 1) + "</body></html>\n")
        parts["viewer (js + css)"] = len(js) + len(css)
        parts["total"] = len(page)
        return page, parts

    def save(self, path: str, standalone: bool = True) -> dict:
        """Write the page; returns bytes per part."""
        page, parts = self.to_html(standalone)
        with open(path, "w") as f:
            f.write(page)
        return parts

    def _repr_html_(self):
        page, _ = self.to_html(True)
        return (f'<iframe srcdoc="{_html.escape(page, quote=True)}" style="width:100%;height:880px;border:0" '
                f'title="{_html.escape(self.title)}"></iframe>')


def view_html(path, *, architecture=None, cre=None, genes=None, signal=None, regions=(), notes=None,
              title="genomeblocks view", subtitle="", standalone=True, **kw) -> dict:
    """One-call version with the default figure (anchor profile, strength, signal, CREs, loops, genes)."""
    v = View(architecture, cre=cre, genes=genes, title=title, subtitle=subtitle, **kw)
    notes = notes or {}
    for r in regions:
        v.region(str(r), str(r), note=notes.get(str(r), ""))
    if architecture is not None:
        v.anchor_profile()
        if "strength" in architecture.vp:
            v.cre_values("node strength", architecture.vp.strength)
    for name, bw in (signal or {}).items():
        v.signal(name, bw)
    v.cres()
    if architecture is not None:
        v.loops()
    if genes is not None:
        v.genes()
    return v.save(path, standalone)


_BODY = """<div class="gbv">
<header><h1>{title}</h1><p class="sub">{sub}</p></header>
<div class="bar">
  <label class="vh" for="gbv-locus">Locus or gene</label>
  <input id="gbv-locus" type="text" spellcheck="false" autocomplete="off" placeholder="a gene name, or chr8:127,700,000-128,100,000">
  <button type="button" id="gbv-go">Go</button>
  <button type="button" id="gbv-out" aria-label="Zoom out">−</button>
  <button type="button" id="gbv-in" aria-label="Zoom in">+</button>
  <div id="gbv-samples" class="chips" aria-label="Samples (click to hide or show)"></div>
  <span class="hint">drag to pan · wheel or double-click to zoom · click a CRE or gene</span>
</div>
<div class="main">
  <aside class="side">
    <section><h2>Regions</h2><div id="gbv-regions" class="list"></div></section>
    <section><h2>Top hubs</h2><div id="gbv-hubs" class="list"></div></section>
  </aside>
  <div class="stage"><div id="gbv-panels" class="panels"></div>
    <section id="gbv-inspect" class="inspect" hidden></section></div>
</div>
<div id="gbv-tip" class="tip" hidden></div>
<p class="foot">Made with genomeblocks · all data is inside this file</p>
</div>"""
