/* genomeblocks view — a small canvas genome browser for genomeblocks tables.
   Data arrives as gzipped typed arrays in <script id="gb-data">; see view.py. */
(() => {
"use strict";
const $ = (id) => document.getElementById(id);
const DATA = JSON.parse($("gb-data").textContent);
const DPR = Math.max(1, Math.min(2, window.devicePixelRatio || 1));
const TYPES = { uint8: Uint8Array, uint16: Uint16Array, int32: Int32Array, uint32: Uint32Array, float32: Float32Array };
const FONT = "12px system-ui, -apple-system, 'Segoe UI', sans-serif";
const SMALL = "11px system-ui, -apple-system, 'Segoe UI', sans-serif";

// ── decoding ─────────────────────────────────────────────────────────────
async function unpack(p) {
  if (!p) return null;
  const s = atob(p.b), u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
  const stream = new Blob([u]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new TYPES[p.t](await new Response(stream).arrayBuffer());
}
function cumsum(a) { const o = new Int32Array(a.length); let s = 0; for (let i = 0; i < a.length; i++) { s += a[i]; o[i] = s; } return o; }
function blocks(chrom, n) {
  const lo = new Int32Array(n), hi = new Int32Array(n), seen = new Uint8Array(n);
  for (let i = 0; i < chrom.length; i++) { const c = chrom[i]; if (!seen[c]) { seen[c] = 1; lo[c] = i; } hi[c] = i + 1; }
  return { lo, hi };
}
function maxLens(chrom, start, end, n) {
  const m = new Int32Array(n);
  for (let i = 0; i < chrom.length; i++) m[chrom[i]] = Math.max(m[chrom[i]], end[i] - start[i]);
  return m;
}
function lowerBound(arr, x, lo, hi) { while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] < x) lo = m + 1; else hi = m; } return lo; }
function rangeOf(set, c, a, b) {          // rows of an interval set that may overlap [a, b) on chrom c
  const lo = set.blk.lo[c], hi = set.blk.hi[c];
  return [lowerBound(set.start, a - set.maxLen[c], lo, hi), lowerBound(set.start, b, lo, hi)];
}
async function intervals(p, nchrom) {
  const chrom = await unpack(p.chrom), start = cumsum(await unpack(p.start)), len = await unpack(p.len);
  const end = new Int32Array(start.length);
  for (let i = 0; i < start.length; i++) end[i] = start[i] + len[i];
  return { chrom, start, end, blk: blocks(chrom, nchrom), maxLen: maxLens(chrom, start, end, nchrom) };
}

// ── model ────────────────────────────────────────────────────────────────
const M = { chroms: DATA.chroms, sizes: DATA.sizes, samples: DATA.samples, labels: DATA.labels || [] };
M.code = new Map(M.chroms.map((c, i) => [c.toLowerCase(), i]));
const state = { panels: [], anchor: null, hidden: new Set(), hover: null };

async function load() {
  const nc = M.chroms.length;
  const C = DATA.cre;
  M.cre = await intervals(C, nc);
  const n = M.cre.start.length;
  M.cre.n = n;
  M.cre.center = new Int32Array(n);
  for (let i = 0; i < n; i++) M.cre.center[i] = (M.cre.start[i] + M.cre.end[i]) >> 1;
  M.cre.label = C.label ? await unpack(C.label) : null;
  M.cre.gene = C.gene ? await unpack(C.gene) : null;
  M.cre.geneNames = DATA.creGeneNames || [];
  M.cre.cols = [];
  for (const c of DATA.creCols || []) M.cre.cols.push({ name: c.name, sample: c.sample, v: await unpack(c.v) });
  M.promoter = M.labels.indexOf("Promoter-TSS");
  if (DATA.edges) {
    const E = DATA.edges, src = cumsum(await unpack(E.src)), dt = await unpack(E.dt);
    const m = src.length, tgt = new Int32Array(m);
    for (let i = 0; i < m; i++) tgt[i] = src[i] + dt[i];
    const scores = [];
    for (const s of E.scores) {
      const v = await unpack(s.v);
      const pos = Array.from(v).filter((x) => x > 0).sort((a, b) => a - b);
      scores.push({ sample: s.sample, v, top: pos.length ? pos[Math.floor(pos.length * 0.95)] : 1 });
    }
    const ptr = new Int32Array(n + 1);
    for (let i = 0; i < m; i++) { ptr[src[i] + 1]++; ptr[tgt[i] + 1]++; }
    for (let i = 0; i < n; i++) ptr[i + 1] += ptr[i];
    const cur = ptr.slice(0, n), nbr = new Int32Array(2 * m), eid = new Int32Array(2 * m);
    for (let i = 0; i < m; i++) {
      nbr[cur[src[i]]] = tgt[i]; eid[cur[src[i]]++] = i;
      nbr[cur[tgt[i]]] = src[i]; eid[cur[tgt[i]]++] = i;
    }
    M.edges = { m, src, tgt, scores, scoreName: E.scoreName, ptr, nbr, eid };
  }
  if (DATA.genes) {
    const G = DATA.genes;
    M.genes = await intervals(G, nc);
    M.genes.strand = await unpack(G.strand);
    M.genes.nex = await unpack(G.nex);
    M.genes.exStart = await unpack(G.exStart);
    M.genes.exLen = await unpack(G.exLen);
    M.genes.names = G.names;
    M.genes.tss = G.tss ? await unpack(G.tss) : null;
    const off = new Int32Array(M.genes.nex.length + 1);
    for (let i = 0; i < M.genes.nex.length; i++) off[i + 1] = off[i] + M.genes.nex[i];
    M.genes.exOff = off;
    M.genes.byName = new Map();
    G.names.forEach((nm, i) => { if (!M.genes.byName.has(nm.toLowerCase())) M.genes.byName.set(nm.toLowerCase(), i); });
  }
  M.tracks = [];
  for (const t of DATA.tracks) M.tracks.push(await loadTrack(t, nc));
}

async function loadTrack(t, nc) {
  const o = Object.assign({}, t);
  if (t.type === "signal") {
    o.series = [];
    for (const s of t.series) {
      const coarse = [];
      for (const w of s.coarse || []) coarse.push({ chrom: w.chrom, start: 0, bin: w.bin, v: await unpack(w.v) });
      const fine = [];
      for (const w of s.fine || []) fine.push({ chrom: w.chrom, start: w.start, bin: w.bin, v: await unpack(w.v) });
      o.series.push({ sample: s.sample, coarse, fine });
    }
  } else if (t.type === "intervals") {
    o.series = [];
    for (const s of t.series) o.series.push(Object.assign({ sample: s.sample, label: s.label }, await intervals(s, nc)));
  } else if (t.type === "points") {
    o.series = [];
    let lo = Infinity, hi = -Infinity;
    for (const s of t.series) {
      const chrom = await unpack(s.chrom), pos = cumsum(await unpack(s.pos)), v = await unpack(s.v);
      const sorted = Array.from(v).filter(Number.isFinite).sort((a, b) => a - b);
      if (sorted.length) { lo = Math.min(lo, sorted[Math.floor(sorted.length * 0.005)]); hi = Math.max(hi, sorted[Math.floor(sorted.length * 0.995)]); }
      o.series.push({ sample: s.sample, chrom, pos, v, blk: blocks(chrom, nc) });
    }
    o.ylo = t.ylim ? t.ylim[0] : Math.min(-2, lo); o.yhi = t.ylim ? t.ylim[1] : Math.max(2, hi);
    o.segments = t.segments || [];
  } else if (t.type === "creval") {
    o.series = t.series.map((s) => ({ sample: s.sample, v: M.cre.cols[s.col].v }));
  }
  return o;
}

// ── helpers ──────────────────────────────────────────────────────────────
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
function sampleColor(i) { return i == null || i < 0 || !M.samples[i] ? css("--gbv-ink") : M.samples[i].color; }
function sampleName(i) { return i == null || !M.samples[i] ? "" : M.samples[i].name; }
function fmtBp(x) {
  const a = Math.abs(x);
  if (a >= 1e6) return (x / 1e6).toFixed(a >= 1e7 ? 1 : 2).replace(/\.?0+$/, "") + " Mb";
  if (a >= 1e3) return (x / 1e3).toFixed(a >= 1e4 ? 0 : 1).replace(/\.0$/, "") + " kb";
  return Math.round(x) + " bp";
}
const fmtNum = (x) => (x === 0 ? "0" : Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : x.toPrecision(2));
const locStr = (p) => `${M.chroms[p.c]}:${Math.round(p.a + 1).toLocaleString()}-${Math.round(p.b).toLocaleString()}`;
function creUid(i) { return `${M.chroms[M.cre.chrom[i]]}:${M.cre.start[i]}-${M.cre.end[i]}`; }
function visibleSamples(series) { return series.filter((s) => s.sample == null || !state.hidden.has(s.sample)); }
function niceStep(span, n) {
  const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  return [1, 2, 5, 10].map((k) => k * mag).find((s) => span / s <= n) || 10 * mag;
}
function smooth(v, sigma) {
  if (sigma <= 0) return v;
  const r = Math.ceil(3 * sigma), k = [];
  let ks = 0;
  for (let i = -r; i <= r; i++) { const w = Math.exp(-(i * i) / (2 * sigma * sigma)); k.push(w); ks += w; }
  const out = new Float64Array(v.length);
  for (let i = 0; i < v.length; i++) {
    let s = 0;
    for (let j = -r; j <= r; j++) { const q = i + j; if (q >= 0 && q < v.length) s += v[q] * k[j + r]; }
    out[i] = s / ks;
  }
  return out;
}

// ── anchors (selected gene / CRE) ────────────────────────────────────────
function anchorFromGene(g) {
  const G = M.genes, c = G.chrom[g];
  const tss = G.tss ? G.tss[g] : (G.strand[g] === 2 ? G.end[g] : G.start[g]);
  const R = DATA.anchorR || 5000, lo = Math.max(0, tss - R), hi = tss + R, rows = [];
  const center = (DATA.anchorMode || "center") === "center";
  const [i0, i1] = rangeOf(M.cre, c, lo, hi);           // same rule as Architecture.support
  for (let i = i0; i < i1; i++) {
    const hit = center ? (M.cre.center[i] >= lo && M.cre.center[i] < hi) : (M.cre.start[i] < hi && M.cre.end[i] > lo);
    if (hit && (!M.edges || M.edges.ptr[i + 1] > M.edges.ptr[i])) rows.push(i);
  }
  return { kind: "gene", gene: g, name: G.names[g], c, rows, tss };
}
function anchorFromCre(i) { return { kind: "cre", name: creUid(i), c: M.cre.chrom[i], rows: [i] }; }
function partners(anchor) {            // [{p, e}] edges from the anchor to non-anchor CREs
  if (!M.edges || !anchor) return [];
  const set = new Set(anchor.rows), out = [];
  for (const v of anchor.rows) {
    for (let k = M.edges.ptr[v]; k < M.edges.ptr[v + 1]; k++) {
      const p = M.edges.nbr[k];
      if (!set.has(p)) out.push({ p, e: M.edges.eid[k] });
    }
  }
  return out;
}
function setAnchor(a) {
  state.anchor = a;
  if (a) { a.partners = partners(a); a.pset = new Map(a.partners.map((x) => [x.p, x])); }
  renderInspector();
  renderAll();
}

// ── panels ───────────────────────────────────────────────────────────────
const HEAD = 16;
const HEIGHTS = { anchor: 76, creval: 54, signal: 54, points: 84, cre: 26, loops: 118, genes: 92 };
function trackHeight(t) {
  if (t.height) return t.height;
  if (t.type === "intervals") return 6 + 12 * t.series.length;
  return HEIGHTS[t.type] || 50;
}
function layout() {
  let y = 34;
  const out = [];
  for (const t of M.tracks) { const h = trackHeight(t); out.push({ t, y: y + HEAD, h }); y += HEAD + h + 6; }
  return { rows: out, height: y + 4 };
}

function makePanel(c, a, b) {
  const el = document.createElement("div");
  el.className = "gbv-panel";
  const bar = document.createElement("div");
  bar.className = "gbv-ptitle";
  const loc = document.createElement("span");
  const close = document.createElement("button");
  close.type = "button"; close.textContent = "×"; close.setAttribute("aria-label", "Close this panel");
  bar.append(loc, close);
  const cv = document.createElement("canvas");
  cv.setAttribute("role", "img");
  el.append(bar, cv);
  const P = { c, a, b, el, cv, loc, hits: [] };
  close.addEventListener("click", () => { state.panels = state.panels.filter((q) => q !== P); el.remove(); syncPanels(); });
  bindPanel(P);
  return P;
}
function setPanels(loci) {
  const host = $("gbv-panels");
  host.replaceChildren();
  state.panels = loci.map(([c, a, b]) => makePanel(c, a, b));
  state.panels.forEach((p) => host.appendChild(p.el));
  syncPanels();
}
function syncPanels() {
  state.panels.forEach((p) => p.el.querySelector("button").hidden = state.panels.length < 2);
  $("gbv-locus").value = state.panels.map(locStr).join("  ");
  renderAll();
}
function clampPanel(P) {
  const size = M.sizes[P.c] || 1;
  let span = Math.max(200, Math.min(P.b - P.a, size));
  let a = Math.max(0, Math.min(P.a, size - span));
  P.a = a; P.b = a + span;
}
function renderAll() { state.panels.forEach(render); }

// ── drawing ──────────────────────────────────────────────────────────────
function render(P) {
  clampPanel(P);
  P.loc.textContent = `${locStr(P)} · ${fmtBp(P.b - P.a)}`;
  const W = Math.max(200, P.el.clientWidth - 2), L = layout();
  P.W = W; P.L = L;
  P.cv.width = W * DPR; P.cv.height = L.height * DPR;
  P.cv.style.width = W + "px"; P.cv.style.height = L.height + "px";
  P.cv.setAttribute("aria-label", `Tracks for ${locStr(P)}`);
  const ctx = P.cv.getContext("2d");
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  ctx.fillStyle = css("--gbv-canvas"); ctx.fillRect(0, 0, W, L.height);
  P.hits = [];
  const X = (pos) => ((pos - P.a) / (P.b - P.a)) * W;
  P.X = X;
  drawRuler(ctx, P, W);
  // anchor band
  if (state.anchor && state.anchor.rows.length && state.anchor.c === P.c) {
    let lo = Infinity, hi = -Infinity;
    for (const i of state.anchor.rows) { lo = Math.min(lo, M.cre.start[i]); hi = Math.max(hi, M.cre.end[i]); }
    const x0 = X(lo), x1 = Math.max(X(hi), x0 + 3);
    ctx.fillStyle = css("--gbv-anchor"); ctx.globalAlpha = 0.28;
    ctx.fillRect(x0, 30, x1 - x0, L.height - 30); ctx.globalAlpha = 1;
  }
  for (const r of L.rows) {
    ctx.save();
    ctx.beginPath(); ctx.rect(0, r.y - HEAD, W, r.h + HEAD); ctx.clip();
    const fn = DRAW[r.t.type];
    if (fn) fn(ctx, P, r, W, X);
    ctx.restore();
    ctx.font = "600 " + SMALL; ctx.fillStyle = css("--gbv-ink2"); ctx.textBaseline = "alphabetic";
    ctx.fillText(r.t.name || r.t.type, 4, r.y - 4);
  }
  drawMarks(ctx, P, L, X);
}

function drawRuler(ctx, P, W) {
  const span = P.b - P.a, step = niceStep(span, Math.max(2, Math.floor(W / 110)));
  ctx.strokeStyle = css("--gbv-grid"); ctx.fillStyle = css("--gbv-muted"); ctx.font = SMALL; ctx.textAlign = "center";
  ctx.beginPath(); ctx.moveTo(0, 28.5); ctx.lineTo(W, 28.5); ctx.stroke();
  for (let t = Math.ceil(P.a / step) * step; t <= P.b; t += step) {
    const x = Math.round(P.X(t)) + 0.5;
    ctx.beginPath(); ctx.moveTo(x, 22); ctx.lineTo(x, 28); ctx.stroke();
    if (x > 30 && x < W - 30) ctx.fillText(fmtBp(t).replace(" ", " "), x, 17);
  }
  ctx.textAlign = "left";
}

function overlay(ctx, r, W, series, label) {    // series: [{color, v (per pixel-bin), bw}]
  let ymax = 0;
  for (const s of series) for (const x of s.v) if (x > ymax) ymax = x;
  const y0 = r.y + r.h;
  ctx.strokeStyle = css("--gbv-grid"); ctx.beginPath(); ctx.moveTo(0, y0 + 0.5); ctx.lineTo(W, y0 + 0.5); ctx.stroke();
  ctx.font = SMALL; ctx.fillStyle = css("--gbv-muted"); ctx.textAlign = "right";
  ctx.fillText(ymax > 0 ? fmtNum(ymax) : (label || "0"), W - 4, r.y + 10); ctx.textAlign = "left";
  if (ymax <= 0) return;
  const order = series.map((s, i) => [s.v.reduce((a, b) => a + b, 0), i]).sort((p, q) => q[0] - p[0]).map((p) => p[1]);
  for (const i of order) {
    const s = series[i], n = s.v.length, bw = W / n;
    ctx.beginPath(); ctx.moveTo(0, y0);
    for (let k = 0; k < n; k++) ctx.lineTo((k + 0.5) * bw, y0 - (s.v[k] / ymax) * (r.h - 4));
    ctx.lineTo(W, y0); ctx.closePath();
    ctx.fillStyle = s.color; ctx.globalAlpha = series.length > 1 ? 0.35 : 0.8; ctx.fill(); ctx.globalAlpha = 1;
    if (series.length > 1) {
      ctx.beginPath();
      for (let k = 0; k < n; k++) { const x = (k + 0.5) * bw, y = y0 - (s.v[k] / ymax) * (r.h - 4); k ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }
      ctx.strokeStyle = s.color; ctx.lineWidth = 1; ctx.stroke();
    }
  }
}
function note(ctx, r, W, text) {
  ctx.font = SMALL; ctx.fillStyle = css("--gbv-muted"); ctx.textAlign = "center";
  ctx.fillText(text, W / 2, r.y + r.h / 2 + 4); ctx.textAlign = "left";
}

const DRAW = {
  signal(ctx, P, r, W) {
    const nb = Math.max(1, Math.floor(W)), bpp = (P.b - P.a) / nb, out = [];
    let any = false;
    for (const s of visibleSamples(r.t.series)) {
      const fine = s.fine.filter((w) => w.chrom === P.c && w.start <= P.a && w.start + w.v.length * w.bin >= P.b);
      const src = fine.length && fine[0].bin <= bpp * 4 ? fine[0] : s.coarse.find((w) => w.chrom === P.c);
      const v = new Float64Array(nb);
      if (src) {
        any = true;
        for (let k = 0; k < nb; k++) {
          const a = P.a + k * bpp, b = a + bpp;
          let i0 = Math.floor((a - src.start) / src.bin), i1 = Math.ceil((b - src.start) / src.bin);
          i0 = Math.max(0, i0); i1 = Math.min(src.v.length, Math.max(i1, i0 + 1));
          let m = 0;
          for (let i = i0; i < i1; i++) if (src.v[i] > m) m = src.v[i];
          v[k] = m;
        }
      }
      out.push({ color: sampleColor(s.sample), v });
    }
    if (!any) return note(ctx, r, W, "no signal exported for this chromosome");
    overlay(ctx, r, W, out);
  },
  creval(ctx, P, r, W) {
    const nb = Math.max(20, Math.floor(W / 3)), bpp = (P.b - P.a) / nb, out = [];
    const [i0, i1] = rangeOf(M.cre, P.c, P.a, P.b);
    for (const s of visibleSamples(r.t.series)) {
      const v = new Float64Array(nb);
      for (let i = i0; i < i1; i++) {
        const k = Math.floor((M.cre.center[i] - P.a) / bpp);
        if (k >= 0 && k < nb) v[k] += s.v[i];
      }
      out.push({ color: sampleColor(s.sample), v: smooth(v, r.t.smooth == null ? 1 : r.t.smooth) });
    }
    overlay(ctx, r, W, out);
  },
  anchor(ctx, P, r, W) {
    const A = state.anchor;
    if (!A || !M.edges) return note(ctx, r, W, "search a gene or click a CRE to see its contact profile");
    const nb = Math.max(20, Math.floor(W / 3)), bpp = (P.b - P.a) / nb, out = [];
    M.edges.scores.forEach((sc, j) => {
      if (state.hidden.has(sc.sample)) return;
      const v = new Float64Array(nb);
      for (const { p, e } of A.partners) {
        if (M.cre.chrom[p] !== P.c) continue;
        const k = Math.floor((M.cre.center[p] - P.a) / bpp);
        if (k >= 0 && k < nb && sc.v[e] > 0) v[k] += sc.v[e];
      }
      out.push({ color: sampleColor(sc.sample), v: smooth(v, r.t.smooth == null ? 1.5 : r.t.smooth) });
    });
    overlay(ctx, r, W, out);
    if (!out.some((s) => s.v.some((x) => x > 0))) note(ctx, r, W, `${A.name}: no contacts with weight in this view`);
  },
  intervals(ctx, P, r, W, X) {
    const vis = r.t.series;
    vis.forEach((s, j) => {
      if (s.sample != null && state.hidden.has(s.sample)) return;
      const y = r.y + 3 + j * 12, [i0, i1] = rangeOf(s, P.c, P.a, P.b);
      ctx.fillStyle = s.sample != null ? sampleColor(s.sample) : (r.t.color || css("--gbv-cre"));
      for (let i = i0; i < i1; i++) {
        if (s.end[i] <= P.a) continue;
        const x0 = X(s.start[i]), x1 = X(s.end[i]);
        ctx.fillRect(x0, y, Math.max(2, x1 - x0), 9);
      }
      if (s.label) { ctx.font = SMALL; ctx.fillStyle = css("--gbv-muted"); ctx.textAlign = "right"; ctx.fillText(s.label, W - 4, y + 8); ctx.textAlign = "left"; }
    });
  },
  points(ctx, P, r, W, X) {
    const t = r.t, y = (v) => r.y + r.h - ((Math.max(t.ylo, Math.min(t.yhi, v)) - t.ylo) / (t.yhi - t.ylo)) * (r.h - 4) - 2;
    let segTop = 0;
    const segs = t.segments.filter((s) => s.chrom === P.c && s.end > P.a && s.start < P.b && !state.hidden.has(s.sample));
    ctx.setLineDash([3, 3]); ctx.strokeStyle = css("--gbv-muted");
    ctx.beginPath(); ctx.moveTo(0, Math.round(y(0)) + 0.5); ctx.lineTo(W, Math.round(y(0)) + 0.5); ctx.stroke(); ctx.setLineDash([]);
    for (const s of visibleSamples(t.series)) {
      const lo = s.blk.lo[P.c], hi = s.blk.hi[P.c];
      const i0 = lowerBound(s.pos, P.a, lo, hi), i1 = lowerBound(s.pos, P.b, lo, hi);
      ctx.fillStyle = sampleColor(s.sample);
      for (let i = i0; i < i1; i++) if (Number.isFinite(s.v[i])) ctx.fillRect(X(s.pos[i]) - 1.2, y(s.v[i]) - 1.2, 2.4, 2.4);
    }
    segs.forEach((s, k) => {
      const yy = r.y + 6 + (k % 4) * 9, x0 = Math.max(0, X(s.start)), x1 = Math.min(W, X(s.end));
      ctx.fillStyle = sampleColor(s.sample); ctx.fillRect(x0, yy, Math.max(2, x1 - x0), 4);
      if (s.label) { ctx.font = SMALL; ctx.fillStyle = css("--gbv-ink"); ctx.fillText(s.label, Math.min(x1 + 3, W - 24), yy + 5); }
      segTop++;
    });
    ctx.font = SMALL; ctx.fillStyle = css("--gbv-muted"); ctx.textAlign = "right";
    ctx.fillText(`${fmtNum(t.yhi)} / ${fmtNum(t.ylo)}`, W - 4, r.y + r.h - 3); ctx.textAlign = "left";
  },
  cre(ctx, P, r, W, X) {
    const [i0, i1] = rangeOf(M.cre, P.c, P.a, P.b), A = state.anchor;
    const aset = A ? new Set(A.rows) : null;
    const colors = { cre: css("--gbv-cre"), prom: css("--gbv-prom"), anchor: css("--gbv-anchor-ink"),
                     pp: css("--gbv-partner-p"), pe: css("--gbv-partner-e") };
    for (let i = i0; i < i1; i++) {
      if (M.cre.end[i] <= P.a) continue;
      const x0 = X(M.cre.start[i]), x1 = X(M.cre.end[i]), w = Math.max(2, x1 - x0);
      const isProm = M.cre.label && M.cre.label[i] === M.promoter;
      let col = isProm ? colors.prom : colors.cre, y = r.y + 7, h = 12;
      if (A) {
        if (aset.has(i)) { col = colors.anchor; y = r.y + 2; h = 22; }
        else if (A.pset.has(i)) { col = isProm ? colors.pp : colors.pe; y = r.y + 4; h = 18; }
        else { ctx.globalAlpha = 0.35; }
      }
      ctx.fillStyle = col; ctx.fillRect(x0, y, w, h); ctx.globalAlpha = 1;
      P.hits.push({ kind: "cre", i, x0: x0 - 2, x1: x0 + w + 2, y0: r.y, y1: r.y + r.h });
    }
  },
  loops(ctx, P, r, W, X) {
    if (!M.edges) return note(ctx, r, W, "no loops");
    const A = state.anchor, base = r.y + r.h - 2, E = M.edges, span = P.b - P.a;
    const [i0, i1] = rangeOf(M.cre, P.c, P.a - span, P.b + span);
    const seen = new Set();
    const edges = [];
    if (A) for (const { e } of A.partners) edges.push(e);
    else for (let i = i0; i < i1; i++) for (let k = E.ptr[i]; k < E.ptr[i + 1]; k++) { const e = E.eid[k]; if (!seen.has(e)) { seen.add(e); edges.push(e); } }
    const cis = [], trans = [];
    for (const e of edges) {
      const cs = M.cre.chrom[E.src[e]], ct = M.cre.chrom[E.tgt[e]];
      if (cs !== ct) trans.push(e); else if (cs === P.c) cis.push(e);
    }
    ctx.lineWidth = A ? 1.6 : 1;
    E.scores.forEach((sc) => {
      if (state.hidden.has(sc.sample)) return;
      ctx.strokeStyle = sampleColor(sc.sample);
      for (const e of cis) {
        const val = sc.v[e];
        if (!(val > 0) && !A) continue;
        const xa = X(M.cre.center[E.src[e]]), xb = X(M.cre.center[E.tgt[e]]);
        if ((xa < 0 && xb < 0) || (xa > W && xb > W)) continue;
        const h = Math.min(r.h - 6, Math.abs(xb - xa) / 2);
        ctx.globalAlpha = Math.max(0.12, Math.min(0.9, val / sc.top));
        ctx.beginPath(); ctx.moveTo(xa, base); ctx.bezierCurveTo(xa, base - h * 1.33, xb, base - h * 1.33, xb, base); ctx.stroke();
      }
    });
    ctx.globalAlpha = 1; ctx.lineWidth = 1;
    ctx.font = SMALL;
    for (const e of trans) {                      // trans: a stub labelled with the partner chromosome
      const s = E.src[e], t = E.tgt[e];
      const here = M.cre.chrom[s] === P.c ? s : M.cre.chrom[t] === P.c ? t : -1;
      if (here < 0) continue;
      const x = X(M.cre.center[here]);
      if (x < 0 || x > W) continue;
      const other = here === s ? t : s;
      ctx.strokeStyle = css("--gbv-trans"); ctx.fillStyle = css("--gbv-trans");
      ctx.beginPath(); ctx.moveTo(x, base); ctx.lineTo(x, r.y + 16); ctx.stroke();
      ctx.fillText(M.chroms[M.cre.chrom[other]], Math.min(x + 3, W - 34), r.y + 14);
      P.hits.push({ kind: "trans", here, other, e, x0: x - 5, x1: x + 30, y0: r.y, y1: r.y + r.h });
    }
  },
  genes(ctx, P, r, W, X) {
    if (!M.genes) return;
    const G = M.genes, [i0, i1] = rangeOf(G, P.c, P.a, P.b), rowsEnd = [], rowH = 22, maxRows = Math.floor(r.h / rowH);
    ctx.font = SMALL;
    let hidden = 0;
    const sel = state.anchor && state.anchor.kind === "gene" ? state.anchor.gene : -1;
    for (let g = i0; g < i1; g++) {
      if (G.end[g] <= P.a) continue;
      const x0 = X(G.start[g]), x1 = X(G.end[g]), name = G.names[g], tw = ctx.measureText(name).width;
      const left = Math.min(x0, (x0 + x1) / 2 - tw / 2), right = Math.max(x1, (x0 + x1) / 2 + tw / 2);
      let row = rowsEnd.findIndex((e) => e < left - 6);
      if (row < 0) { if (rowsEnd.length >= maxRows) { hidden++; continue; } row = rowsEnd.length; rowsEnd.push(0); }
      rowsEnd[row] = right;
      const y = r.y + row * rowH + 6;
      const col = g === sel ? css("--gbv-anchor-ink") : css("--gbv-gene");
      ctx.strokeStyle = col; ctx.fillStyle = col;
      ctx.beginPath(); ctx.moveTo(x0, y + 4.5); ctx.lineTo(x1, y + 4.5); ctx.stroke();
      if (x1 - x0 > 30) {                        // strand chevrons
        const dir = G.strand[g] === 2 ? -1 : 1;
        for (let x = x0 + 10; x < x1 - 6; x += 24) { ctx.beginPath(); ctx.moveTo(x - 2 * dir, y + 2); ctx.lineTo(x + 1 * dir, y + 4.5); ctx.lineTo(x - 2 * dir, y + 7); ctx.stroke(); }
      }
      for (let k = G.exOff[g]; k < G.exOff[g + 1]; k++) {
        const a = G.start[g] + G.exStart[k], b = a + G.exLen[k];
        ctx.fillRect(X(a), y, Math.max(1.5, X(b) - X(a)), 9);
      }
      ctx.font = (g === sel ? "700 " : "") + SMALL;
      ctx.fillText(name, Math.max(2, (x0 + x1) / 2 - tw / 2), y + 20);
      ctx.font = SMALL;
      P.hits.push({ kind: "gene", g, x0: Math.min(x0, left), x1: Math.max(x1, right), y0: y - 2, y1: y + 22 });
    }
    if (hidden) { ctx.fillStyle = css("--gbv-muted"); ctx.textAlign = "right"; ctx.fillText(`+${hidden} more genes: zoom in`, W - 4, r.y + r.h - 2); ctx.textAlign = "left"; }
  },
};

function drawMarks(ctx, P, L, X) {
  const y = L.rows.length ? L.rows[0].y - HEAD + 2 : 30;
  for (const m of DATA.marks || []) {
    if (M.code.get(m.chrom.toLowerCase()) !== P.c || m.pos < P.a || m.pos > P.b) continue;
    const x = X(m.pos);
    ctx.fillStyle = css("--gbv-ink"); ctx.strokeStyle = css("--gbv-ink");
    ctx.beginPath(); ctx.moveTo(x, y + 12); ctx.lineTo(x - 5, y + 4); ctx.lineTo(x + 5, y + 4); ctx.closePath(); ctx.fill();
    ctx.font = "700 " + FONT; ctx.textAlign = "center"; ctx.fillText(m.label, Math.max(30, Math.min(P.W - 30, x)), y + 1); ctx.textAlign = "left";
  }
}

// ── interaction ──────────────────────────────────────────────────────────
function hitAt(P, x, y) {
  for (let k = P.hits.length - 1; k >= 0; k--) {
    const h = P.hits[k];
    if (x >= h.x0 && x <= h.x1 && y >= h.y0 && y <= h.y1) return h;
  }
  return null;
}
function bindPanel(P) {
  const cv = P.cv;
  let drag = null;
  cv.addEventListener("pointerdown", (ev) => {
    drag = { x: ev.clientX, a: P.a, b: P.b, moved: false };
    cv.setPointerCapture(ev.pointerId);
  });
  cv.addEventListener("pointermove", (ev) => {
    const rect = cv.getBoundingClientRect(), x = ev.clientX - rect.left, y = ev.clientY - rect.top;
    if (drag) {
      const dx = ev.clientX - drag.x;
      if (Math.abs(dx) > 3) drag.moved = true;
      if (drag.moved) {
        const shift = (dx / P.W) * (drag.b - drag.a);
        P.a = drag.a - shift; P.b = drag.b - shift;
        requestAnimationFrame(() => { render(P); $("gbv-locus").value = state.panels.map(locStr).join("  "); });
      }
      return;
    }
    showTip(hitAt(P, x, y), ev);
  });
  cv.addEventListener("pointerup", (ev) => {
    const d = drag; drag = null;
    if (d && !d.moved) {
      const rect = cv.getBoundingClientRect(), h = hitAt(P, ev.clientX - rect.left, ev.clientY - rect.top);
      if (h) onHit(P, h);
    }
  });
  cv.addEventListener("pointerleave", () => showTip(null));
  cv.addEventListener("wheel", (ev) => {
    ev.preventDefault();
    const rect = cv.getBoundingClientRect(), f = Math.exp(ev.deltaY * 0.0015);
    zoomAt(P, (ev.clientX - rect.left) / P.W, f);
  }, { passive: false });
  cv.addEventListener("dblclick", (ev) => {
    const rect = cv.getBoundingClientRect();
    zoomAt(P, (ev.clientX - rect.left) / P.W, 0.5);
  });
}
function zoomAt(P, frac, f) {
  const pos = P.a + frac * (P.b - P.a), span = Math.max(200, (P.b - P.a) * f);
  P.a = pos - frac * span; P.b = P.a + span;
  render(P); $("gbv-locus").value = state.panels.map(locStr).join("  ");
}
function onHit(P, h) {
  if (h.kind === "cre") setAnchor(anchorFromCre(h.i));
  else if (h.kind === "gene") setAnchor(anchorFromGene(h.g));
  else if (h.kind === "trans") openPair(h.here, h.other);
}
function openPair(i, j) {
  const win = (k) => [M.cre.chrom[k], M.cre.center[k] - 100000, M.cre.center[k] + 100000];
  setPanels([win(i), win(j)]);
}
function showTip(h, ev) {
  const tip = $("gbv-tip");
  if (!h) { tip.hidden = true; return; }
  const lines = [];
  if (h.kind === "cre" || h.kind === "trans") {
    const i = h.kind === "cre" ? h.i : h.here;
    lines.push(["b", creUid(i)]);
    if (M.cre.label) lines.push(["", M.labels[M.cre.label[i]]]);
    if (M.cre.gene && M.cre.gene[i] >= 0) lines.push(["", "gene: " + M.cre.geneNames[M.cre.gene[i]]]);
    for (const c of M.cre.cols) lines.push(["", `${[c.name, sampleName(c.sample)].filter(Boolean).join(" · ")}: ${fmtNum(c.v[i])}`]);
    if (M.edges) lines.push(["", `${M.edges.ptr[i + 1] - M.edges.ptr[i]} partners · click to inspect`]);
    if (h.kind === "trans") lines.push(["b", `trans partner on ${M.chroms[M.cre.chrom[h.other]]} · click to open side by side`]);
  } else if (h.kind === "gene") {
    const G = M.genes;
    lines.push(["b", G.names[h.g]], ["", `${M.chroms[G.chrom[h.g]]}:${G.start[h.g].toLocaleString()}-${G.end[h.g].toLocaleString()} (${G.strand[h.g] === 2 ? "−" : "+"})`], ["", "click to anchor at its TSS"]);
  }
  tip.replaceChildren(...lines.map(([k, t]) => { const d = document.createElement(k ? "strong" : "span"); d.textContent = t; return d; }));
  tip.hidden = false;
  const x = ev.clientX + 14, y = ev.clientY + 14;
  tip.style.left = Math.min(window.innerWidth - tip.offsetWidth - 8, x) + "px";
  tip.style.top = Math.min(window.innerHeight - tip.offsetHeight - 8, y) + "px";
}

// ── search ───────────────────────────────────────────────────────────────
function parseLocus(s) {
  const m = s.trim().replace(/,/g, "").match(/^([^:\s]+):\s*([\d.]+)\s*(mb|kb)?\s*-\s*([\d.]+)\s*(mb|kb)?$/i);
  if (!m) return null;
  const c = M.code.get(m[1].toLowerCase());
  if (c == null) return null;
  const unit = (u) => (u ? (u.toLowerCase() === "mb" ? 1e6 : 1e3) : 1);
  const u = m[5] || m[3];
  return [c, parseFloat(m[2]) * unit(m[3] || u) - 1, parseFloat(m[4]) * unit(u)];
}
function go(text) {
  const parts = text.trim().split(/\s{2,}|\s+(?=[^\s:]+:)/).filter(Boolean);
  const loci = parts.map(parseLocus);
  if (loci.length && loci.every(Boolean)) { setPanels(loci); return true; }
  if (M.genes) {
    const key = text.trim().toLowerCase();
    let g = M.genes.byName.get(key);
    if (g == null) for (const [k, v] of M.genes.byName) if (k.startsWith(key)) { g = v; break; }
    if (g != null) { openGene(g); return true; }
  }
  const one = M.code.get(text.trim().toLowerCase());
  if (one != null) { setPanels([[one, 0, M.sizes[one]]]); return true; }
  return false;
}
function openGene(g, half) {
  const G = M.genes, tss = G.tss ? G.tss[g] : (G.strand[g] === 2 ? G.end[g] : G.start[g]), h = half || DATA.geneHalf || 500000;
  setPanels([[G.chrom[g], tss - h, tss + h]]);
  setAnchor(anchorFromGene(g));
}

// ── side panel: regions, hubs, inspector ─────────────────────────────────
function listButton(title, note, onClick) {
  const b = document.createElement("button");
  b.type = "button"; b.className = "gbv-item";
  const t = document.createElement("span"); t.className = "t"; t.textContent = title;
  b.append(t);
  if (note) { const n = document.createElement("span"); n.className = "n"; n.textContent = note; b.append(n); }
  b.addEventListener("click", onClick);
  return b;
}
function renderSide() {
  const reg = $("gbv-regions");
  reg.replaceChildren(...(DATA.regions || []).map((r) => listButton(r.label, r.note, () => {
    setPanels(r.loci.map(([c, a, b]) => [M.code.get(c.toLowerCase()), a, b]));
    if (r.gene && M.genes) { const g = M.genes.byName.get(r.gene.toLowerCase()); if (g != null) setAnchor(anchorFromGene(g)); }
    else if (r.anchorRow != null) setAnchor(anchorFromCre(r.anchorRow));
  })));
  const hubs = $("gbv-hubs");
  hubs.replaceChildren(...(DATA.hubs || []).map((h) => listButton(h.label, h.note, () => {
    const i = h.row;
    setPanels([[M.cre.chrom[i], M.cre.center[i] - 250000, M.cre.center[i] + 250000]]);
    setAnchor(anchorFromCre(i));
  })));
  const chips = $("gbv-samples");
  chips.replaceChildren(...M.samples.map((s, i) => {
    const b = document.createElement("button");
    b.type = "button"; b.className = "gbv-chip"; b.setAttribute("aria-pressed", "true");
    const dot = document.createElement("i"); dot.style.background = s.color;
    b.append(dot, document.createTextNode(s.name));
    b.addEventListener("click", () => {
      if (state.hidden.has(i)) state.hidden.delete(i); else state.hidden.add(i);
      b.setAttribute("aria-pressed", String(!state.hidden.has(i)));
      renderAll(); renderInspector();
    });
    return b;
  }));
  chips.hidden = M.samples.length < 1;
}
function renderInspector() {
  const box = $("gbv-inspect"), A = state.anchor;
  if (!A) { box.hidden = true; return; }
  box.hidden = false;
  const head = document.createElement("div"); head.className = "ihead";
  const h = document.createElement("h2");
  h.textContent = A.kind === "gene" ? `${A.name} · anchor = CREs within ${fmtBp(DATA.anchorR || 5000)} of the TSS` : A.name;
  const clear = document.createElement("button"); clear.type = "button"; clear.textContent = "Clear"; clear.className = "gbv-small";
  clear.addEventListener("click", () => setAnchor(null));
  head.append(h, clear);
  const sub = document.createElement("p"); sub.className = "isub";
  const vis = M.edges ? M.edges.scores.map((s, j) => j).filter((j) => !state.hidden.has(M.edges.scores[j].sample)) : [];
  const agg = new Map();                       // one row per partner: O/E summed over its edges to the anchor
  for (const { p, e } of A.partners) {
    let r = agg.get(p);
    if (!r) { r = { p, s: M.edges.scores.map(() => 0), n: 0 }; agg.set(p, r); }
    M.edges.scores.forEach((sc, j) => { r.s[j] += sc.v[e] || 0; });
    r.n++;
  }
  const key = vis.length ? vis[0] : 0;
  const rows = Array.from(agg.values()).sort((x, y) => y.s[key] - x.s[key]);
  sub.textContent = `${A.rows.length} anchor CRE${A.rows.length === 1 ? "" : "s"} · ${rows.length} partners ` +
    `(${rows.filter((r) => M.cre.chrom[r.p] !== A.c).length} on other chromosomes) · ` +
    vis.map((j) => `Σ ${M.edges.scoreName} ${sampleName(M.edges.scores[j].sample)} ${fmtNum(rows.reduce((a, r) => a + r.s[j], 0))}`).join(" · ");
  const tbl = document.createElement("table");
  const thead = document.createElement("tr");
  for (const c of ["partner", "where", "label", ...vis.map((j) => [M.edges.scoreName, sampleName(M.edges.scores[j].sample)].filter(Boolean).join(" · ")), ""]) {
    const th = document.createElement("th"); th.textContent = c; thead.append(th);
  }
  tbl.append(thead);
  const anchorPos = A.rows.length ? M.cre.center[A.rows[0]] : 0;
  for (const { p, s } of rows.slice(0, 40)) {
    const tr = document.createElement("tr");
    const cells = [creUid(p), M.cre.chrom[p] === A.c ? fmtBp(Math.abs(M.cre.center[p] - anchorPos)) : `trans · ${M.chroms[M.cre.chrom[p]]}`,
                   M.cre.label ? M.labels[M.cre.label[p]] + (M.cre.gene && M.cre.gene[p] >= 0 ? ` · ${M.cre.geneNames[M.cre.gene[p]]}` : "") : "",
                   ...vis.map((j) => fmtNum(s[j]))];
    cells.forEach((c, k) => { const td = document.createElement("td"); td.textContent = c; if (k >= 3) td.className = "n"; tr.append(td); });
    const td = document.createElement("td");
    const b = document.createElement("button"); b.type = "button"; b.className = "gbv-small";
    if (M.cre.chrom[p] === A.c) { b.textContent = "Go"; b.addEventListener("click", () => { setPanels([[M.cre.chrom[p], M.cre.center[p] - 100000, M.cre.center[p] + 100000]]); }); }
    else { b.textContent = "Side by side"; b.addEventListener("click", () => openPair(A.rows[0], p)); }
    td.append(b); tr.append(td); tbl.append(tr);
  }
  const wrap = document.createElement("div"); wrap.className = "itbl"; wrap.append(tbl);
  const more = document.createElement("p"); more.className = "isub";
  more.textContent = rows.length > 40 ? `showing the 40 strongest of ${rows.length}` : "";
  box.replaceChildren(head, sub, wrap, more);
}

// ── boot ─────────────────────────────────────────────────────────────────
async function boot() {
  const stage = $("gbv-panels");
  try {
    await load();
  } catch (e) {
    stage.textContent = "This browser could not open the embedded data (" + e + "). Try a current Chrome, Edge, Firefox or Safari.";
    document.body.dataset.ready = "error";
    return;
  }
  renderSide();
  $("gbv-go").addEventListener("click", () => go($("gbv-locus").value));
  $("gbv-locus").addEventListener("keydown", (ev) => { if (ev.key === "Enter") go($("gbv-locus").value); });
  $("gbv-in").addEventListener("click", () => state.panels.forEach((p) => zoomAt(p, 0.5, 0.5)));
  $("gbv-out").addEventListener("click", () => state.panels.forEach((p) => zoomAt(p, 0.5, 2)));
  new ResizeObserver(() => renderAll()).observe(stage);
  const mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.addEventListener) mq.addEventListener("change", renderAll);
  new MutationObserver(renderAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  const hash = decodeURIComponent(location.hash.slice(1));
  const first = (DATA.regions || [])[0];
  if (!(hash && go(hash))) {
    if (first) document.querySelector("#gbv-regions .gbv-item").click();
    else setPanels([[0, 0, Math.min(M.sizes[0], 2e6)]]);
  }
  document.body.dataset.ready = "1";
}
boot();
})();
