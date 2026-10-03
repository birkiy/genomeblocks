/* Tiny SVG chart kit for the benchmark report: horizontal bars (emphasis)
   and multi-series log/linear line charts. Colours come from CSS tokens so
   both themes work; every label goes in via textContent. */
(function () {
  const NS = "http://www.w3.org/2000/svg";
  const tip = document.createElement("div");
  tip.className = "tip"; tip.hidden = true; tip.setAttribute("role", "status");
  document.addEventListener("DOMContentLoaded", () => document.body.appendChild(tip));

  function el(tag, attrs, parent) {
    const n = document.createElementNS(NS, tag);
    for (const k in attrs || {}) n.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(n);
    return n;
  }
  function txt(parent, x, y, s, cls, anchor) {
    const t = el("text", { x, y, class: cls || "", "text-anchor": anchor || "start" }, parent);
    t.textContent = s;
    return t;
  }
  function showTip(evt, value, label) {
    tip.replaceChildren();
    const v = document.createElement("strong"); v.textContent = value;
    const l = document.createElement("span"); l.textContent = label;
    tip.append(v, l);
    tip.hidden = false;
    const r = (evt.target.getBoundingClientRect ? evt.target.getBoundingClientRect() : null);
    const x = evt.clientX || (r ? r.left + r.width / 2 : 0);
    const y = evt.clientY || (r ? r.top : 0);
    const w = tip.offsetWidth;
    tip.style.left = Math.min(window.innerWidth - w - 12, Math.max(12, x + 14)) + "px";
    tip.style.top = (y + window.scrollY - 12 - tip.offsetHeight) + "px";
  }
  function hideTip() { tip.hidden = true; }
  function bindTip(node, value, label) {
    node.setAttribute("tabindex", "0");
    node.setAttribute("aria-label", label + ": " + value);
    node.addEventListener("pointermove", e => showTip(e, value, label));
    node.addEventListener("focus", e => showTip(e, value, label));
    node.addEventListener("pointerleave", hideTip);
    node.addEventListener("blur", hideTip);
  }

  // ── scales & ticks ────────────────────────────────────────────────────────
  function scale(d0, d1, r0, r1, log) {
    if (log) {
      const a = Math.log10(d0), b = Math.log10(d1);
      return v => r0 + (Math.log10(v) - a) / (b - a) * (r1 - r0);
    }
    return v => r0 + (v - d0) / (d1 - d0) * (r1 - r0);
  }
  function logTicks(lo, hi) {
    const out = [];
    for (let p = Math.floor(Math.log10(lo)); p <= Math.ceil(Math.log10(hi)); p++) {
      const v = Math.pow(10, p);
      if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
    }
    return out;
  }
  function niceTicks(lo, hi, n) {
    const span = hi - lo, step0 = span / (n || 5);
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= (n || 5));
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }

  // ── horizontal emphasis bars ──────────────────────────────────────────────
  // spec: {rows:[{label, value, hl, tip?}], log, fmt, axis, sort:'asc'|'desc'|null}
  function hbar(host, spec) {
    const rows = spec.rows.slice();
    if (spec.sort === "asc") rows.sort((a, b) => a.value - b.value);
    if (spec.sort === "desc") rows.sort((a, b) => b.value - a.value);
    const W = 760, labelW = spec.labelW || 290, rowH = 30, top = 6, axisH = 34;
    const H = top + rows.length * rowH + axisH;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, class: "chart", role: "img",
                            "aria-label": spec.aria || "" }, null);
    const vals = rows.map(r => r.value);
    let lo = Math.min(...vals), hi = Math.max(...vals);
    const x0 = labelW, x1 = W - 92;
    let sx, ticks;
    if (spec.log) {
      lo = Math.pow(10, Math.floor(Math.log10(lo) - 0.15));
      hi = Math.pow(10, Math.ceil(Math.log10(hi) + 0.05));
      sx = scale(lo, hi, x0, x1, true); ticks = logTicks(lo, hi);
    } else {
      lo = 0; ticks = niceTicks(0, hi * 1.02, 5); hi = ticks[ticks.length - 1];
      sx = scale(0, hi, x0, x1, false);
    }
    const base = spec.log ? x0 : sx(0);
    const plotB = top + rows.length * rowH;
    for (const t of ticks) {
      el("line", { x1: sx(t), x2: sx(t), y1: top, y2: plotB, class: "grid" }, svg);
      txt(svg, sx(t), plotB + 16, spec.fmt(t), "tick", "middle");
    }
    if (spec.axis) txt(svg, (x0 + x1) / 2, plotB + 31, spec.axis, "axis-label", "middle");
    el("line", { x1: base, x2: base, y1: top, y2: plotB, class: "baseline" }, svg);
    rows.forEach((r, i) => {
      const y = top + i * rowH + rowH / 2;
      const lab = txt(svg, labelW - 12, y + 4, r.label, r.hl ? "rlabel hl" : "rlabel", "end");
      const xe = Math.max(sx(r.value), base + 2);
      const h = 16, rad = 4;
      const w = xe - base;
      const p = `M${base},${y - h / 2} h${Math.max(0, w - rad)} a${rad},${rad} 0 0 1 ${rad},${rad}` +
                ` v${h - 2 * rad} a${rad},${rad} 0 0 1 -${rad},${rad} h-${Math.max(0, w - rad)} z`;
      el("path", { d: p, class: r.hl ? "bar hl" : "bar" }, svg);
      txt(svg, xe + 6, y + 4, spec.fmt(r.value, true), "vlabel");
      const hit = el("rect", { x: 0, y: y - rowH / 2, width: W, height: rowH, class: "hit" }, svg);
      bindTip(hit, spec.fmt(r.value, true), r.tip || r.label);
      void lab;
    });
    host.appendChild(svg);
  }

  // ── multi-series lines ────────────────────────────────────────────────────
  // spec: {series:[{name, color, pts:[[x,y],...]}], xlog, ylog, xfmt, yfmt,
  //        xLabel, yLabel, legend:true}
  function lines(host, spec) {
    const W = spec.W || 760, H = spec.H || 330;
    const m = { l: 74, r: spec.endLabels ? 190 : 24, t: 14, b: 46 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, class: "chart", role: "img",
                            "aria-label": spec.aria || "" }, null);
    const xs = spec.series.flatMap(s => s.pts.map(p => p[0]));
    const ys = spec.series.flatMap(s => s.pts.map(p => p[1]));
    let xlo = Math.min(...xs), xhi = Math.max(...xs), ylo = Math.min(...ys), yhi = Math.max(...ys);
    let xt, yt;
    if (spec.xlog) { xt = logTicks(xlo, xhi); }
    else { xt = spec.xticks || niceTicks(xlo, xhi, 6); }
    if (spec.ylog) {
      ylo = Math.pow(10, Math.floor(Math.log10(ylo))); yhi = Math.pow(10, Math.ceil(Math.log10(yhi)));
      yt = logTicks(ylo, yhi);
    } else {
      ylo = 0; yt = niceTicks(0, yhi * 1.05, 5); yhi = yt[yt.length - 1];
    }
    const sx = scale(xlo, xhi, m.l, W - m.r, !!spec.xlog);
    const sy = scale(ylo, yhi, H - m.b, m.t, !!spec.ylog);
    for (const t of yt) {
      el("line", { x1: m.l, x2: W - m.r, y1: sy(t), y2: sy(t), class: "grid" }, svg);
      txt(svg, m.l - 8, sy(t) + 4, spec.yfmt(t), "tick", "end");
    }
    for (const t of xt) txt(svg, sx(t), H - m.b + 17, spec.xfmt(t), "tick", "middle");
    el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: "baseline" }, svg);
    if (spec.xLabel) txt(svg, (m.l + W - m.r) / 2, H - 8, spec.xLabel, "axis-label", "middle");
    if (spec.yLabel) {
      const t = txt(svg, 0, 0, spec.yLabel, "axis-label", "middle");
      t.setAttribute("transform", `translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)`);
    }
    const ends = [];
    spec.series.forEach(s => {
      const d = s.pts.map((p, i) => (i ? "L" : "M") + sx(p[0]).toFixed(1) + "," + sy(p[1]).toFixed(1)).join(" ");
      el("path", { d, class: "line", style: `stroke:${s.color}` }, svg);
      s.pts.forEach(p => {
        const c = el("circle", { cx: sx(p[0]), cy: sy(p[1]), r: 4, class: "dot", style: `fill:${s.color}` }, svg);
        const hit = el("circle", { cx: sx(p[0]), cy: sy(p[1]), r: 12, class: "hit" }, svg);
        bindTip(hit, spec.yfmt(p[1], true), `${s.name} · ${spec.xfmt(p[0], true)}`);
        void c;
      });
      const last = s.pts[s.pts.length - 1];
      ends.push({ y: sy(last[1]), x: sx(last[0]), name: s.name, color: s.color });
    });
    if (spec.endLabels) {
      // de-collide end labels with leader lines
      ends.sort((a, b) => a.y - b.y);
      let prev = -1e9;
      ends.forEach(e => { e.ly = Math.max(e.y, prev + 15); prev = e.ly; });
      ends.forEach(e => {
        const lx = W - m.r + 14;
        el("line", { x1: e.x + 5, y1: e.y, x2: lx - 3, y2: e.ly, class: "leader" }, svg);
        txt(svg, lx, e.ly + 4, e.name, "endlabel");
      });
    }
    host.appendChild(svg);
    if (spec.legend !== false) {
      const lg = document.createElement("div");
      lg.className = "legend";
      spec.series.forEach(s => {
        const i = document.createElement("span");
        const k = document.createElement("i"); k.style.background = s.color;
        const t = document.createElement("span"); t.textContent = s.name;
        i.append(k, t); lg.appendChild(i);
      });
      host.appendChild(lg);
    }
  }

  window.Charts = { hbar, lines };
})();
