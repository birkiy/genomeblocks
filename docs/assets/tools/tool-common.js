/* ============================================================
   Tools — shared JS utilities.

   Exposes window.GBTools. Every tool page pulls this script and
   then its own inline logic; GBTools provides parsers, set
   algebra, drop-zone wiring, an SRI-pinned CDN loader, and a
   local-backend client stub used by the (future) heatmap tool.

   No dependencies. Runs in every evergreen browser.
   ============================================================ */
(function () {
  'use strict';

  // ── BED / narrowPeak / broadPeak parser ───────────────────
  // Accepts BED3+ (tab-separated: chrom, start, end, ...). Skips
  // blank lines, comments (#), `track` / `browser` header lines,
  // and any row whose first three columns are not a valid interval.
  function parseBED(text) {
    const out = [];
    const lines = text.split(/\r?\n/);
    let bad = 0;
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      if (!line) continue;
      const trimmed = line.trimStart();
      if (!trimmed || trimmed[0] === '#' ||
          trimmed.startsWith('track') ||
          trimmed.startsWith('browser')) continue;
      const parts = line.split(/\t/);
      if (parts.length < 3) { bad++; continue; }
      const chrom = parts[0];
      const start = +parts[1];
      const end   = +parts[2];
      if (!chrom || !Number.isFinite(start) || !Number.isFinite(end) || end <= start) {
        bad++; continue;
      }
      out.push({ chrom: chrom, start: start, end: end });
    }
    return { intervals: out, skipped: bad };
  }

  // Per-set self-merge — sort + collapse overlapping/adjacent
  // intervals. Matches `Loci.make(file).sort().merge()`.
  function mergeSingle(intervals) {
    const sorted = intervals.slice().sort(function (a, b) {
      if (a.chrom !== b.chrom) return a.chrom < b.chrom ? -1 : 1;
      return a.start - b.start;
    });
    let merged = 0;
    let cur = null;
    for (let i = 0; i < sorted.length; i++) {
      const iv = sorted[i];
      if (cur && cur.chrom === iv.chrom && iv.start < cur.end) {
        if (iv.end > cur.end) cur.end = iv.end;
      } else {
        if (cur) merged++;
        cur = { chrom: iv.chrom, start: iv.start, end: iv.end };
      }
    }
    if (cur) merged++;
    return merged;
  }

  // Merged-universe labelling for set-algebra overlap counts.
  // Each input set carries a bit. We sort ALL intervals globally,
  // sweep, merge overlapping runs and OR their source bits. The
  // resulting map (mask → merged-region count) answers every
  // Venn-region count. Mirrors bedtools multiinter.
  function buildRegions(sets) {
    const all = [];
    sets.forEach(function (s) {
      for (let i = 0; i < s.intervals.length; i++) {
        const iv = s.intervals[i];
        all.push({ chrom: iv.chrom, start: iv.start, end: iv.end, src: s.bit });
      }
    });
    all.sort(function (x, y) {
      if (x.chrom !== y.chrom) return x.chrom < y.chrom ? -1 : 1;
      return x.start - y.start;
    });

    const counts = Object.create(null);
    const bpPerMask = Object.create(null);
    let cur = null;
    for (let i = 0; i < all.length; i++) {
      const iv = all[i];
      if (cur && cur.chrom === iv.chrom && iv.start < cur.end) {
        if (iv.end > cur.end) cur.end = iv.end;
        cur.src |= iv.src;
      } else {
        if (cur) {
          counts[cur.src] = (counts[cur.src] || 0) + 1;
          bpPerMask[cur.src] = (bpPerMask[cur.src] || 0) + (cur.end - cur.start);
        }
        cur = { chrom: iv.chrom, start: iv.start, end: iv.end, src: iv.src };
      }
    }
    if (cur) {
      counts[cur.src] = (counts[cur.src] || 0) + 1;
      bpPerMask[cur.src] = (bpPerMask[cur.src] || 0) + (cur.end - cur.start);
    }
    return { counts: counts, bp: bpPerMask };
  }

  // ── Drop-zone wiring ──────────────────────────────────────
  // Attach drag/drop + file-picker handlers to a .gb-tool-drop
  // element. The caller gets one callback: onFile(File).
  function attachDropZone(slot, opts) {
    opts = opts || {};
    const onFile = opts.onFile || function () {};
    const picker = slot.querySelector('.gb-tool-picker');

    if (picker) {
      picker.addEventListener('change', function (e) {
        const f = e.target.files && e.target.files[0];
        if (f) onFile(f);
      });
    }
    ['dragenter', 'dragover'].forEach(function (ev) {
      slot.addEventListener(ev, function (e) {
        e.preventDefault(); e.stopPropagation();
        slot.classList.add('is-dragover');
      });
    });
    ['dragleave', 'dragend', 'drop'].forEach(function (ev) {
      slot.addEventListener(ev, function (e) {
        e.preventDefault(); e.stopPropagation();
        slot.classList.remove('is-dragover');
      });
    });
    slot.addEventListener('drop', function (e) {
      const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      if (f) onFile(f);
    });
  }

  // ── Formatters ────────────────────────────────────────────
  function fmtInt(n) { return (n || 0).toLocaleString('en-US'); }
  function humanBp(n) {
    if (!Number.isFinite(n) || n < 0) return '0 bp';
    if (n >= 1e6) return (n / 1e6).toFixed(2) + ' Mb';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + ' kb';
    return n + ' bp';
  }
  function truncateName(name, max) {
    max = max || 22;
    if (!name) return '';
    if (name.length <= max) return name;
    return name.slice(0, max - 1) + '…';
  }
  function escapeHtml(s) {
    return String(s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#39;');
  }

  // ── CDN loader with SRI ───────────────────────────────────
  // Loads <script>s sequentially (respecting peer-dep order) and
  // resolves when each has executed. Pass { src, integrity }.
  function loadScripts(items) {
    return items.reduce(function (p, item) {
      return p.then(function () { return loadOneScript(item); });
    }, Promise.resolve());
  }

  function loadOneScript(item) {
    return new Promise(function (resolve, reject) {
      const existing = document.querySelector('script[data-gbtools-src="' + item.src + '"]');
      if (existing) {
        if (existing.dataset.loaded === '1') return resolve();
        existing.addEventListener('load', function () { resolve(); }, { once: true });
        existing.addEventListener('error', function () { reject(new Error('Failed to load ' + item.src)); }, { once: true });
        return;
      }
      const s = document.createElement('script');
      s.src = item.src;
      s.async = false;
      s.dataset.gbtoolsSrc = item.src;
      if (item.integrity) {
        s.integrity = item.integrity;
        s.crossOrigin = 'anonymous';
      }
      s.addEventListener('load', function () { s.dataset.loaded = '1'; resolve(); }, { once: true });
      s.addEventListener('error', function () { reject(new Error('Failed to load ' + item.src)); }, { once: true });
      document.head.appendChild(s);
    });
  }

  // ── Backend client (for the future heatmap tool) ──────────
  // Points at a locally-running `genomeblocks serve` process.
  // The docs site stays static/GH-Pages-hostable; tools that
  // need Python compute (Loci.signal on BigWigs) call into this.
  // Unused by the Venn tool today.
  const backend = {
    baseUrl: 'http://localhost:4455',

    setBaseUrl: function (url) { this.baseUrl = String(url || '').replace(/\/$/, ''); },

    probe: function () {
      const url = this.baseUrl + '/healthz';
      return fetch(url, { mode: 'cors' })
        .then(function (r) { return r.ok; })
        .catch(function () { return false; });
    },

    post: function (path, body) {
      return fetch(this.baseUrl + path, {
        method: 'POST',
        mode: 'cors',
        headers: { 'content-type': 'application/json' },
        body: JSON.stringify(body || {}),
      }).then(function (r) {
        if (!r.ok) return r.text().then(function (t) { throw new Error(t || r.statusText); });
        return r.json();
      });
    },
  };

  // ── Expose ────────────────────────────────────────────────
  window.GBTools = {
    parseBED: parseBED,
    mergeSingle: mergeSingle,
    buildRegions: buildRegions,
    attachDropZone: attachDropZone,
    fmtInt: fmtInt,
    humanBp: humanBp,
    truncateName: truncateName,
    escapeHtml: escapeHtml,
    loadScripts: loadScripts,
    backend: backend,
  };
})();
