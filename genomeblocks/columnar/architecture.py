"""Architecture as two tables: vertices = Loci rows, edges = a sorted table.

    vertices   the CRE Loci itself (row i = vertex i), plus vertex columns in ``vp``
    edges      src | tgt | ep columns (w, n, d, ...)    src < tgt, int32

Edges are kept sorted so that the intra-chromosomal (cis) edges of each
chromosome form one contiguous block, in genome order, and every
inter-chromosomal (trans) edge sits in one final block::

    | cis chr1 | cis chr2 | ... | cis chrX | trans |
    ^offsets['chr1']           ^offsets['chrX']  ^trans

Seen as an adjacency matrix this is block-diagonal plus a few off-diagonal
dots — one graph, nothing split apart:

    * ``A.chrom('chr2')``, ``A.cis``, ``A.trans`` are slices (zero-copy views);
      each chromosome can be processed on its own core.
    * neighbours come from one adjacency index over *all* edges, so a CRE's
      trans partners are always there.
    * graph algorithms (components, centrality, ...) use a graph-tool Graph
      built from the arrays on demand (``A.graph()``) — it sees every edge.

Vertex ids are Loci row numbers, so vertex columns, annotations and signal
cubes line up with the Loci by position: no uid -> vertex dictionary.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from .loci import Loci


class Props(dict):
    """dict with attribute access: ``A.ep.w`` == ``A.ep['w']``."""

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(k) from None

    def __setattr__(self, k, v):
        self[k] = v


def _pl_model(x, C, alpha):
    x = np.asarray(x, float)
    out = np.full(x.shape, np.inf)
    pos = x > 0
    out[pos] = C * (x[pos] ** (-alpha))
    return out


def _pl_expect(x, y):
    """Power-law fit of y ~ C * x^-alpha (same procedure as the classic normalize)."""
    from scipy.optimize import curve_fit
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if not m.any():
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}
    try:
        b, a = np.polyfit(np.log10(x[m]), np.log10(y[m]), 1)
        popt, _ = curve_fit(_pl_model, x[m], y[m], p0=[10.0 ** a, -b])
        return _pl_model(x, *popt), {"alpha": float(popt[1]), "C": float(popt[0])}
    except (RuntimeError, np.linalg.LinAlgError, ValueError):
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}


def _cross(k1, i1, k2, j2, n):
    """Per-group cross product: all (i, j) with i from group k in (k1, i1) and
    j from the same group in (k2, j2). Both inputs sorted by group."""
    c2 = np.bincount(k2, minlength=n)
    off2 = np.concatenate([[0], np.cumsum(c2)[:-1]])
    rep = c2[k1]
    src = np.repeat(i1, rep)
    first = np.repeat(off2[k1], rep)
    within = np.arange(rep.sum()) - np.repeat(np.cumsum(rep) - rep, rep)
    return src, j2[first + within]


class Architecture:
    """CRE interaction graph stored as a vertex table (the Loci) + an edge table."""

    def __init__(self, loci: Loci, src=(), tgt=(), *, ep=None, vp=None, name: str = "Architecture",
                 _canonical: bool = False):
        if not loci.is_sorted:
            raise ValueError("Architecture needs genome-sorted Loci (Loci.make sorts by default; "
                             "or call .sort()).")
        self.loci = loci
        self.name = name
        src = np.asarray(src, np.int32)
        tgt = np.asarray(tgt, np.int32)
        ep = {k: np.asarray(v, float) for k, v in (ep or {}).items()}
        if not _canonical:
            src, tgt, order = self._canonical_order(src, tgt)
            ep = {k: v[order] for k, v in ep.items()}
        self.src, self.tgt = src, tgt
        self.ep = Props(ep)
        self.vp = Props(vp or {})
        self._reset()

    def _reset(self):
        self._blocks = None
        self._csr = None
        self._gt = None
        self._keys = None

    def _canonical_order(self, src, tgt):
        lo, hi = np.minimum(src, tgt), np.maximum(src, tgt)
        c = self.loci.codes
        rank = self.loci.genome.rank
        cis = c[lo] == c[hi]
        block = np.where(cis, rank[c[lo]], len(rank) + rank[c[lo]])
        order = np.lexsort((hi, lo, block))
        return lo[order], hi[order], order

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, loci: Loci, bedpe: str, *, name: str = "Skeleton", r: int = 2500, dmax: float = 1e9,
             verbose: bool = True) -> "Architecture":
        """Edges between every pair of CREs within ``r`` of the two anchor
        midpoints of a loop (same rule as the classic ``Architecture.make``)."""
        import pandas as pd
        from . import _intervals as K
        if not loci.is_sorted:
            loci = loci.sort()
        df = pd.read_csv(bedpe, sep="\t", header=None, comment="#", usecols=range(6),
                         dtype={0: str, 3: str})
        m1 = (df[1].to_numpy(np.int64) + df[2].to_numpy(np.int64)) // 2
        m2 = (df[4].to_numpy(np.int64) + df[5].to_numpy(np.int64)) // 2
        keep = np.abs(m1 - m2) <= dmax
        g = loci.genome
        c1, c2 = g.encode(df[0].to_numpy()[keep]), g.encode(df[3].to_numpy()[keep])
        m1, m2 = m1[keep], m2[keep]
        k1, i1 = K.overlap_pairs(c1, m1 - r, m1 + r, loci.codes, loci.starts, loci.ends)
        k2, j2 = K.overlap_pairs(c2, m2 - r, m2 + r, loci.codes, loci.starts, loci.ends)
        o2 = np.argsort(k2, kind="stable")
        src, tgt = _cross(k1, i1, k2[o2], j2[o2], len(m1))
        ok = src != tgt
        lo, hi = np.minimum(src[ok], tgt[ok]), np.maximum(src[ok], tgt[ok])
        key = np.unique(lo.astype(np.int64) * len(loci) + hi)
        A = cls(loci, key // len(loci), key % len(loci), name=name)
        n = len(A.src)
        A.ep.update(w=np.zeros(n), n=np.zeros(n), d=np.zeros(n))
        if verbose:
            mapped = len(np.intersect1d(np.unique(k1), np.unique(k2)))
            print(f"[INFO] {len(df)} loops | {mapped} mapped ({100 * mapped / max(len(df), 1):.1f}%) | "
                  f"loci={A.n_loci}, links={A.n_links} ({A.n_trans} trans)")
        return A

    # ── sizes ─────────────────────────────────────────────────────────────
    @property
    def n_links(self) -> int:
        return len(self.src)

    @property
    def degree(self) -> np.ndarray:
        n = len(self.loci)
        return np.bincount(self.src, minlength=n) + np.bincount(self.tgt, minlength=n)

    @property
    def n_loci(self) -> int:
        """Vertices with at least one link (rows of the Loci with degree > 0)."""
        return int((self.degree > 0).sum())

    @property
    def is_cis(self) -> np.ndarray:
        c = self.loci.codes
        return c[self.src] == c[self.tgt]

    @property
    def n_trans(self) -> int:
        a, b = self.blocks["trans"]
        return b - a

    def __len__(self):
        return self.n_loci

    # ── blocks & views ────────────────────────────────────────────────────
    @property
    def blocks(self) -> dict:
        """{chrom: (lo, hi)} edge ranges for each cis block, plus 'trans'."""
        if self._blocks is None:
            c, names = self.loci.codes, self.loci.genome.names
            cs, ct = c[self.src], c[self.tgt]
            trans = np.flatnonzero(cs != ct)
            t0 = int(trans[0]) if len(trans) else len(cs)
            cut = np.flatnonzero(np.diff(cs[:t0])) + 1
            lo = np.concatenate([[0], cut]) if t0 else np.zeros(0, int)
            hi = np.concatenate([cut, [t0]]) if t0 else np.zeros(0, int)
            b = {names[cs[a]]: (int(a), int(z)) for a, z in zip(lo, hi)}
            b["trans"] = (t0, len(cs))
            self._blocks = b
        return self._blocks

    def _view(self, sel, name) -> "Architecture":
        """Edges ``sel`` (slice → zero-copy) sharing the same Loci and vertex columns."""
        v = Architecture.__new__(Architecture)
        v.loci, v.name = self.loci, name
        v.src, v.tgt = self.src[sel], self.tgt[sel]
        v.ep = Props({k: a[sel] for k, a in self.ep.items()})
        v.vp = self.vp                               # shared: same rows, same columns
        v._reset()
        return v

    def chrom(self, chrom: str) -> "Architecture":
        """Cis edges of one chromosome (a view: no copy)."""
        a, b = self.blocks.get(chrom, (0, 0))
        return self._view(slice(a, b), f"{self.name}:{chrom}")

    @property
    def cis(self) -> "Architecture":
        return self._view(slice(0, self.blocks["trans"][0]), f"{self.name}:cis")

    @property
    def trans(self) -> "Architecture":
        a, b = self.blocks["trans"]
        return self._view(slice(a, b), f"{self.name}:trans")

    def chroms(self):
        """Iterate ``(chrom, view)`` over the cis blocks (e.g. to fan out per chromosome)."""
        for c in self.blocks:
            if c != "trans":
                yield c, self.chrom(c)

    def subgraph(self, rows=None, *, mask=None, vp: Optional[str] = None, values=None,
                 name: Optional[str] = None) -> "Architecture":
        """Edges whose two ends are both selected (rows, bool mask, or vp in values)."""
        n = len(self.loci)
        if mask is None:
            mask = np.zeros(n, bool)
            if rows is not None:
                mask[np.asarray(rows)] = True
            elif vp is not None:
                vals = values if isinstance(values, (list, tuple, set, np.ndarray)) else [values]
                mask = np.isin(self.vp[vp], list(vals))
        sel = np.flatnonzero(mask[self.src] & mask[self.tgt])
        return self._view(sel, name or f"{self.name}_sub")

    def region(self, region, start=None, end=None, *, both: bool = True) -> "Architecture":
        """Edges inside a region ('chr1:1-5,000,000' or chrom, start, end).
        ``both=False`` also keeps edges with one end inside (trans partners too)."""
        mask = np.zeros(len(self.loci), bool)
        mask[self.loci.overlap_rows(region, start, end)] = True
        hit = (mask[self.src] & mask[self.tgt]) if both else (mask[self.src] | mask[self.tgt])
        label = region if start is None else f"{region}:{start}-{end}"
        return self._view(np.flatnonzero(hit), f"{self.name}:{label}")

    # ── adjacency (all edges, both directions) ───────────────────────────
    def _adj(self):
        if self._csr is None:
            n = len(self.loci)
            a = np.concatenate([self.src, self.tgt])
            b = np.concatenate([self.tgt, self.src])
            e = np.concatenate([np.arange(self.n_links), np.arange(self.n_links)])
            o = np.argsort(a, kind="stable")
            indptr = np.concatenate([[0], np.cumsum(np.bincount(a, minlength=n))])
            self._csr = (indptr, b[o], e[o])
        return self._csr

    def _row(self, x) -> int:
        if isinstance(x, str):
            return self.loci.uids[x]
        if hasattr(x, "row"):
            return x.row
        return int(x)

    def neighbor_rows(self, x):
        """(partner rows, edge ids) of one CRE — cis and trans — as arrays."""
        i = self._row(x)
        indptr, nb, eid = self._adj()
        return nb[indptr[i]:indptr[i + 1]], eid[indptr[i]:indptr[i + 1]]

    def neighbors(self, x):
        """All partners of one CRE (row / uid / Locus) — cis and trans — as a DataFrame."""
        import pandas as pd
        i = self._row(x)
        j, e = self.neighbor_rows(i)
        L = self.loci
        d = {"row": j, "uid": L.uid[j], "chrom": L.chroms[j], "cis": L.codes[j] == L.codes[i]}
        d.update({k: v[e] for k, v in self.ep.items()})
        return pd.DataFrame(d)

    def __getitem__(self, key):
        """Classic API: ``A[uid]`` -> {neighbour uid: w}; ``A[uid1, uid2]`` -> edge props."""
        if isinstance(key, tuple) and len(key) == 2:
            i, j = sorted((self._row(key[0]), self._row(key[1])))
            e = self._edge_id(i, j)
            return None if e < 0 else {k: float(v[e]) for k, v in self.ep.items()}
        j, e = self.neighbor_rows(key)
        w = self.ep["w"][e] if "w" in self.ep else np.zeros(len(e))
        return dict(zip(self.loci.uid[j].tolist(), w.tolist()))

    def _edge_id(self, i, j) -> int:
        if self._keys is None:
            k = self.src.astype(np.int64) * len(self.loci) + self.tgt
            o = np.argsort(k)
            self._keys = (k[o], o)
        k, o = self._keys
        q = np.int64(i) * len(self.loci) + j
        p = np.searchsorted(k, q)
        return int(o[p]) if p < len(k) and k[p] == q else -1

    def __contains__(self, uid) -> bool:
        i = self.loci.uids.get(uid)
        return i is not None and self.degree[i] > 0

    # ── Hi-C weights ──────────────────────────────────────────────────────
    def add_mcool(self, mcool: str, *, resolution: Optional[int] = None, name: str = "w",
                  verbose: bool = True) -> "Architecture":
        """Edge weight = Hi-C count of the (bin, bin) pixel holding the two CREs,
        shared equally among the edges that fall in the same pixel.

        Same rule as the classic ``add_mcool``. Pixels are read one chromosome
        block at a time (rows of the lower bin), so memory stays bounded and
        trans pixels are found in the same pass.
        """
        import cooler
        uri = f"{mcool}::resolutions/{resolution}" if resolution else mcool
        clr = cooler.Cooler(uri)
        if clr.binsize is None:
            raise NotImplementedError("variable-size bins")
        L = self.loci
        nb = int(clr.info["nbins"])
        names = L.genome.names
        off = np.full(len(names), -1, np.int64)
        last = np.full(len(names), -1, np.int64)
        for c in clr.chromnames:
            if c in L.genome.code:
                lo, hi = clr.extent(c)
                off[L.genome.code[c]], last[L.genome.code[c]] = lo, hi - 1
        vb = np.where(off[L.codes] >= 0, np.minimum(off[L.codes] + L.starts // clr.binsize,
                                                    last[L.codes]), -1)
        b1, b2 = vb[self.src], vb[self.tgt]
        ok = (b1 >= 0) & (b2 >= 0)
        key = np.minimum(b1, b2) * nb + np.maximum(b1, b2)
        ukey, inv, n_share = np.unique(key[ok], return_inverse=True, return_counts=True)
        count = np.zeros(len(ukey))
        ubin1 = ukey // nb
        with clr.open("r") as h5:
            bin1_off = h5["indexes/bin1_offset"]
            px1, px2, pxc = h5["pixels/bin1_id"], h5["pixels/bin2_id"], h5["pixels/count"]
            for c in clr.chromnames:                     # one chromosome block at a time
                lo, hi = clr.extent(c)
                a, z = np.searchsorted(ubin1, [lo, hi])
                if a == z:
                    continue
                p0, p1 = int(bin1_off[ubin1[a]]), int(bin1_off[ubin1[z - 1] + 1])
                pk = px1[p0:p1].astype(np.int64) * nb + px2[p0:p1]
                hit = np.searchsorted(pk, ukey[a:z])
                hit = np.minimum(hit, len(pk) - 1)
                found = pk[hit] == ukey[a:z]
                count[a:z][found] = pxc[p0:p1][hit[found]]
        w = np.zeros(self.n_links)
        w[ok] = (count / n_share)[inv]
        self.ep[name] = w
        self._gt = None
        if verbose:
            print(f"[INFO] Set distributed weights for {int((w > 0).sum())}/{self.n_links} edges "
                  f"from cooler. [{name}]")
        return self

    # ── distance-decay normalisation ─────────────────────────────────────
    def normalize(self, *, source: str = "w", name: str = "n", verbose: bool = True) -> "Architecture":
        """Observed / expected. Cis: power-law fit on distance (classic rule).
        Trans: no distance, so expected = mean trans weight. ``ep.d`` is inf for trans."""
        if source not in self.ep:
            raise ValueError(f"Source edge property '{source}' not found in the graph.")
        cen = self.loci.centers
        cis = self.is_cis
        d = np.full(self.n_links, np.inf)
        d[cis] = np.abs(cen[self.src[cis]] - cen[self.tgt[cis]])
        w = self.ep[source]
        exp = np.empty(self.n_links)
        exp[cis], fit = _pl_expect(d[cis], w[cis])
        tw = w[~cis]
        exp[~cis] = tw.mean() if len(tw) and tw.mean() > 0 else np.nan
        self.ep.d = d
        self.ep[name] = w / np.maximum(exp, 1e-12)
        self.ep[name][~np.isfinite(self.ep[name])] = 0.0
        self.fit = fit
        self._gt = None
        if verbose:
            print(f"[INFO] Power-law fit: alpha={fit['alpha']:.3f}, C={fit['C']:.3e} on "
                  f"{int(cis.sum())} cis edges; {int((~cis).sum())} trans edges use the mean trans "
                  f"weight → ep.{name}")
        return self

    def prune(self, *, dist_prop: str = "d", verbose: bool = True) -> "Architecture":
        """Drop zero-distance (co-located) cis edges. Trans edges (d = inf) stay."""
        if dist_prop not in self.ep:
            raise ValueError(f"Edge property '{dist_prop}' not found — run normalize() first.")
        keep = self.ep[dist_prop] > 0
        n0 = int((~keep).sum())
        if n0:
            self.src, self.tgt = self.src[keep], self.tgt[keep]
            for k in list(self.ep):
                self.ep[k] = self.ep[k][keep]
            self._reset()
        if verbose:
            print(f"[INFO] prune: removed {n0} zero-distance edges → {self.n_links} edges.")
        return self

    # ── gene annotation ───────────────────────────────────────────────────
    def annotate(self, genes, *, key: str = "n", name: str = "gene", verbose: bool = True) -> "Architecture":
        """Region label per CRE + gene assignment (classic two-stage rule):
        promoter CREs get their nearest TSS gene; every other CRE gets the gene
        of its highest-``key`` promoter neighbour (cis or trans)."""
        from .genes import LABELS
        if key not in self.ep:
            raise ValueError(f"Edge property '{key}' not found — run add_mcool/normalize first.")
        L = self.loci
        lab = genes.labels(L)
        is_prom = lab == 5
        gene = np.full(len(L), "", dtype=object)
        pr = np.flatnonzero(is_prom)
        gene[pr] = genes.nearest_tss(L.take(pr))[0]
        src, tgt, w = self.src.astype(np.int64), self.tgt.astype(np.int64), self.ep[key]
        eligible = is_prom & (gene != "")
        c1 = ~is_prom[src] & eligible[tgt]
        c2 = ~is_prom[tgt] & eligible[src]
        idx = np.concatenate([src[c1], tgt[c2]])
        pg = np.concatenate([gene[tgt[c1]], gene[src[c2]]])
        we = np.concatenate([w[c1], w[c2]])
        n_assigned = 0
        if idx.size:
            order = np.lexsort((-we, idx))
            _, first = np.unique(idx[order], return_index=True)
            sel = order[first]
            gene[idx[sel]] = pg[sel]
            n_assigned = int((pg[sel] != "").sum())
        self.vp.annot = LABELS[lab]
        self.vp[name] = gene
        if verbose:
            on = self.degree > 0
            print(f"[INFO] Annotated {int(on.sum())} loci: {int((is_prom & on).sum())} promoter CREs | "
                  f"{n_assigned}/{int((~is_prom & on).sum())} non-promoter CREs assigned to a "
                  f"top-'{key}' promoter gene → vp.{name}.")
        return self

    # ── hubs ──────────────────────────────────────────────────────────────
    def strength(self, key: str = "n", name: str = "strength", *, verbose: bool = True) -> "Architecture":
        """Sum of incident ``ep[key]`` per vertex (two bincounts)."""
        if key not in self.ep:
            raise ValueError(f"Edge property '{key}' not found.")
        n = len(self.loci)
        w = self.ep[key]
        self.vp[name] = np.bincount(self.src, w, n) + np.bincount(self.tgt, w, n)
        if verbose:
            print(f"[INFO] Summed ep.{key} → vp.{name} (node strength).")
        return self

    def elbow(self, key: str, *, verbose: bool = True):
        """Slope-1 knee on the sorted ``vp[key]`` curve (classic rule).
        Returns (cutoff, uids sorted by value, descending)."""
        from scipy.ndimage import uniform_filter1d
        v = self.vp[key]
        rows = np.flatnonzero(v > 0)
        rows = rows[np.argsort(-v[rows], kind="stable")]
        uids = self.loci.uid[rows].tolist()
        y = v[rows]
        if len(y) < 3:
            return len(y), uids
        ys = y[::-1]
        m = len(ys)
        xn = np.arange(m) / (m - 1)
        yn = (ys - ys[0]) / (ys[-1] - ys[0])
        yn_s = uniform_filter1d(yn, size=max(11, m // 200), mode="nearest")
        crossed = np.flatnonzero(np.gradient(yn_s, xn) >= 1.0)
        i = int(crossed[0]) if len(crossed) else m
        cutoff = m - i
        if verbose:
            print(f"[INFO] Slope-1 on vp.{key} (value≈{(ys[i] if i < m else ys[-1]):.3g} at cutoff): "
                  f"cutoff at {cutoff}/{m} ({100 * cutoff / m:.1f}%)")
        return cutoff, uids

    def prime_hubs(self, key: str = "n", gene: str = "gene", *, verbose: bool = True) -> dict:
        for req in ("annot", gene):
            if req not in self.vp:
                raise ValueError(f"vp.{req} missing — run annotate() first.")
        sname = f"strength_{key.replace('n_', '')}" if key != "n" else "strength"
        if sname not in self.vp:
            self.strength(key=key, name=sname, verbose=verbose)
        cutoff, uids = self.elbow(sname, verbose=verbose)
        hub_uids = uids[:cutoff]
        rows = np.array([self.loci.uids[u] for u in hub_uids], np.int64)
        prom = np.char.find(self.vp.annot[rows].astype(str), "Promoter") >= 0 if len(rows) else \
            np.zeros(0, bool)
        g = self.vp[gene][rows]
        p_genes = set(g[prom][g[prom] != ""])
        e_genes = set(g[~prom][g[~prom] != ""])
        if verbose:
            print(f"[INFO] Prime hubs: {len(hub_uids)} hubs ({int(prom.sum())} promoters, "
                  f"{int((~prom).sum())} enhancers)")
            print(f"[INFO] Prime genes: {len(p_genes | e_genes)} = {len(p_genes)} promoter + "
                  f"{len(e_genes)} enhancer (overlap: {len(p_genes & e_genes)})")
        return {"prime_genes": p_genes | e_genes, "promoter_genes": p_genes, "enhancer_genes": e_genes,
                "hub_uids": hub_uids, "cutoff": cutoff,
                "promoter_uids": [u for u, p in zip(hub_uids, prom) if p],
                "enhancer_uids": [u for u, p in zip(hub_uids, prom) if not p]}

    # ── graph algorithms: graph-tool, built from the arrays when asked ───
    def graph(self, vprops=()):
        """A graph-tool Graph over all edges (cis + trans); vertex i = Loci row i.

        Built once and cached; edge columns are re-synced on every call, so
        it always reflects the current weights. ``vprops`` copies numeric
        vertex columns too.
        """
        import graph_tool.all as gt
        if self._gt is None:
            g = gt.Graph(directed=False)
            g.add_vertex(len(self.loci))
            g.add_edge_list(np.column_stack([self.src, self.tgt]))
            self._gt = g
        g = self._gt
        for k, v in self.ep.items():
            if k not in g.ep:
                g.ep[k] = g.new_ep("double")
            g.ep[k].a[:] = v
        for k in vprops:
            g.vp[k] = g.new_vp("double", vals=np.asarray(self.vp[k], float))
        return g

    def components(self, name: str = "component") -> np.ndarray:
        """Connected component per vertex (isolated rows get -1)."""
        import graph_tool.all as gt
        comp, _ = gt.label_components(self.graph())
        c = np.asarray(comp.a, np.int64).copy()
        c[self.degree == 0] = -1
        self.vp[name] = c
        return c

    # ── tables out ────────────────────────────────────────────────────────
    def edges_frame(self):
        import pandas as pd
        L = self.loci
        d = {"src": self.src, "tgt": self.tgt, "uid1": L.uid[self.src], "uid2": L.uid[self.tgt],
             "chrom1": L.chroms[self.src], "chrom2": L.chroms[self.tgt], "cis": self.is_cis}
        d.update(self.ep)
        return pd.DataFrame(d)

    def to_frame(self):
        """Vertex table for vertices with links (classic layout: uid + vertex props)."""
        import pandas as pd
        rows = np.flatnonzero(self.degree > 0)
        d = {"uid": self.loci.uid[rows]}
        d.update({k: np.asarray(v)[rows] for k, v in self.vp.items()})
        return pd.DataFrame(d)

    def block_counts(self):
        """chrom x chrom edge counts (cis on the diagonal) as a DataFrame."""
        import pandas as pd
        L = self.loci
        c1, c2 = L.codes[self.src], L.codes[self.tgt]
        used = np.unique(np.concatenate([c1, c2]))
        used = used[np.argsort(L.genome.rank[used])]
        pos = np.full(len(L.genome), -1)
        pos[used] = np.arange(len(used))
        M = np.zeros((len(used), len(used)), np.int64)
        np.add.at(M, (pos[c1], pos[c2]), 1)
        M = M + M.T - np.diag(np.diag(M))
        names = [L.genome.names[c] for c in used]
        return pd.DataFrame(M, index=names, columns=names)

    def to_legacy(self):
        """The classic graph-tool ``genomeblocks.Architecture`` (e.g. for its draw helpers)."""
        import contextlib
        import io
        from ..architecture import Architecture as Classic
        with contextlib.redirect_stdout(io.StringIO()):
            G = Classic(name=self.name)
        rows = np.flatnonzero(self.degree > 0)
        pos = np.full(len(self.loci), -1, np.int64)
        pos[rows] = np.arange(len(rows))
        G.add_vertex(len(rows))
        uids = self.loci.uid[rows]
        for v, u in zip(G.vertices(), uids.tolist()):
            G.vp.uid[v] = u
        G.index = {u: G.vertex(k) for k, u in enumerate(uids.tolist())}
        props = [G.ep[k] if k in G.ep else G.new_ep("double") for k in self.ep]
        for k, p in zip(self.ep, props):
            G.ep[k] = p
        G.add_edge_list(np.column_stack([pos[self.src], pos[self.tgt]] + [self.ep[k] for k in self.ep]),
                        eprops=props)
        for k, v in self.vp.items():
            v = np.asarray(v)[rows]
            if v.dtype == object:
                G.vp[k] = G.new_vp("string")
                for vert, s in zip(G.vertices(), v.tolist()):
                    G.vp[k][vert] = s
            else:
                G.vp[k] = G.new_vp("double", vals=v.astype(float))
        return G

    # ── set operations (edge keys) ────────────────────────────────────────
    def _keyset(self):
        return self.src.astype(np.int64) * len(self.loci) + self.tgt

    def __or__(self, other: "Architecture") -> "Architecture":
        if other.loci is not self.loci:
            raise ValueError("union needs both graphs on the same Loci")
        k = np.concatenate([self._keyset(), other._keyset()])
        _, first = np.unique(k, return_index=True)
        ep = {e: np.concatenate([self.ep.get(e, np.zeros(self.n_links)),
                                 other.ep.get(e, np.zeros(other.n_links))])[first]
              for e in set(self.ep) | set(other.ep)}
        src = np.concatenate([self.src, other.src])[first]
        tgt = np.concatenate([self.tgt, other.tgt])[first]
        return Architecture(self.loci, src, tgt, ep=ep, name=f"{self.name}|{other.name}")

    def __and__(self, other: "Architecture") -> "Architecture":
        if other.loci is not self.loci:
            raise ValueError("intersection needs both graphs on the same Loci")
        keep = np.isin(self._keyset(), other._keyset())
        v = self._view(np.flatnonzero(keep), f"{self.name}&{other.name}")
        v.vp = Props()
        return v

    def copy(self) -> "Architecture":
        v = self._view(slice(None), self.name)
        v.src, v.tgt = v.src.copy(), v.tgt.copy()
        v.ep = Props({k: a.copy() for k, a in v.ep.items()})
        v.vp = Props({k: np.array(a, copy=True) for k, a in self.vp.items()})
        return v

    # ── persistence: parquet tables + a small json ───────────────────────
    def save(self, path: str):
        """``path/`` gets vertices.parquet (Loci + vp), edges.parquet (src, tgt, ep), meta.json."""
        import json
        import os
        import pyarrow as pa
        import pyarrow.parquet as pq
        os.makedirs(path, exist_ok=True)
        vt = self.loci.to_arrow()
        for k, v in self.vp.items():
            vt = vt.append_column(f"vp:{k}", pa.array(v))
        pq.write_table(vt, os.path.join(path, "vertices.parquet"))
        et = pa.table({"src": self.src, "tgt": self.tgt, **{f"ep:{k}": v for k, v in self.ep.items()}})
        pq.write_table(et, os.path.join(path, "edges.parquet"))
        with open(os.path.join(path, "meta.json"), "w") as f:
            json.dump({"name": self.name, "format": "genomeblocks.columnar.Architecture/1"}, f)

    @classmethod
    def load(cls, path: str, *, genome=None) -> "Architecture":
        import json
        import os
        import pyarrow.parquet as pq
        vt = pq.read_table(os.path.join(path, "vertices.parquet"))
        vp_cols = [c for c in vt.column_names if c.startswith("vp:")]
        tmp = os.path.join(path, "vertices.parquet")
        L = Loci.load(tmp, genome=genome)
        for c in vp_cols:
            L.cols.pop(c, None)
        L._sorted = True
        vp = {c[3:]: vt.column(c).to_numpy(zero_copy_only=False) for c in vp_cols}
        et = pq.read_table(os.path.join(path, "edges.parquet"))
        ep = {c[3:]: et.column(c).to_numpy() for c in et.column_names if c.startswith("ep:")}
        meta = json.load(open(os.path.join(path, "meta.json")))
        return cls(L, et.column("src").to_numpy(), et.column("tgt").to_numpy(), ep=ep, vp=vp,
                   name=meta.get("name", "Architecture"), _canonical=True)

    # ── display ───────────────────────────────────────────────────────────
    def __repr__(self):
        return (f"Architecture(name='{self.name}', loci={self.n_loci:,}, links={self.n_links:,} "
                f"[{self.n_links - self.n_trans:,} cis · {self.n_trans:,} trans], "
                f"edge_props=[{', '.join(self.ep)}], vertex_props=[{', '.join(self.vp)}])")

    def _repr_html_(self):
        from ._display import table_html
        b = self.blocks
        cis = [(c, hi - lo) for c, (lo, hi) in b.items() if c != "trans"]
        rows = [[c, f"{lo:,}–{hi:,}", f"{hi - lo:,}"] for c, (lo, hi) in b.items() if c != "trans"]
        show = rows[:4] + [["⋮", "", ""]] + rows[-2:] if len(rows) > 6 else rows
        t0, t1 = b["trans"]
        show.append(["trans", f"{t0:,}–{t1:,}", f"{t1 - t0:,}"])
        mb = (self.src.nbytes + self.tgt.nbytes + sum(v.nbytes for v in self.ep.values())) / 1e6
        title = (f"Architecture '{self.name}' · {self.n_loci:,} loci with links · {self.n_links:,} edges "
                 f"({len(cis)} cis blocks + trans) · {mb:.1f} MB")
        note = (f"edge columns: {', '.join(self.ep) or '—'}   ·   vertex columns: "
                f"{', '.join(self.vp) or '—'}   ·   vertex i = Loci row i")
        return table_html(title, ["block", "edge rows", "edges"], show, note=note)
