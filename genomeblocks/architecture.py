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

Graph algorithms run through the graph backend (graph-tool by default;
igraph, networkx or scipy on request), and the graph exports to each of
them, to a scipy sparse matrix, to AnnData and to plain tables.
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
    if m.sum() < 3 or len(np.unique(x[m])) < 2:        # a fit needs a few distances
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}
    try:
        b, a = np.polyfit(np.log10(x[m]), np.log10(y[m]), 1)
        popt, _ = curve_fit(_pl_model, x[m], y[m], p0=[10.0 ** a, -b])
        return _pl_model(x, *popt), {"alpha": float(popt[1]), "C": float(popt[0])}
    except (RuntimeError, TypeError, np.linalg.LinAlgError, ValueError):
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

    def __getstate__(self):
        d = self.__dict__.copy()
        d.update(_blocks=None, _csr=None, _gt=None, _keys=None)
        return d

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
    def make(cls, loci, bedpe, *, name: str = "Skeleton", r: int = 2500, dmax: float = 1e9,
             trans: bool = True, verbose: bool = True, backend: Optional[str] = None) -> "Architecture":
        """Edges between every pair of CREs within ``r`` bp of the two anchor
        midpoints of a loop.

        ``bedpe`` is a BEDPE path, a :class:`~genomeblocks.Pairs` or a frame
        with BEDPE columns. Cis loops longer than ``dmax`` are skipped; trans
        loops are kept unless ``trans=False``. ``backend`` picks the interval
        engine for the anchor-to-CRE mapping. Vertex ``i`` is ``loci`` row ``i``
        (``loci`` must be genome-sorted, as ``Loci.make`` returns it).
        """
        from .backends.intervals import overlap_pairs
        from .bedpe import as_pairs
        from .interop import as_loci
        loci = as_loci(loci)
        if not loci.is_sorted:
            print("[WARN] loci were not in genome order; the Architecture uses a sorted copy "
                  "(A.loci), so align vertex columns to A.loci rather than to the input.")
            loci = loci.sort()
        P = as_pairs(bedpe, genome=loci.genome)
        P = P.take(loci._check(P.a).codes >= 0) if len(P) else P
        a, b = loci._check(P.a), loci._check(P.b)
        m1, m2 = a.centers, b.centers
        cis = a.codes == b.codes
        keep = np.where(cis, np.abs(m1 - m2) <= dmax, trans)
        m1, m2 = m1[keep], m2[keep]
        A1 = Loci(a.codes[keep], np.maximum(m1 - r, 0), m1 + r, genome=loci.genome)
        A2 = Loci(b.codes[keep], np.maximum(m2 - r, 0), m2 + r, genome=loci.genome)
        k1, i1 = overlap_pairs(A1, loci, backend=backend)
        k2, j2 = overlap_pairs(A2, loci, backend=backend)
        o1 = np.argsort(k1, kind="stable")
        o2 = np.argsort(k2, kind="stable")
        src, tgt = _cross(k1[o1], i1[o1], k2[o2], j2[o2], len(A1))
        ok = src != tgt
        lo, hi = np.minimum(src[ok], tgt[ok]), np.maximum(src[ok], tgt[ok])
        key = np.unique(lo.astype(np.int64) * len(loci) + hi)
        A = cls(loci, key // len(loci), key % len(loci), name=name)
        A.ep["w"] = np.zeros(len(A.src))          # add_mcool fills w; normalize adds n and d
        if verbose:
            mapped = len(np.intersect1d(np.unique(k1), np.unique(k2)))
            print(f"[INFO] {len(P)} loops | {mapped} mapped ({100 * mapped / max(len(P), 1):.1f}%) | "
                  f"loci={A.n_loci}, links={A.n_links} ({A.n_trans} trans)")
        return A

    @classmethod
    def from_edges(cls, loci, src, tgt, *, name: str = "Architecture", vp=None, **ep) -> "Architecture":
        """From vertex rows: ``src[k]``–``tgt[k]`` is edge k; keyword arrays become
        edge columns. Duplicate and reversed pairs are kept as given."""
        from .interop import as_loci
        loci = as_loci(loci)
        return cls(loci, src, tgt, ep=ep, vp=vp, name=name)

    @classmethod
    def from_frame(cls, loci, df, *, src="src", tgt="tgt", name: str = "Architecture") -> "Architecture":
        """From an edge table (pandas / polars / arrow) with vertex rows in
        ``src`` / ``tgt`` or uids in ``uid1`` / ``uid2``; other numeric
        columns become edge columns (what :meth:`edges_frame` writes)."""
        from .interop import _as_pandas, as_loci
        loci = as_loci(loci)
        pdf = _as_pandas(df)
        if src in pdf and tgt in pdf:
            s, t = pdf[src].to_numpy(np.int64), pdf[tgt].to_numpy(np.int64)
        elif "uid1" in pdf and "uid2" in pdf:
            u = loci.uids
            bad = [x for x in pdf["uid1"].tolist() + pdf["uid2"].tolist() if x not in u]
            if bad:
                raise ValueError(f"{len(bad)} uid(s) in the edge table are not in the loci, e.g. {bad[:3]}")
            s = np.array([u[x] for x in pdf["uid1"]], np.int64)
            t = np.array([u[x] for x in pdf["uid2"]], np.int64)
        else:
            raise ValueError(f"an edge table needs vertex rows in {src!r} / {tgt!r} or uids in "
                             f"'uid1' / 'uid2' (what edges_frame() writes); got columns {list(pdf.columns)}")
        skip = {src, tgt, "uid1", "uid2", "chrom1", "chrom2", "cis"}
        ep = {k: pdf[k].to_numpy(np.float64) for k in pdf.columns
              if k not in skip and np.issubdtype(pdf[k].dtype, np.number)}
        return cls(loci, s, t, ep=ep, name=name)

    @classmethod
    def from_scipy(cls, loci, M, *, weight: str = "w", name: str = "Architecture") -> "Architecture":
        """From a square sparse matrix over the loci rows (upper triangle read)."""
        from scipy.sparse import triu
        from .interop import as_loci
        loci = as_loci(loci)
        U = triu(M, k=1).tocoo()
        return cls(loci, U.row, U.col, ep={weight: U.data.astype(float)}, name=name)

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
        """``A[uid]`` -> {neighbour uid: w}; ``A[uid1, uid2]`` -> that edge's columns."""
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

    # ── CREs near a position (gene support) ──────────────────────────────
    def near_rows(self, key, start=None, end=None, *, r: int = 0, mode: str = "overlap",
                  linked: bool = True) -> np.ndarray:
        """Row numbers of the CREs near a window.

        ``key`` is a Locus, a region string, or a chromosome with ``start`` /
        ``end``; ``r`` widens the window on each side (a 1-bp TSS gives the
        window TSS ± r). ``mode='overlap'`` keeps CREs that intersect the
        window, ``'center'`` only those whose midpoint is in it. ``linked=True``
        keeps CREs with at least one edge."""
        if mode not in ("overlap", "center"):
            raise ValueError(f"mode must be 'overlap' or 'center', got {mode!r}")
        if hasattr(key, "chrom"):
            key, start, end = key.chrom, key.start, key.end
        elif start is None:
            from .locus import parse_region
            key, start, end = parse_region(key)
        if end is None:
            end = start
        start, end = sorted((int(start), int(end)))
        lo, hi = max(0, start - r), end + r
        rows = self.loci.overlap_rows(key, lo, hi) if hi > lo else np.zeros(0, np.int64)
        if mode == "center":
            c = self.loci.centers[rows]
            rows = rows[(c >= lo) & (c < hi)]
        if linked:
            rows = rows[self.degree[rows] > 0]
        return rows

    def near(self, key, start=None, end=None, *, r: int = 0, mode: str = "overlap") -> Loci:
        """CREs near a window, as a Loci."""
        return self.loci.take(self.near_rows(key, start, end, r=r, mode=mode))

    def support(self, genes, *, r: int = 5000, mode: str = "overlap", uids: bool = True,
                rows: bool = False, linked: bool = True, backend: Optional[str] = None) -> dict:
        """CREs within ``r`` bp of each gene's TSS: ``{gene_name: [uid, ...]}``.

        The window is TSS ± ``r`` (the TSS is the gene's 5'-most base),
        computed for every gene at once. ``rows=True`` returns CRE row arrays
        instead of uids (what you want for array work)."""
        if mode not in ("overlap", "center"):
            raise ValueError(f"mode must be 'overlap' or 'center', got {mode!r}")
        from .backends.intervals import overlap_pairs
        from .genes import tss_base
        G, L = genes.genes, self.loci
        tss = tss_base(G.starts, G.ends, G.strands)
        win = L._check(Loci(G.codes, np.maximum(tss - r, 0), tss + r + 1, genome=G.genome))
        gi, ci = overlap_pairs(win, L, backend=backend)
        if mode == "center":
            c = L.centers[ci]
            ok = (c >= win.starts[gi]) & (c < win.ends[gi])
            gi, ci = gi[ok], ci[ok]
        if linked:
            ok = self.degree[ci] > 0
            gi, ci = gi[ok], ci[ok]
        o = np.lexsort((L.starts[ci], gi))
        gi, ci = gi[o], ci[o]
        cut = np.flatnonzero(np.diff(gi)) + 1
        names = G.cols["gene_name"]
        ids = G.cols["gene_id"]
        out = {}
        for grp in np.split(np.arange(len(gi)), cut) if len(gi) else []:
            k = gi[grp[0]]
            key = names[k] or ids[k]
            sel = ci[grp]
            out[key] = sel if rows else (L.uid[sel].tolist() if uids else [L[i] for i in sel])
        return out

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
        try:
            clr = cooler.Cooler(uri)
        except (KeyError, OSError) as e:
            try:
                groups = cooler.fileops.list_coolers(str(mcool))
            except Exception:                              # noqa: BLE001 — not a cooler file at all
                raise ValueError(f"{mcool}: not a .cool / .mcool file ({e})") from None
            res = [g.rsplit("/", 1)[-1] for g in groups if g.startswith("/resolutions/")]
            if res:
                raise ValueError(f"{mcool} is multi-resolution: pass resolution= one of {', '.join(res)}") from None
            raise ValueError(f"{mcool} is a single-resolution cooler: call add_mcool without resolution=") from None
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
        from .se import knee
        v = np.asarray(self.vp[key], float)
        cutoff, rows = knee(v)
        uids = self.loci.uid[rows].tolist()
        ys = v[rows][::-1]
        m, i = len(ys), len(ys) - cutoff
        if verbose and m:
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

    # ── graph algorithms: through the graph backend ─────────────────────
    def graph(self, vprops=(), *, backend: Optional[str] = None):
        """The graph as an object of the graph backend (graph-tool by default;
        'igraph', 'networkx' or 'scipy'): vertex i = Loci row i, every edge
        (cis + trans), edge columns as edge properties, ``vprops`` copied too.

        The graph-tool graph is built once and cached; its edge columns are
        re-synced on every call, so it reflects the current weights."""
        from .backends import resolve
        from .backends import graph as G
        b = resolve("graph", backend)
        if b != "graph-tool":
            return G.native(self, b, vprops=vprops)
        if self._gt is None:
            self._gt = G.native(self, "graph-tool")
        g = self._gt
        for k, v in self.ep.items():
            if k not in g.ep:
                g.ep[k] = g.new_ep("double")
            g.ep[k].a[:] = v
        for k in vprops:
            g.vp[k] = g.new_vp("double", vals=np.asarray(self.vp[k], float))
        return g

    def components(self, name: str = "component", *, backend: Optional[str] = None) -> np.ndarray:
        """Connected component per vertex, numbered by first row (rows without
        links get -1); stored as ``vp[name]``. Any graph backend gives the same labels."""
        from .backends.graph import components
        c = components(self, backend=backend)
        self.vp[name] = c
        return c

    def pagerank(self, weight: Optional[str] = None, *, damping: float = 0.85, name: str = "pagerank",
                 backend: Optional[str] = None) -> np.ndarray:
        """PageRank per vertex (``weight`` = an edge column or None), stored as ``vp[name]``."""
        from .backends.graph import pagerank
        p = pagerank(self, weight, damping=damping, backend=backend)
        self.vp[name] = p
        return p

    def to_graph_tool(self, vprops=()):
        return self.graph(vprops, backend="graph-tool")

    def to_networkx(self, vprops=(), *, loci_cols: bool = True):
        """networkx Graph (node i = Loci row i; edge columns and ``vprops`` as attributes)."""
        from .backends.graph import native
        return native(self, "networkx", vprops=vprops, loci_cols=loci_cols)

    def to_igraph(self, vprops=(), *, loci_cols: bool = True):
        """igraph Graph (vertex i = Loci row i; edge columns and ``vprops`` as attributes)."""
        from .backends.graph import native
        return native(self, "igraph", vprops=vprops, loci_cols=loci_cols)

    def to_scipy(self, weight: Optional[str] = "w"):
        """Symmetric CSR adjacency over the Loci rows (``weight=None``: 1 per edge)."""
        from .backends.graph import to_scipy
        return to_scipy(self, weight if weight in self.ep or weight is None else None)

    def to_anndata(self, *, weights=None):
        """AnnData with the vertices as ``obs`` (coordinates + vertex columns) and
        one sparse adjacency per edge column in ``obsp`` — the layout scanpy's
        graph tools read (e.g. ``obsp['n']`` as connectivities)."""
        import anndata as ad
        import pandas as pd
        L = self.loci
        obs = L.to_pandas()
        for k, v in self.vp.items():
            obs[k] = np.asarray(v)
        obs.index = pd.Index(L.names, dtype=str)
        a = ad.AnnData(X=np.zeros((len(L), 0), np.float32), obs=obs)
        for k in (weights or list(self.ep)):
            a.obsp[k] = self.to_scipy(k)
        a.uns["genomeblocks"] = {"name": self.name, "edges": int(self.n_links)}
        return a

    def draw(self, region, **kw):
        """Draw the subgraph of a region; see :func:`genomeblocks.architecture_draw.draw`."""
        from .architecture_draw import draw
        return draw(self, region, **kw)

    # ── tables out ────────────────────────────────────────────────────────
    def edges_frame(self):
        """One row per edge: src / tgt rows, uids, chromosomes, cis, edge columns."""
        import pandas as pd
        L = self.loci
        d = {"src": self.src, "tgt": self.tgt, "uid1": L.uid[self.src], "uid2": L.uid[self.tgt],
             "chrom1": L.chroms[self.src], "chrom2": L.chroms[self.tgt], "cis": self.is_cis}
        d.update(self.ep)
        return pd.DataFrame(d)

    def to_pandas(self):
        """The edge table as a DataFrame (same as :meth:`edges_frame`)."""
        return self.edges_frame()

    def to_polars(self):
        import polars as pl
        return pl.from_pandas(self.edges_frame())

    def vertices_frame(self, linked: bool = True):
        """One row per vertex (with links, unless ``linked=False``): coordinates,
        uid and every vertex column."""
        import pandas as pd
        rows = np.flatnonzero(self.degree > 0) if linked else np.arange(len(self.loci))
        df = self.loci.take(rows).to_pandas(uid=True)
        df.index = rows
        for k, v in self.vp.items():
            df[k] = np.asarray(v)[rows]
        return df

    def to_frame(self):
        """Vertex table for vertices with links (uid + vertex columns)."""
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
            json.dump({"name": self.name, "format": "genomeblocks.Architecture/2"}, f)

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
