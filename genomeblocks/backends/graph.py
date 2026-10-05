"""Graph backends for the Architecture: graph-tool (default), igraph, networkx, scipy.

The Architecture is two tables (vertices = Loci rows, edges = src/tgt + edge
columns); a backend turns them into its own graph object when an algorithm
needs one. Every function here returns plain numpy arrays aligned to the
Loci rows, so results do not depend on the engine:

    native(A, backend)         the engine's own graph (vertex i = Loci row i),
                               edge columns as edge attributes
    components(A, backend)     connected component per row (-1 for rows without links),
                               numbered by their first row, so every engine agrees
    pagerank(A, ...)           PageRank per row
    degree / strength          plain numpy (no engine needed)
    layout(A, kind, backend)   2-D positions for drawing
"""
from __future__ import annotations

import numpy as np

from . import resolve


def _weights(A, weight):
    if weight is None:
        return None
    if isinstance(weight, str):
        return np.asarray(A.ep[weight], np.float64)
    return np.asarray(weight, np.float64)


def to_scipy(A, weight="w", *, n=None):
    """Symmetric CSR adjacency (n x n over Loci rows); ``weight=None`` gives 1 per edge.
    Parallel edges are summed."""
    from scipy.sparse import coo_matrix
    n = len(A.loci) if n is None else n
    if weight is None:
        w = np.ones(A.n_links)
    elif isinstance(weight, str):
        if weight not in A.ep:
            raise ValueError(f"edge column {weight!r} not found; have: {', '.join(A.ep) or 'none'}")
        w = np.asarray(A.ep[weight], np.float64)
    else:
        w = np.asarray(weight, np.float64)
    r = np.concatenate([A.src, A.tgt])
    c = np.concatenate([A.tgt, A.src])
    return coo_matrix((np.concatenate([w, w]), (r, c)), shape=(n, n)).tocsr()


def native(A, backend=None, *, vprops=(), loci_cols: bool = False):
    """The Architecture as a graph of the chosen engine (vertex i = Loci row i).

    Edge columns become edge properties/attributes; ``vprops`` names vertex
    columns to copy too (numeric and string). ``loci_cols=True`` also adds the
    chromosome / start / end of every vertex.
    """
    b = resolve("graph", backend)
    n = len(A.loci)
    vcols = {k: np.asarray(A.vp[k]) for k in vprops}
    if loci_cols:
        vcols.update(chrom=A.loci.chroms, start=A.loci.starts, end=A.loci.ends)
    if b == "graph-tool":
        import graph_tool.all as gt
        g = gt.Graph(directed=False)
        g.add_vertex(n)
        g.add_edge_list(np.column_stack([A.src, A.tgt]))
        for k, v in A.ep.items():
            g.ep[k] = g.new_ep("double", vals=np.asarray(v, np.float64))
        for k, v in vcols.items():
            if v.dtype.kind in "biuf":
                g.vp[k] = g.new_vp("double", vals=v.astype(np.float64))
            else:
                p = g.new_vp("string")
                for vert, s in zip(g.vertices(), v.astype(str).tolist()):
                    p[vert] = s
                g.vp[k] = p
        return g
    if b == "igraph":
        import igraph as ig
        g = ig.Graph(n=n, edges=np.column_stack([A.src, A.tgt]).tolist(), directed=False)
        for k, v in A.ep.items():
            g.es[k] = np.asarray(v, np.float64).tolist()
        for k, v in vcols.items():
            g.vs[k] = v.tolist()
        return g
    if b == "networkx":
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(range(n))
        # nx.Graph holds one edge per pair: parallel edges are merged and their
        # numeric columns summed, the same sum to_scipy() puts in the adjacency
        s, t, inv = _merged_pairs(A)
        keys = list(A.ep)
        cols = [np.bincount(inv, weights=np.asarray(A.ep[k], np.float64), minlength=len(s)).tolist()
                for k in keys]
        g.add_edges_from((a, b_, dict(zip(keys, vals))) for a, b_, *vals in zip(s.tolist(), t.tolist(), *cols))
        for k, v in vcols.items():
            nx.set_node_attributes(g, dict(enumerate(v.tolist())), k)
        return g
    return to_scipy(A, "w" if "w" in A.ep else None)


def _merged_pairs(A):
    """Unique undirected (src, tgt) pairs and, per edge, the index of its pair."""
    lo, hi = np.minimum(A.src, A.tgt), np.maximum(A.src, A.tgt)
    n = len(A.loci)
    key, inv = np.unique(lo.astype(np.int64) * max(n, 1) + hi, return_inverse=True)
    return key // max(n, 1), key % max(n, 1), inv.ravel()


def _canonical(labels, active):
    """Renumber component labels by their first row; inactive rows get -1."""
    labels = np.asarray(labels, np.int64)
    out = np.full(len(labels), -1, np.int64)
    if not active.any():
        return out
    lab = labels[active]
    _, first = np.unique(lab, return_index=True)
    order = np.argsort(first)                        # components by their first row
    rank = np.empty(len(order), np.int64)
    rank[order] = np.arange(len(order))
    uniq = np.unique(lab)
    out[active] = rank[np.searchsorted(uniq, lab)]
    return out


def components(A, backend=None) -> np.ndarray:
    """Connected component per Loci row (-1 for rows without links)."""
    b = resolve("graph", backend)
    active = A.degree > 0
    if b == "graph-tool":
        import graph_tool.all as gt
        comp, _ = gt.label_components(A.graph(backend="graph-tool"))
        lab = np.asarray(comp.a, np.int64).copy()
    elif b == "igraph":
        lab = np.asarray(native(A, "igraph").connected_components().membership, np.int64)
    elif b == "networkx":
        import networkx as nx
        lab = np.zeros(len(A.loci), np.int64)
        for k, comp in enumerate(nx.connected_components(native(A, "networkx"))):
            lab[list(comp)] = k
    else:
        from scipy.sparse.csgraph import connected_components
        _, lab = connected_components(to_scipy(A, None), directed=False)
    return _canonical(lab, active)


def pagerank(A, weight=None, *, damping: float = 0.85, backend=None, tol: float = 1e-10) -> np.ndarray:
    """PageRank per Loci row (rows without links included, as every engine does)."""
    b = resolve("graph", backend)
    w = _weights(A, weight)
    if b == "graph-tool":
        import graph_tool.all as gt
        g = A.graph(backend="graph-tool")
        ep = g.new_ep("double", vals=w) if w is not None else None
        return np.asarray(gt.pagerank(g, damping=damping, weight=ep, epsilon=tol).a, np.float64).copy()
    if b == "igraph":
        g = native(A, "igraph")
        return np.asarray(g.pagerank(damping=damping, weights=None if w is None else w.tolist(),
                                     implementation="prpack"), np.float64)
    if b == "networkx":
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(range(len(A.loci)))
        ww = np.ones(A.n_links) if w is None else np.asarray(w, np.float64)
        s, t, inv = _merged_pairs(A)                      # parallel edges: weights summed (as scipy)
        g.add_weighted_edges_from(zip(s.tolist(), t.tolist(),
                                      np.bincount(inv, weights=ww, minlength=len(s)).tolist()))
        pr = nx.pagerank(g, alpha=damping, weight="weight", tol=tol, max_iter=1000)
        return np.array([pr[i] for i in range(len(A.loci))], np.float64)
    # scipy: power iteration on the column-stochastic matrix, dangling mass spread uniformly
    M = to_scipy(A, None if w is None else w)
    n = M.shape[0]
    out_w = np.asarray(M.sum(axis=1)).ravel()
    dangling = out_w == 0
    inv = np.where(dangling, 0.0, 1.0 / np.where(dangling, 1.0, out_w))
    P = M.multiply(inv[:, None]).T.tocsr()
    x = np.full(n, 1.0 / n)
    for _ in range(1000):
        x_new = damping * (P @ x + x[dangling].sum() / n) + (1 - damping) / n
        if np.abs(x_new - x).sum() < tol * n:
            x = x_new
            break
        x = x_new
    return x / x.sum()


def layout(A, rows, *, kind: str = "spring", backend=None, seed: int = 0) -> np.ndarray:
    """2-D positions (len(rows) x 2) for the subgraph induced by ``rows``.

    ``kind``: 'spring' (the engine's force-directed layout), 'circular' or
    'genomic' (x = genomic position, y = 0: a linear browser-like layout).
    """
    rows = np.asarray(rows, np.int64)
    m = len(rows)
    if m == 0:
        return np.zeros((0, 2))
    if kind == "circular":
        t = np.linspace(0, 2 * np.pi, m, endpoint=False)
        return np.column_stack([np.cos(t), np.sin(t)])
    if kind == "genomic":
        return np.column_stack([A.loci.centers[rows].astype(float), np.zeros(m)])
    if kind != "spring":
        raise ValueError(f"kind must be 'spring', 'circular' or 'genomic', got {kind!r}")
    pos_of = np.full(len(A.loci), -1, np.int64)
    pos_of[rows] = np.arange(m)
    keep = (pos_of[A.src] >= 0) & (pos_of[A.tgt] >= 0)
    s, t = pos_of[A.src[keep]], pos_of[A.tgt[keep]]
    b = resolve("graph", backend)
    if b == "graph-tool":
        import graph_tool.all as gt
        gt.seed_rng(seed)
        np.random.seed(seed)
        g = gt.Graph(directed=False)
        g.add_vertex(m)
        g.add_edge_list(np.column_stack([s, t]))
        return np.asarray(gt.sfdp_layout(g).get_2d_array([0, 1]).T, np.float64)
    if b == "igraph":
        import igraph as ig
        import random
        random.seed(seed)
        g = ig.Graph(n=m, edges=np.column_stack([s, t]).tolist(), directed=False)
        return np.asarray(g.layout_fruchterman_reingold().coords, np.float64)
    if b == "networkx":
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(range(m))
        g.add_edges_from(zip(s.tolist(), t.tolist()))
        p = nx.spring_layout(g, seed=seed)
        return np.array([p[i] for i in range(m)], np.float64)
    # scipy: spectral layout from the two smallest non-trivial Laplacian eigenvectors,
    # deterministic (dense solver for small graphs, seeded start vector otherwise)
    import warnings
    from scipy.sparse import coo_matrix, diags
    from scipy.sparse.linalg import eigsh
    circle = np.column_stack([np.cos(np.linspace(0, 2 * np.pi, m, endpoint=False)),
                              np.sin(np.linspace(0, 2 * np.pi, m, endpoint=False))])
    if m < 3 or not len(s):
        return circle
    W = coo_matrix((np.ones(2 * len(s)), (np.r_[s, t], np.r_[t, s])), shape=(m, m)).tocsr()
    L = (diags(np.asarray(W.sum(1)).ravel()) - W).astype(float)
    try:
        if m <= 2000:
            vals, vecs = np.linalg.eigh(L.toarray())
            vecs = vecs[:, 1:3]
        else:
            rng = np.random.default_rng(seed)
            vals, vecs = eigsh(L, k=min(3, m - 1), which="SM", v0=rng.normal(size=m))
            vecs = vecs[:, 1:3]
    except Exception as e:                                 # noqa: BLE001
        warnings.warn(f"spectral layout failed ({e}); using a circular layout", RuntimeWarning, stacklevel=2)
        return circle
    pos = np.zeros((m, 2))
    pos[:, :vecs.shape[1]] = vecs
    for j in range(vecs.shape[1]):                         # fix the sign: largest entry positive
        k = np.argmax(np.abs(pos[:, j]))
        if pos[k, j] < 0:
            pos[:, j] *= -1
    rng = np.random.default_rng(seed)
    scale = float(np.ptp(pos, axis=0).max()) or 1.0
    pos += rng.normal(scale=1e-3 * scale, size=pos.shape)  # separate structurally equivalent rows
    return pos
