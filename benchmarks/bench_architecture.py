#!/usr/bin/env python3
"""Architecture (graph-tool) pipeline: per-step cost on a realistic graph.

Input: 100k CREs, 50k loops (BEDPE), a 5 kb .mcool built from the 5M Hi-C
pairs, and the synthetic GTF. Runs the README pipeline:

    make -> add_mcool -> normalize -> annotate -> strength

and, where a step has a vectorised core, compares it with the per-element
loop it replaced (annotate stage 2) or could use (strength).
"""
from __future__ import annotations

import contextlib
import io
import subprocess
import sys
from pathlib import Path

import numpy as np

from common import DATA, Recorder, timeit

from genomeblocks import Architecture, Genes, Loci

BIN = Path(sys.executable).parent
MCOOL = DATA / "hic_5kb.mcool"


def quiet(fn):
    """Architecture prints progress and tqdm bars; keep timings clean."""
    def run():
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return fn()
    return run


def build_mcool():
    if MCOOL.exists():
        return
    cool = DATA / "hic_5kb.cool"
    subprocess.run([str(BIN / "cooler"), "cload", "pairs", "-c1", "2", "-p1", "3",
                    "-c2", "4", "-p2", "5", f"{DATA / 'hg38.chrom.sizes'}:5000",
                    str(DATA / "hic.pairs"), str(cool)], check=True)
    subprocess.run([str(BIN / "cooler"), "zoomify", "-r", "5000", "-o", str(MCOOL),
                    str(cool)], check=True)


def annotate_stage2_loop(G, key, is_prom, gene_pm):
    """The per-vertex / per-edge loop that annotate() stage 2 replaced."""
    out = {}
    for v in G.vertices():
        if is_prom[int(v)]:
            continue
        best, best_w = "", -np.inf
        for e in v.all_edges():
            u = e.target() if e.source() == v else e.source()
            if not is_prom[int(u)] or not gene_pm[u]:
                continue
            w = G.ep[key][e]
            if w > best_w:
                best, best_w = gene_pm[u], w
        if best:
            out[int(v)] = best
    return out


def annotate_stage2_vec(G, key, is_prom, gene_pm):
    """annotate() stage 2, as implemented in genomeblocks."""
    ew = G.get_edges(eprops=[G.ep[key]])
    src = ew[:, 0].astype(np.int64); tgt = ew[:, 1].astype(np.int64); w = ew[:, 2]
    gene_arr = np.asarray(list(gene_pm), dtype=object)
    eligible = is_prom & (gene_arr != "")
    c1 = ~is_prom[src] & eligible[tgt]
    c2 = ~is_prom[tgt] & eligible[src]
    idx = np.concatenate([src[c1], tgt[c2]])
    pg = np.concatenate([gene_arr[tgt[c1]], gene_arr[src[c2]]])
    we = np.concatenate([w[c1], w[c2]])
    order = np.lexsort((-we, idx))
    _, first = np.unique(idx[order], return_index=True)
    return dict(zip(idx[order][first].tolist(), pg[order][first].tolist()))


if __name__ == "__main__":
    import graph_tool.all as gt
    rec = Recorder("architecture")
    build_mcool()
    cre = Loci.make(str(DATA / "peaks_A_100000.bed"))
    cre.cgr, cre.uids
    genes = Genes.make(str(DATA / "genes.gtf"))
    bedpe = str(DATA / "loops.bedpe")

    t = timeit(quiet(lambda: Architecture.make(cre, bedpe, r=2500, verbose=False)), repeat=3)
    G = quiet(lambda: Architecture.make(cre, bedpe, r=2500, verbose=False))()
    rec.add(step="make (50k loops -> graph)", seconds=t["median"], runs=t["runs"],
            n_vertices=G.n_loci, n_edges=G.n_links)

    t = timeit(quiet(lambda: G.add_mcool(cre, str(MCOOL), resolution=5000, verbose=False)),
               repeat=1, warmup=0)
    rec.add(step="add_mcool (5 kb)", seconds=t["median"], runs=t["runs"], n_edges=G.n_links)

    t = timeit(quiet(lambda: G.normalize(cre, verbose=False)), repeat=3)
    rec.add(step="normalize (power-law O/E)", seconds=t["median"], runs=t["runs"],
            n_edges=G.n_links)

    t = timeit(quiet(lambda: G.annotate(cre, genes, key="n", verbose=False)), repeat=3)
    rec.add(step="annotate (total)", seconds=t["median"], runs=t["runs"], n_vertices=G.n_loci)

    # stage-2 core: vectorised vs the loop it replaced
    is_prom = np.array(["Promoter" in a for a in G.vp.annot], bool)
    gene_pm = G.vp.gene
    res = {}
    t = timeit(lambda: res.__setitem__("v", annotate_stage2_vec(G, "n", is_prom, gene_pm)), repeat=5)
    rec.add(step="annotate stage 2: vectorised (genomeblocks)", seconds=t["median"],
            runs=t["runs"], n_edges=G.n_links)
    t = timeit(lambda: res.__setitem__("l", annotate_stage2_loop(G, "n", is_prom, gene_pm)),
               repeat=3)
    rec.add(step="annotate stage 2: per-vertex loop", seconds=t["median"], runs=t["runs"],
            n_edges=G.n_links, same_result=res["v"] == res["l"])

    t = timeit(quiet(lambda: G.strength(key="n", verbose=False)), repeat=3)
    rec.add(step="strength: per-edge Python loop (genomeblocks)", seconds=t["median"],
            runs=t["runs"], n_edges=G.n_links)
    ref = np.asarray(G.vp.strength.a).copy()
    t = timeit(lambda: res.__setitem__("s", gt.incident_edges_op(G, "out", "sum", G.ep.n)),
               repeat=5)
    rec.add(step="strength: graph-tool incident_edges_op", seconds=t["median"],
            runs=t["runs"], n_edges=G.n_links,
            same_result=bool(np.allclose(np.asarray(res["s"].a), ref)))
    rec.save()
