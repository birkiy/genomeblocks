#!/usr/bin/env python3
"""Columnar prototype (branch columnar-prototype) vs the main branch.

Parts
  steps   every pipeline step, same inputs, both implementations (in-process;
          the classic modules on this branch are byte-identical to main)
  views   per-chromosome views, region / neighbour lookups, copies, graph-tool
  trans   the same pipeline with 5% inter-chromosomal loops
  scale   whole Architecture pipeline at 12.5k → 1M loops
  e2e     fresh Python processes, one per implementation — the "main" process
          imports genomeblocks from a checkout of origin/main — wall time per
          step, total, peak memory; plus reloading the saved result

Inputs: 100k CREs, the synthetic GTF (20k genes), loops.bedpe (50k loops),
the 5 kb Hi-C .mcool; loops_trans.bedpe + hic_trans_5kb.cool for ``trans``.
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import pickle
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from common import DATA, RESULTS, Recorder, timeit

ROOT = Path(__file__).resolve().parent.parent
MAIN = Path(os.environ.get("GB_MAIN_CHECKOUT", ROOT.parent / "genomeblocks-main"))
BED = str(DATA / "peaks_A_100000.bed")
GTF = str(DATA / "genes.gtf")
LOOPS = str(DATA / "loops.bedpe")
MCOOL = str(DATA / "hic_5kb.mcool")
TRANS_LOOPS = str(DATA / "loops_trans.bedpe")
TRANS_COOL = str(DATA / "hic_trans_5kb.cool")


def quiet(fn):
    def run():
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return fn()
    return run


def classic_identical_to_main() -> bool:
    """The classic modules this branch runs are exactly origin/main's."""
    r = subprocess.run(["git", "diff", "--quiet", "origin/main", "--", "genomeblocks",
                        ":(exclude)genomeblocks/columnar"], cwd=ROOT)
    return r.returncode == 0


# ── in-process step timings ──────────────────────────────────────────────────

def part_steps(rec):
    import genomeblocks.columnar as gbc
    from genomeblocks import Architecture, Genes, Loci

    def add(step, impl, t, **kw):
        rec.add(part="steps", step=step, impl=impl, seconds=t["median"], runs=t["runs"], **kw)

    add("Loci.make (100k CREs)", "main", timeit(lambda: Loci.make(BED), repeat=5))
    add("Loci.make (100k CREs)", "columnar", timeit(lambda: gbc.Loci.make(BED), repeat=5))
    add("Genes.make (GTF, 877k lines)", "main", timeit(quiet(lambda: Genes.make(GTF)), repeat=3))
    add("Genes.make (GTF, 877k lines)", "columnar", timeit(lambda: gbc.Genes.make(GTF), repeat=3))

    ol = Loci.make(BED); ol.cgr; ol.uids
    og = quiet(lambda: Genes.make(GTF))()
    L = gbc.Loci.make(BED)
    G = gbc.Genes.make(GTF)

    def reset_main(): og._annot = None
    def reset_col(): G._annot = None
    add("gene annotation index", "main", timeit(lambda: og.annot, setup=reset_main, repeat=3))
    add("gene annotation index", "columnar", timeit(lambda: G.annot, setup=reset_col, repeat=3))
    og.annot; G.annot
    add("Genes.annotations (100k CREs)", "main", timeit(lambda: og.annotations(ol), repeat=3))
    add("Genes.annotations (100k CREs)", "columnar", timeit(lambda: G.annotations(L), repeat=5))

    add("Architecture.make (50k loops)", "main",
        timeit(quiet(lambda: Architecture.make(ol, LOOPS, verbose=False)), repeat=3))
    add("Architecture.make (50k loops)", "columnar",
        timeit(lambda: gbc.Architecture.make(L, LOOPS, verbose=False), repeat=5))
    O = quiet(lambda: Architecture.make(ol, LOOPS, verbose=False))()
    A = gbc.Architecture.make(L, LOOPS, verbose=False)

    add("add_mcool (5 kb)", "main",
        timeit(quiet(lambda: O.add_mcool(ol, MCOOL, resolution=5000, verbose=False)),
               repeat=1, warmup=0))
    add("add_mcool (5 kb)", "columnar",
        timeit(lambda: A.add_mcool(MCOOL, resolution=5000, verbose=False), repeat=5))
    add("normalize (power-law O/E)", "main", timeit(quiet(lambda: O.normalize(ol, verbose=False)), repeat=3))
    add("normalize (power-law O/E)", "columnar", timeit(lambda: A.normalize(verbose=False), repeat=5))
    add("annotate (labels + genes)", "main",
        timeit(quiet(lambda: O.annotate(ol, og, key="n", verbose=False)), repeat=3))
    add("annotate (labels + genes)", "columnar", timeit(lambda: A.annotate(G, verbose=False), repeat=5))
    add("strength", "main", timeit(quiet(lambda: O.strength(key="n", verbose=False)), repeat=3))
    add("strength", "columnar", timeit(lambda: A.strength(verbose=False), repeat=5))
    add("prime_hubs", "main", timeit(quiet(lambda: O.prime_hubs(verbose=False)), repeat=3))
    add("prime_hubs", "columnar", timeit(lambda: A.prime_hubs(verbose=False), repeat=5))

    with tempfile.TemporaryDirectory() as d:
        pk, pq = os.path.join(d, "arch.pkl"), os.path.join(d, "arch")
        add("save", "main", timeit(lambda: pickle.dump(O, open(pk, "wb")), repeat=3))
        add("save", "columnar", timeit(lambda: A.save(pq), repeat=5))
        add("load", "main", timeit(quiet(lambda: pickle.load(open(pk, "rb"))), repeat=3))
        add("load", "columnar", timeit(lambda: gbc.Architecture.load(pq), repeat=5))
        size_main = os.path.getsize(pk)
        size_col = sum(os.path.getsize(os.path.join(pq, f)) for f in os.listdir(pq))
    rec.add(part="steps_meta", n_vertices=O.n_loci, n_edges=O.n_links,
            file_mb_main=size_main / 1e6, file_mb_columnar=size_col / 1e6,
            identical_to_main=classic_identical_to_main())
    return ol, og, L, G, O, A


# ── views & lookups ─────────────────────────────────────────────────────────

def part_views(rec, ol, O, A):
    import graph_tool.all as gt

    def add(op, impl, t, **kw):
        rec.add(part="views", op=op, impl=impl, seconds=t["median"], runs=t["runs"], **kw)

    uid = np.array(list(O.vp.uid), dtype=object)
    on8 = np.array([u.startswith("chr8:") for u in uid])
    add("one chromosome (chr8)", "main",
        timeit(quiet(lambda: O.subgraph(filter_func=lambda v: on8[int(v)])), repeat=3),
        n_edges=A.chrom("chr8").n_links)
    add("one chromosome (chr8)", "columnar", timeit(lambda: A.chrom("chr8"), repeat=50))
    add("all cis / all trans split", "columnar", timeit(lambda: (A.cis, A.trans), repeat=50))
    region = ("chr8", 20_000_000, 40_000_000)
    in_reg = np.array([u.startswith("chr8:") and region[1] <= int(u.split(":")[1].split("-")[0]) < region[2]
                       for u in uid])
    add("region chr8:20-40 Mb", "main",
        timeit(quiet(lambda: O.subgraph(filter_func=lambda v: in_reg[int(v)])), repeat=3))
    add("region chr8:20-40 Mb", "columnar", timeit(lambda: A.region(*region), repeat=20))

    rng = np.random.default_rng(0)
    pick = rng.choice(uid, 1000, replace=False).tolist()
    add("neighbours of 1,000 CREs", "main", timeit(lambda: [O[u] for u in pick], repeat=5))
    A._csr = None
    add("neighbours of 1,000 CREs", "columnar", timeit(lambda: [A[u] for u in pick], repeat=5))
    add("copy", "main", timeit(quiet(lambda: O.copy()), repeat=1, warmup=0))
    add("copy", "columnar", timeit(lambda: A.copy(), repeat=5))

    def build():
        A._gt = None
        return A.graph()
    add("graph-tool Graph from tables", "columnar", timeit(build, repeat=5))
    add("connected components (graph-tool)", "main", timeit(lambda: gt.label_components(O), repeat=5))
    A.graph()
    add("connected components (graph-tool)", "columnar", timeit(lambda: A.components(), repeat=5))
    add("pagerank (graph-tool)", "main", timeit(lambda: gt.pagerank(O, weight=O.ep.n), repeat=3))
    g = A.graph()
    add("pagerank (graph-tool)", "columnar", timeit(lambda: gt.pagerank(g, weight=g.ep.n), repeat=3))


# ── inter-chromosomal edges ─────────────────────────────────────────────────

def part_trans(rec, ol, og, L, G):
    import genomeblocks.columnar as gbc
    from genomeblocks import Architecture

    def add(step, impl, seconds, **kw):
        rec.add(part="trans", step=step, impl=impl, seconds=seconds, **kw)

    A = gbc.Architecture.make(L, TRANS_LOOPS, verbose=False)
    O = quiet(lambda: Architecture.make(ol, TRANS_LOOPS, verbose=False))()
    steps = [
        ("make", lambda: gbc.Architecture.make(L, TRANS_LOOPS, verbose=False),
         quiet(lambda: Architecture.make(ol, TRANS_LOOPS, verbose=False))),
        ("add_mcool", lambda: A.add_mcool(TRANS_COOL, verbose=False),
         quiet(lambda: O.add_mcool(ol, TRANS_COOL, verbose=False))),
        ("normalize", lambda: A.normalize(verbose=False), quiet(lambda: O.normalize(ol, verbose=False))),
        ("annotate", lambda: A.annotate(G, verbose=False), quiet(lambda: O.annotate(ol, og, verbose=False))),
        ("strength", lambda: A.strength(verbose=False), quiet(lambda: O.strength(verbose=False))),
    ]
    main_ok = True
    for name, new, old in steps:
        t = timeit(new, repeat=3)
        add(name, "columnar", t["median"])
        if not main_ok:
            add(name, "main", None, status="not reached")
            continue
        try:
            t0 = time.perf_counter(); old()
            add(name, "main", time.perf_counter() - t0, status="ok")
        except Exception as e:  # noqa: BLE001 — record the failure, keep going
            add(name, "main", None, status=f"{type(e).__name__}: {e}")
            main_ok = False
    t = A.trans
    cis = A.cis
    comp_all = A.components()
    comp_cis = A.cis.components(name="cc")
    nc = lambda c: int(len(np.unique(c[c >= 0])))
    big_all = np.bincount(comp_all[comp_all >= 0]).max()
    big_cis = np.bincount(comp_cis[comp_cis >= 0]).max()
    blocks = A.block_counts()
    rec.add(part="trans_meta", n_edges=A.n_links, n_cis=cis.n_links, n_trans=t.n_links,
            trans_with_contacts=int((t.ep.w > 0).sum()),
            cis_with_contacts=int((cis.ep.w > 0).sum()),
            trans_n_median=float(np.median(t.ep.n[t.ep.w > 0])),
            cis_n_median=float(np.median(cis.ep.n[cis.ep.w > 0])),
            components_all=nc(comp_all), components_cis_only=nc(comp_cis),
            largest_component_all=int(big_all), largest_component_cis_only=int(big_cis),
            genes_via_trans=int(((A.vp.gene != "") & (A.vp.annot != "Promoter-TSS")).sum()),
            block_counts={"chroms": list(blocks.index), "matrix": blocks.to_numpy().tolist()})


# ── scaling ──────────────────────────────────────────────────────────────────

def part_scale(rec, ol, og, L, G, main_upto=200_000):
    import genomeblocks.columnar as gbc
    from genomeblocks import Architecture
    for n in (12_500, 50_000, 200_000, 1_000_000):
        loops = LOOPS if n == 50_000 else str(DATA / f"loops_{n}.bedpe")
        t0 = time.perf_counter()
        A = gbc.Architecture.make(L, loops, verbose=False)
        A.add_mcool(MCOOL, resolution=5000, verbose=False).normalize(verbose=False)
        A.annotate(G, verbose=False).strength(verbose=False)
        rec.add(part="scale", loops=n, impl="columnar", seconds=time.perf_counter() - t0,
                n_edges=A.n_links)
        if n > main_upto:
            continue
        t0 = time.perf_counter()
        O = quiet(lambda: Architecture.make(ol, loops, verbose=False))()
        quiet(lambda: O.add_mcool(ol, MCOOL, resolution=5000, verbose=False))()
        quiet(lambda: O.normalize(ol, verbose=False))()
        quiet(lambda: O.annotate(ol, og, key="n", verbose=False))()
        quiet(lambda: O.strength(key="n", verbose=False))()
        rec.add(part="scale", loops=n, impl="main", seconds=time.perf_counter() - t0, n_edges=O.n_links)


# ── fresh-process end to end ─────────────────────────────────────────────────

E2E_MAIN = r'''
import contextlib, io, json, pickle, sys, time
T = {}; t00 = t = time.perf_counter()
def lap(k):
    global t
    T[k] = time.perf_counter() - t; t = time.perf_counter()
from genomeblocks import Architecture, Genes, Loci
import genomeblocks.architecture, genomeblocks.genes; lap("import")
q = lambda: contextlib.redirect_stdout(io.StringIO())
L = Loci.make(BED); lap("Loci.make")
with q(), contextlib.redirect_stderr(io.StringIO()):
    G = Genes.make(GTF); lap("Genes.make")
    A = Architecture.make(L, LOOPS, verbose=False); lap("make")
    A.add_mcool(L, MCOOL, resolution=5000, verbose=False); lap("add_mcool")
    A.normalize(L, verbose=False); lap("normalize")
    A.annotate(L, G, key="n", verbose=False); lap("annotate")
    A.strength(key="n", verbose=False); lap("strength")
    A.prime_hubs(verbose=False); lap("prime_hubs")
pickle.dump(A, open(OUT, "wb")); lap("save")
T["total"] = time.perf_counter() - t00
import genomeblocks
print(json.dumps({"steps": T, "peak_mb": PEAK(),
                  "module": genomeblocks.__file__}))
'''

E2E_COL = r'''
import json, sys, time
T = {}; t00 = t = time.perf_counter()
def lap(k):
    global t
    T[k] = time.perf_counter() - t; t = time.perf_counter()
import genomeblocks.columnar as gbc
from genomeblocks.columnar import Architecture, Genes, Loci; lap("import")
L = Loci.make(BED); lap("Loci.make")
G = Genes.make(GTF); lap("Genes.make")
A = Architecture.make(L, LOOPS, verbose=False); lap("make")
A.add_mcool(MCOOL, resolution=5000, verbose=False); lap("add_mcool")
A.normalize(verbose=False); lap("normalize")
A.annotate(G, verbose=False); lap("annotate")
A.strength(verbose=False); lap("strength")
A.prime_hubs(verbose=False); lap("prime_hubs")
A.save(OUT); lap("save")
T["total"] = time.perf_counter() - t00
import genomeblocks
print(json.dumps({"steps": T, "peak_mb": PEAK(),
                  "module": genomeblocks.__file__}))
'''

LOAD_MAIN = r'''
import json, pickle, time, contextlib, io
t = time.perf_counter()
with contextlib.redirect_stdout(io.StringIO()):
    import genomeblocks.architecture
    A = pickle.load(open(OUT, "rb"))
print(json.dumps({"seconds": time.perf_counter() - t, "n_edges": A.n_links,
                  "peak_mb": PEAK()}))
'''

LOAD_COL = r'''
import json, time
t = time.perf_counter()
from genomeblocks.columnar import Architecture
A = Architecture.load(OUT)
print(json.dumps({"seconds": time.perf_counter() - t, "n_edges": A.n_links,
                  "peak_mb": PEAK()}))
'''


# peak resident memory of *this* process: VmHWM starts fresh at exec, whereas
# ru_maxrss would carry over the (large) benchmark parent's peak
PEAK_SRC = """
def PEAK():
    for line in open("/proc/self/status"):
        if line.startswith("VmHWM:"):
            return int(line.split()[1]) / 1024
"""


def _run(code, pythonpath, **consts):
    env = dict(os.environ, PYTHONPATH=str(pythonpath))
    head = "".join(f"{k} = {v!r}\n" for k, v in consts.items()) + PEAK_SRC
    r = subprocess.run([sys.executable, "-c", head + code], env=env, cwd=str(DATA),
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def part_e2e(rec):
    if not (MAIN / "genomeblocks").is_dir():
        raise SystemExit(f"main checkout missing: git worktree add {MAIN} origin/main")
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=MAIN, capture_output=True,
                          text=True).stdout.strip()
    with tempfile.TemporaryDirectory() as d:
        for impl, code, load, path, out in (
                ("main", E2E_MAIN, LOAD_MAIN, MAIN, os.path.join(d, "a.pkl")),
                ("columnar", E2E_COL, LOAD_COL, ROOT, os.path.join(d, "a"))):
            r = _run(code, path, BED=BED, GTF=GTF, LOOPS=LOOPS, MCOOL=MCOOL, OUT=out)
            assert r["module"].startswith(str(path)), r["module"]
            for k, v in r["steps"].items():
                rec.add(part="e2e", impl=impl, step=k, seconds=v)
            rec.add(part="e2e_mem", impl=impl, peak_mb=r["peak_mb"], main_commit=head)
            ld = _run(load, path, OUT=out)
            rec.add(part="e2e_load", impl=impl, seconds=ld["seconds"], peak_mb=ld["peak_mb"],
                    n_edges=ld["n_edges"])


if __name__ == "__main__":
    parts = sys.argv[1:] or ["steps", "views", "trans", "scale", "e2e"]
    rec = Recorder("prototype")
    state = None
    if {"steps", "views", "trans", "scale"} & set(parts):
        state = part_steps(rec) if "steps" in parts else None
        if state is None:
            import genomeblocks.columnar as gbc
            from genomeblocks import Architecture, Genes, Loci
            ol = Loci.make(BED); ol.cgr; ol.uids
            og = quiet(lambda: Genes.make(GTF))(); og.annot
            L, G = gbc.Loci.make(BED), gbc.Genes.make(GTF)
            O = quiet(lambda: Architecture.make(ol, LOOPS, verbose=False))()
            for f in (lambda: O.add_mcool(ol, MCOOL, resolution=5000, verbose=False),
                      lambda: O.normalize(ol, verbose=False)):
                quiet(f)()
            A = gbc.Architecture.make(L, LOOPS, verbose=False).add_mcool(MCOOL, verbose=False)
            A.normalize(verbose=False)
            state = (ol, og, L, G, O, A)
    ol, og, L, G, O, A = state if state else (None,) * 6
    if "views" in parts:
        part_views(rec, ol, O, A)
    if "trans" in parts:
        part_trans(rec, ol, og, L, G)
    if "scale" in parts:
        part_scale(rec, ol, og, L, G)
    if "e2e" in parts:
        part_e2e(rec)
    rec.save()
