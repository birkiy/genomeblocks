#!/usr/bin/env python3
"""Architecture pipeline: per-step cost on a realistic graph.

Input: 100k CREs, 50k loops (BEDPE), a 5 kb .mcool built from the 5M Hi-C
pairs, and the synthetic GTF. Runs the README pipeline step by step:

    make -> add_mcool -> normalize -> annotate -> strength -> prime_hubs

plus the graph algorithms (components, pagerank) through each installed
graph backend, and save / load. Every step works on the edge table, so
there is no object-vs-array comparison left to make; the point is the cost
of each step at this size.
"""
from __future__ import annotations

import contextlib
import io
import subprocess
import sys
import tempfile
from pathlib import Path


from common import DATA, Recorder, timeit

import genomeblocks as gb
from genomeblocks import Architecture, Genes, Loci

BIN = Path(sys.executable).parent
MCOOL = DATA / "hic_5kb.mcool"


def quiet(fn):
    """Architecture prints progress; keep timings clean."""
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
    subprocess.run([str(BIN / "cooler"), "zoomify", "-r", "5000", "-o", str(MCOOL), str(cool)], check=True)


if __name__ == "__main__":
    rec = Recorder("architecture")
    build_mcool()
    cre = Loci.make(str(DATA / "peaks_A_100000.bed"))
    genes = Genes.make(str(DATA / "genes.gtf"))
    bedpe = str(DATA / "loops.bedpe")

    t = timeit(quiet(lambda: Architecture.make(cre, bedpe, r=2500, verbose=False)), repeat=3)
    A = quiet(lambda: Architecture.make(cre, bedpe, r=2500, verbose=False))()
    rec.add(step="make (50k loops -> graph)", seconds=t["median"], runs=t["runs"],
            n_vertices=A.n_loci, n_edges=A.n_links)

    t = timeit(quiet(lambda: A.add_mcool(str(MCOOL), resolution=5000, verbose=False)), repeat=1, warmup=0)
    rec.add(step="add_mcool (5 kb)", seconds=t["median"], runs=t["runs"], n_edges=A.n_links)

    t = timeit(quiet(lambda: A.normalize(verbose=False)), repeat=3)
    rec.add(step="normalize (power-law O/E)", seconds=t["median"], runs=t["runs"], n_edges=A.n_links)

    t = timeit(quiet(lambda: A.annotate(genes, key="n", verbose=False)), repeat=3)
    rec.add(step="annotate (labels + genes)", seconds=t["median"], runs=t["runs"], n_vertices=A.n_loci)

    t = timeit(quiet(lambda: A.strength(verbose=False)), repeat=5)
    rec.add(step="strength (numpy bincount)", seconds=t["median"], runs=t["runs"], n_edges=A.n_links)

    t = timeit(quiet(lambda: A.prime_hubs(verbose=False)), repeat=3)
    rec.add(step="prime_hubs (elbow on strength)", seconds=t["median"], runs=t["runs"])

    t = timeit(lambda: A.support(genes, r=5000), repeat=3)
    rec.add(step="support (CREs per TSS ± 5 kb)", seconds=t["median"], runs=t["runs"])

    for b in gb.backends.families()["graph"]:
        if not gb.backends.installed("graph", b):
            continue
        t = timeit(lambda: A.components(backend=b), repeat=3)
        rec.add(step="components", engine=b, seconds=t["median"], runs=t["runs"], n_edges=A.n_links)
        t = timeit(lambda: A.pagerank("n", backend=b), repeat=3)
        rec.add(step="pagerank (O/E weights)", engine=b, seconds=t["median"], runs=t["runs"], n_edges=A.n_links)

    with tempfile.TemporaryDirectory() as td:
        t = timeit(lambda: A.save(f"{td}/A"), repeat=3)
        rec.add(step="save (parquet)", seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: Architecture.load(f"{td}/A"), repeat=3)
        rec.add(step="load (parquet)", seconds=t["median"], runs=t["runs"])
    rec.save()
