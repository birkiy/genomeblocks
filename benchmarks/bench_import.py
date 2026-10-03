#!/usr/bin/env python3
"""Import cost: the lazy ``__getattr__`` top level vs touching each block.

Each statement runs in a fresh interpreter (median of 15); the bare
interpreter start-up is subtracted.
"""
from __future__ import annotations

import statistics
import subprocess
import sys
import time

from common import Recorder

STMTS = [
    ("import genomeblocks", "import genomeblocks"),
    ("+ Locus", "from genomeblocks import Locus"),
    ("+ Loci (pulls signal, motifs, atlas, bedpe)", "from genomeblocks import Loci"),
    ("+ Genes", "from genomeblocks import Genes"),
    ("+ browser (matplotlib)", "from genomeblocks import browser"),
    ("+ Architecture (graph-tool)", "from genomeblocks import Architecture"),
    ("everything", "from genomeblocks import Loci, Genes, Atlas, browser, Architecture, "
                   "scan_motifs, compare_heatmap, coverage"),
]


def wall(stmt, n=15):
    runs = []
    for _ in range(n):
        t0 = time.perf_counter()
        subprocess.run([sys.executable, "-W", "ignore", "-c", stmt], check=True,
                       capture_output=True)
        runs.append(time.perf_counter() - t0)
    return statistics.median(runs), runs


if __name__ == "__main__":
    rec = Recorder("import")
    base, _ = wall("pass")
    rec.add(stmt="python -c pass", seconds=base)
    for label, stmt in STMTS:
        t, runs = wall(stmt)
        rec.add(stmt=label, seconds=t - base, runs=[r - base for r in runs])
    rec.save(baseline=base)
