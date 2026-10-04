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


def importtime_breakdown(stmt, n=7):
    """Median cumulative import time (s) of each genomeblocks submodule."""
    acc = {}
    for _ in range(n):
        err = subprocess.run([sys.executable, "-X", "importtime", "-c", stmt],
                             capture_output=True, text=True).stderr
        for line in err.splitlines():
            parts = [p.strip() for p in line.replace("import time:", "").split("|")]
            if len(parts) == 3 and parts[2].startswith("genomeblocks.") and parts[1].isdigit():
                acc.setdefault(parts[2], []).append(int(parts[1]) / 1e6)
    return {k: statistics.median(v) for k, v in acc.items()}


if __name__ == "__main__":
    rec = Recorder("import")
    base, _ = wall("pass")
    rec.add(stmt="python -c pass", seconds=base)
    for label, stmt in STMTS:
        t, runs = wall(stmt)
        rec.add(stmt=label, seconds=t - base, runs=[r - base for r in runs])
    # what touching Loci pulls in, per submodule (cumulative, -X importtime)
    for mod, s in sorted(importtime_breakdown("from genomeblocks import Loci").items(),
                         key=lambda kv: -kv[1]):
        rec.add(stmt=f"breakdown: {mod}", seconds=s)
    rec.save(baseline=base)
