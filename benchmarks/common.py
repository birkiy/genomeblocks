"""Shared helpers for the genomeblocks benchmark suite.

Every bench script imports this: paths, a small timing harness, result
recording (one JSON file per bench), and environment capture.
"""
from __future__ import annotations

import gc
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("GB_BENCH_DATA", HERE / "data"))
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
for _d in (DATA, RESULTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)


def timeit(fn, *, repeat: int = 5, warmup: int = 1, setup=None) -> dict:
    """Run ``fn`` ``warmup + repeat`` times; return median/min/max seconds.

    ``setup`` (optional) runs before every call, outside the timed region —
    use it to reset caches (e.g. drop a Loci's lazily built index).
    """
    for _ in range(warmup):
        if setup: setup()
        fn()
    runs = []
    for _ in range(repeat):
        if setup: setup()
        gc.collect()
        t0 = time.perf_counter()
        fn()
        runs.append(time.perf_counter() - t0)
    return {"median": statistics.median(runs), "min": min(runs),
            "max": max(runs), "runs": runs}


class Recorder:
    """Collect result rows and write them to ``results/<name>.json``."""

    def __init__(self, name: str):
        self.name = name
        self.rows: list[dict] = []

    def add(self, **row):
        self.rows.append(row)
        t = row.get("seconds")
        extra = {k: v for k, v in row.items() if k not in ("seconds", "runs")}
        print(f"  {extra}  ->  {t:.4f}s" if isinstance(t, float) else f"  {row}",
              flush=True)

    def save(self, **meta):
        out = {"bench": self.name, "env": env_info(), "meta": meta,
               "rows": self.rows}
        p = RESULTS / f"{self.name}.json"
        p.write_text(json.dumps(out, indent=1, default=float))
        print(f"[saved] {p}")


def env_info() -> dict:
    cpu = platform.processor() or ""
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip(); break
    except OSError:
        pass
    try:
        mem_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
    except (ValueError, OSError):
        mem_gb = None
    vers = {}
    for mod in ("numpy", "scipy", "pandas", "pybigtools", "pyBigWig", "cgranges",
                "pyranges", "bioframe", "lightmotif", "Bio", "MOODS",
                "graph_tool", "intervaltree"):
        try:
            m = __import__(mod)
            vers[mod] = getattr(m, "__version__", "installed")
        except Exception:
            pass
    try:
        bt = subprocess.run(["bedtools", "--version"], capture_output=True,
                            text=True).stdout.strip()
        vers["bedtools"] = bt.split()[-1] if bt else None
    except FileNotFoundError:
        pass
    return {"cpu": cpu, "cores": os.cpu_count(), "mem_gb": mem_gb,
            "python": sys.version.split()[0], "platform": platform.platform(),
            "versions": vers}


def read_chromsizes(path=None) -> dict:
    path = path or DATA / "hg38.chrom.sizes"
    out = {}
    for line in open(path):
        c, s = line.split()[:2]
        out[c] = int(s)
    return out
