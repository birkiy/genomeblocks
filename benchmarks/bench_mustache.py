#!/usr/bin/env python3
"""Hi-C loops: Mustache 1.3.3 against genomeblocks' mustache(), on the same balanced cooler.

  synthetic  a 5 kb cooler of chr21 + chr22 built here (into $GB_BENCH_DATA/mustache_syn_5kb.cool, reused):
             Poisson contacts with a power-law distance decay times per-bin visibility, plus a 3x3 blob at each
             loop of data/loops.bedpe on those chromosomes (5 bins to 2 Mb apart), balanced with cooler's ICE.
             Mustache runs when MUSTACHE is set (its command, e.g. the ``mustache`` of an env holding
             mustache-hic; MUSTACHE_EXEC is an optional command prefix such as "singularity exec img.sif"),
             with ``-p`` = GB_BENCH_THREADS (default 4), the threads genomeblocks gets. Without it only
             genomeblocks is timed.
  real       GB_BENCH_MUSTACHE_REAL=samples.tsv, one row per sample with the columns name, cool (a cooler URI
             with a weight column), mustache_tsv (Mustache's output on that cooler with its defaults),
             mustache_seconds (its wall time) and an optional note. genomeblocks is timed and compared with
             that output; no path is stored in the results.

Compared: the loop sets (chromosome, bin1, bin2), the largest relative FDR difference and the DETECTION_SCALE
of the shared loops, and Mustache's output file line for line (sorted: Mustache's row order across blocks
follows its processes). On the synthetic cooler, also the planted loops recovered (within one bin).
"""
from __future__ import annotations

import json
import os
import resource
import shlex
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder, read_chromsizes

COOL = DATA / "mustache_syn_5kb.cool"
LOOPS = DATA / "loops.bedpe"
MUSTACHE = os.environ.get("MUSTACHE")
MUSTACHE_EXEC = shlex.split(os.environ.get("MUSTACHE_EXEC", ""))
THREADS = int(os.environ.get("GB_BENCH_THREADS", "4"))
REAL = os.environ.get("GB_BENCH_MUSTACHE_REAL")
CHROMS = ("chr21", "chr22")
RES = 5000


def make_cool(path: Path):
    """The synthetic cooler (see the module notes); returns the planted loops as (chrom, bin1, bin2)."""
    import cooler
    import pandas as pd
    rng = np.random.default_rng(21)
    cs = read_chromsizes()
    lp = pd.read_csv(LOOPS, sep="\t", header=None, usecols=[0, 1, 2, 4, 5, 7])
    bins, pix, planted, off = [], [], [], 0
    for c in CHROMS:
        n = -(-cs[c] // RES)
        st = np.arange(n) * RES
        bins.append(pd.DataFrame({"chrom": c, "start": st, "end": np.minimum(st + RES, cs[c])}))
        vis = rng.lognormal(0, 0.25, n)
        i, d = np.meshgrid(np.arange(n), np.arange(401), indexing="ij")
        i, j = i.ravel(), (i + d).ravel()
        keep = j < n
        i, j = i[keep], j[keep]
        lam = 60.0 * (j - i + 1.0) ** -0.75 * vis[i] * vis[j]
        L = lp[lp[0] == c]
        a = ((L[1] + L[2]) // 2 // RES).to_numpy()
        b = ((L[4] + L[5]) // 2 // RES).to_numpy()
        ok = (b - a >= 5) & (b - a <= 400)
        a, b, score = a[ok], b[ok], L[7].to_numpy()[ok]
        planted += [(c, int(x), int(y)) for x, y in zip(a, b)]
        boost = np.zeros(len(i))
        key = i * n + j
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                k = np.searchsorted(key, (a + dx) * n + (b + dy))
                k = np.minimum(k, len(key) - 1)
                hit = key[k] == (a + dx) * n + (b + dy)
                np.add.at(boost, k[hit], 0.5 * score[hit] * np.exp(-(dx * dx + dy * dy) / 2.0))
        cnt = rng.poisson(lam * (1 + boost))
        m = cnt > 0
        pix.append(pd.DataFrame({"bin1_id": i[m] + off, "bin2_id": j[m] + off, "count": cnt[m]}))
        off += n
    cooler.create_cooler(str(path), pd.concat(bins, ignore_index=True), pd.concat(pix, ignore_index=True),
                         ordered=True)
    cooler.balance_cooler(cooler.Cooler(str(path)), store=True)
    with open(path.with_suffix(".planted.json"), "w") as f:
        json.dump(planted, f)
    return planted


def time_gb(cool, repeat=3):
    """Median wall time and peak memory (the process's own VmHWM) of mustache() in a fresh process, after a
    warm-up run."""
    code = f"""
import json, statistics, time
from genomeblocks._mustache import mustache
mustache({cool!r}, threads={THREADS}, verbose=False)
ts = []
for _ in range({repeat}):
    t = time.perf_counter(); mustache({cool!r}, threads={THREADS}, verbose=False); ts.append(time.perf_counter() - t)
hwm = next(int(l.split()[1]) for l in open("/proc/self/status") if l.startswith("VmHWM"))
print(json.dumps({{"seconds": statistics.median(ts), "runs": ts, "peak_mb": hwm / 1024}}))
"""
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def run_mustache(cool, out: Path):
    """``mustache -f cool -r 5000 -o out -p THREADS``: wall time and peak memory."""
    t = time.perf_counter()
    with open(out.with_suffix(".log"), "w") as log:
        subprocess.run(MUSTACHE_EXEC + [MUSTACHE, "-f", str(cool), "-r", str(RES), "-o", str(out), "-p", str(THREADS)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    return {"seconds": time.perf_counter() - t,
            "peak_mb": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024}


def compare(res, ref_tsv, tmp: Path, planted=None):
    """Agreement of a mustache() result with a Mustache output file."""
    import pandas as pd
    mine = res.write(str(tmp), "gb")["loops"]
    key = ["BIN1_CHR", "BIN1_START", "BIN2_START"]
    R = pd.read_csv(ref_tsv, sep="\t")
    G = pd.read_csv(mine, sep="\t")
    m = R.merge(G, on=key, how="outer", suffixes=("_ref", "_gb"), indicator=True)
    sh = m[m["_merge"] == "both"]
    rel = ((sh.FDR_gb - sh.FDR_ref).abs() / sh.FDR_ref).max() if len(sh) else 0.0
    a = sorted(open(ref_tsv).read().splitlines()[1:])
    b = sorted(open(mine).read().splitlines()[1:])
    row = {"loops_ref": len(R), "loops_gb": len(G), "loops_shared": len(sh),
           "max_rel_FDR": float(rel), "FDR_identical": int((sh.FDR_gb == sh.FDR_ref).sum()),
           "scales_equal": bool((sh.DETECTION_SCALE_gb == sh.DETECTION_SCALE_ref).all()),
           "lines_identical": a == b, "lines_differing": len(set(a) ^ set(b))}
    if planted is not None:
        near = {(c, x + dx, y + dy) for c, x, y in planted for dx in (-1, 0, 1) for dy in (-1, 0, 1)}
        for name, T in (("ref", R), ("gb", G)):
            got = set(zip(T.BIN1_CHR, T.BIN1_START // RES, T.BIN2_START // RES))
            row[f"planted_found_{name}"] = sum(any((c, x + dx, y + dy) in got for dx in (-1, 0, 1)
                                                   for dy in (-1, 0, 1)) for c, x, y in planted)
            row[f"calls_at_planted_{name}"] = len(got & near)
        row["planted"] = len(planted)
    return row


if __name__ == "__main__":
    from genomeblocks._mustache import mustache
    rec = Recorder("mustache")
    parts = os.environ.get("GB_BENCH_PARTS", "synthetic,real").split(",")
    if "synthetic" in parts:
        if COOL.exists():
            planted = [tuple(p) for p in json.load(open(COOL.with_suffix(".planted.json")))]
        else:
            planted = make_cool(COOL)
        ref = None
        with tempfile.TemporaryDirectory(dir=DATA) as d:
            if MUSTACHE:                              # first, so the children's peak memory is Mustache's
                ref_tsv = Path(d) / "mustache.tsv"
                ref = run_mustache(COOL, ref_tsv)
            else:
                print("[skip] Mustache: set MUSTACHE (and MUSTACHE_EXEC) to compare")
            rec.add(part="synthetic", impl="genomeblocks", threads=THREADS, **time_gb(str(COOL)))
            if ref:
                res = mustache(str(COOL), threads=THREADS)
                rec.add(part="synthetic", impl="Mustache", threads=THREADS, **ref)
                rec.add(part="synthetic", impl="agreement", **compare(res, ref_tsv, Path(d), planted))
    if "real" in parts and REAL:
        import pandas as pd
        for s in pd.read_table(REAL).fillna("").itertuples():
            t = time_gb(s.cool, repeat=1)
            res = mustache(s.cool, threads=THREADS)
            with tempfile.TemporaryDirectory(dir=DATA) as d:
                agree = compare(res, s.mustache_tsv, Path(d))
            rec.add(part="real", sample=s.name, threads=THREADS, seconds=t["seconds"], peak_mb=t["peak_mb"],
                    mustache_seconds=float(s.mustache_seconds), note=str(getattr(s, "note", "")), **agree)
    rec.save(mustache="1.3.3 (github.com/ay-lab/mustache df71aa1)", threads=THREADS)
