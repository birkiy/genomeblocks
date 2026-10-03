#!/usr/bin/env python3
"""bigWig signal extraction: ``Loci.signal()`` vs the usual alternatives.

Workload: the heatmap/profile case — a ±3 kb window around each peak centre,
200 bins, mean aggregation, written into a (loci x tracks x bins) cube.
Tracks are genome-wide variable-span bigWigs (~13.8M records, 165 MB each);
the files sit in the OS page cache, as they do in an interactive session.

Parts:
  engines   one track, same loci: every engine, with a correctness check
  scaling   throughput vs number of loci
  parallel  16 tracks: workers=1..4 processes, threads, deepTools -p
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder, timeit

import genomeblocks.signal as gs
from genomeblocks import Loci

FLANK, NBINS = 3_000, 200
BWS = [str(DATA / f"signal_{k}.bw") for k in range(4)]
BIN = Path(sys.executable).parent


def pick_loci(n: int, seed: int = 0) -> Loci:
    """``n`` random A peaks (kept in genome order) whose ±FLANK window fits
    inside the chromosome, so every engine extracts identical windows."""
    from common import read_chromsizes
    cs = read_chromsizes()
    allp = [l for l in Loci.make(str(DATA / "peaks_A_100000.bed"))
            if FLANK <= (l.start + l.end) // 2 <= cs[l.chrom] - FLANK]
    idx = np.sort(np.random.default_rng(seed).choice(len(allp), n, replace=False))
    return Loci([allp[int(i)] for i in idx])


def windows(L):
    for l in L:
        c = (l.start + l.end) // 2
        yield l.chrom, c - FLANK, c + FLANK


# ── baseline engines (same cube layout as signal()) ──────────────────────────

def pybigwig_stats(L, bw, exact):
    import pyBigWig
    cube = np.zeros((len(L), 1, NBINS), np.float32)
    f = pyBigWig.open(bw)
    for i, (c, s, e) in enumerate(windows(L)):
        xs = f.stats(c, s, e, nBins=NBINS, exact=exact)
        cube[i, 0] = [0.0 if x is None else x for x in xs]
    f.close()
    return cube


def pybigwig_values(L, bw):
    import pyBigWig
    cube = np.zeros((len(L), 1, NBINS), np.float32)
    f = pyBigWig.open(bw)
    for i, (c, s, e) in enumerate(windows(L)):
        v = np.asarray(f.values(c, s, e, numpy=True)).reshape(NBINS, -1)
        with np.errstate(all="ignore"):
            cube[i, 0] = np.nan_to_num(np.nanmean(v, axis=1))
    f.close()
    return cube


def pybigtools_values(L, bw):
    import pybigtools
    cube = np.zeros((len(L), 1, NBINS), np.float32)
    f = pybigtools.open(bw)
    for i, (c, s, e) in enumerate(windows(L)):
        v = f.values(c, s, e, missing=np.nan).reshape(NBINS, -1)
        with np.errstate(all="ignore"):
            cube[i, 0] = np.nan_to_num(np.nanmean(v, axis=1))
    f.close()
    return cube


def deeptools(L, bws, procs):
    """computeMatrix reference-point (centre, ±3 kb, 30 bp bins = 200)."""
    with tempfile.TemporaryDirectory() as td:
        bed = Path(td) / "r.bed"
        bed.write_text("".join(f"{l.chrom}\t{l.start}\t{l.end}\n" for l in L))
        cmd = [str(BIN / "computeMatrix"), "reference-point", "--referencePoint",
               "center", "-a", str(FLANK), "-b", str(FLANK), "--binSize", "30",
               "-R", str(bed), "-S", *bws, "-p", str(procs),
               "-o", str(Path(td) / "m.gz"), "--outFileNameMatrix",
               str(Path(td) / "m.tab"), "--missingDataAsZero"]
        subprocess.run(cmd, check=True, capture_output=True)
        m = np.loadtxt(Path(td) / "m.tab", comments="#", skiprows=3, ndmin=2)
    return m.reshape(len(L), len(bws), NBINS).astype(np.float32)


def gb(L, bws, **kw):
    return gs.signal(L, bws, n_bins=NBINS, flank=FLANK, progress=False,
                     verbose=False, **kw)


# ── parts ────────────────────────────────────────────────────────────────────

def compare(ref, cube):
    a, b = ref.ravel().astype(float), cube.ravel().astype(float)
    return {"max_abs_diff": float(np.abs(a - b).max()),
            "pearson_r": float(np.corrcoef(a, b)[0, 1])}


def part_engines(rec):
    print("\n== engines (1 track) ==")
    bw = BWS[0]
    n_main, n_slow = 5_000, 300
    L, Ls = pick_loci(n_main), pick_loci(n_slow, seed=1)
    ref_s = gb(Ls, [bw])                      # reference for correctness
    eng = [
        ("genomeblocks · pybigtools (exact)", lambda X: gb(X, [bw]), L),
        ("genomeblocks · pybigtools (zoom, exact=False)", lambda X: gb(X, [bw], exact=False), L),
        ("genomeblocks · pure-Python reader", lambda X: gb(X, [bw], backend="bigwig"), L),
        ("pybigtools values() + numpy binning", lambda X: pybigtools_values(X, bw), L),
        ("pyBigWig values() + numpy binning", lambda X: pybigwig_values(X, bw), L),
        ("pyBigWig stats(nBins) exact", lambda X: pybigwig_stats(X, bw, True), Ls),
        ("pyBigWig stats(nBins) zoom", lambda X: pybigwig_stats(X, bw, False), Ls),
        ("deepTools computeMatrix (-p 1)", lambda X: deeptools(X, [bw], 1), L),
    ]
    for name, fn, X in eng:
        rep = 3 if X is L else 1
        t = timeit(lambda: fn(X), repeat=rep, warmup=1)
        chk = compare(ref_s, fn(Ls))
        rec.add(part="engines", engine=name, n_loci=len(X), n_tracks=1,
                seconds=t["median"], runs=t["runs"],
                rate=len(X) / t["median"], **chk)


def part_scaling(rec):
    print("\n== scaling with #loci (1 track) ==")
    for n in (1_000, 5_000, 20_000, 50_000):
        L = pick_loci(n)
        for name, kw in (("genomeblocks · pybigtools (exact)", {}),
                         ("genomeblocks · pybigtools (zoom, exact=False)", {"exact": False}),
                         ("genomeblocks · pure-Python reader", {"backend": "bigwig"})):
            t = timeit(lambda: gb(L, [BWS[0]], **kw), repeat=3 if n <= 20_000 else 1)
            rec.add(part="scaling", engine=name, n_loci=n, n_tracks=1,
                    seconds=t["median"], runs=t["runs"], rate=n / t["median"])
        t = timeit(lambda: pybigwig_values(L, BWS[0]), repeat=1)
        rec.add(part="scaling", engine="pyBigWig values() + numpy binning", n_loci=n,
                n_tracks=1, seconds=t["median"], runs=t["runs"], rate=n / t["median"])


def threaded(L, bws, workers):
    """signal()'s chunking, but on a ThreadPoolExecutor (one shared cube)."""
    from concurrent.futures import ThreadPoolExecutor
    n, T = len(L), len(bws)
    chroms = [l.chrom for l in L]
    st = np.fromiter((l.start for l in L), np.int64, n)
    en = np.fromiter((l.end for l in L), np.int64, n)
    cube = np.zeros((n, T, NBINS), np.float32)
    t_ranges = gs._even_ranges(T, min(workers, T))

    def run(t_lo, t_hi):
        hs = [gs._open_pybigtools(p) for p in bws[t_lo:t_hi]]
        gs._extract_chunk(cube, hs, hs[0].chroms(), chroms, st, en, t_lo, 0, n,
                          NBINS, FLANK, "mean", False, True)
        for h in hs:
            h.close()

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(lambda r: run(*r), t_ranges))
    return cube


def part_parallel(rec):
    print("\n== parallel (16 tracks) ==")
    T, n = 16, 5_000
    bws = [BWS[i % 4] for i in range(T)]
    L = pick_loci(n)
    ref = None
    real_cpu = gs.cpu_count
    for w in (1, 2, 3, 4):
        # signal() caps workers at cpu_count()//2; lift the cap past 2 so the
        # sweep can use all 4 cores of this machine
        gs.cpu_count = (lambda: 8) if w > 2 else real_cpu
        res = {}
        t = timeit(lambda: res.__setitem__("c", gb(L, bws, workers=w)), repeat=3)
        ref = res["c"] if ref is None else ref
        rec.add(part="parallel", engine="genomeblocks · processes" +
                (" (cap lifted)" if w > 2 else ""), workers=w, n_loci=n, n_tracks=T,
                seconds=t["median"], runs=t["runs"], rate=n * T / t["median"],
                identical=bool(np.array_equal(ref, res["c"])))
    gs.cpu_count = real_cpu
    for w in (2, 4):
        t = timeit(lambda: threaded(L, bws, w), repeat=3)
        rec.add(part="parallel", engine="threads (same chunks)", workers=w,
                n_loci=n, n_tracks=T, seconds=t["median"], runs=t["runs"],
                rate=n * T / t["median"])
    for p in (1, 4):
        t = timeit(lambda: deeptools(L, bws, p), repeat=1, warmup=0)
        rec.add(part="parallel", engine="deepTools computeMatrix", workers=p,
                n_loci=n, n_tracks=T, seconds=t["median"], runs=t["runs"],
                rate=n * T / t["median"])


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    parts = sys.argv[1:] or ["engines", "scaling", "parallel"]
    rec = Recorder("signal" if len(parts) == 3 else "signal_" + "_".join(parts))
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save(flank=FLANK, n_bins=NBINS, tracks=BWS)
