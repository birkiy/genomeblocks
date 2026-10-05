#!/usr/bin/env python3
"""bigWig signal extraction: ``Loci.signal()`` vs the usual alternatives.

Workload: the heatmap/profile case — a ±3 kb window around each peak centre,
200 bins, mean aggregation, written into a (loci x tracks x bins) cube.
Tracks are genome-wide variable-span bigWigs (~13.8M records, 165 MB each);
the files sit in the OS page cache, as they do in an interactive session.

Parts:
  engines   one track, same loci: every engine, with a correctness check
  scaling   throughput vs number of loci
  parallel  16 tracks: workers=1..8 processes, threads, deepTools -p
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
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
    allp = Loci.make(str(DATA / "peaks_A_100000.bed"))
    size = np.array([cs[c] for c in allp.chroms], np.int64)
    c = allp.centers
    allp = allp.take(np.flatnonzero((c >= FLANK) & (c <= size - FLANK)))
    idx = np.sort(np.random.default_rng(seed).choice(len(allp), n, replace=False))
    return allp.take(idx)


def windows(L):
    c = L.centers
    yield from zip(L.chroms.tolist(), (c - FLANK).tolist(), (c + FLANK).tolist())


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
        ("genomeblocks · pure-Python reader", lambda X: gb(X, [bw], backend="python"), L),
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
    from genomeblocks.backends.bigwig import open_bigwig
    n, T = len(L), len(bws)
    chroms = L.chroms
    st, en = L.starts, L.ends
    cube = np.zeros((n, T, NBINS), np.float32)
    t_ranges = gs._even_ranges(T, min(workers, T))

    def run(t_lo, t_hi):
        hs = [open_bigwig(p, backend="pybigtools") for p in bws[t_lo:t_hi]]
        gs._extract(cube, hs, chroms, st, en, t_lo, 0, n, NBINS, FLANK, "mean", False, True)
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
    for w in (1, 2, 4, 8):
        res = {}
        t = timeit(lambda: res.__setitem__("c", gb(L, bws, workers=w)), repeat=3)
        ref = res["c"] if ref is None else ref
        rec.add(part="parallel", engine="genomeblocks · processes", workers=w, n_loci=n,
                n_tracks=T, seconds=t["median"], runs=t["runs"], rate=n * T / t["median"],
                identical=bool(np.array_equal(ref, res["c"])))
    for w in (2, 4, 8):
        t = timeit(lambda: threaded(L, bws, w), repeat=3)
        rec.add(part="parallel", engine="threads (same chunks)", workers=w,
                n_loci=n, n_tracks=T, seconds=t["median"], runs=t["runs"],
                rate=n * T / t["median"])
    for p in (1, 4, 8):
        t = timeit(lambda: deeptools(L, bws, p), repeat=1, warmup=0)
        rec.add(part="parallel", engine="deepTools computeMatrix", workers=p,
                n_loci=n, n_tracks=T, seconds=t["median"], runs=t["runs"],
                rate=n * T / t["median"])


def retime_deeptools(workers: int = 4):
    """Re-time one deepTools point in place (e.g. after it ran on a busy machine)."""
    import json
    from common import RESULTS
    p = RESULTS / "signal.json"
    d = json.loads(p.read_text())
    T, n = 16, 5_000
    bws = [BWS[i % 4] for i in range(T)]
    L = pick_loci(n)
    t = timeit(lambda: deeptools(L, bws, workers), repeat=1, warmup=0)
    for r in d["rows"]:
        if r["part"] == "parallel" and r["engine"].startswith("deepTools") and r["workers"] == workers:
            r.update(seconds=t["median"], runs=t["runs"], rate=n * T / t["median"], retimed=True)
    p.write_text(json.dumps(d, indent=1))
    print(f"[retimed] deepTools -p {workers}: {t['median']:.1f}s")


if __name__ == "__main__":
    if sys.argv[1:] == ["retime_deeptools"]:
        import warnings
        warnings.filterwarnings("ignore")
        retime_deeptools()
        sys.exit()
    import warnings
    warnings.filterwarnings("ignore")
    parts = sys.argv[1:] or ["engines", "scaling", "parallel"]
    rec = Recorder("signal" if len(parts) == 3 else "signal_" + "_".join(parts))
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save(flank=FLANK, n_bins=NBINS, tracks=BWS)
