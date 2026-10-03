#!/usr/bin/env python3
"""Atlas (GIGGLE-style) enrichment: one sparse mat-vec vs per-track loops.

Index: 500 peak files (median 12k peaks; 7.4M peaks total) tiled at 1 kb into
a single bin x track CSR matrix. Query: 20k peaks.

Parts:
  build     Atlas.make: time, nnz, resident bytes vs BED bytes
  query     search() time vs number of tracks, against
              - GIGGLE (C, on-disk index; part "giggle")
              - per-track cgranges loop (Loci indexes prebuilt, reused)
              - pyranges overlap() per track (vectorised join)
              - bedtools intersect -C (all files in one call)
  accuracy  bin-level counts vs exact interval counts (Spearman), by bin size
  bootstrap Atlas.bootstrap(n=100) per-iteration cost
"""
from __future__ import annotations

import glob
import os
import subprocess
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder, timeit

from genomeblocks import Atlas, Loci

TRACKS = sorted(glob.glob(str(DATA / "atlas" / "track_*.bed")))
QUERY = str(DATA / "atlas_query.bed")
CS = str(DATA / "hg38.chrom.sizes")
T_SWEEP = [10, 50, 100, 250, 500]


def sub_atlas(a: Atlas, T: int) -> Atlas:
    return Atlas(bin_size=a.bin_size, chrom_names=a.chrom_names,
                 chrom_sizes=a.chrom_sizes, chrom_offsets=a.chrom_offsets,
                 n_bins=a.n_bins, track_names=a.track_names[:T],
                 track_n_peaks=a.track_n_peaks[:T], track_n_bins=a.track_n_bins[:T],
                 M=a.M[:, :T].tocsr())


def exact_counts(Q: Loci, tracks: list) -> np.ndarray:
    """#query intervals overlapping each track (cgranges, prebuilt indexes)."""
    out = np.zeros(len(tracks), np.int64)
    for j, t in enumerate(tracks):
        cg = t.cgr
        out[j] = sum(1 for q in Q if any(True for _ in cg.overlap(q.chrom, q.start, q.end)))
    return out


def part_build(rec):
    print("\n== build ==")
    bed_bytes = sum(__import__("os").path.getsize(p) for p in TRACKS)
    n_peaks = sum(sum(1 for _ in open(p)) for p in TRACKS)
    for bs in (200, 1000, 5000):
        for w in (1, 3):
            if bs != 1000 and w == 1:
                continue
            t0 = time.perf_counter()
            a = Atlas.make(TRACKS, chromsizes=CS, bin_size=bs, workers=w, verbose=False)
            dt = time.perf_counter() - t0
            mem = a.M.data.nbytes + a.M.indices.nbytes + a.M.indptr.nbytes
            rec.add(part="build", bin_size=bs, workers=w, seconds=dt,
                    nnz=int(a.M.nnz), index_bytes=int(mem), bed_bytes=int(bed_bytes),
                    n_peaks=int(n_peaks), n_tracks=len(TRACKS))
            if bs == 1000 and w == 3:
                a.save(str(DATA / "atlas_1kb.npz"))


def part_query(rec):
    import pyranges as pr
    print("\n== query vs #tracks ==")
    A = Atlas.load(str(DATA / "atlas_1kb.npz"))
    Q = Loci.make(QUERY)
    loci_tracks = [Loci.make(p) for p in TRACKS]
    for t in loci_tracks:
        t.cgr                                   # prebuild every index once
    gq = pr.read_bed(QUERY)
    for T in T_SWEEP:
        a = sub_atlas(A, T)
        tq = timeit(lambda: a.search(Q), repeat=5)
        rec.add(part="query", engine="Atlas.search (incl. Fisher + DataFrame)",
                n_tracks=T, seconds=tq["median"], runs=tq["runs"])
        qb = a._bins_union(a._intervals_to_bin_ranges(Q))
        tq = timeit(lambda: a._overlap_counts(qb), repeat=5)
        rec.add(part="query", engine="Atlas sparse row-sum only", n_tracks=T,
                seconds=tq["median"], runs=tq["runs"])
        rep = 3 if T <= 100 else 1
        tq = timeit(lambda: exact_counts(Q, loci_tracks[:T]), repeat=rep,
                    warmup=1 if T <= 100 else 0)
        rec.add(part="query", engine="per-track cgranges loop (prebuilt)",
                n_tracks=T, seconds=tq["median"], runs=tq["runs"])
        # pyranges.count_overlaps fails under pandas 3 (pyranges 0.1.4), so
        # count per track with the vectorised overlap() join instead
        grs = [pr.read_bed(p) for p in TRACKS[:T]]
        tq = timeit(lambda: [len(gq.overlap(g)) for g in grs], repeat=rep,
                    warmup=1 if T <= 100 else 0)
        rec.add(part="query", engine="pyranges overlap(), per track", n_tracks=T,
                seconds=tq["median"], runs=tq["runs"])
        if T <= 250:
            def _bt():
                p = subprocess.run(["bedtools", "intersect", "-C", "-a", QUERY,
                                    "-b", *TRACKS[:T]], capture_output=True)
                return len(p.stdout)
            tq = timeit(_bt, repeat=1, warmup=0)
            rec.add(part="query", engine="bedtools intersect -C (CLI)",
                    n_tracks=T, seconds=tq["median"], runs=tq["runs"])


def part_accuracy(rec):
    from scipy.stats import spearmanr
    print("\n== accuracy by bin size ==")
    Q = Loci.make(QUERY)
    loci_tracks = [Loci.make(p) for p in TRACKS]
    exact = exact_counts(Q, loci_tracks)
    for bs in (200, 1000, 5000):
        a = Atlas.make(TRACKS, chromsizes=CS, bin_size=bs, workers=3, verbose=False)
        df = a.search(Q)
        by = dict(zip(df["name"], df["overlaps"]))
        ov = np.array([by[n] for n in a.track_names])
        rho = spearmanr(ov, exact).statistic
        # do the enrichment rankings agree with an exact interval-level test?
        top_exact = set(np.argsort(-exact / a.track_n_peaks)[:25])
        top_atlas = set(np.argsort(-ov / a.track_n_bins)[:25])
        tq = timeit(lambda: a.search(Q), repeat=5)
        rec.add(part="accuracy", bin_size=bs, spearman_counts=float(rho),
                top25_overlap=len(top_exact & top_atlas), seconds=tq["median"],
                runs=tq["runs"], nnz=int(a.M.nnz))


GIGGLE = os.environ.get("GIGGLE", "giggle")


def _giggle_inputs():
    """bgzipped copies of the tracks + query (GIGGLE reads .bed.gz)."""
    gz = DATA / "atlas_gz"
    gz.mkdir(exist_ok=True)
    for p in TRACKS:
        out = gz / (Path(p).name + ".gz")
        if not out.exists():
            with open(out, "wb") as fh:
                subprocess.run(["bgzip", "-c", p], stdout=fh, check=True)
    q = DATA / "atlas_query.bed.gz"
    if not q.exists():
        with open(q, "wb") as fh:
            subprocess.run(["bgzip", "-c", QUERY], stdout=fh, check=True)
    return sorted(gz.glob("track_*.bed.gz")), q


def _giggle_search(idx, q, cwd):
    out = subprocess.run([GIGGLE, "search", "-i", idx, "-q", str(q), "-s"], cwd=cwd,
                         capture_output=True, text=True, check=True).stdout
    rows = {}
    for line in out.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        rows[Path(f[0]).name.replace(".bed.gz", "")] = (int(f[2]), float(f[7]))
    return rows


def part_giggle(rec):
    """GIGGLE (Layer et al. 2018): build + search at each collection size."""
    from scipy.stats import spearmanr
    print("\n== GIGGLE ==")
    gz, q = _giggle_inputs()
    base = DATA / "giggle"
    base.mkdir(exist_ok=True)
    for T in T_SWEEP:
        tdir, idir = base / f"t{T}", f"i{T}"
        tdir.mkdir(exist_ok=True)
        for p in gz[:T]:
            link = tdir / p.name
            if not link.exists():
                link.symlink_to(p)
        build = ["bash", "-c", f'{GIGGLE} index -i "t{T}/*.gz" -o {idir} -f -s > /dev/null']
        tb = timeit(lambda: subprocess.run(build, cwd=base, check=True), repeat=1, warmup=0)
        size = sum(f.stat().st_size for f in (base / idir).rglob("*") if f.is_file())
        rec.add(part="giggle", kind="build", n_tracks=T, seconds=tb["median"], runs=tb["runs"],
                index_bytes=int(size))
        ts = timeit(lambda: _giggle_search(idir, q, base), repeat=5)
        rec.add(part="giggle", kind="query", engine="GIGGLE search -s (CLI)", n_tracks=T,
                seconds=ts["median"], runs=ts["runs"])
    # agreement at full size: GIGGLE vs exact interval counts and vs Atlas
    res = _giggle_search(f"i{T_SWEEP[-1]}", q, base)
    names = [Path(p).name.replace(".bed", "") for p in TRACKS]
    g_ov = np.array([res[n][0] for n in names])
    g_sc = np.array([res[n][1] for n in names])
    Q = Loci.make(QUERY)
    exact = exact_counts(Q, [Loci.make(p) for p in TRACKS])
    A = Atlas.load(str(DATA / "atlas_1kb.npz"))
    df = A.search(Q)
    a_sc = df.set_index("name").loc[A.track_names, "giggle_score"].to_numpy()
    top_g = set(np.argsort(-g_sc)[:25]); top_a = set(np.argsort(-a_sc)[:25])
    rec.add(part="giggle", kind="agreement", n_tracks=len(TRACKS),
            seconds=0.0, spearman_overlaps_vs_exact=float(spearmanr(g_ov, exact).statistic),
            spearman_score_vs_atlas=float(spearmanr(g_sc, a_sc).statistic),
            top25_shared_with_atlas=len(top_g & top_a))


def part_load(rec):
    """Cold, one-shot use (like a CLI call): load the saved index, then search.

    GIGGLE's search is a CLI that opens its on-disk index each call, so this
    is the like-for-like comparison to it; ``query`` times a loaded Atlas.
    """
    import sys
    print("\n== cold load + search ==")
    npz = str(DATA / "atlas_1kb.npz")
    t = timeit(lambda: Atlas.load(npz), repeat=3)
    rec.add(part="load", kind="Atlas.load (npz)", n_tracks=len(TRACKS), seconds=t["median"],
            runs=t["runs"])
    code = ("from genomeblocks import Atlas, Loci; "
            f"Atlas.load({npz!r}).search(Loci.make({QUERY!r}))")
    t = timeit(lambda: subprocess.run([sys.executable, "-W", "ignore", "-c", code], check=True),
               repeat=3)
    rec.add(part="load", kind="fresh python: import + load + search", n_tracks=len(TRACKS),
            seconds=t["median"], runs=t["runs"])


def part_bootstrap(rec):
    print("\n== bootstrap ==")
    A = Atlas.load(str(DATA / "atlas_1kb.npz"))
    Q = Loci.make(QUERY)
    for n in (10, 100):
        t = timeit(lambda: A.bootstrap(Q, n=n, seed=0, verbose=False), repeat=3)
        rec.add(part="bootstrap", n_iter=n, n_tracks=len(A), seconds=t["median"],
                runs=t["runs"], per_iter=t["median"] / n)


if __name__ == "__main__":
    import sys
    parts = sys.argv[1:] or ["build", "query", "giggle", "load", "accuracy", "bootstrap"]
    rec = Recorder("atlas")
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save(n_tracks=len(TRACKS), query=QUERY)
