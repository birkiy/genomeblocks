#!/usr/bin/env python3
"""Atlas (GIGGLE-style) enrichment: one sparse mat-vec vs per-track loops.

Index: 500 peak files (median 12k peaks; 7.4M peaks total) tiled at 1 kb into
a single bin x track CSR matrix. Query: 20k peaks.

Parts:
  build     Atlas.make: time, nnz, resident bytes vs BED bytes
  query     search() time vs number of tracks, against
              - per-track cgranges loop (Loci indexes prebuilt, reused)
              - pyranges.count_overlaps (all tracks in one call)
              - bedtools intersect -C (all files in one call)
  accuracy  bin-level counts vs exact interval counts (Spearman), by bin size
  bootstrap Atlas.bootstrap(n=100) per-iteration cost
"""
from __future__ import annotations

import glob
import subprocess
import time

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
        grs = {i: pr.read_bed(p) for i, p in enumerate(TRACKS[:T])}
        tq = timeit(lambda: pr.count_overlaps(grs, gq), repeat=rep,
                    warmup=1 if T <= 100 else 0)
        rec.add(part="query", engine="pyranges.count_overlaps", n_tracks=T,
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
    parts = sys.argv[1:] or ["build", "query", "accuracy", "bootstrap"]
    rec = Recorder("atlas")
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save(n_tracks=len(TRACKS), query=QUERY)
