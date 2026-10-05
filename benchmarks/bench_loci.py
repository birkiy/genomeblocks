#!/usr/bin/env python3
"""Loci interval algebra: genomeblocks vs the usual Python/CLI tools.

What is measured (A and B are n-peak sets, n = 1k .. 1M):

  intersect   A & B  -> loci of A overlapping any locus of B   (bedtools -u)
  merge       A.merge()                                         (bedtools merge)
  make        Loci.make(bed)  BED parsing into Locus objects
  latency     one overlap query against an indexed set (interactive use)
  caching     1st call (index build + query) vs. later calls (query only)

Engines: genomeblocks (its numpy kernel, the pip default; cgranges when
installed), pyranges, bioframe, bedtools (CLI, includes file I/O + process
start), intervaltree, and a naive per-chromosome double loop. Every backend
of the intervals family against each other is bench_backends.py.
"""
from __future__ import annotations

import subprocess

import numpy as np
import pandas as pd

from common import DATA, Recorder, timeit

import genomeblocks as gb
from genomeblocks import Loci

SIZES = [1_000, 10_000, 100_000, 1_000_000]
HAS_CGRANGES = gb.backends.installed("intervals", "cgranges")


# ── engines ──────────────────────────────────────────────────────────────────

def fresh(L: Loci) -> Loci:
    """Same loci, no cached index."""
    out = L.copy()
    out._dirty()
    return out


def naive_intersect(A, B):
    by = {}
    for b in B:
        by.setdefault(b.chrom, []).append((b.start, b.end))
    out = []
    for a in A:
        for s, e in by.get(a.chrom, ()):
            if s < a.end and e > a.start:
                out.append(a); break
    return out


def intervaltree_intersect(A, B):
    from intervaltree import IntervalTree
    trees = {}
    for b in B:
        trees.setdefault(b.chrom, IntervalTree()).addi(b.start, b.end)
    return [a for a in A if a.chrom in trees and trees[a.chrom].overlaps(a.start, a.end)]


def to_df(L):
    return L.to_bioframe()[["chrom", "start", "end"]]


# ── benches ──────────────────────────────────────────────────────────────────

def bench_intersect(rec: Recorder):
    import pyranges as pr
    import bioframe as bf
    print("\n== intersect ==")
    for n in SIZES:
        pa, pb = DATA / f"peaks_A_{n}.bed", DATA / f"peaks_B_{n}.bed"
        A, B = Loci.make(str(pa)), Loci.make(str(pb))
        rep = 5 if n <= 100_000 else 3
        truth = None

        def add(engine, t, n_out):
            nonlocal truth
            if truth is None:
                truth = n_out
            rec.add(op="intersect", engine=engine, n=n, seconds=t["median"],
                    runs=t["runs"], n_out=n_out, agrees=(n_out == truth))

        # genomeblocks: the numpy kernel (pip default), then cgranges when installed
        res = {}
        t = timeit(lambda: res.__setitem__("o", A & B), repeat=rep)
        add("genomeblocks", t, len(res["o"]))
        if HAS_CGRANGES:
            Bc = fresh(B)
            t = timeit(lambda: res.__setitem__("o", A.intersect(Bc, backend="cgranges")), repeat=rep,
                       setup=Bc._dirty)
            add("genomeblocks (cgranges)", t, len(res["o"]))
            t = timeit(lambda: res.__setitem__("o", A.intersect(Bc, backend="cgranges")), repeat=rep)
            add("genomeblocks (cgranges, warm index)", t, len(res["o"]))

        # pyranges / bioframe on their native frames
        ga = pr.read_bed(str(pa)); gb_ = pr.read_bed(str(pb))
        t = timeit(lambda: res.__setitem__("o", ga.overlap(gb_)), repeat=rep)
        add("pyranges", t, len(res["o"]))
        da, db = to_df(A), to_df(B)

        def _bf():
            ov = bf.overlap(da, db, how="inner", return_index=True)
            return da.loc[np.unique(ov["index"].to_numpy())]
        t = timeit(lambda: res.__setitem__("o", _bf()), repeat=rep)
        add("bioframe", t, len(res["o"]))

        # bedtools (process + file I/O + text output)
        def _bt():
            out = subprocess.run(["bedtools", "intersect", "-u", "-a", str(pa),
                                  "-b", str(pb)], capture_output=True, text=True).stdout
            return out.count("\n")
        t = timeit(lambda: res.__setitem__("o", _bt()), repeat=rep)
        add("bedtools (CLI)", t, res["o"])

        if n <= 100_000:
            t = timeit(lambda: res.__setitem__("o", intervaltree_intersect(A, B)),
                       repeat=1 if n == 100_000 else 3, warmup=0 if n == 100_000 else 1)
            add("intervaltree", t, len(res["o"]))
        if n <= 10_000:
            t = timeit(lambda: res.__setitem__("o", naive_intersect(A, B)),
                       repeat=1, warmup=0)
            add("naive double loop", t, len(res["o"]))


def bench_merge(rec: Recorder):
    import pyranges as pr
    import bioframe as bf
    print("\n== merge (unsorted union A+B; sort + merge) ==")
    for n in SIZES:
        pa = DATA / f"peaks_AB_{n}.bed"
        A = Loci.make(str(pa))
        rep = 5 if n <= 100_000 else 3
        res = {}
        t = timeit(lambda: res.__setitem__("o", A.merge()), repeat=rep)
        rec.add(op="merge", engine="genomeblocks", n=n, seconds=t["median"],
                runs=t["runs"], n_out=len(res["o"]))
        ga = pr.read_bed(str(pa))
        t = timeit(lambda: res.__setitem__("o", ga.merge()), repeat=rep)
        rec.add(op="merge", engine="pyranges", n=n, seconds=t["median"],
                runs=t["runs"], n_out=len(res["o"]))
        da = to_df(A)
        t = timeit(lambda: res.__setitem__("o", bf.merge(da, min_dist=0)), repeat=rep)
        rec.add(op="merge", engine="bioframe", n=n, seconds=t["median"],
                runs=t["runs"], n_out=len(res["o"]))

        def _bt():  # bedtools merge needs sorted input
            return subprocess.run(f"bedtools sort -i {pa} | bedtools merge -i -",
                                  shell=True, capture_output=True,
                                  text=True).stdout.count("\n")
        t = timeit(lambda: res.__setitem__("o", _bt()), repeat=rep)
        rec.add(op="merge", engine="bedtools (CLI)", n=n, seconds=t["median"],
                runs=t["runs"], n_out=res["o"])


def bench_make(rec: Recorder):
    import pyranges as pr
    print("\n== make (BED parse) ==")
    for n in SIZES:
        pa = str(DATA / f"peaks_A_{n}.bed")
        rep = 5 if n <= 100_000 else 3
        t = timeit(lambda: Loci.make(pa), repeat=rep)
        rec.add(op="make", engine="genomeblocks", n=n, seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: pd.read_csv(pa, sep="\t", header=None, usecols=[0, 1, 2]),
                   repeat=rep)
        rec.add(op="make", engine="pandas.read_csv", n=n, seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: pr.read_bed(pa), repeat=rep)
        rec.add(op="make", engine="pyranges", n=n, seconds=t["median"], runs=t["runs"])


def bench_latency(rec: Recorder):
    """Per-query latency of one overlap lookup against an indexed 100k set."""
    import pyranges as pr
    import bioframe as bf
    from intervaltree import IntervalTree
    print("\n== single-query latency ==")
    n = 100_000
    B = Loci.make(str(DATA / f"peaks_B_{n}.bed"))
    Q = Loci.make(str(DATA / "peaks_A_1000.bed"))[:200]
    for kind, label in (("genomeblocks", "genomeblocks (numpy point index)"),
                        ("cgranges", "genomeblocks (cgranges)"), ("ncls", "genomeblocks (ncls)")):
        if not gb.backends.installed("intervals", kind):
            continue
        Bk = fresh(B)
        Bk.overlap_rows("chr1", 0, 1, backend=kind)                 # build the index once
        t = timeit(lambda: [Bk.overlap_rows(q, backend=kind) for q in Q], repeat=5)
        rec.add(op="latency", engine=label, n=n, seconds=t["median"] / len(Q),
                runs=[r / len(Q) for r in t["runs"]])
    trees = {}
    for b in B:
        trees.setdefault(b.chrom, IntervalTree()).addi(b.start, b.end)
    t = timeit(lambda: [trees[q.chrom].overlap(q.start, q.end) for q in Q], repeat=5)
    rec.add(op="latency", engine="intervaltree", n=n, seconds=t["median"] / len(Q),
            runs=[r / len(Q) for r in t["runs"]])
    gb_ = pr.read_bed(str(DATA / f"peaks_B_{n}.bed"))
    Qs = Q[:20]
    t = timeit(lambda: [gb_[q.chrom, q.start:q.end] for q in Qs], repeat=3)
    rec.add(op="latency", engine="pyranges", n=n, seconds=t["median"] / len(Qs),
            runs=[r / len(Qs) for r in t["runs"]])
    db = to_df(B)
    t = timeit(lambda: [bf.select(db, (q.chrom, q.start, q.end)) for q in Qs], repeat=3)
    rec.add(op="latency", engine="bioframe", n=n, seconds=t["median"] / len(Qs),
            runs=[r / len(Qs) for r in t["runs"]])


def bench_caching(rec: Recorder):
    """The point index is built once and memoised on the Loci: later lookups are free."""
    print("\n== index caching ==")
    for n in SIZES:
        B = Loci.make(str(DATA / f"peaks_B_{n}.bed"))
        Bc = fresh(B)
        t_build = timeit(lambda: Bc.overlap_rows("chr1", 0, 1), repeat=3, setup=Bc._dirty)
        Bc.overlap_rows("chr1", 0, 1)
        t_query = timeit(lambda: Bc.overlap_rows("chr1", 0, 1), repeat=3)
        rec.add(op="caching", engine="index build", n=n, seconds=t_build["median"],
                runs=t_build["runs"])
        rec.add(op="caching", engine="query (warm)", n=n, seconds=t_query["median"],
                runs=t_query["runs"])


if __name__ == "__main__":
    rec = Recorder("loci")
    bench_intersect(rec)
    bench_merge(rec)
    bench_make(rec)
    bench_latency(rec)
    bench_caching(rec)
    rec.save(sizes=SIZES)
