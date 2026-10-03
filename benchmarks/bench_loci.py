#!/usr/bin/env python3
"""Loci interval algebra: genomeblocks vs the usual Python/CLI tools.

What is measured (A and B are n-peak sets, n = 1k .. 1M):

  intersect   A & B  -> loci of A overlapping any locus of B   (bedtools -u)
  merge       A.merge()                                         (bedtools merge)
  make        Loci.make(bed)  BED parsing into Locus objects
  latency     one overlap query against an indexed set (interactive use)
  caching     1st call (index build + query) vs. later calls (query only)

Engines: genomeblocks with cgranges (C index, the conda path), genomeblocks
with its pure-Python fallback (the pip path), a bounded-scan variant of that
fallback (a proposed fix), pyranges, bioframe, bedtools (CLI, includes file
I/O + process start), intervaltree, and a naive per-chromosome double loop.
"""
from __future__ import annotations

import bisect
import subprocess

import numpy as np
import pandas as pd

from common import DATA, Recorder, timeit

import genomeblocks.loci as gl
from genomeblocks import Loci

SIZES = [1_000, 10_000, 100_000, 1_000_000]


# ── engines ──────────────────────────────────────────────────────────────────

class BoundedPyIndex(gl._PyIntervalIndex):
    """The fallback index with the scan bounded on the left.

    Tracks the longest interval per chromosome so a query only walks intervals
    whose start lies in [qs - max_len, qe) — O(log n + k) instead of O(n).
    """
    __slots__ = ("_maxlen",)

    def index(self):
        super().index()
        self._maxlen = {c: max((e - s for s, e, _ in ivs), default=0)
                        for c, (ivs, _) in self._sorted.items()}

    def overlap(self, chrom, start, end):
        data = self._sorted.get(chrom)
        if data is None:
            return
        ivs, starts = data
        qs, qe = int(start), int(end)
        lo = bisect.bisect_left(starts, qs - self._maxlen[chrom])
        hi = bisect.bisect_left(starts, qe)
        for k in range(lo, hi):
            s, e, label = ivs[k]
            if e > qs:
                yield (s, e, label)


def use_index(kind: str):
    """Point genomeblocks at an interval-index implementation."""
    import cgranges
    if kind == "cgranges":
        gl._get_cgranges = lambda: cgranges
        gl._PyIntervalIndex = _ORIG_PY
    elif kind == "python":
        gl._get_cgranges = lambda: None
        gl._PyIntervalIndex = _ORIG_PY
    elif kind == "python-bounded":
        gl._get_cgranges = lambda: None
        gl._PyIntervalIndex = BoundedPyIndex


_ORIG_PY = gl._PyIntervalIndex


def fresh(L: Loci) -> Loci:
    """Same loci, no cached index."""
    out = Loci(L)
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
    return pd.DataFrame({"chrom": [l.chrom for l in L],
                         "start": [l.start for l in L],
                         "end": [l.end for l in L]})


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

        # genomeblocks + cgranges: cold (index built inside the call) and warm
        use_index("cgranges")
        Bc = fresh(B)
        res = {}
        t = timeit(lambda: res.__setitem__("o", A & Bc), repeat=rep,
                   setup=lambda: setattr(Bc, "_cgr", None))
        add("genomeblocks (cgranges)", t, len(res["o"]))
        Bc.cgr
        t = timeit(lambda: res.__setitem__("o", A & Bc), repeat=rep)
        add("genomeblocks (cgranges, warm index)", t, len(res["o"]))

        # pure-Python fallback — quadratic, so only at small n
        if n <= 100_000:
            use_index("python")
            Bp = fresh(B)
            t = timeit(lambda: res.__setitem__("o", A & Bp),
                       repeat=1 if n == 100_000 else rep, warmup=0 if n == 100_000 else 1,
                       setup=lambda: setattr(Bp, "_cgr", None))
            add("genomeblocks (pure-Python fallback)", t, len(res["o"]))
        use_index("python-bounded")
        Bb = fresh(B)
        t = timeit(lambda: res.__setitem__("o", A & Bb), repeat=rep,
                   setup=lambda: setattr(Bb, "_cgr", None))
        add("genomeblocks (fallback, bounded scan)", t, len(res["o"]))
        use_index("cgranges")

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
    for kind, label in (("cgranges", "genomeblocks (cgranges)"),
                        ("python", "genomeblocks (pure-Python fallback)"),
                        ("python-bounded", "genomeblocks (fallback, bounded scan)")):
        use_index(kind)
        Bk = fresh(B); Bk.cgr
        t = timeit(lambda: [Bk.overlaps(q) for q in Q], repeat=5)
        rec.add(op="latency", engine=label, n=n, seconds=t["median"] / len(Q),
                runs=[r / len(Q) for r in t["runs"]])
    use_index("cgranges")
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
    """Index is built once and memoised on the Loci: repeated queries are free."""
    print("\n== index caching ==")
    use_index("cgranges")
    for n in SIZES:
        A = Loci.make(str(DATA / f"peaks_A_{n}.bed"))
        B = Loci.make(str(DATA / f"peaks_B_{n}.bed"))
        Bc = fresh(B)
        t_build = timeit(lambda: Bc._build_cgr(), repeat=3)
        Bc.cgr
        t_query = timeit(lambda: A & Bc, repeat=3)
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
