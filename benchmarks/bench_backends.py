#!/usr/bin/env python3
"""Interval backends: what switching engines costs, and what pays off.

For A & B (A's loci overlapping B) and merge, at n = 1k .. 1M:

  cgranges (today)            genomeblocks' object path
  numpy, from objects         convert Loci -> arrays on every call, then numpy
  numpy, columnar             arrays already stored (no conversion)
  pyranges, from Loci         Loci.to_pyranges() on every call, then pyranges
  pyranges, native            frames already built

Also: conversion costs on their own, the cost of an 'auto' backend decision,
and single-lookup latency for each engine.
"""
from __future__ import annotations

import sys
import time

import numpy as np

from common import DATA, Recorder, timeit
from prototypes import backends as B
from prototypes.npintervals import Intervals, merge, overlaps_any

from genomeblocks import Loci

SIZES = [1_000, 10_000, 100_000, 1_000_000]


def part_ops(rec):
    import pyranges as pr
    for n in SIZES:
        A = Loci.make(str(DATA / f"peaks_A_{n}.bed"))
        Bl = Loci.make(str(DATA / f"peaks_B_{n}.bed"))
        AB = Loci.make(str(DATA / f"peaks_AB_{n}.bed"))
        rep = 5 if n <= 100_000 else 3
        res = {}
        Bl.cgr
        t = timeit(lambda: res.__setitem__("ref", A & Bl), repeat=rep)
        ref = {l.uid for l in res["ref"]}
        rec.add(part="ops", op="A & B", engine="cgranges (today)", n=n, seconds=t["median"], runs=t["runs"])

        def np_obj():
            m = overlaps_any(Intervals.from_loci(A), Intervals.from_loci(Bl))
            return Loci([A[i] for i in np.flatnonzero(m).tolist()])
        t = timeit(lambda: res.__setitem__("o", np_obj()), repeat=rep)
        rec.add(part="ops", op="A & B", engine="numpy, converting from objects", n=n, seconds=t["median"],
                runs=t["runs"], same={l.uid for l in res["o"]} == ref)
        ia, ib = Intervals.from_loci(A), Intervals.from_loci(Bl)
        t = timeit(lambda: res.__setitem__("m", overlaps_any(ia, ib)), repeat=rep)
        rec.add(part="ops", op="A & B", engine="numpy, columnar (no conversion)", n=n, seconds=t["median"],
                runs=t["runs"], same={A[i].uid for i in np.flatnonzero(res["m"]).tolist()} == ref)
        if n <= 100_000:
            t = timeit(lambda: A.to_pyranges().overlap(Bl.to_pyranges()), repeat=3)
            rec.add(part="ops", op="A & B", engine="pyranges, converting from Loci", n=n,
                    seconds=t["median"], runs=t["runs"])
        ga, gb = pr.read_bed(str(DATA / f"peaks_A_{n}.bed")), pr.read_bed(str(DATA / f"peaks_B_{n}.bed"))
        t = timeit(lambda: ga.overlap(gb), repeat=rep)
        rec.add(part="ops", op="A & B", engine="pyranges, native frames", n=n, seconds=t["median"], runs=t["runs"])

        # merge of the unsorted 2n union
        t = timeit(lambda: res.__setitem__("mr", AB.merge()), repeat=3)
        mref = sorted((l.chrom, l.start, l.end) for l in res["mr"])
        rec.add(part="ops", op="merge", engine="cgranges (today)", n=2 * n, seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: res.__setitem__("mo", merge(Intervals.from_loci(AB))), repeat=3)
        mo = res["mo"]
        rec.add(part="ops", op="merge", engine="numpy, converting from objects", n=2 * n, seconds=t["median"],
                runs=t["runs"], same=sorted(zip([mo.names[c] for c in mo.codes], mo.starts.tolist(), mo.ends.tolist())) == mref)
        iab = Intervals.from_loci(AB)
        t = timeit(lambda: merge(iab), repeat=rep)
        rec.add(part="ops", op="merge", engine="numpy, columnar (no conversion)", n=2 * n, seconds=t["median"], runs=t["runs"])
        gab = pr.read_bed(str(DATA / f"peaks_AB_{n}.bed"))
        t = timeit(lambda: gab.merge(), repeat=rep)
        rec.add(part="ops", op="merge", engine="pyranges, native frames", n=2 * n, seconds=t["median"], runs=t["runs"])


def part_convert(rec):
    import pandas as pd
    for n in (100_000, 1_000_000):
        A = Loci.make(str(DATA / f"peaks_A_{n}.bed"))
        ia = Intervals.from_loci(A)
        from genomeblocks.locus import Locus
        steps = [
            ("Loci -> numpy arrays", lambda: Intervals.from_loci(A)),
            ("numpy arrays -> Loci objects", lambda: ia.to_loci(Locus, Loci)),
            ("Loci.to_frame() (pandas)", lambda: A.to_frame()),
            ("Loci.to_pyranges()", lambda: A.to_pyranges()),
            ("numpy arrays -> pandas frame", lambda: pd.DataFrame({"chrom": pd.Categorical.from_codes(ia.codes, ia.names),
                                                                    "start": ia.starts, "end": ia.ends})),
        ]
        for name, fn in steps:
            if n == 1_000_000 and "to_pyranges" in name:
                continue
            t = timeit(fn, repeat=3)
            rec.add(part="convert", step=name, n=n, seconds=t["median"], runs=t["runs"])


def part_dispatch(rec):
    """How expensive is the decision itself?"""
    class Dummy:
        pass
    for name in ("cgranges", "numpy", "pyranges"):
        B.register("intervals", name, Dummy())
    B.set_policy("intervals", lambda op="", n=0: "numpy" if op != "lookup" and n >= 20_000 else "cgranges")
    N = 200_000
    t0 = time.perf_counter()
    for i in range(N):
        B.get("intervals", op="overlaps_any", n=i)
    per_auto = (time.perf_counter() - t0) / N
    with B.use(intervals="numpy"):
        t0 = time.perf_counter()
        for i in range(N):
            B.get("intervals", op="overlaps_any", n=i)
        per_fixed = (time.perf_counter() - t0) / N
    rec.add(part="dispatch", step="auto policy (size + op check)", seconds=per_auto)
    rec.add(part="dispatch", step="user-pinned backend", seconds=per_fixed)


def part_lookup(rec):
    """Single-interval latency: can a numpy engine serve point lookups too?"""
    n = 100_000
    Bl = Loci.make(str(DATA / f"peaks_B_{n}.bed"))
    Q = Loci.make(str(DATA / "peaks_A_1000.bed"))[:200]
    Bl.cgr
    t = timeit(lambda: [any(True for _ in Bl.cgr.overlap(q.chrom, q.start, q.end)) for q in Q], repeat=5)
    rec.add(part="lookup", engine="cgranges", seconds=t["median"] / len(Q), runs=[r / len(Q) for r in t["runs"]])
    ib = Intervals.from_loci(Bl)
    from prototypes.npintervals import STRIDE
    code = {c: i for i, c in enumerate(ib.names)}
    bs = ib.codes * STRIDE + ib.starts
    o = np.argsort(bs, kind="stable")
    bs, rmax = bs[o], np.maximum.accumulate((ib.codes * STRIDE + ib.ends)[o])

    def one(q):
        c = code[q.chrom] * STRIDE
        k = int(np.searchsorted(bs, c + q.end))
        return k > 0 and rmax[k - 1] > c + q.start
    t = timeit(lambda: [one(q) for q in Q], repeat=5)
    rec.add(part="lookup", engine="numpy (one searchsorted per query)", seconds=t["median"] / len(Q),
            runs=[r / len(Q) for r in t["runs"]])
    QI = Intervals.from_loci(Q)
    t = timeit(lambda: overlaps_any(QI, ib), repeat=5)
    rec.add(part="lookup", engine="numpy, all 200 queries in one call", seconds=t["median"] / len(Q),
            runs=[r / len(Q) for r in t["runs"]])


if __name__ == "__main__":
    parts = sys.argv[1:] or ["ops", "convert", "dispatch", "lookup"]
    rec = Recorder("backends")
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save()
