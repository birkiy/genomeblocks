#!/usr/bin/env python3
"""genomeblocks.columnar.Loci on the bench_loci inputs: make, A & B, merge.

Same files and the same identical-output checks as bench_loci, so the rows sit
next to the classic Loci and the other engines on the Benchmarks page.
"""
from __future__ import annotations

from common import DATA, Recorder, timeit

import genomeblocks.columnar as gbc
from genomeblocks import Loci

SIZES = [1_000, 10_000, 100_000, 1_000_000]


def key(L):
    return sorted((l.chrom, l.start, l.end) for l in L)


if __name__ == "__main__":
    rec = Recorder("loci_columnar")
    for n in SIZES:
        pa, pb, pab = (str(DATA / f"peaks_{s}_{n}.bed") for s in ("A", "B", "AB"))
        rep = 5 if n <= 100_000 else 3
        g = gbc.Genome()
        t = timeit(lambda: gbc.Loci.make(pa, genome=g, sort=False), repeat=rep)
        rec.add(op="make", engine="genomeblocks.columnar", n=n, seconds=t["median"], runs=t["runs"])

        A = gbc.Loci.make(pa, genome=g, sort=False)
        B = gbc.Loci.make(pb, genome=g, sort=False)
        res = {}
        t = timeit(lambda: res.__setitem__("o", A & B), repeat=rep)
        same = key(res["o"]) == key(Loci.make(pa) & Loci.make(pb)) if n <= 100_000 else None
        rec.add(op="intersect", engine="genomeblocks.columnar", n=n, seconds=t["median"], runs=t["runs"],
                n_out=len(res["o"]), agrees=same)

        AB = gbc.Loci.make(pab, genome=g, sort=False)
        t = timeit(lambda: res.__setitem__("m", AB.merge()), repeat=rep)
        same = key(res["m"]) == key(Loci.make(pab).merge()) if n <= 100_000 else None
        rec.add(op="merge", engine="genomeblocks.columnar", n=n, seconds=t["median"], runs=t["runs"],
                n_out=len(res["m"]), agrees=same)
    rec.save(sizes=SIZES)
