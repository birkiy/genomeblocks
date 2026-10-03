#!/usr/bin/env python3
"""Measured upside of small fixes for the slow spots the suite found.

Each pair runs the current genomeblocks code path and a drop-in variant with
identical output (checked), so the gain is attributable to that one change.

  sort    Loci.sort(): dataclass __lt__ comparisons vs a (chrom, start) key
  merge   Loci.merge(): current vs key-sorted single pass
"""
from __future__ import annotations

from common import DATA, Recorder, timeit

from genomeblocks import Loci
from genomeblocks.locus import Locus


def sort_key(L):
    out = Loci(sorted(L, key=lambda l: (l.chrom, l.start)))
    out.uids
    return out


def merge_key(L):
    """Loci.merge with the key-based sort; same single pass otherwise."""
    if not L:
        return Loci()
    s = sorted(L, key=lambda l: (l.chrom, l.start))
    out = []
    c, a, b, st = s[0].chrom, s[0].start, s[0].end, s[0].strand
    for l in s[1:]:
        if l.chrom == c and l.start <= b:
            if l.end > b:
                b = l.end
        else:
            out.append(Locus(c, a, b, st))
            c, a, b, st = l.chrom, l.start, l.end, l.strand
    out.append(Locus(c, a, b, st))
    return Loci(out)


if __name__ == "__main__":
    rec = Recorder("fixes")
    for n in (100_000, 1_000_000):
        L = Loci.make(str(DATA / f"peaks_AB_{n}.bed"))     # 2n unsorted intervals
        res = {}
        t = timeit(lambda: res.__setitem__("a", L.sort()), repeat=3)
        rec.add(op="sort", variant="current (__lt__)", n=len(L), seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: res.__setitem__("b", sort_key(L)), repeat=3)
        rec.add(op="sort", variant="key=(chrom, start)", n=len(L), seconds=t["median"],
                runs=t["runs"], same=[x.uid for x in res["a"]] == [x.uid for x in res["b"]])
        t = timeit(lambda: res.__setitem__("a", L.merge()), repeat=3)
        rec.add(op="merge", variant="current", n=len(L), seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: res.__setitem__("b", merge_key(L)), repeat=3)
        rec.add(op="merge", variant="key sort + single pass", n=len(L), seconds=t["median"],
                runs=t["runs"], same=[x.uid for x in res["a"]] == [x.uid for x in res["b"]])
    rec.save()
