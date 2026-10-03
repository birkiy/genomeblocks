#!/usr/bin/env python3
"""Measured upside of small fixes for the slow spots the suite found.

Each pair runs the current genomeblocks code path and a drop-in variant with
identical output (checked), so the gain is attributable to that one change.

  loci    Loci.sort() / merge(): dataclass __lt__ comparisons vs a (chrom, start) key
  motifs  scan_motifs_matrix: one scan() per (window, motif) vs one native
          call per motif over the concatenated windows
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


def scan_concat(seqs, path, thr=13.0, pseudo=0.1):
    """Per-window hit counts with one native call per motif.

    Windows of equal length are concatenated and striped once; each motif is
    scored over the whole block with ``calculate().threshold()``, and hits are
    split back per window with numpy (hits crossing a window edge dropped).
    """
    import numpy as np
    import lightmotif
    W = len(seqs[0])
    striped = lightmotif.stripe("".join(seqs))
    n = len(seqs)
    names, cols = [], []
    for m in lightmotif.load(path, format="jaspar16"):
        pssm = m.counts.normalize(pseudo).log_odds()
        w = len(m.counts)
        pos = np.asarray(pssm.calculate(striped).threshold(thr), dtype=np.int64)
        pos = pos[pos <= n * W - w]
        ok = (pos % W) + w <= W
        names.append(m.name)
        cols.append(np.bincount(pos[ok] // W, minlength=n))
    return names, np.stack(cols, axis=1)


def part_motifs(rec):
    import numpy as np
    from bench_motifs import JASPAR, R, THR, subset_file, windows
    from genomeblocks.motifs import make_genome, scan_motifs_matrix
    genome = make_genome(str(DATA / "genome.fa"))
    for M in (100, 1019):
        path = subset_file(M if M < 1019 else None)
        L = windows(genome, 1000)
        seqs = [l.sequence(genome, r=R).upper() for l in L]
        res = {}
        t = timeit(lambda: res.__setitem__("a", scan_motifs_matrix(
            L, genome, path, motif_format="jaspar16", r=R, threshold=THR, norm=False,
            workers=1, verbose=False)), repeat=3)
        rec.add(part="motifs", op=f"scan 1,000 windows × {M} motifs", variant="current (scan() per window × motif)",
                n=M, seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: res.__setitem__("b", scan_concat(seqs, path, THR)), repeat=3)
        names, mat = res["b"]
        same = bool(np.array_equal(res["a"][names].to_numpy(), mat))
        rec.add(part="motifs", op=f"scan 1,000 windows × {M} motifs", variant="concatenated, 1 call per motif",
                n=M, seconds=t["median"], runs=t["runs"], same=same)


def part_loci(rec):
    for n in (100_000, 1_000_000):
        L = Loci.make(str(DATA / f"peaks_AB_{n}.bed"))     # 2n unsorted intervals
        res = {}
        t = timeit(lambda: res.__setitem__("a", L.sort()), repeat=3)
        rec.add(part="loci", op="sort", variant="current (__lt__)", n=len(L), seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: res.__setitem__("b", sort_key(L)), repeat=3)
        rec.add(part="loci", op="sort", variant="key=(chrom, start)", n=len(L), seconds=t["median"],
                runs=t["runs"], same=[x.uid for x in res["a"]] == [x.uid for x in res["b"]])
        t = timeit(lambda: res.__setitem__("a", L.merge()), repeat=3)
        rec.add(part="loci", op="merge", variant="current", n=len(L), seconds=t["median"], runs=t["runs"])
        t = timeit(lambda: res.__setitem__("b", merge_key(L)), repeat=3)
        rec.add(part="loci", op="merge", variant="key sort + single pass", n=len(L), seconds=t["median"],
                runs=t["runs"], same=[x.uid for x in res["a"]] == [x.uid for x in res["b"]])


if __name__ == "__main__":
    import sys
    rec = Recorder("fixes")
    for part in (sys.argv[1:] or ["loci", "motifs"]):
        globals()[f"part_{part}"](rec)
    rec.save()
