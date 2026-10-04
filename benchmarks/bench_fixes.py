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


def _genome_index(loci):
    """Windows on one genome-wide axis: chrom -> offset, sorted global
    starts/ends, and the original loci index of each sorted window."""
    import numpy as np
    chroms = sorted({l.chrom for l in loci})
    cid = {c: i for i, c in enumerate(chroms)}
    span = np.zeros(len(chroms), np.int64)
    for l in loci:
        span[cid[l.chrom]] = max(span[cid[l.chrom]], l.end)
    off = np.concatenate([[0], np.cumsum(span + 1)[:-1]])
    gs = np.fromiter((off[cid[l.chrom]] + l.start for l in loci), np.int64, len(loci))
    ge = np.fromiter((off[cid[l.chrom]] + l.end for l in loci), np.int64, len(loci))
    order = np.argsort(gs, kind="stable")
    return cid, off, gs[order], ge[order], order


def _locate_global(codes, cats, pos, cid, off, gs, ge, order):
    """Window index (original loci order) for each (chrom code, pos); -1 if none."""
    import numpy as np
    lut = np.array([cid.get(str(c), -1) for c in cats], np.int64)
    c = lut[codes]
    ok = c >= 0
    g = np.where(ok, off[np.maximum(c, 0)] + pos, -1)
    i = np.searchsorted(gs, g, side="right") - 1
    ok &= i >= 0
    ok &= g < ge[np.maximum(i, 0)]
    return np.where(ok, order[np.maximum(i, 0)], -1)


def count_pairs_fast(loci, path, chunksize=2_000_000):
    """count_pairs (all partners) with integer chrom codes and one global
    searchsorted per anchor; same output frame as genomeblocks'."""
    import numpy as np
    import pandas as pd
    from genomeblocks.bedpe import read_pairs_chunks
    cid, off, gs, ge, order = _genome_index(loci)
    names, acc = {}, []
    for ch in read_pairs_chunks(path, format="pairs", chunksize=chunksize):
        for a, b in (("1", "2"), ("2", "1")):
            ca, cb = ch["chrom" + a], ch["chrom" + b]
            w = _locate_global(ca.cat.codes.to_numpy(), ca.cat.categories, ch["pos" + a].to_numpy(),
                               cid, off, gs, ge, order)
            p_lut = np.array([names.setdefault(str(c), len(names)) for c in cb.cat.categories], np.int64)
            p = p_lut[cb.cat.codes.to_numpy()]
            hit = w >= 0
            acc.append((w[hit], p[hit]))
    P = len(names)
    w = np.concatenate([x for x, _ in acc]); p = np.concatenate([y for _, y in acc])
    mat = np.bincount(w * P + p, minlength=len(loci) * P).reshape(len(loci), P)
    out = pd.DataFrame({"chrom": [l.chrom for l in loci], "start": [l.start for l in loci],
                        "end": [l.end for l in loci], "uid": [l.uid for l in loci]})
    for name in sorted(names):
        col = mat[:, names[name]]
        if col.any():
            out[name] = col
    return out


def count_pairs_2d_fast(loci, path, chunksize=2_000_000):
    """count_pairs_2d (loci_b=None) the same way: one locate per anchor,
    then a single sparse COO sum."""
    import numpy as np
    import scipy.sparse as sp
    from genomeblocks.bedpe import read_pairs_chunks
    cid, off, gs, ge, order = _genome_index(loci)
    rows, cols = [], []
    for ch in read_pairs_chunks(path, format="pairs", chunksize=chunksize):
        w1 = _locate_global(ch["chrom1"].cat.codes.to_numpy(), ch["chrom1"].cat.categories,
                            ch["pos1"].to_numpy(), cid, off, gs, ge, order)
        w2 = _locate_global(ch["chrom2"].cat.codes.to_numpy(), ch["chrom2"].cat.categories,
                            ch["pos2"].to_numpy(), cid, off, gs, ge, order)
        ok = (w1 >= 0) & (w2 >= 0)
        rows += [w1[ok], w2[ok]]; cols += [w2[ok], w1[ok]]
    r = np.concatenate(rows); c = np.concatenate(cols)
    n = len(loci)
    m = sp.coo_matrix((np.ones(len(r), np.int64), (r, c)), shape=(n, n))
    m.sum_duplicates()
    return m.tocsr()


def part_pairs(rec):
    import numpy as np
    from common import read_chromsizes
    from genomeblocks.bedpe import count_pairs, count_pairs_2d
    path = str(DATA / "hic.pairs")
    cs = read_chromsizes()
    w50, w500 = Loci.tile_genome(cs, 50_000), Loci.tile_genome(cs, 500_000)
    res = {}
    t = timeit(lambda: res.__setitem__("a", count_pairs(w50, path, format="pairs", verbose=False)), repeat=3)
    rec.add(part="pairs", op="count_pairs (50 kb x partner chrom)", variant="current", n=5_000_000,
            seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: res.__setitem__("b", count_pairs_fast(w50, path)), repeat=3)
    same = list(res["a"].columns) == list(res["b"].columns) and res["a"].equals(res["b"])
    rec.add(part="pairs", op="count_pairs (50 kb x partner chrom)", variant="integer codes + global searchsorted",
            n=5_000_000, seconds=t["median"], runs=t["runs"], same=bool(same))
    t = timeit(lambda: res.__setitem__("a", count_pairs_2d(w500, path, format="pairs", verbose=False)),
               repeat=1, warmup=0)
    rec.add(part="pairs", op="count_pairs_2d (500 kb)", variant="current", n=5_000_000,
            seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: res.__setitem__("b", count_pairs_2d_fast(w500, path)), repeat=3)
    same = (res["a"] != res["b"]).nnz == 0
    rec.add(part="pairs", op="count_pairs_2d (500 kb)", variant="integer codes + global searchsorted",
            n=5_000_000, seconds=t["median"], runs=t["runs"], same=bool(same))


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
    for part in (sys.argv[1:] or ["loci", "motifs", "pairs"]):
        globals()[f"part_{part}"](rec)
    rec.save()
