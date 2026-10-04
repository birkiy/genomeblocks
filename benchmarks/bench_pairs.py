#!/usr/bin/env python3
"""Hi-C pair counting: vectorised streaming vs per-pair loops and cooler.

Input: 5M read pairs (4DN .pairs, 75% cis, P(s) ~ s^-1).

  count_pairs     50 kb windows x partner chromosome, one streaming pass
  count_pairs_2d  500 kb window x window sparse matrix
  baselines       naive per-pair Python loop (cgranges lookup per anchor),
                  pandas parse-only (the I/O floor), cooler cload pairs
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from common import DATA, Recorder, read_chromsizes, timeit

from genomeblocks import Loci
from genomeblocks.bedpe import count_pairs, count_pairs_2d, read_pairs_chunks

PAIRS = str(DATA / "hic.pairs")
N_PAIRS = 5_000_000
BIN = Path(sys.executable).parent


def naive_count(loci, path, limit):
    """Per-pair Python loop: look each anchor up in the window index."""
    from collections import defaultdict
    counts = defaultdict(lambda: np.zeros(len(loci), np.int64))
    cg = loci.cgr
    with open(path) as f:
        n = 0
        for line in f:
            if line[0] == "#":
                continue
            _, c1, p1, c2, p2, *_ = line.split("\t")
            p1, p2 = int(p1), int(p2)
            for c, p, other in ((c1, p1, c2), (c2, p2, c1)):
                for *_, i in cg.overlap(c, p, p + 1):
                    counts[other][i] += 1
            n += 1
            if n == limit:
                break
    return counts


def parse_only(path):
    n = 0
    for ch in read_pairs_chunks(path, format="pairs"):
        n += len(ch)
    return n


def cooler_cload(path, binsize):
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "m.cool"
        subprocess.run([str(BIN / "cooler"), "cload", "pairs", "-c1", "2", "-p1", "3",
                        "-c2", "4", "-p2", "5",
                        f"{DATA / 'hg38.chrom.sizes'}:{binsize}", path, str(out)],
                       check=True, capture_output=True)


if __name__ == "__main__":
    rec = Recorder("pairs")
    cs = read_chromsizes()
    win50 = Loci.tile_genome(cs, 50_000)
    win500 = Loci.tile_genome(cs, 500_000)
    print(f"windows: {len(win50):,} x 50 kb, {len(win500):,} x 500 kb")

    t = timeit(lambda: parse_only(PAIRS), repeat=3)
    rec.add(engine="pandas chunked parse only (I/O floor)", n_pairs=N_PAIRS,
            seconds=t["median"], runs=t["runs"], rate=N_PAIRS / t["median"])

    res = {}
    t = timeit(lambda: res.__setitem__("cp", count_pairs(win50, PAIRS, format="pairs",
                                                          verbose=False)), repeat=3)
    rec.add(engine="genomeblocks count_pairs (50 kb x partner chrom)", n_pairs=N_PAIRS,
            seconds=t["median"], runs=t["runs"], rate=N_PAIRS / t["median"])

    t = timeit(lambda: count_pairs(win50, PAIRS, format="pairs", target_chrom="chr1",
                                   verbose=False), repeat=3)
    rec.add(engine="genomeblocks count_pairs (target_chrom=chr1)", n_pairs=N_PAIRS,
            seconds=t["median"], runs=t["runs"], rate=N_PAIRS / t["median"])

    t = timeit(lambda: res.__setitem__("m2", count_pairs_2d(win500, PAIRS, format="pairs",
                                                             verbose=False)), repeat=3)
    rec.add(engine="genomeblocks count_pairs_2d (500 kb)", n_pairs=N_PAIRS,
            seconds=t["median"], runs=t["runs"], rate=N_PAIRS / t["median"],
            nnz=int(res["m2"].nnz))

    t = timeit(lambda: cooler_cload(PAIRS, 500_000), repeat=3)
    rec.add(engine="cooler cload pairs (500 kb, CLI)", n_pairs=N_PAIRS,
            seconds=t["median"], runs=t["runs"], rate=N_PAIRS / t["median"])

    lim = 200_000
    t = timeit(lambda: res.__setitem__("nv", naive_count(win50, PAIRS, lim)), repeat=1)
    # agreement on the same first `lim` pairs
    head = DATA / "hic_head.pairs"
    with open(PAIRS) as f, open(head, "w") as g:
        k = 0
        for line in f:
            g.write(line)
            k += line[0] != "#"
            if k == lim:
                break
    vec = count_pairs(win50, str(head), format="pairs", verbose=False)
    same = all(np.array_equal(vec[c].to_numpy(), res["nv"][c]) for c in res["nv"])
    rec.add(engine="naive per-pair loop (cgranges lookup)", n_pairs=lim,
            seconds=t["median"], runs=t["runs"], rate=lim / t["median"],
            agrees_with_count_pairs=bool(same))
    rec.save(n_pairs=N_PAIRS)
