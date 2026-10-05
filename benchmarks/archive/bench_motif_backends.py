#!/usr/bin/env python3
"""Motif backends: lightmotif (today), MOODS, and a MOODS-style numba scanner.

Every engine fills the same (windows x motifs) hit-count matrix that
``scan_motifs_matrix(norm=False)`` returns: 1,019 JASPAR motifs, 500 bp
windows, forward strand, log2-odds >= 13 with pseudocount 0.1.
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from bench_motifs import JASPAR, R, THR, logodds, motif_counts, windows
from common import DATA, Recorder, timeit

from genomeblocks.motifs import make_genome, scan_motifs_matrix

_MOODS: dict = {}


def _moods_init(mats):
    import MOODS.scan
    sc = MOODS.scan.Scanner(7)
    sc.set_motifs([logodds(c).T.tolist() for _, c in mats], [0.25] * 4, [THR] * len(mats))
    _MOODS["sc"] = sc


def _moods_chunk(seqs):
    sc = _MOODS["sc"]
    return [[len(h) for h in sc.scan(s)] for s in seqs]


def moods_matrix(seqs, mats, workers=1):
    """MOODS as a scan_motifs_matrix backend: same matrix shape and order."""
    if workers == 1:
        _moods_init(mats)
        return np.array(_moods_chunk(seqs), np.int32)
    chunks = [seqs[i::workers] for i in range(workers)]
    with ProcessPoolExecutor(workers, initializer=_moods_init, initargs=(mats,)) as ex:
        parts = list(ex.map(_moods_chunk, chunks))
    out = np.zeros((len(seqs), len(mats)), np.int32)
    for i, part in enumerate(parts):
        out[i::workers] = part
    return out


def lookahead_matrix(enc, packed, threads):
    import numba
    from prototypes.lookahead import scan_counts
    numba.set_num_threads(threads)
    return scan_counts(*enc, *packed, THR)


if __name__ == "__main__":
    from prototypes.lookahead import encode, prepare, work_fraction
    rec = Recorder("motif_backends")
    genome = make_genome(str(DATA / "genome.fa"))
    mats = motif_counts()
    names = [n for n, _ in mats]
    packed = prepare([logodds(c) for _, c in mats])
    for N in (1000, 5000):
        L = windows(genome, N, seed=1)
        seqs = [l.sequence(genome, r=R).upper() for l in L]
        enc = encode(seqs)
        res = {}
        for w in (1, 4):
            t = timeit(lambda: res.__setitem__("ref", scan_motifs_matrix(
                L, genome, JASPAR, motif_format="jaspar16", r=R, threshold=THR, norm=False,
                workers=w, verbose=False)), repeat=3)
            rec.add(part="scan", engine="lightmotif (scan_motifs_matrix, today)", workers=w, n_seqs=N,
                    seconds=t["median"], runs=t["runs"])
        ref = res["ref"][names].to_numpy()
        for w in (1, 4):
            t = timeit(lambda: res.__setitem__("m", moods_matrix(seqs, mats, w)), repeat=3)
            d = int((res["m"] != ref).sum())
            rec.add(part="scan", engine="MOODS backend", workers=w, n_seqs=N, seconds=t["median"],
                    runs=t["runs"], cells_differing=d, cells=int(ref.size))
        for w in (1, 4):
            t = timeit(lambda: res.__setitem__("l", lookahead_matrix(enc, packed, w)), repeat=3)
            d = int((res["l"] != ref).sum())
            rec.add(part="scan", engine="lookahead scanner (numba prototype)", workers=w, n_seqs=N,
                    seconds=t["median"], runs=t["runs"], cells_differing=d, cells=int(ref.size))
        if N == 1000:
            frac = work_fraction(*enc, *packed, THR)
            rec.add(part="pruning", n_seqs=N, seconds=0.0, columns_scored_fraction=float(frac))
    rec.save(threshold=THR)
