"""Prototype: MOODS-style lookahead-filtered PWM scanning in numba.

For every motif the columns are visited most-informative first. After each
column the scanner knows the best score still reachable (score so far + the
best possible value of every remaining column); once that bound falls below
the threshold the position is abandoned. At a stringent threshold almost all
positions die after one or two columns, which is where MOODS gets its speed.

Output matches ``scan_motifs_matrix(norm=False)``: hit counts per
(sequence, motif), forward strand, score >= threshold.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange

CODE = np.full(256, 4, np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
    CODE[_b + 32] = _i                       # lower case


def encode(seqs):
    """Concatenate sequences into one uint8 array (A=0 C=1 G=2 T=3 other=4)."""
    lens = np.array([len(s) for s in seqs], np.int64)
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)
    seq = CODE[np.frombuffer("".join(seqs).encode(), np.uint8)]
    return seq, starts, lens


def prepare(logodds):
    """Pack (W, 4) log-odds matrices for the kernel.

    order[m, j]  column visited at step j (largest score spread first)
    bound[m, j]  best score the columns order[m, j:] can still add
    """
    M = len(logodds)
    W = np.array([len(x) for x in logodds], np.int64)
    Wmax = int(W.max())
    L = np.zeros((M, Wmax, 4), np.float64)
    order = np.zeros((M, Wmax), np.int64)
    bound = np.zeros((M, Wmax + 1), np.float64)
    for m, x in enumerate(logodds):
        w = len(x)
        L[m, :w] = x
        o = np.argsort(-(x.max(axis=1) - x.min(axis=1)), kind="stable")
        order[m, :w] = o
        colmax = x.max(axis=1)[o]
        bound[m, :w] = np.cumsum(colmax[::-1])[::-1]
    return L, order, bound, W


@njit(parallel=True, cache=True)
def scan_counts(seq, starts, lens, L, order, bound, W, thr):
    n_seq, M = len(starts), L.shape[0]
    out = np.zeros((n_seq, M), np.int32)
    for m in prange(M):
        w = W[m]
        for si in range(n_seq):
            s0, n = starts[si], lens[si]
            c = 0
            for p in range(n - w + 1):
                score = 0.0
                ok = True
                for j in range(w):
                    col = order[m, j]
                    b = seq[s0 + p + col]
                    if b > 3:
                        ok = False
                        break
                    score += L[m, col, b]
                    if score + bound[m, j + 1] < thr:
                        ok = False
                        break
                if ok:
                    c += 1
            out[si, m] = c
    return out


@njit(cache=True)
def work_fraction(seq, starts, lens, L, order, bound, W, thr):
    """Columns actually scored / columns a full scan would score."""
    done, full = 0, 0
    for m in range(L.shape[0]):
        w = W[m]
        for si in range(len(starts)):
            s0, n = starts[si], lens[si]
            for p in range(n - w + 1):
                full += w
                score = 0.0
                for j in range(w):
                    done += 1
                    col = order[m, j]
                    b = seq[s0 + p + col]
                    if b > 3:
                        break
                    score += L[m, col, b]
                    if score + bound[m, j + 1] < thr:
                        break
    return done / full
