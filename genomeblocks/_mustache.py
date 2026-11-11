"""Loops from a Hi-C / HiChIP cooler: Mustache's scale-space caller, in numpy.

:func:`mustache` calls the loops Mustache 1.3.3 (Roayaei Ardakany et al.
2020, Genome Biol 21:256; github.com/ay-lab/mustache at df71aa1) calls with
``mustache -f <cool> -r <res> -o <out>``, its defaults otherwise. Line numbers
below are ``mustache/mustache.py`` at that commit. Per chromosome:

1. **Contacts** (``read_cooler`` :399): the balanced upper-triangle pixels
   ``count * w1 * w2`` (cooler's own product, ``weight`` column) within
   ``distance`` bins, NaN and non-positive values dropped. ``distance``
   defaults to 2 Mb (200 bins at >= 10 kb, 2000 bins at <= 1 kb) and is
   clamped to [200, 10 000] bins (:996). Chromosomes: those longer than 1 Mb.
2. **Normalisation** (``normalize_sparse`` :622), per diagonal: when the
   chromosome's last contact lies more than ``distance`` + 2 Mb from its start,
   each value (+ 0.001) is z-scored against the contacts within a 2 Mb window
   along its diagonal (``np.convolve``; the diagonal's mean / SD where the
   window holds < 30 contacts) and scaled by 1 + log30(1 + diagonal mean);
   otherwise against the whole diagonal, and the last diagonal (``distance``)
   is left as it is.
3. **Blocks** (``regulator`` :893): dense blocks of max(2 ``distance``, 2000)
   bins overlapping by ``distance`` bins; diagonals 0-4 and those beyond
   ``distance`` are set to 2. The tested pixels are the non-zero ones at least
   4 bins off the diagonal.
4. **Scale space** (:697): per octave sigma0 * 2^o, Gaussians of 12 sigmas
   sigma0 2^o 2^((k-1)/10) (radius ceil(2 sigma), reflected edges) and their
   11 differences L_k. A pixel is taken at L_k (k = 2..10) when L_k is a 3x3
   maximum there, larger than the 3x3 maxima of L_k-1 and L_k+1, a 3x3
   maximum of L_k-1 or L_k+1, and larger than any L it was taken at before;
   its P = 1 - F(|L_k|) for the exponential fitted (ML) to |L_k| of the
   block's tested pixels, its scale sigma_(k+1).
5. **Test** (:775): Benjamini-Hochberg over the taken pixels of the block;
   candidates FDR < ``pt``; a candidate is dropped when its row is 0, when
   fewer than ``st`` of the (2s+1)^2 pixels around it (s = ceil(scale)) or
   fewer than 0.6 of the (4s+1)^2 are tested pixels, or when its value is not
   above twice the mean non-zero value of its diagonal.
6. **Loops** (:828): the candidates' 3x3 neighbourhoods are joined into
   8-connected components; each gives the pixel of lowest FDR among all its
   pixels (the neighbours included, ties in row-major order). A block keeps
   the loops not already inside the previous block.

Mustache's quirks are kept because they change the calls: ``-i`` never
reaches the scale count (10 is hard-coded, :711); a block with fewer than
10 000 tested pixels gives no loops (``len(pFound)``, :775); the sparsity
window of a candidate within ``s`` bins of the block's first row or column is
empty (negative slice start) and drops it; the lowest-FDR pixel of a component
may be a neighbour that failed the filters. Here the blocks' Gaussians are
computed only on the band of tested pixels (the two passes of scipy's
``gaussian_filter`` on sheared copies of the band, so every value is the one
scipy gives on the full block), the 3x3 maxima on the band, the filters
vectorised; blocks run on ``threads`` threads (Mustache: 4 processes).

Against Mustache on the same balanced coolers (``benchmarks/bench_mustache.py``:
a synthetic cooler, HiChIP at 5 and 10 kb) the loops are the same and FDRs
agree to 1e-11 (relative), most of them to the bit, 9-15x faster on the same
4 cores. The rest is rounding: numpy 2's ``exp`` (scipy's Gaussian kernels)
and summation differ from numpy 1.26's, which Mustache needs; the box sums
here are cumulative, Mustache's dot products; and Mustache sums each diagonal
in the order of Python sets.
"""
from __future__ import annotations

import math
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from .bedpe import Pairs, open_cooler
from .genome import Genome
from .loci import Loci

_SCALES = 10            # mustache.py:711, ``s = 10`` (the CLI's -i is never used)

HEADER = "BIN1_CHR\tBIN1_START\tBIN1_END\tBIN2_CHROMOSOME\tBIN2_START\tBIN2_END\tFDR\tDETECTION_SCALE\n"


@dataclass
class Mustache:
    """The result of :func:`mustache`.

    ``loops`` is a :class:`~genomeblocks.Pairs` of bin pairs (anchors a / b,
    ``end = start + binsize`` as Mustache writes them) in Mustache's order
    (chromosome, block, component), with ``count`` (raw contacts of the
    pixel), ``balanced`` (count * w1 * w2), ``P`` (the exponential P of the
    pixel's detection scale; NaN for a component pixel that was not detected),
    ``FDR`` (Mustache's FDR column: BH q within the block; 1 for a pixel that
    is not a tested pixel, 2 for one that was never detected) and ``scale``
    (``DETECTION_SCALE``, the sigma of detection; 1 when not detected).
    Mustache has no expected count: its statistic is the scale-space
    difference of the normalised map."""
    loops: Pairs
    binsize: int
    distance: int
    pt: float
    timings: dict = field(default_factory=dict)

    def __repr__(self):
        return (f"Mustache(loops={len(self.loops):,} (FDR < {self.pt}), binsize={self.binsize}, "
                f"distance={self.distance:,})")

    def write(self, outdir: str, prefix: str) -> dict:
        """Mustache's output, ``{outdir}/{prefix}.tsv``: its header and one row per loop with the same text
        (``str`` of the numbers, as Mustache prints them). Returns the paths."""
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / f"{prefix}.tsv"
        L = self.loops
        c1, c2 = L.a.chroms.astype(str), L.b.chroms.astype(str)
        s1, e1, s2, e2 = L.a.starts, L.a.ends, L.b.starts, L.b.ends
        q, sc = L.cols["FDR"], L.cols["scale"]
        with open(path, "w") as f:
            f.write(HEADER)
            f.writelines(f"{c1[i]}\t{s1[i]}\t{e1[i]}\t{c2[i]}\t{s2[i]}\t{e2[i]}\t{float(q[i])!r}\t{float(sc[i])!r}\n"
                         for i in range(len(L)))
        return {"loops": str(path)}


def mustache(cool: str, *, resolution: Optional[int] = None, chromosomes: Optional[Union[str, Sequence[str]]] = None,
             distance: Optional[int] = None, pt: float = 0.2, st: float = 0.88, sigma0: float = 1.6,
             octaves: int = 2, weight: Union[str, bool, np.ndarray, None] = "weight", threads: int = 4,
             verbose: bool = True) -> Mustache:
    """Mustache loops (multi-scale blob detection) from a balanced cooler, with Mustache's defaults.

    Args:
        cool: a .cool, a cooler URI, or an .mcool with ``resolution``.
        chromosomes: the chromosomes to scan (Mustache's ``-ch``); default: every chromosome longer than 1 Mb.
        distance: largest loop distance in bp (``-d``); default 2 Mb (see the module notes for Mustache's
            resolution-dependent default and clamping, which apply).
        pt: FDR threshold of the loops (``-pt``, strict ``<``, at most 1).
        st: sparsity threshold (``-st``).
        sigma0, octaves: first scale and number of octaves (``-sz``, ``-oc``).
        weight: the balancing weights — a bin-table column name (``-norm``; ``True`` = ``"weight"``), an array
            with one weight per bin of the cooler, or ``None`` / ``False`` for raw counts (which Mustache never
            uses on a cooler).
        threads: blocks scanned in parallel (Mustache's ``-p``).
    """
    if not 0 < pt <= 1:
        raise ValueError("pt must be in (0, 1]")
    T = {}
    t0 = time.perf_counter()
    clr = open_cooler(cool, resolution, caller="mustache")
    res = int(clr.binsize)
    dist_bp = _distance_filter(distance, res)
    D = int(dist_bp // res)
    names = list(clr.chromnames)
    sizes = clr.chromsizes
    if chromosomes is None:
        chroms = [c for c in names if sizes[c] > 1_000_000]
    else:
        chroms = [chromosomes] if isinstance(chromosomes, str) else list(chromosomes)
        missing = [c for c in chroms if c not in names]
        if missing:
            raise ValueError(f"not in the cooler: {', '.join(missing)}")
    g = Genome.from_sizes(sizes)
    data = {}
    with clr.open("r") as h5:
        w = _weights(h5, weight)
        chrom_off = h5["indexes/chrom_offset"][:]
        bin1_off = h5["indexes/bin1_offset"]
        for c in chroms:
            ci = names.index(c)
            lo, hi = int(chrom_off[ci]), int(chrom_off[ci + 1])
            p0, p1 = int(bin1_off[lo]), int(bin1_off[hi])
            b1 = h5["pixels/bin1_id"][p0:p1]
            b2 = h5["pixels/bin2_id"][p0:p1]
            cnt = h5["pixels/count"][p0:p1]
            val = (w[b1] * w[b2]) * cnt if w is not None else cnt.astype(np.float64)
            x, y = b1 - lo, b2 - lo
            keep = (b2 < hi) & (y - x <= dist_bp / res) & (val > 0)
            if keep.any():
                data[c] = (x[keep].astype(np.int64), y[keep].astype(np.int64), val[keep], cnt[keep])
    T["read"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    with ThreadPoolExecutor(max(1, threads)) as ex:
        z = dict(zip(data, ex.map(lambda c: _normalize(data[c][0], data[c][1], data[c][2], res, D), data)))
    T["normalize"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    octave_values = [sigma0 * (2 ** i) for i in range(octaves)]
    jobs = []
    for c in chroms:
        if c not in data:
            continue
        x, y = data[c][0], data[c][1]
        starts, ends = _blocks(int(y.max()) + 1, D)
        for i in range(len(starts)):
            mask = -1 if i == 0 else (ends[i - 1] - starts[i] if i == len(starts) - 1 else D)
            jobs.append((c, starts[i], ends[i], mask))

    def run(job):
        c, s, e, mask = job
        x, y = data[c][0], data[c][1]
        return _block(x, y, z[c], s, e, mask, D, octave_values, pt, st)

    with ThreadPoolExecutor(max(1, threads)) as ex:
        found = list(ex.map(run, jobs))
    T["scale space"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    parts = {}
    for (c, *_), f in zip(jobs, found):
        if len(f[0]):
            parts.setdefault(c, []).append(f)
    cols = {k: [] for k in ("chrom", "x", "y", "count", "balanced", "P", "FDR", "scale")}
    for c, fs in parts.items():
        X = np.concatenate([f[0] for f in fs])
        Y = np.concatenate([f[1] for f in fs])
        x, y, val, cnt = data[c]
        n = int(y.max()) + 2
        key = x * n + y                                   # row-major: sorted
        k = np.searchsorted(key, X * n + Y)
        k = np.minimum(k, len(key) - 1)
        hit = key[k] == X * n + Y
        cols["chrom"].append(np.full(len(X), g.code[c], np.int32))
        cols["x"].append(X)
        cols["y"].append(Y)
        cols["count"].append(np.where(hit, cnt[k], 0))
        cols["balanced"].append(np.where(hit, val[k], 0.0))
        for j, name in ((2, "P"), (3, "FDR"), (4, "scale")):
            cols[name].append(np.concatenate([f[j] for f in fs]))
    cat = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in cols.items()}
    code = cat["chrom"].astype(np.int32)
    X, Y = cat["x"].astype(np.int64), cat["y"].astype(np.int64)
    loops = Pairs(Loci(code, X * res, (X + 1) * res, genome=g), Loci(code, Y * res, (Y + 1) * res, genome=g),
                  {"count": cat["count"].astype(np.int64), "balanced": cat["balanced"].astype(float),
                   "P": cat["P"].astype(float), "FDR": cat["FDR"].astype(float), "scale": cat["scale"].astype(float)})
    T["tables"] = time.perf_counter() - t0
    out = Mustache(loops, res, dist_bp, pt, T)
    if verbose:
        print(f"[mustache] {len(data)} chromosomes, {len(jobs)} blocks, {len(loops):,} loops (FDR < {pt}) | "
              + ", ".join(f"{k} {v:.1f}s" for k, v in T.items()))
    return out


# ── pieces ─────────────────────────────────────────────────────────────────

def _distance_filter(distance: Optional[int], res: int) -> int:
    """Mustache's distance limit in bp (``main``, mustache.py:996-1014)."""
    if not distance:
        if 200 * res >= 2_000_000:
            return 200 * res
        if 2000 * res <= 2_000_000:
            return 2000 * res
        return 2_000_000
    distance = int(distance)
    if distance < 200 * res:
        return 200 * res
    if distance > 10000 * res:
        return 10000 * res
    if distance > 10_000_000:
        return 10_000_000
    return distance


def _weights(h5, weight):
    """Per-bin multiplicative weights (None for raw counts)."""
    if weight is None or weight is False:
        return None
    if isinstance(weight, (str, bool)):
        name = "weight" if weight is True else weight
        if name not in h5["bins"]:
            raise ValueError(f"no bins/{name} in the cooler: balance it (cooler balance) or pass weight=")
        return h5["bins"][name][:].astype(np.float64)
    w = np.asarray(weight, np.float64)
    if w.shape != (h5["bins/start"].shape[0],):
        raise ValueError("weight array needs one value per bin of the cooler")
    return w


def _normalize(x, y, v, res, D):
    """``normalize_sparse`` (mustache.py:622) on the diagonals >= 4 (the only ones Mustache reads back):
    the normalised values, in the pixels' order.

    Mustache's local branch convolves each whole diagonal with a 2 Mb box (``np.convolve(..., 'same')``,
    i.e. the window [p - K//2, p + (K-1)//2]; the diagonal is never shorter than the box there) and reads the
    sums back at the pixels. Here the window sums are taken at the pixels only, from extended-precision
    cumulative sums (within an ulp or two of the convolution's dot products), except on a diagonal where, in
    a window of >= 30 contacts, the variance or the pixel's deviation from the window mean is at the level of
    rounding (equal values, e.g. integer counts): there the z-score is the rounding residue of the sums (or
    exactly 0, which drops the pixel), so the diagonal is convolved as Mustache does it."""
    v = v.copy()
    n = int(max(x.max(), y.max())) + 1
    d = y - x
    o = np.argsort(d, kind="stable")
    ds = d[o]
    bounds = np.searchsorted(ds, np.arange(D + 3))
    with np.errstate(all="ignore"):
        if (n - D) * res > 2_000_000:
            K = int(2_000_000 / res)
            h1, h2 = K // 2, (K - 1) // 2
            for dd in range(4, 2 + D):
                idx = o[bounds[dd]:bounds[dd + 1]]
                if len(idx) == 0:                                  # nothing to normalise
                    continue
                xs, vs = x[idx], v[idx]                            # xs ascending
                a = vs + 0.001
                std = np.std(vs)
                mean = np.mean(vs)
                lo = np.searchsorted(xs, xs - h1, "left")
                hi = np.searchsorted(xs, xs + h2, "right")
                counts = (hi - lo).astype(np.float64)
                cs = np.zeros(len(a) + 1, np.longdouble)
                np.cumsum(a, dtype=np.longdouble, out=cs[1:])
                s = (cs[hi] - cs[lo]).astype(np.float64)
                np.cumsum(a ** 2, dtype=np.longdouble, out=cs[1:])
                s2 = (cs[hi] - cs[lo]).astype(np.float64)
                tie = (np.abs(s2 - s ** 2 / counts) <= 1e-6 * s2) | (np.abs(a - s / counts) <= 1e-9 * a)
                if np.any(tie & (counts >= 30)):
                    vals = np.zeros(n - dd)
                    vals[xs] = a
                    kernel = np.ones(K)
                    counts = np.convolve(vals != 0, kernel, mode="same")[xs]
                    s = np.convolve(vals, kernel, mode="same")[xs]
                    s2 = np.convolve(vals ** 2, kernel, mode="same")[xs]
                local_var = (s2 - s ** 2 / counts) / (counts - 1)
                std2 = std ** 2
                np.nan_to_num(local_var, copy=False, neginf=std2, posinf=std2, nan=std2)
                local_mean = s / counts
                local_mean[counts < 30] = mean
                local_var[counts < 30] = std2
                z = np.nan_to_num((a - local_mean) / np.sqrt(local_var), nan=0, posinf=0, neginf=0)
                v[idx] = z * (1 + math.log(1 + mean, 30))
        else:
            for dd in range(4, min(D, n)):
                idx = o[bounds[dd]:bounds[dd + 1]]
                if len(idx) == 0:
                    continue
                vs = v[idx]
                v[idx] = np.nan_to_num((vs - np.mean(vs)) / np.std(vs), nan=0, posinf=0, neginf=0)
    return v


def _blocks(n: int, D: int):
    """Block starts / ends in bins (``regulator``, mustache.py:893-910)."""
    C = max(2 * D, 2000)
    if n <= C:
        return [0], [n]
    start, end = [0], [C]
    while end[-1] < n:
        start.append(end[-1] - D)
        end.append(start[-1] + C)
    end[-1] = n
    start[-1] = end[-1] - C
    return start, end


def _sigmas(o: float):
    """The 12 (sigma, truncate) of an octave, as mustache.py:714-750 computes them."""
    out = []
    for k in range(1, _SCALES + 3):
        sigma = o if k == 1 else o * 2 ** ((k - 1) / _SCALES)
        w = 2 * math.ceil(2 * sigma) + 1
        out.append((sigma, (((w - 1) / 2) - 0.5) / sigma))
    return out


@lru_cache(maxsize=4)
def _off_band(C: int, D: int) -> np.ndarray:
    """The block cells Mustache sets to 2: distance <= 4 (``tril_indices(4)``) or > D."""
    d = np.arange(C)[None, :] - np.arange(C)[:, None]
    return (d <= 4) | (d > D)


def _reflect(i, C):
    """Indices of scipy's 'reflect' extension of 0..C-1 (one reflection; anything further is clipped)."""
    i = np.where(i < 0, -1 - i, i)
    return np.clip(np.where(i >= C, 2 * C - 1 - i, i), 0, C - 1)


def _band_dog(cm, D, R, r, c, sigs):
    """The block's scale space on its band only: for each octave in turn (``sigs``: its 12 (sigma, truncate))
    yields (L_k, 3x3 max of L_k) at the tested pixels (r, c) (row-major, as Mustache's ``c[nz]``),
    k = 1..11, each value the one scipy gives on the full block (``gaussian_filter``, reflected edges;
    ``maximum_filter`` 3x3, constant 0). A Gaussian an octave shares with a later one (the last two of an
    octave are the first two of the next) is computed once.

    ``gaussian_filter`` is a pass along rows i (axis 0) then along columns j. Here the first pass runs on
    the band sheared to (j, j - i) — each of its lines is a column of the block read upwards — and the
    second on its result sheared to (i, j - i); with symmetric weights every sum is the same expression."""
    from numpy.lib.stride_tricks import as_strided
    from scipy.ndimage import gaussian_filter1d
    C = cm.shape[0]
    kg0, kg1 = 2, D + 3                      # distances where G is kept (the 3x3 maxima need 3..D+1)
    kh0 = kg0 - R                            # first distance of the first pass's output in use
    k0, j0 = kg0 - 2 * R, kg0 - R            # the sheared input: distances k0.., columns j0..
    j = np.arange(j0, C + kg1 + R - 1)[:, None]
    k = np.arange(k0, kg1 + 2 * R)[None, :]
    T = cm[_reflect(j - k, C), _reflect(j, C)]
    del j, k
    KH = (kg1 + R) - kh0
    tail = C - kg1                                       # rows reaching beyond the block (maximum_filter's 0)
    out = np.arange(tail, C)[:, None] + np.arange(kg0, kg1)[None, :] >= C
    at_L = r * (kg1 - kg0) + (c - r - kg0)               # the tested pixels in L (distances 2..D+2)
    at_M = r * (D - 3) + (c - r - 4)                     # ... and in the maxima (distances 4..D)

    def gauss(s, tr):
        H = gaussian_filter1d(T, s, axis=1, truncate=tr)
        rs, cs = H.strides
        S = as_strided(H[kh0 - j0:, kh0 - k0:], shape=(C, KH), strides=(rs, rs + cs))
        return gaussian_filter1d(S, s, axis=1, truncate=tr)[:, R:R + kg1 - kg0]

    shared = {}
    for o, sig in enumerate(sigs):
        later = {key for s2 in sigs[o + 1:] for key in s2}
        prev = None
        for key in sig:
            G = shared[key] if key in shared else gauss(*key)
            if key in later:
                shared[key] = G
            if prev is not None:
                L = prev - G
                L[tail:][out] = 0
                A = np.maximum(L[:, :-2], L[:, 1:-1])
                np.maximum(A, L[:, 2:], out=A)                                          # distances 3..D+1
                Mx = A[:, 1:-1].copy()                                                  # distances 4..D
                np.maximum(Mx[1:], A[:-1, 2:], out=Mx[1:])                              # row above, k + 1
                np.maximum(Mx[:-1], A[1:, :-2], out=Mx[:-1])                            # row below, k - 1
                Mx[0] = np.maximum(Mx[0], 0)                                            # rows -1 and C: 0
                Mx[-1] = np.maximum(Mx[-1], 0)
                yield np.take(L, at_L), np.take(Mx, at_M)
            prev = G


def _bh(p):
    """statsmodels ``multipletests(p, method='fdr_bh')[1]`` (the same arithmetic)."""
    o = np.argsort(p)
    ps = p[o]
    n = len(ps)
    q = ps / (np.arange(1, n + 1) / float(n))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q[q > 1] = 1
    out = np.empty_like(q)
    out[o] = q
    return out


def _block(x, y, v, start, end, mask, D, octave_values, pt, st):
    """``mustache()`` (mustache.py:697) on one block, then ``process_block``'s overlap mask: loop rows
    (X, Y, P, FDR, scale) in chromosome bins."""
    from scipy.ndimage import label
    from scipy.special import expm1
    empty = (np.zeros(0, np.int64),) * 2 + (np.zeros(0),) * 3
    C = max(2 * D, 2000)
    i0, i1 = np.searchsorted(x, [start, end])
    xs, ys, vs = x[i0:i1], y[i0:i1], v[i0:i1]
    f = ys < end
    xs, ys, vs = xs[f] - start, ys[f] - start, vs[f]
    t = (ys - xs >= 4) & (vs != 0)                     # mustache's nz, row-major
    r, c = xs[t], ys[t]
    nnz = len(r)
    if nnz < 10000:                  # mustache.py:701 (< 50) and :775 (``len(pFound)`` < 10000): no loops
        return empty
    cm = np.zeros((C, C))
    cm[xs, ys] = vs
    cm[_off_band(C, D)] = 2
    pAll = np.full(nnz, 2.0)
    praw = np.full(nnz, np.nan)
    Scales = np.ones(nnz)
    vAll = np.zeros(nnz)
    sigs = [_sigmas(o) for o in octave_values]
    R = max(int(t * s + 0.5) for sig in sigs for s, t in sig)     # scipy's radius for (sigma, truncate)
    Ls = _band_dog(cm, D, R, r, c, sigs)
    for sig in sigs:                                   # mustache.py:714-773
        (Lp, mp), (Lc, mc) = next(Ls), next(Ls)
        for i in range(3, _SCALES + 2):
            Ln, mn = next(Ls)
            a = np.abs(Lc)
            loc = float(a.min())                       # expon.fit (ML) and 1 - expon.cdf, as scipy does them
            scale = float(a.mean()) - loc
            pval = 1 - (-expm1(-((a - loc) / scale))) if scale > 0 else np.full(nnz, np.nan)
            up = (Lc > vAll) & (Lc == mc) & ((Lp == mp) | (Ln == mn)) & (Lc > mp) & (Lc > mn)
            vAll[up] = Lc[up]
            Scales[up] = sig[i - 1][0]
            pAll[up] = pval[up]
            praw[up] = pval[up]
            (Lp, mp), (Lc, mc) = (Lc, mc), (Ln, mn)
    found = pAll != 2
    pAll[found] = _bh(pAll[found])
    cand = np.flatnonzero(pAll < pt)
    X, Y = r[cand], c[cand]
    sc = Scales[cand]
    # sparsity (mustache.py:800-810): tested pixels in the (2s+1)^2 and (4s+1)^2 squares around each candidate,
    # Python's slicing kept: the square is cut at the block's end, and empty when it starts before row/col 0
    W = C + 2
    key = r * W + c                                    # row-major: sorted
    s = np.ceil(sc).astype(np.int64)
    keep = X != 0
    for s_, thr in ((s, st), (2 * s, 0.6)):
        off = np.arange(-int(s_.max(initial=0)), int(s_.max(initial=0)) + 1)
        rows = X[:, None] + off[None, :]
        use = (np.abs(off)[None, :] <= s_[:, None]) & (rows < C)
        lo = rows * W + (Y - s_)[:, None]
        hi = rows * W + np.minimum(Y + s_, C - 1)[:, None]
        cnt = ((np.searchsorted(key, hi, "right") - np.searchsorted(key, lo, "left")) * use).sum(1)
        cnt[(X - s_ < 0) | (Y - s_ < 0)] = 0
        keep &= ~(cnt / ((2 * s_ + 1) ** 2) < thr)
    X, Y = X[keep], Y[keep]
    if len(X) == 0:
        return empty
    # twice the diagonal's mean non-zero value (mustache.py:822)
    k = Y - X
    uk, inv = np.unique(k, return_inverse=True)
    nzd = (dg[dg != 0] for dg in (np.diagonal(cm, kk) for kk in uk))    # never empty: holds the candidate
    means = np.array([np.mean(a) for a in nzd])[inv]
    ok = cm[X, Y] > 2 * means
    X, Y = X[ok], Y[ok]
    if len(X) == 0:
        return empty
    # components of the candidates' 3x3 neighbourhoods; the lowest-FDR pixel of each (mustache.py:828-850)
    m = int(Y.max()) + 2
    lm = np.zeros((m, m), bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            lm[X + dx, Y + dy] = True
    lab, nf = label(lm, structure=np.ones((3, 3)))
    pr, pc = np.nonzero(lab)
    lb = lab[pr, pc]
    q = pr * W + pc
    j = np.minimum(np.searchsorted(key, q), nnz - 1)
    j[key[j] != q] = -1                                # not a tested pixel
    hit = j >= 0
    oq = np.where(hit, pAll[np.maximum(j, 0)], 1.0)
    order = np.lexsort((np.arange(len(pr)), oq, lb))
    first = order[np.r_[0, np.flatnonzero(np.diff(lb[order])) + 1]]
    jj = j[first]
    hj = jj >= 0
    lx, ly = pr[first] + start, pc[first] + start
    P = np.where(hj, praw[np.maximum(jj, 0)], np.nan)
    Q = oq[first]
    Sc = np.where(hj, Scales[np.maximum(jj, 0)], 1.0)
    keep = (lx >= start + mask) | (ly >= start + mask)
    return lx[keep].astype(np.int64), ly[keep].astype(np.int64), P[keep], Q[keep], Sc[keep]
