"""Hi-C loops: Juicer's HiCCUPS (``hiccups --cpu``) on a cooler, in numpy.

:func:`hiccups` calls the loops ``juicer_tools hiccups --cpu`` calls (Juicer
Tools 2.13.06, the version of the ENCODE Hi-C pipeline; the HiCCUPS code is the
same in Juicer Tools 3.0.0) with its default settings — 5, 10 and 25 kb, peak
width 4 / 2 / 1, window 7 / 5 / 3, 10% FDR, cluster radius 20 / 20 / 50 kb, the
thresholds 0.02 / 1.5 / 1.75 / 2 — on a cooler instead of a .hic. The steps,
with Juicer's source (``juicebox/tools/...``) for each:

1. **Normalised matrix** (``DatasetReaderV2.readNormalizedBlock``):
   c = float32(count / (v_i v_j)) with Juicer's divisive vector v; pixels whose
   value is NaN are dropped (so a masked bin is 0, not NaN, in the windows
   below). Cooler weights w become v = 1/w, scaled like Juicer's vectors so the
   normalised and raw matrices have the same sum (``NormalizationCalculations
   .getSumFactor``) — the scale matters to HiCCUPS through step 4's "< 16".
2. **Expected by distance** (``ExpectedValueCalculation``): normalised
   contacts per distance over all chromosomes / the possible pixels (len //
   res bins per chromosome), averaged over a window that grows until it holds
   400 contacts, divided by each chromosome's observed-to-expected factor;
   stored as float32 like a v9 .hic. Or Juicer's own (``dump expected``).
3. **Search region** (``HiCCUPSRegionHandler``): tiles of ``matrix_size`` - 40
   bins (20 bins margin); in the CPU path only tile pairs at most 8 Mb apart,
   and never the first or last 20 bins of a chromosome nor pixels within peak
   width + 2 bins of the diagonal.
4. **Local backgrounds** (``GPUController.runCPUVersion``), for pixel (i, j) at
   distance dd = j - i, window w = min(window, (dd - 1) // 2), at least peak
   width + 1, at most 20: *lower left* = rows i+1..i+w x columns j-w..j-1
   minus the peak box (rows i+1..i+p, columns j-p..j-1), grown one bin at a
   time while its normalised sum is < 16 (until w = 20 or 2w >= dd);
   *donut* = the (2w+1)^2 square minus the (2p+1)^2 peak square minus the
   cross through (i, j); *vertical* = rows i-w..i-p-1 and i+p+1..i+w in
   columns j-1..j+1; *horizontal* the transpose. The lower-left and donut
   boxes keep only cells above the diagonal; vertical / horizontal do not
   check, and reach the diagonal cell when w = dd - 1. Each background's
   expected count is (sum of c) x E(dd) / (sum of E over the same cells) x
   v_i v_j, with E summed for every cell (masked bins included). Juicer sums
   these in float32 in a fixed order; here window sums come from float64
   summed-area tables, and every pixel whose lambda chunk or window could
   change by that rounding (rigorous float32 error bound) is recomputed in
   Juicer's own order, so the outcome is Juicer's bit for bit.
5. **Lambda chunks**: chunk = floor(ln e / ln 2^0.33) (0 for e <= 1, at most
   39) per background — note 0.33, not 1/3 — while chunk k is tested against
   Poisson(2^((k+1)/3)) (``ArrayTools.generatePoissonPMF``, commons-math 2.2's
   saddle-point pmf, ported here).
6. **FDR thresholds** (``HiCCUPSUtils.calculateThresholdAndFDR``): per chunk,
   histograms of the observed counts (0..9999) of every searched pixel; the
   threshold is the last count before the first where expected tail x 1/FDR
   <= observed tail; the q-value of (chunk, count) is expected tail / observed
   tail (float32).
7. **Enriched pixels**: observed - max of the four thresholds > 0; FDR of each
   background from its table.
8. **Post-processing** (``HiCCUPSUtils.postProcess``): drop pixels with a
   NaN norm within 5 bins of either anchor; cluster greedily from the highest
   observed count (ties: larger start1, then start2) — a pixel joins when within
   the radius of the running integer centroid, the radius growing by the
   cluster's rounded spread; keep the seed if observed > 1.75 x lower-left and
   donut, > 1.5 x horizontal and vertical, > 2 x lower-left or donut, and
   either it absorbed others or its four FDRs sum to <= 0.02.
9. **Merging resolutions** (``mergeAllResolutions``): 5 kb loops with a 10 kb
   loop within 20 kb, within 110 kb of the diagonal or with observed > 100,
   plus 10 kb loops with no 5 kb loop within 20 kb, plus 25 kb loops with no
   merged loop within 50 kb. Rows come in Juicer's order — including the
   order a ``java.util.HashSet`` gives them — so the files match byte for byte.

``write`` produces Juicer's files (``merged_loops.bedpe``,
``postprocessed_pixels_<res>.bedpe``, ``enriched_pixels_<res>.bedpe``,
``fdr_thresholds_<res>``) with Java's ``Float.toString`` numbers.

Against ``juicer_tools hiccups --cpu --threads 1`` on the same matrix and SCALE vectors the four files are
byte-identical (synthetic Hi-C). With several threads Juicer's own q-values (fdrBL ... fdrV) vary between
runs by up to a few percent — ``GPUOutputContainer.updateHistograms`` is ``synchronized`` on a per-tile
object, so concurrent tiles lose histogram counts — while pixels, loops, expected values, thresholds and
clusters stay the same; on a 114M-pair HiChIP sample this port gives those same 21,413 loops.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np

from .bedpe import Pairs, open_cooler
from .genome import Genome
from .loci import Loci

MARGIN = 20                 # HiCCUPS.regionMargin
W1, W2 = 40, 10_000         # lambda chunks, largest observed count tested
KR_NEIGHBORHOOD = 5         # HiCCUPS.krNeighborhood
DEFAULTS = {5000: (4, 7, 20000), 10000: (2, 5, 20000), 25000: (1, 3, 50000)}   # peak width, window, radius
JUICER_VERSION = "2.13.06"
_F32 = np.float32
_U = 2.0 ** -24             # float32 unit roundoff
_LOGNORM = float(_F32(math.log(2.0 ** 0.33)))
_ENRICHED = ["observed", "expectedBL", "expectedDonut", "expectedH", "expectedV", "binBL", "binDonut", "binH",
             "binV", "fdrBL", "fdrDonut", "fdrH", "fdrV"]
_FINAL = ["observed", "expectedBL", "expectedDonut", "expectedH", "expectedV", "fdrBL", "fdrDonut", "fdrH", "fdrV",
          "numCollapsed", "centroid1", "centroid2", "radius"]
_CHUNK = 1 << 16          # pixels per vectorised step
_HEADER = "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor"


@dataclass
class HiCCUPS:
    """The result of :func:`hiccups`.

    ``loops`` (Juicer's ``merged_loops.bedpe``), ``pixels[res]``
    (``postprocessed_pixels_<res>.bedpe``) and ``enriched[res]``
    (``enriched_pixels_<res>.bedpe``) are :class:`~genomeblocks.Pairs` of bins,
    in Juicer's row order, with the columns ``count`` (observed raw count),
    ``expected`` (= ``expectedDonut``), ``peak`` (observed - largest threshold),
    ``expectedBL`` / ``expectedDonut`` / ``expectedH`` / ``expectedV`` (expected
    raw counts of the four backgrounds), ``binBL`` ... (lambda chunks),
    ``fdrBL`` ... (q-values), ``fdr`` (their float32 sum, the quantity tested
    against 0.02) and ``resolution``; ``pixels`` and ``loops`` add
    ``numCollapsed``, ``centroid1``, ``centroid2`` and ``radius``.
    ``thresholds[res]`` is the observed-count threshold per lambda chunk and
    background (``fdr_thresholds_<res>``); ``timings`` seconds per step."""
    loops: Pairs
    pixels: Dict[int, Pairs]
    enriched: Dict[int, Pairs]
    thresholds: Dict[int, "object"]
    params: Dict[int, dict]
    timings: dict = field(default_factory=dict)

    def __repr__(self):
        per = ", ".join(f"{r // 1000} kb: {len(p):,}" for r, p in self.pixels.items())
        return f"HiCCUPS(loops={len(self.loops):,}; postprocessed pixels {per})"

    def write(self, outdir: str, prefix: str = "") -> dict:
        """Juicer's output files in ``outdir`` (names prefixed with ``prefix``): ``merged_loops.bedpe``,
        ``postprocessed_pixels_<res>.bedpe``, ``enriched_pixels_<res>.bedpe`` and ``fdr_thresholds_<res>``, in
        juicer_tools 2.13.06's format. Returns the paths."""
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        paths = {}
        for res in self.pixels:
            p = paths[f"enriched_{res}"] = out / f"{prefix}enriched_pixels_{res}.bedpe"
            _write_bedpe(self.enriched[res], p, _ENRICHED, "0,0,0")
            p = paths[f"thresholds_{res}"] = out / f"{prefix}fdr_thresholds_{res}"
            t = self.thresholds[res]
            cols = [_jfloats(t[k].to_numpy()) for k in ("BL", "Donut", "H", "V")]
            p.write_text("".join(f"{i}\t{a}\t{b}\t{c}\t{d}\n" for i, (a, b, c, d) in enumerate(zip(*cols))))
            p = paths[f"postprocessed_{res}"] = out / f"{prefix}postprocessed_pixels_{res}.bedpe"
            _write_bedpe(self.pixels[res], p, _FINAL, "0,255,255")
        paths["merged"] = out / f"{prefix}merged_loops.bedpe"
        _write_bedpe(self.loops, paths["merged"], _FINAL, "0,255,255")
        return {k: str(v) for k, v in paths.items()}


def hiccups(cool, *, resolutions: Sequence[int] = (5000, 10000, 25000), weights="weight", divisive: bool = False,
            expected=None, chroms: Optional[Sequence[str]] = None, fdr=0.1, peak_width=None, window=None,
            cluster_radius=None, thresholds=(0.02, 1.5, 1.75, 2.0), restrict: bool = True, matrix_size: int = 512,
            workers: int = 1, verbose: bool = True) -> HiCCUPS:
    """Juicer HiCCUPS loops (``juicer_tools hiccups --cpu``) from a cooler.

    Args:
        cool: an .mcool holding ``resolutions``, or a dict {resolution: cooler URI}. Raw counts
            (``pixels/count``) are read. Juicer orders rows by its chromosome index, taken here as the
            cooler's chromosome order (that of the chrom.sizes both were built with).
        resolutions: the resolutions to call; 5, 10 and 25 kb have Juicer's defaults, others need
            ``peak_width``, ``window`` and ``cluster_radius``.
        weights: the normalisation: a bins column name (default the cooler's ``weight``), an array over the
            cooler's bins, or a dict {resolution: either}. Cooler-style (balanced = count x w_i x w_j) unless
            ``divisive``.
        divisive: the arrays are Juicer normalisation vectors (balanced = count / (v_i v_j)), e.g. from
            ``juicer_tools dump norm SCALE file.hic <chr> BP <res>``, used as they are. Cooler weights are
            inverted and scaled the way Juicer scales its vectors (same total as the raw matrix).
        expected: None to compute Juicer's normalised expected from the cooler, or a dict {resolution: {chrom:
            array}} of per-chromosome expected values by distance in bins (``juicer_tools dump expected``).
        chroms: chromosomes to search (Juicer's ``-c``); default all. The expected uses all chromosomes.
        fdr, peak_width, window, cluster_radius: Juicer's ``-f``, ``-p``, ``-i``, ``-d``; a scalar, or a dict /
            sequence per resolution.
        thresholds: Juicer's ``-t``: the FDR sum, observed / expected for horizontal and vertical, for donut and
            lower left, and for donut or lower left.
        restrict: the CPU path's search restriction (tile pairs at most 8 Mb apart); False searches whole
            chromosomes, like the GPU path.
        matrix_size: Juicer's ``-m`` (the tile size, hence the searched pixels when ``restrict``).
        workers: processes for the background sums (tiles in parallel, forked: they share the matrix); the
            result does not depend on it.
    """
    import pandas as pd

    T = {}
    res_list = [int(r) for r in resolutions]
    params = {}
    for k, res in enumerate(res_list):
        dpw, dwin, drad = DEFAULTS.get(res, (None, None, None))
        pw, win, rad = (_per_res(x, res, k, d) for x, d in ((peak_width, dpw), (window, dwin), (cluster_radius, drad)))
        if pw is None or win is None or rad is None or not 0 < pw < win < MARGIN:
            raise ValueError(f"resolution {res}: needs 0 < peak_width < window < {MARGIN} and cluster_radius "
                             f"(Juicer has defaults only for 5, 10 and 25 kb)")
        params[res] = {"peak_width": int(pw), "window": int(win), "cluster_radius": int(rad),
                       "fdr": float(_per_res(fdr, res, k, 0.1)), "region_width": matrix_size - 2 * MARGIN,
                       "restrict": bool(restrict)}
    enriched, pixels, thr_tables = {}, {}, {}
    names = None
    for res in res_list:
        P = params[res]
        t0 = time.perf_counter()
        clr = _open(cool, res)
        names = list(clr.chromnames)
        data = _read_matrix(clr, weights.get(res, "weight") if isinstance(weights, dict) else weights, divisive)
        T[f"read {res}"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        exp = _expected_from(expected[res], data) if expected is not None else _juicer_expected(data, res)
        T[f"expected {res}"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        use = [c for c in (chroms if chroms is not None else names)
               if data["v"].get(c) is not None and c in exp and len(data["v"][c]) > 2 * MARGIN]
        out = _search(data, exp, use, res, P, workers)
        cand, nflag, npix = out["cand"], out["flagged"], out["pixels"]
        T[f"backgrounds {res}"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        thr, fdr_log = _fdr_tables(out["hist"].reshape(4, W1, W2), 1.0 / P["fdr"])
        thr_tables[res] = pd.DataFrame(thr.T, columns=["BL", "Donut", "H", "V"])
        en = _enriched(data, exp, cand, thr, fdr_log, res, P)
        T[f"enriched {res}"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        pp = _postprocess(en, data, res, P["cluster_radius"], thresholds)
        T[f"postprocess {res}"] = time.perf_counter() - t0
        enriched[res], pixels[res] = en, pp
        if verbose:
            print(f"[hiccups] {res} bp: {npix:,} pixels searched ({nflag:,} re-summed in float32), "
                  f"{len(en['observed']):,} enriched, {len(pp['observed']):,} postprocessed")
        del data, cand

    t0 = time.perf_counter()
    merged = _merge_resolutions(pixels, names)
    T["merge"] = time.perf_counter() - t0
    clr = _open(cool, res_list[0])
    g = Genome.from_sizes(clr.chromsizes)
    result = HiCCUPS(_as_pairs(merged, names, g), {r: _as_pairs(p, names, g) for r, p in pixels.items()},
                     {r: _as_pairs(p, names, g) for r, p in enriched.items()}, thr_tables, params, T)
    if verbose:
        print(f"[hiccups] {len(result.loops):,} merged loops | " + ", ".join(f"{k} {v:.1f}s" for k, v in T.items()))
    return result


# ── inputs ─────────────────────────────────────────────────────────────────

def _per_res(x, res, k, default):
    if x is None:
        return default
    if isinstance(x, dict):
        return x.get(res, default)
    if isinstance(x, (list, tuple)):
        return x[k]
    return x


def _open(cool, res):
    if isinstance(cool, dict):
        clr = open_cooler(cool[res], caller="hiccups")
    else:
        try:
            clr = open_cooler(cool, res, caller="hiccups")
        except ValueError:
            clr = open_cooler(cool, caller="hiccups")
    if int(clr.binsize) != res:
        raise ValueError(f"{cool}: bin size {clr.binsize}, not {res}")
    return clr


def _read_matrix(clr, weights, divisive):
    """Intra-chromosomal pixels per chromosome (local bin ids, raw counts) and Juicer's normalisation vector per
    chromosome (float32 values in float64, NaN where invalid; None when no bin is valid — Juicer has no vector
    then and skips the chromosome)."""
    names = list(clr.chromnames)
    with clr.open("r") as h5:
        chrom_off = h5["indexes/chrom_offset"][:].astype(np.int64)
        b1_off = h5["indexes/bin1_offset"][:].astype(np.int64)
        if isinstance(weights, str):
            if weights not in h5["bins"]:
                raise ValueError(f"no bins column {weights!r} in the cooler (balance it, or pass weights=)")
            wall = h5["bins"][weights][:].astype(np.float64)
        else:
            wall = np.asarray(weights, np.float64)
            if len(wall) != chrom_off[-1]:
                raise ValueError(f"weights: {len(wall)} values for {chrom_off[-1]} bins")
        pix, v = {}, {}
        for k, c in enumerate(names):
            lo, hi = chrom_off[k], chrom_off[k + 1]
            p0, p1 = b1_off[lo], b1_off[hi]
            b1 = h5["pixels/bin1_id"][p0:p1]
            b2 = h5["pixels/bin2_id"][p0:p1]
            cnt = h5["pixels/count"][p0:p1]
            m = b2 < hi
            b1, b2, cnt = (b1[m] - lo).astype(np.int64), (b2[m] - lo).astype(np.int64), cnt[m].astype(np.float64)
            pix[c] = (b1, b2, cnt)
            w = wall[lo:hi]
            if divisive:
                vec = w.astype(_F32).astype(np.float64)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    vec = 1.0 / w
                vec[~np.isfinite(vec) | (vec <= 0)] = np.nan
                ok = ~np.isnan(vec[b1]) & ~np.isnan(vec[b2])
                if ok.any():                                    # NormalizationCalculations.getSumFactor
                    two = np.where(b1[ok] == b2[ok], 1.0, 2.0)
                    norm_sum = np.sum(two * cnt[ok] / (vec[b1[ok]] * vec[b2[ok]]))
                    vec = vec * math.sqrt(norm_sum / np.sum(two * cnt[ok]))
                vec = vec.astype(_F32).astype(np.float64)
            v[c] = vec if (np.isfinite(vec) & (vec > 0)).any() else None
    sizes = {c: int(s) for c, s in zip(names, clr.chromsizes.values)}
    return {"pix": pix, "v": v, "sizes": sizes, "names": names}


def _juicer_expected(data, res):
    """Juicer's normalised expected per chromosome (``ExpectedValueCalculation`` as ``pre`` / ``addNorm`` run
    it), float32 by distance in bins, as hiccups reads it from a v9 .hic."""
    names, sizes = data["names"], data["sizes"]
    nbins_all = max(sizes.values()) // res + 1
    actual = np.zeros(nbins_all)
    possible = np.zeros(nbins_all)
    counts = {}
    for c in names:
        vec = data["v"][c]
        if vec is None:
            continue
        b1, b2, cnt = data["pix"][c]
        x = vec.astype(_F32)[b1]
        y = vec.astype(_F32)[b2]
        ok = (x > 0) & ~np.isnan(x) & (y > 0) & ~np.isnan(y)
        if not ok.any():
            continue
        val = (cnt[ok].astype(_F32) / (x[ok] * y[ok])).astype(np.float64)    # float math, as Juicer
        actual += np.bincount(b2[ok] - b1[ok], weights=val, minlength=nbins_all)[:nbins_all]
        counts[c] = float(val.sum())
        n = sizes[c] // res
        possible[:n] += n - np.arange(n)
    if not counts:
        return {}
    nmax = max(sizes[c] // res for c in counts)
    a, p = actual.tolist(), possible.tolist()
    dens = [0.0] * nmax
    num, den, lo, hi = a[0], p[0], 0, 0
    for ii in range(nmax):                      # ExpectedValueCalculation.computeDensity, as written
        if num < 400:
            while num < 400 and hi < nmax:
                hi += 1
                num += a[hi]
                den += p[hi]
        elif num >= 400 and hi - lo > 0:
            while hi - lo > 0 and hi < nbins_all and lo < nbins_all and num - a[lo] - a[hi] >= 400:
                num = num - a[lo] - a[hi]
                den = den - p[lo] - p[hi]
                lo += 1
                hi -= 1
        dens[ii] = num / den
        if hi + 2 < nmax:
            num += a[hi + 1] + a[hi + 2]
            den += p[hi + 1] + p[hi + 2]
            hi += 2
        elif hi + 1 < nmax:
            num += a[hi + 1]
            den += p[hi + 1]
            hi += 1
    dens = np.array(dens)
    stored = dens.astype(_F32).astype(np.float64)
    out = {}
    for c in counts:
        n = sizes[c] // res
        k = min(n, nmax)
        f = float(_F32(np.sum((n - np.arange(k)) * dens[:k]) / counts[c]))
        W = len(data["v"][c])
        out[c] = (stored[np.minimum(np.arange(W), nmax - 1)] / f).astype(_F32)
    return out


def _expected_from(given, data):
    out = {}
    for c, vec in data["v"].items():
        if c in given and vec is not None:
            e = np.asarray(given[c], np.float64).astype(_F32)
            if len(e) < len(vec):
                e = np.concatenate([e, np.full(len(vec) - len(e), e[-1], _F32)])
            out[c] = e
    return out


# ── the search (GPUController.runCPUVersion) ───────────────────────────────

def _row_tiles(W, res, P):
    """The searched row tiles of a chromosome of W bins: (first row, end row, end column) each
    (``HiCCUPSRegionHandler``)."""
    rw = P["region_width"]
    ntile = -(-W // rw) + 1
    maxdiff = 8_000_000 // (rw * res) if P["restrict"] else ntile
    out = []
    for i in range(ntile):
        r0 = min(MARGIN + i * rw, W - MARGIN)
        if r0 >= W - MARGIN:
            continue
        r1 = min(r0 + rw, W - MARGIN)
        j = min(i + maxdiff, ntile - 1)
        cj = min(MARGIN + j * rw, W - MARGIN)
        while cj >= W - MARGIN and j > i:                       # the last tiles may be empty
            j -= 1
            cj = min(MARGIN + j * rw, W - MARGIN)
        out.append((r0, r1, min(cj + rw, W - MARGIN)))
    return out


class _Block:
    """Dense normalised matrix over rows [R0, R1) x columns [C0, C1): float32(count / (v_i v_j)), upper
    triangle, NaN dropped (``readNormalizedBlock`` + ``extractLocalBoundedRegion``)."""

    def __init__(self, pix, v, R0, R1, C0, C1):
        b1, b2, cnt = pix
        s, e = np.searchsorted(b1, [R0, R1])
        r, c, n = b1[s:e], b2[s:e], cnt[s:e]
        m = (c >= C0) & (c < C1)
        r, c, n = r[m], c[m], n[m]
        with np.errstate(invalid="ignore", divide="ignore"):
            val = (n / (v[r] * v[c])).astype(_F32).astype(np.float64)
        ok = ~np.isnan(val)
        self.B = np.zeros((R1 - R0, C1 - C0))
        self.B[r[ok] - R0, c[ok] - C0] = val[ok]
        self.R0, self.R1, self.C0, self.C1 = R0, R1, C0, C1


def _sat(M):
    S = np.zeros((M.shape[0] + 1, M.shape[1] + 1))
    np.cumsum(np.cumsum(M, axis=0), axis=1, out=S[1:, 1:])
    return S.ravel()


class _Strip(_Block):
    """A row tile with its 20-bin margins and the summed-area tables of its normalised counts: above the diagonal
    (``s``, for the lower-left and donut squares) and all of it (``a``: Juicer's c matrix, upper triangle with the
    diagonal, for the cross and the vertical / horizontal strips)."""

    def __init__(self, pix, v, r0, r1, c1, W):
        super().__init__(pix, v, r0 - MARGIN, min(r1 + MARGIN, W), r0 - MARGIN, min(c1 + MARGIN, W))
        self.tables()

    def tables(self):
        k = np.arange(self.C0, self.C1)[None, :] - np.arange(self.R0, self.R1)[:, None]
        self.n = self.C1 - self.C0 + 1
        self.S = {"s": _sat(np.where(k > 0, self.B, 0.0)), "a": _sat(self.B)}
        # rounding of a rectangle sum: 4 table entries, each a sum along rows then columns
        self.eps = 4 * (self.B.shape[0] + self.B.shape[1]) * np.finfo(np.float64).eps * max(self.S["a"][-1], 1.0)

    def rect(self, kind, base, a, b, e, f):
        """Sums over rows r+a .. r+b-1 x columns c+e .. c+f-1 of the pixels at ``base``."""
        S, n = self.S[kind], self.n
        return S[base + (b * n + f)] - S[base + (a * n + f)] - S[base + (b * n + e)] + S[base + (a * n + e)]


def _boxes(w, pw):
    """The four backgrounds for window ``w``, as (added rectangles, removed rectangles, bound on the float32
    operations); a rectangle is (rows a..b-1, columns e..f-1 relative to the pixel, kind): ``s`` keeps the cells
    above the diagonal (Juicer's ``i + diff - j < 0``), ``a`` all of them (no check in Juicer)."""
    return {"BL": ([(1, w + 1, -w, 0, "s")], [(1, pw + 1, -pw, 0, "s")], w * w + pw * pw),
            "Donut": ([(-w, w + 1, -w, w + 1, "s")],
                      [(-pw, pw + 1, -pw, pw + 1, "s"), (-w, -pw, 0, 1, "a"), (pw + 1, w + 1, 0, 1, "a"),
                       (0, 1, -w, -pw, "a"), (0, 1, pw + 1, w + 1, "a")], (2 * w + 1) ** 2 + (2 * pw + 1) ** 2 + 4 * w),
            "H": ([(-1, 2, -w, -pw, "a"), (-1, 2, pw + 1, w + 1, "a")], [], 6 * w),
            "V": ([(-w, -pw, -1, 2, "a"), (pw + 1, w + 1, -1, 2, "a")], [], 6 * w)}


def _dtables(d32, pw, ddmax):
    """Sums of the expected over each background, by window and distance: the expected matrix is Toeplitz, so
    they depend on (w, dd) only. Returns {background: (sum, sum of |terms|)} as (MARGIN + 1) x (ddmax + 1)."""
    d = d32.astype(np.float64)
    need = ddmax + 2 * MARGIN + 2
    if len(d) < need:
        d = np.concatenate([d, np.full(need - len(d), d[-1])])
    dd = np.arange(ddmax + 1)

    def rect_sum(a, b, e, f, kind):
        t = (np.arange(e, f)[None, :] - np.arange(a, b)[:, None]).ravel()
        if not len(t):
            return np.zeros(len(dd))
        mult = np.bincount(t - t.min())
        acc = np.zeros(len(dd))
        for i in np.flatnonzero(mult):
            k = dd + (t.min() + i)
            acc += mult[i] * (np.where(k > 0, d[np.clip(k, 0, None)], 0.0) if kind == "s" else d[np.abs(k)])
        return acc

    out = {k: (np.zeros((MARGIN + 1, ddmax + 1)), np.zeros((MARGIN + 1, ddmax + 1))) for k in ("BL", "Donut", "H", "V")}
    for w in range(pw + 1, MARGIN + 1):
        for key, (plus, minus, _) in _boxes(w, pw).items():
            p = sum((rect_sum(*r) for r in plus), np.zeros(len(dd)))
            m = sum((rect_sum(*r) for r in minus), np.zeros(len(dd)))
            out[key][0][w], out[key][1][w] = p - m, p + m
    return out


def _search(data, exp, chroms, res, P, workers=1):
    """Every searched pixel of the chromosomes: (chunk, observed) histograms, and the pixels with a count as
    candidates (chrom, row, column, observed, chunks). Tiles of all chromosomes are shared out to ``workers``
    processes."""
    pw = P["peak_width"]
    tiles = {c: _row_tiles(len(data["v"][c]), res, P) for c in chroms}
    v32 = {c: data["v"][c].astype(_F32) for c in chroms}

    dtab = {c: _dtables(exp[c], pw, max((c1 - r0 for r0, _, c1 in tiles[c]), default=0)) for c in chroms}

    def work(tasks):
        hist = np.zeros((4, W1 * W2), np.int64)
        cand = {c: [] for c in chroms}
        later = {c: [] for c in chroms}
        npix = 0
        for chrom, (r0, r1, c1) in tasks:
            pix, v, d32 = data["pix"][chrom], data["v"][chrom], exp[chrom]
            rr = np.arange(r0, r1)
            ncol = np.maximum(c1 - (rr + pw + 3), 0)
            if not ncol.any():
                continue
            R = np.repeat(rr, ncol)
            C = np.arange(len(R)) - np.repeat(np.cumsum(ncol) - ncol, ncol) + np.repeat(rr + pw + 3, ncol)
            npix += len(R)
            keep = ~np.isnan(v[R]) & ~np.isnan(v[C])    # a NaN norm makes every expected NaN: not tested
            R, C = R[keep], C[keep]
            if not len(R):
                continue
            st = _Strip(pix, v, r0, r1, c1, len(v))
            o, bins, flag = (np.concatenate(x, axis=-1) for x in zip(*(        # chunks that fit the caches
                _fast(st, R[lo:lo + _CHUNK], C[lo:lo + _CHUNK], d32, v32[chrom], dtab[chrom], pw, P["window"])
                for lo in range(0, len(R), _CHUNK))))
            if flag.any():                              # Juicer's float32 sums decide these, after the tiles
                fi = np.flatnonzero(flag)
                later[chrom].append((R[fi], C[fi], o[fi], _patches(st, R[fi], C[fi])))
                R, C, o, bins = R[~flag], C[~flag], o[~flag], bins[:, ~flag]
            _count(hist, cand[chrom], R, C, o, bins)
        flagged = 0
        for chrom, parts in later.items():
            if parts:
                R, C, o, patches = (np.concatenate(x) for x in zip(*parts))
                bins = _exact(patches, C - R, v32[chrom][R], v32[chrom][C], exp[chrom], pw, P["window"])["bins"]
                _count(hist, cand[chrom], R, C, o, bins)
                flagged += len(R)
        return hist, cand, flagged, npix

    tasks = [(c, t) for c in chroms for t in tiles[c]]
    n = max(1, min(workers, len(tasks)))
    shares = [tasks[i::n] for i in range(n)]
    parts = _fork_map(work, shares) if n > 1 else [work(shares[0])]
    cand = []
    for c in chroms:
        got = [x for p in parts for x in p[1][c]]
        if got:
            cand.append((c,) + tuple(np.concatenate(x, axis=-1) for x in zip(*got)))
    return {"cand": cand, "hist": sum(p[0] for p in parts), "flagged": sum(p[2] for p in parts),
            "pixels": sum(p[3] for p in parts)}


_FORKED = {}


def _run_forked(i):
    fn, args = _FORKED["job"]
    return fn(args[i])


def _fork_map(fn, args):
    """``[fn(a) for a in args]`` in forked worker processes (they share the parent's arrays copy-on-write; only
    the results travel back). Threads instead where fork does not exist."""
    import multiprocessing as mp
    if "fork" not in mp.get_all_start_methods():
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(len(args)) as ex:
            return list(ex.map(fn, args))
    _FORKED["job"] = (fn, args)
    try:
        with mp.get_context("fork").Pool(len(args)) as pool:
            return pool.map(_run_forked, range(len(args)))
    finally:
        _FORKED.clear()


def _count(hist, cand, R, C, o, bins):
    """Histogram the pixels with four finite expected values (``cleanUpBinNans``); keep those with a count as
    candidates (the thresholds are >= 0: no peak without a count)."""
    valid = (bins >= 0).all(axis=0)
    for t in range(4):
        hist[t] += _hist_add(bins[t, valid], o[valid])
    keep = valid & (o >= 1) & (o < W2)
    if keep.any():
        cand.append((R[keep], C[keep], o[keep], bins[:, keep].astype(np.int8)))


def _hist_add(b, o):
    """Counts of (chunk, observed) for observed < W2 (``GPUOutputContainer.updateHistograms``)."""
    m = o < W2
    return np.bincount(b[m].astype(np.int64) * W2 + o[m], minlength=W1 * W2)


def _start_window(dd, pw, window):
    w = np.minimum(window, (dd - 1) // 2)
    w = np.where(w <= pw, pw + 1, w)
    return np.minimum(w, MARGIN)


def _groups(w):
    """(window, indices) for each window size present."""
    order = np.argsort(w, kind="stable")
    vals, first = np.unique(w[order], return_index=True)
    return zip(vals.tolist(), np.split(order, first[1:]))


def _fast(st, R, C, d32, v32, dtab, pw, window):
    """Window sums for every pixel from float64 summed-area tables (counts) and Toeplitz tables (expected); the
    lambda chunks they give, and a flag where Juicer's float32 sums could give another window or chunk (those
    pixels go through :func:`_exact`)."""
    n = len(R)
    dd = C - R
    base = (R - st.R0) * st.n + (C - st.C0)
    w = _start_window(dd, pw, window)
    eps = st.eps

    def lower_left(idx, wv):
        out = st.rect("s", base[idx], 1, wv + 1, -wv, 0)
        pk = st.rect("s", base[idx], 1, pw + 1, -pw, 0)
        return out - pk, out + pk

    def near16(idx, wv):
        E, A = lower_left(idx, wv)
        return E, np.abs(E - 16) <= 1.01 * (wv * wv + pw * pw) * _U * A + eps

    E, flag = near16(np.arange(n), w)
    grow = np.flatnonzero(E < 16)
    if len(grow):
        # Juicer grows the box one bin at a time while its sum is < 16, stopping at the first w with w >= 20 or
        # 2w >= dd (``stop``). The sum only grows with w, so bisect for the first w with a sum >= 16; Juicer's
        # decisions all follow from the two around it, which are checked against the float32 error.
        w0 = w[grow]
        stop = np.maximum(w0 + 1, np.minimum(MARGIN, (dd[grow] + 1) // 2))
        lo, hi = w0.copy(), stop.copy()
        Ehi, _ = lower_left(grow, hi)
        todo = np.flatnonzero((Ehi >= 16) & (hi - lo > 1))
        while len(todo):
            mid = (lo[todo] + hi[todo]) // 2
            Em, _ = lower_left(grow[todo], mid)
            up = Em >= 16
            hi[todo[up]] = mid[up]
            lo[todo[~up]] = mid[~up]
            todo = todo[hi[todo] - lo[todo] > 1]
        w[grow] = hi
        _, f1 = near16(grow, hi)
        _, f0 = near16(grow, np.maximum(hi - 1, w0))
        flag[grow] |= f1 | f0

    k1, k2 = v32[R], v32[C]
    ddv = d32[dd]
    F = ddv.astype(np.float64) * k1.astype(np.float64) * k2.astype(np.float64)
    bins = np.empty((4, n), np.int16)
    Ec, Ac = {k: np.empty(n) for k in dtab}, {k: np.empty(n) for k in dtab}
    nops = {k: np.empty(n) for k in dtab}
    for wv, idx in _groups(w):
        b = base[idx]
        for key, (plus, minus, nop) in _boxes(wv, pw).items():
            p = sum(st.rect(r[4], b, *r[:4]) for r in plus)
            m = sum((st.rect(r[4], b, *r[:4]) for r in minus), 0.0)
            Ec[key][idx], Ac[key][idx], nops[key][idx] = p - m, p + m, nop
    for t, key in enumerate(("BL", "Donut", "H", "V")):
        Ed, Ad = dtab[key][0][w, dd], dtab[key][1][w, dd]
        e = _expected_count(Ec[key].astype(_F32), Ed.astype(_F32), ddv, k1, k2).astype(np.float64)
        # how far Juicer's float32 e can be from this one: its sums' rounding, then 5 float32 operations each
        dE = 1.01 * nops[key] * _U * Ac[key] + np.where(Ac[key] > 0, eps, 0.0)
        dEd = 1.01 * nops[key] * _U * Ad + 1e-12 * Ad
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            de = np.where(Ed == 0, 0.0, np.abs(F / Ed) * dE + np.abs(e) * (dEd / np.abs(Ed) + 12 * _U))
            x = np.log(e) / _LOGNORM                    # chunk edges at x = 1 .. 39 (e <= 1 and x < 1: chunk 0)
            pos = e - de > 0
            r = de / e
            dx = 1.01 * r / (1 - r) / _LOGNORM + 1e-12      # >= -ln(1 - r) / ln 2^0.33
            lo, hi = np.floor(x - dx), np.floor(x + dx)
            chunk = np.where(np.isfinite(e) & (e > 1), np.minimum(np.floor(x), W1 - 1), 0)
            near = np.where(pos, (hi > lo) & (hi >= 1) & (lo <= W1 - 2), e + de >= math.exp(_LOGNORM))
        bins[t] = np.where(np.isnan(e), -1, chunk)
        flag |= np.isfinite(e) & near
    with np.errstate(invalid="ignore"):
        o = np.floor(st.B.ravel()[(R - st.R0) * st.B.shape[1] + (C - st.C0)] * k1.astype(np.float64)
                     * k2.astype(np.float64) + 0.5)
    return o.astype(np.int64), bins, flag


def _expected_count(E, Ed, ddv, k1, k2):
    """``((E * d[dd]) / Ed) * kr1 * kr2`` in float32, Juicer's order."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return ((E * ddv) / Ed) * k1 * k2


def _chunk64(x):
    out = np.zeros(len(x), np.int16)
    big = np.isfinite(x) & (x > 1)
    out[big] = np.minimum(np.floor(np.log(x[big]) / _LOGNORM), W1 - 1)
    return out


def _chunk(e):
    """Lambda chunk of float32 expected values: floor(ln e / ln 2^0.33), 0 for e <= 1 or infinite, at most 39;
    -1 for NaN."""
    out = _chunk64(e.astype(np.float64))
    out[np.isnan(e)] = -1
    return out


def _patches(st, R, C):
    """The (2 x 20 + 1)^2 normalised cells around each pixel, float32 as Juicer holds them."""
    ncb = st.B.shape[1]
    k = np.arange(-MARGIN, MARGIN + 1)
    base = (R - st.R0) * ncb + (C - st.C0)
    cells = st.B.ravel()[base[:, None] + (k[:, None] * ncb + k[None, :]).ravel()]
    return cells.astype(_F32).reshape(len(R), 2 * MARGIN + 1, 2 * MARGIN + 1)


def _exact_at(st, R, C, d32, v32, pw, window):
    """:func:`_exact` for pixels of a block."""
    return _exact(_patches(st, R, C), C - R, v32[R], v32[C], d32, pw, window)


def _exact(patches, dd, k1, k2, d32, pw, window):
    """Juicer's CPU kernel for these pixels (their cells in ``patches``, their distance ``dd`` and norm values
    ``k1``, ``k2``): the same cells added in the same order in float32 (``float += double``), vectorised over the
    pixels with one numpy step per window cell."""
    n = len(dd)
    w = _start_window(dd, pw, window)

    def add(E, D, s, di, dj, didx, sign):
        val = patches[s, MARGIN + di, MARGIN + dj].astype(np.float64)
        E[s] = (E[s].astype(np.float64) + sign * val).astype(_F32)
        D[s] = D[s] + d32[didx] if sign > 0 else D[s] - d32[didx]

    Ebl, Dbl = np.zeros(n, _F32), np.zeros(n, _F32)
    wm = int(w.max())
    for di in range(1, wm + 1):                         # lower left, then its peak box out
        for dj in range(-wm, 0):
            s = np.flatnonzero((di <= w) & (-dj <= w) & (dd + dj - di > 0))
            if len(s):
                add(Ebl, Dbl, s, di, dj, dd[s] + dj - di, 1)
    for di in range(1, pw + 1):
        for dj in range(-pw, 0):
            s = np.flatnonzero(dd + dj - di > 0)
            if len(s):
                add(Ebl, Dbl, s, di, dj, dd[s] + dj - di, -1)
    active = np.flatnonzero(Ebl < 16)
    while len(active):                                  # one pass per size, the peak box removed cell by cell
        w[active] += 1
        Ebl[active] = 0
        Dbl[active] = 0
        wa = w[active]
        wm = int(wa.max())
        for di in range(1, wm + 1):
            for dj in range(-wm, 0):
                s = active[(di <= wa) & (-dj <= wa) & (dd[active] + dj - di > 0)]
                if not len(s):
                    continue
                add(Ebl, Dbl, s, di, dj, dd[s] + dj - di, 1)
                if di <= pw and -dj <= pw:
                    add(Ebl, Dbl, s, di, dj, dd[s] + dj - di, -1)
        stop = (wa >= MARGIN) | (2 * wa >= dd[active])
        active = active[~stop & (Ebl[active] < 16)]

    Edo, Ddo = np.zeros(n, _F32), np.zeros(n, _F32)
    wm = int(w.max())
    for di in range(-wm, wm + 1):                       # donut square, then the peak square out
        for dj in range(-wm, wm + 1):
            s = np.flatnonzero((abs(di) <= w) & (abs(dj) <= w) & (dd + dj - di > 0))
            if len(s):
                add(Edo, Ddo, s, di, dj, dd[s] + dj - di, 1)
    for di in range(-pw, pw + 1):
        for dj in range(-pw, pw + 1):
            s = np.flatnonzero(dd + dj - di > 0)
            if len(s):
                add(Edo, Ddo, s, di, dj, dd[s] + dj - di, -1)
    Ev, Dv, Eh, Dh = (np.zeros(n, _F32) for _ in range(4))
    for di in list(range(-wm, -pw)) + list(range(pw + 1, wm + 1)):   # the column through the pixel: vertical
        s = np.flatnonzero(abs(di) <= w)
        if len(s):
            add(Edo, Ddo, s, di, 0, np.abs(di - dd[s]), -1)
            for j in (-1, 0, 1):
                add(Ev, Dv, s, di, j, np.abs(di - dd[s] - j), 1)
    for dj in list(range(-wm, -pw)) + list(range(pw + 1, wm + 1)):   # the row through the pixel: horizontal
        s = np.flatnonzero(abs(dj) <= w)
        if len(s):
            add(Edo, Ddo, s, 0, dj, np.abs(dd[s] + dj), -1)
            for i in (-1, 0, 1):
                add(Eh, Dh, s, i, dj, np.abs(i - dd[s] - dj), 1)
    ddv = d32[dd]
    e = np.vstack([_expected_count(E, D, ddv, k1, k2) for E, D in ((Ebl, Dbl), (Edo, Ddo), (Eh, Dh), (Ev, Dv))])
    return {"e": e, "bins": np.vstack([_chunk(x) for x in e]), "w": w}


# ── thresholds ─────────────────────────────────────────────────────────────

_STIRLING = [0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
             0.0413406959554092940938221, 0.03316287351993628748511048, 0.02767792568499833914878929,
             0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
             0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
             0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
             0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637,
             0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100,
             0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
             0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867,
             0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690]


def _poisson_pmf(index: int, width: int) -> np.ndarray:
    """``ArrayTools.generatePoissonPMF``: P(X = k) for k < width, X ~ Poisson(2^((index + 1) / 3)), as
    commons-math 2.2's ``PoissonDistributionImpl.probability`` computes it (saddle-point expansion)."""
    mu = math.pow(2.0, (index + 1.0) / 3.0)
    x = np.arange(width, dtype=np.float64)
    z2 = x * x
    with np.errstate(divide="ignore", invalid="ignore"):
        stir = (0.083333333333333333333 - (0.00277777777777777777778 - (0.00079365079365079365079365 - (
            0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2) / z2) / z2) / x
    small = x < 15
    stir[small] = np.array(_STIRLING)[(2 * x[small]).astype(int)]
    dev = np.empty(width)
    ser = np.abs(x - mu) < 0.1 * (x + mu)
    with np.errstate(divide="ignore", invalid="ignore"):
        dev[~ser] = x[~ser] * np.log(x[~ser] / mu) + mu - x[~ser]
    xs = x[ser]
    d = xs - mu
    v = d / (xs + mu)
    s1 = v * d
    s = np.full(len(xs), np.nan)
    ej = 2.0 * xs * v
    v = v * v
    j = 1
    live = np.ones(len(xs), bool)
    while live.any():                                   # the series, each value until it stops changing
        s[live] = s1[live]
        ej[live] *= v[live]
        s1[live] = s[live] + ej[live] / (j * 2 + 1)
        live &= s1 != s
        j += 1
    dev[ser] = s1
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.exp(-stir - dev) / np.sqrt(2 * math.pi * x)
    p[0] = math.exp(-mu)
    return p


def _fdr_tables(hist, inv_fdr):
    """Observed-count thresholds (4 x W1, float32) and q-value tables (4 x W1 x W2, float32) from the (chunk,
    observed) histograms (``HiCCUPSUtils.calculateThresholdAndFDR``)."""
    rcs = np.cumsum(hist[:, :, ::-1], axis=2)[:, :, ::-1]
    thr = np.full((4, W1), W2, _F32)
    fdr_log = np.zeros((4, W1, W2), _F32)
    for i in range(W1):
        pmf = _poisson_pmf(i, W2).astype(_F32)
        for t in range(4):
            tot = rcs[t, i, 0]
            if tot <= 0:
                continue
            rexp = np.cumsum((_F32(tot) * pmf)[::-1], dtype=_F32)[::-1]
            ok = inv_fdr * rexp.astype(np.float64) <= rcs[t, i]
            if ok.any():
                j = int(np.argmax(ok))
                thr[t, i] = (W2 - 2) if j == 0 else (j - 1)
            pos = rcs[t, i] > 0
            stop = W2 if pos.all() else int(np.argmin(pos))
            fdr_log[t, i, :stop] = rexp[:stop] / rcs[t, i, :stop].astype(_F32)
    return thr, fdr_log


# ── enriched pixels ────────────────────────────────────────────────────────

def _enriched(data, exp, cand, thr, fdr_log, res, P):
    """Pixels whose observed count exceeds every threshold of their chunks, with Juicer's float32 expected values,
    q-values and peak, in juicer_tools' single-thread order."""
    names = data["names"]
    parts = []
    for (chrom, R, C, o, b) in cand:
        if not len(R):
            continue
        bi = b.astype(np.int64)
        peak = o - thr[np.arange(4)[:, None], bi].max(axis=0).astype(np.float64)
        m = peak > 0
        if not m.any():
            continue
        R, C, o, bi = R[m], C[m], o[m], bi[:, m]
        v = data["v"][chrom]
        v32 = v.astype(_F32)
        W = len(v)
        patches = np.zeros((len(R), 2 * MARGIN + 1, 2 * MARGIN + 1), _F32)
        for r0, r1, c1 in _row_tiles(W, res, P):        # the cells around them, tile by tile
            sel = np.flatnonzero((R >= r0) & (R < r1))
            if len(sel):
                blk = _Block(data["pix"][chrom], v, r0 - MARGIN, min(r1 + MARGIN, W), r0 - MARGIN,
                             min(c1 + MARGIN, W))
                patches[sel] = _patches(blk, R[sel], C[sel])
        ex = _exact(patches, C - R, v32[R], v32[C], exp[chrom], P["peak_width"], P["window"])
        e = ex["e"]
        if not np.array_equal(ex["bins"], bi):
            raise RuntimeError(f"hiccups: lambda chunks of {chrom} pixels changed when re-summed in float32 "
                               f"(a float error bound is too tight; please report)")
        rw = P["region_width"]
        order = np.lexsort((C, R, (C - MARGIN) // rw, (R - MARGIN) // rw))
        parts.append({"chrom": np.full(len(R), names.index(chrom)), "start1": R[order] * res,
                      "start2": C[order] * res, "res": np.full(len(R), res), "observed": o[order].astype(_F32),
                      "e": e[:, order], "bins": bi[:, order]})
    tbl = _concat(parts)
    b = tbl["bins"]
    o = tbl["observed"].astype(np.int64)
    tbl["fdr"] = np.vstack([fdr_log[t][b[t], o] for t in range(4)]) if len(o) else np.zeros((4, 0), _F32)
    tbl["peak"] = (tbl["observed"] - thr[np.arange(4)[:, None], b].max(axis=0)).astype(_F32) if len(o) \
        else np.zeros(0, _F32)
    return _take(tbl, _key_order(tbl["chrom"]))


def _concat(parts):
    keys = ("chrom", "start1", "start2", "res", "observed", "e", "bins")
    if not parts:
        return {"chrom": np.zeros(0, np.int64), "start1": np.zeros(0, np.int64), "start2": np.zeros(0, np.int64),
                "res": np.zeros(0, np.int64), "observed": np.zeros(0, _F32), "e": np.zeros((4, 0), _F32),
                "bins": np.zeros((4, 0), np.int64)}
    return {k: np.concatenate([p[k] for p in parts], axis=-1) for k in keys}


def _take(tbl, idx):
    return {k: (v[..., idx] if isinstance(v, np.ndarray) else v) for k, v in tbl.items()}


def _key_order(chrom):
    """Stable order of Juicer's Feature2DList keys ``"<i>_<i>"`` (i = chromosome index + 1) sorted as strings."""
    keys = np.array([f"{i + 1}_{i + 1}" for i in range(int(chrom.max()) + 1)] if len(chrom) else [], dtype=object)
    if not len(chrom):
        return np.zeros(0, np.int64)
    rank = np.argsort(np.argsort(keys.astype(str), kind="stable"), kind="stable")
    return np.argsort(rank[chrom], kind="stable")


# ── post-processing (HiCCUPSUtils.postProcess) ─────────────────────────────

def _postprocess(tbl, data, res, radius, thresholds):
    """Low-mappability filter, clustering to centroids, the observed / expected and FDR filters."""
    names = data["names"]
    f, t1, t2, t3 = thresholds
    keep = np.ones(len(tbl["chrom"]), bool)
    for k in np.unique(tbl["chrom"]):                   # removeLowMapQ: no NaN norm within 5 bins
        m = tbl["chrom"] == k
        v = data["v"][names[k]]
        nan = np.isnan(v)
        cnan = np.concatenate([[0], np.cumsum(nan)])
        for col in ("start1", "start2"):
            i = tbl[col][m] // res
            lo, hi = i - KR_NEIGHBORHOOD, i + KR_NEIGHBORHOOD + 1
            keep[np.flatnonzero(m)[(cnan[np.clip(hi, 0, len(v))] - cnan[np.clip(lo, 0, len(v))]) > 0]] = False
    tbl = _take(tbl, np.flatnonzero(keep))
    parts = [_coalesce(_take(tbl, np.flatnonzero(tbl["chrom"] == k)), res, radius)
             for k in _unique_in_order(tbl["chrom"])]
    out = _concat_full(parts)
    o = np.round(out["observed"].astype(np.float64))     # Math.round(float), int
    eBL, eDo, eH, eV = (out["e"][t].astype(np.float64) for t in range(4))
    fsum = ((out["fdr"][0] + out["fdr"][1]) + out["fdr"][2]) + out["fdr"][3]       # float32, left to right
    ok = ((o > t2 * eBL) & (o > t2 * eDo) & (o > t1 * eH) & (o > t1 * eV) & ((o > t3 * eBL) | (o > t3 * eDo))
          & ((out["numCollapsed"] > 1) | (fsum.astype(np.float64) <= f)))
    return _take(out, np.flatnonzero(ok))


def _unique_in_order(x):
    _, i = np.unique(x, return_index=True)
    return x[np.sort(i)]


def _concat_full(parts):
    if not parts:
        empty = _concat([])
        empty.update(fdr=np.zeros((4, 0), _F32), peak=np.zeros(0, _F32),
                     **{k: np.zeros(0, np.int64) for k in ("numCollapsed", "centroid1", "centroid2", "radius")})
        return empty
    return {k: np.concatenate([p[k] for p in parts], axis=-1) for k in parts[0]}


def _coalesce(tbl, res, radius):
    """``coalescePixelsToCentroid`` on one chromosome: seeds by decreasing observed (ties: larger start1, then
    start2); a pixel joins the seed's cluster when within the radius of the running integer centroid, visited in
    the same order; the radius grows by the cluster's rounded spread."""
    n = len(tbl["observed"])
    s1, s2 = tbl["start1"].astype(np.int64), tbl["start2"].astype(np.int64)
    order = np.lexsort((-s2, -s1, -tbl["observed"].astype(np.float64)))
    rank = np.empty(n, np.int64)
    rank[order] = np.arange(n)
    by1 = np.argsort(s1, kind="stable")
    s1_sorted = s1[by1]
    alive = np.ones(n, bool)
    seeds, ncol, cx_out, cy_out, rad_out = [], [], [], [], []
    for seed in order:
        if not alive[seed]:
            continue
        alive[seed] = False
        members = [seed]
        cx, cy = int(s1[seed]), int(s2[seed])
        r = 0.0
        R = float(radius)
        cur = rank[seed]
        while True:                                     # the next remaining pixel (in order) within R
            lo = np.searchsorted(s1_sorted, cx - R, side="left")
            hi = np.searchsorted(s1_sorted, cx + R, side="right")
            cand = by1[lo:hi]
            cand = cand[alive[cand] & (rank[cand] > cur)]
            if not len(cand):
                break
            dx = (cx - s1[cand]).astype(np.float64)
            dy = (cy - s2[cand]).astype(np.float64)
            cand = cand[np.sqrt(dx * dx + dy * dy) <= R]
            if not len(cand):
                break
            px = cand[np.argmin(rank[cand])]
            members.append(px)
            alive[px] = False
            cur = rank[px]
            m = np.array(members)
            cx = int(float(np.sum(s1[m].astype(np.float64))) / len(m))
            cy = int(float(np.sum(s2[m].astype(np.float64))) / len(m))
            dx = (cx - s1[m]).astype(np.float64)
            dy = (cy - s2[m]).astype(np.float64)
            r = float(np.floor(np.sqrt(dx * dx + dy * dy).max() + 0.5))
            R = radius + r
        seeds.append(seed)
        ncol.append(len(members))
        cx_out.append(cx + res // 2)
        cy_out.append(cy + res // 2)
        rad_out.append(int(np.floor(r + 0.5)))
    out = _take(tbl, np.array(seeds, np.int64))
    out.update(numCollapsed=np.array(ncol, np.int64), centroid1=np.array(cx_out, np.int64),
               centroid2=np.array(cy_out, np.int64), radius=np.array(rad_out, np.int64))
    return out


# ── merging resolutions (HiCCUPSUtils.mergeAllResolutions) ─────────────────

def _near_any(a, b, radius):
    """For each pixel of ``a``: is a pixel of ``b`` on the same chromosome within ``radius`` (starts, Euclidean,
    inclusive)? (``extractReproducibleCentroids``)."""
    from scipy.spatial import cKDTree
    out = np.zeros(len(a["chrom"]), bool)
    for k in np.unique(a["chrom"]):
        ia = np.flatnonzero(a["chrom"] == k)
        ib = np.flatnonzero(b["chrom"] == k)
        if not len(ib):
            continue
        pb = np.column_stack([b["start1"][ib], b["start2"][ib]]).astype(np.float64)
        pa = np.column_stack([a["start1"][ia], a["start2"][ia]]).astype(np.float64)
        dist, j = cKDTree(pb).query(pa, k=1, distance_upper_bound=radius * (1 + 1e-9) + 1)
        hit = np.isfinite(dist)
        dx = pb[j[hit], 0] - pa[hit, 0]
        dy = pb[j[hit], 1] - pa[hit, 1]
        ok = np.zeros(len(ia), bool)
        ok[hit] = np.sqrt(dx * dx + dy * dy) <= radius
        out[ia] = ok
    return out


def _merge_resolutions(lists, names):
    """5 + 10 kb, then 25 kb, as Juicer merges them; per chromosome the rows in Juicer's order, which a Java
    HashSet sets (``Feature2DList.removeDuplicates``)."""
    have5, have10, have25 = 5000 in lists, 10000 in lists, 25000 in lists
    merged = None
    if have5 or have10:
        if have5 and have10:
            five, ten = lists[5000], lists[10000]
            near5 = _near_any(five, ten, 2 * 10000)
            near10 = _near_any(ten, five, 2 * 10000)
            d = np.abs(five["start1"] - five["start2"])
            parts = [_take(five, np.flatnonzero(near5)), _take(ten, np.flatnonzero(~near10)),
                     _take(five, np.flatnonzero(d < 110000)), _take(five, np.flatnonzero(five["observed"] > 100))]
            merged = _hashset(_append_by_key(parts), names)
        else:
            merged = lists[5000] if have5 else lists[10000]
        merged = _hashset(merged, names)
    if have25:
        if merged is not None:
            t = lists[25000]
            far = ~_near_any(t, merged, 2 * 25000)
            merged = _append_by_key([merged, _take(t, np.flatnonzero(far))])
        else:
            merged = lists[25000]
        merged = _hashset(merged, names)
    if merged is None:                                  # no 5, 10 or 25 kb: all lists, unfiltered
        merged = _hashset(_append_by_key(list(lists.values())), names)
    return merged


def _append_by_key(parts):
    """``Feature2DList.add``: per chromosome, the rows of each part after those of the previous ones."""
    tbl = _concat_full([p for p in parts if len(p["chrom"])] or parts[:1])
    src = np.concatenate([np.full(len(p["chrom"]), i) for i, p in enumerate(parts)]) if parts else np.zeros(0)
    pos = np.concatenate([np.arange(len(p["chrom"])) for p in parts]) if parts else np.zeros(0)
    if not len(tbl["chrom"]):
        return tbl
    return _take(tbl, np.lexsort((pos, src, tbl["chrom"])))


def _java_string_hash(s: str) -> int:
    h = 0
    for ch in s:
        h = (31 * h + ord(ch)) & 0xFFFFFFFF
    return h


def _hashset(tbl, names):
    """``new ArrayList<>(new HashSet<>(list))`` per chromosome key, keys in Juicer's (string) order: duplicates
    dropped, rows in the HashSet's iteration order (Feature2D.hashCode = Objects.hash(chr1, end1, start1, chr2,
    end2, start2); HashMap buckets, resized while a bucket of a table smaller than 64 reaches 9 entries)."""
    import warnings
    if not len(tbl["chrom"]):
        return tbl
    out = []
    for k in np.unique(tbl["chrom"]):
        idx = np.flatnonzero(tbl["chrom"] == k)
        s1, s2 = tbl["start1"][idx].astype(np.int64), tbl["start2"][idx].astype(np.int64)
        e1, e2 = s1 + tbl["res"][idx], s2 + tbl["res"][idx]
        _, first = np.unique(np.column_stack([s1, s2, e1, e2]), axis=0, return_index=True)
        uniq = np.sort(first)
        hc = _java_string_hash(names[k])
        h = np.full(len(idx), 1, np.int64)
        for x in (np.full(len(idx), hc, np.int64), e1, s1, np.full(len(idx), hc, np.int64), e2, s2):
            xh = (x ^ (x >> 32)) & 0xFFFFFFFF                   # Long.hashCode / String.hashCode
            h = (31 * h + xh) & 0xFFFFFFFF
        spread = (h ^ (h >> 16)) & 0xFFFFFFFF
        n = len(idx)
        cap = 1 << (max(int(np.float32(n) / np.float32(0.75)) + 1, 16) - 1).bit_length()   # tableSizeFor
        while True:
            b = spread[uniq] & (cap - 1)
            fullest = np.bincount(b).max()
            if fullest < 9:
                break
            if cap < 64:
                cap *= 2
                continue
            warnings.warn("hiccups: a Java HashMap bucket would be treeified; the order of these merged loops "
                          "may differ from Juicer's (their content does not)")
            break
        out.append(idx[uniq[np.argsort(b, kind="stable")]])
    order = np.concatenate(out)
    tbl = _take(tbl, order)
    return _take(tbl, _key_order(tbl["chrom"]))


# ── output ─────────────────────────────────────────────────────────────────

def _as_pairs(tbl, names, g) -> Pairs:
    codes = np.array([g.code[c] for c in names], np.int32)[tbl["chrom"]] if len(tbl["chrom"]) \
        else np.zeros(0, np.int32)
    a = Loci(codes, tbl["start1"], tbl["start1"] + tbl["res"], genome=g)
    b = Loci(codes, tbl["start2"], tbl["start2"] + tbl["res"], genome=g)
    e, fdr = tbl["e"], tbl["fdr"]
    cols = {"count": tbl["observed"].astype(np.int64), "expected": e[1], "peak": tbl["peak"],
            "expectedBL": e[0], "expectedDonut": e[1], "expectedH": e[2], "expectedV": e[3],
            "binBL": tbl["bins"][0], "binDonut": tbl["bins"][1], "binH": tbl["bins"][2], "binV": tbl["bins"][3],
            "fdrBL": fdr[0], "fdrDonut": fdr[1], "fdrH": fdr[2], "fdrV": fdr[3],
            "fdr": ((fdr[0] + fdr[1]) + fdr[2]) + fdr[3], "resolution": tbl["res"]}
    for k in ("numCollapsed", "centroid1", "centroid2", "radius"):
        if k in tbl:
            cols[k] = tbl[k]
    return Pairs(a, b, cols)


def _write_bedpe(P: Pairs, path, keys, color):
    """A Juicer Feature2DList file (``exportFeatureList(file, true, ...)``): header, version line, rows."""
    cols = {"observed": _jfloats(P.cols["count"].astype(_F32))}
    for k in ("expectedBL", "expectedDonut", "expectedH", "expectedV", "fdrBL", "fdrDonut", "fdrH", "fdrV"):
        cols[k] = _jfloats(P.cols[k])
    for k in ("binBL", "binDonut", "binH", "binV"):
        cols[k] = _jfloats(P.cols[k].astype(_F32))
    for k in ("numCollapsed", "centroid1", "centroid2", "radius"):
        if k in P.cols:
            cols[k] = [str(int(x)) for x in P.cols[k]]
    chroms = P.a.chroms
    with open(path, "w") as fh:
        fh.write(_HEADER + "".join("\t" + k for k in keys) + "\n")
        fh.write(f"# juicer_tools version {JUICER_VERSION}\n")
        for i in range(len(P)):
            fh.write(f"{chroms[i]}\t{P.a.starts[i]}\t{P.a.ends[i]}\t{chroms[i]}\t{P.b.starts[i]}\t{P.b.ends[i]}"
                     f"\t.\t.\t.\t.\t{color}" + "".join("\t" + cols[k][i] for k in keys) + "\n")


def _jfloats(x) -> list:
    """Java's ``String.valueOf(float)`` of each value: the shortest digits that round-trip, as d.ddd for
    1e-3 <= |x| < 1e7, else d.dddE<exp>."""
    out = []
    for v in np.asarray(x, _F32):
        if np.isnan(v):
            out.append("NaN")
            continue
        if np.isinf(v):
            out.append("Infinity" if v > 0 else "-Infinity")
            continue
        if v == 0:
            out.append("-0.0" if np.signbit(v) else "0.0")
            continue
        s = np.format_float_scientific(v, unique=True, trim="-")
        mant, ex = s.split("e")
        ex = int(ex)
        sign = "-" if mant.startswith("-") else ""
        digits = mant.lstrip("-").replace(".", "")
        if -3 <= ex < 7:
            point = ex + 1
            if point <= 0:
                body = "0." + "0" * (-point) + digits
            elif point >= len(digits):
                body = digits + "0" * (point - len(digits)) + ".0"
            else:
                body = digits[:point] + "." + digits[point:]
        else:
            body = digits[0] + "." + (digits[1:] or "0") + f"E{ex}"
        out.append(sign + body)
    return out
