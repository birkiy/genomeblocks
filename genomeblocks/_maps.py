"""HiChIP / PLAC-seq loops: MAPS (Juric et al. 2019, PLoS Comput Biol 15:e1006982) in numpy.

:func:`maps` calls the bin pairs MAPS calls with its default path
(``run_pipeline.sh``: ``model=pospoisson``, ``binning_range=1000000``,
``length_cutoff=1000``, ``sex_chroms_to_process=X``, ``filter_file=None``) from
a pairs file instead of FASTQ → feather → per-chromosome text files → pandas →
R. The steps, with the lines of github.com/ijuric/MAPS (master 21222de, same
logic as github.com/HuMingLab/MAPS) they follow:

1. **Short and long pairs** (``feather_split_rongxin.py`` L104-130, L256-266):
   intra-chromosomal pairs on ``chroms`` (chr1-22 + chrX). An end's position is
   its 5' end as feather has it (pysam ``reference_start`` for a forward read,
   ``reference_end`` for a reverse one), i.e. ``pos - 1`` / ``pos`` for the
   1-based 5' positions of a pairs file. |TLEN| = the distance of the two 5'
   ends (what ``samtools fixmate`` writes). *Short*: |TLEN| <= ``short_cutoff``,
   opposite strands, distinct ends; it counts once in the bin of the midpoint of
   its two ends (``MAPS.py`` L124). *Long*: |TLEN| > ``short_cutoff``, any
   orientation; each end is binned at its 5' end (feather bins the midpoint of
   the aligned read, ``bedtools bamtobed -bedpe`` — a pairs file has no read
   span; this is the one input-level difference, see ``benchmarks/bench_maps.py``).
2. **Peak bins** (``MAPS.py`` L106-117): a peak covers bins
   ``floor(start/B) .. ceil(end/B) - 1``; MAPS's "special case" is kept: when the
   first peak of a chromosome spans exactly 4 bins, the bin before it is *not*
   a peak bin, even if another peak covers it.
3. **Bin pairs** (``MAPS.py`` L137-158): long pairs counted per (bin1 <= bin2),
   kept at 1 <= bin2 - bin1 <= ``binning_range / B``; AND = both bins are peak
   bins, XOR = exactly one is (the NOT set is dropped).
4. **Covariates** (``MAPS.py`` L168-212), per bin pair, from MAPS's genomic
   features file (effective fragment length f, GC, mappability m per bin) and
   the short-pair count s of each bin: ``logl = log((f1+1)(f2+1)/B^2)``,
   ``loggc = log(gc1 gc2)``, ``logm = log(m1 m2)``,
   ``logdist = log((1 + d) / (range/B))``, ``logShortCount = log((s1+1)(s2+1) /
   ((max s1 + 1)(max s2 + 1)))`` (maxima over the set's pairs of that
   chromosome, distance-1 pairs included); pairs with f = 0 or without a
   features row are dropped.
5. **Regression** (``MAPS_regression_and_peak_caller.r`` L88, L104-128): pairs
   at distance > 1 bin; per chromosome and set, VGAM ``vglm(count ~ logl +
   loggc + logm + logdist + logShortCount, pospoisson())`` — zero-truncated
   Poisson, log link — fitted by VGAM's own Fisher scoring (same start
   ``0.95 (mean y + 1/16) + 0.05 y``, working weights, log-likelihood
   convergence test and half-stepping; ``vglm.fit``, ``pospoisson``, identical
   in VGAM 1.1-2, the version MAPS lists, and 1.1-14). The log-likelihood is
   computed in log space, as VGAM 1.1-2's ``dpospois(log = TRUE)`` does; VGAM
   >= 1.1-6 computes it with ``dgaitdpois``, which underflows to -Inf when a
   count is far above its starting mean (counts of ~400 at the start's ~30),
   and ``vglm`` then stops with "missing value where TRUE/FALSE needed" — MAPS
   fails that way on deep samples with a current VGAM. Expected = the
   truncated mean ``l / (1 - e^-l)``; p =
   ``ppois(c, E, lower=FALSE) / ppois(0, E, lower=FALSE)`` — P(Y > c) with the
   *truncated mean* as Poisson rate, as coded. Pairs with p <= 1/n are dropped
   and the model refitted; its coefficients rounded to 10 decimals give
   ``expected`` (rounded to 10 decimals, then truncated-mean transformed),
   ``ratio``, ``P`` and BH ``fdr`` within the chromosome x set (NA p, from an
   expected that rounds to 0, are left out as ``p.adjust`` does).
6. **Calls** (L73-75, L309-310): count >= 12, ratio >= 2, -log10 fdr > 2.
7. **Clusters** (``label_peaks`` L171-227, ``classify_peaks`` L229-261): calls
   of a chromosome within ``gap`` (15 kb) in both bin starts are neighbours;
   calls without one are singletons, the rest are connected components
   numbered by their first row; a cluster's -log10 P is the sum over its calls
   (p = 0 replaced by the chromosome's smallest p); clusters above the knee of
   the sorted sums are ``BroadPeak``, below ``SharpPeak``.
8. **Output** (``MAPS_peak_formatting.r``): singletons need fdr < 1e-4;
   the summit of a cluster is its call(s) of smallest fdr (as 15-digit text);
   :meth:`MAPS.write` writes ``{prefix}.2.sig3Dinteractions.bedpe`` in R's
   format (``write.table``, ``scipen = 999``, 15 significant digits).

Row order (and with it cluster numbering) follows MAPS.py under pandas >= 2.2
(inner merges keep the left order): AND pairs by (bin1, bin2); XOR pairs with
the peak at bin1 first, then those with the peak at bin2.

Not ported: feather itself (mapping, MAPQ, duplicates), MAPS's negative-binomial
model, the blacklist (``filter_file``) and the summary/QC text files.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .bedpe import PAIRS_FORMAT_COLUMNS, Pairs, _detect_pairs_format, open_cooler
from .genome import Genome
from .loci import Loci
from .loops import _bh, _q15

# strand columns (0-based) of the pairs layouts genomeblocks.bedpe knows
PAIRS_STRAND_COLUMNS = {"allvalidpairs": (3, 6), "pairs": (5, 6), "juicer": (1, 5)}
_SETS = ("AND", "XOR")


@dataclass
class MAPS:
    """The result of :func:`maps`.

    ``tested``, ``significant`` and ``loops`` are :class:`~genomeblocks.Pairs`
    of bins (anchors ``[bin start, bin start + B)``).

    ``tested``: every bin pair of the regressions (distance > 1 bin), per
    chromosome AND then XOR, with ``set`` ("AND" / "XOR"), ``peak1``,
    ``peak2``, ``count``, ``dist`` (bins), ``effective_length1/2``, ``gc1/2``,
    ``mappability1/2``, ``short_count1/2``, the covariates ``logl``, ``loggc``,
    ``logm``, ``logdist``, ``logShortCount``, ``expected1`` and ``P1`` (first
    fit), ``expected`` (MAPS ``expected2``), ``ratio`` (``ratio2``), ``P``
    (``p_val_reg2``), ``fdr`` and ``p_bonferroni``.

    ``significant``: the calls (MAPS's ``.peaks`` table) with ``count``,
    ``expected``, ``ratio``, ``P``, ``fdr``, ``set``, ``CountNei``, ``label``
    (0 = singleton), ``NegLog10P``, ``ClusterSize``, ``ClusterType``.

    ``loops``: the rows of MAPS's ``sig3Dinteractions.bedpe`` — the calls minus
    singletons with fdr >= 1e-4 — with ``count``, ``expected``, ``fdr``,
    ``P``, ``ratio``, ``ClusterLabel``, ``ClusterSize``, ``ClusterType``,
    ``ClusterNegLog10P``, ``ClusterSummit``.

    ``fits``: one row per chromosome x set (rows, refit rows, IRLS iterations,
    first and second-pass coefficients). ``qc``: MAPS's ``maps.qc`` sums
    (AND / XOR / NOT long pairs per chromosome)."""
    tested: Pairs
    significant: Pairs
    loops: Pairs
    fits: "object"
    qc: "object"
    binsize: int
    fdr_label: str = "2"
    timings: dict = field(default_factory=dict)

    def __repr__(self):
        return (f"MAPS(tested={len(self.tested):,}, significant={len(self.significant):,}, "
                f"loops={len(self.loops):,}, binsize={self.binsize})")

    def write(self, outdir: str, prefix: str, *, tested: bool = False) -> dict:
        """MAPS's main output ``{prefix}.{fdr}.sig3Dinteractions.bedpe`` (``MAPS_peak_formatting.r``; MAPS
        names it ``{dataset}.{B/1000}k.2.sig3Dinteractions.bedpe``), byte for byte as R writes it, and with
        ``tested=True`` every tested bin pair as parquet. Returns the paths."""
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        paths = {"bedpe": out / f"{prefix}.{self.fdr_label}.sig3Dinteractions.bedpe"}
        _write_bedpe(self.loops, self.binsize, paths["bedpe"])
        if tested:
            paths["tested"] = out / f"{prefix}.tested.parquet"
            self.tested.save(str(paths["tested"]))
        return {k: str(v) for k, v in paths.items()}


def maps(pairs_or_cool, peaks, *, features, resolution: int = 5000, binning_range: int = 1_000_000,
         chroms: Optional[Sequence[str]] = None, short=None, short_cutoff: int = 1000, format: str = "auto",
         columns: Optional[Tuple[int, int, int, int]] = None, strand_columns: Optional[Tuple[int, int]] = None,
         chunksize: int = 5_000_000, count_cutoff: int = 12, ratio_cutoff: float = 2.0,
         fdr_cutoff: float = 2.0, gap: int = 15_000, singleton_fdr: float = 1e-4,
         verbose: bool = True) -> MAPS:
    """MAPS loops (AND / XOR zero-truncated Poisson regression) from a pairs file or a cooler and 1D peaks.

    Args:
        pairs_or_cool: a pairs file — HiC-Pro allValidPairs, 4DN / pairtools ``.pairs``, or the Juicer-like
            layout ``readID strand1 chr1 pos1 frag1 strand2 chr2 pos2 frag2`` (strands ``+``/``-`` or
            ``0``/``16``); 1-based 5' positions — or a cooler (.cool, URI, .mcool with ``resolution``), whose
            raw counts are then the long-pair counts (it holds short pairs too, which only matters for
            distance-1 pairs, i.e. for the short-count maxima of step 4).
        peaks: the 1D peaks (anything :func:`~genomeblocks.as_loci` reads; file order matters for MAPS's
            first-peak quirk).
        features: MAPS's genomic features table for this genome, enzyme and resolution (tab-separated
            ``chr start end effective_length gc mappability``, no header), as a path or a frame; e.g.
            ``MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_5Kb_el.GRCh38.txt`` of github.com/ijuric/MAPS.
        resolution: bin size B in bp.
        binning_range: largest bin-pair distance (bp).
        chroms: chromosomes analysed (MAPS: ``chr1..chr{N_CHROMS}`` + sex chromosomes); default chr1-22, chrX.
        short: the short-pair count of each bin when it is not derived from ``pairs_or_cool``: a pairs file
            (its short pairs are counted), or ``{chrom: counts per bin}``. Required for a cooler.
        short_cutoff: largest |TLEN| of a short pair.
        format, columns: pairs layout (``"auto"``, ``"allvalidpairs"``, ``"pairs"``, ``"juicer"``) or the
            0-based (chrom1, pos1, chrom2, pos2) columns; ``strand_columns``: the 0-based strand columns
            (needed with ``columns``).
        count_cutoff, ratio_cutoff, fdr_cutoff: calls need count >= ``count_cutoff``, count / expected >=
            ``ratio_cutoff`` and -log10 fdr > ``fdr_cutoff``.
        gap: clustering distance (bp, both anchors).
        singleton_fdr: singletons are reported when their fdr is below this.
    """
    import pandas as pd

    T = {}
    B = int(resolution)
    chroms = list(chroms) if chroms is not None else [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    cidx = {c: i for i, c in enumerate(chroms)}
    rng = float(binning_range) / float(B)                 # MAPS's BIN_RANGE (a float)
    max_d = int(np.floor(rng))

    # ── genomic features ──
    t0 = time.perf_counter()
    feat = _read_features(features, chroms, B)
    T["features"] = time.perf_counter() - t0

    # ── long pairs per bin pair, short pairs per bin ──
    t0 = time.perf_counter()
    if _is_cooler(pairs_or_cool):
        if short is None:
            raise ValueError("a cooler holds no short pairs: pass short= (a pairs file or {chrom: counts})")
        lng, totals = _cooler_long(pairs_or_cool, resolution, chroms, max_d)
        shc = None
    else:
        lng, shc, totals = _scan_pairs(pairs_or_cool, chroms, B, short_cutoff, max_d, format, columns,
                                       strand_columns, chunksize, want_short=short is None, verbose=verbose)
    if short is not None:
        if isinstance(short, dict):
            shc = {c: np.asarray(short[c], np.float64) for c in chroms if c in short}
        else:
            _, shc, _ = _scan_pairs(short, chroms, B, short_cutoff, max_d, format, columns, strand_columns,
                                    chunksize, want_short=True, want_long=False, verbose=verbose)
    T["pairs"] = time.perf_counter() - t0

    # ── peak bins, AND / XOR sets, covariates ──
    t0 = time.perf_counter()
    from .interop import as_loci
    g = Genome(chroms)
    pk = as_loci(peaks, genome=g)
    pk_chrom = np.asarray(pk.chroms, dtype=object)
    rows, qc = [], []
    for c in chroms:
        b1, dd, cnt = lng.get(c, (np.zeros(0, np.int64),) * 3)
        m = dd >= 1                                         # MAPS.py L145-151: 1 <= |bin2 - bin1| <= range
        b1, dd, cnt = b1[m], dd[m], cnt[m]
        if not len(b1):
            continue
        b2 = b1 + dd
        sel = pk_chrom == c
        nb = int(max(b2.max(), feat[c][0].size if c in feat else 0, 0)) + 2
        is_peak = _peak_bins(pk.starts[sel], pk.ends[sel], B, nb)
        p1, p2 = is_peak[b1], is_peak[b2]
        s_and = p1 & p2
        and_i = np.flatnonzero(s_and)
        xor_i = np.r_[np.flatnonzero(p1 & ~p2), np.flatnonzero(~p1 & p2)]
        qc.append((c, int(cnt[s_and].sum()), int(cnt[xor_i].sum()),
                   int(totals.get(c, 0)) - int(cnt[s_and].sum()) - int(cnt[xor_i].sum())))
        s = shc.get(c, np.zeros(0)) if shc is not None else np.zeros(0)
        for name, idx in (("AND", and_i), ("XOR", xor_i)):
            r = _covariates(b1[idx], b2[idx], cnt[idx], p1[idx], p2[idx], feat.get(c), s, B, rng)
            if r is not None:
                r["chrom"] = c
                r["set"] = name
                rows.append(r)
    T["covariates"] = time.perf_counter() - t0

    # ── per chromosome x set: positive-Poisson fit, refit, test ──
    t0 = time.perf_counter()
    keep = [r for r in rows if len(r["count"])]
    for r in keep:                                          # R: subset(mm, dist > 1)
        m = r["dist"] > 1
        for k in list(r):
            if isinstance(r[k], np.ndarray):
                r[k] = r[k][m]
    keep = [r for r in keep if len(r["count"])]
    dataset_length = sum(len(r["count"]) for r in keep)
    fits = []
    for r in keep:
        fits.append(_regress(r, dataset_length))
    T["regression"] = time.perf_counter() - t0

    # ── tables, calls, clusters ──
    t0 = time.perf_counter()
    cols = ["peak1", "peak2", "count", "dist", "effective_length1", "gc1", "mappability1", "short_count1",
            "effective_length2", "gc2", "mappability2", "short_count2", "logl", "loggc", "logm", "logdist",
            "logShortCount", "expected1", "P1", "expected", "ratio", "P", "p_bonferroni", "fdr"]
    cat = {k: (np.concatenate([r[k] for r in keep]) if keep else np.zeros(0)) for k in cols}
    chrom_a = np.concatenate([np.full(len(r["count"]), cidx[r["chrom"]], np.int32) for r in keep]) if keep \
        else np.zeros(0, np.int32)
    set_a = np.concatenate([np.full(len(r["count"]), r["set"], dtype=object) for r in keep]) if keep \
        else np.zeros(0, object)
    s1 = np.concatenate([r["bin1"] for r in keep]) * B if keep else np.zeros(0, np.int64)
    s2 = np.concatenate([r["bin2"] for r in keep]) * B if keep else np.zeros(0, np.int64)
    cat["count"] = cat["count"].astype(np.int64)
    cat["dist"] = cat["dist"].astype(np.int64)
    tested = Pairs(Loci(chrom_a, s1, s1 + B, genome=g), Loci(chrom_a, s2, s2 + B, genome=g),
                   {"set": set_a, **cat})
    with np.errstate(divide="ignore"):
        call = (cat["count"] >= count_cutoff) & (cat["ratio"] >= ratio_cutoff) & (-np.log10(cat["fdr"]) > fdr_cutoff)
    significant = _label(tested, np.flatnonzero(call), chroms, B, gap)
    loops = _format(significant, singleton_fdr)
    T["calls"] = time.perf_counter() - t0

    fits_df = pd.DataFrame(fits)
    qc_df = pd.DataFrame(qc, columns=["chrom", "AND_set", "XOR_set", "NOT_set"])
    fl = f"{fdr_cutoff:g}"
    res = MAPS(tested, significant, loops, fits_df, qc_df, B, fl, T)
    if verbose:
        print(f"[maps] {len(tested):,} bin pairs tested, {len(significant):,} significant, {len(loops):,} reported | "
              + ", ".join(f"{k} {v:.1f}s" for k, v in T.items()))
    return res


# ── inputs ─────────────────────────────────────────────────────────────────

def _is_cooler(x) -> bool:
    s = str(x)
    return "::" in s or s.endswith((".cool", ".mcool"))


def _read_features(features, chroms, B) -> Dict[str, tuple]:
    """{chrom: (has_row, effective_length, gc, mappability)} indexed by bin (``start // B``)."""
    import pandas as pd
    if isinstance(features, (str, Path)):
        df = pd.read_csv(features, sep="\t", header=None)
    else:
        df = pd.DataFrame(features).copy()
    df = df.iloc[:, :6]
    df.columns = ["chr", "start", "end", "effective_length", "gc", "mappability"]
    out = {}
    for c, d in df[df["chr"].isin(chroms)].groupby("chr", sort=False):
        b = d["start"].to_numpy(np.int64) // B
        n = int(b.max()) + 1
        has = np.zeros(n, bool)
        el, gc, mp = np.zeros(n, np.int64), np.zeros(n), np.zeros(n)
        has[b] = True
        el[b] = d["effective_length"].to_numpy(np.int64)
        gc[b] = d["gc"].to_numpy(np.float64)
        mp[b] = d["mappability"].to_numpy(np.float64)
        out[c] = (has, el, gc, mp)
    return out


def _pairs_reader(filename, format, columns, strand_columns, chunksize):
    """Chunks of (chrom1, pos1, chrom2, pos2, strand1, strand2) as pandas frames (categorical chroms / strands)."""
    import pandas as pd
    if columns is None:
        if format == "auto":
            format = _detect_pairs_format(filename)
        if format not in PAIRS_FORMAT_COLUMNS:
            raise ValueError(f"Unknown pairs format: {format!r}. Known: {sorted(PAIRS_FORMAT_COLUMNS)}, "
                             f"or pass columns=(c1, p1, c2, p2) and strand_columns=(s1, s2).")
        columns = PAIRS_FORMAT_COLUMNS[format]
        if strand_columns is None:
            strand_columns = PAIRS_STRAND_COLUMNS[format]
    elif strand_columns is None:
        raise ValueError("MAPS needs strands: pass strand_columns=(s1, s2) with columns=")
    use = list(columns) + list(strand_columns)
    names = ["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
    dtype = {use[0]: "category", use[1]: "int64", use[2]: "category", use[3]: "int64",
             use[4]: "category", use[5]: "category"}
    rename = dict(zip(use, names))
    for ch in pd.read_csv(filename, sep="\t", header=None, comment="#", usecols=use, dtype=dtype,
                          chunksize=chunksize, engine="c"):
        yield ch.rename(columns=rename)


def _codes(cat_series, lut: dict, default=-1) -> np.ndarray:
    cats = cat_series.cat.categories
    m = np.array([lut.get(str(x), default) for x in cats], np.int64)
    codes = cat_series.cat.codes.to_numpy()
    out = np.full(len(codes), default, np.int64)
    ok = codes >= 0
    out[ok] = m[codes[ok]]
    return out


def _scan_pairs(filename, chroms, B, short_cutoff, max_d, format, columns, strand_columns, chunksize,
                want_short=True, want_long=True, verbose=True):
    """Long pairs per (bin1, distance) and short pairs per bin, per chromosome, streaming the file."""
    cidx = {c: i for i, c in enumerate(chroms)}
    rev_lut = {"-": 1, "16": 1, "+": 0, "0": 0}
    nbmax = 1 << 22                                        # bins per chromosome in the key
    keys, kcnt, short_keys = [], [], []
    totals = np.zeros(len(chroms), np.int64)
    n = 0
    for ch in _pairs_reader(filename, format, columns, strand_columns, chunksize):
        n += len(ch)
        c1 = _codes(ch["chrom1"], cidx)
        c2 = _codes(ch["chrom2"], cidx)
        ok = (c1 >= 0) & (c1 == c2)
        r1 = _codes(ch["strand1"], rev_lut)[ok]
        r2 = _codes(ch["strand2"], rev_lut)[ok]
        if (r1 < 0).any() or (r2 < 0).any():
            raise ValueError("strands must be +/- or 0/16")
        c = c1[ok]
        p1 = ch["pos1"].to_numpy()[ok] - (r1 == 0)          # feather's 5' ends (0-based start / end)
        p2 = ch["pos2"].to_numpy()[ok] - (r2 == 0)
        tlen = np.abs(p1 - p2)
        if want_short:
            s = (tlen <= short_cutoff) & (r1 != r2) & (p1 != p2)
            short_keys.append(c[s] * nbmax + (p1[s] + p2[s]) // (2 * B))
        if want_long:
            lo = tlen > short_cutoff
            a, b, cl = p1[lo] // B, p2[lo] // B, c[lo]
            totals += np.bincount(cl, minlength=len(chroms))
            b1 = np.minimum(a, b)
            d = np.abs(a - b)
            k = d <= max_d
            u, cn = np.unique((cl[k] * nbmax + b1[k]) * (max_d + 1) + d[k], return_counts=True)
            keys.append(u)
            kcnt.append(cn)
        if verbose:
            print(f"\r[maps] {n:,} pairs read", end="", flush=True)
    if verbose:
        print()
    lng = {}
    if want_long and keys:
        u, inv = np.unique(np.concatenate(keys), return_inverse=True)
        cnt = np.bincount(inv, weights=np.concatenate(kcnt)).astype(np.int64)
        d = u % (max_d + 1)
        rest = u // (max_d + 1)
        cc, b1 = rest // nbmax, rest % nbmax
        for i, cname in enumerate(chroms):
            m = cc == i
            if m.any():
                lng[cname] = (b1[m], d[m], cnt[m])
    shc = None
    if want_short:
        sk = np.concatenate(short_keys) if short_keys else np.zeros(0, np.int64)
        shc = {}
        for i, cname in enumerate(chroms):
            v = sk[sk // nbmax == i] % nbmax
            if len(v):
                shc[cname] = np.bincount(v).astype(np.float64)
    return lng, shc, {c: int(totals[i]) for i, c in enumerate(chroms)}


def _cooler_long(cool, resolution, chroms, max_d):
    """Raw pixel counts per (bin1, distance) within each chromosome (distance <= max_d bins)."""
    clr = open_cooler(cool, resolution if not str(cool).endswith(".cool") and "::" not in str(cool) else None,
                      caller="maps")
    lng, totals = {}, {}
    for c in chroms:
        if c not in clr.chromnames:
            continue
        off = clr.offset(c)
        px = clr.matrix(balance=False, as_pixels=True, join=False).fetch(c)
        b1 = px["bin1_id"].to_numpy(np.int64) - off
        b2 = px["bin2_id"].to_numpy(np.int64) - off
        cnt = px["count"].to_numpy(np.int64)
        totals[c] = int(cnt.sum())
        d = b2 - b1
        k = (d <= max_d) & (cnt > 0)
        o = np.lexsort((b2[k], b1[k]))
        lng[c] = (b1[k][o], d[k][o], cnt[k][o])
    return lng, totals


def _peak_bins(starts, ends, B, nb) -> np.ndarray:
    """MAPS.py L106-117: bins floor(s/B) .. ceil(e/B) - 1 of every peak; the first-peak quirk."""
    is_peak = np.zeros(nb, bool)
    if not len(starts):
        return is_peak
    sb = np.floor(np.asarray(starts) / B).astype(np.int64)
    eb = np.ceil(np.asarray(ends) / B).astype(np.int64)
    special = eb[0] - sb[0] == 4                           # MACS2.shape[1] - 1 with 5 columns
    if special:
        sb = sb.copy()
        sb[0] -= 1
    ln = np.maximum(eb - sb, 0)
    idx = np.repeat(sb - np.r_[0, np.cumsum(ln)[:-1]], ln) + np.arange(ln.sum())
    idx = idx[(idx >= 0) & (idx < nb)]
    is_peak[idx] = True
    if special and 0 <= sb[0] < nb:
        is_peak[sb[0]] = False
    return is_peak


def _covariates(b1, b2, cnt, p1, p2, feat, short, B, rng):
    """MAPS.py L168-212 for one chromosome x set (rows in MAPS order)."""
    if feat is None:
        return None
    has, el, gc, mp = feat
    nf = len(has)
    ok = (b1 < nf) & (b2 < nf)
    ok[ok] = has[b1[ok]] & has[b2[ok]]
    ok[ok] = (el[b1[ok]] > 0) & (el[b2[ok]] > 0)
    b1, b2, cnt, p1, p2 = b1[ok], b2[ok], cnt[ok], p1[ok], p2[ok]

    def sc(b):
        out = np.zeros(len(b))
        m = b < len(short)
        out[m] = short[b[m]]
        return out
    s1, s2 = sc(b1), sc(b2)
    el1, el2 = el[b1], el[b2]
    gc1, gc2, m1, m2 = gc[b1], gc[b2], mp[b1], mp[b2]
    dist = (b2 - b1).astype(np.float64)
    r = {"bin1": b1, "bin2": b2, "count": cnt, "peak1": p1.astype(np.int8), "peak2": p2.astype(np.int8),
         "dist": dist, "effective_length1": el1, "gc1": gc1, "mappability1": m1, "short_count1": s1,
         "effective_length2": el2, "gc2": gc2, "mappability2": m2, "short_count2": s2}
    with np.errstate(divide="ignore"):
        r["logl"] = np.log((el1 + 1.0) * (el2 + 1.0) / (B * B))
        r["loggc"] = np.log(gc1 * gc2)
        r["logm"] = np.log(m1 * m2)
        r["logdist"] = np.log((1.0 + dist) / rng)
        mx = ((s1.max() + 1.0) * (s2.max() + 1.0)) if len(s1) else 1.0
        r["logShortCount"] = np.log((s1 + 1.0) * (s2 + 1.0) / mx)
    return r


# ── model ──────────────────────────────────────────────────────────────────

_TERMS = ("logl", "loggc", "logm", "logdist", "logShortCount")


def _ztp_loglik(y, eta, lgy):
    lam = np.exp(eta)
    return float(np.sum(y * eta - lam - lgy - np.log(-np.expm1(-lam))))


def pospoisson_fit(X, y, *, epsilon: float = 1e-7, maxit: int = 30):
    """Zero-truncated Poisson regression, log link, as VGAM ``vglm(y ~ X, pospoisson())`` fits it.

    VGAM ``vglm.fit`` with ``pospoisson()`` (1.1-2 and 1.1-14): start ``l = 0.95 (mean y + 1/16) + 0.05 y``; Fisher
    scoring with score ``y - l - l / expm1(l)`` and weight ``l^2 (1 + 1/expm1(l)) (1/l - 1/expm1(l))``
    (floored at eps^0.75); stop when ``sqrt(n) |dlogL| / (|logL| + epsilon) <= epsilon`` or after ``maxit``
    iterations; half steps while the log-likelihood falls. ``X`` holds the intercept column.
    Returns (coefficients, eta, iterations)."""
    from scipy.special import gammaln
    y = np.asarray(y, np.float64)
    n, p = X.shape
    if n < p:
        raise ValueError(f"There are {p} parameters but only {n} observations")
    lgy = gammaln(y + 1)
    locn = float(np.sum(y)) / n + 1 / 16
    eta = np.log(0.95 * locn + (1 - 0.95) * y)
    wzeps = np.finfo(float).eps ** 0.75
    sq = np.sqrt(n)

    def working(eta):
        lam = np.exp(eta)
        with np.errstate(over="ignore"):                  # expm1(l) = Inf for l > 709: 1 / Inf = 0, as in R
            t6 = np.expm1(lam)
        deriv = (y / lam - 1 - 1 / t6) * lam
        wz = np.maximum((1 + 1 / t6) * (1 / lam - 1 / t6) * lam ** 2, wzeps)
        U = np.sqrt(wz)
        return U, eta + (deriv / U) / U

    def valid(eta):
        lam = np.exp(eta)
        return bool(np.all(np.isfinite(lam)) and np.all(lam > 0))

    new_crit = _ztp_loglik(y, eta, lgy)
    U, z = working(eta)
    old_coef = None
    it = 1
    while True:
        coef = np.linalg.lstsq(X * U[:, None], z * U, rcond=None)[0]
        eta = X @ coef
        old_crit, new_crit = new_crit, _ztp_loglik(y, eta, lgy)
        half = old_coef is not None and (not np.isfinite(new_crit) or new_crit < old_crit or not valid(eta))
        stop = False
        if half:
            saved, step = coef, 1.0
            while True:
                step /= 2
                if step < 1e-6:
                    stop = True
                    break
                coef = (1 - step) * old_coef + step * saved
                eta = X @ coef
                new_crit = _ztp_loglik(y, eta, lgy)
                if valid(eta) and np.isfinite(new_crit) and new_crit > old_crit:
                    break
        more = (not stop) and it < maxit and sq * abs(old_crit - new_crit) / (abs(old_crit) + epsilon) > epsilon
        if not more:
            break
        it += 1
        U, z = working(eta)
        old_coef = coef
    return coef, eta, it


def _ztp_p(count, E):
    """``ppois(count, E, lower.tail = FALSE) / ppois(0, E, lower.tail = FALSE)`` (R, as MAPS codes it)."""
    from scipy.special import pdtrc
    return pdtrc(count, E) / -np.expm1(-E)


def _regress(r, dataset_length) -> dict:
    """MAPS_regression_and_peak_caller.r ``pospoisson_regression`` on one chromosome x set (in place)."""
    y = r["count"].astype(np.float64)
    n = len(y)
    X = np.column_stack([np.ones(n)] + [r[k] for k in _TERMS])
    info = {"chrom": r["chrom"], "set": r["set"], "n": n}
    try:
        c1, eta, it1 = pospoisson_fit(X, y)
        lam = np.exp(eta)
        E1 = -lam / np.expm1(-lam)
        P1 = _ztp_p(y, E1)
        keep = P1 > 1 / n
        c2, _, it2 = pospoisson_fit(X[keep], y[keep])
    except (ValueError, np.linalg.LinAlgError) as e:
        nan = np.full(n, np.nan)
        r.update(expected1=nan, P1=nan, expected=nan, ratio=nan, P=nan, p_bonferroni=nan, fdr=nan)
        info["error"] = str(e)
        return info
    co = np.round(c2, 10)
    lin = co[0] + co[1] * r["logl"] + co[2] * r["loggc"] + co[3] * r["logm"] + co[4] * r["logdist"] \
        + co[5] * r["logShortCount"]
    E2 = np.round(np.exp(lin), 10)
    with np.errstate(invalid="ignore", divide="ignore"):
        E2 = E2 / (1 - np.exp(-E2))                         # 0 / 0 = NaN when exp(lin) rounds to 0, as in R
        P = _ztp_p(y, E2)
        ratio = y / E2
    fdr = np.full(n, np.nan)
    ok = ~np.isnan(P)                                       # p.adjust drops NA and counts only the rest
    fdr[ok] = _bh(P[ok])
    r.update(expected1=E1, P1=P1, expected=E2, ratio=ratio, P=P, p_bonferroni=P * dataset_length, fdr=fdr)
    info["n_na"] = int((~ok).sum())
    info.update(n_refit=int(keep.sum()), iter1=it1, iter2=it2,
                **{f"pass1_{k}": v for k, v in zip(("Intercept",) + _TERMS, c1)},
                **{f"pass2_{k}": v for k, v in zip(("Intercept",) + _TERMS, c2)})
    return info


# ── calls: clusters, classes, formatting ───────────────────────────────────

def _neighbours(b1, b2, k):
    """Index pairs (i, j) with |b1_i - b1_j| <= k and |b2_i - b2_j| <= k (i != j)."""
    key = b1.astype(np.int64) * (1 << 32) + b2
    o = np.argsort(key, kind="stable")
    sk = key[o]
    I, J = [], []
    for dx in range(-k, k + 1):
        for dy in range(-k, k + 1):
            if dx == 0 and dy == 0:
                continue
            q = (b1 + dx).astype(np.int64) * (1 << 32) + (b2 + dy)
            pos = np.searchsorted(sk, q)
            hit = pos < len(sk)
            hit[hit] = sk[pos[hit]] == q[hit]
            I.append(np.flatnonzero(hit))
            J.append(o[pos[hit]])
    return np.concatenate(I), np.concatenate(J)


def _label(tested: Pairs, sig_idx, chroms, B, gap) -> Pairs:
    """``label_peaks`` + ``classify_peaks`` on the calls (rows of ``tested``), in MAPS's row order."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    S = tested.take(sig_idx)
    k = int(gap // B)
    codes = S.a.codes
    # R: unique(df$chr) of rbind(AND calls, XOR calls)
    isand = S.cols["set"] == "AND"
    order_c = list(dict.fromkeys(np.r_[codes[isand], codes[~isand]].tolist()))
    out_idx, label, nei, nlp, size = [], [], [], [], []
    for c in order_c:
        rows = np.flatnonzero(codes == c)
        rows = np.r_[rows[isand[rows]], rows[~isand[rows]]]          # AND calls, then XOR calls
        b1 = S.a.starts[rows] // B
        b2 = S.b.starts[rows] // B
        p = S.cols["P"][rows].copy()
        if (p == 0).any():
            nz = p[p != 0]
            p[p == 0] = nz.min() if len(nz) else 1111111
        I, J = _neighbours(b1, b2, k)
        cn = 1 + np.bincount(I, minlength=len(rows))
        u = np.flatnonzero(cn == 1)
        v = np.flatnonzero(cn >= 2)
        lab = np.zeros(len(rows), np.int64)
        nl = np.zeros(len(rows))
        cs = np.ones(len(rows), np.int64)
        nl[u] = -np.log10(p[u])
        if len(v):
            ncomp, comp = connected_components(coo_matrix((np.ones(len(I)), (I, J)), shape=(len(rows),) * 2),
                                               directed=False)
            cv = comp[v]
            first = np.full(ncomp, len(rows))
            np.minimum.at(first, cv, v)
            rank = np.empty(ncomp, np.int64)
            used = np.unique(cv)
            rank[used[np.argsort(first[used], kind="stable")]] = np.arange(1, len(used) + 1)
            lab[v] = rank[cv]
            # R sums in long double, in row order
            ov = v[np.argsort(lab[v], kind="stable")]
            lv = lab[ov]
            starts = np.r_[0, np.flatnonzero(np.diff(lv)) + 1]
            sums = np.add.reduceat((-np.log10(p[ov])).astype(np.longdouble), starts).astype(np.float64)
            nl[v] = sums[lab[v] - 1]
            cs[v] = np.bincount(lab[v])[lab[v]]
        sel = np.r_[u, v]
        out_idx.append(rows[sel])
        label.append(lab[sel])
        nei.append(cn[sel])
        nlp.append(nl[sel])
        size.append(cs[sel])
    if out_idx:
        oi = np.concatenate(out_idx)
        label, nei, nlp, size = map(np.concatenate, (label, nei, nlp, size))
    else:
        oi = np.zeros(0, np.int64)
        label = nei = size = np.zeros(0, np.int64)
        nlp = np.zeros(0)
    R = S.take(oi)
    ctype = np.where(label == 0, "Singleton", "SharpPeak").astype(object)
    cl = label != 0
    if cl.any():
        # classify_peaks: knee of the sorted cluster sums
        ccode = R.a.codes[cl]
        _, first = np.unique(np.c_[ccode, label[cl]], axis=0, return_index=True)
        y = np.sort(nlp[cl][np.sort(first)], kind="stable")
        z1 = np.arange(1, len(y) + 1) / len(y)
        z2 = y / y.max()
        u2 = (-1 / np.sqrt(2)) * z1 + (1 / np.sqrt(2)) * z2
        ref = y[np.flatnonzero(u2 == u2.min())[0]]
        ctype[cl] = np.where(nlp[cl] < ref, "SharpPeak", "BroadPeak")
    keep = {k: R.cols[k] for k in ("set", "peak1", "peak2", "count", "dist", "expected", "ratio", "P", "fdr")}
    return Pairs(R.a, R.b, {**keep, "CountNei": nei, "label": label, "NegLog10P": nlp, "ClusterSize": size,
                            "ClusterType": ctype})


def _format(S: Pairs, singleton_fdr: float) -> Pairs:
    """MAPS_peak_formatting.r: drop weak singletons, summits, singleton labels."""
    fdr = _q15(S.cols["fdr"])                                   # read back from the .peaks text file
    single = S.cols["label"] == 0
    keep = np.flatnonzero(~single | (fdr < singleton_fdr))
    P = S.take(keep)
    fdr, single = fdr[keep], single[keep]
    chrom = np.asarray(P.a.chroms, dtype=object)
    lab = np.array([f"{c}_{l}" for c, l in zip(chrom, P.cols["label"])], dtype=object)
    summit = np.zeros(len(P))
    if len(P):
        import pandas as pd
        mn = pd.Series(fdr).groupby(lab, sort=False).transform("min").to_numpy()
        summit = (fdr == mn).astype(np.float64)
    summit[single] = 1
    lab[single] = [f"{s}{i}" for i, s in enumerate(lab[single], start=1)]
    return Pairs(P.a, P.b, {"count": P.cols["count"], "expected": P.cols["expected"], "fdr": P.cols["fdr"],
                            "P": P.cols["P"], "ratio": P.cols["ratio"], "set": P.cols["set"], "ClusterLabel": lab,
                            "ClusterSize": P.cols["ClusterSize"], "ClusterType": P.cols["ClusterType"],
                            "ClusterNegLog10P": P.cols["NegLog10P"], "ClusterSummit": summit})


def _r_num(v) -> str:
    """A number as R ``write.table`` prints it with ``options(scipen = 999)`` (15 significant digits)."""
    v = float(v)
    if np.isnan(v):
        return "NA"
    if np.isinf(v):
        return "Inf" if v > 0 else "-Inf"
    s = f"{v:.15g}"
    if "e" in s:
        s = format(Decimal(s), "f")
    return s


def _write_bedpe(L: Pairs, B: int, path) -> None:
    """``{SET}.{FDR}.sig3Dinteractions.bedpe`` as MAPS_peak_formatting.r writes it."""
    head = ["chr1", "start1", "end1", "chr2", "start2", "end2", "count", "expected", "fdr", "ClusterLabel",
            "ClusterSize", "ClusterType", "ClusterNegLog10P", "ClusterSummit"]
    c = L.cols
    chrom = np.asarray(L.a.chroms, dtype=object)
    with open(path, "w") as f:
        f.write("\t".join(head) + "\n")
        for i in range(len(L)):
            s1, s2 = int(L.a.starts[i]), int(L.b.starts[i])
            f.write(f"{chrom[i]}\t{s1}\t{s1 + B}\t{chrom[i]}\t{s2}\t{s2 + B}\t{int(c['count'][i])}\t"
                    f"{_r_num(c['expected'][i])}\t{_r_num(c['fdr'][i])}\t{c['ClusterLabel'][i]}\t"
                    f"{int(c['ClusterSize'][i])}\t{c['ClusterType'][i]}\t{_r_num(c['ClusterNegLog10P'][i])}\t"
                    f"{_r_num(c['ClusterSummit'][i])}\n")
