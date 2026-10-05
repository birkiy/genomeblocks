"""bigWig signal over loci, as a cube of (rows x tracks x bins).

    S = cre.signal(["atac.bw", "h3k27ac.bw"], n_bins=200, flank=3000)

Row ``i`` of the cube is row ``i`` of the loci (rows on chromosomes a
bigWig lacks stay 0), so it lines up with every other table built from
them. Files are read through the bigwig backend: pybigtools by default,
pyBigWig or the pure-Python reader on request (``backend=``); open
pyBigWig / pybigtools handles are accepted as tracks too.

Convert a cube with :func:`genomeblocks.interop.cube_to_xarray`,
``cube_to_anndata`` or ``cube_to_pandas``.
"""
from __future__ import annotations

import os
import shutil
from os import cpu_count
from typing import List, Optional, Sequence, Tuple

import numpy as np


def _available_ram_bytes() -> int:
    """Best-effort available physical RAM in bytes."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (ValueError, AttributeError, OSError):
        return shutil.disk_usage("/").free


def _even_ranges(total: int, n_chunks: int) -> List[Tuple[int, int]]:
    """Split [0, total) into at most ``n_chunks`` non-empty contiguous ranges."""
    n_chunks = max(1, min(n_chunks, total))
    base, rem = divmod(total, n_chunks)
    ranges, off = [], 0
    for i in range(n_chunks):
        step = base + (1 if i < rem else 0)
        if step == 0:
            break
        ranges.append((off, off + step))
        off += step
    return ranges


def _windows(starts, ends, sizes, n_bins, flank, span):
    """Per row: (left, right, first bin, number of bins) of the window to read
    and whether the chromosome is in the file. Windows are not clipped: the
    bigwig handles bin a window that leaves the chromosome over the full grid
    (bases outside have no data), so edge bins keep the signal they have and
    every backend agrees."""
    if span:
        L, R = starts.astype(np.int64), ends.astype(np.int64)
    else:
        c = (starts + ends) // 2
        L, R = c - flank, c + flank
    pre = np.zeros(len(L), np.int64)
    core = np.full(len(L), max(1, n_bins), np.int64)
    ok = (R > L) & (sizes > 0)
    return L, R, pre, core, ok


def _extract(cube, handles, chroms, starts, ends, t_lo, l_lo, l_hi, n_bins, flank, agg, span, exact):
    """Fill cube[l_lo:l_hi, t_lo:t_lo + len(handles), :]."""
    for off, h in enumerate(handles):
        sizes_d = h.chroms()
        rows = np.arange(l_lo, l_hi)
        sizes = np.array([sizes_d.get(c, 0) for c in chroms[l_lo:l_hi]], np.int64)
        L, R, pre, core, ok = _windows(starts[l_lo:l_hi], ends[l_lo:l_hi], sizes, n_bins, flank, span)
        t = t_lo + off
        stats = h.stats_array
        for r, c, a, b, p, k in zip(rows[ok].tolist(), chroms[l_lo:l_hi][ok].tolist(), L[ok].tolist(),
                                    R[ok].tolist(), pre[ok].tolist(), core[ok].tolist()):
            cube[r, t, p:p + k] = stats(c, a, b, n_bins=k, stat=agg, exact=exact, missing=0.0)


def _mp_worker(paths, shm_name, shape, dtype_str, chroms, starts_b, ends_b, t_lo, t_hi, l_lo, l_hi,
               n_bins, flank, agg, span, exact, backend):
    """Process-pool worker: attach the shared cube, extract, detach."""
    from multiprocessing import shared_memory
    from .backends.bigwig import open_bigwig
    starts = np.frombuffer(starts_b, dtype=np.int64)
    ends = np.frombuffer(ends_b, dtype=np.int64)
    shm = shared_memory.SharedMemory(name=shm_name)
    try:
        cube = np.ndarray(shape, dtype=np.dtype(dtype_str), buffer=shm.buf)
        hs = [open_bigwig(p, backend=backend) for p in paths[t_lo:t_hi]]
        try:
            _extract(cube, hs, np.asarray(chroms, dtype=object), starts, ends, t_lo, l_lo, l_hi,
                     n_bins, flank, agg, span, exact)
        finally:
            for h in hs:
                h.close()
    finally:
        shm.close()


def _tracks(bigwigs) -> list:
    if isinstance(bigwigs, dict):
        return list(bigwigs.values())
    if isinstance(bigwigs, (str, os.PathLike)) or not hasattr(bigwigs, "__iter__"):
        return [bigwigs]
    return list(bigwigs)


def signal(loci, bigwigs, *, n_bins: int = 200, flank: int = 3_000, agg: str = "mean",
           dtype=np.float32, progress: bool = True, workers: Optional[int] = 1, span: bool = False,
           verbose: bool = True, backend: Optional[str] = None, exact: bool = True) -> np.ndarray:
    """bigWig signal of every locus: a ``(rows, tracks, bins)`` array.

    Args:
        loci: anything :func:`~genomeblocks.as_loci` takes.
        bigwigs: one bigWig or several (paths, open pyBigWig / pybigtools
            handles, or a ``{name: path}`` dict — tracks in its order).
        n_bins: bins per locus.
        flank: bp each side of the locus centre (ignored with ``span=True``).
        agg: 'mean' | 'min' | 'max' | 'sum' | 'std' | 'coverage'.
        span: bin the whole locus instead of centre ± flank.
        workers: 1 (default) reads sequentially — one native call per region
            and track; >1 uses that many processes writing into a shared-memory
            cube (needs paths); ``None`` = half the cores.
        backend: 'pybigtools' (default), 'pybigwig' or 'python'.
        exact: base-accurate binning (default); False lets pybigtools /
            pyBigWig use zoom levels (~3x faster, approximate).
    """
    from .backends import resolve
    from .backends.bigwig import open_bigwig
    from .interop import as_loci
    L = as_loci(loci)
    tracks = _tracks(bigwigs)
    n_loci, n_tracks = len(L), len(tracks)
    if n_loci == 0:
        raise ValueError("No loci provided.")
    if n_tracks == 0:
        raise ValueError("No bigWig provided.")
    need = n_loci * n_tracks * n_bins * np.dtype(dtype).itemsize
    if need > 0.5 * _available_ram_bytes():
        raise MemoryError(f"The signal cube needs ~{need / 1e9:.1f} GB, over half of the available RAM; "
                          f"extract in chunks of loci (L[a:b].signal(...)).")
    paths = all(isinstance(t, (str, os.PathLike)) for t in tracks)
    if workers is None:
        workers = max(1, (cpu_count() or 2) // 2)
    workers = max(1, min(workers, n_tracks * max(1, (n_loci + 999) // 1000), cpu_count() or 1))
    if workers > 1 and not paths:
        raise ValueError("workers > 1 needs bigWig paths (open handles cannot be shared between "
                         "processes); pass the file paths, or workers=1.")
    name = resolve("bigwig", backend) if paths or backend else "handle"
    if verbose:
        print(f"[INFO] Extracting {n_tracks} bigwigs for {n_loci} loci into {n_bins} bins (span={span}, "
              f"agg='{agg}', backend='{name}', exact={exact}, workers={workers}).")
    chroms = L.chroms
    starts, ends = L.starts, L.ends
    cube = np.zeros((n_loci, n_tracks, n_bins), dtype=dtype)
    if workers <= 1:
        hs = [open_bigwig(t, backend=backend) for t in tracks]
        try:
            from tqdm import tqdm
            step = 512
            bar = tqdm(total=n_loci * n_tracks, dynamic_ncols=True, disable=not progress)
            for lo in range(0, n_loci, step):
                hi = min(n_loci, lo + step)
                _extract(cube, hs, chroms, starts, ends, 0, lo, hi, n_bins, flank, agg, span, exact)
                bar.update((hi - lo) * n_tracks)
            bar.close()
        finally:
            for h in hs:
                h.close()
    else:
        _run_multiprocess(cube, [str(t) for t in tracks], chroms, starts, ends, n_bins, flank, agg, span,
                          exact, resolve("bigwig", backend), workers, progress)
    return cube


def _run_multiprocess(cube, paths, chroms, starts, ends, n_bins, flank, agg, span, exact, backend,
                      workers, progress):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from multiprocessing import shared_memory
    from tqdm import tqdm
    n_loci, n_tracks = len(chroms), len(paths)
    t_splits = min(workers, n_tracks)
    l_per_t = max(1, workers // t_splits)
    shm = shared_memory.SharedMemory(create=True, size=cube.nbytes)
    try:
        shared = np.ndarray(cube.shape, dtype=cube.dtype, buffer=shm.buf)
        shared[:] = 0
        sb, eb = starts.tobytes(), ends.tobytes()
        chrom_list = chroms.tolist()
        tasks = [(paths, shm.name, cube.shape, cube.dtype.name, chrom_list, sb, eb, t_lo, t_hi, l_lo, l_hi,
                  n_bins, flank, agg, span, exact, backend)
                 for t_lo, t_hi in _even_ranges(n_tracks, t_splits)
                 for l_lo, l_hi in _even_ranges(n_loci, l_per_t)]
        bar = tqdm(total=len(tasks), dynamic_ncols=True, desc="chunks", disable=not progress)
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for f in as_completed([ex.submit(_mp_worker, *t) for t in tasks]):
                f.result()
                bar.update()
        bar.close()
        cube[:] = shared
    finally:
        shm.close()
        shm.unlink()


def _tmm_norm_factors(data, trim_lfc=0.3, trim_mag=0.05, index_ref=None):
    """edgeR Trimmed-Mean-of-M-values factors, one per column (Robinson & Oshlack 2010),
    scaled to a geometric mean of 1. Rows are features, columns samples."""
    x = np.asarray(data, dtype=float).T
    lib_size = x.sum(axis=1)
    mask = x == 0
    if index_ref is None:
        xr = x.copy()
        xr[:, np.all(mask, axis=0)] = np.nan
        p75 = np.nanpercentile(xr, 75, axis=1)
        index_ref = int(np.argmin(np.abs(p75 - p75.mean())))
    mask = mask.copy()
    mask[:, mask[index_ref]] = True
    x = x.copy()
    x[mask] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        norm_x = x / lib_size[:, None]
        logs = np.log2(norm_x)
        m_g = logs - logs[index_ref]
        a_g = (logs + logs[index_ref]) / 2
        pm = np.nanquantile(m_g, [trim_lfc, 1 - trim_lfc], axis=1, method="nearest")[..., None]
        pa = np.nanquantile(a_g, [trim_mag, 1 - trim_mag], axis=1, method="nearest")[..., None]
        mask = mask | (m_g < pm[0]) | (m_g > pm[1])
        mask = mask | (a_g < pa[0]) | (a_g > pa[1])
        w = (1 - norm_x) / x
        w = 1 / (w + w[index_ref])
    w[mask] = 0
    m_g[mask] = 0
    w /= w.sum(axis=1)[:, None]
    f = np.sum(w * m_g, axis=1)
    f -= f.mean()
    return 2 ** f


def tmm(cube: np.ndarray) -> np.ndarray:
    """TMM-normalise a (regions x tracks x bins) cube: per-track TMM factors and
    library-size scaling (counts per million), so tracks become comparable."""
    means = np.nanmean(cube, axis=2)
    lib_size = np.nansum(means, axis=0)
    has = lib_size > 0
    if not has.any():
        raise ValueError("tmm: every track is empty (no signal in any region)")
    if not has.all():
        import warnings
        warnings.warn(f"tmm: track(s) {np.flatnonzero(~has).tolist()} have no signal and are left at 0",
                      RuntimeWarning, stacklevel=2)
    scale = np.zeros(cube.shape[1])
    factors = _tmm_norm_factors(means[:, has])
    scale[has] = 1.0 / (factors * lib_size[has] / 1_000_000)
    return cube * scale[None, :, None]
