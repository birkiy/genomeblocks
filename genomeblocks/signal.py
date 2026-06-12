from __future__ import annotations
import shutil
from os import cpu_count
from typing import List, Tuple, Sequence, Dict

import numpy as np
from tqdm import tqdm

from .loci import Loci


# ── backend selection ────────────────────────────────────────────────────────
# Prefer pybigtools (Rust, fastest) → fall back to our pure-Python reader.
# Both are thread-safe, so the threading architecture works with either.

def _open_pybigtools(path):
    """Adapter wrapping pybigtools handle to match our API."""
    import pybigtools
    return _PyBigToolsHandle(pybigtools.open(path, "r"))


# pybigtools' native summary kinds — anything else falls back to a values() reduction.
_NATIVE_SUMMARY = {'mean', 'min', 'max'}


class _PyBigToolsHandle:
    """Thin wrapper so pybigtools matches BigWigReader's interface.

    The hot path is ``stats_array`` → one FFI call into Rust that does I/O,
    decompression, and binning without releasing a single Python list.
    """
    __slots__ = ('_h',)
    def __init__(self, h):
        self._h = h

    def chroms(self):
        return dict(self._h.chroms())

    def stats_array(self, chrom, start, end, *, n_bins=1,
                    stat='mean', exact=True, missing=0.0):
        if stat in _NATIVE_SUMMARY:
            return self._h.values(chrom, start, end,
                                  bins=n_bins, summary=stat,
                                  exact=exact, missing=missing)
        # sum / std / coverage: reduce from per-base values in numpy.
        vals = self._h.values(chrom, start, end, missing=np.nan)
        n = vals.size
        if n == 0 or n_bins <= 0:
            return np.full(max(n_bins, 1), missing, dtype=np.float64)
        if n % n_bins:
            pad = n_bins - (n % n_bins)
            vals = np.concatenate([vals, np.full(pad, np.nan)])
        chunks = vals.reshape(n_bins, -1)
        with np.errstate(all='ignore'):
            if stat == 'sum':
                out = np.nansum(chunks, axis=1)
            elif stat == 'std':
                out = np.nanstd(chunks, axis=1)
            elif stat in ('coverage', 'cov'):
                out = (~np.isnan(chunks)).mean(axis=1)
            else:
                raise ValueError(f"Unknown stat: {stat!r}")
        return np.where(np.isnan(out), missing, out).astype(np.float64)

    def stats(self, chrom, start, end, *, n_bins=1, nBins=None,
              stat='mean', type=None, **_):
        """Legacy list-returning API (used by bench_bigwig.py)."""
        if nBins is not None:
            n_bins = nBins
        if type is not None:
            stat = type
        arr = self.stats_array(chrom, start, end, n_bins=n_bins,
                               stat=stat, missing=np.nan)
        return [None if np.isnan(v) else float(v) for v in arr]

    def values(self, chrom, start, end):
        return self._h.values(chrom, start, end, missing=0.0).astype(np.float64, copy=False)

    def close(self):
        self._h.close()


def _detect_backend():
    """Return (opener_func, backend_name).

    Prefer pybigtools (Rust) — releases GIL, scales with threads,
    fastest for base-pair resolution.  Falls back to pure-python
    reader (zero compiled deps, exact pyBigWig match for binned stats).
    """
    try:
        import pybigtools  # noqa: F401
        return _open_pybigtools, 'pybigtools'
    except ImportError:
        from . import bigwig
        return bigwig.open, 'bigwig'


_bw_open, _bw_backend = _detect_backend()


def _even_ranges(total: int, n_chunks: int) -> List[Tuple[int, int]]:
    """Split [0, total) into *n_chunks* non-empty contiguous ranges.

    Unlike naive ceil-div slicing, this distributes the remainder across
    the early chunks so no chunk is empty even when total % n_chunks != 0.
    """
    n_chunks = max(1, min(n_chunks, total))
    base, rem = divmod(total, n_chunks)
    ranges = []
    off = 0
    for i in range(n_chunks):
        step = base + (1 if i < rem else 0)
        if step == 0:
            break
        ranges.append((off, off + step))
        off += step
    return ranges


def plan_workers(n_tracks: int, n_loci: int, *, cores: int | None = None,
                max_bw_parallel: int = 6) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """
    Plan the distribution of work across multiple threads.

    Args:
        n_tracks: Number of bigwig tracks
        n_loci: Number of genomic loci
        cores: Number of CPU cores to use (default: all available)
        max_bw_parallel: Maximum number of bigwig files to process in parallel

    Returns:
        List of work chunks as [((track_start,track_end), (loci_start,loci_end)), ...]
    """
    c = cores or cpu_count() or 1
    t_chunks = min(n_tracks, max_bw_parallel, c)
    l_per_t = max(1, c // max(1, t_chunks))
    t_ranges = _even_ranges(n_tracks, t_chunks)
    l_ranges = _even_ranges(n_loci, l_per_t)
    return [(tr, lr) for tr in t_ranges for lr in l_ranges]


def _extract_chunk(cube: np.ndarray, hs: list, ch: dict,
                   chroms: Sequence[str],
                   starts: np.ndarray, ends: np.ndarray,
                   t_lo: int, l_lo: int, l_hi: int,
                   n_bins: int, flank: int, agg: str,
                   span: bool, exact: bool):
    """Write stats into ``cube[l_lo:l_hi, t_lo:t_lo+len(hs), :]``.

    Called on the main thread for the sequential path and inside a child
    process (attached to SharedMemory) for the multiprocessing path.
    """
    native = all(hasattr(h, 'stats_array') for h in hs)
    for r in range(l_lo, l_hi):
        chrom = chroms[r]
        if chrom not in ch:
            continue
        cs = int(starts[r]); ce = int(ends[r])
        chrom_size = ch[chrom]
        if span:
            L, R = max(0, cs), min(chrom_size, ce)
            pre, core = 0, max(1, n_bins)
        else:
            c = (cs + ce) // 2
            L, R = c - flank, c + flank
            pre = max(0, -L)
            post = max(0, R - chrom_size)
            L = max(0, L); R = min(chrom_size, R)
            core = n_bins - (pre + post)
        if core <= 0 or R <= L:
            continue
        for off, h in enumerate(hs):
            t = t_lo + off
            if native:
                xs = h.stats_array(chrom, L, R, n_bins=core,
                                   stat=agg, exact=exact, missing=0.0)
            else:
                raw = h.stats(chrom, L, R, n_bins=core, stat=agg)
                xs = np.fromiter(
                    (0.0 if x is None else x for x in raw),
                    dtype=np.float64, count=len(raw),
                )
            cube[r, t, pre:pre + core] = xs


def _mp_worker(bw_paths, shm_name, shape, dtype_str,
               chrom_list, starts_bytes, ends_bytes,
               t_lo, t_hi, l_lo, l_hi,
               n_bins, flank, agg, span, exact, backend_name):
    """Process-pool worker: attach SharedMemory, extract, detach."""
    from multiprocessing import shared_memory
    import numpy as _np

    if backend_name == 'pybigtools':
        opener = _open_pybigtools
    else:
        from . import bigwig as _bw_mod
        opener = _bw_mod.open

    starts = _np.frombuffer(starts_bytes, dtype=_np.int64)
    ends = _np.frombuffer(ends_bytes, dtype=_np.int64)
    shm = shared_memory.SharedMemory(name=shm_name)
    try:
        cube = _np.ndarray(shape, dtype=_np.dtype(dtype_str), buffer=shm.buf)
        hs = [opener(p) for p in bw_paths[t_lo:t_hi]]
        try:
            ch = hs[0].chroms()
            _extract_chunk(cube, hs, ch, chrom_list, starts, ends,
                           t_lo, l_lo, l_hi, n_bins, flank, agg, span, exact)
        finally:
            for h in hs:
                h.close()
    finally:
        shm.close()


def signal(
    loci: Loci,
    bigwigs: Sequence[str],
    *,
    n_bins: int = 200,
    flank: int = 3_000,
    agg: str = "mean",
    dtype: str | np.dtype = np.float32,
    progress: bool = True,
    workers: int = 1,
    span: bool = False,
    verbose: bool = True,
    backend: str | None = None,
    exact: bool = True,
    ) -> np.ndarray:
    """Extract signal from bigwig files for given genomic loci.

    Default path is **sequential** — a single native-Rust pass with
    pybigtools reaches ~50k region-tracks/s, which is fastest for
    typical heatmap / browser / per-locus-profile workloads.

    For scale (many bigwigs × many loci) pass ``workers > 1``: extraction
    runs in a ``ProcessPoolExecutor`` writing into a shared-memory cube.
    Multiprocessing (not threading) is used because pybigtools serialises
    concurrent Python threads; processes give linear scale up to CPU count.

    Args:
        loci: Genomic loci to extract signal from
        bigwigs: BigWig file paths
        n_bins: Bins per region
        flank: bp each side of locus center (ignored when ``span=True``)
        agg: 'mean' | 'min' | 'max' | 'sum' | 'std' | 'coverage'
        dtype: Output cube dtype (float32 default)
        progress: Show a tqdm bar
        workers: 1 (default) = sequential. >1 = multiprocessing with that
            many processes. Capped at ``min(workers, n_tracks, cpu_count())``.
        span: Use full locus span instead of center±flank
        backend: 'pybigtools' | 'bigwig' | None (auto-detect)
        exact: pybigtools base-accurate binning when True (default); zoom
            interpolation when False (~3× faster, approximate). Ignored by
            the pure-python backend.

    Returns:
        ndarray of shape ``(n_loci, n_tracks, n_bins)``.
    """
    global _bw_open, _bw_backend
    if backend is not None:
        if backend == 'pybigtools':
            _bw_open = _open_pybigtools
            backend_name = 'pybigtools'
        elif backend == 'bigwig':
            from . import bigwig as _bw_mod
            _bw_open = _bw_mod.open
            backend_name = 'bigwig'
        else:
            raise ValueError(f"Unknown backend: {backend!r}")
    else:
        backend_name = _bw_backend

    n_loci, n_tracks = len(loci), len(bigwigs)
    if n_loci == 0:
        raise ValueError("No loci provided.")

    bytes_need = n_loci * n_tracks * n_bins * np.dtype(dtype).itemsize
    if bytes_need > 0.5 * shutil.disk_usage("/").free:
        raise MemoryError("Cube may exceed safe RAM; try disk-chunk mode.")

    # flatten loci into pickle-friendly columns (needed for MP, cheap for seq)
    chrom_list = [loc.chrom for loc in loci]
    starts = np.fromiter((loc.start for loc in loci),
                         dtype=np.int64, count=n_loci)
    ends = np.fromiter((loc.end for loc in loci),
                       dtype=np.int64, count=n_loci)

    # Use at most half the cores: each worker spins up pybigtools' own tokio
    # pool, so we leave the other half free for the OS / user processes.
    _core_cap = max(1, (cpu_count() or 2) // 2)
    workers = max(1, min(workers, n_tracks * max(1, (n_loci + 999) // 1000),
                         _core_cap))

    if verbose:
        print(f"[INFO] Extracting {n_tracks} bigwigs for {n_loci} loci "
              f"into {n_bins} bins (span={span}, agg='{agg}', "
              f"backend='{backend_name}', exact={exact}, workers={workers}).")

    cube = np.zeros((n_loci, n_tracks, n_bins), dtype=dtype)

    if workers <= 1:
        _run_sequential(cube, bigwigs, chrom_list, starts, ends,
                        n_bins, flank, agg, span, exact,
                        progress=progress)
    else:
        _run_multiprocess(cube, list(bigwigs), chrom_list, starts, ends,
                          n_bins, flank, agg, span, exact,
                          backend_name, workers, progress=progress)

    return cube


def _run_sequential(cube, bigwigs, chrom_list, starts, ends,
                    n_bins, flank, agg, span, exact, *, progress):
    if not bigwigs:
        return
    hs = [_bw_open(p) for p in bigwigs]
    try:
        ch = hs[0].chroms()
        n_loci = len(chrom_list)
        n_tracks = len(hs)
        STEP = 256  # progress granularity — keeps tqdm overhead negligible
        bar = tqdm(total=n_loci * n_tracks, dynamic_ncols=True,
                   disable=not progress)
        try:
            for l_lo in range(0, n_loci, STEP):
                l_hi = min(n_loci, l_lo + STEP)
                _extract_chunk(cube, hs, ch, chrom_list, starts, ends,
                               0, l_lo, l_hi, n_bins, flank, agg,
                               span, exact)
                bar.update((l_hi - l_lo) * n_tracks)
        finally:
            bar.close()
    finally:
        for h in hs:
            h.close()


def _run_multiprocess(cube, bigwigs, chrom_list, starts, ends,
                      n_bins, flank, agg, span, exact,
                      backend_name, workers, *, progress):
    from multiprocessing import shared_memory
    from concurrent.futures import ProcessPoolExecutor, as_completed

    n_loci = len(chrom_list)
    n_tracks = len(bigwigs)

    # Prefer splitting by track (each process owns its own bigwig handles).
    # If workers > n_tracks, sub-split the loci axis.
    t_splits = min(workers, n_tracks)
    l_per_t = max(1, workers // t_splits)
    t_ranges = _even_ranges(n_tracks, t_splits)
    l_ranges = _even_ranges(n_loci, l_per_t)

    shm = shared_memory.SharedMemory(create=True, size=cube.nbytes)
    try:
        shm_cube = np.ndarray(cube.shape, dtype=cube.dtype, buffer=shm.buf)
        shm_cube[:] = 0

        # bytes views of starts/ends avoid re-pickling numpy in each task
        starts_bytes = starts.tobytes()
        ends_bytes = ends.tobytes()

        tasks = []
        for t_lo, t_hi in t_ranges:
            for l_lo, l_hi in l_ranges:
                tasks.append((bigwigs, shm.name, cube.shape, cube.dtype.name,
                              chrom_list, starts_bytes, ends_bytes,
                              t_lo, t_hi, l_lo, l_hi,
                              n_bins, flank, agg, span, exact, backend_name))

        bar = tqdm(total=len(tasks), dynamic_ncols=True,
                   desc='chunks', disable=not progress)
        try:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futures = [ex.submit(_mp_worker, *t) for t in tasks]
                for f in as_completed(futures):
                    f.result()
                    bar.update()
        finally:
            bar.close()

        cube[:] = shm_cube  # copy out before unlink
    finally:
        shm.close()
        shm.unlink()


def tmm(cube: np.ndarray) -> np.ndarray:
    """TMM-normalize a signal cube (regions × tracks × bins).

    Computes per-track TMM normalization factors and library-size
    scaling so that tracks become comparable.

    Args:
        cube: Signal array of shape (regions, tracks, bins)

    Returns:
        Normalized copy of the cube (same shape).
    """
    import conorm
    means = np.nanmean(cube, axis=2)                        # (regions, tracks)
    factors = conorm.tmm_norm_factors(means)                 # (tracks,)
    lib_size = means.sum(0)                                  # (tracks,)
    scale = 1.0 / (factors * lib_size / 1_000_000)
    return cube * scale[None, :, None]


Loci.signal = signal
