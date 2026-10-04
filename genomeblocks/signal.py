from __future__ import annotations
import os
import shutil
from os import cpu_count
from typing import List, Tuple, Sequence

import numpy as np
from tqdm import tqdm

from .loci import Loci


# ── backend selection ────────────────────────────────────────────────────────
# Prefer pybigtools (Rust, fastest) → fall back to our pure-Python reader.
# Both are thread-safe, so the threading architecture works with either.

def _open_pybigtools(path):
    """Adapter wrapping pybigtools handle to match our API."""
    global _FILL_KW
    import pybigtools
    h = pybigtools.open(path, "r")
    if _FILL_KW is None:
        _FILL_KW = _fill_kwarg(h)
    return _PyBigToolsHandle(h)


# pybigtools' native summary kinds — anything else falls back to a values() reduction.
_NATIVE_SUMMARY = {'mean', 'min', 'max'}


def _fill_kwarg(h) -> str:
    """Name of pybigtools' fill-value argument for ``values()``.

    pybigtools 0.3 renamed ``missing=`` to ``fillna=`` and warns on every
    ``missing=`` call; a later release drops it. Older releases only know
    ``missing=``.
    """
    try:
        from importlib.metadata import version
        major, minor = (int(x) for x in version('pybigtools').split('.')[:2])
        return 'fillna' if (major, minor) >= (0, 3) else 'missing'
    except Exception:          # no metadata / odd version string: probe the handle
        import inspect
        try:
            params = inspect.signature(h.values).parameters
        except (TypeError, ValueError):
            return 'missing'
        return 'fillna' if 'fillna' in params else 'missing'


_FILL_KW: str | None = None    # resolved on the first pybigtools open


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
                                  exact=exact, **{_FILL_KW: missing})
        # sum / std / coverage: reduce from per-base values in numpy.
        vals = self._h.values(chrom, start, end, **{_FILL_KW: np.nan})
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
        return self._h.values(chrom, start, end, **{_FILL_KW: 0.0}).astype(np.float64, copy=False)

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


def _available_ram_bytes() -> int:
    """Best-effort available physical RAM in bytes.

    Uses POSIX ``sysconf`` (Linux/macOS); if the keys are unavailable, falls
    back to free space on ``/`` so the guard still returns *something* rather
    than crashing the extraction.
    """
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (ValueError, AttributeError, OSError):
        return shutil.disk_usage("/").free


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


def _resolve_opener(backend: str | None):
    """Return ``(opener, backend_name)`` for a requested backend.

    Pure resolution — does **not** mutate module state, so concurrent
    ``signal()`` calls with different backends don't race on globals.
    ``backend=None`` uses the auto-detected default.
    """
    if backend is None:
        return _bw_open, _bw_backend
    if backend == 'pybigtools':
        return _open_pybigtools, 'pybigtools'
    if backend == 'bigwig':
        from . import bigwig as _bw_mod
        return _bw_mod.open, 'bigwig'
    raise ValueError(f"Unknown backend: {backend!r}")


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
    workers: int | None = 1,
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
            many processes, capped at ``cpu_count()`` and at the available
            work. ``None`` = half the cores, the polite choice on a shared
            machine.
        span: Use full locus span instead of center±flank
        backend: 'pybigtools' | 'bigwig' | None (auto-detect)
        exact: pybigtools base-accurate binning when True (default); zoom
            interpolation when False (~3× faster, approximate). Ignored by
            the pure-python backend.

    Returns:
        ndarray of shape ``(n_loci, n_tracks, n_bins)``.
    """
    opener, backend_name = _resolve_opener(backend)

    n_loci, n_tracks = len(loci), len(bigwigs)
    if n_loci == 0:
        raise ValueError("No loci provided.")

    bytes_need = n_loci * n_tracks * n_bins * np.dtype(dtype).itemsize
    if bytes_need > 0.5 * _available_ram_bytes():
        raise MemoryError(
            f"Signal cube needs ~{bytes_need / 1e9:.1f} GB, over half of "
            f"available RAM; extract in loci chunks and stream to disk.")

    # flatten loci into pickle-friendly columns (needed for MP, cheap for seq)
    chrom_list = [loc.chrom for loc in loci]
    starts = np.fromiter((loc.start for loc in loci),
                         dtype=np.int64, count=n_loci)
    ends = np.fromiter((loc.end for loc in loci),
                       dtype=np.int64, count=n_loci)

    # workers=None defaults to half the cores: each worker spins up
    # pybigtools' own tokio pool, so the other half stays free for the OS /
    # other users. An explicit count is honoured up to the core count.
    if workers is None:
        workers = max(1, (cpu_count() or 2) // 2)
    workers = max(1, min(workers, n_tracks * max(1, (n_loci + 999) // 1000),
                         cpu_count() or 1))

    if verbose:
        print(f"[INFO] Extracting {n_tracks} bigwigs for {n_loci} loci "
              f"into {n_bins} bins (span={span}, agg='{agg}', "
              f"backend='{backend_name}', exact={exact}, workers={workers}).")

    cube = np.zeros((n_loci, n_tracks, n_bins), dtype=dtype)

    if workers <= 1:
        _run_sequential(cube, bigwigs, chrom_list, starts, ends,
                        n_bins, flank, agg, span, exact,
                        opener=opener, progress=progress)
    else:
        _run_multiprocess(cube, list(bigwigs), chrom_list, starts, ends,
                          n_bins, flank, agg, span, exact,
                          backend_name, workers, progress=progress)

    return cube


def _run_sequential(cube, bigwigs, chrom_list, starts, ends,
                    n_bins, flank, agg, span, exact, *, opener, progress):
    if not bigwigs:
        return
    hs = [opener(p) for p in bigwigs]
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


def _tmm_norm_factors(data, trim_lfc=0.3, trim_mag=0.05, index_ref=None):
    """edgeR Trimmed-Mean-of-M-values normalization factors, one per column.

    A self-contained implementation of the standard TMM algorithm (Robinson &
    Oshlack, *Genome Biology* 2010) so genomeblocks carries no extra dependency
    for :func:`tmm`. Rows are features, columns are samples; returns a factor
    per sample, scaled to a geometric mean of 1.
    """
    x = np.asarray(data, dtype=float).T                     # (samples, features)
    lib_size = x.sum(axis=1)
    mask = x == 0
    if index_ref is None:
        xr = x.copy()
        xr[:, np.all(mask, axis=0)] = np.nan                # drop all-zero features
        p75 = np.nanpercentile(xr, 75, axis=1)
        index_ref = int(np.argmin(np.abs(p75 - p75.mean())))
    mask = mask.copy()
    mask[:, mask[index_ref]] = True                          # mask where the ref is 0
    x = x.copy(); x[mask] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        norm_x = x / lib_size[:, None]
        logs = np.log2(norm_x)
        m_g = logs - logs[index_ref]                         # log fold-change vs ref
        a_g = (logs + logs[index_ref]) / 2                   # average abundance
        pm = np.nanquantile(m_g, [trim_lfc, 1 - trim_lfc], axis=1, method="nearest")[..., None]
        pa = np.nanquantile(a_g, [trim_mag, 1 - trim_mag], axis=1, method="nearest")[..., None]
        mask = mask | (m_g < pm[0]) | (m_g > pm[1])
        mask = mask | (a_g < pa[0]) | (a_g > pa[1])
        w = (1 - norm_x) / x                                 # asymptotic variance
        w = 1 / (w + w[index_ref])
    w[mask] = 0
    m_g[mask] = 0
    w /= w.sum(axis=1)[:, None]
    f = np.sum(w * m_g, axis=1)
    f -= f.mean()                                            # geometric mean -> 1
    return 2 ** f


def tmm(cube: np.ndarray) -> np.ndarray:
    """TMM-normalize a signal cube (regions × tracks × bins).

    Computes per-track TMM normalization factors and library-size
    scaling so that tracks become comparable.

    Args:
        cube: Signal array of shape (regions, tracks, bins)

    Returns:
        Normalized copy of the cube (same shape).
    """
    means = np.nanmean(cube, axis=2)                        # (regions, tracks)
    factors = _tmm_norm_factors(means)                       # (tracks,)
    lib_size = means.sum(0)                                  # (tracks,)
    scale = 1.0 / (factors * lib_size / 1_000_000)
    return cube * scale[None, :, None]


Loci.signal = signal
