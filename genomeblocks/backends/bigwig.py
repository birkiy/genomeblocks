"""bigWig backends: one handle interface over pybigtools, pyBigWig and a
pure-Python reader.

    h = open_bigwig("signal.bw", backend=None)     # path or an open handle
    h.chroms()                                     # {chrom: size}
    h.stats_array(chrom, start, end, n_bins=100, stat="mean", exact=True, missing=0.0)
    h.values(chrom, start, end)                    # per base, NaN where no data
    h.close()

``stat`` is one of mean, min, max, sum, std, coverage. pybigtools (Rust) is
the default; pyBigWig (C, libBigWig) and the pure-Python reader (numpy +
mmap, no compiled dependency) give the same numbers. An already-open
``pyBigWig`` or ``pybigtools`` handle is accepted too, so signal extraction
works on files a caller has opened elsewhere.
"""
from __future__ import annotations

import numpy as np

from . import resolve

_NATIVE = {"mean", "min", "max"}
_FILL_KW = None


def _fill_kwarg(h) -> str:
    """pybigtools renamed ``missing=`` to ``fillna=`` in 0.3 (and warns on the old one)."""
    try:
        from importlib.metadata import version
        major, minor = (int(x) for x in version("pybigtools").split(".")[:2])
        return "fillna" if (major, minor) >= (0, 3) else "missing"
    except Exception:
        import inspect
        try:
            return "fillna" if "fillna" in inspect.signature(h.values).parameters else "missing"
        except (TypeError, ValueError):
            return "missing"


def _reduce(vals, n_bins, stat, missing):
    """Bin per-base values (NaN = no data) for stats the engines lack natively."""
    n = vals.size
    if n == 0 or n_bins <= 0:
        return np.full(max(n_bins, 1), missing, dtype=np.float64)
    if n % n_bins:
        vals = np.concatenate([vals, np.full(n_bins - n % n_bins, np.nan)])
    chunks = vals.reshape(n_bins, -1)
    with np.errstate(all="ignore"):
        if stat == "sum":
            out = np.nansum(chunks, axis=1)
        elif stat == "std":
            out = np.nanstd(chunks, axis=1)
        elif stat in ("coverage", "cov"):
            out = (~np.isnan(chunks)).mean(axis=1)
        elif stat == "mean":
            out = np.nanmean(chunks, axis=1)
        elif stat == "min":
            out = np.nanmin(chunks, axis=1)
        elif stat == "max":
            out = np.nanmax(chunks, axis=1)
        else:
            raise ValueError(f"unknown stat {stat!r}")
    return np.where(np.isnan(out), missing, out).astype(np.float64)


class PyBigToolsHandle:
    """pybigtools: one FFI call per region does I/O, decompression and binning."""
    backend = "pybigtools"
    __slots__ = ("_h", "_own")

    def __init__(self, src):
        global _FILL_KW
        import pybigtools
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self._h = pybigtools.open(str(src), "r") if self._own else src
        if _FILL_KW is None:
            _FILL_KW = _fill_kwarg(self._h)

    def chroms(self):
        return dict(self._h.chroms())

    def stats_array(self, chrom, start, end, *, n_bins=1, stat="mean", exact=True, missing=0.0):
        if stat in _NATIVE:
            return np.asarray(self._h.values(chrom, start, end, bins=n_bins, summary=stat,
                                             exact=exact, **{_FILL_KW: missing}), np.float64)
        return _reduce(self.values(chrom, start, end), n_bins, stat, missing)

    def values(self, chrom, start, end):
        return np.asarray(self._h.values(chrom, start, end, **{_FILL_KW: np.nan}), np.float64)

    def close(self):
        if self._own:
            self._h.close()


class PyBigWigHandle:
    """pyBigWig (libBigWig, C)."""
    backend = "pybigwig"
    __slots__ = ("_h", "_own")

    def __init__(self, src):
        import pyBigWig
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self._h = pyBigWig.open(str(src)) if self._own else src

    def chroms(self):
        return dict(self._h.chroms())

    def stats_array(self, chrom, start, end, *, n_bins=1, stat="mean", exact=True, missing=0.0):
        if stat in _NATIVE or stat in ("std", "coverage", "sum"):
            kind = {"sum": "sum", "coverage": "coverage", "cov": "coverage"}.get(stat, stat)
            try:
                v = self._h.stats(chrom, int(start), int(end), type=kind, nBins=int(n_bins), exact=exact)
            except RuntimeError:                       # region outside the chromosome
                v = [None] * n_bins
            return np.array([missing if x is None else x for x in v], np.float64)
        return _reduce(self.values(chrom, start, end), n_bins, stat, missing)

    def values(self, chrom, start, end):
        try:
            return np.asarray(self._h.values(chrom, int(start), int(end), numpy=True), np.float64)
        except RuntimeError:
            return np.full(max(int(end) - int(start), 0), np.nan)

    def close(self):
        if self._own:
            self._h.close()


class PythonHandle:
    """The pure-Python reader (struct + zlib + numpy + mmap)."""
    backend = "python"
    __slots__ = ("_h", "_own")

    def __init__(self, src):
        from ._bbi import BigWigReader
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self._h = BigWigReader(src) if self._own else src

    def chroms(self):
        return self._h.chroms()

    def stats_array(self, chrom, start, end, *, n_bins=1, stat="mean", exact=True, missing=0.0):
        if stat in ("mean", "min", "max", "std", "sum", "coverage"):
            return self._h.stats_array(chrom, start, end, n_bins=n_bins, stat=stat, exact=exact,
                                       missing=missing)
        return _reduce(self.values(chrom, start, end), n_bins, stat, missing)

    def values(self, chrom, start, end):
        return self._h.values(chrom, start, end)

    def close(self):
        if self._own:
            self._h.close()


_HANDLES = {"pybigtools": PyBigToolsHandle, "pybigwig": PyBigWigHandle, "python": PythonHandle}


def handle_backend(src):
    """The backend an already-open handle belongs to (None for a path)."""
    if isinstance(src, (PyBigToolsHandle, PyBigWigHandle, PythonHandle)):
        return src.backend
    mod = type(src).__module__ or ""
    if mod.startswith("pyBigWig"):
        return "pybigwig"
    if mod.startswith("pybigtools"):
        return "pybigtools"
    if type(src).__name__ == "BigWigReader":
        return "python"
    return None


def open_bigwig(src, *, backend=None):
    """A uniform handle for a bigWig path or an open pyBigWig / pybigtools handle."""
    if isinstance(src, (PyBigToolsHandle, PyBigWigHandle, PythonHandle)):
        return src
    own = handle_backend(src)
    if own is None and hasattr(src, "stats_array"):     # anything speaking our handle interface
        return src
    if own is not None:                                  # someone else's open handle
        return _HANDLES[own](src)
    return _HANDLES[resolve("bigwig", backend)](src)
