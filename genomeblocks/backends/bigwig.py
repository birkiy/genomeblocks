"""bigWig backends: one handle interface over pybigtools, pyBigWig and a
pure-Python reader.

    h = open_bigwig("signal.bw", backend=None)     # path or an open handle
    h.chroms()                                     # {chrom: size}
    h.stats_array(chrom, start, end, n_bins=100, stat="mean", exact=True, missing=0.0)
    h.values(chrom, start, end)                    # per base, NaN where no data
    h.close()

``stat`` is one of mean, min, max, sum, std (population), coverage. pybigtools
(Rust) is the default; pyBigWig (C, libBigWig) and the pure-Python reader
(numpy + mmap, no compiled dependency) give the same numbers: bin ``b`` of a
window of ``n`` bases is ``[floor(n*b/n_bins), floor(n*(b+1)/n_bins))`` for
every engine, and bins outside the chromosome are ``missing``. An already-open
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
    """Bin per-base values (NaN = no data) for stats the engines lack natively.

    Bin ``b`` covers bases ``[floor(n*b/n_bins), floor(n*(b+1)/n_bins))`` — the
    integer edges libBigWig and pybigtools use — so every backend bins alike,
    whatever the bin size. ``std`` is the population standard deviation;
    ``coverage`` is the fraction of the bin's bases with data; a bin without
    data is ``missing``.
    """
    n = vals.size
    if n == 0 or n_bins <= 0:
        return np.full(max(n_bins, 1), missing, dtype=np.float64)
    edges = (np.arange(n_bins + 1, dtype=np.int64) * n) // n_bins
    lo, hi = edges[:-1], edges[1:]
    width = hi - lo
    finite = np.isfinite(vals)
    v = np.where(finite, vals, 0.0)
    cs = np.concatenate([[0.0], np.cumsum(v)])
    cc = np.concatenate([[0], np.cumsum(finite)])
    k = cc[hi] - cc[lo]                                        # bases with data per bin
    s = cs[hi] - cs[lo]
    with np.errstate(all="ignore"):
        if stat == "sum":
            out = s
        elif stat == "mean":
            out = s / k
        elif stat == "std":
            cs2 = np.concatenate([[0.0], np.cumsum(v * v)])
            m = s / k
            out = np.sqrt(np.maximum((cs2[hi] - cs2[lo]) / k - m * m, 0.0))
        elif stat in ("coverage", "cov"):
            out = k / width
        elif stat in ("min", "max"):
            fill = np.inf if stat == "min" else -np.inf
            w = np.where(finite, vals, fill)
            red = np.minimum.reduceat if stat == "min" else np.maximum.reduceat
            out = red(w, np.minimum(lo, n - 1))
        else:
            raise ValueError(f"unknown stat {stat!r}; use mean, min, max, sum, std or coverage")
    out = np.asarray(out, np.float64)
    return np.where((k > 0) & (width > 0) & np.isfinite(out), out, missing)


class _Handle:
    """What the three handles share: a clipped, padded path for windows that
    leave the chromosome (or name one the file lacks), so every engine answers
    the same way — ``missing`` per bin, NaN per base."""
    __slots__ = ("_h", "_own", "_sizes")

    def chroms(self):
        if self._sizes is None:
            self._sizes = dict(self._h.chroms())
        return self._sizes

    def _inside(self, chrom, start, end):
        size = self.chroms().get(chrom)
        return size is not None and 0 <= start and end <= size

    def stats_array(self, chrom, start, end, *, n_bins=1, stat="mean", exact=True, missing=0.0):
        start, end = int(start), int(end)
        if self._inside(chrom, start, end) and n_bins <= end - start:
            return self._stats(chrom, start, end, n_bins, stat, exact, missing)
        return _reduce(self.values(chrom, start, end), n_bins, stat, missing)   # edges / tiny windows

    def values(self, chrom, start, end):
        start, end = int(start), int(end)
        n = max(end - start, 0)
        if self._inside(chrom, start, end):
            return self._values(chrom, start, end)
        out = np.full(n, np.nan)
        size = self.chroms().get(chrom)
        if size is not None and n:
            a, b = max(start, 0), min(end, size)
            if b > a:
                out[a - start:b - start] = self._values(chrom, a, b)
        return out

    def close(self):
        if self._own:
            self._h.close()


class PyBigToolsHandle(_Handle):
    """pybigtools: one FFI call per region does I/O, decompression and binning."""
    backend = "pybigtools"
    __slots__ = ()

    def __init__(self, src):
        global _FILL_KW
        import pybigtools
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self._h = pybigtools.open(str(src), "r") if self._own else src
        self._sizes = None
        if _FILL_KW is None:
            _FILL_KW = _fill_kwarg(self._h)

    def _stats(self, chrom, start, end, n_bins, stat, exact, missing):
        if stat in _NATIVE:
            return np.asarray(self._h.values(chrom, start, end, bins=n_bins, summary=stat,
                                             exact=exact, **{_FILL_KW: missing}), np.float64)
        return _reduce(self._values(chrom, start, end), n_bins, stat, missing)

    def _values(self, chrom, start, end):
        return np.asarray(self._h.values(chrom, start, end, **{_FILL_KW: np.nan}), np.float64)


class PyBigWigHandle(_Handle):
    """pyBigWig (libBigWig, C)."""
    backend = "pybigwig"
    __slots__ = ()

    def __init__(self, src):
        import pyBigWig
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self._h = pyBigWig.open(str(src)) if self._own else src
        self._sizes = None

    def _stats(self, chrom, start, end, n_bins, stat, exact, missing):
        if stat in _NATIVE or stat in ("coverage", "cov", "sum"):
            kind = {"sum": "sum", "coverage": "coverage", "cov": "coverage"}.get(stat, stat)
            v = self._h.stats(chrom, start, end, type=kind, nBins=int(n_bins), exact=exact)
            return np.array([missing if x is None else x for x in v], np.float64)
        return _reduce(self._values(chrom, start, end), n_bins, stat, missing)   # std: population

    def _values(self, chrom, start, end):
        return np.asarray(self._h.values(chrom, start, end, numpy=True), np.float64)


class PythonHandle(_Handle):
    """The pure-Python reader (struct + zlib + numpy + mmap)."""
    backend = "python"
    __slots__ = ()

    def __init__(self, src):
        from ._bbi import BigWigReader
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self._h = BigWigReader(src) if self._own else src
        self._sizes = None

    def _stats(self, chrom, start, end, n_bins, stat, exact, missing):
        if (end - start) % max(n_bins, 1) == 0 and stat in ("mean", "min", "max", "std", "sum", "coverage"):
            return self._h.stats_array(chrom, start, end, n_bins=n_bins, stat=stat, exact=exact,
                                       missing=missing)
        return _reduce(self._values(chrom, start, end), n_bins, stat, missing)   # fractional bins

    def _values(self, chrom, start, end):
        return self._h.values(chrom, start, end)


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
        if backend is not None and resolve("bigwig", backend) != own:
            raise ValueError(f"the track is an open {own} handle, so backend={backend!r} cannot apply: "
                             f"pass the file path instead, or drop backend=")
        return _HANDLES[own](src)
    return _HANDLES[resolve("bigwig", backend)](src)
