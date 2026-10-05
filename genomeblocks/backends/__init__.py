"""Backends: which library does the work.

genomeblocks keeps its data as columnar tables; the heavy lifting on those
tables runs through a *backend*, one per family of operations:

    family      what it does                               default (automatic)        others on request
    ─────────   ─────────────────────────────────────────  ─────────────────────────  ──────────────────────────────────
    intervals   overlap, nearest, merge, point lookups     genomeblocks (numpy)       cgranges, ncls, bioframe, pyranges,
                                                                                      bedtools
    bigwig      reading bigWig signal                      pybigtools                 pybigwig, python
    motifs      scoring motif matrices along sequences     lightmotif                 moods, biopython
    fasta       fetching sequence                          genomeblocks (indexed)     pysam, pyfaidx, memory, biopython
    tables      parsing text tables (BED, GTF, pairs)      polars, else pandas        pandas
    graph       graph algorithms on an Architecture        graph-tool, else scipy     igraph, networkx, scipy

Every default installs with ``pip install genomeblocks``; an engine is only
picked automatically where all candidates give the same answer (polars and
pandas parse to the same columns; graph-tool and scipy give the same
components). Anything else is a request: ``backend=`` on the call, or a
``with`` block for a stretch of code::

    import genomeblocks as gb
    gb.backends()                                  # what is installed, what is used
    cre.intersect(other, backend="pyranges")       # this call
    with gb.use_backend(intervals="bioframe", graph="networkx"):
        ...                                        # this block

A requested backend that is not installed raises an error with the install
command — there is no silent switch to another engine. The test suite runs
every backend against the default and checks they agree;
``benchmarks/bench_backends.py`` times them against each other.
"""
from __future__ import annotations

import importlib.util
from contextlib import contextmanager
from typing import Dict, Optional

# family -> backend -> (module to probe, install hint); None = always available
_FAMILIES: Dict[str, Dict[str, tuple]] = {
    "intervals": {
        "genomeblocks": (None, ""),
        "cgranges": ("cgranges", "conda install -c bioconda cgranges  (or pip install git+https://github.com/lh3/cgranges)"),
        "ncls": ("ncls", "pip install ncls"),
        "bioframe": ("bioframe", "pip install bioframe"),
        "pyranges": ("pyranges", "pip install pyranges"),
        "bedtools": ("pybedtools", "pip install pybedtools  (and bedtools on PATH)"),
    },
    "bigwig": {
        "pybigtools": ("pybigtools", "pip install pybigtools"),
        "pybigwig": ("pyBigWig", "pip install pyBigWig"),
        "python": (None, ""),
    },
    "motifs": {
        "lightmotif": ("lightmotif", "pip install lightmotif"),
        "moods": ("MOODS", "pip install MOODS-python"),
        "biopython": ("Bio", "pip install biopython"),
    },
    "fasta": {
        "genomeblocks": (None, ""),
        "pysam": ("pysam", "pip install pysam"),
        "pyfaidx": ("pyfaidx", "pip install pyfaidx"),
        "memory": (None, ""),
        "biopython": ("Bio", "pip install biopython"),
    },
    "tables": {
        "polars": ("polars", "pip install polars"),
        "pandas": ("pandas", "pip install pandas"),
    },
    "graph": {
        "graph-tool": ("graph_tool", "conda install -c conda-forge graph-tool"),
        "scipy": ("scipy", "pip install scipy"),
        "igraph": ("igraph", "pip install igraph"),
        "networkx": ("networkx", "pip install networkx"),
    },
}

#: what "auto" may pick, in order: only engines that give identical answers
AUTO = {
    "intervals": ["genomeblocks"],
    "bigwig": ["pybigtools", "python"],
    "motifs": ["lightmotif"],
    "fasta": ["genomeblocks"],
    "tables": ["polars", "pandas"],
    "graph": ["graph-tool", "scipy"],
}

_scoped: Dict[str, str] = {}
_installed_cache: Dict[str, bool] = {}

#: aliases accepted anywhere a backend name is
_ALIASES = {"numpy": "genomeblocks", "graph_tool": "graph-tool", "graphtool": "graph-tool",
            "gt": "graph-tool", "nx": "networkx", "pybedtools": "bedtools", "pyBigWig": "pybigwig",
            "bio": "biopython", "Bio": "biopython", "MOODS": "moods", "python-igraph": "igraph",
            "dict": "memory"}


def families():
    """The backend families and their backends, default first."""
    return {f: list(v) for f, v in _FAMILIES.items()}


def installed(family: str, backend: str) -> bool:
    """True when ``backend`` of ``family`` can be imported here."""
    backend = _ALIASES.get(backend, backend)
    mod, _ = _FAMILIES[family][backend]
    if mod is None:
        return True
    if mod not in _installed_cache:
        try:
            _installed_cache[mod] = importlib.util.find_spec(mod.split(".")[0]) is not None
        except (ImportError, ValueError):
            _installed_cache[mod] = False
    return _installed_cache[mod]


def _check_family(family: str):
    if family not in _FAMILIES:
        raise ValueError(f"unknown backend family {family!r}; families: {', '.join(_FAMILIES)}")


def resolve(family: str, backend: Optional[str] = None) -> str:
    """The backend a call uses: ``backend`` if given, else the one chosen by an
    enclosing :func:`use_backend` block, else the automatic default. A named
    backend that is unknown or not installed raises (never a silent switch)."""
    _check_family(family)
    name = backend or _scoped.get(family) or "auto"
    name = _ALIASES.get(name, name)
    if name == "auto":
        for cand in AUTO[family]:
            if installed(family, cand):
                return cand
        need = " or ".join(_FAMILIES[family][c][1] for c in AUTO[family] if _FAMILIES[family][c][1])
        raise ImportError(f"no {family} backend is installed: {need}")
    if name not in _FAMILIES[family]:
        raise ValueError(f"unknown {family} backend {name!r}; choose from: {', '.join(_FAMILIES[family])}")
    if not installed(family, name):
        raise ImportError(f"the {name!r} {family} backend is not installed: {_FAMILIES[family][name][1]}")
    return name


@contextmanager
def use_backend(**choices):
    """Use other backends inside a ``with`` block (``None`` or "auto" = default)::

        with gb.use_backend(intervals="pyranges", bigwig="pybigwig"):
            ...

    Names are checked on entry, so a missing backend fails before any work."""
    previous = {}
    for family, name in choices.items():
        _check_family(family)
        if name not in (None, "auto"):
            resolve(family, name)
    try:
        for family, name in choices.items():
            previous[family] = _scoped.get(family)
            if name in (None, "auto"):
                _scoped.pop(family, None)
            else:
                _scoped[family] = _ALIASES.get(name, name)
        yield
    finally:
        for family, name in previous.items():
            if name is None:
                _scoped.pop(family, None)
            else:
                _scoped[family] = name


def backends():
    """A table of every family and backend: installed, the automatic default,
    the one in use right now, and how to install the missing ones."""
    import pandas as pd
    rows = []
    for family, names in _FAMILIES.items():
        try:
            active = resolve(family)
        except ImportError:
            active = None
        default = next((n for n in AUTO[family] if installed(family, n)), None)
        for name in names:
            rows.append({"family": family, "backend": name, "installed": installed(family, name),
                         "default": name == default, "in use": name == active,
                         "install": "" if installed(family, name) else _FAMILIES[family][name][1]})
    return pd.DataFrame(rows)


def unsupported(family: str, backend: str, op: str, supported) -> NotImplementedError:
    """The error raised when a backend cannot do an operation."""
    return NotImplementedError(f"the {backend!r} {family} backend has no {op}; "
                               f"use one of: {', '.join(supported)}")


# ``genomeblocks.backends`` is this module (it holds the engines); calling it
# lists them, so ``gb.backends()`` works although the name is a subpackage.
import sys as _sys
import types as _types


class _Listing(_types.ModuleType):
    def __call__(self):
        return backends()


_sys.modules[__name__].__class__ = _Listing
