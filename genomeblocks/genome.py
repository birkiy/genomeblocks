"""Genome: chromosome names <-> integer codes, plus sizes.

Each table stores chromosomes as small integer *codes* instead of strings,
relative to its :class:`Genome`. Tables made from one another (``take``,
``slop``, an Architecture's vertices, ...) share their Genome; tables read
separately each get their own, and an operation between two of them
re-codes the second one onto the first (a small lookup table, O(n)). Pass
one Genome to several constructors (``genome=g``) to skip even that.

Codes are append-only: a new chromosome name gets the next free code, and
existing codes never change. Ordering (for ``sort``) uses a separate
natural-sort rank (chr1, chr2, ..., chr10, ..., chrX), so results never
depend on the order in which names were seen.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional

import numpy as np


def _natural_key(name: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


class Genome:
    """Chromosome names <-> integer codes (+ optional sizes)."""

    def __init__(self, names: Iterable[str] = (), sizes: Optional[Dict[str, int]] = None,
                 name: Optional[str] = None):
        self.name = name
        self.names: List[str] = []
        self.code: Dict[str, int] = {}
        self._sizes: Dict[str, int] = dict(sizes or {})
        self._rank = None
        for n in list(names) + [n for n in self._sizes if n not in set(names)]:
            self._add(n)

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def from_sizes(cls, sizes, name: Optional[str] = None) -> "Genome":
        """From chromosome sizes: a ``{chrom: length}`` dict, a ``.chrom.sizes``
        path, a pandas Series, a cooler, or anything with a ``chromsizes`` mapping."""
        return cls(sizes=read_sizes(sizes), name=name)

    @classmethod
    def from_fasta(cls, fasta, name: Optional[str] = None) -> "Genome":
        """Chromosome names and sizes of a FASTA (its .fai is used when present)."""
        from .backends.fasta import open_fasta
        return cls(sizes=open_fasta(fasta).sizes(), name=name)

    @property
    def sizes(self) -> Dict[str, int]:
        """{chrom: length} for the chromosomes whose size is known."""
        return dict(self._sizes)

    def set_sizes(self, sizes) -> "Genome":
        """Add or update chromosome sizes (adds unseen chromosomes)."""
        for k, v in read_sizes(sizes).items():
            self._add(k)
            self._sizes[k] = int(v)
        return self

    def _add(self, n: str) -> int:
        c = self.code.get(n)
        if c is None:
            c = self.code[n] = len(self.names)
            self.names.append(n)
            self._rank = None
        return c

    # ── encode / decode ───────────────────────────────────────────────────
    def encode(self, values) -> np.ndarray:
        """Chromosome names -> int32 codes (unseen names are added)."""
        import pandas as pd
        arr = np.asarray(values, dtype=object)
        if not len(arr):
            return np.zeros(0, np.int32)
        inv, uniq = pd.factorize(arr, sort=False)
        if (inv < 0).any():
            n = int((inv < 0).sum())
            raise ValueError(f"{n} of {len(arr)} rows have a missing chromosome (None / NaN); "
                             f"drop those rows or fill the chrom column first")
        lut = np.fromiter((self._add(str(u)) for u in uniq), np.int32, len(uniq))
        return lut[inv]

    def decode(self, codes) -> np.ndarray:
        """int codes -> object array of names."""
        return np.asarray(self.names, dtype=object)[np.asarray(codes)]

    def __getitem__(self, name: str) -> int:
        return self.code[name]

    def __contains__(self, name: str) -> bool:
        return name in self.code

    def __len__(self) -> int:
        return len(self.names)

    @property
    def rank(self) -> np.ndarray:
        """rank[code] = position of that chromosome in natural sort order."""
        if self._rank is None or len(self._rank) != len(self.names):
            order = sorted(range(len(self.names)), key=lambda i: _natural_key(self.names[i]))
            r = np.empty(len(self.names), np.int64)
            r[order] = np.arange(len(order))
            self._rank = r
        return self._rank

    def size(self, name: str) -> Optional[int]:
        return self._sizes.get(name)

    def __repr__(self) -> str:
        head = ", ".join(self.names[:5]) + (", ..." if len(self.names) > 5 else "")
        return f"Genome({self.name or ''}{' · ' if self.name else ''}{len(self)} chroms: {head})"


def read_sizes(chromsizes) -> Dict[str, int]:
    """{chrom: length} from a dict, a .chrom.sizes path, a pandas Series, a
    Genome, a cooler, or anything with a ``chromsizes`` mapping."""
    if isinstance(chromsizes, Genome):
        return dict(chromsizes._sizes)
    if isinstance(chromsizes, str) or hasattr(chromsizes, "__fspath__"):
        out = {}
        with open(chromsizes) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    c, n = line.split()[:2]
                    out[c] = int(n)
        return out
    if hasattr(chromsizes, "chromsizes"):
        chromsizes = chromsizes.chromsizes
    if hasattr(chromsizes, "items"):
        return {str(k): int(v) for k, v in chromsizes.items()}
    raise TypeError(f"cannot read chromosome sizes from {type(chromsizes).__name__}")
