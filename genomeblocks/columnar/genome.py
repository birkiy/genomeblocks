"""Genome: the one chromosome dictionary every table shares.

Each table stores chromosomes as small integer *codes* instead of strings.
Codes only mean something relative to a name list, so all tables in a
session point at the same :class:`Genome` — then a code in ``Loci`` and the
same code in ``Genes`` or ``Architecture`` are the same chromosome, and no
table ever has to re-code another one before comparing.

Codes are append-only: a new chromosome name gets the next free code, and
existing codes never change. Ordering (for ``sort``) uses a separate
natural-sort rank (chr1, chr2, ..., chr10, ..., chrX), recomputed lazily.
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
        """From a ``.chrom.sizes`` path or a ``{chrom: length}`` mapping."""
        if isinstance(sizes, str):
            d = {}
            with open(sizes) as f:
                for line in f:
                    if line.strip() and not line.startswith("#"):
                        c, n = line.split()[:2]
                        d[c] = int(n)
            sizes = d
        return cls(sizes=dict(sizes), name=name)

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
        inv, uniq = pd.factorize(np.asarray(values, dtype=object), sort=False)
        lut = np.fromiter((self._add(str(u)) for u in uniq), np.int32, len(uniq))
        return lut[inv] if len(inv) else np.zeros(0, np.int32)

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


#: the session-wide default; every table uses it unless told otherwise
DEFAULT = Genome(name="default")


def default_genome() -> Genome:
    return DEFAULT


def set_default_genome(g: Genome) -> Genome:
    global DEFAULT
    DEFAULT = g
    return g
