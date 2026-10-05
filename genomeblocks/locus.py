"""Locus: one genomic interval.

Tables (:class:`~genomeblocks.Loci`, :class:`~genomeblocks.Genes`) hold
intervals as columns; a ``Locus`` is what you get when you look at one row
(``L[i]`` is a :class:`~genomeblocks.loci.LocusView`, a ``Locus`` whose fields
read the columns) and what you pass to say "this region".
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Tuple

_REGION = re.compile(r"^\s*([^:\s]+)\s*:\s*([\d,._]+)\s*([kKmM][bB])?\s*[-–]\s*([\d,._]+)\s*([kKmM][bB])?\s*$")


@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    strand: str = "."

    @property
    def uid(s) -> str:
        """``chrom:start-end(strand)``; :meth:`from_uid` reverses it."""
        return f"{s.chrom}:{s.start}-{s.end}({s.strand})"

    @classmethod
    def from_uid(cls, uid: str) -> "Locus":
        """Rebuild a Locus from ``chrom:start-end(strand)`` or ``chrom:start-end``."""
        body, _, strand = uid.partition("(")
        chrom, _, span = body.rpartition(":")
        start, _, end = span.partition("-")
        if not chrom or not end:
            raise ValueError(f"not a locus uid: {uid!r}")
        return cls(chrom, int(start), int(end), strand.rstrip(")") or ".")

    @classmethod
    def parse(cls, region) -> "Locus":
        """A Locus from a region in any common spelling: ``'chr1:1,000-2,000'``,
        ``'chr8:127.7-128.1 Mb'``, ``'chr1:5kb-10kb'``, a ``(chrom, start, end)``
        tuple, or anything with ``chrom/start/end`` attributes."""
        chrom, start, end = parse_region(region)
        return cls(chrom, start, end, getattr(region, "strand", ".") or ".")

    @property
    def length(l) -> int:
        return l.end - l.start

    @property
    def center(l) -> int:
        return (l.start + l.end) // 2

    def __eq__(s, o: object) -> bool:
        return isinstance(o, Locus) and s.uid == o.uid

    def __ne__(s, o: object) -> bool:
        return not s.__eq__(o)

    def __hash__(s) -> int:
        return hash(s.uid)

    def __lt__(s, o: "Locus") -> bool:
        return NotImplemented if not isinstance(o, Locus) else \
            s.start < o.start if s.chrom == o.chrom else s.chrom < o.chrom

    def __le__(s, o: "Locus") -> bool:
        return NotImplemented if not isinstance(o, Locus) else \
            s.start <= o.start if s.chrom == o.chrom else s.chrom < o.chrom

    def __gt__(s, o: "Locus") -> bool:
        return NotImplemented if not isinstance(o, Locus) else \
            s.start > o.start if s.chrom == o.chrom else s.chrom > o.chrom

    def __ge__(s, o: "Locus") -> bool:
        return NotImplemented if not isinstance(o, Locus) else \
            s.start >= o.start if s.chrom == o.chrom else s.chrom > o.chrom

    def distance_to(s, o: "Locus") -> int:
        """|center - other.center| on the same chromosome (NotImplemented across)."""
        return abs(s.center - o.center) if s.chrom == o.chrom else NotImplemented

    def overlaps(s, o: "Locus") -> bool:
        """Same chromosome and intersecting half-open spans."""
        return (s.chrom == o.chrom) and not (s.end <= o.start or s.start >= o.end)

    def sequence(s, fasta, r=None) -> str:
        """The bases of this locus (or ``center ± r``) from a FASTA path, a
        ``{chrom: str}`` dict, or an open pyfaidx / pysam / Biopython handle."""
        from .backends.fasta import open_fasta
        a, b = (s.center - r, s.center + r) if r is not None else (s.start, s.end)
        return open_fasta(fasta).fetch(s.chrom, max(0, a), b)

    def copy(s) -> "Locus":
        return Locus(s.chrom, s.start, s.end, s.strand)


_REGIONS_RE = re.compile(r"[^\s:]+:\s*[\d.,_]+\s*(?:[kKmM][bB])?\s*[-–]\s*[\d.,_]+\s*(?:[kKmM][bB])?")


def parse_regions(text) -> list:
    """Every region in ``text`` as (chrom, start, end): ``'chr8:127.7-128.0 Mb chr1:1-2 Mb'``
    holds two (a split view). Tuples, Locus-like objects and lists of any of
    these are accepted too."""
    if isinstance(text, (list, tuple)) and not (len(text) == 3 and isinstance(text[0], str)
                                                and not isinstance(text[1], str)
                                                and str(text[1]).lstrip("-").isdigit()):
        return [r for x in text for r in parse_regions(x)]
    if not isinstance(text, str):
        return [parse_region(text)]
    found = _REGIONS_RE.findall(text)
    if not found:
        return [parse_region(text)]
    return [parse_region(x) for x in found]


def parse_region(region) -> Tuple[str, int, int]:
    """(chrom, start, end) from a region string, a tuple, or a Locus-like object.

    Strings accept thousands separators and kb / Mb units on either number:
    ``'chr1:1,000-2,000'``, ``'chr8:127.7-128.1 Mb'``, ``'chr2:5kb-12kb'``.
    """
    if hasattr(region, "chrom") and hasattr(region, "start"):
        return str(region.chrom), int(region.start), int(region.end)
    if isinstance(region, (tuple, list)) and len(region) == 3:
        return str(region[0]), int(region[1]), int(region[2])
    if isinstance(region, str):
        m = _REGION.match(region)
        if m:
            chrom, a, ua, b, ub = m.groups()
            unit = ub or ua
            scale = {"kb": 1_000, "mb": 1_000_000}.get((unit or "").lower(), 1)
            sa = {"kb": 1_000, "mb": 1_000_000}.get((ua or unit or "").lower(), 1)

            def num(x, k):
                return int(round(float(x.replace(",", "").replace("_", "")) * k))
            return chrom, num(a, sa), num(b, scale)
    raise ValueError(f"cannot parse region {region!r}; use 'chr1:1,000-2,000', 'chr1:1.2-1.5 Mb' "
                     f"or (chrom, start, end)")
