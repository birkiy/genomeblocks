"""3-set Venn overlap of BED-like region files.

Inputs (JSON object, all paths absolute — Node resolves them against the
session drive before calling us):

    {"a": "/path/to/A.bed", "b": "/path/to/B.bed", "c": "/path/to/C.bed"}

Output:

    {
      "only_a": int, "only_b": int, "only_c": int,
      "ab_only": int, "ac_only": int, "bc_only": int,
      "abc":     int,
      "a_total": int, "b_total": int, "c_total": int
    }

Overlap rule: two regions overlap when they share ≥1 bp on the same
chromosome (matches ``bedtools intersect`` default). We only read the first
three columns so both ``.bed`` and ``.narrowPeak`` work.

Implementation is pure Python with a per-chromosome sweep — no pandas
dependency for this small job. Swap in an interval tree if region counts
grow past ~1e5 per set.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Iterable


Interval = tuple[int, int]


def _parse_bed(path: str) -> dict[str, list[Interval]]:
    """Read chrom/start/end from a BED-ish file, skipping comments and malformed rows.

    Returns intervals grouped by chromosome, each chromosome's list **sorted
    by start**. Rows where start >= end are discarded (the demo files contain
    a few such malformed rows — matches bedtools' behaviour).
    """
    by_chrom: dict[str, list[Interval]] = defaultdict(list)
    p = Path(path)
    with p.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if start >= end:
                continue
            by_chrom[chrom].append((start, end))

    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def _any_overlap(
    query: Interval, sorted_intervals: list[Interval], starts: list[int]
) -> bool:
    """True if ``query`` overlaps any interval in the sorted list.

    Intervals are sorted by start; ``starts`` is the parallel list of start
    coordinates (hoisted out of this hot path). We binary-search for the
    first interval whose start is ≥ query.end — every candidate that can
    overlap must appear before that index — and then scan those candidates
    for one whose end > query.start.

    We can't short-circuit the scan because intervals aren't sorted by end,
    but every candidate does have start < query.end, so the scan cost is
    proportional to the number of intervals starting before the query ends,
    not to the whole set.
    """
    if not sorted_intervals:
        return False
    qs, qe = query
    idx = bisect_right(starts, qe - 1)  # last candidate has start <= qe-1
    for i in range(idx - 1, -1, -1):
        if sorted_intervals[i][1] > qs:
            return True
    return False


def _overlap_mask(
    set_a: dict[str, list[Interval]], set_b: dict[str, list[Interval]]
) -> dict[str, list[bool]]:
    """For each chromosome in A, mark which A-intervals overlap any B-interval."""
    mask: dict[str, list[bool]] = {}
    for chrom, intervals in set_a.items():
        b_intervals = set_b.get(chrom, [])
        b_starts = [iv[0] for iv in b_intervals]
        mask[chrom] = [_any_overlap(iv, b_intervals, b_starts) for iv in intervals]
    return mask


def _iter_masks(*masks: dict[str, list[bool]]) -> Iterable[tuple[bool, ...]]:
    """Iterate aligned mask entries across multiple per-chromosome masks.

    All masks must be derived from the same A set, so keys/lengths line up.
    """
    base = masks[0]
    for chrom, values in base.items():
        others = [m.get(chrom, [False] * len(values)) for m in masks[1:]]
        for i, v in enumerate(values):
            yield (v, *(o[i] for o in others))


def run(inputs: dict) -> dict:
    a_path = inputs["a"]
    b_path = inputs["b"]
    c_path = inputs["c"]

    A = _parse_bed(a_path)
    B = _parse_bed(b_path)
    C = _parse_bed(c_path)

    a_total = sum(len(v) for v in A.values())
    b_total = sum(len(v) for v in B.values())
    c_total = sum(len(v) for v in C.values())

    # For each set, ask "do I overlap X?" for every other set. Those boolean
    # vectors give us the 7 disjoint Venn regions.
    a_hits_b = _overlap_mask(A, B)
    a_hits_c = _overlap_mask(A, C)
    b_hits_a = _overlap_mask(B, A)
    b_hits_c = _overlap_mask(B, C)
    c_hits_a = _overlap_mask(C, A)
    c_hits_b = _overlap_mask(C, B)

    only_a = 0
    only_b = 0
    only_c = 0
    ab_only = 0
    ac_only = 0
    bc_only = 0
    # For the "all three" region we count the A-side (standard convention —
    # region counts are reported relative to one reference set).
    abc = 0

    for hb, hc in _iter_masks(a_hits_b, a_hits_c):
        if hb and hc:
            abc += 1
        elif hb and not hc:
            ab_only += 1
        elif hc and not hb:
            ac_only += 1
        else:
            only_a += 1

    for ha, hc in _iter_masks(b_hits_a, b_hits_c):
        if not ha and not hc:
            only_b += 1
        elif hc and not ha:
            bc_only += 1  # B∩C regions from B's perspective
        # B∩A (covered by ab_only above from A's side) and A∩B∩C (covered
        # by abc from A's side) are intentionally skipped here to keep the
        # Venn regions disjoint in the reported totals.

    for ha, hb in _iter_masks(c_hits_a, c_hits_b):
        if not ha and not hb:
            only_c += 1
        # A∩C and B∩C and A∩B∩C are counted from the other sets' sides.

    return {
        "only_a": only_a,
        "only_b": only_b,
        "only_c": only_c,
        "ab_only": ab_only,
        "ac_only": ac_only,
        "bc_only": bc_only,
        "abc": abc,
        "a_total": a_total,
        "b_total": b_total,
        "c_total": c_total,
    }
