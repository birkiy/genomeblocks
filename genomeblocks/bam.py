"""BAM pileup extraction for the region browser.

Kept matplotlib-free (drawing lives in ``browser.py``) and ``pysam`` is
imported lazily so ``import genomeblocks`` stays cheap and pysam remains an
optional dependency — only code paths that actually touch a ``.bam`` pull it
in, mirroring how ``signal.py`` defers ``pybigtools``.

The one hot call is :func:`pileup_counts` → ``pysam``'s ``count_coverage``,
one C-level pass that returns per-base A/C/G/T depth over the window. Total
coverage is the column sum; the reference base (from :func:`reference_seq`)
tells the drawer which allele is the match (gray) and which are mismatches
(colored), exactly the IGV coverage-track look.
"""
from __future__ import annotations
from typing import Tuple, Union

import numpy as np

# A/C/G/T in the order ``pysam.count_coverage`` returns them.
BASES = ("A", "C", "G", "T")


# ── lazy handles ─────────────────────────────────────────────────────────────

def _pysam():
    try:
        import pysam
    except ImportError as e:  # pragma: no cover - trivial guard
        raise ImportError(
            "BAM support needs pysam — install it with "
            "`pip install genomeblocks[bam]` (or `conda install -c bioconda pysam`)."
        ) from e
    return pysam


def _open_bam(path: str):
    return _pysam().AlignmentFile(path, "rb")


def _open_fasta(path: str):
    return _pysam().FastaFile(path)


def _resolve_contig(names, chrom: str) -> str:
    """Return the contig name present in ``names`` matching ``chrom``.

    Tolerates the common ``chr1`` vs ``1`` mismatch between a BAM/FASTA and the
    requested region so the viewer doesn't silently render an empty track.
    """
    names = set(names)
    if chrom in names:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    if alt in names:
        return alt
    raise KeyError(f"contig {chrom!r} not found (have e.g. {sorted(names)[:3]}…)")


# ── extraction ───────────────────────────────────────────────────────────────

def pileup_counts(bam: Union[str, "object"], chrom: str, start: int, end: int,
                  *, min_baseq: int = 15) -> np.ndarray:
    """Per-base A/C/G/T read depth over ``[start, end)``.

    Returns an ``int64`` array of shape ``(4, end - start)`` (rows = A, C, G, T).
    Bases with quality < ``min_baseq`` are skipped, and unmapped / secondary /
    QC-fail / duplicate reads are excluded (``count_coverage``'s ``'all'``
    read filter) — the same reads IGV drops from its coverage track.

    ``bam`` may be a path (opened and closed here) or an already-open
    ``pysam.AlignmentFile``.
    """
    opened = isinstance(bam, str)
    h = _open_bam(bam) if opened else bam
    try:
        contig = _resolve_contig(h.references, chrom)
        span = end - start
        clen = h.get_reference_length(contig)
        lo, hi = max(0, start), min(end, clen)
        if hi <= lo:
            return np.zeros((4, span), dtype=np.int64)
        cov = h.count_coverage(contig, lo, hi, quality_threshold=min_baseq,
                               read_callback="all")
        arr = np.vstack([np.asarray(c, dtype=np.int64) for c in cov])
    finally:
        if opened:
            h.close()
    # Pad back to the full requested window if we clamped to the contig edges.
    if lo != start or hi != end:
        out = np.zeros((4, span), dtype=np.int64)
        out[:, lo - start:hi - start] = arr
        return out
    return arr


def reference_seq(fasta: Union[str, "object"], chrom: str, start: int,
                  end: int) -> str:
    """Uppercased reference bases over ``[start, end)``.

    Shorter than ``end - start`` only near a contig edge; the drawer pads with
    ``N`` in that case. ``fasta`` may be a path or an open ``pysam.FastaFile``.
    """
    opened = isinstance(fasta, str)
    h = _open_fasta(fasta) if opened else fasta
    try:
        contig = _resolve_contig(h.references, chrom)
        clen = h.get_reference_length(contig)
        lo, hi = max(0, start), min(end, clen)
        if hi <= lo:
            return ""
        pre = "N" * (lo - start)
        post = "N" * (end - hi)
        return pre + h.fetch(contig, lo, hi).upper() + post
    finally:
        if opened:
            h.close()


def coverage(bam: Union[str, "object"], region, *, min_baseq: int = 15
             ) -> np.ndarray:
    """Total per-base coverage over ``region`` (an ``int64`` array).

    ``region`` is anything the browser accepts — a ``Locus``, ``(chrom, start,
    end)`` tuple, or ``'chr1:1,000-2,000'`` string. Convenience wrapper around
    :func:`pileup_counts` for scripting outside the browser.
    """
    from .browser import _parse_region
    chrom, start, end = _parse_region(region)
    return pileup_counts(bam, chrom, start, end, min_baseq=min_baseq).sum(axis=0)
