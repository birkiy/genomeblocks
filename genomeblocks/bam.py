"""BAM pileup extraction for the region browser.

Kept matplotlib-free (drawing lives in ``browserview.py``) and ``pysam`` is
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
from typing import Optional, Union

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
    from .locus import parse_region
    chrom, start, end = parse_region(region)
    return pileup_counts(bam, chrom, start, end, min_baseq=min_baseq).sum(axis=0)


# ── BAM -> bigWig ─────────────────────────────────────────────────────────────

EXCLUDE = 0x4 | 0x100 | 0x200 | 0x800       # unmapped, secondary, failed QC, supplementary
NORMALIZE = (None, "CPM", "RPKM", "BPM", "RPGC")


def _runs(s, e):
    """Depth over half-open intervals [s, e) on one axis: the runs with depth > 0 as (starts, ends, depth),
    touching runs of equal depth merged (``bedtools genomecov -bg`` on one chromosome). O(n log n) in the
    intervals, whatever the axis length."""
    if len(s) == 0:
        z = np.zeros(0, np.int64)
        return z, z, z
    s, e = np.sort(s), np.sort(e)
    b = np.unique(np.concatenate([s, e]))
    d = np.searchsorted(s, b, side="right") - np.searchsorted(e, b, side="right")
    a, z, d = b[:-1], b[1:], d[:-1]
    keep = d > 0
    a, z, d = a[keep], z[keep], d[keep]
    if len(d) == 0:
        return a, z, d
    brk = np.r_[True, (a[1:] != z[:-1]) | (d[1:] != d[:-1])]
    first = np.flatnonzero(brk)
    last = np.r_[first[1:], len(d)] - 1
    return a[first], z[last], d[first]


def _contig_coverage(bam, contig, length, bin_size, exclude, min_mapq, extend):
    """One contig's coverage: (contig, run starts, run ends, counts, reads counted, bases they cover, the sum
    of the counts over bases or bins). Reads are taken as their aligned blocks (split at deletions and
    introns), or extended to fragments (``extend``)."""
    import itertools
    pysam = _pysam()
    f = pysam.AlignmentFile(bam, "rb")
    try:
        it = f.fetch(contig)
        q = int(min_mapq)
        # the filters are written out in each comprehension: a filter function costs ~0.4 s per million reads
        if extend is None:
            blocks = ([r.get_blocks() for r in it if not (r.flag & exclude) and r.mapping_quality >= q] if q > 0
                      else [r.get_blocks() for r in it if not (r.flag & exclude)])
            n = len(blocks)
            lens = np.fromiter(map(len, blocks), np.int64, n)
            se = np.fromiter(itertools.chain.from_iterable(itertools.chain.from_iterable(blocks)), np.int64,
                             2 * int(lens.sum())).reshape(-1, 2)
            s, e = se[:, 0], se[:, 1]
            rid = np.repeat(np.arange(n), lens)
        elif extend == "pairs":
            rec = np.array([(r.flag, r.reference_start, r.reference_end, r.next_reference_start, r.template_length)
                            for r in it if not (r.flag & exclude) and r.mapping_quality >= q], np.int64).reshape(-1, 5)
            fl, s, e, ms, tl = rec.T
            proper = (fl & 0x1 > 0) & (fl & 0x2 > 0) & (tl != 0)
            first = proper & (fl & 0x40 > 0)
            take = first | ~proper                                  # one interval per fragment, mates of proper pairs once
            left = np.minimum(s, ms)
            s, e = np.where(first, left, s)[take], np.where(first, left + np.abs(tl), e)[take]
            n = len(s)
            rid = np.arange(n)
        else:
            L = int(extend)
            rec = np.array([(r.reference_start, r.reference_end, r.is_reverse) for r in it
                            if not (r.flag & exclude) and r.mapping_quality >= q], np.int64).reshape(-1, 3)
            s, e, rev = rec[:, 0], rec[:, 1], rec[:, 2].astype(bool)
            s, e = np.where(rev, np.minimum(s, e - L), s), np.where(rev, e, np.maximum(e, s + L))
            n = len(s)
            rid = np.arange(n)
    finally:
        f.close()
    s, e = np.clip(s, 0, length), np.clip(e, 0, length)
    good = e > s
    s, e, rid = s[good], e[good], rid[good]
    nbases = int((e - s).sum())
    if bin_size == 1:
        a, z, d = _runs(s, e)
        total = int((d * (z - a)).sum())
    else:
        b0, b1 = s // bin_size, (e - 1) // bin_size + 1                    # bins [b0, b1) each interval touches
        if extend is None and len(rid) > 1:                                 # a read counts once per bin: merge its blocks
            new = np.r_[True, (rid[1:] != rid[:-1]) | (b0[1:] >= b1[:-1])]
            first = np.flatnonzero(new)
            last = np.r_[first[1:], len(new)] - 1
            b0, b1 = b0[first], b1[last]
        a, z, d = _runs(b0, b1)
        total = int((d * (z - a)).sum())
        a, z = a * bin_size, np.minimum(z * bin_size, length)
    return contig, a, z, d, n, nbases, total


def bam_to_bigwig(bam: str, path: str, *, bin_size: int = 1, normalize: Optional[str] = None,
                  extend=None, effective_genome_size: Optional[int] = None, min_mapq: int = 0,
                  ignore_duplicates: bool = False, exclude: int = EXCLUDE, chroms=None,
                  workers: Optional[int] = 1) -> str:
    """A coverage bigWig from a BAM, ``deepTools bamCoverage``'s job.

    Args:
        bam: an indexed BAM (``samtools index``).
        path: the bigWig to write.
        bin_size: 1 (default) writes the depth of every base (``bedtools genomecov -bg``, megadepth);
            ``b`` > 1 writes, per bin of ``b`` bases from each chromosome's start, the number of reads that
            overlap it (deepTools' value; a read counts once per bin).
        normalize: None (raw counts), "CPM" (per million reads counted), "RPKM" (per million reads and
            kilobase of bin), "BPM" (per million of the summed bin values) or "RPGC" (1x average coverage:
            ``effective_genome_size / (reads counted x fragment length)``) — deepTools' formulas.
        extend: None counts each read over its aligned blocks (split at deletions and introns, as deepTools
            does); an int extends every read to that many bases in its direction (reads already longer are
            kept); "pairs" counts properly paired reads once, over the fragment between the mates, and other
            reads over their aligned span.
        effective_genome_size: for "RPGC" (e.g. 2,913,022,398 for GRCh38 reads of 50 bp).
        min_mapq: reads below this mapping quality are left out.
        ignore_duplicates: leave out reads flagged as duplicates (0x400).
        exclude: SAM flags of the reads to leave out; default unmapped, secondary, failed QC and
            supplementary.
        chroms: the chromosomes to write (default every chromosome of the BAM header).
        workers: processes, one chromosome at a time each; ``None`` = half the cores.

    The reads counted are those that pass the filters on the chromosomes written; bases with no read hold
    no data. Returns ``path``.
    """
    import os
    from .signal import _write_bigwig
    if normalize is not None and normalize not in NORMALIZE:
        raise ValueError(f"normalize must be one of {', '.join(map(str, NORMALIZE))}")
    if normalize == "RPGC" and not effective_genome_size:
        raise ValueError('normalize="RPGC" needs effective_genome_size')
    if extend is not None and extend != "pairs" and int(extend) <= 0:
        raise ValueError('extend must be None, a positive fragment length or "pairs"')
    bin_size = int(bin_size)
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")
    pysam = _pysam()
    with pysam.AlignmentFile(str(bam), "rb") as f:
        if not f.has_index():
            raise ValueError(f"{bam}: no index; run `samtools index {bam}` first")
        sizes = dict(zip(f.references, map(int, f.lengths)))
        mapped = {s.contig: s.mapped for s in f.get_index_statistics()}
    if chroms is not None:
        want = [chroms] if isinstance(chroms, str) else list(chroms)
        missing = [c for c in want if c not in sizes]
        if missing:
            raise ValueError(f"not in the BAM: {', '.join(missing[:5])}")
    else:
        want = list(sizes)
    todo = sorted((c for c in want if mapped.get(c, 1) > 0), key=lambda c: -sizes[c])    # longest first
    flags = int(exclude) | (0x400 if ignore_duplicates else 0)
    args = [(str(bam), c, sizes[c], bin_size, flags, int(min_mapq), extend) for c in todo]
    if workers is None:
        workers = max(1, (os.cpu_count() or 2) // 2)
    workers = max(1, min(int(workers), len(args) or 1))
    if workers == 1:
        res = [_contig_coverage(*a) for a in args]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as ex:
            res = list(ex.map(_contig_coverage, *zip(*args)))
    n = sum(r[4] for r in res)
    scale = 1.0
    if normalize in ("CPM", "RPKM", "RPGC") and n == 0:
        raise ValueError(f"no reads counted in {bam}: nothing to normalise by")
    if normalize == "CPM":
        scale = 1e6 / n
    elif normalize == "RPKM":
        scale = 1e9 / (n * bin_size)
    elif normalize == "BPM":
        tot = sum(r[6] for r in res)
        scale = 1e6 / tot if tot else 1.0
    elif normalize == "RPGC":
        frag = int(extend) if extend is not None and extend != "pairs" else sum(r[5] for r in res) / n
        scale = float(effective_genome_size) / (n * frag)
    chunks = [(c, a, z, d * scale) for c, a, z, d, *_ in res]
    return _write_bigwig(path, {c: sizes[c] for c in want}, chunks)
