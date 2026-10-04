"""HiChIP short-range tracks: the ChIP-like signal hidden in HiChIP pairs.

Ligation pairs closer than ~1 kb are mostly undigested ChIP fragments, so
their 5' read ends behave like ChIP-seq reads. This module does, in one
process, what the usual shell recipe does with awk, sort, bedtools and a
bedGraph converter:

    allValidPairs ─▶ cis pairs ≤ max_dist ─▶ both 5' ends (stranded)
                 ─▶ BED6 for MACS3                       (write_bed / macs3)
                 ─▶ ends extended to 147 bp fragments ─▶ coverage ─▶ bigWig

Peak calling stays with MACS3 (``macs3()`` is a thin wrapper around the CLI).
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .genome import Genome, default_genome
from .loci import Loci

# HiC-Pro allValidPairs: readID chr1 pos1 strand1 chr2 pos2 strand2 size [frag1 frag2 mapq1 mapq2 ...]
_AVP = {"chr1": 1, "pos1": 2, "strand1": 3, "chr2": 4, "pos2": 5, "strand2": 6}


def shortrange_ends(pairs: str, max_dist: int = 1000, *, columns: Optional[Dict[str, int]] = None,
                    genome: Optional[Genome] = None) -> Loci:
    """5' read ends of every cis pair with |pos2 - pos1| <= ``max_dist``.

    Both ends of each pair become one row of a stranded, 1-bp Loci
    (BED start = pos - 1). ``columns`` maps chr1/pos1/strand1/chr2/pos2/strand2
    to 0-based column numbers (defaults: HiC-Pro allValidPairs; for 4DN
    .pairs use {chr1: 1, pos1: 2, chr2: 3, pos2: 4, strand1: 5, strand2: 6}).
    The file is scanned lazily, so only the short-range rows are held in memory.
    """
    import polars as pl
    c = dict(_AVP, **(columns or {}))
    names = {i: f"column_{i + 1}" for i in c.values()}
    lf = pl.scan_csv(pairs, separator="\t", has_header=False, comment_prefix="#", quote_char=None,
                     infer_schema_length=0)
    col = lambda k: pl.col(names[c[k]])
    lf = (lf.select(col("chr1").alias("c1"), col("pos1").cast(pl.Int64).alias("p1"), col("strand1").alias("s1"),
                    col("chr2").alias("c2"), col("pos2").cast(pl.Int64).alias("p2"), col("strand2").alias("s2"))
            .filter((pl.col("c1") == pl.col("c2")) & ((pl.col("p2") - pl.col("p1")).abs() <= max_dist)))
    df = lf.collect(engine="streaming")
    g = genome or default_genome()
    uniq = pl.concat([df["c1"], df["c2"]]).unique().to_list()
    code = lambda col: df[col].replace_strict(uniq, [g._add(u) for u in uniq], return_dtype=pl.Int32).to_numpy()
    codes = np.concatenate([code("c1"), code("c2")])
    pos = np.concatenate([df["p1"].to_numpy(), df["p2"].to_numpy()])
    strand = np.concatenate([(df["s1"] == "-").to_numpy(), (df["s2"] == "-").to_numpy()])
    return Loci(codes, pos - 1, pos, np.where(strand, 2, 1).astype(np.int8), genome=g).sort()


def write_bed(ends: Loci, path: str) -> str:
    """BED6 of the ends (the MACS3 input)."""
    import polars as pl
    names = pl.Series(ends.genome.names)
    pl.DataFrame({"c": names.gather(ends.codes), "s": ends.starts, "e": ends.ends, "n": ".", "v": 0,
                  "d": pl.Series([".", "+", "-"]).gather(ends.strands.astype(np.int64))}).write_csv(
        path, separator="\t", include_header=False)
    return path


def fragments(ends: Loci, extsize: int = 147, chrom_sizes: Optional[Dict[str, int]] = None) -> Loci:
    """Each 5' end extended ``extsize`` bp in its read direction (clipped to the chromosome)."""
    plus = ends.strands != 2
    s = np.where(plus, ends.starts, ends.ends - extsize)
    e = np.where(plus, ends.starts + extsize, ends.ends)
    s = np.maximum(s, 0)
    if chrom_sizes:
        lim = np.array([chrom_sizes.get(n, np.iinfo(np.int64).max) for n in ends.genome.names], np.int64)
        e = np.minimum(e, lim[ends.codes])
    return Loci(ends.codes, s, e, ends.strands, genome=ends.genome)


def coverage(L: Loci):
    """Per-base coverage as runs: yields (chrom, start, end, depth) for depth > 0,
    adjacent equal depths merged (the same records as ``bedtools genomecov -bg``)."""
    names = L.genome.names
    for c in sorted(np.unique(L.codes), key=lambda x: names[x]):
        m = L.codes == c
        st, en = np.sort(L.starts[m]), np.sort(L.ends[m])
        b = np.unique(np.concatenate([st, en]))
        depth = np.searchsorted(st, b, side="right") - np.searchsorted(en, b, side="right")
        a, z, d = b[:-1], b[1:], depth[:-1]
        keep = d > 0
        a, z, d = a[keep], z[keep], d[keep]
        if len(d) == 0:
            continue
        brk = np.r_[True, (a[1:] != z[:-1]) | (d[1:] != d[:-1])]           # merge touching equal runs
        first = np.flatnonzero(brk)
        last = np.r_[first[1:], len(d)] - 1
        yield names[c], a[first], z[last], d[first]


def to_bigwig(L: Loci, path: str, chrom_sizes: Dict[str, int]) -> str:
    """Write the coverage of ``L`` (e.g. ``fragments(ends)``) as a bigWig."""
    import pybigtools

    def records():
        for chrom, a, z, d in coverage(L):
            for x in zip([chrom] * len(a), a.tolist(), z.tolist(), d.astype(float).tolist()):
                yield x
    out = pybigtools.open(path, "w")
    out.write({k: int(v) for k, v in sorted(chrom_sizes.items())}, records())
    return path


def to_bedgraph(L: Loci, path: str) -> str:
    import polars as pl
    with open(path, "wb") as f:
        for chrom, a, z, d in coverage(L):
            pl.DataFrame({"c": [chrom] * len(a), "a": a, "z": z, "d": d}).write_csv(
                f, separator="\t", include_header=False)
    return path


def macs3(bed: str, name: str, outdir: str, *, gsize: str = "hs", extsize: int = 147, q: float = 0.01,
          exe: str = "macs3") -> str:
    """``macs3 callpeak`` with the HiChIP short-range settings; returns the narrowPeak path."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    subprocess.run([exe, "callpeak", "-t", bed, "-f", "BED", "-g", gsize, "-n", name, "--outdir", outdir,
                    "--nomodel", "--extsize", str(extsize), "-q", str(q), "--keep-dup", "all"],
                   check=True, capture_output=True)
    return str(Path(outdir) / f"{name}_peaks.narrowPeak")


def shortrange_track(pairs: str, out_prefix: str, chrom_sizes: Dict[str, int], *, max_dist: int = 1000,
                     extsize: int = 147, genome: Optional[Genome] = None) -> dict:
    """The whole recipe minus peak calling: ends BED (for MACS3) + coverage bigWig."""
    ends = shortrange_ends(pairs, max_dist, genome=genome)
    bed = write_bed(ends, f"{out_prefix}_shortrange_ends.bed")
    bw = to_bigwig(fragments(ends, extsize, chrom_sizes), f"{out_prefix}_shortrange.bw", chrom_sizes)
    return {"ends": ends, "bed": bed, "bigwig": bw}
