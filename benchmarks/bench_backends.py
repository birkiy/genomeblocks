#!/usr/bin/env python3
"""Every backend of every family, timed through the genomeblocks API.

The same call, the same data, a different engine behind it — so the numbers
say what ``backend=`` buys or costs, and the agreement column says every
engine gives the same answer (the suite's parity rule). Engines that are not
installed are skipped and listed in the results.

Parts (one per family):
  intervals   A & B (overlap_pairs), merge, nearest at n = 1k .. 1M peaks
  bigwig      Loci.signal: 10k peaks x 1 track, 200 bins, ±3 kb
  motifs      scan_motifs_matrix: 1k windows x 100 JASPAR motifs
  fasta       Loci.sequences: 10k windows of 500 bp
  tables      Loci.make on 1M peaks; Genes.make on the GTF
  graph       Architecture.components / pagerank on the loops graph

Data from ``make_data.py`` (``$GB_BENCH_DATA``). ``GB_BENCH_SIZES=1000,10000``
limits the interval sizes (handy for a smoke run on small data).
"""
from __future__ import annotations

import os

import numpy as np

from common import DATA, Recorder, timeit

import genomeblocks as gb
from genomeblocks import Architecture, Genes, Loci, Pairs
from genomeblocks.motifs import scan_motifs_matrix

SIZES = [int(x) for x in os.environ.get("GB_BENCH_SIZES", "1000,10000,100000,1000000").split(",")]
FLANK, NBINS = 3_000, 200


def installed(family):
    return [b for b in gb.backends.families()[family] if gb.backends.installed(family, b)]


def skipped(rec, family):
    for b in gb.backends.families()[family]:
        if not gb.backends.installed(family, b):
            rec.add(part=family, op="skipped", engine=b, n=0, seconds=None,
                    note=gb.backends._FAMILIES[family][b][1])


def pairs_key(qi, ri):
    return np.lexsort((ri, qi)).shape[0], int(qi.sum()), int(ri.sum())


# ── intervals ────────────────────────────────────────────────────────────────

def part_intervals(rec):
    from genomeblocks.backends.intervals import _MERGE, _NEAREST
    for n in SIZES:
        A = Loci.make(str(DATA / f"peaks_A_{n}.bed"))
        B = Loci.make(str(DATA / f"peaks_B_{n}.bed"))
        rep = 5 if n <= 100_000 else 3
        ref = {}
        for b in installed("intervals"):
            res = {}
            t = timeit(lambda: res.__setitem__("p", A.overlap_pairs(B, backend=b)), repeat=rep,
                       setup=lambda: (A._dirty(), B._dirty()))
            k = pairs_key(*res["p"])
            ref.setdefault("pairs", k)
            rec.add(part="intervals", op="A & B", engine=b, n=n, seconds=t["median"], runs=t["runs"],
                    n_out=k[0], agrees=k == ref["pairs"])
            t = timeit(lambda: res.__setitem__("q", A.overlap_pairs(B, backend=b)), repeat=rep)   # warm index
            rec.add(part="intervals", op="A & B (warm index)", engine=b, n=n, seconds=t["median"], runs=t["runs"])
            if b in _MERGE:
                t = timeit(lambda: res.__setitem__("m", A.merge(backend=b)), repeat=rep)
                k = (len(res["m"]), int(res["m"].starts.sum()))
                ref.setdefault("merge", k)
                rec.add(part="intervals", op="merge", engine=b, n=n, seconds=t["median"], runs=t["runs"],
                        n_out=k[0], agrees=k == ref["merge"])
            if b in _NEAREST and n <= 100_000:
                t = timeit(lambda: res.__setitem__("d", A.nearest(B, backend=b)[1]), repeat=rep)
                k = int(np.abs(res["d"]).sum())
                ref.setdefault("nearest", k)
                rec.add(part="intervals", op="nearest (distances)", engine=b, n=n, seconds=t["median"],
                        runs=t["runs"], agrees=k == ref["nearest"])
    skipped(rec, "intervals")


# ── bigwig ───────────────────────────────────────────────────────────────────

def part_bigwig(rec):
    n = min(10_000, SIZES[-1])
    L = Loci.make(str(DATA / f"peaks_A_{max(s for s in SIZES if s <= 10_000) if any(s <= 10_000 for s in SIZES) else SIZES[0]}.bed"))
    L = L.take(slice(0, n))
    bw = str(DATA / "signal_0.bw")
    ref = None
    for b in installed("bigwig"):
        res = {}
        t = timeit(lambda: res.__setitem__("S", L.signal([bw], n_bins=NBINS, flank=FLANK, backend=b, verbose=False,
                                                           progress=False, dtype=np.float64)), repeat=3)
        S = res["S"]
        if ref is None:
            ref = S
        rec.add(part="bigwig", op=f"signal {len(L):,} x 1 track x {NBINS} bins", engine=b, n=len(L),
                seconds=t["median"], runs=t["runs"], agrees=bool(np.allclose(S, ref, atol=1e-3)),
                max_abs_diff=float(np.abs(S - ref).max()))
    skipped(rec, "bigwig")


# ── motifs ───────────────────────────────────────────────────────────────────

def part_motifs(rec, n_windows=1_000, n_motifs=100):
    fa = str(DATA / "genome.fa")
    jaspar = str(DATA / "jaspar.txt")
    lib = gb.load_motifs(jaspar, format="jaspar16").take(list(range(n_motifs)))
    sizes = gb.Genome.from_fasta(fa).sizes
    rng = np.random.default_rng(0)
    chroms = sorted(sizes)
    c = rng.choice(chroms, n_windows)
    pos = np.array([rng.integers(20_000, sizes[x] - 20_000) for x in c])
    L = Loci.from_frame({"chrom": c, "start": pos, "end": pos + 1}).sort()
    ref = None
    for b in installed("motifs"):
        res = {}
        t = timeit(lambda: res.__setitem__("M", scan_motifs_matrix(L, fa, lib, r=250, threshold=13.0, norm=False,
                                                                     workers=1, backend=b, verbose=False)), repeat=3)
        M = res["M"].to_numpy()
        if ref is None:
            ref = M
        rec.add(part="motifs", op=f"scan {n_windows:,} windows x {n_motifs} motifs", engine=b, n=n_windows,
                seconds=t["median"], runs=t["runs"], agrees=bool((M == ref).all()), hits=int(M.sum()))
    skipped(rec, "motifs")


# ── fasta ────────────────────────────────────────────────────────────────────

def part_fasta(rec, n_windows=10_000):
    fa = str(DATA / "genome.fa")
    sizes = gb.Genome.from_fasta(fa).sizes
    rng = np.random.default_rng(1)
    chroms = sorted(sizes)
    c = rng.choice(chroms, n_windows)
    pos = np.array([rng.integers(1_000, sizes[x] - 1_000) for x in c])
    L = Loci.from_frame({"chrom": c, "start": pos, "end": pos + 500}).sort()
    ref = None
    for b in installed("fasta"):
        if b == "memory":
            src = gb.read_fasta(fa)                       # the whole genome as a dict
            t = timeit(lambda: L.sequences(src), repeat=3)
            seqs = L.sequences(src)
        else:
            t = timeit(lambda: L.sequences(fa, backend=b), repeat=3)
            seqs = L.sequences(fa, backend=b)
        if ref is None:
            ref = seqs
        rec.add(part="fasta", op=f"sequences {n_windows:,} x 500 bp", engine=b, n=n_windows, seconds=t["median"],
                runs=t["runs"], agrees=seqs == ref)
    skipped(rec, "fasta")


# ── tables ───────────────────────────────────────────────────────────────────

def part_tables(rec):
    n = SIZES[-1]
    bed = str(DATA / f"peaks_A_{n}.bed")
    gtf = str(DATA / "genes.gtf")
    ref = {}
    for b in installed("tables"):
        res = {}
        t = timeit(lambda: res.__setitem__("L", Loci.make(bed, keep=True, backend=b)), repeat=3)
        ref.setdefault("bed", res["L"])
        rec.add(part="tables", op=f"Loci.make {n:,} peaks (keep=True)", engine=b, n=n, seconds=t["median"],
                runs=t["runs"], agrees=res["L"].equals(ref["bed"], cols=True))
        t = timeit(lambda: res.__setitem__("G", Genes.make(gtf, backend=b)), repeat=3)
        ref.setdefault("gtf", res["G"])
        same = all(getattr(res["G"], k).equals(getattr(ref["gtf"], k), cols=True)
                   for k in ("genes", "transcripts", "features"))
        rec.add(part="tables", op="Genes.make (GTF)", engine=b, n=len(res["G"].transcripts), seconds=t["median"],
                runs=t["runs"], agrees=same)
    skipped(rec, "tables")


# ── graph ────────────────────────────────────────────────────────────────────

def part_graph(rec):
    n = max(s for s in SIZES if s <= 100_000) if any(s <= 100_000 for s in SIZES) else SIZES[0]
    cre = Loci.make(str(DATA / f"peaks_A_{n}.bed"))
    A = Architecture.make(cre, str(DATA / "loops.bedpe"), r=2500, verbose=False)
    A.ep["w"] = np.random.default_rng(0).uniform(1, 10, A.n_links)
    ref = {}
    for b in installed("graph"):
        res = {}
        t = timeit(lambda: res.__setitem__("c", A.components(backend=b)), repeat=3)
        ref.setdefault("components", res["c"])
        rec.add(part="graph", op=f"components ({A.n_links:,} edges)", engine=b, n=A.n_links, seconds=t["median"],
                runs=t["runs"], agrees=bool((res["c"] == ref["components"]).all()))
        t = timeit(lambda: res.__setitem__("p", A.pagerank("w", backend=b)), repeat=3)
        ref.setdefault("pagerank", res["p"])
        rec.add(part="graph", op=f"pagerank ({A.n_links:,} edges)", engine=b, n=A.n_links, seconds=t["median"],
                runs=t["runs"], agrees=bool(np.allclose(res["p"], ref["pagerank"], atol=1e-6)))
        t = timeit(lambda: A.graph(backend=b), repeat=3)
        rec.add(part="graph", op="build the engine's graph", engine=b, n=A.n_links, seconds=t["median"], runs=t["runs"])
    skipped(rec, "graph")


PARTS = {"intervals": part_intervals, "bigwig": part_bigwig, "motifs": part_motifs, "fasta": part_fasta,
         "tables": part_tables, "graph": part_graph}

if __name__ == "__main__":
    import sys
    rec = Recorder("backends")
    for name in (sys.argv[1:] or PARTS):
        print(f"--- {name} ---", flush=True)
        PARTS[name](rec)
    rec.save(sizes=SIZES, backends={f: installed(f) for f in PARTS})
