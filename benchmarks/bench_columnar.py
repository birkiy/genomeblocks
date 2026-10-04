#!/usr/bin/env python3
"""Objects vs columns: run genomeblocks' own functions on a columnar Loci.

For each function: the same input as today's object Loci and as the
prototype ColumnarLoci (prototypes/columnar.py). The genomeblocks code is
NOT changed; the columnar Loci reaches it through LocusView objects. Where a
function has an obvious vectorised form, that variant is timed too. Every
pair of outputs is compared for equality.

Parts: memory, access, functions
"""
from __future__ import annotations

import contextlib
import gc
import io
import sys
import tracemalloc

import numpy as np

from common import DATA, Recorder, read_chromsizes, timeit
from prototypes.columnar import ColumnarLoci, annotations_columnar, bin_ranges_columnar

from genomeblocks import Atlas, Genes, Loci


def quiet(fn):
    def run():
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return fn()
    return run


def peak_bytes(fn):
    gc.collect()
    tracemalloc.start()
    obj = fn()
    cur, _ = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return obj, cur


def part_memory(rec):
    print("\n== memory ==")
    path = str(DATA / "peaks_A_1000000.bed")
    o, b_obj = peak_bytes(lambda: Loci.make(path))
    rec.add(part="memory", storage="objects (Loci of Locus)", n=len(o), bytes=int(b_obj),
            per_interval=b_obj / len(o), seconds=0.0)
    del o
    c, b_col = peak_bytes(lambda: ColumnarLoci.make(path))
    rec.add(part="memory", storage="columns (ColumnarLoci)", n=len(c), bytes=int(b_col),
            per_interval=b_col / len(c), seconds=0.0)
    del c
    t = timeit(lambda: Loci.make(path), repeat=3)
    rec.add(part="memory", storage="objects: Loci.make time", n=1_000_000, seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: ColumnarLoci.make(path), repeat=3)
    rec.add(part="memory", storage="columns: make time", n=1_000_000, seconds=t["median"], runs=t["runs"])


def part_access(rec):
    print("\n== access patterns (100k) ==")
    path = str(DATA / "peaks_A_100000.bed")
    O, C = Loci.make(path), ColumnarLoci.make(path)
    n = len(O)
    idx = np.random.default_rng(0).integers(0, n, n).tolist()
    for name, fn_o, fn_c in (
        ("loop: for l in loci: total += l.start", lambda: sum(l.start for l in O), lambda: sum(l.start for l in C)),
        ("random access: loci[i].end", lambda: [O[i].end for i in idx], lambda: [C[i].end for i in idx]),
        ("build the uid -> index map", lambda: {l.uid: i for i, l in enumerate(O)}, lambda: (setattr(C, "_uids", None), C.uids)),
    ):
        t = timeit(fn_o, repeat=3)
        rec.add(part="access", pattern=name, storage="objects", seconds=t["median"], runs=t["runs"])
        t = timeit(fn_c, repeat=3)
        rec.add(part="access", pattern=name, storage="columns via LocusView", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: int(C.starts.sum()), repeat=5)
    rec.add(part="access", pattern="loop: for l in loci: total += l.start", storage="columns, vectorised",
            seconds=t["median"], runs=t["runs"])


def same_loci(a, b):
    return sorted((l.chrom, l.start, l.end, l.strand) for l in a) == sorted((l.chrom, l.start, l.end, l.strand) for l in b)


def add(rec, fn, storage, t, **kw):
    rec.add(part="functions", function=fn, storage=storage, seconds=t["median"], runs=t["runs"], **kw)


def part_functions(rec):
    print("\n== genomeblocks functions ==")
    pa, pb, pab = (str(DATA / f"peaks_{x}_100000.bed") for x in ("A", "B", "AB"))
    OA, OB, OAB = Loci.make(pa), Loci.make(pb), Loci.make(pab)
    CA, CB, CAB = ColumnarLoci.make(pa), ColumnarLoci.make(pb), ColumnarLoci.make(pab)
    OB.cgr; CB.cgr
    res = {}

    # set algebra
    for name, fo, fc in (("A & B", lambda: OA & OB, lambda: CA & CB),
                         ("A - B", lambda: OA - OB, lambda: CA - CB),
                         ("merge (200k)", lambda: OAB.merge(), lambda: CAB.merge()),
                         ("make → slop → sort → merge", lambda: Loci.make(pab).slop(100).sort().merge(),
                          lambda: ColumnarLoci.make(pab).slop(100).sort().merge())):
        t = timeit(lambda: res.__setitem__("o", fo()), repeat=3)
        add(rec, name, "objects", t)
        t = timeit(lambda: res.__setitem__("c", fc()), repeat=3)
        add(rec, name, "columns, vectorised", t, same=same_loci(res["o"], res["c"]))

    # Genes
    genes = quiet(lambda: Genes.make(str(DATA / "genes.gtf")))()
    genes.annot
    t = timeit(lambda: res.__setitem__("o", genes.annotations(OA)), repeat=3)
    add(rec, "Genes.annotations (100k CREs)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", genes.annotations(CA)), repeat=3)
    add(rec, "Genes.annotations (100k CREs)", "columns, today's code", t, same=res["o"].equals(res["c"]))
    t = timeit(lambda: res.__setitem__("v", annotations_columnar(genes, CA)), repeat=3)
    add(rec, "Genes.annotations (100k CREs)", "columns, vectorised", t, same=res["o"].equals(res["v"]))
    t = timeit(lambda: res.__setitem__("o", genes.nearest_genes(OA)), repeat=3)
    add(rec, "Genes.nearest_genes (100k)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", genes.nearest_genes(CA)), repeat=3)
    key = lambda d: d.sort_values(["Name", "Name_b"]).reset_index(drop=True)
    add(rec, "Genes.nearest_genes (100k)", "columns, today's code", t, same=key(res["o"]).equals(key(res["c"])))

    # Atlas
    A = Atlas.load(str(DATA / "atlas_1kb.npz"))
    OQ, CQ = Loci.make(str(DATA / "atlas_query.bed")), ColumnarLoci.make(str(DATA / "atlas_query.bed"))
    t = timeit(lambda: res.__setitem__("o", A.search(OQ)), repeat=5)
    add(rec, "Atlas.search (20k peaks, 500 tracks)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", A.search(CQ)), repeat=5)
    add(rec, "Atlas.search (20k peaks, 500 tracks)", "columns, today's code", t, same=res["o"].equals(res["c"]))
    t = timeit(lambda: res.__setitem__("o", A._intervals_to_bin_ranges(OQ)), repeat=5)
    add(rec, "Atlas: peaks → bin ranges", "objects", t)
    t = timeit(lambda: res.__setitem__("c", bin_ranges_columnar(A, CQ)), repeat=5)
    add(rec, "Atlas: peaks → bin ranges", "columns, vectorised", t, same=bool(np.array_equal(res["o"], res["c"])))

    # signal
    from genomeblocks.signal import signal
    O5 = Loci(OA[i] for i in range(0, 100_000, 20))
    C5 = ColumnarLoci(O5)
    bw = [str(DATA / "signal_0.bw")]
    kw = dict(n_bins=200, flank=3000, progress=False, verbose=False)
    t = timeit(lambda: res.__setitem__("o", signal(O5, bw, **kw)), repeat=3)
    add(rec, "signal (5k loci, 1 bigWig)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", signal(C5, bw, **kw)), repeat=3)
    add(rec, "signal (5k loci, 1 bigWig)", "columns, today's code", t, same=bool(np.array_equal(res["o"], res["c"])))

    # motifs
    from bench_motifs import JASPAR, windows, subset_file
    from genomeblocks.motifs import make_genome, scan_motifs_matrix
    genome = make_genome(str(DATA / "genome.fa"))
    OW = windows(genome, 1000)
    CW = ColumnarLoci(OW)
    path = subset_file(100)
    mk = dict(motif_format="jaspar16", r=250, threshold=13.0, norm=False, workers=1, verbose=False)
    t = timeit(lambda: res.__setitem__("o", scan_motifs_matrix(OW, genome, path, **mk)), repeat=3)
    add(rec, "scan_motifs_matrix (1k windows × 100)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", scan_motifs_matrix(CW, genome, path, **mk)), repeat=3)
    add(rec, "scan_motifs_matrix (1k windows × 100)", "columns, today's code", t, same=res["o"].equals(res["c"]))

    # bedpe
    from genomeblocks.bedpe import count_pairs, pair_to_bed, read_bedpe
    loops = read_bedpe(str(DATA / "loops.bedpe"), verbose=False)
    t = timeit(lambda: res.__setitem__("o", pair_to_bed(OA, loops, verbose=False)), repeat=3)
    add(rec, "pair_to_bed (50k loops)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", pair_to_bed(CA, loops, verbose=False)), repeat=3)
    add(rec, "pair_to_bed (50k loops)", "columns, today's code", t,
        same=[str(p) for p in res["o"]] == [str(p) for p in res["c"]])
    cs = read_chromsizes()
    OWIN = Loci.tile_genome(cs, 50_000)
    CWIN = ColumnarLoci(OWIN)
    head = str(DATA / "hic_head.pairs")
    t = timeit(lambda: res.__setitem__("o", count_pairs(OWIN, head, format="pairs", verbose=False)), repeat=3)
    add(rec, "count_pairs (200k pairs, 62k windows)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", count_pairs(CWIN, head, format="pairs", verbose=False)), repeat=3)
    add(rec, "count_pairs (200k pairs, 62k windows)", "columns, today's code", t, same=res["o"].equals(res["c"]))

    # Architecture: make + normalize-free annotate (weights set identically)
    from genomeblocks import Architecture

    def build(L):
        G = quiet(lambda: Architecture.make(L, str(DATA / "loops.bedpe"), r=2500, verbose=False))()
        G.ep["n"] = G.new_edge_property("float")
        G.ep["n"].a = (np.arange(G.num_edges()) * 7919 % 1000) / 1000.0
        return G
    t = timeit(lambda: res.__setitem__("o", build(OA)), repeat=3)
    add(rec, "Architecture.make (50k loops)", "objects", t)
    t = timeit(lambda: res.__setitem__("c", build(CA)), repeat=3)
    edges = lambda G: sorted(tuple(sorted((G.vp.uid[e.source()], G.vp.uid[e.target()]))) for e in G.edges())
    add(rec, "Architecture.make (50k loops)", "columns, today's code", t, same=edges(res["o"]) == edges(res["c"]))
    GO, GC = res["o"], res["c"]
    t = timeit(quiet(lambda: GO.annotate(OA, genes, key="n", verbose=False)), repeat=3)
    add(rec, "Architecture.annotate", "objects", t)
    t = timeit(quiet(lambda: GC.annotate(CA, genes, key="n", verbose=False)), repeat=3)
    lab = lambda G: {G.vp.uid[v]: (G.vp.annot[v], G.vp.gene[v]) for v in G.vertices()}
    add(rec, "Architecture.annotate", "columns, today's code", t, same=lab(GO) == lab(GC))


def part_row_access(rec):
    """One row at a time: objects vs views vs the usual table libraries."""
    import pandas as pd
    import polars as pl
    import pyarrow as pa
    p = str(DATA / "peaks_A_100000.bed")
    O, C = Loci.make(p), ColumnarLoci.make(p)
    df = C.to_frame()
    pdf, tb = pl.from_pandas(df), pa.Table.from_pandas(df)
    idx = np.random.default_rng(0).integers(0, len(C), 10_000).tolist()
    for name, fn in (("Locus objects: loci[i].start", lambda: [O[i].start for i in idx]),
                     ("columns via LocusView: loci[i].start", lambda: [C[i].start for i in idx]),
                     ("raw numpy column: starts[i]", lambda: [C.starts[i] for i in idx]),
                     ("polars: df.row(i)", lambda: [pdf.row(i) for i in idx]),
                     ("Arrow: table['Start'][i]", lambda: [tb["Start"][i] for i in idx]),
                     ("pandas: df.at[i, 'Start']", lambda: [df.at[i, "Start"] for i in idx]),
                     ("pandas: df.iloc[i]['Start']", lambda: [df.iloc[i]["Start"] for i in idx])):
        t = timeit(fn, repeat=3)
        rec.add(part="row_access", way=name, seconds=t["median"] / len(idx), runs=[r / len(idx) for r in t["runs"]])
    for name, fn in (("columns → Arrow table", lambda: pa.Table.from_arrays([pa.array(C.starts), pa.array(C.ends)], names=["s", "e"])),
                     ("columns → polars frame", lambda: pl.DataFrame({"s": C.starts, "e": C.ends})),
                     ("columns → pandas frame", lambda: pd.DataFrame({"s": C.starts, "e": C.ends}))):
        t = timeit(fn, repeat=5)
        rec.add(part="row_access", way=name, seconds=t["median"], runs=t["runs"], n=len(C))


def part_annot_breakdown(rec):
    """Where vectorised annotations spends its time when Genes stays objects."""
    from prototypes.npintervals import Intervals, overlaps_any
    genes = quiet(lambda: Genes.make(str(DATA / "genes.gtf")))()
    genes.annot
    CA = ColumnarLoci.make(str(DATA / "peaks_A_100000.bed"))
    res = {}
    t = timeit(lambda: res.__setitem__("iv", {k: Intervals.from_loci(list(v)) for k, v in genes.annot.items()}), repeat=3)
    rec.add(part="annot_breakdown", step="convert the gene index (objects) to arrays", seconds=t["median"], runs=t["runs"])
    a, ivs = CA._iv(), res["iv"]

    def label():
        out = np.full(len(CA), "Intergenic", dtype=object)
        for key, name in (("body", "Intronic"), ("exon", "Exonic"), ("utr3", "3UTR"),
                          ("utr5", "5UTR"), ("prom", "Promoter-TSS")):
            out[overlaps_any(a, ivs[key])] = name
        return out
    t = timeit(label, repeat=5)
    rec.add(part="annot_breakdown", step="label 100k CREs (5 vectorised overlap tests)", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: (setattr(CA, "_uids", None), list(CA.uids)), repeat=5)
    rec.add(part="annot_breakdown", step="uid strings for the output table", seconds=t["median"], runs=t["runs"])


if __name__ == "__main__":
    parts = sys.argv[1:] or ["memory", "access", "functions"]
    rec = Recorder("columnar")
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save()
