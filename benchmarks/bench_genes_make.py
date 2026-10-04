#!/usr/bin/env python3
"""Where Genes.make spends its time, and what would make it faster.

GTF: 20k genes, 70k transcripts, 877k lines (synthetic, GENCODE-shaped).

  components  read lines / split / parse attributes / full make
  parse       current make vs. a faster pure-Python parser (same objects)
              vs. columnar parsers (pandas, pyranges.read_gtf, polars)
  annot       the lazy promoter/exon/UTR index: current vs. key-sort fix
              vs. numpy merge (from objects, and from a columnar table)
"""
from __future__ import annotations

import re
import sys

import numpy as np

from common import DATA, Recorder, timeit
from prototypes.npintervals import Intervals, merge, slop

from genomeblocks import Genes
from genomeblocks.genes import Gene, Transcript, _parse_attributes
from genomeblocks.locus import CDS, UTR, Exon

GTF = str(DATA / "genes.gtf")


# ── a faster pure-Python parser that builds the very same objects ───────────

_ATTR = re.compile(r'(\S+) "?([^";]*)"?;?')


def make_fast(filename, promoter_r=1000):
    """Genes.make's logic with a regex attribute parser and no tqdm."""
    genes = Genes(filename=filename, _promoter_r=promoter_r)
    findall = _ATTR.findall
    with open(filename) as f:
        for line in f:
            if line[0] == "#" or line == "\n":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9:
                chrom, _, ft, start, end, _, strand, _, attributes = fields
            elif len(fields) == 8:
                chrom, _, ft, start, end, _, strand, attributes = fields
            else:
                continue
            start, end = int(start), int(end)
            attrs = dict(findall(attributes))
            gid = attrs["gene_id"]
            if ft == "exon":
                genes[gid].transcripts[attrs["transcript_id"]].exons.append(
                    Exon(chrom, start, end, strand, exon_number=int(attrs.get("exon_number", 0))))
            elif ft == "CDS":
                genes[gid].transcripts[attrs["transcript_id"]].cds.append(
                    CDS(chrom, start, end, strand, exon_number=int(attrs.get("exon_number", 0))))
            elif ft in ("five_prime_UTR", "three_prime_UTR"):
                genes[gid].transcripts[attrs["transcript_id"]].utr.append(
                    UTR(chrom, start, end, strand, exon_number=int(attrs.get("exon_number", 0)),
                        type="5'" if ft == "five_prime_UTR" else "3'"))
            elif ft == "transcript":
                tid = attrs["transcript_id"]
                if gid not in genes:
                    genes[gid] = Gene(chrom, start, end, strand, gene_id=gid)
                genes[gid].transcripts[tid] = Transcript(chrom, start, end, strand, tid)
            elif ft == "gene":
                genes[gid] = Gene(chrom, start, end, strand, gene_id=gid,
                                  gene_name=attrs.get("gene_name", gid),
                                  gene_type=attrs.get("gene_type"))
    return genes


def same_genes(a, b) -> bool:
    if a.keys() != b.keys():
        return False
    for gid, g in a.items():
        h = b[gid]
        if (g.chrom, g.start, g.end, g.strand, g.gene_name, g.gene_type) != \
           (h.chrom, h.start, h.end, h.strand, h.gene_name, h.gene_type):
            return False
        if g.transcripts.keys() != h.transcripts.keys():
            return False
        for tid, t in g.transcripts.items():
            u = h.transcripts[tid]
            for x, y in ((t.exons, u.exons), (t.cds, u.cds), (t.utr, u.utr)):
                if [(e.start, e.end, e.exon_number) for e in x] != [(e.start, e.end, e.exon_number) for e in y]:
                    return False
    return True


# ── columnar parsers (one table, no per-feature objects) ────────────────────

KEYS = ("gene_id", "transcript_id", "gene_name", "gene_type", "exon_number")


def parse_pandas(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attr"],
                     usecols=["chrom", "feature", "start", "end", "strand", "attr"],
                     dtype={"chrom": "category", "feature": "category", "strand": "category"})
    for k in KEYS:
        df[k] = df["attr"].str.extract(f'{k} "?([^";]*)"?;', expand=False)
    return df.drop(columns="attr")


def parse_polars(path):
    import polars as pl
    df = pl.read_csv(path, separator="\t", comment_prefix="#", has_header=False,
                     columns=[0, 2, 3, 4, 6, 8], quote_char=None,
                     new_columns=["chrom", "feature", "start", "end", "strand", "attr"])
    return df.with_columns([pl.col("attr").str.extract(f'{k} "?([^";]*)"?;', 1).alias(k) for k in KEYS]).drop("attr")


def objects_from_table(df, filename=GTF, promoter_r=1000):
    """Build Genes.make's exact object tree from the polars table."""
    import polars as pl
    genes = Genes(filename=filename, _promoter_r=promoter_r)
    cols = ["chrom", "start", "end", "strand", "gene_id", "transcript_id", "gene_name", "gene_type", "exon_number"]
    for ft in ("gene", "transcript", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
        sub = df.filter(pl.col("feature") == ft).select(cols)
        rows = zip(*(sub[c].to_list() for c in cols))
        if ft == "gene":
            for c, s, e, st, gid, _, gn, gt, _ in rows:
                genes[gid] = Gene(c, s, e, st, gene_id=gid, gene_name=gn if gn is not None else gid, gene_type=gt)
        elif ft == "transcript":
            for c, s, e, st, gid, tid, *_ in rows:
                if gid not in genes:
                    genes[gid] = Gene(c, s, e, st, gene_id=gid)
                genes[gid].transcripts[tid] = Transcript(c, s, e, st, tid)
        else:
            for c, s, e, st, gid, tid, _, _, en in rows:
                n = int(en) if en is not None else 0
                t = genes[gid].transcripts[tid]
                if ft == "exon":
                    t.exons.append(Exon(c, s, e, st, exon_number=n))
                elif ft == "CDS":
                    t.cds.append(CDS(c, s, e, st, exon_number=n))
                else:
                    t.utr.append(UTR(c, s, e, st, exon_number=n, type="5'" if ft == "five_prime_UTR" else "3'"))
    return genes


def parse_pyranges(path):
    import pyranges as pr
    return pr.read_gtf(path)


# ── annot builds ────────────────────────────────────────────────────────────

def annot_current(g):
    g._annot = None
    return g.annot


def annot_keysort(g):
    """Same Loci pipeline, with the key-based sort and no eager uid dict."""
    from genomeblocks.loci import Loci
    from genomeblocks.locus import Locus

    def sm(L):
        s = sorted(L, key=lambda l: (l.chrom, l.start))
        out = []
        if s:
            c, a, b = s[0].chrom, s[0].start, s[0].end
            for l in s[1:]:
                if l.chrom == c and l.start <= b:
                    b = max(b, l.end)
                else:
                    out.append(Locus(c, a, b)); c, a, b = l.chrom, l.start, l.end
            out.append(Locus(c, a, b))
        return Loci(out)
    T = [t for x in g.values() for t in x.transcripts.values()]
    return {"body": Loci(g.values()),
            "prom": sm(Loci(g.get_tss().values()).slop(g._promoter_r)),
            "exon": sm([e for t in T for e in t.exons]),
            "utr5": sm([u for t in T for u in t.utr if u.type == "5'"]),
            "utr3": sm([u for t in T for u in t.utr if u.type == "3'"])}


def annot_numpy_objects(g):
    T = [t for x in g.values() for t in x.transcripts.values()]
    tss = Intervals.from_loci(list(g.get_tss().values()))      # '-' TSS is Locus(end, end-1), as in genomeblocks
    return {"prom": merge(slop(tss, g._promoter_r)),
            "exon": merge(Intervals.from_loci([e for t in T for e in t.exons])),
            "utr5": merge(Intervals.from_loci([u for t in T for u in t.utr if u.type == "5'"])),
            "utr3": merge(Intervals.from_loci([u for t in T for u in t.utr if u.type == "3'"]))}


def annot_numpy_table(df, promoter_r=1000):
    """From the polars feature table directly: no objects at all."""
    import polars as pl

    def iv(sub):
        return Intervals.from_frame(sub.select(["chrom", "start", "end"]).to_pandas())
    genes = df.filter(pl.col("feature") == "gene")
    tss = genes.with_columns(
        pl.when(pl.col("strand") == "+").then(pl.col("start")).otherwise(pl.col("end")).alias("start"),
        pl.when(pl.col("strand") == "+").then(pl.col("start") + 1).otherwise(pl.col("end") - 1).alias("end"))
    return {"prom": merge(slop(iv(tss), promoter_r)),
            "exon": merge(iv(df.filter(pl.col("feature") == "exon"))),
            "utr5": merge(iv(df.filter(pl.col("feature") == "five_prime_UTR"))),
            "utr3": merge(iv(df.filter(pl.col("feature") == "three_prime_UTR")))}


def as_tuples(x):
    if isinstance(x, Intervals):
        return sorted(zip([x.names[c] for c in x.codes], x.starts.tolist(), x.ends.tolist()))
    return sorted((l.chrom, l.start, l.end) for l in x)


def part_components(rec):
    print("\n== components ==")
    t = timeit(lambda: sum(1 for _ in open(GTF)), repeat=3)
    rec.add(part="components", step="read lines", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: [l.strip().split("\t") for l in open(GTF) if l[0] != "#"], repeat=3)
    rec.add(part="components", step="read + split", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: [_parse_attributes(l.strip().split("\t")[8]) for l in open(GTF) if l[0] != "#"], repeat=3)
    rec.add(part="components", step="read + split + _parse_attributes", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: [dict(_ATTR.findall(l.rstrip("\n").split("\t")[8])) for l in open(GTF) if l[0] != "#"], repeat=3)
    rec.add(part="components", step="read + split + regex attributes", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: Genes.make(GTF), repeat=3)
    rec.add(part="components", step="Genes.make (all)", seconds=t["median"], runs=t["runs"])


def part_parse(rec):
    print("\n== parse ==")
    res = {}
    t = timeit(lambda: res.__setitem__("a", Genes.make(GTF)), repeat=3)
    rec.add(part="parse", variant="Genes.make (current)", objects=True, seconds=t["median"], runs=t["runs"])
    import gc

    def no_gc(fn):
        def run():
            gc.disable()
            try:
                return fn()
            finally:
                gc.enable()
        return run
    t = timeit(no_gc(lambda: res.__setitem__("g", Genes.make(GTF))), repeat=3)
    rec.add(part="parse", variant="Genes.make, garbage collector paused", objects=True,
            seconds=t["median"], runs=t["runs"], same=same_genes(res["a"], res["g"]))
    t = timeit(lambda: res.__setitem__("b", make_fast(GTF)), repeat=3)
    rec.add(part="parse", variant="fast Python parser, same objects", objects=True, seconds=t["median"],
            runs=t["runs"], same=same_genes(res["a"], res["b"]))
    for name, fn in (("pandas read_csv + str.extract (table)", parse_pandas),
                     ("pyranges.read_gtf (table)", parse_pyranges),
                     ("polars read_csv + str.extract (table)", parse_polars)):
        t = timeit(lambda: res.__setitem__("c", fn(GTF)), repeat=3)
        rec.add(part="parse", variant=name, objects=False, seconds=t["median"], runs=t["runs"],
                rows=int(len(res["c"])))
    t = timeit(lambda: res.__setitem__("d", objects_from_table(parse_polars(GTF))), repeat=3)
    rec.add(part="parse", variant="polars table, then the same objects", objects=True,
            seconds=t["median"], runs=t["runs"], same=same_genes(res["a"], res["d"]))


def part_annot(rec):
    print("\n== annot ==")
    g = Genes.make(GTF)
    res = {}
    t = timeit(lambda: res.__setitem__("a", annot_current(g)), repeat=3)
    rec.add(part="annot", variant="current (Loci.sort/merge)", seconds=t["median"], runs=t["runs"])
    ref = {k: as_tuples(v) for k, v in res["a"].items() if k != "body"}
    t = timeit(lambda: res.__setitem__("b", annot_keysort(g)), repeat=3)
    rec.add(part="annot", variant="key-sort fix (same Loci objects)", seconds=t["median"], runs=t["runs"],
            same=all(as_tuples(res["b"][k]) == ref[k] for k in ref))
    t = timeit(lambda: res.__setitem__("c", annot_numpy_objects(g)), repeat=3)
    rec.add(part="annot", variant="numpy merge (from Gene objects)", seconds=t["median"], runs=t["runs"],
            same=all(as_tuples(res["c"][k]) == ref[k] for k in ref))
    df = parse_polars(GTF)
    t = timeit(lambda: res.__setitem__("d", annot_numpy_table(df)), repeat=3)
    rec.add(part="annot", variant="numpy merge (from columnar table)", seconds=t["median"], runs=t["runs"],
            same=all(as_tuples(res["d"][k]) == ref[k] for k in ref))


if __name__ == "__main__":
    parts = sys.argv[1:] or ["components", "parse", "annot"]
    rec = Recorder("genes_make")
    for p in parts:
        globals()[f"part_{p}"](rec)
    rec.save(gtf=GTF)
