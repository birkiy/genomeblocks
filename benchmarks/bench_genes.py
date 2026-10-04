#!/usr/bin/env python3
"""Gene models: GTF parsing, region annotation, nearest gene, isoform choice.

GTF: 20k genes, ~70k transcripts, 877k lines (about 1/4 of GENCODE).

  make            Genes.make(gtf)  — pure-Python line parser into dataclasses
  annot           lazy promoter / exon / UTR / body index build (first access)
  annotations     label N CREs (Promoter-TSS / UTR / Exonic / Intronic / ...)
                  with cgranges vs the pure-Python fallback index
  nearest_genes   nearest TSS per CRE (pyranges)
  select_isoforms ATAC-supported isoform choice from peaks + a bigWig
"""
from __future__ import annotations

import gc
import time

from common import DATA, Recorder, timeit

import genomeblocks.loci as gl
from genomeblocks import Genes, Loci

GTF = str(DATA / "genes.gtf")


def build_annot(g):
    g._annot = None
    return g.annot


if __name__ == "__main__":
    import cgranges
    from bench_loci import _ORIG_PY
    rec = Recorder("genes")
    n_lines = sum(1 for _ in open(GTF))

    t = timeit(lambda: Genes.make(GTF), repeat=3)
    rec.add(step="Genes.make (GTF parse)", n=n_lines, seconds=t["median"], runs=t["runs"],
            rate=n_lines / t["median"])
    g = Genes.make(GTF)
    n_tx = sum(len(x.transcripts) for x in g.values())

    t = timeit(lambda: build_annot(g), repeat=3)
    rec.add(step="annot index build (prom/exon/UTR merge)", n=len(g), seconds=t["median"],
            runs=t["runs"])

    for n in (10_000, 100_000):
        L = Loci.make(str(DATA / f"peaks_A_{n}.bed"))
        for label, getcg, py in (("cgranges", lambda: cgranges, _ORIG_PY),
                                 ("pure-Python fallback", lambda: None, _ORIG_PY)):
            gl._get_cgranges, gl._PyIntervalIndex = getcg, py
            build_annot(g)
            t = timeit(lambda: g.annotations(L), repeat=3 if n <= 10_000 else 1)
            rec.add(step=f"annotations ({label})", n=n, seconds=t["median"], runs=t["runs"],
                    rate=n / t["median"])
        gl._get_cgranges, gl._PyIntervalIndex = (lambda: cgranges), _ORIG_PY
        build_annot(g)
        t = timeit(lambda: g.nearest_genes(L), repeat=3)
        rec.add(step="nearest_genes (pyranges)", n=n, seconds=t["median"], runs=t["runs"],
                rate=n / t["median"])

    # isoform selection: peaks only, and peaks + bigWig TSS scoring
    peaks = Loci.make(str(DATA / "peaks_A_100000.bed"))
    peaks.cgr
    bw = str(DATA / "signal_0.bw")
    for label, kw in (("select_isoforms (peaks)", dict(cre=peaks)),
                      ("select_isoforms (peaks + bigWig)", dict(cre=peaks, bw=bw)),
                      ("select_isoforms (bigWig only)", dict(bw=bw))):
        t = timeit(lambda: g.select_isoforms(verbose=False, **kw), repeat=3)
        rec.add(step=label, n=n_tx, seconds=t["median"], runs=t["runs"],
                rate=n_tx / t["median"])
    rec.save(gtf_lines=n_lines, genes=len(g), transcripts=n_tx)
