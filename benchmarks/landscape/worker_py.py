#!/usr/bin/env python3
"""The Python competitors that live in their own environment, one task per call.

    python worker_py.py <task> <outdir> <args...>

Same contract as worker.R: one JSON object per line per (engine, op) with the
timed runs and the normalised answer file; inputs are read outside the clock.
Run it with the interpreter that has the competitors (``GB_LANDSCAPE_PY``), so
nothing here shadows the packages genomeblocks' own backends use.
"""
from __future__ import annotations

import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPEAT = int(os.environ.get("GB_REPEAT", "3"))
task, outdir, *rest = sys.argv[1:]
OUT = Path(outdir)
OUT.mkdir(parents=True, exist_ok=True)


def emit(engine, op, secs, n_out=None, out=None, **kw):
    print(json.dumps({"engine": engine, "op": op, "seconds": secs, "n_out": n_out, "out": out, **kw}), flush=True)


def run(engine, op, fn, times=REPEAT):
    try:
        res = fn()
        secs = []
        for _ in range(times):
            gc.collect()
            t0 = time.perf_counter()
            res = fn()
            secs.append(time.perf_counter() - t0)
        emit(engine, op, secs, n_out=res.get("n"), out=res.get("file"))
    except Exception as e:                                  # one engine failing must not kill the task
        emit(engine, op, None, error=f"{type(e).__name__}: {e}")


def read_bed(p):
    return pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])


def save(df, name):
    p = OUT / name
    df.to_csv(p, sep="\t", header=False, index=False)
    return {"n": len(df), "file": str(p)}


def _rnanorm_tmm(rnanorm, counts):
    t = rnanorm.TMM().fit(counts.T)                         # samples x genes
    f = np.asarray(t.get_norm_factors(counts.T), float).ravel()
    return save(pd.DataFrame({"sample": counts.columns, "factor": f}), "rnanorm_tmm.tsv")


def _rnanorm_cpm(rnanorm, counts):
    c = rnanorm.CPM().fit_transform(counts.T)
    df = pd.DataFrame(np.asarray(c).T, index=counts.index, columns=counts.columns)
    p = OUT / "rnanorm_cpm.tsv"
    df.to_csv(p, sep="\t", header=False)
    return {"n": len(df), "file": str(p)}


# ── intervals ────────────────────────────────────────────────────────────────
if task == "intervals":
    # ops: intersect, merge, nearest, and the jobs added 2026-10-23 (overlap_count, coverage, genome_depth,
    # subtract, complement); rest[3] is the chrom.sizes file complement needs
    A, B = read_bed(rest[0]), read_bed(rest[1])
    ops = rest[2].split(",")
    sizes = pd.read_csv(rest[3], sep="\t", header=None, names=["chrom", "size"]) if len(rest) > 3 else None
    A["idx"] = np.arange(len(A))

    # pyranges 1.x (the "pyranges1" package), else whatever pyranges the environment has
    try:
        try:
            import pyranges1 as pr
        except ImportError:
            import pyranges as pr
        v = getattr(pr, "__version__", "?")
        name = f"pyranges {v}"
        pa = pr.PyRanges(A.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
        pb = pr.PyRanges(B.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
        if "intersect" in ops:
            def _i():
                r = pa.overlap(pb)
                d = r.df if hasattr(r, "df") else pd.DataFrame(r)
                return save(d[["Chromosome", "Start", "End"]], "pyranges_intersect.tsv")
            run(name, "intersect", _i)
        if "merge" in ops:
            def _m():
                r = pa.merge_overlaps() if hasattr(pa, "merge_overlaps") else pa.merge()
                d = r.df if hasattr(r, "df") else pd.DataFrame(r)
                return save(d[["Chromosome", "Start", "End"]], "pyranges_merge.tsv")
            run(name, "merge", _m)
        if "nearest" in ops:
            def _n():
                r = pa.nearest_ranges(pb) if hasattr(pa, "nearest_ranges") else pa.nearest(pb)
                d = r.df if hasattr(r, "df") else pd.DataFrame(r)
                bs = "Start_b" if "Start_b" in d.columns else [c for c in d.columns if c.startswith("Start") and c != "Start"][0]
                be = "End_b" if "End_b" in d.columns else [c for c in d.columns if c.startswith("End") and c != "End"][0]
                return save(pd.DataFrame({"idx": d["idx"], "bstart": d[bs], "bend": d[be]}), "pyranges_nearest.tsv")
            run(name, "nearest", _n)
        # pyranges 1.x has no coverage of A by B; coverage_ranges is the genome-wide depth track (it needs the
        # runlength package)
        tab = lambda r: r.df if hasattr(r, "df") else pd.DataFrame(r)
        if "overlap_count" in ops:
            run(name, "overlap_count", lambda: save(tab(pa.count_overlaps(pb))[["Chromosome", "Start", "End", "Count"]], "pyranges_count.tsv"))
        if "genome_depth" in ops:
            run(name, "genome_depth", lambda: save(tab(pa.coverage_ranges())[["Chromosome", "Start", "End", "Coverage"]], "pyranges_depth.tsv"))
        if "subtract" in ops:
            run(name, "subtract", lambda: save(tab(pa.subtract_overlaps(pb))[["Chromosome", "Start", "End"]], "pyranges_subtract.tsv"))
        if "complement" in ops and sizes is not None:
            cs = dict(zip(sizes["chrom"], sizes["size"].astype(int)))
            # include_first_interval: the gap before a chromosome's first interval is part of the complement
            run(name, "complement", lambda: save(tab(pa.complement_ranges(include_first_interval=True, chromsizes=cs))[["Chromosome", "Start", "End"]],
                                                 "pyranges_complement.tsv"))
    except ImportError:
        pass

    # polars-bio
    try:
        import polars as pl
        import polars_bio as pb_
        # polars-bio reads intervals as 1-based closed unless told otherwise; the BED files here are 0-based
        # half-open (coverage of [100, 200) by [180, 260) is 21 bases in the default mode, 20 in this one)
        pb_.set_option(pb_.POLARS_BIO_COORDINATE_SYSTEM_ZERO_BASED, True)
        v = getattr(pb_, "__version__", "?")
        name = f"polars-bio {v}"
        la = pl.from_pandas(A)
        lb = pl.from_pandas(B[["chrom", "start", "end"]])
        if "intersect" in ops:
            def _i():
                r = pb_.overlap(la, lb, cols1=("chrom", "start", "end"), cols2=("chrom", "start", "end"), output_type="polars.DataFrame")
                c = [c for c in r.columns if c.startswith("chrom")][0]
                s = [c for c in r.columns if c.startswith("start")][0]
                e = [c for c in r.columns if c.startswith("end")][0]
                d = r.select([c, s, e]).unique().to_pandas()
                return save(d, "polarsbio_intersect.tsv")
            run(name, "intersect", _i)
        if "merge" in ops:
            def _m():
                r = pb_.merge(la.select(["chrom", "start", "end"]), cols=("chrom", "start", "end"), output_type="polars.DataFrame")
                return save(r.select(["chrom", "start", "end"]).to_pandas(), "polarsbio_merge.tsv")
            run(name, "merge", _m)
        if "nearest" in ops:
            def _n():
                r = pb_.nearest(la, lb, cols1=("chrom", "start", "end"), cols2=("chrom", "start", "end"), output_type="polars.DataFrame")
                idx = "idx_1" if "idx_1" in r.columns else "idx"
                d = r.select([idx, "start_2", "end_2"]).unique(subset=[idx]).to_pandas()
                d.columns = ["idx", "bstart", "bend"]
                return save(d, "polarsbio_nearest.tsv")
            run(name, "nearest", _n)
        # the jobs added 2026-10-23 (polars-bio's depth() is for alignments, so genome_depth is not run)
        c3 = ("chrom", "start", "end")
        la3, lb3 = la.select(list(c3)), lb

        def _first(r, *keys):
            return [next(c for c in r.columns if c.startswith(k)) for k in keys]
        if "overlap_count" in ops:
            def _c():
                r = pb_.count_overlaps(la3, lb3, cols1=list(c3), cols2=list(c3), output_type="polars.DataFrame")
                return save(r.select(_first(r, "chrom", "start", "end", "count")).to_pandas(), "polarsbio_count.tsv")
            run(name, "overlap_count", _c)
        if "coverage" in ops:
            def _v():
                r = pb_.coverage(la3, lb3, cols1=list(c3), cols2=list(c3), output_type="polars.DataFrame")
                return save(r.select(_first(r, "chrom", "start", "end", "coverage")).to_pandas(), "polarsbio_coverage.tsv")
            run(name, "coverage", _v)
        if "subtract" in ops:
            def _s():
                r = pb_.subtract(la3, lb3, cols1=list(c3), cols2=list(c3), output_type="polars.DataFrame")
                return save(r.select(_first(r, "chrom", "start", "end")).to_pandas(), "polarsbio_subtract.tsv")
            run(name, "subtract", _s)
        if "complement" in ops and sizes is not None:
            view = pl.DataFrame({"chrom": sizes["chrom"].tolist(), "start": [0] * len(sizes), "end": sizes["size"].astype(int).tolist()})

            def _x():
                r = pb_.complement(la3, view_df=view, output_type="polars.DataFrame")
                return save(r.select(_first(r, "chrom", "start", "end")).to_pandas(), "polarsbio_complement.tsv")
            run(name, "complement", _x)
    except ImportError:
        pass

    # BiocPy genomicranges
    try:
        import genomicranges as grp
        from iranges import IRanges
        v = getattr(grp, "__version__", "?")
        name = f"BiocPy genomicranges {v}"
        ga = grp.GenomicRanges(seqnames=A["chrom"].tolist(), ranges=IRanges(A["start"].tolist(), (A["end"] - A["start"]).tolist()))
        gb_ = grp.GenomicRanges(seqnames=B["chrom"].tolist(), ranges=IRanges(B["start"].tolist(), (B["end"] - B["start"]).tolist()))
        if "intersect" in ops:
            def _i():
                r = ga.subset_by_overlaps(gb_)
                d = pd.DataFrame({"chrom": r.get_seqnames(), "start": r.get_start(), "end": np.asarray(r.get_end()) + 1})   # inclusive end -> BED
                return save(d, "biocpy_intersect.tsv")
            run(name, "intersect", _i)
        if "merge" in ops:
            def _m():
                r = ga.reduce()
                d = pd.DataFrame({"chrom": r.get_seqnames(), "start": r.get_start(), "end": np.asarray(r.get_end()) + 1})
                return save(d, "biocpy_merge.tsv")
            run(name, "merge", _m)
        if "nearest" in ops:
            def _n():
                # BiocPy's nearest(query) searches *self* for each range of the argument: subject.nearest(query)
                h = np.asarray(gb_.nearest(ga, select="arbitrary"), int)      # one B index per A row, -1 = none
                st, en = np.asarray(gb_.get_start()), np.asarray(gb_.get_end()) + 1
                ok = h >= 0
                idx = np.flatnonzero(ok)
                return save(pd.DataFrame({"idx": idx, "bstart": st[h[ok]], "bend": en[h[ok]]}), "biocpy_nearest.tsv")
            run(name, "nearest", _n)
        # the jobs added 2026-10-23; BiocPy has no coverage of A by B
        if "overlap_count" in ops:
            def _c():
                cnt = np.asarray(gb_.count_overlaps(ga), np.int64)       # subject.count_overlaps(query): one count per A row
                return save(pd.DataFrame({"chrom": A["chrom"], "start": A["start"], "end": A["end"], "count": cnt}), "biocpy_count.tsv")
            run(name, "overlap_count", _c)
        if "subtract" in ops and len(A) > 100_000:
            # 13 s at 10k, 181 s at 100k, 3,000 s per run at 1M (and Jaccard 0.24 against bedtools) on 2026-10-23
            emit(name, "subtract", None, error=f"not run at n = {len(A):,}: 3,000 s per run at 1 M on 2026-10-23 "
                                                f"(superlinear; Jaccard 0.24); timed at 10k and 100k")
        elif "subtract" in ops:
            def _s():
                r = ga.subtract(gb_).unlist()                           # one range list per A row, flattened
                return save(pd.DataFrame({"chrom": r.get_seqnames(), "start": r.get_start(), "end": np.asarray(r.get_end()) + 1}), "biocpy_subtract.tsv")
            run(name, "subtract", _s)
        if "complement" in ops and sizes is not None:
            # A's 0-based starts went in as BiocPy starts, so a chromosome [0, L) is the closed range 0 .. L-1
            last = {c: int(s) - 1 for c, s in zip(sizes["chrom"], sizes["size"])}

            def _g():
                r = ga.gaps(start=0, end=last)
                keep = np.asarray(r.get_strand()) == 0                   # gaps() also returns whole chromosomes on + and -
                return save(pd.DataFrame({"chrom": np.asarray(r.get_seqnames())[keep], "start": np.asarray(r.get_start())[keep],
                                          "end": (np.asarray(r.get_end()) + 1)[keep]}), "biocpy_complement.tsv")
            run(name, "complement", _g)
        if "genome_depth" in ops:
            emit(name, "genome_depth", None, error="not run: GenomicRanges.coverage() returns a dense per-base float64 array per "
                                                     "chromosome (about 25 GB for hg38)")
    except ImportError:
        pass

    # superintervals / ailist: indexed overlap only (A rows with any hit in B)
    for mod, label in (("superintervals", "superintervals"), ("ailist", "AIList (ailist)")):
        try:
            m = __import__(mod)
        except ImportError:
            continue
        if "intersect" not in ops:
            continue
        from importlib.metadata import version as _ver
        try:
            v = _ver(mod)
        except Exception:
            v = getattr(m, "__version__", "?")

        def _si(m=m, mod=mod):
            keep = np.zeros(len(A), bool)
            for chrom, bsub in B.groupby("chrom"):
                asub = A[A["chrom"] == chrom]
                if mod == "superintervals":
                    idx = m.IntervalMap()
                    for s, e in zip(bsub["start"].to_numpy(), bsub["end"].to_numpy() - 1):      # closed ends
                        idx.add(int(s), int(e), 0)
                    idx.build()
                    qs = np.array(asub["start"].to_numpy(), dtype=np.int32)       # count_batch wants C ints
                    qe = np.array(asub["end"].to_numpy() - 1, dtype=np.int32)
                    hit = np.asarray(idx.count_batch(qs, qe)) > 0
                else:
                    al = m.AIList()
                    for s, e in zip(bsub["start"].to_numpy(), bsub["end"].to_numpy()):
                        al.add(int(s), int(e))
                    al.construct()
                    hit = np.array([len(al.intersect(int(s), int(e))) > 0 for s, e in zip(asub["start"], asub["end"])], bool)
                keep[asub.index[hit]] = True
            return save(A.loc[keep, ["chrom", "start", "end"]], f"{mod}_intersect.tsv")
        run(f"{label} {v}", "intersect", _si, 1)

# ── expression: Python ports of tximport / edgeR / DESeq2 normalisation ─────
if task == "expression":
    counts_file, quant_dir, tx2gene_file = rest[0], rest[1], rest[2]
    counts = pd.read_csv(counts_file, sep="\t", index_col=0)
    X = counts.to_numpy(np.float64)
    try:
        import rnanorm
        from importlib.metadata import version as _ver
        v = _ver("rnanorm")
        run(f"rnanorm {v}", "tmm", lambda: _rnanorm_tmm(rnanorm, counts))
        run(f"rnanorm {v}", "cpm", lambda: _rnanorm_cpm(rnanorm, counts))
    except ImportError:
        pass
    try:
        import pydeseq2
        from pydeseq2.preprocessing import deseq2_norm
        v = getattr(pydeseq2, "__version__", "?")

        def _pyd():
            _, sf = deseq2_norm(counts.T)                   # samples x genes
            return save(pd.DataFrame({"sample": counts.columns, "factor": np.asarray(sf, float).ravel()}), "pydeseq2_sf.tsv")
        run(f"pydeseq2 {v}", "size_factors", _pyd)
    except ImportError:
        pass
    try:
        import pytximport
        v = getattr(pytximport, "__version__", "?")
        files = sorted(Path(quant_dir).glob("*/quant.sf"))
        t2g = pd.read_csv(tx2gene_file, sep="\t", header=None, names=["transcript_id", "gene_id"])

        def _ptx():
            r = pytximport.tximport([str(f) for f in files], data_type="salmon", transcript_gene_map=t2g,
                                    counts_from_abundance="length_scaled_tpm", output_type="xarray", ignore_transcript_version=True)
            c = r["counts"]
            genes_dim = next(d for d in c.dims if d != "file_path" and d != "file")
            df = pd.DataFrame(np.asarray(c.transpose(genes_dim, ...)), index=np.asarray(c.coords[genes_dim]).astype(str))
            p = OUT / "pytximport_genes.tsv"
            df.to_csv(p, sep="\t", header=False)
            return {"n": len(df), "file": str(p)}
        run(f"pytximport {v}", "tx_to_gene", _ptx, 1)
    except ImportError:
        pass



def save_matrix(m, name):
    """A gene x sample matrix (index = genes) without header."""
    p = OUT / name
    m.to_csv(p, sep="\t", header=False)
    return {"n": len(m), "file": str(p)}


# ── expression_de: the DE family on simulated counts (args: <counts.tsv> <coldata.tsv>; ~ batch + group) ──
if task == "expression_de":
    import copy
    cts = pd.read_csv(rest[0], sep="\t", index_col=0)
    cd = pd.read_csv(rest[1], sep="\t", index_col=0, dtype=str)
    X = cts.T                                                   # samples x genes, the layout both ports take
    try:
        import pydeseq2
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        name = f"PyDESeq2 {pydeseq2.__version__}"

        def _dds():
            return DeseqDataSet(counts=X, metadata=cd, design="~batch + group", quiet=True, n_cpus=1)
        fitted = []

        def _de():
            dds = _dds()
            dds.deseq2()
            s = DeseqStats(dds, contrast=["group", "B", "A"], quiet=True)
            s.summary()
            fitted.append(s)
            r = s.results_df
            return save(pd.DataFrame({"gene": r.index, "p": r["pvalue"].to_numpy(), "lfc": r["log2FoldChange"].to_numpy()}), "pydeseq2_de.tsv")
        run(name, "de_test", _de, 1)

        def _disp():
            dds = _dds()
            dds.fit_size_factors()
            dds.fit_genewise_dispersions()
            dds.fit_dispersion_trend()
            dds.fit_dispersion_prior()
            dds.fit_MAP_dispersions()
            return save(pd.DataFrame({"gene": dds.var_names, "disp": dds.var["dispersions"].to_numpy()}), "pydeseq2_disp.tsv")
        run(name, "dispersion", _disp, 1)
        if fitted:
            # lfc_shrink rewrites the results in place: one fresh copy per call, made outside the clock
            copies = [copy.deepcopy(fitted[-1]) for _ in range(2)]

            def _shr():
                s = copies.pop()
                s.lfc_shrink(coeff="group[T.B]")
                r = s.results_df
                return save(pd.DataFrame({"gene": r.index, "lfc": r["log2FoldChange"].to_numpy()}), "pydeseq2_shrink.tsv")
            run(name, "lfc_shrinkage", _shr, 1)

        def _vst():
            dds = _dds()
            dds.vst(use_design=False)
            return save_matrix(pd.DataFrame(np.asarray(dds.layers["vst_counts"]).T, index=dds.var_names, columns=dds.obs_names), "pydeseq2_vst.tsv")
        run(name, "vst", _vst, 1)
    except ImportError:
        pass
    try:
        import inmoose
        from inmoose import deseq2 as ide
        name = f"InMoose {inmoose.__version__} deseq2"

        def _ide():
            d = ide.DESeq(ide.DESeqDataSet(X, clinicalData=cd, design="~ batch + group"), quiet=True)
            r = d.results(contrast=["group", "B", "A"])
            return save(pd.DataFrame({"gene": r.index, "p": r["pvalue"].to_numpy(), "lfc": r["log2FoldChange"].to_numpy()}), "inmoose_deseq2_de.tsv")
        run(name, "de_test", _ide, 1)

        def _idisp():
            d = ide.DESeqDataSet(X, clinicalData=cd, design="~ batch + group")
            d = d.estimateSizeFactors()
            d = d.estimateDispersions(quiet=True)
            return save(pd.DataFrame({"gene": d.var_names, "disp": d.var["dispersion"].to_numpy()}), "inmoose_deseq2_disp.tsv")
        run(name, "dispersion", _idisp, 1)

        def _ivst():
            v = ide.varianceStabilizingTransformation(ide.DESeqDataSet(X, clinicalData=cd, design="~ batch + group"), blind=True)
            return save_matrix(pd.DataFrame(np.asarray(v.X).T, index=v.var_names, columns=v.obs_names), "inmoose_deseq2_vst.tsv")
        run(name, "vst", _ivst, 1)
    except ImportError:
        pass
    try:
        import inmoose
        from inmoose import edgepy as ie
        from patsy import dmatrix
        name = f"InMoose {inmoose.__version__} edgepy"
        dm = dmatrix("~ batch + group", cd)                     # edgepy reads the coefficient names from a patsy matrix

        def _ql():
            y = ie.DGEList(counts=cts.to_numpy(float))
            y.estimateGLMCommonDisp(design=dm)
            y.estimateGLMTagwiseDisp(design=dm)
            t = ie.glmQLFTest(y.glmQLFit(design=dm), coef=2)          # a DataFrame-like DGELRT (DESeq2-style column names)
            return save(pd.DataFrame({"gene": cts.index, "p": np.asarray(t["pvalue"]), "lfc": np.asarray(t["log2FoldChange"])}), "inmoose_edgepy_de.tsv")
        run(name, "de_test", _ql, 1)
    except ImportError:
        pass
    try:
        import inmoose
        from inmoose.pycombat import pycombat_seq
        run(f"InMoose {inmoose.__version__} pycombat_seq", "batch_correction",
            lambda: save_matrix(pd.DataFrame(np.asarray(pycombat_seq(cts, cd["batch"].tolist(), covar_mod=cd[["group"]])), index=cts.index),
                                "inmoose_combatseq_counts.tsv"), 1)
    except ImportError:
        pass


# ── track_tmm: TMM-normalised CPM of a region x sample count table (args: <counts.tsv>) ──
if task == "track_tmm":
    m = pd.read_csv(rest[0], sep="\t")
    try:
        import rnanorm
        from importlib.metadata import version as _ver

        def _t():
            x = rnanorm.TMM().fit_transform(m.T)                 # samples x regions in, the same out
            p = OUT / "rnanorm_tmm_cpm.tsv"
            pd.DataFrame(np.asarray(x).T).to_csv(p, sep="\t", header=False, index=False)
            return {"n": len(m), "file": str(p)}
        run(f"rnanorm {_ver('rnanorm')} TMM", "track_tmm", _t)
    except ImportError:
        pass


# ── genes_window: (peak, gene) pairs with the TSS within w (args: <peaks.bed> <tss.bed> <w>) ──
if task == "genes_window":
    P = read_bed(rest[0])
    T = pd.read_csv(rest[1], sep="\t", header=None, names=["chrom", "start", "end", "gene"])
    w = int(rest[2])
    try:
        import pyranges1 as pr
        pa = pr.PyRanges(P.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
        pt = pr.PyRanges(T.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))

        def _j():
            r = pa.join_overlaps(pt, slack=w)
            d = r.df if hasattr(r, "df") else pd.DataFrame(r)
            return save(d[["Chromosome", "Start", "End", "gene"]], "pyranges1_window.tsv")
        run(f"pyranges {pr.__version__} join_overlaps (slack)", "genes_in_window", _j)
    except ImportError:
        pass


# ── chromvar: per-sample motif deviations with pychromVAR (args: <peaks.fa> <counts.tsv> <motifs.jaspar> <groups>) ──
if task == "chromvar":
    try:
        import anndata as ad
        import pychromvar as pc
        from Bio import motifs as bmotifs
        from importlib.metadata import version as _ver
        cnt = pd.read_csv(rest[1], sep="\t", index_col=0)
        width = len(open(rest[0]).read().split("\n")[1])
        mot = list(bmotifs.parse(open(rest[2]), "jaspar"))
        # pychromVAR 0.0.4 still passes AnnData(dtype=...), which anndata >= 0.11 removed: drop the argument
        _init = ad.AnnData.__init__
        if "dtype" not in _init.__code__.co_varnames:
            def _compat(self, *x, dtype=None, **k):
                _init(self, *x, **k)
            ad.AnnData.__init__ = _compat

        def _cv():
            a = ad.AnnData(X=cnt.T.to_numpy(float), obs=pd.DataFrame(index=cnt.columns), var=pd.DataFrame(index=[f"{p}-0-{width}" for p in cnt.index]))
            pc.add_peak_seq(a, genome_file=rest[0], delimiter="-")
            pc.add_gc_bias(a)
            pc.get_bg_peaks(a, niterations=50, n_jobs=1)
            pc.match_motif(a, motifs=mot, p_value=5e-5, background="even")
            dev = pc.compute_deviations(a, n_jobs=1)
            names = [m.matrix_id for m in mot] if len(mot) == dev.shape[1] else list(dev.var_names)
            p = OUT / "pychromvar_z.tsv"
            pd.DataFrame(np.asarray(dev.X).T, index=names).to_csv(p, sep="\t", header=False)
            return {"n": dev.shape[1], "file": str(p)}
        run(f"pychromVAR {_ver('pychromvar')}", "motif_activity", _cv, 1)
    except ImportError:
        pass
