#!/usr/bin/env python3
"""Expression tables: reading quantifications, lining them up with the gene
tables, tximport-style summarising and normalising — the two table backends
against each other, and genomeblocks against R (tximport, edgeR, DESeq2) for
exactness and time.

Data (generated here, seeded): GENCODE-shaped GTFs from make_data.gtf (20k
genes / ~70k transcripts, and 60k genes / ~210k transcripts), salmon-style
quant.sf files for 6 / 24 (/ 96) samples with version-less transcript ids,
and a featureCounts matrix of the summed gene counts.

  read           Expression.make over N quant.sf files: polars vs pandas tables
                 backend; R tximport(txOut=TRUE) (readr) — same numbers?
  align          E.align(genes=): the id rule over every transcript
  to_genes       to_genes(how) x 3 vs R tximport::summarizeToGene — max |Δ|
  normalise      tmm / size_factors / log_cpm / filter_by_expr vs edgeR
                 calcNormFactors / DESeq2 estimateSizeFactorsForMatrix / cpm /
                 filterByExpr — max |Δ|, disagreements
  featurecounts  one gene x sample file: polars vs pandas vs R read.delim /
                 data.table::fread

R: GB_RSCRIPT="/path/Rscript:/other/Rscript" (default: Rscript on PATH); each
R step runs on the first Rscript that has its packages, and is skipped (and
listed) otherwise. GB_BENCH_EXPR="20000x6,60000x24" picks the sizes.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder, timeit

import genomeblocks as gb
from genomeblocks import Expression, Genes

EXPR = DATA / "expression"
EXPR.mkdir(parents=True, exist_ok=True)
SIZES = [tuple(int(v) for v in s.split("x"))
         for s in os.environ.get("GB_BENCH_EXPR", "20000x6,60000x24,60000x96").split(",")]
RSCRIPTS = [p for p in os.environ.get("GB_RSCRIPT", "Rscript").split(":") if p]

# ── R side ──────────────────────────────────────────────────────────────────

R_COMMON = r'''
args <- commandArgs(TRUE)
tm <- function(f, n = 3) { ts <- numeric(n); for (i in seq_len(n)) ts[i] <- system.time(f())[["elapsed"]]; ts }
out <- list(versions = list(R = R.version.string))
'''

R_TXIMPORT = R_COMMON + r'''
suppressMessages({ library(tximport); library(data.table) })
dir <- args[1]; outdir <- args[2]; reps <- as.integer(args[3])
samples <- readLines(file.path(dir, "samples.txt"))
files <- file.path(dir, samples, "quant.sf"); names(files) <- samples
tx2gene <- read.delim(file.path(dir, "tx2gene.tsv"), header = FALSE, col.names = c("tx", "gene"))
out$versions$tximport <- as.character(packageVersion("tximport"))
out$versions$readr <- if (requireNamespace("readr", quietly = TRUE)) as.character(packageVersion("readr")) else NA
# read (txOut): the files into transcript-level matrices
out$read <- tm(function() tximport(files, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE), reps)
txi <- tximport(files, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)
fwrite(as.data.table(txi$counts), file.path(outdir, "tx_counts.tsv"), sep = "\t")
fwrite(as.data.table(txi$abundance), file.path(outdir, "tx_abundance.tsv"), sep = "\t")
fwrite(as.data.table(txi$length), file.path(outdir, "tx_length.tsv"), sep = "\t")
writeLines(rownames(txi$counts), file.path(outdir, "tx_ids.txt"))
# summarise: the three countsFromAbundance modes, from the transcript-level object
for (cfa in c("no", "scaledTPM", "lengthScaledTPM")) {
  out[[paste0("summarize_", cfa)]] <- tm(function() summarizeToGene(txi, tx2gene, countsFromAbundance = cfa), reps)
  g <- summarizeToGene(txi, tx2gene, countsFromAbundance = cfa)
  fwrite(as.data.table(g$counts), file.path(outdir, paste0("gene_counts_", cfa, ".tsv")), sep = "\t")
  if (cfa == "no") {
    fwrite(as.data.table(g$length), file.path(outdir, "gene_length.tsv"), sep = "\t")
    fwrite(as.data.table(g$abundance), file.path(outdir, "gene_abundance.tsv"), sep = "\t")
    writeLines(rownames(g$counts), file.path(outdir, "gene_ids.txt"))
  }
}
# the whole thing end to end, as people call it
out$tximport_end_to_end <- tm(function() tximport(files, type = "salmon", tx2gene = tx2gene,
                                                  countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE), reps)
if (requireNamespace("DESeq2", quietly = TRUE)) {
  suppressMessages(library(DESeq2))
  out$versions$DESeq2 <- as.character(packageVersion("DESeq2"))
  x <- as.matrix(fread(file.path(dir, "gene_counts.tsv"), header = TRUE), rownames = 1)
  out$size_factors <- tm(function() estimateSizeFactorsForMatrix(x), reps)
  writeLines(format(estimateSizeFactorsForMatrix(x), digits = 15), file.path(outdir, "size_factors.txt"))
}
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA), file.path(outdir, "tximport.json"))
'''

R_EDGER = R_COMMON + r'''
suppressMessages({ library(edgeR); library(data.table) })
dir <- args[1]; outdir <- args[2]; reps <- as.integer(args[3])
out$versions$edgeR <- as.character(packageVersion("edgeR"))
out$versions$data.table <- as.character(packageVersion("data.table"))
x <- as.matrix(fread(file.path(dir, "gene_counts.tsv"), header = TRUE), rownames = 1)
group <- factor(readLines(file.path(dir, "group.txt")))
y <- DGEList(x, group = group)
out$tmm <- tm(function() calcNormFactors(y), reps)
yn <- calcNormFactors(y)
writeLines(format(yn$samples$norm.factors, digits = 15), file.path(outdir, "tmm.txt"))
out$logcpm <- tm(function() cpm(yn, log = TRUE, prior.count = 2), reps)
fwrite(as.data.table(cpm(yn, log = TRUE, prior.count = 2)), file.path(outdir, "logcpm.tsv"), sep = "\t")
out$filter <- tm(function() filterByExpr(y, group = group), reps)
writeLines(as.character(as.integer(filterByExpr(y, group = group))), file.path(outdir, "filter.txt"))
# reading a featureCounts file
fc <- file.path(dir, "featureCounts.txt")
out$fc_read_delim <- tm(function() read.delim(fc, comment.char = "#", check.names = FALSE), reps)
out$fc_fread <- tm(function() fread(fc, skip = 1), reps)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA), file.path(outdir, "edger.json"))
'''


def _rscript_with(packages):
    """The first Rscript on GB_RSCRIPT (or PATH) that has every package; None otherwise."""
    for r in RSCRIPTS:
        exe = shutil.which(r) or (r if os.path.exists(r) else None)
        if not exe:
            continue
        code = "cat(all(sapply(c(%s), requireNamespace, quietly=TRUE)))" % ", ".join(f'"{p}"' for p in packages)
        try:
            ok = subprocess.run([exe, "-e", code], capture_output=True, text=True, timeout=120).stdout.strip()
        except (subprocess.SubprocessError, OSError):
            continue
        if ok.endswith("TRUE"):
            return exe
    return None


def run_r(script_text, name, args, rscript):
    path = EXPR / f"{name}.R"
    path.write_text(script_text)
    t0 = time.perf_counter()
    r = subprocess.run([rscript] + [str(path)] + [str(a) for a in args], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{name}.R failed:\n{r.stderr[-3000:]}")
    return time.perf_counter() - t0


def _read_tsv_matrix(path):
    import polars as pl
    return pl.read_csv(path, separator="\t").to_numpy().astype(float)


# ── data ────────────────────────────────────────────────────────────────────

def gtf(n_genes: int) -> Path:
    out = EXPR / f"genes_{n_genes}.gtf"
    if out.exists():
        return out
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import make_data
    if not (DATA / "hg38.chrom.sizes").exists():
        make_data.chromsizes()
    keep = make_data.DATA
    make_data.DATA = EXPR
    try:
        make_data.gtf(n_genes)
    finally:
        make_data.DATA = keep
    (EXPR / "genes.gtf").rename(out)
    return out


def quants(n_genes: int, n_samples: int, G: Genes) -> Path:
    """salmon-style quant.sf per sample (ids without versions), tx2gene.tsv,
    samples.txt, group.txt, gene_counts.tsv (summed) and featureCounts.txt."""
    d = EXPR / f"salmon_{n_genes}_{n_samples}"
    if (d / "featureCounts.txt").exists():
        return d
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    T, Gt = G.transcripts, G.genes
    tx_ids = np.array([t.rsplit(".", 1)[0] for t in T.cols["transcript_id"]], dtype=object)
    gene_of = T.cols["gene"]
    n_tx, n_g = len(T), len(Gt)
    length = np.clip(rng.lognormal(np.log(1800), 0.7, n_tx), 200, 30_000).round()
    gene_mu = rng.lognormal(3.0, 1.6, n_g)
    share = rng.gamma(1.0, 1.0, n_tx)
    tot = np.zeros(n_g)
    np.add.at(tot, gene_of, share)
    mu_tx = gene_mu[gene_of] * share / tot[gene_of]
    lib = rng.lognormal(0, 0.35, n_samples)
    samples = [f"{'WT' if j < n_samples // 2 else 'KO'}_{j + 1}" for j in range(n_samples)]
    counts = np.empty((n_tx, n_samples))
    for j in range(n_samples):
        m = mu_tx * lib[j]
        counts[:, j] = rng.negative_binomial(8, 8 / (8 + m))
    eff = np.maximum(length - 180, 1.0)[:, None] + rng.uniform(-5, 5, (1, n_samples)).round(2)
    rate = counts / eff
    tpm = rate / rate.sum(axis=0) * 1e6
    for j, s in enumerate(samples):
        (d / s).mkdir(exist_ok=True)
        with open(d / s / "quant.sf", "w") as f:
            f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
            f.writelines(f"{tx_ids[i]}\t{int(length[i])}\t{eff[i, j]:.3f}\t{tpm[i, j]:.6f}\t{counts[i, j]:.3f}\n"
                         for i in range(n_tx))
    (d / "samples.txt").write_text("\n".join(samples) + "\n")
    (d / "group.txt").write_text("\n".join(s.split("_")[0] for s in samples) + "\n")
    with open(d / "tx2gene.tsv", "w") as f:
        f.writelines(f"{tx_ids[i]}\t{Gt.cols['gene_id'][gene_of[i]]}\n" for i in range(n_tx))
    gc = np.zeros((n_g, n_samples))
    np.add.at(gc, gene_of, counts)
    with open(d / "gene_counts.tsv", "w") as f:
        f.write("gene\t" + "\t".join(samples) + "\n")
        f.writelines(f"{Gt.cols['gene_id'][i]}\t" + "\t".join(str(int(v)) for v in gc[i]) + "\n" for i in range(n_g))
    with open(d / "featureCounts.txt", "w") as f:
        f.write('# Program:featureCounts v2.0.6; Command:"featureCounts" "-a" "genes.gtf"\n')
        f.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + "\t".join(f"bam/{s}.sorted.bam" for s in samples) + "\n")
        chroms, starts, ends, strands = Gt.chroms, Gt.starts, Gt.ends, Gt.strand
        f.writelines(f"{Gt.cols['gene_id'][i]}\t{chroms[i]}\t{starts[i] + 1}\t{ends[i]}\t{strands[i]}\t"
                     f"{ends[i] - starts[i]}\t" + "\t".join(str(int(v)) for v in gc[i]) + "\n" for i in range(n_g))
    return d


def diff(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if a.shape != b.shape or not ok.any():
        return {"max_abs_diff": None, "max_rel_diff": None, "agrees": False}
    d = np.abs(a[ok] - b[ok])
    rel = d / np.maximum(np.abs(b[ok]), 1e-300)
    return {"max_abs_diff": float(d.max()), "max_rel_diff": float(rel.max()), "agrees": bool(d.max() < 1e-6)}


# ── the bench ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    rec = Recorder("expression")
    r_versions = {}
    r_tx = _rscript_with(["tximport", "data.table", "jsonlite"])
    r_edger = _rscript_with(["edgeR", "data.table", "jsonlite"])
    print(f"[R] tximport/DESeq2: {r_tx or 'not found'} · edgeR: {r_edger or 'not found'}")
    tables = [b for b in gb.backends.families()["tables"] if gb.backends.installed("tables", b)]
    for n_genes, n_samples in SIZES:
        reps = 3 if n_samples <= 24 else 1
        G = Genes.make(str(gtf(n_genes)))
        d = quants(n_genes, n_samples, G)
        samples = [ln.strip() for ln in open(d / "samples.txt")]
        glob = [str(d / s / "quant.sf") for s in samples]          # a list keeps this order (a glob sorts)
        n_tx = len(G.transcripts)
        size = dict(n_genes=len(G), n_transcripts=n_tx, n_samples=n_samples)
        print(f"=== {len(G):,} genes · {n_tx:,} transcripts · {n_samples} samples ===")

        # read: the two table backends
        ref = None
        for b in tables:
            t = timeit(lambda: Expression.make(glob, backend=b, verbose=False), repeat=reps, warmup=1)
            E = Expression.make(glob, backend=b, verbose=False)
            ex = {"agrees": True} if ref is None else diff(E.counts, ref.counts) | {
                "agrees": all(np.array_equal(E.layers[k], ref.layers[k]) for k in ref.layers)}
            ref = ref or E
            rec.add(part="read", step="read quant.sf files", engine=f"genomeblocks ({b})", seconds=t["median"],
                    runs=t["runs"], n=n_tx * n_samples, **size, **ex)
        E_tx = ref

        # align + summarise in genomeblocks
        t = timeit(lambda: E_tx.align(genes=G, level="transcript", verbose=False), repeat=reps)
        A = E_tx.align(genes=G, level="transcript", verbose=False)
        rec.add(part="align", step="align to genes.transcripts", engine="genomeblocks", seconds=t["median"],
                runs=t["runs"], n=n_tx, rule=A.meta["key_rule"], **size)
        gb_genes = {}
        for how in ("counts", "scaled_tpm", "length_scaled_tpm"):
            t = timeit(lambda: A.to_genes(how=how, verbose=False), repeat=reps)
            gb_genes[how] = A.to_genes(how=how, verbose=False)
            rec.add(part="to_genes", step=f"to_genes ({how})", engine="genomeblocks", seconds=t["median"],
                    runs=t["runs"], n=n_tx, **size)
        t = timeit(lambda: Expression.make(glob, genes=G, verbose=False), repeat=reps)
        rec.add(part="to_genes", step="read + align + to_genes (length_scaled_tpm)", engine="genomeblocks",
                seconds=t["median"], runs=t["runs"], n=n_tx * n_samples, **size)
        Eg = gb_genes["counts"]

        # normalisation in genomeblocks (on the summed gene counts, as R gets them)
        group = np.array([s.split("_")[0] for s in Eg.sample_names])
        norm = {}
        for step, fn in (("tmm", lambda: Eg.tmm()), ("size_factors", lambda: Eg.size_factors()),
                         ("log_cpm (tmm)", lambda: Eg.log_cpm(norm="tmm")),
                         ("filter_by_expr", lambda: Eg.filter_by_expr(group))):
            t = timeit(fn, repeat=reps)
            norm[step] = fn()
            rec.add(part="normalise", step=step, engine="genomeblocks", seconds=t["median"], runs=t["runs"],
                    n=len(Eg) * n_samples, **size)

        # featureCounts file: the two backends
        fc = str(d / "featureCounts.txt")
        fref = None
        for b in tables:
            t = timeit(lambda: Expression.make(fc, backend=b, verbose=False), repeat=reps)
            F = Expression.make(fc, backend=b, verbose=False)
            ex = {"agrees": True} if fref is None else {"agrees": np.array_equal(F.counts, fref.counts)}
            fref = fref or F
            rec.add(part="featurecounts", step="read featureCounts matrix", engine=f"genomeblocks ({b})",
                    seconds=t["median"], runs=t["runs"], n=len(F) * n_samples, **size, **ex)
        assert np.array_equal(fref.take_samples(Eg.sample_names).counts, Eg.counts), \
            "featureCounts matrix differs from the summed counts"

        # R: tximport (+ DESeq2 size factors)
        rout = d / "r_out"
        rout.mkdir(exist_ok=True)
        if r_tx:
            print(f"  [R] tximport on {n_samples} files ...", flush=True)
            run_r(R_TXIMPORT, "expr_tximport", [d, rout, reps], r_tx)
            j = json.loads((rout / "tximport.json").read_text())
            r_versions.update(j["versions"])
            tx_ids = [ln.strip() for ln in open(rout / "tx_ids.txt")]
            order = np.array([E_tx.row(i) for i in tx_ids])          # R's rows onto ours (same order, checked)
            ex = diff(E_tx.counts[order], _read_tsv_matrix(rout / "tx_counts.tsv"))
            ex_a = diff(E_tx.tpm[order], _read_tsv_matrix(rout / "tx_abundance.tsv"))
            ex_l = diff(E_tx.length[order], _read_tsv_matrix(rout / "tx_length.tsv"))
            rec.add(part="read", step="read quant.sf files", engine="R tximport (txOut)", seconds=float(np.median(j["read"])),
                    runs=j["read"], n=n_tx * n_samples, **size, **ex,
                    max_abs_diff_abundance=ex_a["max_abs_diff"], max_abs_diff_length=ex_l["max_abs_diff"])
            gids = [ln.strip() for ln in open(rout / "gene_ids.txt")]
            gorder = G.rows(gids)
            for how, cfa in (("counts", "no"), ("scaled_tpm", "scaledTPM"), ("length_scaled_tpm", "lengthScaledTPM")):
                rc = _read_tsv_matrix(rout / f"gene_counts_{cfa}.tsv")
                ex = diff(gb_genes[how].counts[gorder], rc)
                extra = {}
                if how == "counts":
                    extra["max_abs_diff_length"] = diff(gb_genes[how].length[gorder], _read_tsv_matrix(rout / "gene_length.tsv"))["max_abs_diff"]
                    extra["max_abs_diff_abundance"] = diff(gb_genes[how].tpm[gorder], _read_tsv_matrix(rout / "gene_abundance.tsv"))["max_abs_diff"]
                k = f"summarize_{cfa}"
                rec.add(part="to_genes", step=f"to_genes ({how})", engine="R tximport::summarizeToGene",
                        seconds=float(np.median(j[k])), runs=j[k], n=n_tx, **size, **ex, **extra)
            rec.add(part="to_genes", step="read + align + to_genes (length_scaled_tpm)", engine="R tximport",
                    seconds=float(np.median(j["tximport_end_to_end"])), runs=j["tximport_end_to_end"],
                    n=n_tx * n_samples, **size)
            if "size_factors" in j:
                sf = np.array([float(v) for v in open(rout / "size_factors.txt").read().split()])
                rec.add(part="normalise", step="size_factors", engine="R DESeq2", seconds=float(np.median(j["size_factors"])),
                        runs=j["size_factors"], n=len(Eg) * n_samples, **size, **diff(norm["size_factors"], sf))
        else:
            rec.add(part="read", step="skipped", engine="R tximport", reason="no Rscript with tximport", **size)

        # R: edgeR + file reading
        if r_edger:
            print(f"  [R] edgeR on {len(Eg):,} x {n_samples} ...", flush=True)
            run_r(R_EDGER, "expr_edger", [d, rout, reps], r_edger)
            j = json.loads((rout / "edger.json").read_text())
            r_versions.update(j["versions"])
            tmm = np.array([float(v) for v in open(rout / "tmm.txt").read().split()])
            rec.add(part="normalise", step="tmm", engine="R edgeR calcNormFactors", seconds=float(np.median(j["tmm"])),
                    runs=j["tmm"], n=len(Eg) * n_samples, **size, **diff(norm["tmm"], tmm))
            lc = _read_tsv_matrix(rout / "logcpm.tsv")
            rec.add(part="normalise", step="log_cpm (tmm)", engine="R edgeR cpm(log=TRUE)",
                    seconds=float(np.median(j["logcpm"])), runs=j["logcpm"], n=len(Eg) * n_samples, **size,
                    **diff(norm["log_cpm (tmm)"], lc))
            flt = np.array([int(v) for v in open(rout / "filter.txt").read().split()], bool)
            rec.add(part="normalise", step="filter_by_expr", engine="R edgeR filterByExpr",
                    seconds=float(np.median(j["filter"])), runs=j["filter"], n=len(Eg) * n_samples, **size,
                    agrees=bool((flt == norm["filter_by_expr"]).all()),
                    disagreements=int((flt != norm["filter_by_expr"]).sum()))
            for step, key in (("R read.delim", "fc_read_delim"), ("R data.table::fread", "fc_fread")):
                rec.add(part="featurecounts", step="read featureCounts matrix", engine=step,
                        seconds=float(np.median(j[key])), runs=j[key], n=len(Eg) * n_samples, **size)
        else:
            rec.add(part="normalise", step="skipped", engine="R edgeR", reason="no Rscript with edgeR", **size)

    rec.save(sizes=[{"n_genes": a, "n_samples": b} for a, b in SIZES], r_versions=r_versions,
             rscript={"tximport": r_tx, "edgeR": r_edger})
