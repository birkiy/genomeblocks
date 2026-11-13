#!/usr/bin/env python3
"""The landscape benchmark: genomeblocks against the tools of every theme that
the main suite does not cover, scored the way dynverse scored trajectory
methods (accuracy · scalability · stability · usability).

Themes (one part each) and the engines added here — the main suite's engines
(bedtools, pyranges 0.x, bioframe, deepTools, FIMO, MOODS, cooler, tximport,
edgeR, DESeq2, ...) are read from their own results by build_landscape.py:

  intervals     R: GenomicRanges, plyranges, valr, data.table::foverlaps;
                Python: pyranges 1.x, polars-bio, BiocPy genomicranges,
                superintervals, AIList, and (for the jobs added 2026-10-23:
                overlap counts, coverage, genome depth, subtract, complement)
                bioframe and pyranges 0.x; CLI: BEDOPS, bedtools (reference)
  annotation    R: txdbmaker, rtracklayer, GenomicFeatures (DIY labels),
                ChIPseeker, ChIPpeakAnno; CLI: HOMER annotatePeaks.pl, UROPA,
                bedtools closest (reference)
  signal        R: rtracklayer + EnrichedHeatmap, genomation;
                CLI: bigWigAverageOverBed; pyBigWig per-base (reference)
  motifs        R: universalmotif, Biostrings, motifmatchr; CLI: HOMER homer2
  cooccupancy   R: LOLA, regioneR; CLI: GIGGLE, bedtools fisher
  graph         R: GenomicInteractions / InteractionSet, igraph; CLI: bedtools pairtobed
  expression    Python: pytximport, pydeseq2, rnanorm
  signal_bigwig genomeblocks' bam_to_bigwig and Loci.to_bigwig on signal C1 / C9's inputs
  hic_loops     gb.mustache / gb.hiccups next to Mustache 1.3.3 and Juicer HiCCUPS (MUSTACHE, JUICER_TOOLS)
  stability     every theme's main task on the same data presented differently
                (rows shuffled; chromosomes in another order): stable = identical answer

Every row records the median time, the runs, and ``accuracy`` in [0, 1]: 1 when
the normalised answer equals the reference (sets, labels, gaps), else the
agreement (Pearson / Spearman / fraction of rows agreeing).

Environment: GB_LANDSCAPE_R (an Rscript with the Bioconductor packages),
GB_LANDSCAPE_PY (a python with the Python competitors; its bin/ also holds
BEDOPS and bigWigAverageOverBed), GB_LANDSCAPE_CLI (a bin/ with featureCounts,
htseq-count, megadepth, wiggletools, intervene, bedGraphToBigWig),
GB_LANDSCAPE_HIC (an env with cooler, cooltools, hictk, HiCExplorer), GIGGLE
(binary), GIGGLE_LD (dir with libhts), HOMER_BIN (dir with annotatePeaks.pl /
homer2), MUSTACHE / JUICER_TOOLS / JAVA (the hic_loops part's originals);
bedtools, bgzip on PATH.
GB_BENCH_PARTS=intervals,signal runs some parts (other parts' rows are kept);
GB_BENCH_SIZES=10000,100000 limits the interval sizes.
"""
from __future__ import annotations

import contextlib
import glob
import io
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from common import DATA, Recorder, timeit, tool

import genomeblocks as gb
from genomeblocks import Loci, Genes, Pairs, Atlas, Architecture, Expression

HERE = Path(__file__).resolve().parent
WORK = Path(os.environ.get("GB_LANDSCAPE_WORK", DATA / "landscape_work"))
WORK.mkdir(parents=True, exist_ok=True)
RSCRIPT = os.environ.get("GB_LANDSCAPE_R", shutil.which("Rscript"))
PYTHON = os.environ.get("GB_LANDSCAPE_PY")
PYBIN = Path(PYTHON).parent if PYTHON else None
GIGGLE = os.environ.get("GIGGLE")
GIGGLE_LD = os.environ.get("GIGGLE_LD")
HOMER_BIN = Path(os.environ.get("HOMER_BIN", Path(tool("homer2")).parent))
SIZES = [int(x) for x in os.environ.get("GB_BENCH_SIZES", "10000,100000,1000000").split(",")]
CS = str(DATA / "hg38.chrom.sizes")
REPEAT = int(os.environ.get("GB_REPEAT", "3"))


# ── helpers ──────────────────────────────────────────────────────────────────

def quiet(fn):
    def run():
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return fn()
    return run


def worker(kind, task, outdir, *args, env=None):
    """Run worker.R / worker_py.py; return its JSON rows (an unusable worker -> [])."""
    if kind == "R":
        if not RSCRIPT:
            return []
        cmd = [RSCRIPT, str(HERE / "landscape" / "worker.R"), task, str(outdir), *map(str, args)]
    else:
        if not PYTHON:
            return []
        cmd = [PYTHON, str(HERE / "landscape" / "worker_py.py"), task, str(outdir), *map(str, args)]
    e = {**os.environ, "GB_REPEAT": str(REPEAT), **(env or {})}
    p = subprocess.run(cmd, capture_output=True, text=True, env=e)
    rows = []
    for line in p.stdout.splitlines():
        if line.startswith("{"):
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    if p.returncode != 0 and not rows:
        print(f"  [{kind} worker {task} failed] {p.stderr[-800:]}", flush=True)
    return rows


def record_worker(rec, part, rows, n=None, score=None, **extra):
    """Add worker rows; ``score(row) -> (accuracy, detail)`` reads the engine's answer file."""
    for r in rows:
        if r.get("error"):
            rec.add(part=part, task=r["op"], engine=r["engine"], n=n, seconds=None, runs=None,
                    accuracy=None, error=r["error"][:300], **extra)
            continue
        try:
            acc, detail = (None, "") if score is None else score(r)
        except Exception as e:                                   # a scoring failure is recorded, not fatal
            acc, detail = None, f"scoring failed: {type(e).__name__}: {str(e)[:200]}"
        rec.add(part=part, task=r["op"], engine=r["engine"], n=n, seconds=float(np.median(r["seconds"])),
                runs=r["seconds"], n_out=r.get("n_out"), accuracy=acc, detail=detail, **extra)


def bed_set(df_or_path):
    df = pd.read_csv(df_or_path, sep="\t", header=None, usecols=[0, 1, 2]) if isinstance(df_or_path, (str, Path)) else df_or_path
    df = df.iloc[:, :3]
    return set(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(int), df.iloc[:, 2].astype(int)))


def loci_set(L):
    return set(zip(L.chroms.astype(str), L.starts.astype(int), L.ends.astype(int)))


def set_score(ref, got):
    """1 when identical, else the Jaccard similarity."""
    if ref == got:
        return 1.0, f"identical ({len(ref):,})"
    j = len(ref & got) / max(1, len(ref | got))
    return j, f"differs: {len(got):,} vs {len(ref):,} (Jaccard {j:.4f})"


def gap(a_s, a_e, b_s, b_e):
    """Bases between two intervals: 0 when they overlap or are book-ended."""
    return np.maximum(0, np.maximum(a_s, b_s) - np.minimum(a_e, b_e))


def gap_score(ref_gap, idx, bstart, bend, A):
    """Per-row nearest distance agreement against the reference gaps (-1 = none)."""
    idx = np.asarray(idx, int)
    g = gap(A.starts[idx], A.ends[idx], np.asarray(bstart, np.int64), np.asarray(bend, np.int64))
    ok = ref_gap[idx] >= 0
    same = (g[ok] == ref_gap[idx][ok])
    frac = float(same.mean()) if ok.any() else 0.0
    return (1.0 if frac == 1.0 else frac), f"{same.sum():,}/{ok.sum():,} rows agree"


def coord_rows(L):
    """coordinate -> every row of ``L`` with it (BED sets can repeat a coordinate)."""
    out = {}
    for i, k in enumerate(zip(L.chroms.astype(str), L.starts.astype(int), L.ends.astype(int))):
        out.setdefault(k, []).append(i)
    return out


def scatter_by_coords(key, frame, values, n, fill=-1):
    """Values of a sorted-copy output back onto the rows of the original table, duplicates included."""
    out = np.full(n, fill, np.int64)
    for r, v in zip(frame.itertuples(index=False), values):
        for i in key[(str(r[0]), int(r[1]), int(r[2]))]:
            out[i] = v
    return out


def read_tsv(path, names):
    return pd.read_csv(path, sep="\t", header=None, names=names)


def corr_score(ref, got, tol=1e-6, rtol=None):
    """1 when every value agrees within ``tol`` × the largest reference value (or, with
    ``rtol``, element-wise within that relative tolerance); else Pearson r."""
    ref, got = np.asarray(ref, float).ravel(), np.asarray(got, float).ravel()
    if ref.shape != got.shape:
        return 0.0, f"shape {got.shape} vs {ref.shape}"
    mad = float(np.nanmax(np.abs(ref - got)))
    scale = float(np.nanmax(np.abs(ref))) or 1.0
    if rtol is not None and float(np.nanmax(np.abs(ref - got) / (np.abs(ref) + 1e-9))) <= rtol:
        return 1.0, f"max |Δ| {mad:.2e} (within {rtol:g} relative)"
    if mad <= tol * scale:
        return 1.0, f"max |Δ| {mad:.2e}"
    r = float(np.corrcoef(np.nan_to_num(ref), np.nan_to_num(got))[0, 1])
    return max(0.0, r), f"Pearson {r:.5f}, max |Δ| {mad:.3g}"


def spearman_score(ref, got):
    from scipy.stats import spearmanr
    rho = float(spearmanr(ref, got).statistic)
    return max(0.0, rho), f"Spearman {rho:.4f}"


def coord_values(df, col):
    """{(chrom, start, end): value} from a table whose first three columns are the interval."""
    return dict(zip(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(int), df.iloc[:, 2].astype(int)),
                    df.iloc[:, col].astype(float).round().astype(np.int64)))


def value_score(ref, got):
    """Share of the reference intervals whose value the engine reproduces (a missing row counts as wrong)."""
    same = sum(1 for k, v in ref.items() if got.get(k) == v)
    frac = same / max(1, len(ref))
    return (1.0 if same == len(ref) else frac), f"{same:,}/{len(ref):,} intervals agree"


def merged_set(df):
    """The bases an interval table covers, as merged (chrom, start, end) runs: subtract and complement are
    compared on the bases they leave, so a tool that merges its answer and one that keeps every A row agree."""
    d = df.iloc[:, :3].copy()
    d.columns = ["chrom", "start", "end"]
    d = d[d["end"] > d["start"]].sort_values(["chrom", "start", "end"])
    out = set()
    for chrom, g in d.groupby("chrom", sort=False):
        s, e = g["start"].to_numpy(np.int64), np.maximum.accumulate(g["end"].to_numpy(np.int64))
        first = np.flatnonzero(np.r_[True, s[1:] > e[:-1]])           # book-ended rows fuse
        last = np.r_[first[1:], len(s)] - 1
        out.update(zip([str(chrom)] * len(first), s[first].tolist(), e[last].tolist()))
    return out


def depth_runs(df):
    """A depth track as (chrom, start, end, depth) runs with depth > 0, touching runs of equal depth merged,
    so two tools that split the same depth function differently agree."""
    d = df.iloc[:, :4].copy()
    d.columns = ["chrom", "start", "end", "depth"]
    d["depth"] = d["depth"].astype(float).round().astype(np.int64)
    d = d[d["depth"] > 0].sort_values(["chrom", "start"])
    out = set()
    for chrom, g in d.groupby("chrom", sort=False):
        s, e, v = (g[c].to_numpy(np.int64) for c in ("start", "end", "depth"))
        first = np.flatnonzero(np.r_[True, (s[1:] != e[:-1]) | (v[1:] != v[:-1])])
        last = np.r_[first[1:], len(s)] - 1
        out.update(zip([str(chrom)] * len(first), s[first].tolist(), e[last].tolist(), v[first].tolist()))
    return out


def write_bed(L, path, names=None):
    with open(path, "w") as f:
        for i, (c, s, e) in enumerate(zip(L.chroms, L.starts, L.ends)):
            f.write(f"{c}\t{s}\t{e}\t{names[i] if names is not None else i}\n")
    return str(path)


def cli(cmd, **kw):
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        raise RuntimeError(f"{cmd[0]}: {p.stderr[-500:]}")
    return p.stdout


def time_cli(fn, times=REPEAT, warmup=1):
    t = timeit(fn, repeat=times, warmup=warmup)
    return t["median"], t["runs"]


# ── intervals ────────────────────────────────────────────────────────────────

def part_intervals(rec):
    print("\n== intervals ==")
    for n in SIZES:
        a_bed, b_bed = str(DATA / f"peaks_A_{n}.bed"), str(DATA / f"peaks_B_{n}.bed")
        # sort=False keeps the file's row order, which is the row index the workers report (Loci.make would
        # reorder rows that share a start)
        A, B = Loci.make(a_bed, sort=False), Loci.make(b_bed, sort=False)
        rep = REPEAT if n <= 100_000 else max(1, REPEAT - 1)
        # bedtools: the independent reference for the sets
        bt = tool("bedtools")
        out = {}
        t, runs = time_cli(lambda: out.__setitem__("i", cli([bt, "intersect", "-u", "-a", a_bed, "-b", b_bed])), rep)
        ref_int = bed_set(pd.read_csv(io.StringIO(out["i"]), sep="\t", header=None, usecols=[0, 1, 2]))
        rec.add(part="intervals", task="intersect", engine="bedtools intersect -u (CLI)", n=n, seconds=t, runs=runs,
                accuracy=1.0, detail="reference", n_out=len(ref_int))
        t, runs = time_cli(lambda: out.__setitem__("m", cli([bt, "merge", "-i", a_bed])), rep)
        ref_mrg = bed_set(pd.read_csv(io.StringIO(out["m"]), sep="\t", header=None))
        rec.add(part="intervals", task="merge", engine="bedtools merge (CLI)", n=n, seconds=t, runs=runs,
                accuracy=1.0, detail="reference", n_out=len(ref_mrg))
        # genomeblocks (cold index, like the backends bench)
        res = {}
        t = timeit(lambda: res.__setitem__("i", A.take(A.overlap_any(B))), repeat=rep, setup=lambda: (A._dirty(), B._dirty()))
        acc, d = set_score(ref_int, loci_set(res["i"]))
        rec.add(part="intervals", task="intersect", engine="genomeblocks", n=n, seconds=t["median"], runs=t["runs"],
                accuracy=acc, detail=d, n_out=len(res["i"]))
        t = timeit(lambda: res.__setitem__("m", A.merge()), repeat=rep)
        acc, d = set_score(ref_mrg, loci_set(res["m"]))
        rec.add(part="intervals", task="merge", engine="genomeblocks", n=n, seconds=t["median"], runs=t["runs"],
                accuracy=acc, detail=d, n_out=len(res["m"]))
        t = timeit(lambda: res.__setitem__("d", A.nearest(B)[1]), repeat=rep, setup=lambda: (A._dirty(), B._dirty()))
        ref_gap = np.asarray(res["d"], np.int64)
        # bedtools closest as the independent reference for the gaps (it needs both files sorted the same
        # way: sorted copies made outside the clock, rows mapped back by coordinates)
        od = WORK / f"intervals_{n}"
        od.mkdir(parents=True, exist_ok=True)
        sa, sb = str(od / "A.sorted.bed"), str(od / "B.sorted.bed")
        for src, dst in ((a_bed, sa), (b_bed, sb)):
            with open(dst, "w") as f:
                subprocess.run(["sort", "-k1,1", "-k2,2n", src], stdout=f, check=True)
        key = coord_rows(A)
        t2, runs2 = time_cli(lambda: out.__setitem__("c", cli([bt, "closest", "-t", "first", "-a", sa, "-b", sb])), rep)
        c = pd.read_csv(io.StringIO(out["c"]), sep="\t", header=None)
        nb = c.shape[1] // 2
        hit = (c.iloc[:, nb + 1] >= 0).to_numpy()
        g = gap(c.iloc[:, 1].to_numpy(), c.iloc[:, 2].to_numpy(), c.iloc[:, nb + 1].to_numpy(), c.iloc[:, nb + 2].to_numpy())
        g = np.where(hit, g, -1)
        bt_gap = scatter_by_coords(key, c, g, n)
        same = (bt_gap == ref_gap)
        rec.add(part="intervals", task="nearest", engine="bedtools closest (CLI)", n=n, seconds=t2, runs=runs2,
                accuracy=1.0, detail="reference", n_out=int(hit.sum()))
        rec.add(part="intervals", task="nearest", engine="genomeblocks", n=n, seconds=t["median"], runs=t["runs"],
                accuracy=float(same.mean()), detail=f"{same.sum():,}/{len(same):,} rows agree with bedtools closest")
        # the confirmed core jobs added on 2026-10-23 (landscape/jobs.json): overlap counts (C2), coverage of
        # A by B (C4), genome-wide depth (C5), subtract (C7), complement (C8). bedtools is the reference for
        # each; only a tool's own function for the job is timed, never a recipe composed from other
        # primitives. bedtools genomecov / complement need sorted input: the sorted copy is made outside the
        # clock, like BEDOPS' sort-bed copies
        refs = {}
        a_sorted, g_lex = str(od / "A.sorted.bed"), str(od / "genome.lex.txt")      # the genome file in the same order
        with open(g_lex, "w") as f:
            subprocess.run(["sort", "-k1,1", CS], stdout=f, check=True)
        bt_new = [("overlap_count", "bedtools intersect -c (CLI)", [bt, "intersect", "-c", "-a", a_bed, "-b", b_bed]),
                  ("coverage", "bedtools coverage (CLI)", [bt, "coverage", "-a", a_bed, "-b", b_bed]),
                  ("genome_depth", "bedtools genomecov -bg (CLI)", [bt, "genomecov", "-bg", "-i", a_sorted, "-g", g_lex]),
                  ("subtract", "bedtools subtract (CLI)", [bt, "subtract", "-a", a_bed, "-b", b_bed]),
                  ("complement", "bedtools complement (CLI)", [bt, "complement", "-i", a_sorted, "-g", g_lex])]
        for task, eng, cmd in bt_new:
            t, runs = time_cli(lambda cmd=cmd: out.__setitem__(task, cli(cmd)), rep)
            df = pd.read_csv(io.StringIO(out[task]), sep="\t", header=None)
            refs[task] = {"overlap_count": lambda: coord_values(df, -1), "coverage": lambda: coord_values(df, -3),
                          "genome_depth": lambda: depth_runs(df)}.get(task, lambda: merged_set(df))()
            rec.add(part="intervals", task=task, engine=eng, n=n, seconds=t, runs=runs, accuracy=1.0, detail="reference", n_out=len(df))

        def score_new(r, df=None):
            df = pd.read_csv(r["out"], sep="\t", header=None) if df is None else df
            if r["op"] in ("overlap_count", "coverage"):
                return value_score(refs[r["op"]], coord_values(df, 3))
            if r["op"] == "genome_depth":
                return set_score(refs[r["op"]], depth_runs(df))
            return set_score(refs[r["op"]], merged_set(df))

        def run_inproc(task, engine, fn):
            try:
                t = timeit(lambda: res.__setitem__(task, fn()), repeat=rep)
                acc, d = score_new({"op": task}, res[task])
                rec.add(part="intervals", task=task, engine=engine, n=n, seconds=t["median"], runs=t["runs"], accuracy=acc, detail=d,
                        n_out=len(res[task]))
            except Exception as e:                               # one engine failing must not lose the others
                rec.add(part="intervals", task=task, engine=engine, n=n, seconds=None, runs=None, accuracy=None, error=f"{type(e).__name__}: {e}"[:300])

        # genomeblocks: genome-wide depth (genomeblocks.hichip.coverage takes any intervals) and, since 2026-10-24,
        # Loci.count_overlaps / coverage / subtract / complement on its own numpy kernels
        from genomeblocks import hichip
        # the native answer is timed (an array aligned to A's rows, or a Loci); the table for scoring is built after
        def run_gb(task, engine, call, to_frame):
            try:
                t = timeit(lambda: res.__setitem__(task, call()), repeat=rep)
                df = to_frame(res[task])
                acc, d = score_new({"op": task}, df)
                rec.add(part="intervals", task=task, engine=engine, n=n, seconds=t["median"], runs=t["runs"], accuracy=acc, detail=d, n_out=len(df))
            except Exception as e:
                rec.add(part="intervals", task=task, engine=engine, n=n, seconds=None, runs=None, accuracy=None, error=f"{type(e).__name__}: {e}"[:300])
        per_row = lambda v: pd.DataFrame({"c": A.chroms, "s": A.starts, "e": A.ends, "n": v})
        as_frame = lambda L: pd.DataFrame({"c": L.chroms, "s": L.starts, "e": L.ends})
        run_gb("overlap_count", "genomeblocks Loci.count_overlaps", lambda: A.count_overlaps(B), per_row)
        run_gb("coverage", "genomeblocks Loci.coverage", lambda: A.coverage(B), per_row)
        run_gb("subtract", "genomeblocks Loci.subtract", lambda: A.subtract(B), as_frame)
        run_gb("complement", "genomeblocks Loci.complement", lambda: A.complement(CS), as_frame)

        def _runs_frame(runs):                                   # (chrom, starts, ends, depths) per chromosome
            return pd.DataFrame({"chrom": np.repeat([r[0] for r in runs], [len(r[1]) for r in runs]),
                                 "start": np.concatenate([r[1] for r in runs]), "end": np.concatenate([r[2] for r in runs]),
                                 "depth": np.concatenate([r[3] for r in runs])})
        run_gb("genome_depth", "genomeblocks hichip.coverage", lambda: list(hichip.coverage(A)), _runs_frame)
        # bioframe and pyranges 0.x live in this environment (the main suite runs them through genomeblocks'
        # backends; here they are called directly)
        Adf = pd.read_csv(a_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
        Bdf = pd.read_csv(b_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
        sizes = pd.read_csv(CS, sep="\t", header=None, names=["chrom", "end"])
        view = pd.DataFrame({"chrom": sizes["chrom"], "start": 0, "end": sizes["end"], "name": sizes["chrom"]})
        try:
            import bioframe as bf
            name = f"bioframe {bf.__version__}"
            run_inproc("overlap_count", name, lambda: bf.count_overlaps(Adf, Bdf)[["chrom", "start", "end", "count"]])
            run_inproc("coverage", name, lambda: bf.coverage(Adf, Bdf)[["chrom", "start", "end", "coverage"]])
            run_inproc("subtract", name, lambda: bf.subtract(Adf, Bdf)[["chrom", "start", "end"]])
            run_inproc("complement", name, lambda: bf.complement(Adf, view_df=view)[["chrom", "start", "end"]])
        except ImportError:
            pass
        try:
            import pyranges as pr0
            if not pr0.__version__.startswith("1."):
                name = f"pyranges {pr0.__version__}"
                ren = {"chrom": "Chromosome", "start": "Start", "end": "End"}
                pa0, pb0 = pr0.PyRanges(Adf.rename(columns=ren)), pr0.PyRanges(Bdf.rename(columns=ren))

                def _cov0():
                    d = pa0.coverage(pb0).df
                    return pd.DataFrame({"c": d["Chromosome"], "s": d["Start"], "e": d["End"],
                                         "bases": (d["FractionOverlaps"] * (d["End"] - d["Start"])).round()})
                run_inproc("overlap_count", name, lambda: pa0.count_overlaps(pb0).df[["Chromosome", "Start", "End", "NumberOverlaps"]])
                run_inproc("coverage", name, _cov0)
                run_inproc("subtract", name, lambda: pa0.subtract(pb0).df[["Chromosome", "Start", "End"]])
        except ImportError:
            pass
        ops = "intersect,merge,nearest,overlap_count,coverage,genome_depth,subtract,complement"

        def score(r):
            if r["op"] in refs:
                return score_new(r)
            if r["op"] == "intersect":
                return set_score(ref_int, bed_set(r["out"]))
            if r["op"] == "merge":
                return set_score(ref_mrg, bed_set(r["out"]))
            d = read_tsv(r["out"], ["idx", "bstart", "bend"])
            return gap_score(ref_gap, d["idx"], d["bstart"], d["bend"], A)

        od = WORK / f"intervals_{n}"
        record_worker(rec, "intervals", worker("R", "intervals", od / "r", a_bed, b_bed, ops, CS), n=n, score=score)
        record_worker(rec, "intervals", worker("py", "intervals", od / "py", a_bed, b_bed, ops, CS), n=n, score=score)
        # BEDOPS (its inputs must be sort-bed sorted: sorted once, outside the clock)
        if PYBIN and (PYBIN / "bedops").exists():
            sa, sb = str(od / "A.sortbed.bed"), str(od / "B.sortbed.bed")
            with open(sa, "w") as f:
                subprocess.run([str(PYBIN / "sort-bed"), a_bed], stdout=f, check=True)
            with open(sb, "w") as f:
                subprocess.run([str(PYBIN / "sort-bed"), b_bed], stdout=f, check=True)
            bo = str(PYBIN / "bedops")
            t, runs = time_cli(lambda: out.__setitem__("bi", cli([bo, "--element-of", "1", sa, sb])), rep)
            acc, d = set_score(ref_int, bed_set(pd.read_csv(io.StringIO(out["bi"]), sep="\t", header=None)))
            rec.add(part="intervals", task="intersect", engine="BEDOPS bedops --element-of (CLI)", n=n, seconds=t, runs=runs, accuracy=acc, detail=d)
            t, runs = time_cli(lambda: out.__setitem__("bm", cli([bo, "--merge", sa])), rep)
            acc, d = set_score(ref_mrg, bed_set(pd.read_csv(io.StringIO(out["bm"]), sep="\t", header=None)))
            rec.add(part="intervals", task="merge", engine="BEDOPS bedops --merge (CLI)", n=n, seconds=t, runs=runs, accuracy=acc, detail=d)
            cf = str(PYBIN / "closest-features")
            t, runs = time_cli(lambda: out.__setitem__("bc", cli([cf, "--closest", "--delim", "\t", sa, sb])), rep)
            c = pd.read_csv(io.StringIO(out["bc"]), sep="\t", header=None)
            na = c.notna().all(axis=1)
            nca = sum(1 for _ in open(sa).readline().split("\t"))        # columns of A in the output
            g = gap(c.iloc[:, 1].to_numpy(), c.iloc[:, 2].to_numpy(), pd.to_numeric(c.iloc[:, nca + 1], errors="coerce").to_numpy(),
                    pd.to_numeric(c.iloc[:, nca + 2], errors="coerce").to_numpy())
            # rows of the sorted A map to rows of A by coordinates (duplicates included)
            bo_gap = scatter_by_coords(key, c, np.where(na.to_numpy(), g, -1), n)
            ok = ref_gap >= 0
            same = bo_gap[ok] == ref_gap[ok]
            rec.add(part="intervals", task="nearest", engine="BEDOPS closest-features (CLI)", n=n, seconds=t, runs=runs,
                    accuracy=float(same.mean()), detail=f"{same.sum():,}/{ok.sum():,} rows agree")
            # the 2026-10-23 core jobs: bedmap for counts and covered bases, bedops --difference for subtract and
            # (against the whole genome) complement; no single BEDOPS command gives a depth track
            graw, gbed = str(od / "genome.bed"), str(od / "genome.sortbed.bed")
            Path(graw).write_text("".join(f"{c}\t0\t{e}\n" for c, e in zip(sizes["chrom"], sizes["end"])))
            with open(gbed, "w") as f:
                subprocess.run([str(PYBIN / "sort-bed"), graw], stdout=f, check=True)
            bm = str(PYBIN / "bedmap")
            for task, eng, cmd in (("overlap_count", "BEDOPS bedmap --count (CLI)", [bm, "--echo", "--count", "--delim", "\t", sa, sb]),
                                   ("coverage", "BEDOPS bedmap --bases-uniq (CLI)", [bm, "--echo", "--bases-uniq", "--delim", "\t", sa, sb]),
                                   ("subtract", "BEDOPS bedops --difference (CLI)", [bo, "--difference", sa, sb]),
                                   ("complement", "BEDOPS bedops --difference genome (CLI)", [bo, "--difference", gbed, sa])):
                try:
                    t, runs = time_cli(lambda cmd=cmd: out.__setitem__(task, cli(cmd)), rep)
                    df = pd.read_csv(io.StringIO(out[task]), sep="\t", header=None)
                    if task in ("overlap_count", "coverage"):
                        df = pd.concat([df.iloc[:, :3], df.iloc[:, -1]], axis=1)
                    acc, d = score_new({"op": task}, df)
                    rec.add(part="intervals", task=task, engine=eng, n=n, seconds=t, runs=runs, accuracy=acc, detail=d, n_out=len(df))
                except Exception as e:
                    rec.add(part="intervals", task=task, engine=eng, n=n, seconds=None, runs=None, accuracy=None,
                            error=f"{type(e).__name__}: {e}"[:300])


# ── annotation ───────────────────────────────────────────────────────────────

LABEL_MAP = {  # every tool's vocabulary -> the shared six classes
    "promoter-tss": "Promoter", "promoter": "Promoter", "5' utr": "5UTR", "5utr": "5UTR", "3' utr": "3UTR", "3utr": "3UTR",
    "exon": "Exon", "exonic": "Exon", "1st exon": "Exon", "non-coding": "Exon", "intron": "Intron", "intronic": "Intron",
    "1st intron": "Intron", "intergenic": "Intergenic", "distal intergenic": "Intergenic", "downstream": "Intergenic",
    "tts": "Intergenic",
}


def collapse(label):
    s = str(label).lower()
    for k in sorted(LABEL_MAP, key=len, reverse=True):
        if s.startswith(k):
            return LABEL_MAP[k]
    return "Intergenic" if "intergenic" in s else s


def tss_windows_from_gtf(gtf, r, level="transcript"):
    """TSS ± r windows straight from the GTF (independent of genomeblocks' parser):
    one per transcript, or one per gene (the gene line's start / end by strand)."""
    rows = []
    key = 'transcript_id "' if level == "transcript" else 'gene_id "'
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != level:
                continue
            s, e = int(p[3]) - 1, int(p[4])
            tss = s if p[6] == "+" else e - 1
            tid = p[8].split(key)[1].split('"')[0]
            rows.append((p[0], max(0, tss - r), tss + r + 1, tid, tss))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "tid", "tss"])


def tx_level(engine: str) -> bool:
    """Tools whose nearest TSS / promoter is per transcript, not per gene."""
    return any(k in engine for k in ("ChIPseeker", "ChIPpeakAnno", "HOMER", "UROPA", "transcript-level"))


def part_annotation(rec):
    print("\n== annotation ==")
    gtf = str(DATA / "genes.gtf")
    n = min(100_000, max(SIZES))
    pk = str(DATA / f"peaks_A_{n}.bed")
    P = Loci.make(pk, sort=False)                       # file order = the row index the workers report
    r = 1000
    res = {}
    t = timeit(lambda: res.__setitem__("G", Genes.make(gtf)), repeat=REPEAT)
    G = res["G"]
    rec.add(part="annotation", task="parse_gtf", engine="genomeblocks Genes.make", n=sum(1 for _ in open(gtf)), seconds=t["median"],
            runs=t["runs"], accuracy=None, n_out=len(G.transcripts))
    t = timeit(lambda: res.__setitem__("lab", G.labels(P)), repeat=REPEAT, setup=lambda: setattr(G, "_annot", None))
    from genomeblocks.genes import LABELS
    gb_lab = np.array([collapse(x) for x in LABELS[res["lab"]]])
    t2 = timeit(lambda: res.__setitem__("nt", G.nearest_tss(P)), repeat=REPEAT)
    gb_gap = np.asarray(res["nt"][1], np.int64)
    # independent references: bedtools closest to TSS windows read straight from the GTF, at both
    # levels (one TSS per gene = genomeblocks' definition; one per transcript = ChIPseeker's, HOMER's,
    # UROPA's). Both files sorted the same way, as bedtools requires; rows mapped back by coordinates.
    od = WORK / "annotation"
    od.mkdir(parents=True, exist_ok=True)
    pks = str(od / "peaks.sorted.bed")
    with open(pks, "w") as fh:
        subprocess.run(["sort", "-k1,1", "-k2,2n", pk], stdout=fh, check=True)
    key = coord_rows(P)
    nb = sum(1 for _ in open(pks).readline().split("\t"))
    refs, wins = {}, {}
    for level in ("gene", "transcript"):
        win = tss_windows_from_gtf(gtf, r, level)
        wins[level] = win
        wb = str(od / f"tss_windows_{level}.bed")
        win[["chrom", "start", "end", "tid"]].to_csv(wb + ".unsorted", sep="\t", header=False, index=False)
        with open(wb, "w") as fh:
            subprocess.run(["sort", "-k1,1", "-k2,2n", wb + ".unsorted"], stdout=fh, check=True)
        out = {}
        tb, runs = time_cli(lambda: out.__setitem__("c", cli([tool("bedtools"), "closest", "-t", "first", "-a", pks, "-b", wb])), REPEAT)
        c = pd.read_csv(io.StringIO(out["c"]), sep="\t", header=None)
        hit = (c.iloc[:, nb + 1] >= 0).to_numpy()
        g = np.where(hit, gap(c.iloc[:, 1].to_numpy(), c.iloc[:, 2].to_numpy(), c.iloc[:, nb + 1].to_numpy(), c.iloc[:, nb + 2].to_numpy()), -1)
        refs[level] = scatter_by_coords(key, c, g, n)
        rec.add(part="annotation", task="nearest_tss", engine=f"bedtools closest ({level}-level TSS windows, CLI)", n=n, seconds=tb, runs=runs,
                accuracy=1.0, detail=f"reference ({len(win):,} windows)")
    same = refs["gene"] == gb_gap
    rec.add(part="annotation", task="nearest_tss", engine="genomeblocks nearest_tss", n=n, seconds=t2["median"], runs=t2["runs"],
            accuracy=float(same.mean()), detail=f"{same.sum():,}/{n:,} rows agree with bedtools closest (gene-level TSS)")
    tx_tss = dict(zip(wins["transcript"]["tid"], zip(wins["transcript"]["chrom"], wins["transcript"]["tss"])))

    def gap_from_tss(idx, chrom, tss1, level="transcript"):
        """gap between peak rows and TSS ± r windows (tss1 is 1-based), against the reference at ``level``."""
        ref_gap = refs[level]
        idx = np.asarray(idx, int)
        tss0 = np.asarray(tss1, np.int64) - 1
        g = gap(P.starts[idx], P.ends[idx], np.maximum(0, tss0 - r), tss0 + r + 1)
        g = np.where(np.asarray(chrom).astype(str) == P.chroms[idx].astype(str), g, 10 ** 12)
        ok = ref_gap[idx] >= 0
        same = g[ok] == ref_gap[idx][ok]
        return float(same.mean()), f"{same.sum():,}/{ok.sum():,} rows: same distance to the nearest TSS window ({level}-level)"

    rows = worker("R", "annotation", od / "r", pk, gtf, r, env={"GB_REPEAT_SLOW": "1"})
    ref_lab = {}
    for level in ("gene", "transcript"):
        gf = next((x for x in rows if x["op"] == "labels" and x["engine"].startswith("GenomicFeatures") and f"{level}-level" in x["engine"] and x.get("out")), None)
        if gf:
            d = read_tsv(gf["out"], ["idx", "label"]).sort_values("idx")
            ref_lab[level] = np.array([collapse(x) for x in d["label"]])
    if "gene" in ref_lab:
        same = ref_lab["gene"] == gb_lab
        rec.add(part="annotation", task="labels", engine="genomeblocks labels", n=n, seconds=t["median"], runs=t["runs"],
                accuracy=float(same.mean()), detail=f"{same.sum():,}/{n:,} rows: same class as GenomicFeatures (gene-level TSS windows, same precedence)")
    else:
        rec.add(part="annotation", task="labels", engine="genomeblocks labels", n=n, seconds=t["median"], runs=t["runs"], accuracy=None)

    def score(rr):
        if rr["op"] == "parse_gtf":
            return None, ""
        level = "transcript" if tx_level(rr["engine"]) else "gene"
        if rr["op"] == "labels":
            d = read_tsv(rr["out"], ["idx", "label", "tss1"] if "ChIPseeker" in rr["engine"] else ["idx", "label"])
            lab = np.array([collapse(x) for x in d["label"]])
            if rr["engine"].startswith("GenomicFeatures"):
                return 1.0, f"reference for the classes ({level}-level TSS windows)"
            if level not in ref_lab:
                return None, ""
            idx = d["idx"].to_numpy(int)
            same = lab == ref_lab[level][idx]
            det = f"{same.sum():,}/{len(idx):,} rows: same class as GenomicFeatures ({level}-level TSS windows)"
            if "tss1" in d:
                a2, d2 = gap_from_tss(idx, P.chroms[idx], d["tss1"], level)
                rec.add(part="annotation", task="nearest_tss", engine=rr["engine"], n=n, seconds=float(np.median(rr["seconds"])),
                        runs=rr["seconds"], accuracy=a2, detail=d2)
            return float(same.mean()), det
        if rr["op"] == "nearest_tss":
            if "ChIPpeakAnno" in rr["engine"]:
                d = read_tsv(rr["out"], ["idx", "tss1"])
                return gap_from_tss(d["idx"], P.chroms[d["idx"].to_numpy(int)], d["tss1"], level)
            d = read_tsv(rr["out"], ["idx", "gap"])
            idx = d["idx"].to_numpy(int)
            ref_gap = refs[level]
            ok = ref_gap[idx] >= 0
            same = d["gap"].to_numpy()[ok] == ref_gap[idx][ok]
            return float(same.mean()), f"{same.sum():,}/{ok.sum():,} rows agree ({level}-level)"
        return None, ""
    record_worker(rec, "annotation", rows, n=n, score=score)
    # HOMER annotatePeaks.pl with the GTF (its promoter is -1 kb..+100 bp: classes differ by design)
    ap = HOMER_BIN / "annotatePeaks.pl"
    if ap.exists():
        def _homer():
            out["h"] = cli([str(ap), pk, "none", "-gtf", gtf, "-cpu", "1"])
        try:
            th, runs = time_cli(_homer, 1)
            h = pd.read_csv(io.StringIO(out["h"]), sep="\t")
            h = h.rename(columns={h.columns[0]: "PeakID"})
            key = {f"{c}:{s}-{e}": i for i, (c, s, e) in enumerate(zip(P.chroms, P.starts, P.ends))}
            # HOMER reports 1-based starts; map back by coordinates
            idx = np.array([key.get(f"{c}:{s - 1}-{e}", -1) for c, s, e in zip(h["Chr"], h["Start"], h["End"])])
            okm = idx >= 0
            lab = np.array([collapse(str(x).split(" (")[0]) for x in h["Annotation"].fillna("Intergenic")])
            if "transcript" in ref_lab:
                same = lab[okm] == ref_lab["transcript"][idx[okm]]
                rec.add(part="annotation", task="labels", engine="HOMER annotatePeaks.pl -gtf (CLI)", n=n, seconds=th, runs=runs,
                        accuracy=float(same.mean()), detail=f"{same.sum():,}/{okm.sum():,} rows: same class as GenomicFeatures (transcript-level; HOMER's promoter is -1 kb..+100 bp)")
            pid = h["Nearest PromoterID"].astype(str)
            ch = np.array([tx_tss.get(p, ("", -1))[0] for p in pid])
            ts = np.array([tx_tss.get(p, ("", -1))[1] for p in pid])
            ok2 = okm & (ts >= 0)
            a2, d2 = gap_from_tss(idx[ok2], ch[ok2], ts[ok2] + 1, "transcript")
            rec.add(part="annotation", task="nearest_tss", engine="HOMER annotatePeaks.pl -gtf (CLI)", n=n, seconds=th, runs=runs, accuracy=a2, detail=d2)
        except Exception as e:
            rec.add(part="annotation", task="labels", engine="HOMER annotatePeaks.pl -gtf (CLI)", n=n, seconds=None, runs=None, accuracy=None, error=str(e)[:300])
    # UROPA: nearest transcript start within 1 Mb
    if PYBIN and (PYBIN / "uropa").exists():
        ud = od / "uropa"
        shutil.rmtree(ud, ignore_errors=True)
        ud.mkdir(parents=True)
        def _uropa():
            cli([str(PYBIN / "uropa"), "--bed", pk, "--gtf", gtf, "--feature", "transcript", "--feature_anchor", "start",
                 "--distance", "1000000", "--threads", "1", "--outdir", str(ud), "--prefix", "u"], cwd=str(ud))
        try:
            tu, runs = time_cli(_uropa, 1)
            fh = next(iter(glob.glob(str(ud / "*finalhits.txt"))))
            u = pd.read_csv(fh, sep="\t")
            key = {f"{c}:{s}-{e}": i for i, (c, s, e) in enumerate(zip(P.chroms, P.starts, P.ends))}
            idx = np.array([key.get(f"{c}:{s}-{e}", -1) for c, s, e in zip(u["peak_chr"], u["peak_start"], u["peak_end"])])
            fs, fe = pd.to_numeric(u["feat_start"], errors="coerce"), pd.to_numeric(u["feat_end"], errors="coerce")
            # UROPA writes feat_start 0-based and feat_end as the GTF end: TSS (1-based) = start + 1 on '+', end on '-'
            tss1 = np.where(u["feat_strand"].astype(str) == "+", fs + 1, fe).astype(float)
            ok = (idx >= 0) & np.isfinite(tss1)
            a2, d2 = gap_from_tss(idx[ok], u["peak_chr"].to_numpy()[ok], tss1[ok].astype(int), "transcript")
            rec.add(part="annotation", task="nearest_tss", engine="UROPA (CLI)", n=n, seconds=tu, runs=runs, accuracy=a2, detail=d2)
        except Exception as e:
            rec.add(part="annotation", task="nearest_tss", engine="UROPA (CLI)", n=n, seconds=None, runs=None, accuracy=None, error=str(e)[:300])
    annotation_window_jobs(rec, pk, gtf, n)


def annotation_window_jobs(rec, pk, gtf, n, w=50_000):
    """C5 (added 2026-10-23): every gene whose TSS (one per gene) lies within w of a peak, as (peak, gene) pairs;
    bedtools window -w is the reference, and each tool is given the parameter that gives the same rule (the TSS
    base falls in [start - w, end + w)). genomeblocks only reports the nearest gene, so it is not timed."""
    od = WORK / "annotation"
    tss = tss_windows_from_gtf(gtf, 0, "gene")
    tb = str(od / "tss_gene.bed")
    tss[["chrom", "start", "end", "tid"]].to_csv(tb, sep="\t", header=False, index=False)
    out = {}
    t, runs = time_cli(lambda: out.__setitem__("w", cli([tool("bedtools"), "window", "-a", pk, "-b", tb, "-w", str(w)])), REPEAT)
    c = pd.read_csv(io.StringIO(out["w"]), sep="\t", header=None)
    ref = set(zip(c.iloc[:, 0].astype(str), c.iloc[:, 1].astype(int), c.iloc[:, 2].astype(int), c.iloc[:, -1].astype(str)))
    rec.add(part="annotation", task="genes_in_window", engine=f"bedtools window -w {w // 1000} kb (CLI)", n=n, seconds=t, runs=runs, accuracy=1.0,
            detail=f"reference ({len(ref):,} peak-gene pairs)")

    def pair_score(df):
        got = set(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(int), df.iloc[:, 2].astype(int), df.iloc[:, 3].astype(str)))
        return set_score(ref, got)
    try:
        import pyranges as pr0
        if not pr0.__version__.startswith("1."):
            Adf = pd.read_csv(pk, sep="\t", header=None, usecols=[0, 1, 2], names=["Chromosome", "Start", "End"])
            pa0 = pr0.PyRanges(Adf)
            pt0 = pr0.PyRanges(tss.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End", "tid": "gene"})[["Chromosome", "Start", "End", "gene"]])
            rec_py(rec, "annotation", "genes_in_window", f"pyranges {pr0.__version__} join (slack)", n,
                   lambda: pa0.join(pt0, slack=w).df[["Chromosome", "Start", "End", "gene"]], pair_score)
    except ImportError:
        pass
    sc = lambda r: pair_score(pd.read_csv(r["out"], sep="\t", header=None))
    record_worker(rec, "annotation", worker("R", "genes_window", od / "r_window", pk, tb, w), n=n, score=sc)
    record_worker(rec, "annotation", worker("py", "genes_window", od / "py_window", pk, tb, w), n=n, score=sc)


# ── signal ───────────────────────────────────────────────────────────────────

FLANK, NBINS = 3_000, 200


def signal_loci(n, seed=0):
    from common import read_chromsizes
    cs = read_chromsizes()
    allp = Loci.make(str(DATA / "peaks_A_100000.bed"))
    size = np.array([cs[c] for c in allp.chroms], np.int64)
    c = allp.centers
    allp = allp.take(np.flatnonzero((c >= FLANK) & (c <= size - FLANK)))
    rng = np.random.default_rng(seed)
    return allp.take(np.sort(rng.choice(len(allp), n, replace=False)))


def part_signal(rec):
    print("\n== signal ==")
    bw = str(DATA / "signal_0.bw")
    n = 10_000
    L = signal_loci(n)
    od = WORK / "signal"
    od.mkdir(parents=True, exist_ok=True)
    bed = write_bed(L, od / "loci.bed")
    res = {}
    t = timeit(lambda: res.__setitem__("S", L.signal([bw], n_bins=NBINS, flank=FLANK, verbose=False, progress=False, dtype=np.float64)), repeat=REPEAT)
    S = res["S"][:, 0, :]
    # independent reference: pyBigWig per-base values, binned by the same edge rule (1,000 rows)
    import pyBigWig
    f = pyBigWig.open(bw)
    k = 1000
    ref = np.zeros((k, NBINS))
    for i in range(k):
        c0 = L.centers[i]
        v = np.asarray(f.values(str(L.chroms[i]), int(c0 - FLANK), int(c0 + FLANK), numpy=True))
        edges = (np.arange(NBINS + 1) * (2 * FLANK)) // NBINS
        with np.errstate(all="ignore"):
            ref[i] = np.nan_to_num([np.nanmean(v[edges[b]:edges[b + 1]]) for b in range(NBINS)])
    f.close()
    acc, d = corr_score(ref, S[:k], tol=1e-5)
    rec.add(part="signal", task="region_matrix", engine="genomeblocks signal (pybigtools)", n=n, seconds=t["median"], runs=t["runs"], accuracy=acc,
            detail=d + " vs per-base pyBigWig means")
    rows = worker("R", "signal", od / "r", bed, bw, FLANK, NBINS)

    def score(r):
        if r["op"] != "region_matrix":
            return None, ""
        M = pd.read_csv(r["out"], sep="\t", header=None).to_numpy(float)
        return corr_score(S[:len(M)] if len(M) <= len(S) else S, M[:len(S)], tol=1e-4)
    record_worker(rec, "signal", rows, n=n, score=score)
    # bigWigAverageOverBed: one mean per region over the ±3 kb window
    if PYBIN and (PYBIN / "bigWigAverageOverBed").exists():
        wbed = str(od / "windows.bed")
        with open(wbed, "w") as fh:
            for i, (c, c0) in enumerate(zip(L.chroms, L.centers)):
                fh.write(f"{c}\t{c0 - FLANK}\t{c0 + FLANK}\tw{i}\n")
        # Loci.make sorts by start, and windows around centres are not in the peaks' order: map rows back by name
        W = Loci.make(wbed, keep=["name"])
        wi = np.array([int(x[1:]) for x in W["name"]])
        t = timeit(lambda: res.__setitem__("m", W.signal([bw], n_bins=1, span=True, verbose=False, progress=False, dtype=np.float64)), repeat=REPEAT)
        gm = np.zeros(n); gm[wi] = res["m"][:, 0, 0]
        # the reference: the exact per-base mean over covered bases of every window
        f = pyBigWig.open(bw)
        ref_m = np.zeros(n)
        for i, (c, c0) in enumerate(zip(L.chroms, L.centers)):
            v = np.asarray(f.values(str(c), int(c0 - FLANK), int(c0 + FLANK), numpy=True))
            with np.errstate(all="ignore"):
                ref_m[i] = np.nan_to_num(np.nanmean(v))
        f.close()
        out = {}
        tb, runs = time_cli(lambda: out.__setitem__("b", cli([str(PYBIN / "bigWigAverageOverBed"), bw, wbed, "/dev/stdout"])), REPEAT)
        b = pd.read_csv(io.StringIO(out["b"]), sep="\t", header=None, names=["name", "size", "covered", "sum", "mean0", "mean"])
        order = b["name"].str[1:].astype(int).to_numpy()
        bm = np.zeros(n); bm[order] = b["mean"].to_numpy()
        acc, d = corr_score(ref_m, bm, rtol=1e-4)           # its text output carries ~6 significant digits
        rec.add(part="signal", task="region_mean", engine="UCSC bigWigAverageOverBed (CLI)", n=n, seconds=tb, runs=runs, accuracy=acc,
                detail=d + " vs per-base means (its mean column: covered bases only)")
        acc, d = corr_score(ref_m, gm, rtol=1e-4)
        rec.add(part="signal", task="region_mean", engine="genomeblocks signal (span=True, 1 bin)", n=n, seconds=t["median"], runs=t["runs"], accuracy=acc,
                detail=d + " vs per-base means")
    signal_bam_jobs(rec)
    signal_track_jobs(rec)


CLIBIN = Path(os.environ["GB_LANDSCAPE_CLI"]) if os.environ.get("GB_LANDSCAPE_CLI") else None    # gb_bench_cli/bin


def clitool(name):
    """A command from the CLI competitors' environment (or PATH); None when absent."""
    if CLIBIN and (CLIBIN / name).exists():
        return str(CLIBIN / name)
    return shutil.which(name)


def rec_cli(rec, part, task, engine, n, fn, score, rep=REPEAT, warmup=1, **extra):
    """Time a CLI run, score its answer, record the row (an error is recorded, not raised). warmup=0 for engines
    that take most of an hour a run."""
    try:
        t, runs = time_cli(fn, rep, warmup)
        acc, d = score()
        rec.add(part=part, task=task, engine=engine, n=n, seconds=t, runs=runs, accuracy=acc, detail=d, **extra)
    except Exception as e:
        rec.add(part=part, task=task, engine=engine, n=n, seconds=None, runs=None, accuracy=None, error=f"{type(e).__name__}: {e}"[:300], **extra)


def rec_py(rec, part, task, engine, n, fn, score, rep=REPEAT, **extra):
    """Time an in-process call; ``score(answer)`` -> (accuracy, detail)."""
    res = {}
    try:
        t = timeit(lambda: res.__setitem__("a", fn()), repeat=rep)
        acc, d = score(res["a"])
        rec.add(part=part, task=task, engine=engine, n=n, seconds=t["median"], runs=t["runs"], accuracy=acc, detail=d, **extra)
    except Exception as e:
        rec.add(part=part, task=task, engine=engine, n=n, seconds=None, runs=None, accuracy=None, error=f"{type(e).__name__}: {e}"[:300], **extra)


BW_RECORDS = [int(x) for x in os.environ.get("GB_BENCH_BW_RECORDS", "1000000,10000000").split(",")]


def track_bedgraph(od, nrec, lex):
    """C9's input: a sorted bedGraph of about ``nrec`` records over the genome (seeded by ``nrec``), written to
    ``od`` once; returns its path and the records as written (values rounded as in the file)."""
    bg = od / f"track_{nrec}.bg"
    rng = np.random.default_rng(nrec)
    per = np.maximum(1, (np.array([s for _, s in lex]) / sum(s for _, s in lex) * nrec).astype(int))
    frames = []
    for (c, s), k in zip(lex, per):
        gaps = rng.integers(1, max(2, int(s // k) - 20), k)
        st = np.cumsum(gaps) + np.arange(k) * 10
        st = st[st + 10 < s]
        frames.append(pd.DataFrame({"c": c, "s": st, "e": st + rng.integers(1, 10, len(st)), "v": rng.random(len(st)).astype(np.float32)}))
    B = pd.concat(frames, ignore_index=True)
    if not bg.exists():
        B.to_csv(bg, sep="\t", header=False, index=False, float_format="%.6g")
    return bg, pd.read_csv(bg, sep="\t", header=None, names=["c", "s", "e", "v"])


def bw_windows(chroms, sizes, n=400, width=2000, seed=5):
    """Sample windows for checking a coverage track: half around simulated peak summits, half anywhere."""
    rng = np.random.default_rng(seed)
    pk = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None, usecols=[0, 1, 2]).sample(n // 2, random_state=seed)
    code = {c: i for i, c in enumerate(chroms)}
    out = [(code[c], max(0, (s + e) // 2 - width // 2)) for c, s, e in pk.itertuples(index=False)]
    ci = rng.integers(0, len(chroms), n - n // 2)
    out += [(int(c), int(rng.integers(0, sizes[c] - width))) for c in ci]
    return [(c, a, a + width) for c, a in out]


def coverage_score(path, truth, chroms, wins, bin_size=1):
    """A coverage bigWig against the simulated reads over the sample windows: per-base depth (bin_size 1), or
    reads overlapping each genome-aligned bin (deepTools' value for a bin wider than one base)."""
    import pyBigWig
    f = pyBigWig.open(str(path))
    ref, got = [], []
    for c, a, b in wins:
        st = truth.get(c, np.zeros(0, np.int64))
        a -= a % bin_size
        pos = np.arange(a, b - bin_size + 1, bin_size)
        if bin_size == 1:
            ref.append(np.searchsorted(st, pos, side="right") - np.searchsorted(st, pos - READ_LEN, side="right"))
        else:
            ref.append(reads_over(st, pos, pos + bin_size))
        v = np.asarray(f.values(chroms[c], int(a), int(pos[-1] + bin_size), numpy=True), float) if chroms[c] in f.chroms() else np.zeros(pos[-1] + bin_size - a)
        got.append(np.nan_to_num(v)[pos - a])
    f.close()
    return corr_score(np.concatenate(ref), np.concatenate(got), tol=1e-6)


def signal_bam_jobs(rec):
    """The signal jobs on BAM files added 2026-10-23: BAM -> per-base coverage bigWig (C1), read counts per
    region across BAMs (C4), TMM normalisation across samples (C8). genomeblocks has no general BAM -> bigWig and
    no multi-BAM counter, so it is timed only for C8 (signal.tmm)."""
    chroms, sizes = chrom_sizes()
    od = WORK / "signal_bam"
    od.mkdir(parents=True, exist_ok=True)
    gsize = str(od / "genome.lex.txt")
    with open(gsize, "w") as f:
        subprocess.run(["sort", "-k1,1", CS], stdout=f, check=True)
    regions = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    rbed, rsaf = str(od / "regions.bed"), str(od / "regions.saf")
    regions.assign(name=[f"r{i}" for i in range(len(regions))]).to_csv(rbed, sep="\t", header=False, index=False)
    pd.DataFrame({"GeneID": [f"r{i}" for i in range(len(regions))], "Chr": regions["chrom"], "Start": regions["start"] + 1,
                  "End": regions["end"], "Strand": "+"}).to_csv(rsaf, sep="\t", index=False)
    code = {c: i for i, c in enumerate(chroms)}
    rc = regions["chrom"].map(code).to_numpy()
    wins = bw_windows(chroms, sizes)
    for nr in BAM_READS:
        bams = chip_bams(nr)
        truths = [read_truth(str(b) + ".npz") for b in bams]
        rep = REPEAT if nr <= 1_000_000 else max(1, REPEAT - 1)
        # C1: BAM -> per-base coverage bigWig (raw counts, no extension)
        b0, t0 = str(bams[0]), truths[0]
        dt, mg, bg2bw = shutil.which("bamCoverage"), clitool("megadepth"), clitool("bedGraphToBigWig")
        if dt:
            # per base (the same output as megadepth and genomecov), and at deepTools' default 50 bp bins
            for bs in (1, 50):
                out = str(od / f"deeptools_{nr}_{bs}.bw")
                rec_cli(rec, "signal", "bam_to_bigwig", f"deepTools bamCoverage --binSize {bs} (CLI)", nr,
                        lambda bs=bs, out=out: cli([dt, "-b", b0, "-o", out, "--binSize", str(bs), "--normalizeUsing", "None", "-p", "1"]),
                        lambda bs=bs, out=out: coverage_score(out, t0, chroms, wins, bs), 1 if bs == 1 else rep,
                        warmup=0 if bs == 1 else 1)          # per base: 42 min a run at 200k reads (2026-10-23 smoke test)
        if mg:
            pre = str(od / f"megadepth_{nr}")
            rec_cli(rec, "signal", "bam_to_bigwig", "megadepth --bigwig (CLI)", nr,
                    lambda: cli([mg, b0, "--bigwig", "--prefix", pre, "--threads", "1"]),
                    lambda: coverage_score(pre + ".all.bw", t0, chroms, wins), rep)
        if bg2bw:
            bg, out = str(od / f"genomecov_{nr}.bg"), str(od / f"genomecov_{nr}.bw")

            def _gc():
                with open(bg, "w") as fh:
                    p1 = subprocess.Popen([tool("bedtools"), "genomecov", "-bg", "-ibam", b0], stdout=subprocess.PIPE)
                    subprocess.run(["sort", "-k1,1", "-k2,2n"], stdin=p1.stdout, stdout=fh, check=True, env={**os.environ, "LC_ALL": "C"})
                    p1.wait()
                cli([bg2bw, bg, gsize, out])
            rec_cli(rec, "signal", "bam_to_bigwig", "bedtools genomecov -bg + UCSC bedGraphToBigWig (CLI)", nr, _gc,
                    lambda: coverage_score(out, t0, chroms, wins), rep)
        # C4: reads overlapping each of 100k regions, in each of the 4 BAMs (the truth: every read overlapping by >= 1 base)
        truth_counts = np.zeros((len(regions), len(bams)), np.int64)
        for j, tr in enumerate(truths):
            for c in np.unique(rc):
                m = rc == c
                truth_counts[m, j] = reads_over(tr.get(int(c), np.zeros(0, np.int64)), regions["start"].to_numpy()[m], regions["end"].to_numpy()[m])
        tkey = dict(zip(zip(regions["chrom"], regions["start"], regions["end"]), map(tuple, truth_counts)))

        def counts_score(df):
            """df: chrom, start, end, one column per BAM (rows in any order)."""
            got = {(str(r[0]), int(r[1]), int(r[2])): tuple(int(v) for v in r[3:]) for r in df.itertuples(index=False)}
            same = sum(1 for k, v in tkey.items() if got.get(k) == v)
            return (1.0 if same == len(tkey) else same / len(tkey)), f"{same:,}/{len(tkey):,} regions agree in all {len(bams)} BAMs"
        bstr = [str(b) for b in bams]
        mbs = shutil.which("multiBamSummary")
        if mbs:
            tab = str(od / f"multibam_{nr}.tab")
            rec_cli(rec, "signal", "multi_bam_counts", "deepTools multiBamSummary BED-file (CLI)", nr,
                    lambda: cli([mbs, "BED-file", "--BED", rbed, "--bamfiles", *bstr, "-o", str(od / "mbs.npz"), "--outRawCounts", tab, "-p", "1"]),
                    lambda: counts_score(pd.read_csv(tab, sep="\t", comment=None, skiprows=1, header=None)), 1, warmup=0)   # 57 min a run at 200k reads
        out = {}
        rec_cli(rec, "signal", "multi_bam_counts", "bedtools multicov (CLI)", nr,
                lambda: out.__setitem__("mc", cli([tool("bedtools"), "multicov", "-bams", *bstr, "-bed", rbed])),
                lambda: counts_score(pd.read_csv(io.StringIO(out["mc"]), sep="\t", header=None).drop(columns=[3])), rep)
        fc = clitool("featureCounts")
        if fc:
            fco = str(od / f"featurecounts_{nr}.txt")

            def _fc_score():
                d = pd.read_csv(fco, sep="\t", comment="#")
                return counts_score(pd.concat([d[["Chr", "Start", "End"]].assign(Start=d["Start"] - 1), d.iloc[:, 6:]], axis=1))
            rec_cli(rec, "signal", "multi_bam_counts", "Subread featureCounts -F SAF -O (CLI)", nr,
                    lambda: cli([fc, "-F", "SAF", "-a", rsaf, "-O", "-T", "1", "-o", fco, *bstr]), _fc_score, rep)
        try:
            import pysam
            hs = [pysam.AlignmentFile(b) for b in bstr]

            def _pysam():
                return np.array([[h.count(c, s, e) for h in hs] for c, s, e in zip(regions["chrom"], regions["start"], regions["end"])])
            rec_py(rec, "signal", "multi_bam_counts", f"pysam {pysam.__version__} AlignmentFile.count", nr, _pysam,
                   lambda a: counts_score(pd.concat([regions, pd.DataFrame(a)], axis=1)), rep=1)
        except ImportError:
            pass
        rows = worker("R", "bam_counts", od / f"r_{nr}", rsaf, *bstr)
        record_worker(rec, "signal", rows, n=nr, score=lambda r: counts_score(pd.read_csv(r["out"], sep="\t", header=None)))
        # C8: TMM normalisation across the 4 samples, on the true region counts; the answer is the normalised
        # matrix (counts per million after TMM), edgeR the reference
        cm = str(od / f"region_counts_{nr}.tsv")
        pd.DataFrame(truth_counts, columns=[f"s{j}" for j in range(len(bams))]).to_csv(cm, sep="\t", index=False)
        rows = worker("R", "track_tmm", od / f"tmm_{nr}", cm) + worker("py", "track_tmm", od / f"tmm_py_{nr}", cm)
        ref_path = next((r["out"] for r in rows if r["engine"] == "edgeR calcNormFactors + cpm" and not r.get("error")), None)
        ref = pd.read_csv(ref_path, sep="\t", header=None).to_numpy(float) if ref_path else None

        def tmm_score(got):
            if ref is None:
                return None, "no edgeR reference"
            return corr_score(ref, np.asarray(got, float), rtol=1e-6)
        record_worker(rec, "signal", rows, n=nr,
                      score=lambda r: (1.0, "reference") if r["out"] == ref_path else tmm_score(pd.read_csv(r["out"], sep="\t", header=None)))
        from genomeblocks import signal as gsig
        cube = truth_counts.astype(float)[:, :, None]
        rec_py(rec, "signal", "track_tmm", "genomeblocks signal.tmm", nr, lambda: gsig.tmm(cube)[:, :, 0], tmm_score)


def signal_track_jobs(rec):
    """The signal jobs on tracks added 2026-10-23: log2 comparison of two tracks (C7) and writing a bigWig from
    a bedGraph (C9). genomeblocks has neither as a function of its own, so it is not timed."""
    import pyBigWig
    od = WORK / "signal_tracks"
    od.mkdir(parents=True, exist_ok=True)
    # C7 on chr1 copies of two synthetic tracks (both tools read the same files): log2((mean a + 1) / (mean b + 1))
    # per 50 bp genome-aligned bin, deepTools' default bin (per base, bigwigCompare did not finish one run on chr1 in
    # 90 min on 2026-10-24); WiggleTools' bin operator sums, so it is scaled by 1/50 to the mean
    B_ = 50
    src = [str(DATA / "signal_0.bw"), str(DATA / "signal_1.bw")]
    chr1 = [str(od / f"chr1_{j}.bw") for j in range(2)]
    if not all(Path(p).exists() for p in chr1):
        for s_, d_ in zip(src, chr1):
            fi, fo = pyBigWig.open(s_), pyBigWig.open(d_, "w")
            fo.addHeader([("chr1", fi.chroms()["chr1"])])
            iv = fi.intervals("chr1")
            fo.addEntries(["chr1"] * len(iv), [a for a, _, _ in iv], ends=[b for _, b, _ in iv], values=[float(v) for _, _, v in iv])
            fi.close(); fo.close()
    rng = np.random.default_rng(9)
    L1 = 248_956_422
    starts = np.sort(rng.integers(0, (L1 - 4000) // B_, 300)) * B_       # 300 windows of 40 bins, bin-aligned
    bins_at = np.concatenate([a + np.arange(40) * B_ for a in starts])

    def bin_means(path):
        f = pyBigWig.open(path)
        out = np.concatenate([np.nan_to_num(np.asarray(f.values("chr1", int(a), int(a) + 40 * B_, numpy=True))).reshape(40, B_).mean(1) for a in starts])
        f.close()
        return out
    ref = np.log2((bin_means(chr1[0]) + 1) / (bin_means(chr1[1]) + 1))

    def ratio_score(path, bedgraph=False):
        if bedgraph:
            d = pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e", "v"], comment="#")
            d = d[d["c"] == "chr1"].sort_values("s")
            s0, e0, v0 = d["s"].to_numpy(), d["e"].to_numpy(), d["v"].to_numpy(float)
            k = np.searchsorted(s0, bins_at, side="right") - 1
            inside = (k >= 0) & (bins_at < e0[np.maximum(k, 0)])
            got = np.where(inside, v0[np.maximum(k, 0)], 0.0)
        else:
            f = pyBigWig.open(path)
            got = np.concatenate([np.nan_to_num(np.asarray(f.values("chr1", int(a), int(a) + 40 * B_, numpy=True)))[::B_] for a in starts])
            f.close()
        return corr_score(ref, got, tol=1e-4)
    bwc = shutil.which("bigwigCompare")
    if bwc:
        out = str(od / "bigwigcompare.bw")
        rec_cli(rec, "signal", "track_comparison", f"deepTools bigwigCompare --binSize {B_} (CLI)", L1,
                lambda: cli([bwc, "-b1", chr1[0], "-b2", chr1[1], "--operation", "log2", "--pseudocount", "1", "--binSize", str(B_), "-p", "1", "-o", out]),
                lambda: ratio_score(out), max(1, REPEAT - 1))
    wt = clitool("wiggletools")
    if wt:
        out = str(od / "wiggletools.bg")
        sc_ = str(1 / B_)
        rec_cli(rec, "signal", "track_comparison", f"WiggleTools log 2 ratio of {B_} bp bin means (CLI)", L1,
                lambda: (Path(out).unlink(missing_ok=True), cli([wt, "write_bg", out, "log", "2", "ratio", "offset", "1", "scale", sc_, "bin", str(B_), "default", "0", chr1[0],
                             "offset", "1", "scale", sc_, "bin", str(B_), "default", "0", chr1[1]])),         # it will not overwrite
                lambda: ratio_score(out, bedgraph=True), max(1, REPEAT - 1))
    # C9: write a bigWig from a sorted bedGraph of n records (values read back over sample windows)
    chroms, sizes = chrom_sizes()
    lex = sorted(zip(chroms, sizes))
    for nrec in BW_RECORDS:
        bg, B = track_bedgraph(od, nrec, lex)
        n_real = len(B)
        sample = B.sample(5000, random_state=1)

        def wscore(path):
            f = pyBigWig.open(str(path))
            got = np.array([f.values(c, int(s), int(s) + 1)[0] for c, s in zip(sample["c"], sample["s"])], float)
            f.close()
            return corr_score(sample["v"].to_numpy(float), got, rtol=1e-6)
        header = [(c, int(s)) for c, s in lex]
        arrs = {c: (g["s"].to_numpy(np.int64), g["e"].to_numpy(np.int64), g["v"].to_numpy(np.float64)) for c, g in B.groupby("c", sort=False)}
        out = od / f"pybigwig_{nrec}.bw"

        def _pbw():
            f = pyBigWig.open(str(out), "w")
            f.addHeader(header)
            for c, _ in header:
                if c in arrs:
                    s, e, v = arrs[c]
                    f.addEntries([c] * len(s), s.tolist(), ends=e.tolist(), values=v.tolist())
            f.close()
            return out
        rec_py(rec, "signal", "write_bigwig", f"pyBigWig {pyBigWig.__version__} addEntries", n_real, _pbw, wscore, rep=max(1, REPEAT - 1))
        try:
            import pybigtools
            out2 = od / f"pybigtools_{nrec}.bw"

            def _pbt():
                f = pybigtools.open(str(out2), "w")                  # write() closes the file
                f.write(dict(header), ((c, int(s), int(e), float(v)) for c, s, e, v in zip(B["c"], B["s"], B["e"], B["v"])))
                return out2
            rec_py(rec, "signal", "write_bigwig", f"pybigtools {getattr(pybigtools, '__version__', '')} write".replace("  ", " "), n_real, _pbt, wscore, rep=max(1, REPEAT - 1))
        except ImportError:
            pass
        bg2bw = clitool("bedGraphToBigWig")
        gsize = od / "genome.lex.txt"
        gsize.write_text("".join(f"{c}\t{s}\n" for c, s in lex))
        if bg2bw:
            out3 = str(od / f"bg2bw_{nrec}.bw")
            rec_cli(rec, "signal", "write_bigwig", "UCSC bedGraphToBigWig (CLI)", n_real, lambda: cli([bg2bw, str(bg), str(gsize), out3]),
                    lambda: wscore(out3), max(1, REPEAT - 1))
        rows = worker("R", "write_bigwig", od / f"r_{nrec}", str(bg), str(gsize))
        record_worker(rec, "signal", rows, n=n_real, score=lambda r: wscore(r["out"]))


def part_signal_bigwig(rec):
    """genomeblocks' bigWig writers (2026-11) on the inputs the competitors of C1 and C9 ran on: the seeded ChIP
    BAMs (chip_bams; per base and 50 bp bins, one process) scored against the simulated reads, and the seeded
    bedGraph tracks (track_bedgraph) written from memory and read back over 5,000 sample positions. A part of
    its own so it runs without the competitors' environments and keeps their rows."""
    print("\n== signal_bigwig ==")
    chroms, sizes = chrom_sizes()
    od = WORK / "signal_bigwig"
    od.mkdir(parents=True, exist_ok=True)
    wins = bw_windows(chroms, sizes)
    for nr in BAM_READS:
        b0 = str(chip_bams(nr)[0])
        t0 = read_truth(b0 + ".npz")
        rep = REPEAT if nr <= 1_000_000 else max(1, REPEAT - 1)
        for bs in (1, 50):
            out = str(od / f"bam_{nr}_{bs}.bw")
            rec_py(rec, "signal_bigwig", "bam_to_bigwig", f"genomeblocks bam_to_bigwig(bin_size={bs})", nr,
                   lambda bs=bs, out=out: gb.bam_to_bigwig(b0, out, bin_size=bs, workers=1),
                   lambda p, bs=bs: coverage_score(p, t0, chroms, wins, bs), rep=rep)
    lex = sorted(zip(chroms, sizes))
    td = WORK / "signal_tracks"                       # the competitors' bedGraphs, regenerated there if missing
    td.mkdir(parents=True, exist_ok=True)
    for nrec in BW_RECORDS:
        bg, B = track_bedgraph(td, nrec, lex)
        sample = B.sample(5000, random_state=1)

        def wscore(path, sample=sample):
            import pyBigWig
            f = pyBigWig.open(str(path))
            got = np.array([f.values(c, int(s), int(s) + 1)[0] for c, s in zip(sample["c"], sample["s"])], float)
            f.close()
            return corr_score(sample["v"].to_numpy(float), got, rtol=1e-6)
        L = Loci.from_frame(B.rename(columns={"c": "chrom", "s": "start", "e": "end"}))
        out = str(od / f"values_{nrec}.bw")
        rec_py(rec, "signal_bigwig", "write_bigwig", "genomeblocks Loci.to_bigwig", len(B),
               lambda L=L, out=out: L.to_bigwig(out, "v", chrom_sizes=dict(lex)), wscore, rep=max(1, REPEAT - 1))


# ── motifs ───────────────────────────────────────────────────────────────────

def part_motifs(rec):
    print("\n== motifs ==")
    import bench_motifs as bm
    from genomeblocks.motifs import scan_motifs_matrix
    genome = str(DATA / "genome.fa")
    gdict = gb.read_fasta(genome)                    # {chrom: sequence}: what bench_motifs' helpers take
    N, M = 1000, 100
    THR, PVAL = 13.0, 1e-4
    L = bm.windows(gdict, N)
    seqs = L.sequences(gdict, r=bm.R, upper=True)
    path = bm.subset_file(M)
    od = WORK / "motifs"
    od.mkdir(parents=True, exist_ok=True)
    fa = str(od / "seqs.fa")
    bm.write_fasta(seqs, fa)
    lib = gb.load_motifs(path, format="jaspar16")
    pssm = str(od / "pssm.txt")
    with open(pssm, "w") as f:
        for i in range(len(lib)):
            lo = np.asarray(lib.logodds(i), float)
            lo = lo.T if lo.shape[1] == 4 else lo               # (4, W), rows A C G T
            f.write(f">{lib.names[i]}\n" + "\n".join(" ".join(f"{x:.6f}" for x in row) for row in lo) + "\n")
    # independent references (2026-10-23; genomeblocks' own scan used to be the reference): an exact numpy
    # sliding-window scan of the same log-odds matrices for the log-odds task, MEME FIMO (same probabilities, no
    # pseudocount, uniform background, both strands) for the p-value task
    mats = []
    for i in range(len(lib)):
        lo = np.asarray(lib.logodds(i), float)
        mats.append(lo if lo.shape[1] == 4 else lo.T)          # (W, 4), columns A C G T
    code = np.zeros(256, np.uint8)
    for i, b_ in enumerate(b"ACGT"):
        code[b_] = i
    X = np.stack([code[np.frombuffer(sq.encode(), np.uint8)] for sq in seqs])
    ref_lo = {}
    for nm, lo in zip(lib.names, mats):
        win = np.lib.stride_tricks.sliding_window_view(X, lo.shape[0], axis=1)
        ref_lo[str(nm)] = int((lo[np.arange(lo.shape[0]), win].sum(axis=-1) >= THR).sum())
    rec.add(part="motifs", task="scan_logodds", engine="exact numpy sliding-window scan (reference)", n=N, n_motifs=M, seconds=None, runs=None,
            accuracy=1.0, detail=f"reference: {int(sum(ref_lo.values())):,} hits at log2-odds ≥ {THR}")
    # p-value reference: each motif's exact score distribution under the uniform background (scores rounded to 1e-4
    # bits), the smallest score with P(S >= s) <= p, and a scan of both strands with the same rounded scores. A motif
    # too short to reach p (6 bp: the rarest score has p = 4^-6 = 2.4e-4) has no hits
    ref_p = {}
    Xr = 3 - X[:, ::-1]
    for nm, lo in zip(lib.names, mats):
        q, thr = exact_pvalue_threshold(lo, PVAL)
        hits = 0
        if np.isfinite(thr):
            for XX in (X, Xr):
                win = np.lib.stride_tricks.sliding_window_view(XX, len(q), axis=1)
                hits += int((q[np.arange(len(q)), win].sum(axis=-1) >= thr).sum())
        ref_p[str(nm)] = hits
    rec.add(part="motifs", task="scan_pvalue", engine="exact score distribution + numpy scan (reference)", n=N, n_motifs=M, seconds=None, runs=None,
            accuracy=1.0, detail=f"reference: {int(sum(ref_p.values())):,} hits at p ≤ {PVAL}, both strands")
    meme = str(od / "motifs.meme")
    write_meme_probs(lib.names, mats, meme)
    fimo = tool("fimo")
    out = {}
    tf, runs = time_cli(lambda: out.__setitem__("f", cli([fimo, "--text", "--thresh", str(PVAL), "--bfile", "--uniform--", "--motif-pseudo", "0",
                                                          "--skip-matched-sequence", "--verbosity", "1", meme, fa])), REPEAT)
    fimo_p = {str(nm): 0 for nm in lib.names}
    for line in out["f"].splitlines():
        if line and not line.startswith(("motif_id", "#")):
            fimo_p[line.split("\t")[0]] = fimo_p.get(line.split("\t")[0], 0) + 1
    res = {}
    t = timeit(lambda: res.__setitem__("a", scan_motifs_matrix(L, genome, path, format="jaspar16", r=bm.R, threshold=THR, norm=False, workers=1, verbose=False)), repeat=REPEAT)
    gb_lo = {str(k): int(v) for k, v in res["a"].sum(axis=0).to_dict().items()}
    # the p-value task scans both strands: motifmatchr always does
    t2 = timeit(lambda: res.__setitem__("p", scan_motifs_matrix(L, genome, path, format="jaspar16", r=bm.R, pvalue=PVAL, norm=False, both_strands=True,
                                                                workers=1, verbose=False)), repeat=REPEAT)
    gb_p = {str(k): int(v) for k, v in res["p"].sum(axis=0).to_dict().items()}

    def hits_score(ref, got):
        keys = sorted(ref)
        a = np.array([ref[k] for k in keys], float)
        b = np.array([got.get(k, 0) for k in keys], float)
        if np.array_equal(a, b):
            return 1.0, f"identical: {int(b.sum()):,} hits"
        r = float(np.corrcoef(a, b)[0, 1]) if a.std() and b.std() else 0.0
        return max(0.0, r), f"{int(b.sum()):,} vs {int(a.sum()):,} hits, Pearson {r:.4f} over motifs"

    acc, d = hits_score(ref_lo, gb_lo)
    eng = gb.backends.resolve("motifs")                     # the backend genomeblocks picks by itself (fastest installed)
    rec.add(part="motifs", task="scan_logodds", engine=f"genomeblocks scan_motifs_matrix ({eng})", n=N, n_motifs=M, seconds=t["median"], runs=t["runs"],
            accuracy=acc, detail=d + " vs the exact numpy scan")
    acc, d = hits_score(ref_p, gb_p)
    rec.add(part="motifs", task="scan_pvalue", engine=f"genomeblocks scan_motifs_matrix ({eng}, both strands)", n=N, n_motifs=M, seconds=t2["median"],
            runs=t2["runs"], accuracy=acc, detail=d + " vs the exact p-value scan")
    acc, d = hits_score(ref_p, fimo_p)
    rec.add(part="motifs", task="scan_pvalue", engine="MEME FIMO --text, both strands (CLI)", n=N, n_motifs=M, seconds=tf, runs=runs,
            accuracy=acc, detail=d + " vs the exact p-value scan")

    def score(r):
        d = read_tsv(r["out"], ["motif", "hits"])
        got = dict(zip(d["motif"].astype(str), d["hits"]))
        return hits_score(ref_lo if r["op"] == "scan_logodds" else ref_p, got)
    rows = worker("R", "motifs", od / "r", fa, pssm, THR, PVAL, env={"GB_REPEAT_SLOW": "1"})
    record_worker(rec, "motifs", rows, n=N, score=score, n_motifs=M)
    # HOMER: the same PSSMs as probabilities, its natural-log threshold = THR * ln 2
    h2 = HOMER_BIN / "homer2"
    if h2.exists():
        mf = str(od / "motifs.homer")
        with open(mf, "w") as f:
            for i in range(len(lib)):
                lo = np.asarray(lib.logodds(i), float)
                lo = lo if lo.shape[1] == 4 else lo.T        # (W, 4)
                p = 0.25 * np.power(2.0, lo)
                p = p / p.sum(axis=1, keepdims=True)
                cons = "".join("ACGT"[j] for j in p.argmax(axis=1))
                f.write(f">{cons}\t{lib.names[i]}\t{THR * np.log(2):.6f}\n")
                f.write("\n".join("\t".join(f"{x:.6f}" for x in row) for row in p) + "\n")
        out = {}
        def _homer():
            out["h"] = cli([str(h2), "find", "-i", fa, "-m", mf, "-strand", "+", "-p", "1"])
        try:
            th, runs = time_cli(_homer, REPEAT)
            cnt = {}
            for line in out["h"].splitlines():
                p = line.split("\t")
                if len(p) >= 4:
                    cnt[p[3]] = cnt.get(p[3], 0) + 1
            acc, d = hits_score(ref_lo, cnt)
            rec.add(part="motifs", task="scan_logodds", engine="HOMER homer2 find (CLI)", n=N, n_motifs=M, seconds=th, runs=runs, accuracy=acc, detail=d)
        except Exception as e:
            rec.add(part="motifs", task="scan_logodds", engine="HOMER homer2 find (CLI)", n=N, n_motifs=M, seconds=None, runs=None, accuracy=None, error=str(e)[:300])

    motif_new_jobs(rec, lib, mats, od, THR, PVAL, path)



def exact_pvalue_threshold(lo, p, step=1e-4):
    """(integer score matrix, integer threshold) for P(S >= threshold) <= p under a uniform background, from the exact
    distribution of the motif's scores rounded to ``step``; threshold = inf when no score is that rare."""
    q = np.round(np.asarray(lo, float) / step).astype(np.int64)
    dist, base = np.array([1.0]), 0
    for row in q:
        sh = row - row.min()
        new = np.zeros(len(dist) + int(sh.max()))
        for k in sh:
            new[k:k + len(dist)] += dist * 0.25
        dist, base = new, base + int(row.min())
    tail = np.cumsum(dist[::-1])[::-1]                  # P(S >= base + i)
    ok = np.flatnonzero(tail <= p * (1 + 1e-9))
    return q, (base + int(ok[0]) if len(ok) else np.inf)

def write_meme_probs(names, mats, path):
    """MEME-format motifs whose probabilities give back the log2-odds matrices exactly (uniform background)."""
    with open(path, "w") as f:
        f.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for nm, lo in zip(names, mats):
            pr = 0.25 * np.power(2.0, np.asarray(lo, float))
            pr = pr / pr.sum(axis=1, keepdims=True)
            f.write(f"MOTIF {nm}\nletter-probability matrix: alength= 4 w= {len(pr)} nsites= 20 E= 0\n")
            f.write("\n".join(" ".join(f"{x:.8f}" for x in row) for row in pr) + "\n\n")


def read_meme_probs(path):
    """{name: (W, 4) probabilities} from a MEME-format file (STREME / MEME output)."""
    out, name, rows = {}, None, None
    for line in open(path):
        t = line.split()
        if line.startswith("MOTIF"):
            name = t[1]
        elif line.startswith("letter-probability"):
            rows = []
        elif rows is not None:
            if len(t) == 4 and all(x.replace(".", "", 1).replace("e-", "", 1).isdigit() for x in t):
                rows.append([float(x) for x in t])
            elif rows:
                out[name] = np.array(rows)
                rows = None
    if rows:
        out[name] = np.array(rows)
    return out


def read_homer_motifs(path):
    """{name: (W, 4) probabilities} from a HOMER motif file."""
    out, name, rows = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None and rows:
                out[name] = np.array(rows)
            name, rows = line[1:].split("\t")[1] if "\t" in line else line[1:].strip(), []
        elif line.strip():
            rows.append([float(x) for x in line.split()[:4]])
    if name is not None and rows:
        out[name] = np.array(rows)
    return out


def pwm_similarity(a, b, min_overlap=5):
    """Best alignment of two probability matrices over offsets and both strands: the summed per-column Pearson r
    of the overlap divided by the longer motif's width (1 = identical; partial overlaps count less)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    best = 0.0
    for bb in (b, b[::-1, ::-1]):                        # forward, reverse complement (ACGT columns reverse to TGCA)
        for off in range(-(len(bb) - min_overlap), len(a) - min_overlap + 1):
            i0, j0 = max(0, off), max(0, -off)
            k = min(len(a) - i0, len(bb) - j0)
            if k < min_overlap:
                continue
            x, y = a[i0:i0 + k] - 0.25, bb[j0:j0 + k] - 0.25
            num = (x * y).sum(axis=1)
            den = np.sqrt((x * x).sum(axis=1) * (y * y).sum(axis=1)) + 1e-12
            best = max(best, float((num / den).sum() / max(len(a), len(bb))))
    return best


def motif_sets(gdict, names, probs, n=1000, width=200, seed=13):
    """Planted-motif FASTA sets (one sequence per record, also usable as a 'genome' of contigs): T has motif X within
    10 bp of the centre in half its sequences and Z anywhere in 40 %; B has Y in half and Z in 40 %; N has none.
    X, Y, Z are high-information motifs of the library that do not resemble each other. Background: random windows
    of the synthetic genome."""
    d = WORK / "motif_sets"
    meta = d / "planted.json"
    if meta.exists():
        m = json.loads(meta.read_text())
        return {k: d / f"{k}.fa" for k in "TBN"}, m
    d.mkdir(parents=True, exist_ok=True)
    ic = np.array([float((2 + (pr * np.log2(np.maximum(pr, 1e-12))).sum(axis=1)).sum()) if 8 <= len(pr) <= 20 else 0.0 for pr in probs])
    pick = []
    for i in np.argsort(-ic):
        if all(pwm_similarity(probs[i], probs[j]) < 0.35 for j in pick):
            pick.append(int(i))
        if len(pick) == 3:
            break
    X, Y, Z = pick
    rng = np.random.default_rng(seed)
    chroms = [c for c in gdict if len(gdict[c]) > 10 * width]
    comp = str.maketrans("ACGT", "TGCA")

    def background():
        while True:
            c = chroms[rng.integers(len(chroms))]
            a = int(rng.integers(0, len(gdict[c]) - width))
            sq = gdict[c][a:a + width].upper()
            if "N" not in sq:
                return list(sq)

    def site(k):
        sq = "".join("ACGT"[rng.choice(4, p=row / row.sum())] for row in probs[k])
        return sq.translate(comp)[::-1] if rng.random() < 0.5 else sq

    def plant(sq, k, centre):
        w = len(probs[k])
        a = width // 2 - w // 2 + int(rng.integers(-10, 11)) if centre else int(rng.integers(0, width - w))
        sq[a:a + w] = site(k)
    for name, (first, rate1, second) in {"T": (X, 0.5, Z), "B": (Y, 0.5, Z), "N": (None, 0, None)}.items():
        with open(d / f"{name}.fa", "w") as f:
            for i in range(n):
                sq = background()
                if first is not None and rng.random() < rate1:
                    plant(sq, first, centre=(name == "T"))
                if second is not None and rng.random() < 0.4:
                    plant(sq, second, centre=False)
                f.write(f">{name.lower()}{i}\n{''.join(sq)}\n")
    m = {"X": str(names[X]), "Y": str(names[Y]), "Z": str(names[Z]), "width": width, "n": n}
    meta.write_text(json.dumps(m))
    return {k: d / f"{k}.fa" for k in "TBN"}, m


def motif_new_jobs(rec, lib, mats, od, THR, PVAL, path):
    """The motif jobs added 2026-10-23, on planted-motif sets: known-motif enrichment vs background (C3), de novo
    discovery (C5), comparison to a motif database (C6) and differential enrichment between two sets (C8)."""
    names = [str(x) for x in lib.names]
    probs = [0.25 * np.power(2.0, lo) / (0.25 * np.power(2.0, lo)).sum(axis=1, keepdims=True) for lo in mats]
    gdict = gb.read_fasta(str(DATA / "genome.fa"))
    sets, pl = motif_sets(gdict, names, probs)
    T, B, N = (str(sets[k]) for k in "TBN")
    meme = str(od / "motifs.meme")                    # written by part_motifs (the same 100 motifs)
    width, n = pl["width"], pl["n"]
    top = 5

    def rank_score(ranked, want, avoid=()):
        """Share of the checks passed: each motif in ``want`` within the top 5, each in ``avoid`` outside it."""
        head = [str(x) for x in ranked[:top]]
        checks = [w in head for w in want] + [a not in head for a in avoid]
        where = {w: (ranked.index(w) + 1 if w in ranked else None) for w in list(want) + list(avoid)}
        return float(np.mean(checks)), "ranks: " + ", ".join(f"{k} #{v}" for k, v in where.items()) + f" (top {top} wanted for {', '.join(want)}" + (
            f"; not wanted for {', '.join(avoid)})" if avoid else ")")
    jobs = {"known_enrichment": (T, N, [pl["X"], pl["Z"]], []), "differential_enrichment": (T, B, [pl["X"]], [pl["Z"]])}
    # HOMER thresholds: each motif's log-odds at p = 1e-4 (MOODS' exact threshold), in HOMER's natural-log units
    import MOODS.tools
    thr = [MOODS.tools.threshold_from_p([list(col) for col in np.asarray(lo).T], [0.25] * 4, PVAL) for lo in mats]
    hm = str(od / "motifs_p.homer")
    with open(hm, "w") as f:
        for nm, pr, t in zip(names, probs, thr):
            f.write(f">{''.join('ACGT'[j] for j in pr.argmax(axis=1))}\t{nm}\t{t * np.log(2):.6f}\n" + "\n".join("\t".join(f"{x:.6f}" for x in row) for row in pr) + "\n")
    h2, sea, ame = HOMER_BIN / "homer2", tool("sea"), tool("ame")
    from genomeblocks.motifs import scan_motifs_matrix, compare_motifs, compare_motifs_to_ref
    for task, (fg, bg, want, avoid) in jobs.items():
        if h2.exists():
            out = {}

            def _hk(fg=fg, bg=bg):
                out["h"] = cli([str(h2), "known", "-i", fg, "-b", bg, "-m", hm, "-p", "1"])

            def _hk_score():
                rows = [l.split("\t") for l in out["h"].splitlines() if l and not l.startswith("Motif Name")]
                ranked = [r[0].split("\t")[0] for r in rows]
                ranked = [x.split("/")[0] if x not in names else x for x in ranked]
                return rank_score(ranked, want, avoid)
            rec_cli(rec, "motifs", task, "HOMER homer2 known (CLI)", n, _hk, _hk_score)
        if sea:
            out2 = {}

            def _sea(fg=fg, bg=bg):
                out2["s"] = cli([sea, "--text", "--p", fg, "--n", bg, "--m", meme])

            def _sea_score():
                d = pd.read_csv(io.StringIO(out2["s"]), sep="\t", comment="#")
                return rank_score(d.sort_values("LOG_PVALUE")["ID"].astype(str).tolist(), want, avoid)
            rec_cli(rec, "motifs", task, "MEME SEA (CLI)", n, _sea, _sea_score)
        if ame:
            out3 = {}

            def _ame(fg=fg, bg=bg):
                out3["a"] = cli([ame, "--text", "--control", bg, fg, meme])

            def _ame_score():
                d = pd.read_csv(io.StringIO(out3["a"]), sep="\t", comment="#")
                return rank_score(d.sort_values("rank")["motif_ID"].astype(str).tolist(), want, avoid)
            rec_cli(rec, "motifs", task, "MEME AME (CLI)", n, _ame, _ame_score)
        # genomeblocks: scan both sets (p = 1e-4, both strands), then its rank test (enrichment over a reference set,
        # or the two-sample differential test)
        Lf, Lb = (Loci.from_frame(pd.DataFrame({"chrom": [f"{k}{i}" for i in range(n)], "start": 0, "end": width})) for k in (Path(fg).stem.lower(), Path(bg).stem.lower()))

        def _gb(fg=fg, bg=bg, Lf=Lf, Lb=Lb, task=task):
            ma = scan_motifs_matrix(Lf, fg, path, format="jaspar16", r=width // 2, pvalue=PVAL, norm=False, both_strands=True, workers=1, verbose=False)
            mb = scan_motifs_matrix(Lb, bg, path, format="jaspar16", r=width // 2, pvalue=PVAL, norm=False, both_strands=True, workers=1, verbose=False)
            f = compare_motifs_to_ref if task == "known_enrichment" else compare_motifs
            return f(ma, mb, alternative="greater", verbose=False)

        def _gb_score(r):
            col = "Factor" if "Factor" in r.columns else r.columns[0]
            pcol = "p" if "p" in r.columns else next(c for c in r.columns if c.startswith("p_") and not c.startswith("p_adj"))   # p_<group>
            return rank_score(r.sort_values(pcol)[col].astype(str).tolist(), want, avoid)
        rec_py(rec, "motifs", task, "genomeblocks scan_motifs_matrix + " + ("compare_motifs_to_ref" if task == "known_enrichment" else "compare_motifs"),
               n, _gb, _gb_score)
    # C5 de novo: STREME, MEME, HOMER on T against N; score = how well the best discovered motif matches X and Z
    planted = [probs[names.index(pl["X"])], probs[names.index(pl["Z"])]]

    def denovo_score(found):
        if not found:
            return 0.0, "no motif found"
        best = [max(pwm_similarity(pw, f) for f in found.values()) for pw in planted]
        return float(np.mean(best)), f"best match to X {best[0]:.3f}, to Z {best[1]:.3f} ({len(found)} motifs)"
    streme, mm = tool("streme"), tool("meme")
    dn = od / "denovo"
    dn.mkdir(exist_ok=True)
    if streme:
        rec_cli(rec, "motifs", "de_novo", "MEME STREME (CLI)", n,
                lambda: cli([streme, "--p", T, "--n", N, "--dna", "--nmotifs", "5", "--oc", str(dn / "streme"), "--verbosity", "1"]),
                lambda: denovo_score(read_meme_probs(dn / "streme" / "streme.txt")), 1)
    if mm:
        mt = dn / "meme.txt"
        rec_cli(rec, "motifs", "de_novo", "MEME meme -mod zoops (CLI)", n,
                lambda: mt.write_text(cli([mm, T, "-dna", "-text", "-nmotifs", "3", "-revcomp", "-mod", "zoops", "-minw", "6", "-maxw", "20",
                                           "-nostatus", "-searchsize", "100000"])),
                lambda: denovo_score(read_meme_probs(mt)), 1)
    if h2.exists():
        rec_cli(rec, "motifs", "de_novo", "HOMER homer2 denovo (CLI)", n,
                lambda: cli([str(h2), "denovo", "-i", T, "-b", N, "-len", "10", "-S", "5", "-p", "1", "-o", str(dn / "homer.motifs")]),
                lambda: denovo_score(read_homer_motifs(dn / "homer.motifs")), 1)
    # C6 comparison to a database: 20 JASPAR motifs perturbed (Dirichlet noise, a column trimmed at each end) against
    # all 1,019; score = share of queries whose source motif is among the top 3 hits. genomeblocks only compares
    # motifs within one library (pwm_distance_matrix), so it is not timed
    full = gb.load_motifs(str(DATA / "jaspar.txt"), format="jaspar16")
    fnames = [str(x) for x in full.names]
    fm = []
    for i in range(len(full)):
        lo = np.asarray(full.logodds(i), float)
        fm.append(lo if lo.shape[1] == 4 else lo.T)
    dbm = str(od / "jaspar_all.meme")
    write_meme_probs(fnames, fm, dbm)
    rng = np.random.default_rng(3)
    qi = rng.choice([i for i in range(len(fm)) if len(fm[i]) >= 10], 20, replace=False)
    qn, qm = [], []
    for i in qi:
        pr = 0.25 * np.power(2.0, fm[i])
        pr = pr / pr.sum(axis=1, keepdims=True)
        noisy = np.array([rng.dirichlet(row * 30 + 0.5) for row in pr])[1:-1]
        qn.append(f"q_{fnames[i]}")
        qm.append(np.log2(noisy / 0.25))
    qf = str(od / "queries.meme")
    write_meme_probs(qn, qm, qf)
    src = {q: q[2:] for q in qn}

    def topk_score(hits):
        """hits: {query: [targets best first]}."""
        ok = [src[q] in hits.get(q, [])[:3] for q in qn]
        return float(np.mean(ok)), f"{sum(ok)}/{len(qn)} queries have their source motif in the top 3"
    tt = tool("tomtom")
    if tt:
        out4 = {}

        def _tt_score():
            d = pd.read_csv(io.StringIO(out4["t"]), sep="\t", comment="#")
            return topk_score({q: g.sort_values("p-value")["Target_ID"].astype(str).tolist() for q, g in d.groupby("Query_ID")})
        rec_cli(rec, "motifs", "db_comparison", "MEME Tomtom (CLI)", len(qn),
                lambda: out4.__setitem__("t", cli([tt, "-no-ssc", "-text", "-verbosity", "1", "-min-overlap", "5", "-dist", "pearson", "-thresh", "1", qf, dbm])),
                _tt_score)
    rows = worker("R", "motif_compare", od / "r_compare", qf, dbm)

    def _uc_score(r):
        d = read_tsv(r["out"], ["query", "target", "score"])
        return topk_score({q: g.sort_values("score", ascending=False)["target"].astype(str).tolist() for q, g in d.groupby("query")})
    record_worker(rec, "motifs", rows, n=len(qn), score=_uc_score)
    # C9 per-sample motif activity (chromVAR deviations): 2,000 peaks (set T, X planted; set B, Y planted) x 12
    # samples, group B twice as accessible at the T peaks, so X's activity rises and Y's falls in B. genomeblocks has
    # no per-sample activity score, so it is not timed
    md = WORK / "motif_sets"
    tb_fa, cnt_f, grp_f = md / "TB.fa", md / "chromvar_counts.tsv", md / "chromvar_groups.txt"
    if not cnt_f.exists():
        tb_fa.write_text(Path(T).read_text() + Path(B).read_text())
        rng = np.random.default_rng(17)
        n_s = 12
        grp = np.repeat(["A", "B"], n_s // 2)
        is_t = np.r_[np.ones(n, bool), np.zeros(n, bool)]
        mu = np.where(is_t[:, None] & (grp == "B")[None, :], 60.0, 30.0)
        cnt = rng.negative_binomial(10, 10 / (10 + mu))
        pd.DataFrame(cnt, index=[f"t{i}" for i in range(n)] + [f"b{i}" for i in range(n)], columns=[f"S{j + 1}" for j in range(n_s)]).to_csv(cnt_f, sep="\t")
        grp_f.write_text("\n".join(grp) + "\n")
    grp = np.array(grp_f.read_text().split())

    def activity_score(r):
        d = pd.read_csv(r["out"], sep="\t", header=None, index_col=0)
        d.index = [str(x).split()[0] for x in d.index]           # chromVAR names rows by the whole JASPAR header line
        z = d.to_numpy(float)
        a, b = z[:, grp == "A"], z[:, grp == "B"]
        t = (b.mean(1) - a.mean(1)) / np.sqrt(a.var(1, ddof=1) / a.shape[1] + b.var(1, ddof=1) / b.shape[1] + 1e-12)
        order = [str(x) for x in d.index[np.argsort(-np.nan_to_num(t))]]
        up, down = order[:3], order[::-1][:3]
        ok = [pl["X"] in up, pl["Y"] in down]
        return float(np.mean(ok)), (f"X ({pl['X']}) #{order.index(pl['X']) + 1} of {len(order)} by B - A t; Y ({pl['Y']}) #{order[::-1].index(pl['Y']) + 1} from the bottom "
                                    f"(top / bottom 3 wanted)")
    for kind in ("R", "py"):
        rows = worker(kind, "chromvar", od / f"{kind}_chromvar", str(tb_fa), str(cnt_f), path, str(grp_f))
        record_worker(rec, "motifs", rows, n=2 * n, score=activity_score)


# ── cooccupancy ──────────────────────────────────────────────────────────────

def part_cooccupancy(rec):
    print("\n== cooccupancy ==")
    tracks = sorted(glob.glob(str(DATA / "atlas" / "track_*.bed")))
    query = str(DATA / "atlas_query.bed")
    Q = Loci.make(query)
    T = len(tracks)
    names = [Path(p).stem for p in tracks]
    od = WORK / "cooccupancy"
    od.mkdir(parents=True, exist_ok=True)
    # the exact interval-level Fisher test: the reference ranking
    from scipy.stats import fisher_exact
    from genomeblocks.atlas import _fisher_vec
    loci_tracks = [Loci.make(p) for p in tracks]
    a = np.array([int(Q.overlap_any(t).sum()) for t in loci_tracks], np.int64)
    n_t = np.array([len(t) for t in loci_tracks], np.int64)
    from common import read_chromsizes
    G = sum(read_chromsizes().values())
    # genome-wide expectation at the query's footprint: hypergeometric on 1 kb bins of the genome
    nbins = G // 1000
    q_bins = len(Atlas.make([query], chromsizes=CS, bin_size=1000, workers=1, verbose=False).track_n_bins) and int(Atlas.make([query], chromsizes=CS, bin_size=1000, workers=1, verbose=False).track_n_bins[0])
    ref_p = np.array([fisher_exact([[a[j], len(Q) - a[j]], [max(0, int(n_t[j] - a[j])), max(0, nbins - len(Q) - int(n_t[j] - a[j]))]], alternative="greater")[1] for j in range(T)])
    ref_rank = -np.log10(np.maximum(ref_p, 1e-300))
    # designed truth: the fraction of each track drawn from the query's universe
    from make_data import make_peaks
    universe = Loci.from_frame(make_peaks(300_000, read_chromsizes(), np.random.default_rng(7)))
    truth = np.array([t.overlap_any(universe).mean() for t in loci_tracks])
    rec.add(part="cooccupancy", task="enrichment_search", engine="exact interval Fisher (reference)", n=T, seconds=None, runs=None, accuracy=1.0,
            detail=f"Spearman with the designed truth {spearman_score(truth, ref_rank)[0]:.4f}")
    res = {}
    t = timeit(lambda: res.__setitem__("A", Atlas.make(tracks, chromsizes=CS, bin_size=1000, workers=4, verbose=False)), repeat=1)
    A = res["A"]
    rec.add(part="cooccupancy", task="index_collection", engine="genomeblocks Atlas.make (1 kb bins, 4 workers)", n=T, seconds=t["median"], runs=t["runs"], accuracy=None)
    t = timeit(lambda: res.__setitem__("s", A.search(Q)), repeat=REPEAT)
    s = res["s"].set_index("name").loc[names]
    acc, d = spearman_score(ref_rank, s["giggle_score"].to_numpy())
    rec.add(part="cooccupancy", task="enrichment_search", engine="genomeblocks Atlas.search", n=T, seconds=t["median"], runs=t["runs"], accuracy=acc,
            detail=d + f" vs exact Fisher; vs designed truth {spearman_score(truth, s['giggle_score'].to_numpy())[0]:.4f}")

    def score(r):
        if r["op"] == "enrichment_search":
            d = read_tsv(r["out"], ["name", "support", "pValueLog", "oddsRatio"]).set_index("name")
            d = d.reindex(names)
            acc, det = spearman_score(ref_rank, d["pValueLog"].fillna(0).to_numpy())
            return acc, det + f" vs exact Fisher; vs designed truth {spearman_score(truth, d['pValueLog'].fillna(0).to_numpy())[0]:.4f}"
        if r["op"] == "permutation_null":
            d = read_tsv(r["out"], ["name", "observed", "z", "p"]).set_index("name")
            k = len(d)
            idx = [names.index(x) for x in d.index]
            acc, det = spearman_score(ref_rank[idx], d["z"].to_numpy())
            return acc, det + f" over {k} tracks"
        return None, ""
    rows = worker("R", "cooccupancy", od / "r", query, str(DATA / "atlas"), CS, 10, 100, env={"GB_REPEAT_SLOW": "1"})
    record_worker(rec, "cooccupancy", rows, n=T, score=score)
    # GIGGLE (C, on-disk index)
    if GIGGLE and Path(GIGGLE).exists():
        gd = od / "giggle"
        gd.mkdir(exist_ok=True)
        env = {**os.environ, "LD_LIBRARY_PATH": (GIGGLE_LD or "") + ":" + os.environ.get("LD_LIBRARY_PATH", "")}
        bgz = tool("bgzip")
        for p in tracks + [query]:
            dst = gd / (Path(p).name + ".gz")
            if not dst.exists():
                srt = subprocess.run(["sort", "-k1,1", "-k2,2n", p], capture_output=True, text=True, check=True).stdout
                with open(dst, "wb") as fh:
                    subprocess.run([bgz, "-c"], input=srt.encode(), stdout=fh, check=True)
        idx = od / "giggle_index"                 # GIGGLE wants the index beside the input directory, not inside it
        def _index():
            shutil.rmtree(idx, ignore_errors=True)
            subprocess.run([GIGGLE, "index", "-i", str(gd / "track_*.bed.gz"), "-o", str(idx), "-s", "-f"], capture_output=True, check=True, env=env)
        try:
            ti, runs = time_cli(_index, 1)
            rec.add(part="cooccupancy", task="index_collection", engine="GIGGLE index (CLI)", n=T, seconds=ti, runs=runs, accuracy=None)
            out = {}
            def _search():
                out["g"] = subprocess.run([GIGGLE, "search", "-i", str(idx), "-q", str(gd / (Path(query).name + ".gz")), "-s"],
                                          capture_output=True, text=True, check=True, env=env).stdout
            ts, runs = time_cli(_search, REPEAT)
            sc = {}
            for line in out["g"].splitlines():            # "#file ... combo_score" header, then one row per file
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\t\n").split("\t")
                sc[Path(p[0].strip()).name.replace(".bed.gz", "")] = float(p[7])
            v = np.array([sc.get(nm, 0.0) for nm in names])
            acc, d = spearman_score(ref_rank, v)
            rec.add(part="cooccupancy", task="enrichment_search", engine="GIGGLE search -s (CLI)", n=T, seconds=ts, runs=runs, accuracy=acc,
                    detail=d + f" vs exact Fisher; vs designed truth {spearman_score(truth, v)[0]:.4f}")
        except Exception as e:
            rec.add(part="cooccupancy", task="enrichment_search", engine="GIGGLE search -s (CLI)", n=T, seconds=None, runs=None, accuracy=None, error=str(e)[:300])
    # bedtools fisher, one call per track
    bt = tool("bedtools")
    srt_cs = str(od / "chrom.sizes.sorted")
    with open(srt_cs, "w") as fh:
        subprocess.run(["sort", "-k1,1", CS], stdout=fh, check=True)
    qs = str(od / "query.sorted.bed")
    with open(qs, "w") as fh:
        subprocess.run(["sort", "-k1,1", "-k2,2n", query], stdout=fh, check=True)
    def _fisher():
        pv = []
        for p in tracks:
            o = cli([bt, "fisher", "-a", qs, "-b", p, "-g", srt_cs])
            vals = [l for l in o.splitlines() if l and not l.startswith("#")]
            pv.append(float(vals[-1].split("\t")[1]))      # right-tail p
        return np.array(pv)
    try:
        t = timeit(lambda: res.__setitem__("f", _fisher()), repeat=1, warmup=0)
        v = -np.log10(np.maximum(res["f"], 1e-300))
        acc, d = spearman_score(ref_rank, v)
        rec.add(part="cooccupancy", task="enrichment_search", engine="bedtools fisher ×500 (CLI)", n=T, seconds=t["median"], runs=t["runs"], accuracy=acc,
                detail=d + f" vs exact Fisher; vs designed truth {spearman_score(truth, v)[0]:.4f}")
    except Exception as e:
        rec.add(part="cooccupancy", task="enrichment_search", engine="bedtools fisher ×500 (CLI)", n=T, seconds=None, runs=None, accuracy=None, error=str(e)[:300])

    cooccupancy_universe_jobs(rec, A, Q, tracks, names, universe, od)
    cooccupancy_pairwise_jobs(rec, tracks[:20], od)


def bp_merged(df):
    """Merged (chrom, start, end) arrays of an interval table, as one sorted frame."""
    d = df.sort_values(["chrom", "start"])
    out = []
    for c, g in d.groupby("chrom", sort=False):
        s, e = g["start"].to_numpy(np.int64), np.maximum.accumulate(g["end"].to_numpy(np.int64))
        first = np.flatnonzero(np.r_[True, s[1:] > e[:-1]])
        last = np.r_[first[1:], len(s)] - 1
        out.append(pd.DataFrame({"chrom": c, "start": s[first], "end": e[last]}))
    return pd.concat(out, ignore_index=True)


def universe_reference(U, Q, tracks):
    """C5's input and reference: the query as universe regions, and -log10 P of the exact Fisher test of every
    track's universe-region counts (query against the rest of the universe)."""
    from scipy.stats import fisher_exact
    Qu = U.take(np.flatnonzero(U.overlap_any(Q)))            # the query as universe regions (all engines get this)
    nU, nQ = len(U), len(Qu)
    ref = []
    for p in tracks:
        t = Loci.make(p)
        a = int(Qu.overlap_any(t).sum())
        c = int(U.overlap_any(t).sum()) - a
        ref.append(fisher_exact([[a, nQ - a], [c, nU - nQ - c]], alternative="greater")[1])
    return Qu, -np.log10(np.maximum(np.array(ref), 1e-300))


def gb_universe_bootstrap(rec, part, A, Qu, U, names, ref_rank, nperm=100):
    """genomeblocks' C5 row: Atlas.bootstrap with the universe as its pool, ranked by z against the reference."""
    res = {}
    try:
        t = timeit(lambda: res.__setitem__("b", A.bootstrap(Qu, n=nperm, pool=U, seed=1, verbose=False)), repeat=1)
        b = res["b"].set_index("name").reindex(names)
        acc, d = spearman_score(ref_rank, b["z"].fillna(0).to_numpy())
        rec.add(part=part, task="universe_enrichment", engine=f"genomeblocks Atlas.bootstrap (pool = universe, {nperm} draws)", n=len(names),
                seconds=t["median"], runs=t["runs"], accuracy=acc, detail=d + " vs the exact universe Fisher test")
    except Exception as e:
        rec.add(part=part, task="universe_enrichment", engine="genomeblocks Atlas.bootstrap (pool = universe)", n=len(names), seconds=None,
                runs=None, accuracy=None, error=f"{type(e).__name__}: {e}"[:300])


def part_cooccupancy_bootstrap(rec):
    """C5 for genomeblocks alone (2026-11-12, after Atlas.bootstrap counts regions instead of distinct bins): the
    inputs and reference of cooccupancy_universe_jobs, rebuilt from the seeded data, and Atlas.bootstrap's row. A
    part of its own so it runs without LOLA / regioneR's R environment and keeps their rows."""
    print("\n== cooccupancy_bootstrap ==")
    from common import read_chromsizes
    from make_data import make_peaks
    tracks = sorted(glob.glob(str(DATA / "atlas" / "track_*.bed")))
    names = [Path(p).stem for p in tracks]
    Q = Loci.make(str(DATA / "atlas_query.bed"))
    U = Loci.from_frame(make_peaks(300_000, read_chromsizes(), np.random.default_rng(7)))
    A = Atlas.make(tracks, chromsizes=CS, bin_size=1000, workers=4, verbose=False)
    Qu, ref_rank = universe_reference(U, Q, tracks)
    print(f"  {len(Qu):,} of {len(U):,} universe regions in the query; {len(tracks)} tracks")
    rec.add(part="cooccupancy_bootstrap", task="universe_enrichment", engine="exact Fisher in universe regions (reference)", n=len(tracks),
            seconds=None, runs=None, accuracy=1.0, detail=f"{len(Qu):,} of {len(U):,} universe regions in the query")
    gb_universe_bootstrap(rec, "cooccupancy_bootstrap", A, Qu, U, names, ref_rank)


def cooccupancy_universe_jobs(rec, A, Q, tracks, names, universe, od, n_regioner=10, nperm=100):
    """C5 (added 2026-10-23): enrichment against a custom universe. The query and every track are counted in
    universe regions (300k simulated CREs); the reference is the exact Fisher test on those counts. LOLA takes the
    universe as userUniverse and regioneR resamples within it; genomeblocks' Atlas.bootstrap (pool = the universe)
    runs in part_cooccupancy_bootstrap."""
    U = universe
    Qu, ref_rank = universe_reference(U, Q, tracks)
    ub, qb = write_bed(U, od / "universe.bed"), write_bed(Qu, od / "query_in_universe.bed")
    nU, nQ = len(U), len(Qu)
    rec.add(part="cooccupancy", task="universe_enrichment", engine="exact Fisher in universe regions (reference)", n=len(tracks), seconds=None, runs=None,
            accuracy=1.0, detail=f"{nQ:,} of {nU:,} universe regions in the query")

    def score(r):
        if r["op"] != "universe_enrichment":
            return None, ""
        d = read_tsv(r["out"], ["name", "score"]).set_index("name")
        idx = [names.index(x) for x in d.index if x in names]
        acc, det = spearman_score(ref_rank[idx], d["score"].reindex([names[i] for i in idx]).fillna(0).to_numpy())
        return acc, det + f" over {len(idx)} tracks vs the exact universe Fisher test"
    rows = worker("R", "universe", od / "r_universe", qb, ub, str(DATA / "atlas"), n_regioner, nperm, env={"GB_REPEAT_SLOW": "1"})
    record_worker(rec, "cooccupancy", rows, n=len(tracks), score=score)


def cooccupancy_pairwise_jobs(rec, tracks, od):
    """C7 (added 2026-10-23): the pairwise co-occurrence matrix of 20 tracks as base-pair Jaccard (intersection /
    union of the merged sets), against an exact numpy computation. genomeblocks has no pairwise matrix."""
    k = len(tracks)
    frames = [pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"]) for p in tracks]
    merged = [bp_merged(f) for f in frames]
    bp = np.array([int((m["end"] - m["start"]).sum()) for m in merged])
    ref = np.eye(k)
    for i in range(k):
        for j in range(i + 1, k):
            union = int((lambda m: (m["end"] - m["start"]).sum())(bp_merged(pd.concat([merged[i], merged[j]]))))
            ref[i, j] = ref[j, i] = (bp[i] + bp[j] - union) / union
    iu = np.triu_indices(k, 1)
    names = [Path(p).stem for p in tracks]
    pd_ = od / "pairwise"
    pd_.mkdir(exist_ok=True)
    srt = []
    for p in tracks:                                          # bedtools jaccard wants sorted input (made outside the clock)
        q = pd_ / Path(p).name
        with open(q, "w") as fh:
            subprocess.run(["sort", "-k1,1", "-k2,2n", p], stdout=fh, check=True)
        srt.append(str(q))

    def mscore(m):
        return corr_score(ref[iu], np.asarray(m, float)[iu], tol=1e-6)
    bt = tool("bedtools")

    def _bt():
        m = np.eye(k)
        for i in range(k):
            for j in range(i + 1, k):
                o = cli([bt, "jaccard", "-a", srt[i], "-b", srt[j]]).splitlines()
                m[i, j] = m[j, i] = float(o[1].split("\t")[2])
        return m
    rec_py(rec, "cooccupancy", "pairwise_cooccurrence", f"bedtools jaccard x {len(iu[0])} pairs (CLI)", k, _bt, mscore, rep=1)
    iv = clitool("intervene")
    if iv:
        ivd = pd_ / "intervene"

        def _iv():
            shutil.rmtree(ivd, ignore_errors=True)
            env = {**os.environ, "PATH": f"{Path(iv).parent}:{os.environ['PATH']}"}
            pr = subprocess.run([iv, "pairwise", "-i", *srt, "--names", ",".join(names), "--compute", "jaccard", "--htype", "color", "-o", str(ivd)],
                                capture_output=True, text=True, env=env)
            if pr.returncode != 0:
                raise RuntimeError("intervene: " + pr.stderr.strip().splitlines()[-1][:250])
            f = next(iter(glob.glob(str(ivd / "**" / "*jaccard*matrix*.txt"), recursive=True)))
            return pd.read_csv(f, sep="\t", index_col=0).reindex(index=names, columns=names).to_numpy(float)
        rec_py(rec, "cooccupancy", "pairwise_cooccurrence", "intervene pairwise --compute jaccard (CLI)", k, _iv, mscore, rep=1)
    rows = worker("R", "pairwise", pd_ / "r", *srt)
    record_worker(rec, "cooccupancy", rows, n=k, score=lambda r: mscore(pd.read_csv(r["out"], sep="\t", header=None).to_numpy(float)))


# ── graph ────────────────────────────────────────────────────────────────────

def part_graph(rec):
    print("\n== graph ==")
    bedpe = str(DATA / "loops.bedpe")
    cre_bed = str(DATA / "peaks_A_100000.bed")          # the loops are anchored on these: the graph's vertices
    reg_bed = str(DATA / "peaks_B_100000.bed")          # partially co-located: a non-trivial overlap test
    cre, reg = Loci.make(cre_bed), Loci.make(reg_bed)
    od = WORK / "graph"
    od.mkdir(parents=True, exist_ok=True)
    res = {}
    t = timeit(lambda: res.__setitem__("P", Pairs.make(bedpe)), repeat=REPEAT)
    P = res["P"]
    rec.add(part="graph", task="read_bedpe", engine="genomeblocks Pairs.make", n=len(P), seconds=t["median"], runs=t["runs"], accuracy=None)
    bt = tool("bedtools")
    out = {}
    ref = {}
    for typ in ("either", "both"):
        tb, runs = time_cli(lambda typ=typ: out.__setitem__(typ, cli([bt, "pairtobed", "-a", bedpe, "-b", reg_bed, "-type", typ])), REPEAT)
        ids = {line.split("\t")[6] for line in out[typ].splitlines() if line}
        ref[typ] = {int(x[1:]) for x in ids}
        rec.add(part="graph", task=f"loop_overlap_{typ}", engine="bedtools pairtobed (CLI)", n=len(P), seconds=tb, runs=runs, accuracy=1.0, detail="reference", n_out=len(ids))
        t = timeit(lambda typ=typ: res.__setitem__(typ, P.overlapping(reg, both=(typ == "both"))), repeat=REPEAT, setup=lambda: (P.a._dirty(), P.b._dirty(), reg._dirty()))
        got = {int(x[1:]) for x in res[typ]["name"]} if "name" in res[typ].cols else set()
        acc, d = set_score(ref[typ], got)
        rec.add(part="graph", task=f"loop_overlap_{typ}", engine="genomeblocks Pairs.overlapping", n=len(P), seconds=t["median"], runs=t["runs"], accuracy=acc, detail=d, n_out=len(got))
    # the contact graph and its algorithms
    A = quiet(lambda: Architecture.make(cre, bedpe, r=2500, verbose=False))()
    edges = str(od / "edges.tsv")
    pd.DataFrame({"src": A.src, "tgt": A.tgt, "w": np.ones(A.n_links)}).to_csv(edges, sep="\t", header=False, index=False)
    A.ep["one"] = np.ones(A.n_links)
    t = timeit(lambda: res.__setitem__("c", A.components()), repeat=REPEAT)
    comp = np.asarray(res["c"])

    def comp_sizes(member):
        """Multiset of component sizes, singletons left out (isolated vertices are -1 in some engines)."""
        m = np.asarray(member)
        s = np.bincount(m[m >= 0] - m[m >= 0].min()) if (m >= 0).any() else np.array([], int)
        return np.sort(s[s >= 2])
    ref_sizes = comp_sizes(comp)
    rec.add(part="graph", task="components", engine=f"genomeblocks ({gb.backends.resolve('graph', None)})", n=A.n_links, seconds=t["median"], runs=t["runs"], accuracy=1.0, detail=f"reference: {len(ref_sizes):,} components")
    t = timeit(lambda: res.__setitem__("pr", A.pagerank("one")), repeat=REPEAT)
    ref_pr = np.asarray(res["pr"], float)
    rec.add(part="graph", task="pagerank", engine=f"genomeblocks ({gb.backends.resolve('graph', None)})", n=A.n_links, seconds=t["median"], runs=t["runs"], accuracy=1.0, detail="reference")

    def score(r):
        if r["op"].startswith("loop_overlap"):
            d = read_tsv(r["out"], ["idx"])
            return set_score(ref[r["op"].split("_")[-1]], set(d["idx"].astype(int)))
        if r["op"] == "components":
            d = read_tsv(r["out"], ["member"])
            sizes = comp_sizes(d["member"].to_numpy(int))
            same = len(sizes) == len(ref_sizes) and np.array_equal(sizes, ref_sizes)
            return (1.0 if same else 0.0), f"{len(sizes):,} components of 2+ vertices{'' if same else ' (differs)'}"
        if r["op"] == "pagerank":
            d = read_tsv(r["out"], ["pr"])
            return corr_score(ref_pr, d["pr"].to_numpy(), tol=1e-6)
        return None, ""
    rows = worker("R", "graph", od / "r", bedpe, reg_bed, edges, len(cre))
    record_worker(rec, "graph", rows, n=len(P), score=score)


# ── expression ───────────────────────────────────────────────────────────────

def part_expression(rec):
    print("\n== expression ==")
    size = os.environ.get("GB_BENCH_EXPR_SIZE", "60000_24")
    qdir = DATA / "expression" / f"salmon_{size}"
    gtf = str(DATA / "expression" / f"genes_{size.split('_')[0]}.gtf")
    G = Genes.make(gtf)
    files = sorted(glob.glob(str(qdir / "*" / "quant.sf")))
    od = WORK / "expression"
    od.mkdir(parents=True, exist_ok=True)
    res = {}
    t = timeit(lambda: res.__setitem__("E", quiet(lambda: Expression.make(files, genes=G, level="transcript", verbose=False).to_genes(verbose=False))()), repeat=REPEAT)
    E = res["E"]
    rec.add(part="expression", task="tx_to_gene", engine="genomeblocks Expression.make + to_genes (length_scaled_tpm)", n=E.shape[0], n_samples=E.n_samples,
            seconds=t["median"], runs=t["runs"], accuracy=1.0, detail="reference (equals R tximport to 1e-9, see the Expression bench)")
    counts = pd.DataFrame(E.counts, index=G.genes["gene_id"] if "gene_id" in G.genes.cols else G.genes.names, columns=E.sample_names)
    counts = counts.fillna(0.0)
    cf = str(od / "gene_counts.tsv")
    counts.to_csv(cf, sep="\t")
    t2g = str(od / "tx2gene.tsv")
    tx = G.transcripts
    gid = G.genes["gene_id"] if "gene_id" in G.genes.cols else G.genes.names
    pd.DataFrame({"tx": tx["transcript_id"] if "transcript_id" in tx.cols else tx.names, "gene": np.asarray(gid)[tx["gene"]]}).to_csv(t2g, sep="\t", header=False, index=False)
    t = timeit(lambda: res.__setitem__("tmm", E.tmm()), repeat=REPEAT)
    ref_tmm = np.asarray(res["tmm"], float)
    rec.add(part="expression", task="tmm", engine="genomeblocks Expression.tmm", n=E.shape[0], n_samples=E.n_samples, seconds=t["median"], runs=t["runs"],
            accuracy=1.0, detail="reference (equals edgeR calcNormFactors to 4e-16, see the Expression bench)")
    t = timeit(lambda: res.__setitem__("sf", E.size_factors()), repeat=REPEAT)
    ref_sf = np.asarray(res["sf"], float)
    rec.add(part="expression", task="size_factors", engine="genomeblocks Expression.size_factors", n=E.shape[0], n_samples=E.n_samples, seconds=t["median"], runs=t["runs"],
            accuracy=1.0, detail="reference (equals DESeq2 estimateSizeFactors, see the Expression bench)")
    t = timeit(lambda: res.__setitem__("cpm", E.cpm()), repeat=REPEAT)
    ref_cpm = np.nan_to_num(np.asarray(res["cpm"], float))
    rec.add(part="expression", task="cpm", engine="genomeblocks Expression.cpm", n=E.shape[0], n_samples=E.n_samples, seconds=t["median"], runs=t["runs"], accuracy=1.0, detail="reference")
    gene_counts = np.nan_to_num(np.asarray(E.counts, float))
    gene_ids = np.asarray(gid).astype(str)

    def score(r):
        if r["op"] in ("tmm", "size_factors"):
            d = read_tsv(r["out"], ["sample", "factor"]).set_index("sample").reindex(E.sample_names)
            ref = ref_tmm if r["op"] == "tmm" else ref_sf
            return corr_score(ref, d["factor"].to_numpy(), tol=1e-9)
        if r["op"] == "cpm":
            d = pd.read_csv(r["out"], sep="\t", header=None, index_col=0).reindex(counts.index)
            return corr_score(ref_cpm, np.nan_to_num(d.to_numpy(float)), tol=1e-9)
        if r["op"] == "tx_to_gene":
            d = pd.read_csv(r["out"], sep="\t", header=None, index_col=0)
            d.columns = E.sample_names[:d.shape[1]]
            d = d.reindex(gene_ids)
            got = np.nan_to_num(d.to_numpy(float))
            return corr_score(gene_counts, got, tol=1e-6)
        return None, ""
    rows = worker("py", "expression", od / "py", cf, str(qdir), t2g)
    record_worker(rec, "expression", rows, n=E.shape[0], score=score, n_samples=E.n_samples)
    # the confirmed core jobs added 2026-10-23 (landscape/jobs.json): DE test (C7), dispersion estimation (C6),
    # LFC shrinkage (C8), variance-stabilising transform (C9), batch correction (C10), on simulated counts with
    # a known truth. genomeblocks has none of them (DE is planned for 2.1), so it is not timed here
    for ng, ns in DE_SIZES:
        f = de_dataset(ng, ns)
        dd = od / f"de_{ng}_{ns}"
        rows = worker("R", "expression_de", dd / "r", f["counts"], f["coldata"]) + worker("py", "expression_de", dd / "py", f["counts"], f["coldata"])
        record_worker(rec, "expression", rows, n=ng, score=de_scorer(f, rows), n_samples=ns)
    # C1 (added 2026-10-23): reads per gene from an RNA-seq BAM. Every simulated read lies in exons of exactly one
    # gene, so the simulated counts are the answer; genomeblocks reads count tables but does not count BAMs
    gtf_bam = str(DATA / "genes.gtf")
    for nr in BAM_READS:
        bam, truth_f = rna_bam(nr)
        tr = pd.read_csv(truth_f, sep="\t").set_index("gene")["count"]

        def gscore(df, tr=tr):
            """df: gene, count (any extra summary rows are ignored)."""
            df = df[~df.iloc[:, 0].astype(str).str.startswith("__")]
            got = pd.Series(df.iloc[:, 1].to_numpy(), index=df.iloc[:, 0].astype(str)).reindex(tr.index).fillna(-1).astype(int)
            same = int((got == tr).sum())
            return (1.0 if same == len(tr) else same / len(tr)), f"{same:,}/{len(tr):,} genes equal the simulated counts"
        rep = REPEAT if nr <= 1_000_000 else 1
        fc = clitool("featureCounts")
        if fc:
            fco = str(od / f"rna_featurecounts_{nr}.txt")
            rec_cli(rec, "expression", "read_counts_bam", "Subread featureCounts (CLI)", nr,
                    lambda: cli([fc, "-a", gtf_bam, "-t", "exon", "-g", "gene_id", "-T", "1", "-o", fco, str(bam)]),
                    lambda: gscore(pd.read_csv(fco, sep="\t", comment="#").iloc[:, [0, 6]]), rep)
        ht = clitool("htseq-count")
        if ht:
            out = {}
            rec_cli(rec, "expression", "read_counts_bam", "HTSeq htseq-count (CLI)", nr,
                    lambda: out.__setitem__("h", cli([ht, "-f", "bam", "-r", "pos", "-s", "no", "-t", "exon", "-i", "gene_id", str(bam), gtf_bam])),
                    lambda: gscore(pd.read_csv(io.StringIO(out["h"]), sep="\t", header=None)), 1)
        rows = worker("R", "count_bam", od / f"count_bam_{nr}", str(bam), gtf_bam)
        record_worker(rec, "expression", rows, n=nr, score=lambda r, g=gscore: g(pd.read_csv(r["out"], sep="\t", header=None)))


# ── simulated BAM files (signal and expression jobs added 2026-10-23) ───────
BAMS = Path(os.environ.get("GB_LANDSCAPE_BAMS", DATA / "landscape_bams"))
READ_LEN = 50
BAM_READS = [int(x) for x in os.environ.get("GB_BENCH_BAM_READS", "1000000,10000000").split(",")]


def chrom_sizes():
    s = pd.read_csv(CS, sep="\t", header=None, names=["chrom", "size"])
    return s["chrom"].tolist(), s["size"].to_numpy(np.int64)


def write_bam(path, chroms, sizes, ci, starts, rev):
    """Single-end READ_LEN-bp reads (chrom index, 0-based start, reverse?) as a sorted, indexed BAM. The SAM
    text is built with numpy and converted by samtools, which is far faster than pysam record by record."""
    o = np.lexsort((starts, ci))
    ci, starts, rev = ci[o], starts[o], rev[o]
    sam = Path(str(path) + ".sam")
    with open(sam, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.writelines(f"@SQ\tSN:{c}\tLN:{int(s)}\n" for c, s in zip(chroms, sizes))
        names = np.asarray(chroms, dtype=object)
        step = 2_000_000
        for a in range(0, len(starts), step):
            sl = slice(a, a + step)
            q = np.char.add("r", np.arange(a, a + len(starts[sl])).astype(str))
            flag = np.where(rev[sl], "16", "0")
            pos = (starts[sl] + 1).astype(str)
            cols = [q, flag, names[ci[sl]].astype(str), pos, np.full(len(q), f"60\t{READ_LEN}M\t*\t0\t0\t*\t*")]
            f.write("\n".join("\t".join(r) for r in zip(*cols)) + "\n")
    samtools = tool("samtools")
    subprocess.run([samtools, "view", "-b", "-o", str(path), str(sam)], check=True)
    subprocess.run([samtools, "index", str(path)], check=True)
    sam.unlink()


def chip_bams(n_reads, n_bams=4, seed=21):
    """n_bams ChIP-like BAMs of n_reads each: 30 % of the reads within +/- 200 bp of 20,000 peak summits
    (peaks_A_100000.bed), the rest uniform over the genome. Each BAM's read starts are kept beside it
    (.npz) so the truth (counts per region, per-base coverage) is computed exactly."""
    d = BAMS / f"chip_{n_reads}"
    paths = [d / f"chip_{j}.bam" for j in range(n_bams)]
    if all(p.exists() and Path(str(p) + ".npz").exists() for p in paths):
        return paths
    d.mkdir(parents=True, exist_ok=True)
    chroms, sizes = chrom_sizes()
    code = {c: i for i, c in enumerate(chroms)}
    pk = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None, usecols=[0, 1, 2])
    rng = np.random.default_rng(seed)
    pk = pk.iloc[rng.choice(len(pk), 20_000, replace=False)]
    pc, pmid = pk[0].map(code).to_numpy(), ((pk[1] + pk[2]) // 2).to_numpy(np.int64)
    w = sizes / sizes.sum()
    for j, p in enumerate(paths):
        rng = np.random.default_rng(seed + 1 + j)
        n_pk = int(n_reads * 0.3)
        k = rng.integers(0, len(pc), n_pk) if j < n_bams // 2 else rng.integers(0, len(pc) // 2, n_pk)   # the second half differ
        ci = np.concatenate([pc[k], rng.choice(len(chroms), n_reads - n_pk, p=w)])
        st = np.concatenate([pmid[k] + rng.integers(-200, 200, n_pk), (rng.random(n_reads - n_pk) * (sizes[ci[n_pk:]] - READ_LEN)).astype(np.int64)])
        st = np.clip(st, 0, sizes[ci] - READ_LEN)
        write_bam(p, chroms, sizes, ci, st, rng.random(n_reads) < 0.5)
        np.savez(str(p) + ".npz", ci=ci, st=st)
    return paths


def read_truth(npz):
    """{chrom index: sorted read starts} of a simulated BAM."""
    z = np.load(npz)
    o = np.lexsort((z["st"], z["ci"]))
    ci, st = z["ci"][o], z["st"][o]
    cut = np.r_[0, np.flatnonzero(np.diff(ci)) + 1, len(ci)]
    return {int(ci[a]): st[a:b] for a, b in zip(cut[:-1], cut[1:])}


def reads_over(starts, a, b):
    """Reads [s, s + READ_LEN) overlapping [a, b) by at least one base, per region (sorted starts)."""
    return np.searchsorted(starts, b, side="left") - np.searchsorted(starts, a - READ_LEN, side="right")


def rna_bam(n_reads, seed=31):
    """An RNA-seq BAM over data/genes.gtf: reads lie wholly inside exonic stretches that belong to one gene
    (merged exons that touch no other gene's exon), so every counting tool should return the true per-gene
    counts exactly. Returns (bam, truth.tsv)."""
    d = BAMS / f"rna_{n_reads}"
    bam, truth = d / "rna.bam", d / "truth.tsv"
    if bam.exists() and truth.exists():
        return bam, truth
    d.mkdir(parents=True, exist_ok=True)
    chroms, sizes = chrom_sizes()
    code = {c: i for i, c in enumerate(chroms)}
    g = pd.read_csv(DATA / "genes.gtf", sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 8])
    g = g[g[2] == "exon"]
    ex = pd.DataFrame({"c": g[0].map(code), "s": g[3] - 1, "e": g[4], "gene": g[8].str.extract(r'gene_id "([^"]+)"')[0]}).dropna()
    # merge each gene's exons, then drop stretches that touch another gene's exon
    ex = ex.sort_values(["gene", "c", "s"])
    rows = []
    for (gene, c), t in ex.groupby(["gene", "c"], sort=False):
        s, e = t["s"].to_numpy(), np.maximum.accumulate(t["e"].to_numpy())
        first = np.flatnonzero(np.r_[True, s[1:] > e[:-1]])
        last = np.r_[first[1:], len(s)] - 1
        rows.append(pd.DataFrame({"c": int(c), "s": s[first], "e": e[last], "gene": gene}))
    m = pd.concat(rows).sort_values(["c", "s"]).reset_index(drop=True)
    cm = m.groupby("c")["e"].cummax().shift(1)
    touch_prev = (m["c"] == m["c"].shift(1)) & (m["s"] < cm)
    touch_next = (m["c"] == m["c"].shift(-1)) & (m["e"] > m["s"].shift(-1))
    m = m[~(touch_prev | touch_next) & (m["e"] - m["s"] >= READ_LEN + 10)]
    rng = np.random.default_rng(seed)
    genes = m["gene"].unique()
    level = pd.Series(rng.lognormal(0, 1.5, len(genes)), index=genes)
    wt = (m["e"] - m["s"] - READ_LEN) * m["gene"].map(level).to_numpy()
    k = rng.choice(len(m), n_reads, p=(wt / wt.sum()).to_numpy())
    mc, ms, me = m["c"].to_numpy()[k], m["s"].to_numpy()[k], m["e"].to_numpy()[k]
    st = ms + (rng.random(n_reads) * (me - ms - READ_LEN)).astype(np.int64)
    write_bam(bam, chroms, sizes, mc.astype(np.int64), st.astype(np.int64), rng.random(n_reads) < 0.5)
    cnt = pd.Series(m["gene"].to_numpy()[k]).value_counts()
    allg = ex["gene"].unique()
    pd.DataFrame({"gene": allg, "count": cnt.reindex(allg).fillna(0).astype(int).to_numpy()}).to_csv(truth, sep="\t", index=False)
    return bam, truth


DE_SIZES = [tuple(int(v) for v in s.split("x")) for s in os.environ.get("GB_BENCH_DE_SIZES", "20000x12,60000x48").split(",")]


def de_dataset(n_genes, n_samples, seed=7):
    """Counts with a known truth for the DE family: two groups x two batches (balanced; n_samples is a multiple
    of 4), 10 % of genes differential (|log2 FC| ~ 1.5), a per-gene batch shift (sd 0.5 log2), and
    negative-binomial counts with a known dispersion per gene (trend 4 / mean + 0.05, log-normal scatter).
    Written once under WORK/expression_de; returns the three file paths."""
    d = WORK / "expression_de" / f"{n_genes}_{n_samples}"
    files = {k: d / f"{k}.tsv" for k in ("counts", "coldata", "truth")}
    if all(p.exists() for p in files.values()):
        return files
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    samples = [f"S{j + 1}" for j in range(n_samples)]
    group = np.repeat(["A", "B"], n_samples // 2)
    batch = np.tile(np.repeat(["1", "2"], n_samples // 4), 2)
    mu0 = rng.lognormal(np.log(200), 1.8, n_genes)
    de = rng.random(n_genes) < 0.10
    lfc = np.where(de, rng.choice([-1.0, 1.0], n_genes) * rng.normal(1.5, 0.4, n_genes), 0.0)
    blfc = rng.normal(0, 0.5, n_genes)
    disp = (4 / mu0 + 0.05) * rng.lognormal(0, 0.5, n_genes)
    lib = rng.lognormal(0, 0.2, n_samples)
    mu = 2 ** (np.log2(mu0)[:, None] + lfc[:, None] * (group == "B") + blfc[:, None] * (batch == "2")) * lib
    r = 1 / disp[:, None]
    counts = rng.negative_binomial(r, r / (r + mu))
    genes = [f"G{i + 1}" for i in range(n_genes)]
    pd.DataFrame(counts, index=genes, columns=samples).to_csv(files["counts"], sep="\t")
    pd.DataFrame({"group": group, "batch": batch}, index=samples).to_csv(files["coldata"], sep="\t")
    pd.DataFrame({"de": de.astype(int), "lfc": lfc, "batch_lfc": blfc, "disp": disp, "mu": mu0}, index=genes).to_csv(files["truth"], sep="\t")
    return files


def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1.0)
    return out


def log2cpm(m):
    m = np.asarray(m, float)
    return np.log2(m / m.sum(axis=0, keepdims=True) * 1e6 + 1)


def de_scorer(files, rows):
    """score(row) for the DE family against the simulated truth (vst against DESeq2's own vst)."""
    from scipy.stats import rankdata, spearmanr
    truth = pd.read_csv(files["truth"], sep="\t", index_col=0)
    cd = pd.read_csv(files["coldata"], sep="\t", index_col=0, dtype=str)
    raw = pd.read_csv(files["counts"], sep="\t", index_col=0).reindex(truth.index).to_numpy(float)
    g, b = (cd["group"] == "B").to_numpy(), (cd["batch"] == "2").to_numpy()
    expressed = (truth["mu"] > 20).to_numpy()
    vst_ref = next((r["out"] for r in rows if r["op"] == "vst" and r["engine"] == "DESeq2 vst" and not r.get("error")), None)

    def effects(lg):
        """(within-group batch difference, within-batch group difference), per gene, averaged over the strata."""
        bd = np.mean([lg[:, (g == s) & b].mean(1) - lg[:, (g == s) & ~b].mean(1) for s in (True, False)], axis=0)
        gd = np.mean([lg[:, g & (b == s)].mean(1) - lg[:, ~g & (b == s)].mean(1) for s in (True, False)], axis=0)
        return bd, gd
    bd0, gd0 = effects(log2cpm(raw))
    c0 = float(np.corrcoef(bd0[expressed], truth["batch_lfc"].to_numpy()[expressed])[0, 1])
    k0 = float(np.corrcoef(gd0[expressed], truth["lfc"].to_numpy()[expressed])[0, 1])     # the group signal before

    def table(path, names):
        return pd.read_csv(path, sep="\t", header=None, names=names, index_col=0).reindex(truth.index)

    def score(r):
        op, path = r["op"], r["out"]
        if op == "de_test":
            p = table(path, ["gene", "p", "lfc"])["p"].fillna(1.0).to_numpy()
            y = truth["de"].to_numpy() == 1
            rk = rankdata(-np.log10(np.clip(p, 1e-300, 1.0)))
            auc = float((rk[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum()))
            called = bh(p) < 0.05
            tp = int((called & y).sum())
            return auc, (f"AUROC {auc:.4f} against the simulated DE genes; at BH 5 %: {int(called.sum()):,} called, "
                         f"{tp / max(1, called.sum()):.1%} true, {tp / y.sum():.1%} of the true DE genes found")
        if op == "dispersion":
            d = table(path, ["gene", "disp"])["disp"].to_numpy(float)
            ok = np.isfinite(d) & (d > 0)
            rho = float(spearmanr(np.log(d[ok]), np.log(truth["disp"].to_numpy()[ok])).statistic)
            return max(0.0, rho), f"Spearman {rho:.4f} with the simulated dispersions ({ok.sum():,} genes)"
        if op == "lfc_shrinkage":
            d = table(path, ["gene", "lfc"])["lfc"].to_numpy(float)
            ok = np.isfinite(d)
            t = truth["lfc"].to_numpy()[ok]
            rr = float(np.corrcoef(d[ok], t)[0, 1])
            return max(0.0, rr), f"Pearson {rr:.4f} with the simulated log2 FC; RMSE {np.sqrt(np.mean((d[ok] - t) ** 2)):.3f}"
        if op == "vst":
            if vst_ref is None:
                return None, "no DESeq2 vst to compare with"
            if path == vst_ref:
                return 1.0, "reference"
            ref = pd.read_csv(vst_ref, sep="\t", header=None, index_col=0).reindex(truth.index).to_numpy(float)
            got = pd.read_csv(path, sep="\t", header=None, index_col=0).reindex(truth.index).to_numpy(float)
            return corr_score(ref, got, tol=1e-6)
        if op == "batch_correction":
            m = pd.read_csv(path, sep="\t", header=None, index_col=0).reindex(truth.index).to_numpy(float)
            lg = log2cpm(m) if "_counts" in Path(path).name else m
            bd, gd = effects(lg)
            c1 = float(np.corrcoef(bd[expressed], truth["batch_lfc"].to_numpy()[expressed])[0, 1])
            k1 = float(np.corrcoef(gd[expressed], truth["lfc"].to_numpy()[expressed])[0, 1])
            kept = min(1.0, max(0.0, k1) / k0)
            removed = 1 - max(0.0, c1) / c0
            return removed * kept, (f"batch removed {removed:.3f} (residual batch effect vs the simulated one: r {c1:.3f}, "
                                    f"{c0:.3f} before) x group effect kept {kept:.3f} (r with the simulated log2 FC {k1:.3f}, {k0:.3f} before)")
        return None, ""
    return score


# ── stability: the same data, presented differently ──────────────────────────

def perturbed_bed(src, dst, how, rng):
    df = pd.read_csv(src, sep="\t", header=None)
    if how == "shuffled":
        df = df.sample(frac=1.0, random_state=int(rng.integers(1 << 31)))
    elif how == "chrom-order":        # chromosomes grouped, but in another order; sorted within
        order = {c: i for i, c in enumerate(sorted(df[0].unique(), reverse=True))}
        df = df.assign(_k=df[0].map(order)).sort_values(["_k", 1]).drop(columns="_k")
    df.to_csv(dst, sep="\t", header=False, index=False)
    return dst


def part_stability(rec):
    print("\n== stability ==")
    rng = np.random.default_rng(0)
    od = WORK / "stability"
    od.mkdir(parents=True, exist_ok=True)
    n = 10_000
    a_bed, b_bed = str(DATA / f"peaks_A_{n}.bed"), str(DATA / f"peaks_B_{n}.bed")
    A, B = Loci.make(a_bed), Loci.make(b_bed)
    base = {"intersect": loci_set(A.take(A.overlap_any(B))), "merge": loci_set(A.merge())}
    # per-row answers are compared through the interval itself (the row order changes)
    nb, dist = A.nearest(B)
    base_near = {(str(c), int(s), int(e)): int(d) for c, s, e, d in zip(A.chroms, A.starts, A.ends, dist)}
    bt = tool("bedtools")
    for how in ("shuffled", "chrom-order"):
        pa = perturbed_bed(a_bed, str(od / f"A_{how}.bed"), how, rng)
        pb = perturbed_bed(b_bed, str(od / f"B_{how}.bed"), how, rng)
        A2, B2 = Loci.make(pa, sort=False), Loci.make(pb, sort=False)
        try:
            i2 = loci_set(A2.take(A2.overlap_any(B2))); m2 = loci_set(A2.merge())
            nb2, d2 = A2.nearest(B2)
            near2 = {(str(c), int(s), int(e)): int(d) for c, s, e, d in zip(A2.chroms, A2.starts, A2.ends, d2)}
            for task, ok in (("intersect", i2 == base["intersect"]), ("merge", m2 == base["merge"]), ("nearest", near2 == base_near)):
                rec.add(part="stability", theme="intervals", task=task, engine="genomeblocks", variant=how, stable=bool(ok))
        except Exception as e:
            for task in ("intersect", "merge", "nearest"):
                rec.add(part="stability", theme="intervals", task=task, engine="genomeblocks", variant=how, stable=False, error=str(e)[:200])
        # bedtools: intersect works on any order; merge requires sorted input
        for task, cmd in (("intersect", [bt, "intersect", "-u", "-a", pa, "-b", pb]), ("merge", [bt, "merge", "-i", pa])):
            try:
                o = cli(cmd)
                got = bed_set(pd.read_csv(io.StringIO(o), sep="\t", header=None))
                rec.add(part="stability", theme="intervals", task=task, engine="bedtools (CLI)", variant=how, stable=bool(got == base[task]))
            except Exception as e:
                rec.add(part="stability", theme="intervals", task=task, engine="bedtools (CLI)", variant=how, stable=False, error=str(e)[:200])
        try:
            o = cli([bt, "closest", "-t", "first", "-a", pa, "-b", pb])
            c = pd.read_csv(io.StringIO(o), sep="\t", header=None)
            k = c.shape[1] // 2
            g = np.where(c.iloc[:, k + 1] >= 0, gap(c.iloc[:, 1].to_numpy(), c.iloc[:, 2].to_numpy(), c.iloc[:, k + 1].to_numpy(), c.iloc[:, k + 2].to_numpy()), -1)
            near = {(str(r[0]), int(r[1]), int(r[2])): int(gg) for r, gg in zip(c.itertuples(index=False), g)}
            rec.add(part="stability", theme="intervals", task="nearest", engine="bedtools (CLI)", variant=how, stable=bool(near == base_near))
        except Exception as e:
            rec.add(part="stability", theme="intervals", task="nearest", engine="bedtools (CLI)", variant=how, stable=False, error=str(e)[:200])
        A2df = pd.read_csv(pa, sep="\t", header=None, usecols=[0, 1, 2])
        for kind in ("R", "py"):
            rows = worker(kind, "intervals", od / f"{kind}_{how}", pa, pb, "intersect,merge,nearest", env={"GB_REPEAT": "1"})
            for r in rows:
                if r.get("error") or not r.get("out"):
                    rec.add(part="stability", theme="intervals", task=r["op"], engine=r["engine"], variant=how, stable=False, error=str(r.get("error"))[:200])
                    continue
                if r["op"] in ("intersect", "merge"):
                    ok = bed_set(r["out"]) == base[r["op"]]
                else:
                    d = read_tsv(r["out"], ["idx", "bstart", "bend"])
                    idx = d["idx"].to_numpy(int)
                    g = gap(A2df.iloc[idx, 1].to_numpy(), A2df.iloc[idx, 2].to_numpy(), d["bstart"].to_numpy(), d["bend"].to_numpy())
                    near = {(str(A2df.iat[i, 0]), int(A2df.iat[i, 1]), int(A2df.iat[i, 2])): int(gg) for i, gg in zip(idx, g)}
                    ok = all(base_near.get(k) == v for k, v in near.items()) and len(near) == sum(1 for v in base_near.values() if v >= 0)
                rec.add(part="stability", theme="intervals", task=r["op"], engine=r["engine"], variant=how, stable=bool(ok))
        if PYBIN and (PYBIN / "bedops").exists():
            for task, cmd in (("intersect", [str(PYBIN / "bedops"), "--element-of", "1", pa, pb]), ("merge", [str(PYBIN / "bedops"), "--merge", pa])):
                try:
                    o = cli(cmd)
                    got = bed_set(pd.read_csv(io.StringIO(o), sep="\t", header=None))
                    rec.add(part="stability", theme="intervals", task=task, engine="BEDOPS (CLI)", variant=how, stable=bool(got == base[task]))
                except Exception as e:
                    rec.add(part="stability", theme="intervals", task=task, engine="BEDOPS (CLI)", variant=how, stable=False, error=str(e)[:200])
    # signal: shuffled regions must give the same rows
    bw = str(DATA / "signal_0.bw")
    L = signal_loci(2_000)
    S0 = L.signal([bw], n_bins=NBINS, flank=FLANK, verbose=False, progress=False, dtype=np.float64)[:, 0, :]
    perm = rng.permutation(len(L))
    Lp = L.take(perm)
    Sp = Lp.signal([bw], n_bins=NBINS, flank=FLANK, verbose=False, progress=False, dtype=np.float64)[:, 0, :]
    rec.add(part="stability", theme="signal", task="region_matrix", engine="genomeblocks", variant="shuffled", stable=bool(np.allclose(Sp, S0[perm])))
    pbed = write_bed(Lp, od / "signal_shuffled.bed")
    rows = worker("R", "signal", od / "r_signal", pbed, bw, FLANK, NBINS, env={"GB_REPEAT": "1"})
    for r in rows:
        if r["op"] != "region_matrix":
            continue
        if r.get("error") or not r.get("out"):
            rec.add(part="stability", theme="signal", task="region_matrix", engine=r["engine"], variant="shuffled", stable=False, error=str(r.get("error"))[:200])
            continue
        M = pd.read_csv(r["out"], sep="\t", header=None).to_numpy(float)
        # the canonical run of the same engine, for the engine's own convention
        base_rows = worker("R", "signal", od / "r_signal_base", write_bed(L, od / "signal_base.bed"), bw, FLANK, NBINS, env={"GB_REPEAT": "1"})
        b = next((x for x in base_rows if x["engine"] == r["engine"] and x.get("out")), None)
        if b is None:
            continue
        M0 = pd.read_csv(b["out"], sep="\t", header=None).to_numpy(float)
        rec.add(part="stability", theme="signal", task="region_matrix", engine=r["engine"], variant="shuffled",
                stable=bool(M.shape == M0.shape and np.allclose(M, M0[perm], equal_nan=True)))
    # annotation: shuffled peaks keep their classes and distances
    gtf = str(DATA / "genes.gtf")
    G = Genes.make(gtf)
    P = Loci.make(str(DATA / "peaks_A_10000.bed"))
    lab0 = G.labels(P); _, g0 = G.nearest_tss(P)
    perm = rng.permutation(len(P))
    Pp = P.take(perm)
    lab1 = G.labels(Pp); _, g1 = G.nearest_tss(Pp)
    rec.add(part="stability", theme="annotation", task="labels", engine="genomeblocks", variant="shuffled", stable=bool(np.array_equal(lab1, lab0[perm])))
    rec.add(part="stability", theme="annotation", task="nearest_tss", engine="genomeblocks", variant="shuffled", stable=bool(np.array_equal(g1, g0[perm])))
    pbed = write_bed(Pp, od / "peaks_shuffled.bed")
    obed = write_bed(P, od / "peaks_base.bed")
    r_base = worker("R", "annotation", od / "r_annot_base", obed, gtf, 1000, env={"GB_REPEAT": "1", "GB_REPEAT_SLOW": "1"})
    r_pert = worker("R", "annotation", od / "r_annot_pert", pbed, gtf, 1000, env={"GB_REPEAT": "1", "GB_REPEAT_SLOW": "1"})
    for r in r_pert:
        if r["op"] == "parse_gtf" or not r.get("out"):
            continue
        b = next((x for x in r_base if x["engine"] == r["engine"] and x["op"] == r["op"] and x.get("out")), None)
        if b is None:
            continue
        cols = ["idx", "label", "tss1"] if "ChIPseeker" in r["engine"] else (["idx", "tss1"] if "ChIPpeakAnno" in r["engine"] else (["idx", "label"] if r["op"] == "labels" else ["idx", "gap"]))
        d1 = read_tsv(r["out"], cols).set_index("idx").sort_index()
        d0 = read_tsv(b["out"], cols).set_index("idx").sort_index()
        # row i of the perturbed input is row perm[i] of the base input
        d1.index = perm[d1.index.to_numpy(int)]
        d1 = d1.sort_index()
        same = len(d0) == len(d1) and d0.index.equals(d1.index) and all((d0[c].to_numpy() == d1[c].to_numpy()).all() for c in cols[1:])
        rec.add(part="stability", theme="annotation", task=r["op"], engine=r["engine"], variant="shuffled", stable=bool(same))
    # graph: shuffled loops give the same overlapping set and components
    bedpe = str(DATA / "loops.bedpe")
    cre = Loci.make(str(DATA / "peaks_A_100000.bed"))
    reg = Loci.make(str(DATA / "peaks_B_100000.bed"))
    P0 = Pairs.make(bedpe)
    e0 = set(P0.overlapping(reg)["name"]); b0 = set(P0.overlapping(reg, both=True)["name"])
    lines = open(bedpe).read().splitlines()
    rng.shuffle(lines)
    sb = str(od / "loops_shuffled.bedpe")
    open(sb, "w").write("\n".join(lines) + "\n")
    P1 = Pairs.make(sb)
    rec.add(part="stability", theme="graph", task="loop_overlap_either", engine="genomeblocks", variant="shuffled", stable=bool(set(P1.overlapping(reg)["name"]) == e0))
    rec.add(part="stability", theme="graph", task="loop_overlap_both", engine="genomeblocks", variant="shuffled", stable=bool(set(P1.overlapping(reg, both=True)["name"]) == b0))
    A0 = quiet(lambda: Architecture.make(cre, bedpe, r=2500, verbose=False))()
    A1 = quiet(lambda: Architecture.make(cre, sb, r=2500, verbose=False))()

    def sizes2(m):
        m = np.asarray(m); s = np.bincount(m[m >= 0]); return np.sort(s[s >= 2])
    rec.add(part="stability", theme="graph", task="components", engine="genomeblocks", variant="shuffled",
            stable=bool(np.array_equal(sizes2(A0.components()), sizes2(A1.components()))))
    for typ in ("either", "both"):
        o = cli([bt, "pairtobed", "-a", sb, "-b", str(DATA / "peaks_B_100000.bed"), "-type", typ])
        ids = {line.split("\t")[6] for line in o.splitlines() if line}
        rec.add(part="stability", theme="graph", task=f"loop_overlap_{typ}", engine="bedtools pairtobed (CLI)", variant="shuffled", stable=bool(ids == (e0 if typ == "either" else b0)))
    names1 = list(P1["name"])
    rows = worker("R", "graph", od / "r_graph", sb, str(DATA / "peaks_B_100000.bed"), str(WORK / "graph" / "edges.tsv"), len(cre), env={"GB_REPEAT": "1"})
    for r in rows:
        if r["op"].startswith("loop_overlap") and r.get("out"):
            d = read_tsv(r["out"], ["idx"])
            got = {names1[i] for i in d["idx"].to_numpy(int)}
            rec.add(part="stability", theme="graph", task=r["op"], engine=r["engine"], variant="shuffled", stable=bool(got == (e0 if r["op"].endswith("either") else b0)))
    # expression: sample order and gene order must not matter
    size = os.environ.get("GB_BENCH_EXPR_SIZE", "60000_24")
    qdir = DATA / "expression" / f"salmon_{size}"
    Gx = Genes.make(str(DATA / "expression" / f"genes_{size.split('_')[0]}.gtf"))
    files = sorted(glob.glob(str(qdir / "*" / "quant.sf")))
    E0 = quiet(lambda: Expression.make(files, genes=Gx, level="transcript", verbose=False).to_genes(verbose=False))()
    sperm = rng.permutation(len(files))
    E1 = quiet(lambda: Expression.make([files[i] for i in sperm], genes=Gx, level="transcript", verbose=False).to_genes(verbose=False))()
    t0, t1 = E0.tmm(), E1.tmm()
    s0, s1 = E0.size_factors(), E1.size_factors()
    rec.add(part="stability", theme="expression", task="tmm", engine="genomeblocks", variant="samples shuffled", stable=bool(np.allclose(t1, t0[sperm])))
    rec.add(part="stability", theme="expression", task="size_factors", engine="genomeblocks", variant="samples shuffled", stable=bool(np.allclose(s1, s0[sperm])))
    c0 = pd.DataFrame(np.nan_to_num(E0.counts), columns=E0.sample_names)
    gperm = rng.permutation(len(c0))
    ex = od / "expression"
    ex.mkdir(exist_ok=True)
    c0.to_csv(ex / "base.tsv", sep="\t"); c0.iloc[gperm].to_csv(ex / "genes_shuffled.tsv", sep="\t")
    c0.iloc[:, sperm].to_csv(ex / "samples_shuffled.tsv", sep="\t")
    base_rows = worker("py", "expression", ex / "base", str(ex / "base.tsv"), str(qdir), str(WORK / "expression" / "tx2gene.tsv"), env={"GB_REPEAT": "1"})
    for variant, fn in (("genes shuffled", "genes_shuffled.tsv"), ("samples shuffled", "samples_shuffled.tsv")):
        rows = worker("py", "expression", ex / variant.replace(" ", "_"), str(ex / fn), str(qdir), str(WORK / "expression" / "tx2gene.tsv"), env={"GB_REPEAT": "1"})
        for r in rows:
            if r["op"] not in ("tmm", "size_factors") or not r.get("out"):
                continue
            b = next((x for x in base_rows if x["engine"] == r["engine"] and x["op"] == r["op"] and x.get("out")), None)
            if b is None:
                continue
            d0 = read_tsv(b["out"], ["sample", "factor"]).set_index("sample")["factor"]
            d1 = read_tsv(r["out"], ["sample", "factor"]).set_index("sample")["factor"].reindex(d0.index)
            rec.add(part="stability", theme="expression", task=r["op"], engine=r["engine"], variant=variant, stable=bool(np.allclose(d0.to_numpy(), d1.to_numpy())))


# ── hic: Hi-C map jobs added 2026-10-23 (theme graph3d) ──────────────────────
HICBIN = Path(os.environ["GB_LANDSCAPE_HIC"]) / "bin" if os.environ.get("GB_LANDSCAPE_HIC") else None
HIC_RES = 10_000


def hic_sim(n_chrom=2, length=40_000_000, total=15_000_000, seed=41):
    """A 10 kb Hi-C map with a known truth, as a .cool (cis only): P(s) ~ s^-1, TADs of 300 kb - 1.5 Mb (x2.5
    inside), A/B compartments of 1 - 5 Mb (x1.4 same, x0.7 opposite), a loop at the corner of every second TAD of
    at least 200 kb (+8x, Gaussian over ~1 bin), a log-normal bias per bin (2 % of bins empty). Returns the .cool
    path and the truth (boundaries, loops, compartment sign and bias per bin)."""
    import cooler
    d = WORK / "hic_sim"
    cool, truth_f = d / "sim.cool", d / "truth.json"
    if cool.exists() and truth_f.exists():
        return cool, json.loads(truth_f.read_text())
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    nb = length // HIC_RES
    bins, pix, truth = [], [], {"boundaries": [], "loops": [], "comp": [], "bias": [], "chroms": []}
    off = 0
    for k in range(n_chrom):
        chrom = f"chr{k + 1}"
        tad_edges = np.cumsum(np.r_[0, rng.integers(30, 151, nb)])
        tad_edges = tad_edges[tad_edges < nb]
        tad = np.searchsorted(tad_edges, np.arange(nb), side="right") - 1
        comp_edges = np.cumsum(np.r_[0, rng.integers(100, 501, nb)])
        comp = np.where((np.searchsorted(comp_edges, np.arange(nb), side="right") - 1) % 2 == 0, 1.0, -1.0)
        bias = rng.lognormal(0, 0.3, nb)
        bias[rng.random(nb) < 0.02] = 0.0
        i, j = np.meshgrid(np.arange(nb), np.arange(nb), indexing="ij")
        lam = (np.abs(i - j) + 1.0) ** -1.0
        lam *= np.where(tad[:, None] == tad[None, :], 2.5, 1.0)
        lam *= np.where(comp[:, None] == comp[None, :], 1.4, 0.7)
        loops = []
        ends = np.r_[tad_edges[1:], nb]
        for t, (a, b) in enumerate(zip(tad_edges, ends)):
            if t % 2 == 0 and b - a >= 20:
                la, lb = a + 1, b - 2
                loops.append((int(la), int(lb)))
                ii, jj = np.meshgrid(np.arange(la - 3, la + 4), np.arange(lb - 3, lb + 4), indexing="ij")
                ok = (ii >= 0) & (jj < nb)
                bump = 8 * np.exp(-((ii - la) ** 2 + (jj - lb) ** 2) / 2.0)
                lam[ii[ok], jj[ok]] *= 1 + bump[ok]
                lam[jj[ok], ii[ok]] *= 1 + bump[ok]
        lam *= bias[:, None] * bias[None, :]
        iu = np.triu_indices(nb)
        v = lam[iu]
        v *= (total / n_chrom) / v.sum()
        cnt = rng.poisson(v)
        keep = cnt > 0
        pix.append(pd.DataFrame({"bin1_id": iu[0][keep] + off, "bin2_id": iu[1][keep] + off, "count": cnt[keep].astype(np.int32)}))
        st = np.arange(nb) * HIC_RES
        bins.append(pd.DataFrame({"chrom": chrom, "start": st, "end": st + HIC_RES}))
        truth["chroms"].append(chrom)
        truth["boundaries"].append([int(x) for x in tad_edges[1:]])
        truth["loops"].append(loops)
        truth["comp"].append(comp.tolist())
        truth["bias"].append(bias.tolist())
        off += nb
        del lam, i, j
    B = pd.concat(bins, ignore_index=True)
    P = pd.concat(pix, ignore_index=True)
    cooler.create_cooler(str(cool), B, P, ordered=True, dtypes={"count": np.int32})
    # a phasing track for eigenvector sign (GC-like: the compartment sign plus noise)
    ph = B.assign(GC=np.concatenate([np.array(c) * 0.05 + 0.45 + rng.normal(0, 0.02, len(c)) for c in truth["comp"]]))
    ph.to_csv(d / "phasing.tsv", sep="\t", index=False)
    truth_f.write_text(json.dumps(truth))
    return cool, truth


def match_f1(called, true, tol=2):
    """F1 of called against true items (1-D bins or 2-D bin pairs), a match within tol bins (each true item once)."""
    called, true = [tuple(np.atleast_1d(c)) for c in called], [tuple(np.atleast_1d(t)) for t in true]
    used, tp = set(), 0
    for c in called:
        for k, t in enumerate(true):
            if k not in used and max(abs(a - b) for a, b in zip(c, t)) <= tol:
                used.add(k)
                tp += 1
                break
    prec = tp / max(1, len(called))
    rec_ = tp / max(1, len(true))
    f1 = 2 * prec * rec_ / max(1e-12, prec + rec_)
    return f1, f"F1 {f1:.3f} (precision {prec:.3f}, recall {rec_:.3f}; {len(called):,} called, {len(true):,} planted, within {tol} bins)"


def part_hic(rec):
    """Matrix balancing (C3), loop calling (C6), TAD / insulation boundaries (C7) and compartments (C8) on a
    simulated map with a known truth, competitors only: genomeblocks' loop callers run in part_hic_loops on the same map."""
    print("\n== hic ==")
    import cooler
    if not HICBIN:
        rec.add(part="hic", task="part", engine="—", seconds=None, runs=None, accuracy=None, error="GB_LANDSCAPE_HIC not set")
        return
    cool, truth = hic_sim()
    od = WORK / "hic"
    od.mkdir(parents=True, exist_ok=True)
    nbins = sum(len(b) for b in truth["bias"])
    bias = np.concatenate(truth["bias"])
    comp = np.concatenate(truth["comp"])
    hb = lambda x: str(HICBIN / x)
    env_hic = {**os.environ, "PATH": f"{HICBIN}:{os.environ['PATH']}"}

    def hcli(cmd):
        p = subprocess.run(cmd, capture_output=True, text=True, env=env_hic)
        if p.returncode != 0:
            raise RuntimeError(f"{Path(cmd[0]).name}: {p.stderr[-500:]}")
        return p.stdout

    def bias_score(w):
        w = np.asarray(w, float)
        ok = np.isfinite(w) & (w > 0) & (bias > 0)
        r = float(np.corrcoef(np.log(w[ok]), -np.log(bias[ok]))[0, 1])
        return max(0.0, r), f"Pearson {r:.4f} between log weight and minus log simulated bias ({ok.sum():,} of {nbins:,} bins)"
    # C3 balancing: each engine gets a fresh copy (copying is inside the clock for all three and takes < 0.1 s)
    work = str(od / "bal.cool")
    rec_cli(rec, "hic", "matrix_balancing", "cooler balance (CLI)", nbins,
            lambda: (shutil.copy(cool, work), hcli([hb("cooler"), "balance", "--nproc", "1", "--force", work])),
            lambda: bias_score(cooler.Cooler(work).bins()[:]["weight"].to_numpy()))
    if (HICBIN / "hictk").exists():
        work2 = str(od / "bal_hictk.cool")
        rec_cli(rec, "hic", "matrix_balancing", "hictk balance ice (CLI)", nbins,
                lambda: (shutil.copy(cool, work2), hcli([hb("hictk"), "balance", "ice", work2, "--threads", "1", "--force"])),
                lambda: bias_score(cooler.Cooler(work2).bins()[:]["weight"].to_numpy()))
    if (HICBIN / "hicCorrectMatrix").exists():
        out3 = str(od / "bal_hicexplorer.cool")

        def _hx_score():
            c = cooler.Cooler(out3)
            b = c.bins()[:]
            if "weight" in b.columns:
                # HiCExplorer stores divisive correction factors (corrected = raw / (w_i w_j)); cooler's convention multiplies
                return bias_score(1.0 / b["weight"].to_numpy())
            # corrected counts only: per bin, the median of log(raw / corrected) over its row is log b_i + const
            raw = cooler.Cooler(str(cool)).pixels()[:]
            cor = c.pixels()[:].merge(raw, on=["bin1_id", "bin2_id"], suffixes=("_c", "_r"))
            cor = cor[(cor["count_c"] > 0) & (cor["count_r"] > 0)]
            lr = np.log(cor["count_r"] / cor["count_c"])
            est = pd.concat([pd.Series(lr.to_numpy(), index=cor["bin1_id"]), pd.Series(lr.to_numpy(), index=cor["bin2_id"])]).groupby(level=0).median()
            w = np.full(nbins, np.nan)
            w[est.index.to_numpy()] = np.exp(-est.to_numpy())
            return bias_score(w)
        rec_cli(rec, "hic", "matrix_balancing", "HiCExplorer hicCorrectMatrix correct --correctionMethod ICE (CLI)", nbins,
                lambda: hcli([hb("hicCorrectMatrix"), "correct", "--matrix", str(cool), "--correctionMethod", "ICE", "--filterThreshold", "-3", "5",
                              "--outFileName", out3]), _hx_score, 1)
    # the callers below get a balanced map (balancing is the job above, done here outside the clock): cooltools a
    # cooler-balanced copy, HiCExplorer its own hicCorrectMatrix output (it reads weights as divisive factors)
    bal = str(od / "balanced.cool")
    shutil.copy(cool, bal)
    hcli([hb("cooler"), "balance", "--nproc", "1", "--force", bal])
    hx_bal = str(od / "balanced_hicexplorer.cool")
    if (HICBIN / "hicCorrectMatrix").exists() and not Path(hx_bal).exists():
        hcli([hb("hicCorrectMatrix"), "correct", "--matrix", str(cool), "--correctionMethod", "ICE", "--filterThreshold", "-3", "5", "--outFileName", hx_bal])
    true_loops = [(c, a, b) for c, ls in zip(truth["chroms"], truth["loops"]) for a, b in ls]
    true_bnd = [(c, x) for c, bs in zip(truth["chroms"], truth["boundaries"]) for x in bs]
    cidx = {c: k for k, c in enumerate(truth["chroms"])}

    def loop_score(df):
        """df: chrom1, start1, chrom2, start2 (bp) -> F1 against the planted loops, per chromosome (offset keeps them apart)."""
        called = [(cidx[str(c)] * 10 ** 6 + int(s1) // HIC_RES, int(s2) // HIC_RES) for c, s1, s2 in zip(df["c"], df["s1"], df["s2"]) if str(c) in cidx]
        return match_f1(called, [(cidx[c] * 10 ** 6 + a, b) for c, a, b in true_loops])

    def bnd_score(df):
        called = [cidx[str(c)] * 10 ** 6 + int(p) // HIC_RES for c, p in zip(df["c"], df["p"]) if str(c) in cidx]
        return match_f1(called, [cidx[c] * 10 ** 6 + x for c, x in true_bnd])
    if (HICBIN / "cooltools").exists():
        exp = str(od / "expected.tsv")
        dots = str(od / "dots.bedpe")
        rec_cli(rec, "hic", "loop_calling", "cooltools expected-cis + dots (CLI)", nbins,
                lambda: (hcli([hb("cooltools"), "expected-cis", "-p", "1", "-o", exp, bal]), hcli([hb("cooltools"), "dots", "-p", "1", "-o", dots, bal, exp])),
                lambda: loop_score(pd.read_csv(dots, sep="\t").rename(columns={"chrom1": "c", "start1": "s1", "start2": "s2"})), 1)
        ins = str(od / "insulation.tsv")
        rec_cli(rec, "hic", "tad_boundaries", "cooltools insulation 100 kb (CLI)", nbins,
                lambda: hcli([hb("cooltools"), "insulation", "-p", "1", "--threshold", "Li", "-o", ins, bal, "100000"]),
                lambda: (lambda d: bnd_score(d[d["is_boundary_100000"].fillna(False).astype(bool)].assign(p=lambda x: x["start"], c=lambda x: x["chrom"])))(pd.read_csv(ins, sep="\t")))
        eig = str(od / "eigs")
        rec_cli(rec, "hic", "compartments", "cooltools eigs-cis (CLI)", nbins,
                lambda: hcli([hb("cooltools"), "eigs-cis", "-o", eig, "--phasing-track", str(WORK / "hic_sim" / "phasing.tsv"), bal]),
                lambda: comp_score(pd.read_csv(eig + ".cis.vecs.tsv", sep="\t")["E1"].to_numpy(), comp))
    if (HICBIN / "hicDetectLoops").exists():
        lp = str(od / "hicexplorer_loops.bedgraph")
        # the parameters of HiCExplorer's own documentation example for a 10 kb map; no loops -> no file -> F1 0
        rec_cli(rec, "hic", "loop_calling", "HiCExplorer hicDetectLoops (CLI)", nbins,
                lambda: (Path(lp).unlink(missing_ok=True),
                         hcli([hb("hicDetectLoops"), "-m", hx_bal, "-o", lp, "--maxLoopDistance", "2000000", "--windowSize", "10", "--peakWidth", "6",
                               "--pValuePreselection", "0.05", "--pValue", "0.05", "-t", "1", "-tpc", "1"])),
                lambda: loop_score(pd.read_csv(lp, sep="\t", header=None).rename(columns={0: "c", 1: "s1", 4: "s2"}) if Path(lp).exists()
                                   else pd.DataFrame(columns=["c", "s1", "s2"])), 1)
    if (HICBIN / "hicFindTADs").exists():
        pre = str(od / "hicexplorer_tads")
        rec_cli(rec, "hic", "tad_boundaries", "HiCExplorer hicFindTADs (CLI)", nbins,
                lambda: hcli([hb("hicFindTADs"), "-m", hx_bal, "--outPrefix", pre, "--correctForMultipleTesting", "fdr", "--thresholdComparisons", "0.05",
                              "--minDepth", "30000", "--maxDepth", "100000", "--step", "10000", "-p", "1"]),
                lambda: bnd_score(pd.read_csv(pre + "_boundaries.bed", sep="\t", header=None).rename(columns={0: "c", 1: "p"}).assign(p=lambda x: x["p"] + HIC_RES // 2)), 1)
    if (HICBIN / "hicPCA").exists():
        pc = str(od / "hicexplorer_pc1.bedgraph")

        def _pca_score():
            d = pd.read_csv(pc, sep="\t", header=None, names=["c", "s", "e", "v"])
            v = np.full(nbins, np.nan)
            starts = {c: k for k, c in enumerate(truth["chroms"])}
            offs = np.cumsum([0] + [len(b) for b in truth["bias"]])
            for c, s_, x in zip(d["c"], d["s"], d["v"]):
                if str(c) in starts:
                    v[offs[starts[str(c)]] + int(s_) // HIC_RES] = x
            return comp_score(v, comp)
        rec_cli(rec, "hic", "compartments", "HiCExplorer hicPCA (CLI)", nbins,
                lambda: hcli([hb("hicPCA"), "-m", hx_bal, "-o", pc, "--whichEigenvectors", "1", "--format", "bedgraph"]), _pca_score, 1)


def comp_score(e1, comp, n_per_chrom=4000):
    """Mean over chromosomes of |Pearson| between an eigenvector and the simulated compartment sign (an eigenvector's
    sign is arbitrary per chromosome; only callers given a phasing track align them)."""
    e1 = np.asarray(e1, float)
    rs = []
    for a in range(0, len(e1), n_per_chrom):
        x, y = e1[a:a + n_per_chrom], comp[a:a + n_per_chrom]
        ok = np.isfinite(x)
        rs.append(abs(float(np.corrcoef(x[ok], y[ok])[0, 1])) if ok.sum() > 2 else 0.0)
    r = float(np.mean(rs))
    return r, f"|Pearson| {r:.4f} with the simulated A/B compartments, mean over {len(rs)} chromosomes (" + ", ".join(f"{x:.3f}" for x in rs) + ")"


def part_hic_loops(rec):
    """Loop calling (C6) on part_hic's simulated map by the callers genomeblocks has ports of: gb.mustache and
    gb.hiccups, and Mustache 1.3.3 (MUSTACHE: its command) and Juicer HiCCUPS (JUICER_TOOLS: juicer_tools_2.13.06.jar,
    JAVA: its java) when set. A part of its own, so it runs without part_hic's CLI environment and without replacing
    its rows; the map is the seeded hic_sim, ICE-balanced with cooler's defaults (part_hic's `cooler balance`)."""
    print("\n== hic_loops ==")
    import cooler
    cool, truth = hic_sim()
    od = WORK / "hic_loops"
    od.mkdir(parents=True, exist_ok=True)
    nbins = sum(len(b) for b in truth["bias"])
    bal = str(od / "balanced.cool")
    if not Path(bal).exists():
        shutil.copy(cool, bal)
        cooler.balance_cooler(cooler.Cooler(bal), store=True)
    cidx = {c: i for i, c in enumerate(truth["chroms"])}
    true_loops = [(c, a, b) for c, ls in zip(truth["chroms"], truth["loops"]) for a, b in ls]
    print(f"  map: {nbins:,} bins, {len(true_loops)} planted loops")

    def loop_score(df):
        called = [(cidx[str(c)] * 10 ** 6 + int(s1) // HIC_RES, int(s2) // HIC_RES) for c, s1, s2 in zip(df["c"], df["s1"], df["s2"]) if str(c) in cidx]
        return match_f1(called, [(cidx[c] * 10 ** 6 + a, b) for c, a, b in true_loops])

    # the two ports and the two originals on the same balanced 10 kb map, one core each like every engine of
    # part_hic (Mustache -p 1, Juicer --threads 1). Juicer reads a .hic built from the same pixels (pre, outside
    # the clock) and normalises with its own SCALE; the ports use the cooler's ICE weights.
    def loops_df(P):
        return pd.DataFrame({"c": P.a.chroms.astype(str), "s1": P.a.starts, "s2": P.b.starts})
    rec_py(rec, "hic_loops", "loop_calling", "genomeblocks mustache() (port of Mustache 1.3.3)", nbins,
           lambda: gb.mustache(bal, threads=1, verbose=False), lambda m: loop_score(loops_df(m.loops)), rep=1)
    rec_py(rec, "hic_loops", "loop_calling", "genomeblocks hiccups() (port of Juicer HiCCUPS, cooler weights)", nbins,
           lambda: gb.hiccups({HIC_RES: bal}, resolutions=(HIC_RES,), workers=1, verbose=False),
           lambda h: loop_score(loops_df(h.loops)), rep=1)
    if os.environ.get("MUSTACHE"):
        mo = str(od / "mustache_loops.tsv")
        rec_cli(rec, "hic_loops", "loop_calling", "Mustache 1.3.3 (CLI)", nbins,
                lambda: cli(shlex.split(os.environ["MUSTACHE"]) + ["-f", bal, "-r", str(HIC_RES), "-o", mo, "-p", "1"]),
                lambda: loop_score(pd.read_csv(mo, sep="\t").rename(columns={"BIN1_CHR": "c", "BIN1_START": "s1", "BIN2_START": "s2"})), 1)
    if os.environ.get("JUICER_TOOLS"):
        jt = shlex.split(os.environ.get("JAVA", "java")) + ["-Xmx16g", "-Djava.awt.headless=true", "-jar", os.environ["JUICER_TOOLS"]]
        hic, short, sizes_f = od / "sim.hic", od / "sim.short.txt", od / "sim.chrom.sizes"
        if not hic.exists():
            c = cooler.Cooler(str(cool))
            c.chromsizes.to_csv(sizes_f, sep="\t", header=False)
            px = c.pixels(join=True)[:]
            # Juicer's "short with score" format: str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 score (a position inside each bin)
            pd.DataFrame({"s1": 0, "c1": px["chrom1"], "p1": px["start1"] + HIC_RES // 2, "f1": 0, "s2": 0, "c2": px["chrom2"],
                          "p2": px["start2"] + HIC_RES // 2, "f2": 1, "n": px["count"]}).to_csv(short, sep=" ", header=False, index=False)
            subprocess.run(jt + ["pre", "-r", str(HIC_RES), "-k", "SCALE", str(short), str(hic), str(sizes_f)], check=True, capture_output=True)
        jdir = od / "juicer_loops"

        def _juicer_loops():
            ml = jdir / "merged_loops.bedpe"
            d = pd.read_csv(ml, sep="\t", comment="#", header=None) if ml.exists() and ml.stat().st_size else pd.DataFrame(columns=[0, 1, 4])
            # Juicer may drop the "chr" prefix of the names it was given
            d[0] = [x if str(x) in cidx else f"chr{x}" for x in d[0]]
            return loop_score(d.rename(columns={0: "c", 1: "s1", 4: "s2"}))
        rec_cli(rec, "hic_loops", "loop_calling", "Juicer juicer_tools 2.13.06 hiccups --cpu (CLI)", nbins,
                lambda: (shutil.rmtree(jdir, ignore_errors=True),
                         cli(jt + ["hiccups", "--cpu", "--threads", "1", "-r", str(HIC_RES), "-k", "SCALE", str(hic), str(jdir)])),
                _juicer_loops, 1)


PARTS = {"intervals": part_intervals, "annotation": part_annotation, "signal": part_signal, "signal_bigwig": part_signal_bigwig, "motifs": part_motifs,
         "cooccupancy": part_cooccupancy, "cooccupancy_bootstrap": part_cooccupancy_bootstrap, "graph": part_graph, "hic": part_hic, "hic_loops": part_hic_loops, "expression": part_expression, "stability": part_stability}

if __name__ == "__main__":
    todo = os.environ.get("GB_BENCH_PARTS")
    todo = todo.split(",") if todo else list(PARTS)
    # GB_BENCH_OUT=landscape_intervals writes results/landscape_intervals.json instead (one file per
    # SLURM job); landscape/merge_results.py folds those back into results/landscape.json
    rec = Recorder(os.environ.get("GB_BENCH_OUT", "landscape"))
    for name in todo:
        try:
            PARTS[name](rec)
        except Exception as e:                           # a part failing must not lose the others
            import traceback
            traceback.print_exc()
            rec.add(part=name, task="part", engine="—", seconds=None, runs=None, accuracy=None, error=f"{type(e).__name__}: {e}"[:400])
    rec.save(sizes=SIZES, rscript=RSCRIPT, python=PYTHON, giggle=GIGGLE, homer=str(HOMER_BIN))
