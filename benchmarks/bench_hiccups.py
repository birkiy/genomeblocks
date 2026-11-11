#!/usr/bin/env python3
"""Hi-C loops: Juicer's HiCCUPS (juicer_tools 2.13.06 ``hiccups --cpu``) against genomeblocks' port, same matrix.

  synthetic  two chromosomes (40 + 26 Mb) of Hi-C pairs built here: distance decay, per-kb bias, a 300 kb gap
             without contacts and 190 planted loops; written as Juicer short format and as an .mcool (5/10/25 kb,
             Juicer's binning pos // res, cooler-balanced), cached in $GB_BENCH_DATA/hiccups_synth. With
             JUICER_TOOLS (the juicer_tools_2.13.06.jar; JAVA for the java to run it, default ``java``) the
             pairs also go through ``pre -r 5000,10000,25000 -k SCALE``, ``hiccups --cpu --threads
             $GB_BENCH_THREADS -k SCALE`` (Juicer's defaults otherwise) and ``dump norm / expected SCALE``;
             genomeblocks then runs on Juicer's SCALE vectors and expected and is compared file by file.
             Without Juicer only genomeblocks is timed (cooler weights).
  real       GB_BENCH_HICCUPS_REAL=samples.tsv, one row per sample with the columns name, mcool (5/10/25 kb raw
             counts binned like the .hic), juicer_dir (holding loops/ — Juicer's output directory — and
             vectors/norm.<chrom>.<res>, vectors/expected.<chrom>.<res> from ``dump``), juicer_seconds (wall time
             of the hiccups run), juicer_threads and an optional note. genomeblocks runs on the dumped vectors;
             no path is stored in the results.

Compared, per resolution and for the merged loops: the pixel sets, every number column (largest relative
difference), and the files byte for byte (Juicer writes the enriched pixels of a multi-threaded run in thread order:
their rows are compared as sets). With more than one thread Juicer's own fdr* columns change from run to run:
GPUOutputContainer.updateHistograms is synchronized on a per-tile object, so threads lose histogram counts
(the thresholds survive, the q-values move by up to a few %); ``*_columns_differing`` lists what differs.
"""
from __future__ import annotations

import filecmp
import json
import os
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder

SYN = DATA / "hiccups_synth"
JT = os.environ.get("JUICER_TOOLS")
JAVA = shlex.split(os.environ.get("JAVA", "java"))
THREADS = int(os.environ.get("GB_BENCH_THREADS", "1"))
REAL = os.environ.get("GB_BENCH_HICCUPS_REAL")
RES = (5000, 10000, 25000)


def make_synthetic(out: Path, seed: int = 7):
    """Juicer short-format pairs + chrom sizes + an .mcool of the same pairs (pos // res, like Juicer's pre)."""
    import cooler
    import pandas as pd
    out.mkdir(parents=True, exist_ok=True)
    if (out / "synth.mcool").exists():
        return
    rng = np.random.default_rng(seed)
    sizes = {"chr1": 40_000_000, "chr2": 26_000_001}
    (out / "synth.sizes").write_text("".join(f"{c}\t{s}\n" for c, s in sizes.items()))
    frames = []
    for c, L in sizes.items():
        n = int(8e6 * L / sum(sizes.values()))
        nb = L // 1000 + 1
        bias = rng.lognormal(0, 0.35, nb)
        gap = rng.integers(nb // 4, nb // 2)
        bias[gap:gap + 300] = 0
        b1 = rng.choice(nb, size=2 * n, p=bias / bias.sum())
        p1 = b1 * 1000 + rng.integers(0, 1000, 2 * n)
        d = (np.exp(rng.uniform(np.log(500), np.log(L), 2 * n)) * rng.random(2 * n) ** 0.25).astype(np.int64)
        p2 = p1 + d
        ok = (p2 < L) & (rng.random(2 * n) < bias[np.minimum(p2 // 1000, nb - 1)] / bias.max())
        frames.append(pd.DataFrame({"c": c, "p1": p1[ok][:n], "p2": p2[ok][:n]}))
        na = 120 if c == "chr1" else 70
        a1 = rng.integers(300_000, L - 4_000_000, na)
        a2 = a1 + np.exp(rng.uniform(np.log(60_000), np.log(3_000_000), na)).astype(np.int64)
        k = rng.integers(25, 150, na)
        q1 = np.repeat(a1, k) + rng.normal(0, 3000, k.sum()).astype(np.int64)
        q2 = np.repeat(a2, k) + rng.normal(0, 3000, k.sum()).astype(np.int64)
        m = (q1 > 0) & (q2 < L) & (q1 < q2)
        frames.append(pd.DataFrame({"c": c, "p1": q1[m], "p2": q2[m]}))
    df = pd.concat(frames).sort_values(["c", "p1"], kind="stable")
    df["p1"] = df["p1"].clip(lower=1)
    short = pd.DataFrame({"s1": rng.choice([0, 16], len(df)), "c1": df.c, "p1": df.p1, "f1": 0,
                          "s2": rng.choice([0, 16], len(df)), "c2": df.c, "p2": df.p2, "f2": 1})
    short.to_csv(out / "synth.short.txt", sep=" ", header=False, index=False)
    bins = cooler.binnify(pd.Series(sizes), 5000)
    px = pd.DataFrame({"bin1_id": _bin_id(bins, df.c.to_numpy(), df.p1.to_numpy() // 5000),
                       "bin2_id": _bin_id(bins, df.c.to_numpy(), df.p2.to_numpy() // 5000)})
    px = px.value_counts().rename("count").reset_index().sort_values(["bin1_id", "bin2_id"])
    cooler.create_cooler(str(out / "synth.5000.cool"), bins, px)
    cooler.zoomify_cooler(str(out / "synth.5000.cool"), str(out / "synth.mcool"), list(RES), chunksize=10_000_000)
    for r in RES:
        cooler.balance_cooler(cooler.Cooler(f"{out / 'synth.mcool'}::resolutions/{r}"), store=True)


def _bin_id(bins, chrom, idx):
    first = {c: int(i) for c, i in zip(*np.unique(bins.chrom.to_numpy(), return_index=True))}
    names, inv = np.unique(chrom, return_inverse=True)
    return np.array([first[c] for c in names])[inv] + idx


def run_juicer(out: Path, hic: Path, sizes: Path, threads: int):
    """juicer_tools hiccups --cpu (defaults, SCALE) and the vectors it used; wall time of hiccups."""
    jt = JAVA + ["-Xmx32g", "-Djava.awt.headless=true", "-Ddevelopment=false", "-jar", JT]
    if not hic.exists():
        subprocess.run(jt + ["pre", "-r", "5000,10000,25000", "-k", "SCALE", "--threads", str(threads),
                             str(out / "synth.short.txt"), str(hic), str(sizes)], check=True, capture_output=True)
    jdir = out / f"juicer_t{threads}"
    (jdir / "vectors").mkdir(parents=True, exist_ok=True)
    t = time.perf_counter()
    subprocess.run(jt + ["hiccups", "--cpu", "--threads", str(threads), "-k", "SCALE", str(hic), str(jdir / "loops")],
                   check=True, capture_output=True)
    sec = time.perf_counter() - t
    for res in RES:
        for c in [l.split()[0] for l in open(sizes)]:
            for what in ("norm", "expected"):
                subprocess.run(jt + ["dump", what, "SCALE", str(hic), c, "BP", str(res),
                                     str(jdir / "vectors" / f"{what}.{c}.{res}")], check=True, capture_output=True)
    return jdir, sec


def juicer_vectors(mcool, jdir):
    """Juicer's SCALE vectors over the cooler's bins and its expected per chromosome, per resolution."""
    import cooler
    read = lambda p: np.array([float(x) for x in open(p).read().split()]) if p.exists() else None
    norms, exps = {}, {}
    for res in RES:
        clr = cooler.Cooler(f"{mcool}::resolutions/{res}")
        parts, e = [], {}
        for c, n in zip(clr.chromnames, np.diff(clr._load_dset("indexes/chrom_offset"))):
            v = read(Path(jdir) / "vectors" / f"norm.{c}.{res}")
            parts.append(v[:n] if v is not None and len(v) >= n else np.full(n, np.nan))
            x = read(Path(jdir) / "vectors" / f"expected.{c}.{res}")
            if x is not None:
                e[c] = x
        norms[res], exps[res] = np.concatenate(parts), e
    return norms, exps


def time_gb(mcool, kwargs_code: str, threads: int):
    """Wall time and peak memory of genomeblocks' hiccups in a fresh process (vectors loaded outside the timing)."""
    code = f"""
import json, sys, time
sys.path.insert(0, {str(Path(__file__).parent)!r})
from bench_hiccups import juicer_vectors
from genomeblocks._hiccups import hiccups
{kwargs_code}
t = time.perf_counter(); hiccups({str(mcool)!r}, workers={threads}, verbose=False, **kw); s = time.perf_counter() - t
hwm = next(int(l.split()[1]) for l in open("/proc/self/status") if l.startswith("VmHWM"))
print(json.dumps({{"seconds": s, "peak_mb": hwm / 1024}}))
"""
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def _rows(path):
    lines = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]
    head = open(path).readline().lstrip("#").rstrip("\n").split("\t")
    return {(r[0], int(r[1]), int(r[4]), int(r[2]) - int(r[1])): np.array(r[11:], float) for r in lines}, head[11:]


def compare(ref_dir: Path, got_dir: Path) -> dict:
    """Agreement of genomeblocks' files with Juicer's output directory."""
    row = {}
    names = ["merged_loops.bedpe"] + [f"{k}_{r}.bedpe" for r in RES for k in ("postprocessed_pixels",
                                                                            "enriched_pixels")]
    for name in names:
        ref, got = ref_dir / name, got_dir / name
        if not ref.exists():
            continue
        a, cols = _rows(ref)
        b, _ = _rows(got)
        shared = sorted(set(a) & set(b))
        key = name.replace(".bedpe", "")
        row[f"{key}_juicer"], row[f"{key}_gb"], row[f"{key}_shared"] = len(a), len(b), len(shared)
        if shared:
            x, y = np.array([a[k] for k in shared]), np.array([b[k] for k in shared])
            with np.errstate(divide="ignore", invalid="ignore"):
                d = np.where(x == y, 0.0, np.abs(x - y) / np.maximum(np.abs(x), np.abs(y)))
            row[f"{key}_max_rel_diff"] = float(np.nanmax(d))
            # which value columns differ at all: a multi-threaded Juicer run differs from itself in the fdr*
            # columns only (its (chunk, observed) histograms are updated without a lock, so counts get lost)
            row[f"{key}_columns_differing"] = [c for c, dc in zip(cols, (d > 0).any(axis=0)) if dc]
        row[f"{key}_identical"] = filecmp.cmp(ref, got, shallow=False)
        if not row[f"{key}_identical"] and name.startswith("enriched"):
            ra = sorted(open(ref).read().splitlines()[2:])
            rb = sorted(open(got).read().splitlines()[2:])
            row[f"{key}_identical_as_sets"] = ra == rb
    for r in RES:
        f = f"fdr_thresholds_{r}"
        if (ref_dir / f).exists():
            row[f"{f}_identical"] = filecmp.cmp(ref_dir / f, got_dir / f, shallow=False)
    return row


def exact_run(mcool, jdir, out: Path):
    from genomeblocks._hiccups import hiccups
    norms, exps = juicer_vectors(mcool, jdir)
    res = hiccups(str(mcool), weights=norms, divisive=True, expected=exps, workers=THREADS, verbose=True)
    res.write(str(out))
    return compare(Path(jdir) / "loops", out)


KW_JUICER = "norms, exps = juicer_vectors({m!r}, {j!r}); kw = dict(weights=norms, divisive=True, expected=exps)"

if __name__ == "__main__":
    rec = Recorder("hiccups")
    parts = os.environ.get("GB_BENCH_PARTS", "synthetic,real").split(",")
    if "synthetic" in parts:
        make_synthetic(SYN)
        mcool = SYN / "synth.mcool"
        rec.add(part="synthetic", impl="genomeblocks (cooler weights)", threads=THREADS,
                **time_gb(mcool, "kw = {}", THREADS))
        if JT:
            jdir, sec = run_juicer(SYN, SYN / "synth.hic", SYN / "synth.sizes", THREADS)
            rec.add(part="synthetic", impl="juicer_tools 2.13.06 hiccups --cpu", threads=THREADS, seconds=sec)
            rec.add(part="synthetic", impl="genomeblocks (Juicer's SCALE + expected)", threads=THREADS,
                    **time_gb(mcool, KW_JUICER.format(m=str(mcool), j=str(jdir)), THREADS))
            with tempfile.TemporaryDirectory(dir=DATA) as d:
                rec.add(part="synthetic", impl="agreement", **exact_run(mcool, jdir, Path(d)))
        else:
            print("[skip] Juicer: set JUICER_TOOLS=juicer_tools_2.13.06.jar (and JAVA) to compare")
    if "real" in parts and REAL:
        import pandas as pd
        for s in pd.read_table(REAL).fillna("").itertuples():
            t = time_gb(s.mcool, KW_JUICER.format(m=s.mcool, j=s.juicer_dir), THREADS)
            with tempfile.TemporaryDirectory(dir=DATA) as d:
                agree = exact_run(s.mcool, s.juicer_dir, Path(d))
            rec.add(part="real", sample=s.name, threads=THREADS, seconds=t["seconds"], peak_mb=t["peak_mb"],
                    juicer_seconds=float(s.juicer_seconds), juicer_threads=int(s.juicer_threads),
                    note=str(getattr(s, "note", "")), **agree)
    rec.save()
