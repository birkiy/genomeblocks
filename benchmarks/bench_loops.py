#!/usr/bin/env python3
"""HiChIP loops: FitHiChIP 11.0 against genomeblocks.fithichip, on the same cooler and peaks.

  synthetic  data/hichip_5kb.cool (make_data.py hichip_cool: the HiChIP pairs + contacts at the 50k loops of
             loops.bedpe) with data/peaks_A_100000.bed as anchors. FitHiChIP runs when FITHICHIP_HOME is set
             (the directory holding FitHiChIP_HiCPro.sh, as its runner sees it); FITHICHIP_EXEC is a command
             prefix for the runner (e.g. "singularity exec -B /data fithichip.sif"), else FitHiChIP runs on the
             host with its R / python packages. Without FitHiChIP only genomeblocks is timed.
  real       GB_BENCH_LOOPS_REAL=samples.tsv, one row per sample with the columns name, cool (a cooler URI),
             peaks, chromsizes (or -), fithichip_dir (the FitHiC_BiasCorr directory of a FitHiChIP run on the
             same cooler and peaks), fithichip_seconds (its wall time, e.g. from the pipeline trace) and an
             optional note on that run (e.g. its IntType). genomeblocks is timed and compared with that output;
             no path is stored in the results.

Compared, wherever the FitHiChIP output has them: the tested pairs (set, counts), the largest relative
difference of P and q, the loop and merged-loop sets, and FitHiChIP's two IGV bedpe files byte for byte.
FitHiChIP's config: IntType=3 BINSIZE=5000 LowDistThr=5000 UppDistThr=10000000 QVALUE=0.01 BiasType=1
UseP2PBackgrnd=0 MergeInt=1 (FitHiChIP's defaults otherwise), the input as nf-core/hic's FITHICHIP_COOL2MATRIX
writes it (cooler dump, 1-based bin ids).
"""
from __future__ import annotations

import filecmp
import glob
import json
import os
import resource
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder

COOL = DATA / "hichip_5kb.cool"
PEAKS = DATA / "peaks_A_100000.bed"
SIZES = DATA / "hg38.chrom.sizes"
FH_HOME = os.environ.get("FITHICHIP_HOME")
FH_EXEC = shlex.split(os.environ.get("FITHICHIP_EXEC", ""))
REAL = os.environ.get("GB_BENCH_LOOPS_REAL")
SUB = "FitHiChIP_Peak2ALL_b5000_L5000_U10000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr"


def time_gb(cool, peaks, sizes, repeat=3):
    """Median wall time and peak memory of gb.fithichip in a fresh process (one warm-up run). Peak memory is
    the process's own VmHWM: ru_maxrss would carry over the parent's resident size across fork + exec."""
    code = f"""
import json, statistics, time
import genomeblocks as gb
args = dict(chromsizes={sizes!r}, verbose=False)
gb.fithichip({cool!r}, {peaks!r}, **args)
ts = []
for _ in range({repeat}):
    t = time.perf_counter(); gb.fithichip({cool!r}, {peaks!r}, **args); ts.append(time.perf_counter() - t)
hwm = next(int(l.split()[1]) for l in open("/proc/self/status") if l.startswith("VmHWM"))
print(json.dumps({{"seconds": statistics.median(ts), "runs": ts, "peak_mb": hwm / 1024}}))
"""
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def run_fithichip(cool, peaks, sizes, out):
    """FitHiChIP_HiCPro.sh on the cooler's raw 5 kb counts; wall time and peak memory of the run."""
    import cooler
    clr = cooler.Cooler(str(cool))
    out.mkdir(parents=True, exist_ok=True)
    b = clr.bins()[:][["chrom", "start", "end"]]
    b["id"] = np.arange(1, len(b) + 1)
    b.to_csv(out / "abs.bed", sep="\t", header=False, index=False)
    with open(out / "matrix", "w") as f:
        for lo in range(0, clr.info["nnz"], 20_000_000):
            px = clr.pixels()[lo:lo + 20_000_000]
            np.savetxt(f, np.column_stack([px.bin1_id + 1, px.bin2_id + 1, px["count"]]), fmt="%d", delimiter="\t")
    cfg = out / "fithichip.config"
    cfg.write_text(f"Interval={out / 'abs.bed'}\nMatrix={out / 'matrix'}\nPeakFile={peaks}\nChrSizeFile={sizes}\n"
                   f"OutDir={out / 'fithichip'}\nPREFIX=bench\nIntType=3\nBINSIZE=5000\nLowDistThr=5000\n"
                   f"UppDistThr=10000000\nQVALUE=0.01\nUseP2PBackgrnd=0\nBiasType=1\nMergeInt=1\nOverWrite=1\n")
    t = time.perf_counter()
    with open(out / "fithichip.log", "w") as log:
        subprocess.run(FH_EXEC + ["bash", "-c", f"cd {shlex.quote(FH_HOME)} && bash FitHiChIP_HiCPro.sh -C {cfg}"],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    sec = time.perf_counter() - t
    return {"seconds": sec, "peak_mb": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024}, \
        out / "fithichip" / SUB


def compare(res, ref, tmp):
    """Agreement of a gb.fithichip result with a FitHiChIP FitHiC_BiasCorr directory."""
    import pandas as pd
    paths = res.write(str(tmp), "gb")
    one = lambda pat: (glob.glob(f"{ref}/{pat}") or [None])[0]
    key = lambda P: set(zip(P.a.chroms.astype(str).tolist(), P.a.starts.tolist(), P.b.starts.tolist()))
    row = {"tested_gb": len(res.tested), "loops_gb": len(res.loops), "merged_gb": len(res.merged)}
    igv, migv = one("*_FitHiC_Q0.01_IGV.bedpe"), one("Merge_Nearby_Interactions/*_MergeNearContacts_IGV.bedpe")
    S = pd.read_table(one("*_FitHiC_Q0.01.bed"), usecols=[0, 1, 4])
    sr = set(map(tuple, S.to_numpy().tolist()))
    M = pd.read_table(one("Merge_Nearby_Interactions/*_MergeNearContacts.bed"), usecols=[0, 1, 4])
    mr = set(map(tuple, M.to_numpy().tolist()))
    row.update(loops_ref=len(sr), loops_shared=len(sr & key(res.loops)), merged_ref=len(mr),
               merged_shared=len(mr & key(res.merged)))
    for name, mine, theirs in (("igv", paths["igv"], igv), ("merged_igv", paths["merged_igv"], migv)):
        a, b = open(mine).read().splitlines(), open(theirs).read().splitlines()
        row[f"{name}_identical"] = filecmp.cmp(mine, theirs, shallow=False)
        row[f"{name}_lines"] = len(b)
        row[f"{name}_lines_differing"] = sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))
    tested = one("*.interactions_FitHiC.bed")
    if tested:
        R = pd.read_table(tested, usecols=["chr1", "s1", "s2", "cc", "P-Value_Bias", "Q-Value_Bias"],
                          float_precision="round_trip")
        G = pd.DataFrame({"chr1": res.tested.a.chroms.astype(str), "s1": res.tested.a.starts,
                          "s2": res.tested.b.starts, "cc_gb": res.tested.cols["count"], "P_gb": res.tested.cols["P"],
                          "Q_gb": res.tested.cols["Q"]})
        m = R.merge(G, on=["chr1", "s1", "s2"], how="inner")
        rel = lambda a, b: float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-300)))
        row.update(tested_ref=len(R), tested_shared=len(m), counts_equal=bool((m.cc == m.cc_gb).all()),
                   max_rel_P=rel(m["P-Value_Bias"].to_numpy(), m.P_gb.to_numpy()),
                   max_rel_Q=rel(m["Q-Value_Bias"].to_numpy(), m.Q_gb.to_numpy()))
    return row


if __name__ == "__main__":
    import genomeblocks as gb
    rec = Recorder("loops")
    parts = os.environ.get("GB_BENCH_PARTS", "synthetic,real").split(",")
    if "synthetic" in parts and COOL.exists():
        fh = None
        if FH_HOME:                                   # first, so the children's peak memory is FitHiChIP's
            keep = tempfile.TemporaryDirectory(dir=DATA)
            fh, ref = run_fithichip(COOL, PEAKS, SIZES, Path(keep.name))
        else:
            print("[skip] FitHiChIP: set FITHICHIP_HOME (and FITHICHIP_EXEC) to compare")
        rec.add(part="synthetic", impl="genomeblocks", **time_gb(str(COOL), str(PEAKS), str(SIZES)))
        if fh:
            res = gb.fithichip(str(COOL), str(PEAKS), chromsizes=str(SIZES), verbose=True)
            rec.add(part="synthetic", impl="FitHiChIP", **fh)
            rec.add(part="synthetic", impl="agreement", **compare(res, ref, Path(keep.name) / "gb"))
            keep.cleanup()
    elif "synthetic" in parts:
        print(f"[skip] {COOL} missing: python make_data.py hichip_cool")
    if "real" in parts and REAL:
        import pandas as pd
        for s in pd.read_table(REAL).fillna("").itertuples():
            sizes = None if s.chromsizes == "-" else s.chromsizes
            t = time_gb(s.cool, s.peaks, sizes, repeat=1)
            res = gb.fithichip(s.cool, s.peaks, chromsizes=sizes, verbose=True)
            with tempfile.TemporaryDirectory(dir=DATA) as d:
                agree = compare(res, s.fithichip_dir, Path(d))
            rec.add(part="real", sample=s.name, seconds=t["seconds"], peak_mb=t["peak_mb"],
                    fithichip_seconds=float(s.fithichip_seconds), note=str(getattr(s, "note", "")), **agree)
    rec.save()
