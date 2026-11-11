#!/usr/bin/env python3
"""HiChIP loops: MAPS (github.com/ijuric/MAPS) against genomeblocks.maps, on the same pairs and peaks.

  synthetic  data/hichip.allValidPairs (make_data.py hichip_avp) plus contacts at the loops of data/loops.bedpe,
             drawn as make_data.py hichip_cool draws them (same seed; the pairs of data/hichip_5kb.cool), written
             to a temporary allValidPairs; data/peaks_A_100000.bed as 1D peaks; 5 kb.
             MAPS runs when MAPS_HOME is set (a clone of github.com/ijuric/MAPS: bin/MAPS/*.py|r and
             MAPS_data_files/); MAPS_PYTHON / MAPS_RSCRIPT name a Python (pandas, numpy) and an Rscript (VGAM,
             MASS, data.table) for it (default: python3 / Rscript on PATH); give R the VGAM MAPS lists (1.1-2,
             e.g. through R_LIBS): VGAM >= 1.1-6 stops on deep samples (see genomeblocks/_maps.py). Its features file
             (MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_5Kb_el.GRCh38.txt) is used by both; without
             MAPS_HOME genomeblocks alone is timed on a made-up features table.
  real       GB_BENCH_MAPS_REAL=samples.tsv, one row per sample with the columns name, pairs, peaks, features,
             maps_dir (the output directory of a MAPS run on the same pairs, from feather-style files written
             as below), maps_name (that run's NAME, if not the sample name), maps_seconds (its wall time:
             MAPS.py + regression + formatting) and an optional note.
             genomeblocks is timed and compared with that output; no path is stored in the results.

MAPS starts from FASTQ (feather: bwa/chromap, MAPQ, duplicates, short/long split). Here both start from the same
pairs: MAPS gets the per-chromosome NAME.CHR.shrt.vip.bed / NAME.CHR.long.intra.bedpe feather would write, made
with feather's rule (see genomeblocks/_maps.py; long ends as 1-bp intervals at their 5' end), and runs
run_pipeline.sh's MAPS part (MAPS.py, MAPS_regression_and_peak_caller.r, MAPS_peak_formatting.r; hg38, chr1-22 +
chrX, binning_range 1 Mb, pospoisson). Writing those files is timed separately ("split"); genomeblocks' time
includes reading the pairs file. Compared: the tested bin pairs (set, counts, row order), the largest relative
difference of expected, P and fdr, the calls, and MAPS's sig3Dinteractions.bedpe line by line.
"""
from __future__ import annotations

import filecmp
import glob
import json
import os
import re
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from common import DATA, Recorder

AVP = DATA / "hichip.allValidPairs"
LOOPS = DATA / "loops.bedpe"
PEAKS = DATA / "peaks_A_100000.bed"
SIZES = DATA / "hg38.chrom.sizes"
MAPS_HOME = os.environ.get("MAPS_HOME")
MAPS_PY = os.environ.get("MAPS_PYTHON", "python3")
MAPS_R = os.environ.get("MAPS_RSCRIPT", "Rscript")
REAL = os.environ.get("GB_BENCH_MAPS_REAL")
FEATURES = "MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_5Kb_el.GRCh38.txt"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
LAYOUT = {"pairs": [0, 1, 2, 3, 4, 5, 6], "allvalidpairs": [0, 1, 2, 4, 5, 3, 6], "juicer": [0, 2, 3, 6, 7, 1, 5]}


def loop_pairs(out):
    """hichip.allValidPairs + Poisson(score) contacts at each loop of loops.bedpe, ends ~ anchor centre +-
    N(0, 1.5 kb) (make_data.py hichip_cool, seed 10), as an allValidPairs file."""
    import pandas as pd
    from common import read_chromsizes
    rng = np.random.default_rng(10)
    lp = pd.read_csv(LOOPS, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5, 7])
    k = rng.poisson(lp[7].to_numpy())
    i = np.repeat(np.arange(len(lp)), k)
    x = ((lp[1].to_numpy() + lp[2].to_numpy()) // 2)[i] + rng.normal(0, 1500, len(i)).astype(np.int64)
    y = ((lp[4].to_numpy() + lp[5].to_numpy()) // 2)[i] + rng.normal(0, 1500, len(i)).astype(np.int64)
    extra = pd.DataFrame({"id": [f"l{n}" for n in range(len(i))], "c1": lp[0].to_numpy()[i], "p1": np.maximum(x, 1),
                          "s1": "+", "c2": lp[3].to_numpy()[i], "p2": np.maximum(y, 1), "s2": "-"})
    avp = pd.read_csv(AVP, sep="\t", header=None, usecols=range(7), names=["id", "c1", "p1", "s1", "c2", "p2", "s2"])
    pairs = pd.concat([avp, extra], ignore_index=True)
    cs = read_chromsizes()
    pairs["p1"] = np.minimum(pairs.p1, pairs.c1.map(cs))
    pairs["p2"] = np.minimum(pairs.p2, pairs.c2.map(cs))
    pairs.to_csv(out, sep="\t", header=False, index=False)
    return out


def made_up_features(path):
    """A MAPS-style features table (5 kb bins of hg38) when MAPS's own is not at hand."""
    from common import read_chromsizes
    with open(path, "w") as f:
        for c, n in read_chromsizes(SIZES).items():
            if c not in CHROMS:
                continue
            b = np.arange(-(-n // 5000))
            el = 1000 + (b * 37) % 3000
            gc = 0.35 + 0.003 * ((b * 13) % 100)
            mp = 0.6 + 0.008 * ((b * 7) % 50)
            for i in range(len(b)):
                f.write(f"{c}\t{b[i] * 5000}\t{b[i] * 5000 + 5000}\t{el[i]}\t{gc[i]:.4f}\t{mp[i]:.4f}\n")


def time_gb(pairs, peaks, features, repeat=3):
    """Median wall time and peak memory of gb.maps in a fresh process (one warm-up run when repeat > 1). Peak
    memory is the process's own VmHWM."""
    code = f"""
import json, statistics, time
from genomeblocks._maps import maps
args = dict(features={features!r}, verbose=False)
if {repeat} > 1:
    maps({pairs!r}, {peaks!r}, **args)
ts = []
for _ in range({repeat}):
    t = time.perf_counter(); maps({pairs!r}, {peaks!r}, **args); ts.append(time.perf_counter() - t)
hwm = next(int(l.split()[1]) for l in open("/proc/self/status") if l.startswith("VmHWM"))
print(json.dumps({{"seconds": statistics.median(ts), "runs": ts, "peak_mb": hwm / 1024}}))
"""
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def split_pairs(pairs, layout, outdir, name):
    """feather's per-chromosome short / long files from a pairs file (the rule of genomeblocks._maps, written
    with plain pandas): NAME.CHR.shrt.vip.bed and NAME.CHR.long.intra.bedpe."""
    import pandas as pd
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    cols = LAYOUT[layout]
    names = ["id", "c1", "p1", "c2", "p2", "s1", "s2"]
    fs = {c: open(out / f"{name}.{c}.shrt.vip.bed", "w") for c in CHROMS}
    fl = {c: open(out / f"{name}.{c}.long.intra.bedpe", "w") for c in CHROMS}
    for ch in pd.read_csv(pairs, sep="\t", header=None, comment="#", usecols=cols, chunksize=5_000_000,
                          dtype={cols[0]: str, cols[1]: str, cols[3]: str, cols[2]: "int64", cols[4]: "int64",
                                 cols[5]: str, cols[6]: str}):
        ch = ch[cols]
        ch.columns = names
        ch = ch[(ch.c1 == ch.c2) & ch.c1.isin(CHROMS)]
        r1, r2 = ch.s1.isin(["-", "16"]), ch.s2.isin(["-", "16"])
        q1, q2 = ch.p1 - (~r1).astype("int64"), ch.p2 - (~r2).astype("int64")
        t = (q1 - q2).abs()
        s = (t <= 1000) & (r1 != r2) & (q1 != q2)
        lo, hi = np.minimum(q1[s], q2[s]), np.maximum(q1[s], q2[s])
        sh = pd.DataFrame({"c": ch.c1[s], "a": lo, "b": hi, "id": ch.id[s], "e": ""})
        g = t > 1000
        lg = pd.DataFrame({"c": ch.c1[g], "a1": q1[g], "b1": q1[g] + 1, "c2": ch.c1[g], "a2": q2[g], "b2": q2[g] + 1,
                           "id": ch.id[g]})
        for c, d in sh.groupby("c", sort=False):
            d.to_csv(fs[c], sep="\t", header=False, index=False)
        for c, d in lg.groupby("c", sort=False):
            d.to_csv(fl[c], sep="\t", header=False, index=False)
    for f in list(fs.values()) + list(fl.values()):
        f.close()


def run_maps(feather_dir, name, peaks, features, out):
    """run_pipeline.sh's MAPS part on feather-style files; wall time and peak memory of the three steps."""
    out.mkdir(parents=True, exist_ok=True)
    o, fd = str(out) + "/", str(feather_dir) + "/"
    b = Path(MAPS_HOME) / "bin" / "MAPS"
    subprocess.run([MAPS_PY, str(b / "make_maps_runfile.py"), name, o, str(peaks), str(features), fd, fd, "5000",
                    "22", o, "X", "--BINNING_RANGE", "1000000"], check=True, capture_output=True)
    t = time.perf_counter()
    with open(out / "maps.log", "w") as log:
        for cmd in ([MAPS_PY, str(b / "MAPS.py"), f"{o}maps_{name}.maps"],
                    [MAPS_R, str(b / "MAPS_regression_and_peak_caller.r"), o, f"{name}.5k", "5000", "22X", "None",
                     "pospoisson"],
                    [MAPS_R, str(b / "MAPS_peak_formatting.r"), o, f"{name}.5k", "2", "5000"]):
            if "formatting" in cmd[1] and not (out / f"{name}.5k.2.peaks").exists():
                break                                   # MAPS stops (quit()) when no bin pair passes the cutoffs
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True, cwd=str(b))
    return {"seconds": time.perf_counter() - t, "peak_mb": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024}


def compare(res, ref, name, tmp):
    """Agreement of a gb.maps result with a MAPS output directory."""
    import pandas as pd
    T = res.tested
    G = pd.DataFrame({"chr": T.a.chroms.astype(str), "bin1_mid": T.a.starts, "bin2_mid": T.b.starts,
                      "set": T.cols["set"], "count_gb": T.cols["count"], "expected_gb": T.cols["expected"],
                      "P_gb": T.cols["P"], "fdr_gb": T.cols["fdr"], "ratio_gb": T.cols["ratio"]})
    parts = []
    for f in sorted(glob.glob(f"{ref}/reg_raw.*.{name}.5k.*.MAPS2_pospoisson")):
        d = pd.read_csv(f, sep=" ", float_precision="round_trip",
                        usecols=["bin1_mid", "bin2_mid", "count", "chr", "expected2", "p_val_reg2", "fdr"])
        d["set"] = re.search(r"\.(and|xor)\.MAPS2", f).group(1).upper()
        parts.append(d)
    R = pd.concat(parts, ignore_index=True)
    R[["bin1_mid", "bin2_mid"]] = R[["bin1_mid", "bin2_mid"]].astype(np.int64)
    row = {"tested_ref": len(R), "tested_gb": len(G)}
    same = all(np.array_equal(d[["bin1_mid", "bin2_mid"]].to_numpy(),
                              G[(G.chr == c) & (G.set == s)][["bin1_mid", "bin2_mid"]].to_numpy())
               for (c, s), d in R.groupby(["chr", "set"], sort=False))
    m = R.merge(G, on=["chr", "bin1_mid", "bin2_mid", "set"], how="inner")
    rel = lambda a, b: float(np.max(np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)))  # noqa
    with np.errstate(divide="ignore", invalid="ignore"):
        call = lambda c, r, q: (c >= 12) & (r >= 2) & (-np.log10(q) > 2)       # noqa: E731
        cr = call(m["count"], m["count"] / m["expected2"], m["fdr"])
        cg = call(m["count_gb"], m["ratio_gb"], m["fdr_gb"])
    row.update(tested_shared=len(m), same_row_order=bool(same), counts_equal=bool((m["count"] == m.count_gb).all()),
               max_rel_expected=rel(m.expected2.to_numpy(), m.expected_gb.to_numpy()),
               max_rel_P=rel(m.p_val_reg2.to_numpy(), m.P_gb.to_numpy()),
               max_rel_fdr=rel(m.fdr.to_numpy(), m.fdr_gb.to_numpy()),
               calls_ref=int(cr.sum()), calls_gb=int(cg.sum()), calls_shared=int((cr & cg).sum()))
    mine = res.write(str(tmp), f"{name}.5k")["bedpe"]
    theirs = f"{ref}/{name}.5k.2.sig3Dinteractions.bedpe"
    a = open(mine).read().splitlines()
    b = open(theirs).read().splitlines() if os.path.exists(theirs) else a[:1]     # MAPS stopped: no calls
    key = lambda L: {tuple(x.split("\t")[:6]) for x in L[1:]}                   # noqa: E731
    row.update(loops_ref=len(b) - 1, loops_gb=len(a) - 1, loops_shared=len(key(a) & key(b)),
               bedpe_identical=os.path.exists(theirs) and filecmp.cmp(mine, theirs, shallow=False),
               bedpe_lines_differing=sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b)),
               bedpe_lines_differing_beyond_fdr_digits=sum(
                   [v for k, v in enumerate(x.split("\t")) if k not in (8, 12)]
                   != [v for k, v in enumerate(y.split("\t")) if k not in (8, 12)] for x, y in zip(a, b)))
    return row


if __name__ == "__main__":
    from genomeblocks._maps import maps
    rec = Recorder("maps")
    parts = os.environ.get("GB_BENCH_PARTS", "synthetic,real").split(",")
    if "synthetic" in parts and AVP.exists() and LOOPS.exists():
        keep = tempfile.TemporaryDirectory(dir=DATA)
        tmp = Path(keep.name)
        PAIRS = loop_pairs(tmp / "hichip_loops.allValidPairs")
        if MAPS_HOME:
            feats = str(Path(MAPS_HOME) / FEATURES)
            t = time.perf_counter()
            split_pairs(PAIRS, "allvalidpairs", tmp / "feather", "bench")
            rec.add(part="synthetic", impl="split (feather files for MAPS)", seconds=time.perf_counter() - t)
            mp = run_maps(tmp / "feather", "bench", PEAKS, feats, tmp / "maps")
        else:
            print("[skip] MAPS: set MAPS_HOME (and MAPS_PYTHON / MAPS_RSCRIPT) to compare")
            feats = str(tmp / "features.txt")
            made_up_features(feats)
        rec.add(part="synthetic", impl="genomeblocks", **time_gb(str(PAIRS), str(PEAKS), feats))
        if MAPS_HOME:
            rec.add(part="synthetic", impl="MAPS", **mp)
            res = maps(str(PAIRS), str(PEAKS), features=feats, verbose=True)
            rec.add(part="synthetic", impl="agreement", **compare(res, tmp / "maps", "bench", tmp / "gb"))
        keep.cleanup()
    elif "synthetic" in parts:
        print(f"[skip] {AVP} or {LOOPS} missing: python make_data.py loops hichip_avp")
    if "real" in parts and REAL:
        import pandas as pd
        for s in pd.read_table(REAL).fillna("").itertuples():
            t = time_gb(s.pairs, s.peaks, s.features, repeat=1)
            res = maps(s.pairs, s.peaks, features=s.features, verbose=True)
            with tempfile.TemporaryDirectory(dir=DATA) as d:
                agree = compare(res, s.maps_dir, getattr(s, "maps_name", "") or s.name, Path(d))
            rec.add(part="real", sample=s.name, seconds=t["seconds"], peak_mb=t["peak_mb"],
                    maps_seconds=float(s.maps_seconds), note=str(getattr(s, "note", "")), **agree)
    rec.save()
