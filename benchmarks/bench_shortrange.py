#!/usr/bin/env python3
"""HiChIP short-range track: shell recipe vs genomeblocks, same input, same output.

Input: data/hichip.allValidPairs (4M HiC-Pro pairs: 3M H3K27ac-like short-range
+ 1M Hi-C background). Steps compared:

  1. cis pairs <= 1 kb -> both 5' ends (BED6)      awk          | hichip.shortrange_ends
  2. ends -> 147 bp fragments -> sorted -> coverage  awk|sort|bedtools genomecov -bg
                                                                 | hichip.fragments + coverage
  3. bedGraph -> bigWig                              (same writer for both: pybigtools)
  4. MACS3 on the ends (same tool either way; timed once)

Outputs are compared byte-for-byte after sorting (ends) and as records (bedGraph).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time

from common import DATA, Recorder

PY = sys.executable
AVP = str(DATA / "hichip.allValidPairs")
CS = str(DATA / "hg38.chrom.sizes")
MACS3 = os.path.join(os.path.dirname(PY), "macs3")

# extract_shortrange.sh, as described: cis pairs with |pos2 - pos1| <= 1 kb -> both 5' ends as BED6
EXTRACT = r"""awk 'BEGIN{OFS="\t"} $2==$5 { d=$6-$3; if (d<0) d=-d;
  if (d<=1000) { print $2,$3-1,$3,$1"/1",0,$4; print $5,$6-1,$6,$1"/2",0,$7 } }' %s > %s"""
# step 3 of make_tracks.sh, unchanged except sort memory/threads for a 4-core, 16 GB machine
COVER = r"""awk 'BEGIN{OFS="\t"} { if($6=="+"){s=$2; e=$2+147} else {e=$3; s=$3-147}
   if(s<0)s=0; print $1,s,e }' %s \
 | sort -k1,1 -k2,2n -S 4G --parallel=4 -T %s \
 | bedtools genomecov -bg -i - -g <(sort -k1,1 %s) > %s"""

PEAK = """
def PEAK():
    for line in open("/proc/self/status"):
        if line.startswith("VmHWM:"):
            return int(line.split()[1]) / 1024
"""


def run_shell(cmd):
    code = (f"import resource, subprocess, time, json\nt=time.perf_counter()\n"
            f"subprocess.run(['bash','-c',{cmd!r}], check=True)\n"
            f"print(json.dumps({{'s': time.perf_counter()-t, "
            f"'mb': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024}}))")
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def run_gb(tmp):
    code = PEAK + f"""
import json, time
T = {{}}; t0 = t = time.perf_counter()
from genomeblocks.columnar import hichip
sizes = {{l.split()[0]: int(l.split()[1]) for l in open({CS!r}) if l.strip()}}
ends = hichip.shortrange_ends({AVP!r}, 1000); T["ends"] = time.perf_counter() - t; t = time.perf_counter()
hichip.write_bed(ends, {tmp!r} + "/gb_ends.bed"); T["write BED"] = time.perf_counter() - t; t = time.perf_counter()
fr = hichip.fragments(ends, 147, sizes)
hichip.to_bedgraph(fr, {tmp!r} + "/gb.bdg"); T["coverage (bedGraph)"] = time.perf_counter() - t; t = time.perf_counter()
hichip.to_bigwig(fr, {tmp!r} + "/gb.bw", sizes); T["bigWig"] = time.perf_counter() - t
T["total"] = time.perf_counter() - t0
print(json.dumps({{"steps": T, "mb": PEAK(), "n_ends": len(ends)}}))
"""
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True, check=True,
                       env=dict(os.environ, PYTHONPATH=str(DATA.parent.parent)))
    return json.loads(r.stdout.strip().splitlines()[-1])


def bdg_to_bw(bdg, bw):
    code = f"""
import pybigtools, time, json
t = time.perf_counter()
sizes = {{l.split()[0]: int(l.split()[1]) for l in open({CS!r}) if l.strip()}}
def rec():
    for line in open({bdg!r}):
        c, s, e, v = line.split()
        yield c, int(s), int(e), float(v)
pybigtools.open({bw!r}, "w").write(dict(sorted(sizes.items())), rec())
print(json.dumps({{"s": time.perf_counter() - t}}))
"""
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


if __name__ == "__main__":
    rec = Recorder("shortrange")
    with tempfile.TemporaryDirectory(dir=str(DATA)) as tmp:
        sh = {}
        sh["extract ends (awk)"] = run_shell(EXTRACT % (AVP, f"{tmp}/sh_ends.bed"))
        sh["fragments + sort + genomecov"] = run_shell(COVER % (f"{tmp}/sh_ends.bed", tmp, CS, f"{tmp}/sh.bdg"))
        sh["bedGraph -> bigWig"] = bdg_to_bw(f"{tmp}/sh.bdg", f"{tmp}/sh.bw")
        for k, v in sh.items():
            rec.add(part="steps", impl="shell", step=k, seconds=v["s"], peak_mb=v.get("mb"))
        rec.add(part="total", impl="shell", seconds=sum(v["s"] for v in sh.values()),
                peak_mb=max(v.get("mb") or 0 for v in sh.values()))
        gb = run_gb(tmp)
        for k, v in gb["steps"].items():
            if k != "total":
                rec.add(part="steps", impl="genomeblocks", step=k, seconds=v)
        rec.add(part="total", impl="genomeblocks", seconds=gb["steps"]["total"], peak_mb=gb["mb"])
        # same output?
        same_ends = subprocess.run(
            ["bash", "-c", f"cmp -s <(cut -f1-3,6 {tmp}/sh_ends.bed | LC_ALL=C sort) "
                           f"<(cut -f1-3,6 {tmp}/gb_ends.bed | LC_ALL=C sort) && echo same || echo differ"],
            capture_output=True, text=True).stdout.strip()
        r1 = subprocess.run(["bash", "-c", f"cmp -s {tmp}/sh.bdg {tmp}/gb.bdg && echo same || echo differ"],
                            capture_output=True, text=True).stdout.strip()
        n_bdg = int(subprocess.run(["bash", "-c", f"wc -l < {tmp}/sh.bdg"], capture_output=True, text=True).stdout)
        t = time.perf_counter()
        subprocess.run([MACS3, "callpeak", "-t", f"{tmp}/gb_ends.bed", "-f", "BED", "-g", "hs", "-n", "x",
                        "--outdir", tmp, "--nomodel", "--extsize", "147", "-q", "0.01", "--keep-dup", "all"],
                       check=True, capture_output=True)
        macs = time.perf_counter() - t
        n_peaks = sum(1 for _ in open(f"{tmp}/x_peaks.narrowPeak"))
        rec.add(part="check", ends=same_ends, bedgraph=r1, bedgraph_records=n_bdg, n_ends=gb["n_ends"],
                macs3_seconds=macs, macs3_peaks=n_peaks, seconds=macs)
    rec.save()
