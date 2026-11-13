#!/usr/bin/env python3
"""Fold results/landscape_<part>.json files (one per SLURM job) into results/landscape.json.

    python landscape/merge_results.py

Rows of a part present in a part file replace that part's rows in landscape.json;
every part file's environment (the node it ran on) is kept under meta.envs, and
machine paths in error texts and metadata are cut to portable names.
"""
import glob
import json
import re
from pathlib import Path

# machine paths in error texts and run metadata -> short, portable names (results hold no machine paths)
_PATHS = [(re.compile(r"/[^\s'\"]*?/micromamba/envs/([^/\s'\"]+)/"), r"\1/"),
          (re.compile(r"/[^\s'\"]*?/benchmarks/data/"), "data/"),
          (re.compile(r"/(?:home|groups|scratch|tmp)/[^\s'\":]*/"), "")]


def scrub(x):
    if isinstance(x, str):
        for pat, rep in _PATHS:
            x = pat.sub(rep, x)
        return x
    if isinstance(x, dict):
        return {k: scrub(v) for k, v in x.items()}
    if isinstance(x, list):
        return [scrub(v) for v in x]
    return x

RES = Path(__file__).resolve().parent.parent / "results"
main = RES / "landscape.json"
out = json.loads(main.read_text()) if main.exists() else {"bench": "landscape", "env": None, "meta": {}, "rows": []}
envs = out["meta"].get("envs", {})
for f in sorted(glob.glob(str(RES / "landscape_*.json"))):
    if f.endswith("landscape_scores.json"):
        continue
    d = json.loads(Path(f).read_text())
    parts = {r.get("part") for r in d["rows"]}
    out["rows"] = [r for r in out["rows"] if r.get("part") not in parts] + d["rows"]
    for p in parts:
        envs[p] = {k: d["env"].get(k) for k in ("cpu", "cores", "mem_gb", "python")}
    out["env"] = out["env"] or d["env"]
    out["meta"].update({k: v for k, v in d.get("meta", {}).items() if k != "envs"})
    print(f"{Path(f).name}: {len(d['rows'])} rows, parts {sorted(p for p in parts if p)}")
out["meta"]["envs"] = envs
main.write_text(json.dumps(scrub(out), indent=1, default=float))
print(f"[merged] {main}: {len(out['rows'])} rows")
