#!/usr/bin/env python3
"""Generate the synthetic-but-realistic datasets used by the benchmarks.

Everything is seeded and laid out on the real hg38 assembly (chrom sizes ship
with ``bioframe``). Motifs are the real JASPAR CORE vertebrate collection
(shipped with ``pyjaspar``). Nothing is downloaded at run time.

Outputs (in ``$GB_BENCH_DATA`` or ``benchmarks/data``):

  hg38.chrom.sizes          chr1-22, chrX
  peaks_A_<n>.bed           "ATAC-like" peak sets (n = 1k..1M), clustered
  peaks_B_<n>.bed           "ChIP-like" peak sets, partially co-located with A
  signal_<k>.bw             genome-wide bigWigs, ~13M variable-span records each
  atlas/track_<i>.bed       500 peak files for the Atlas (GIGGLE-style) index
  atlas_query.bed           20k-peak query set for the Atlas
  genes.gtf                 GENCODE-shaped GTF (20k genes, ~70k transcripts)
  hic.pairs                 5M Hi-C read pairs (4DN .pairs) with distance decay
  loops.bedpe               50k chromatin loops anchored near A peaks
  loops_<n>.bedpe           12.5k / 200k / 1M loops (scaling)
  loops_trans.bedpe         loops.bedpe + 2.5k inter-chromosomal loops (5%)
  hic_trans_5kb.cool        hic.pairs + ~40 read pairs at each inter-chromosomal loop
  hichip.allValidPairs      HiC-Pro style H3K27ac HiChIP pairs: 3M short-range pairs over
                            H3K27ac-like CREs (with SE-like clusters) + 1M Hi-C background
  genome.fa                 chr21 + chr22, random sequence at 41% GC
  jaspar.txt                JASPAR CORE vertebrates (jaspar16 format)

Usage:
    python make_data.py            # everything
    python make_data.py peaks bw   # just some parts
"""
from __future__ import annotations

import sys
from pathlib import Path
import time

import numpy as np
import pandas as pd

from common import DATA, read_chromsizes

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
PEAK_SIZES = [1_000, 10_000, 100_000, 1_000_000]


def _tic(msg):
    print(f"[make_data] {msg} ...", flush=True)
    return time.perf_counter()


def _toc(t0):
    print(f"            done in {time.perf_counter() - t0:.1f}s", flush=True)


# ── chrom sizes ──────────────────────────────────────────────────────────────

def chromsizes():
    import bioframe
    info = bioframe.assembly_info("hg38").seqinfo.set_index("name")["length"]
    with open(DATA / "hg38.chrom.sizes", "w") as f:
        for c in CHROMS:
            f.write(f"{c}\t{int(info[c])}\n")


# ── peaks ────────────────────────────────────────────────────────────────────

def _hotspots(cs: dict, n: int, rng):
    """Shared 'regulatory hotspot' centres: peaks cluster in gene-dense areas."""
    names = list(cs)
    sizes = np.array([cs[c] for c in names], dtype=float)
    ci = rng.choice(len(names), size=n, p=sizes / sizes.sum())
    pos = (rng.random(n) * sizes[ci]).astype(np.int64)
    return ci, pos


def make_peaks(n: int, cs: dict, rng, *, hot=None, hot_frac=0.6,
               scatter=40_000, median_w=400):
    """``n`` peaks: ``hot_frac`` scattered around hotspots, rest uniform.

    Widths are lognormal (median ``median_w``), clipped to 100-5000 bp —
    the shape of MACS2 narrowPeak / ATAC calls.
    """
    names = list(cs)
    sizes = np.array([cs[c] for c in names], dtype=np.int64)
    n_hot = int(n * hot_frac)
    if hot is None:
        hot = _hotspots(cs, 4000, rng)
    h_ci, h_pos = hot
    pick = rng.integers(0, len(h_ci), n_hot)
    ci_hot = h_ci[pick]
    pos_hot = h_pos[pick] + rng.normal(0, scatter, n_hot).astype(np.int64)
    ci_uni = rng.choice(len(names), size=n - n_hot,
                        p=sizes / sizes.sum())
    pos_uni = (rng.random(n - n_hot) * sizes[ci_uni]).astype(np.int64)
    ci = np.concatenate([ci_hot, ci_uni])
    centre = np.concatenate([pos_hot, pos_uni])
    w = np.clip(rng.lognormal(np.log(median_w), 0.6, n), 100, 5000).astype(np.int64)
    start = np.clip(centre - w // 2, 0, sizes[ci] - w - 1)
    df = pd.DataFrame({"chrom": np.array(names)[ci], "start": start,
                       "end": start + w})
    order = np.lexsort((df["start"].to_numpy(), ci))
    return df.iloc[order].reset_index(drop=True)


def write_bed(df: pd.DataFrame, path, name_prefix="p"):
    out = df[["chrom", "start", "end"]].copy()
    out["name"] = [f"{name_prefix}{i}" for i in range(len(out))]
    out["score"] = 0
    out["strand"] = "."
    out.to_csv(path, sep="\t", header=False, index=False)


def peaks():
    cs = read_chromsizes()
    rng = np.random.default_rng(1)
    hot = _hotspots(cs, 4000, rng)
    for n in PEAK_SIZES:
        a = make_peaks(n, cs, np.random.default_rng(10 + n), hot=hot)
        # a TF binds mostly in open chromatin: 40% of B re-uses A sites
        # (jittered, own width), the rest is B-specific
        rng = np.random.default_rng(20 + n)
        n_sh = int(0.4 * n)
        sh = a.sample(n_sh, random_state=30 + n).copy()
        mid = (sh["start"] + sh["end"]) // 2 + rng.integers(-150, 150, n_sh)
        w = np.clip(rng.lognormal(np.log(600), 0.6, n_sh), 100, 5000).astype(np.int64)
        sh["start"] = np.maximum(mid - w // 2, 0)
        sh["end"] = sh["start"] + w
        b = pd.concat([sh, make_peaks(n - n_sh, cs, rng, hot=hot, median_w=600)])
        b = b.sort_values(["chrom", "start"], kind="stable").reset_index(drop=True)
        write_bed(a, DATA / f"peaks_A_{n}.bed", "a")
        write_bed(b, DATA / f"peaks_B_{n}.bed", "b")
        # unsorted union A+B, the input of the merge benchmark
        ab = pd.concat([a, b]).sample(frac=1, random_state=40 + n)
        write_bed(ab, DATA / f"peaks_AB_{n}.bed", "ab")


# ── bigWigs ──────────────────────────────────────────────────────────────────

def _merge_sorted(s, e):
    """Merge overlapping sorted intervals (numpy)."""
    if len(s) == 0:
        return s, e
    run_end = np.maximum.accumulate(e)
    new = np.empty(len(s), bool)
    new[0] = True
    new[1:] = s[1:] > run_end[:-1]
    grp = np.cumsum(new) - 1
    ms = s[new]
    me = np.zeros(len(ms), dtype=e.dtype)
    np.maximum.at(me, grp, e)
    return ms, me


def _explode(starts, lengths):
    total = int(lengths.sum())
    offs = np.repeat(np.cumsum(lengths) - lengths, lengths)
    return np.repeat(starts, lengths) + (np.arange(total) - offs)


def _chrom_records(size: int, pk_s, pk_e, rng):
    """Variable-span bedGraph records for one chromosome.

    Background: read 'islands' (~1 per 2.5 kb, mean 250 bp) split into
    segments of mean 25 bp with low RPM-like values — what bedtools
    genomecov produces from sparse background reads. Peaks: merged peak
    intervals rendered as 10 bp steps of a Gaussian bump. Zero-coverage
    stretches are omitted, as in real ChIP/ATAC bigWigs.
    """
    n_isl = size // 2500
    isl_s = np.unique(rng.integers(0, size - 2000, n_isl))
    isl_len = np.clip(rng.geometric(1 / 250, len(isl_s)), 20, 1500)
    isl_e = np.minimum(isl_s + isl_len, np.append(isl_s[1:], size))
    isl_len = isl_e - isl_s
    keep = isl_len > 0
    isl_s, isl_e, isl_len = isl_s[keep], isl_e[keep], isl_len[keep]
    # segment islands: breakpoints in the concatenated island coordinate space
    cum = np.cumsum(isl_len)
    total = int(cum[-1])
    bps = np.unique(np.concatenate([
        rng.integers(1, total, total // 25), cum[:-1], [0]]))
    seg_cs = bps
    seg_ce = np.append(bps[1:], total)
    isl_idx = np.searchsorted(cum, seg_cs, side="right")
    base = isl_s[isl_idx] - (cum[isl_idx] - isl_len[isl_idx])
    bg_s, bg_e = seg_cs + base, seg_ce + base
    bg_v = rng.gamma(1.5, 0.25, len(bg_s))
    # drop background segments that touch a peak
    if len(pk_s):
        j = np.searchsorted(pk_s, bg_e, side="left") - 1
        hit = (j >= 0) & (pk_e[np.maximum(j, 0)] > bg_s)
        bg_s, bg_e, bg_v = bg_s[~hit], bg_e[~hit], bg_v[~hit]
        # peaks as 10 bp steps of a Gaussian bump
        step = 10
        n_steps = np.maximum((pk_e - pk_s) // step, 1)
        idx = np.repeat(np.arange(len(pk_s)), n_steps)
        k = np.arange(int(n_steps.sum())) - np.repeat(np.cumsum(n_steps) - n_steps, n_steps)
        ps = pk_s[idx] + k * step
        pe = np.minimum(ps + step, pk_e[idx])
        pe[np.r_[np.cumsum(n_steps) - 1]] = pk_e  # last step reaches peak end
        amp = rng.lognormal(1.5, 0.8, len(pk_s))
        mid = (pk_s + pk_e) / 2
        sd = np.maximum((pk_e - pk_s) / 4, 1)
        z = ((ps + pe) / 2 - mid[idx]) / sd[idx]
        pv = amp[idx] * np.exp(-0.5 * z * z) + 0.3
        s = np.concatenate([bg_s, ps]); e = np.concatenate([bg_e, pe])
        v = np.concatenate([bg_v, pv])
    else:
        s, e, v = bg_s, bg_e, bg_v
    o = np.argsort(s, kind="stable")
    s, e, v = s[o], e[o], v[o]
    ok = e > s
    return s[ok], e[ok], v[ok]


def bigwigs(n_files=4):
    import pyBigWig
    cs = read_chromsizes()
    universe = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None,
                           usecols=[0, 1, 2], names=["chrom", "start", "end"])
    for k in range(n_files):
        rng = np.random.default_rng(100 + k)
        # each track: 40k peaks — 75% from the shared universe, 25% private
        shared = universe.sample(30_000, random_state=100 + k)
        private = make_peaks(10_000, cs, rng)
        pk = pd.concat([shared, private]).sort_values(["chrom", "start"])
        path = DATA / f"signal_{k}.bw"
        bw = pyBigWig.open(str(path), "w")
        bw.addHeader([(c, cs[c]) for c in CHROMS], maxZooms=10)
        n_rec = 0
        for c in CHROMS:
            sub = pk[pk["chrom"] == c]
            ms, me = _merge_sorted(sub["start"].to_numpy(np.int64),
                                   sub["end"].to_numpy(np.int64))
            s, e, v = _chrom_records(cs[c], ms, me, rng)
            for lo in range(0, len(s), 1_000_000):
                hi = lo + 1_000_000
                bw.addEntries([c] * len(s[lo:hi]), s[lo:hi], ends=e[lo:hi],
                              values=v[lo:hi].astype(np.float64))
            n_rec += len(s)
        bw.close()
        print(f"            {path.name}: {n_rec:,} records, "
              f"{path.stat().st_size / 1e6:.0f} MB", flush=True)


# ── Atlas tracks ─────────────────────────────────────────────────────────────

def atlas(n_tracks=500):
    cs = read_chromsizes()
    d = DATA / "atlas"
    d.mkdir(exist_ok=True)
    universe = make_peaks(300_000, cs, np.random.default_rng(7))
    rng = np.random.default_rng(8)
    sizes = np.clip(rng.lognormal(np.log(12_000), 0.7, n_tracks), 1_000, 60_000).astype(int)
    for i, n in enumerate(sizes):
        frac = rng.uniform(0.2, 0.8)       # how "CRE-like" this factor is
        n_u = int(n * frac)
        u = universe.sample(n_u, random_state=int(rng.integers(1 << 31)))
        r = make_peaks(n - n_u, cs, rng)
        df = pd.concat([u, r]).sort_values(["chrom", "start"])
        write_bed(df, d / f"track_{i:04d}.bed", f"t{i}_")
    q = universe.sample(20_000, random_state=9).sort_values(["chrom", "start"])
    write_bed(q, DATA / "atlas_query.bed", "q")


# ── GTF ──────────────────────────────────────────────────────────────────────

def gtf(n_genes=20_000):
    cs = read_chromsizes()
    rng = np.random.default_rng(3)
    names = list(cs)
    sizes = np.array([cs[c] for c in names], dtype=float)
    lines = []
    for g in range(n_genes):
        c = names[rng.choice(len(names), p=sizes / sizes.sum())]
        glen = int(np.clip(rng.lognormal(np.log(25_000), 1.0), 2_000, 1_000_000))
        gs = int(rng.integers(1, cs[c] - glen - 1))
        ge = gs + glen
        strand = "+" if rng.random() < 0.5 else "-"
        gid = f"ENSG{g:011d}.1"
        gname = f"GENE{g}"
        gtype = "protein_coding" if rng.random() < 0.6 else "lncRNA"
        ga = f'gene_id "{gid}"; gene_type "{gtype}"; gene_name "{gname}";'
        lines.append(f"{c}\tSYN\tgene\t{gs}\t{ge}\t.\t{strand}\t.\t{ga}")
        for t in range(1 + rng.poisson(2.5)):
            # isoforms: alternative TSS / TES inside the gene span
            ts = gs + int(rng.integers(0, max(1, glen // 3))) if t else gs
            te = ge - int(rng.integers(0, max(1, glen // 3))) if t else ge
            if te - ts < 500:
                ts, te = gs, ge
            tid = f"ENST{g:09d}{t:02d}.1"
            ta = f'gene_id "{gid}"; transcript_id "{tid}"; gene_type "{gtype}"; gene_name "{gname}";'
            lines.append(f"{c}\tSYN\ttranscript\t{ts}\t{te}\t.\t{strand}\t.\t{ta}")
            n_ex = 1 + rng.poisson(6)
            inner = np.sort(rng.integers(ts + 1, te - 1, max(0, 2 * (n_ex - 1))))
            bounds = np.concatenate([[ts], inner, [te]])
            ex = [(int(bounds[2 * i]), int(bounds[2 * i + 1])) for i in range(len(bounds) // 2)]
            ex = [(a, min(b, a + int(rng.lognormal(np.log(150), 0.6)) + 50)) for a, b in ex[:-1]] + [ex[-1]]
            order = ex if strand == "+" else ex[::-1]
            coding = gtype == "protein_coding" and len(ex) > 1
            for k, (a, b) in enumerate(order, 1):
                ea = f'{ta} exon_number {k};'
                lines.append(f"{c}\tSYN\texon\t{a}\t{b}\t.\t{strand}\t.\t{ea}")
                if not coding:
                    continue
                if k == 1:
                    lines.append(f"{c}\tSYN\tfive_prime_UTR\t{a}\t{b}\t.\t{strand}\t.\t{ea}")
                elif k == len(order):
                    lines.append(f"{c}\tSYN\tthree_prime_UTR\t{a}\t{b}\t.\t{strand}\t.\t{ea}")
                else:
                    lines.append(f"{c}\tSYN\tCDS\t{a}\t{b}\t.\t{strand}\t0\t{ea}")
    with open(DATA / "genes.gtf", "w") as f:
        f.write("##description: synthetic GENCODE-shaped annotation\n")
        f.write("\n".join(lines) + "\n")
    print(f"            genes.gtf: {len(lines):,} lines", flush=True)


# ── Hi-C pairs + loops ───────────────────────────────────────────────────────

def pairs(n=5_000_000):
    cs = read_chromsizes()
    rng = np.random.default_rng(4)
    names = np.array(list(cs))
    sizes = np.array([cs[c] for c in names], dtype=np.int64)
    p = sizes / sizes.sum()
    path = DATA / "hic.pairs"
    with open(path, "w") as f:
        f.write("## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        step = 1_000_000
        for lo in range(0, n, step):
            m = min(step, n - lo)
            c1 = rng.choice(len(names), size=m, p=p)
            p1 = (rng.random(m) * sizes[c1]).astype(np.int64)
            cis = rng.random(m) < 0.75
            # distance decay P(s) ~ s^-1 between 1 kb and the chrom length
            d = (1_000 * rng.random(m) ** (-1 / 0.9)).astype(np.int64)
            sign = np.where(rng.random(m) < 0.5, -1, 1)
            p2_cis = np.clip(p1 + sign * d, 0, sizes[c1] - 1)
            c2 = np.where(cis, c1, rng.choice(len(names), size=m, p=p))
            p2 = np.where(cis, p2_cis, (rng.random(m) * sizes[c2]).astype(np.int64))
            st = np.array(["+", "-"])
            df = pd.DataFrame({"id": [f"r{lo + i}" for i in range(m)],
                               "c1": names[c1], "p1": p1 + 1,
                               "c2": names[c2], "p2": p2 + 1,
                               "s1": st[rng.integers(0, 2, m)],
                               "s2": st[rng.integers(0, 2, m)]})
            df.to_csv(f, sep="\t", header=False, index=False)


def loops(n=50_000, out="loops.bedpe"):
    a = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None,
                    usecols=[0, 1, 2], names=["chrom", "start", "end"])
    rng = np.random.default_rng(5)
    rows = []
    centres = ((a["start"] + a["end"]) // 2).to_numpy()
    chroms = a["chrom"].to_numpy()
    for _ in range(n):
        i = int(rng.integers(0, len(a) - 1))
        j = min(len(a) - 1, i + 1 + int(rng.geometric(0.08)))
        if chroms[j] != chroms[i]:
            j = i - 1 - int(rng.geometric(0.08))
            if j < 0 or chroms[j] != chroms[i]:
                continue
        x, y = sorted((centres[i], centres[j]))
        jit = rng.integers(-1500, 1500, 2)
        x, y = max(0, x + jit[0]), max(0, y + jit[1])
        rows.append(f"{chroms[i]}\t{max(0, x - 2500)}\t{x + 2500}\t{chroms[i]}"
                    f"\t{max(0, y - 2500)}\t{y + 2500}\tL{len(rows)}\t{rng.integers(2, 60)}")
    (DATA / out).write_text("\n".join(rows) + "\n")


def loops_scale():
    """Loop sets of other sizes for the Architecture scaling benchmark."""
    for n in (12_500, 200_000, 1_000_000):
        loops(n, f"loops_{n}.bedpe")


def loops_trans(n=2_500):
    """loops.bedpe plus ``n`` inter-chromosomal loops between A peaks."""
    a = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None,
                    usecols=[0, 1, 2], names=["chrom", "start", "end"])
    rng = np.random.default_rng(7)
    centres = ((a["start"] + a["end"]) // 2).to_numpy()
    chroms = a["chrom"].to_numpy()
    rows = []
    while len(rows) < n:
        i, j = rng.integers(0, len(a), 2)
        if chroms[i] == chroms[j]:
            continue
        x, y = centres[i], centres[j]
        rows.append(f"{chroms[i]}\t{max(0, x - 2500)}\t{x + 2500}\t{chroms[j]}"
                    f"\t{max(0, y - 2500)}\t{y + 2500}\tT{len(rows)}\t{rng.integers(2, 20)}")
    cis = (DATA / "loops.bedpe").read_text()
    (DATA / "loops_trans.bedpe").write_text(cis + "\n".join(rows) + "\n")


def hic_trans():
    """hic.pairs plus ~40 contacts at each inter-chromosomal loop, binned at 5 kb,
    so the trans loops carry real Hi-C weight (the random trans background is
    far too sparse to hit any one pixel)."""
    import shutil
    import subprocess
    rng = np.random.default_rng(8)
    t = pd.read_csv(DATA / "loops_trans.bedpe", sep="\t", header=None, usecols=range(6))
    t = t[t[0] != t[3]]
    k = rng.poisson(40, len(t))
    i = np.repeat(np.arange(len(t)), k)
    x = ((t[1].to_numpy() + t[2].to_numpy()) // 2)[i] + rng.normal(0, 1500, len(i)).astype(np.int64)
    y = ((t[4].to_numpy() + t[5].to_numpy()) // 2)[i] + rng.normal(0, 1500, len(i)).astype(np.int64)
    extra = pd.DataFrame({"id": [f"t{n}" for n in range(len(i))], "c1": t[0].to_numpy()[i],
                          "p1": np.maximum(x, 1), "c2": t[3].to_numpy()[i], "p2": np.maximum(y, 1),
                          "s1": "+", "s2": "-"})
    out = DATA / "hic_trans.pairs"
    shutil.copyfile(DATA / "hic.pairs", out)
    extra.to_csv(out, sep="\t", header=False, index=False, mode="a")
    cooler = Path(sys.executable).parent / "cooler"
    subprocess.run([str(cooler), "cload", "pairs", "-c1", "2", "-p1", "3", "-c2", "4", "-p2", "5",
                    f"{DATA / 'hg38.chrom.sizes'}:5000", str(out), str(DATA / "hic_trans_5kb.cool")],
                   check=True)
    out.unlink()
    print(f"            hic_trans_5kb.cool: +{len(i):,} trans read pairs", flush=True)


def hichip_avp(n_short=3_000_000, n_bg=1_000_000):
    """HiC-Pro allValidPairs with a ChIP-like short-range component.

    H3K27ac weight per A peak: lognormal, x15 inside 400 SE-like clusters
    (CREs within 15 kb of a cluster centre). Short-range pairs sit on a CRE
    (pos1 ~ centre +- 250 bp) with a gap of 150 bp + Exp(450 bp), so most fall
    within 1 kb and some do not; the background is the first ``n_bg`` Hi-C pairs.
    """
    rng = np.random.default_rng(9)
    a = pd.read_csv(DATA / "peaks_A_100000.bed", sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"])
    cen = ((a.start + a.end) // 2).to_numpy()
    ch = a.chrom.to_numpy()
    w = rng.lognormal(0, 1, len(a))
    for c in rng.choice(len(a), 400, replace=False):
        near = (ch == ch[c]) & (np.abs(cen - cen[c]) < 15_000)
        w[near] *= 15
    i = rng.choice(len(a), n_short, p=w / w.sum())
    p1 = cen[i] + rng.normal(0, 250, n_short).astype(np.int64)
    gap = 150 + rng.exponential(450, n_short).astype(np.int64)
    p2 = p1 + gap
    s1 = np.where(rng.random(n_short) < 0.9, "+", "-")
    s2 = np.where(rng.random(n_short) < 0.9, "-", "+")
    short = pd.DataFrame({"id": [f"s{k}" for k in range(n_short)], "c1": ch[i], "p1": np.maximum(p1, 1), "s1": s1,
                          "c2": ch[i], "p2": np.maximum(p2, 2), "s2": s2, "size": gap})
    bg = pd.read_csv(DATA / "hic.pairs", sep="\t", comment="#", header=None, nrows=n_bg,
                     names=["id", "c1", "p1", "c2", "p2", "s1", "s2"])
    bg = bg[["id", "c1", "p1", "s1", "c2", "p2", "s2"]].assign(size=300)
    out = pd.concat([short, bg], ignore_index=True).sample(frac=1, random_state=1)
    out = out.assign(f1="HIC_x_1", f2="HIC_x_2", q1=42, q2=42)
    out.to_csv(DATA / "hichip.allValidPairs", sep="\t", header=False, index=False)
    print(f"            hichip.allValidPairs: {len(out):,} pairs", flush=True)


# ── genome FASTA + motifs ────────────────────────────────────────────────────

def genome():
    cs = read_chromsizes()
    rng = np.random.default_rng(6)
    alphabet = np.frombuffer(b"ACGT", dtype=np.uint8)
    with open(DATA / "genome.fa", "w") as f:
        for c in ("chr21", "chr22"):
            n = cs[c]
            seq = alphabet[rng.choice(4, size=n, p=[0.295, 0.205, 0.205, 0.295])]
            seq[:10_000] = ord("N")                      # telomeric N run
            s = seq.tobytes().decode()
            f.write(f">{c}\n")
            f.write("\n".join(s[i:i + 60] for i in range(0, n, 60)) + "\n")


def motifs():
    from pyjaspar import jaspardb
    ms = jaspardb().fetch_motifs(collection="CORE", tax_group=["vertebrates"])
    with open(DATA / "jaspar.txt", "w") as f:
        for m in ms:
            f.write(f">{m.matrix_id} {m.name}\n")
            for b in "ACGT":
                row = " ".join(str(int(round(x))) for x in m.counts[b])
                f.write(f"{b} [ {row} ]\n")
    print(f"            jaspar.txt: {len(ms)} motifs", flush=True)


STEPS = {"chromsizes": chromsizes, "peaks": peaks, "bw": bigwigs,
         "atlas": atlas, "gtf": gtf, "pairs": pairs, "loops": loops,
         "loops_trans": loops_trans, "loops_scale": loops_scale, "hic_trans": hic_trans,
         "hichip_avp": hichip_avp,
         "genome": genome, "motifs": motifs}

if __name__ == "__main__":
    todo = sys.argv[1:] or list(STEPS)
    for name in todo:
        t0 = _tic(name)
        STEPS[name]()
        _toc(t0)
