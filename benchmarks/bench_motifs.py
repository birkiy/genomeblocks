#!/usr/bin/env python3
"""Motif scanning: ``scan_motifs_matrix`` (lightmotif SIMD) vs alternatives.

Engines: genomeblocks (lightmotif), MEME FIMO, MOODS, Biopython, numpy.

Task: count forward-strand hits (log2-odds >= 13, pseudocount 0.1, uniform
background — genomeblocks' defaults) of JASPAR CORE vertebrate motifs in
500 bp windows around N loci. Every engine scores the identical PSSMs.

Parts:
  engines   N=1000 windows x 100 motifs, every engine + hit-count agreement
  library   N windows x all 1019 motifs for the fast engines
  workers   scan_motifs_matrix process-pool scaling (N=5000, 1019 motifs)
  stripe    cost of striping once vs. once per motif
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np

from common import DATA, Recorder, timeit

from genomeblocks import Loci, read_fasta
from genomeblocks.locus import Locus
from genomeblocks.motifs import scan_motifs_matrix

R, THR, PSEUDO = 250, 13.0, 0.1
JASPAR = str(DATA / "jaspar.txt")


def windows(genome, n, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        c = "chr21" if i % 2 else "chr22"
        pos = int(rng.integers(20_000, len(genome[c]) - 20_000))
        out.append(Locus(c, pos - 100, pos + 100))
    return Loci.from_records(out)


def motif_counts(limit=None):
    """(name, (W,4) count matrix in ACGT order) — the shared PSSM source."""
    import lightmotif
    out = []
    for m in lightmotif.load(JASPAR, format="jaspar16"):
        arr = np.asarray(m.counts, dtype=float)[:, [0, 1, 3, 2]]   # A C T G N -> ACGT
        out.append((m.name, arr))
        if limit and len(out) == limit:
            break
    return out


def logodds(counts):
    p = (counts + PSEUDO) / (counts + PSEUDO).sum(axis=1, keepdims=True)
    return np.log2(p / 0.25)


# ── engines: each returns {motif: total forward hits over all sequences} ─────

def eng_genomeblocks(L, genome, limit, path):
    df = scan_motifs_matrix(L, genome, path, format="jaspar16", r=R,
                            threshold=THR, norm=False, workers=1, verbose=False)
    return df.sum(axis=0).to_dict()


def eng_lightmotif_restripe(seqs, path):
    """Same scanner, but the sequence is re-striped for every motif."""
    import lightmotif
    out = {}
    for m in lightmotif.load(path, format="jaspar16"):
        pssm = m.counts.normalize(PSEUDO).log_odds()
        out[m.name] = sum(sum(1 for _ in lightmotif.scan(pssm, lightmotif.stripe(s), threshold=THR))
                          for s in seqs)
    return out


def eng_biopython(seqs, mats):
    from Bio.Seq import Seq
    from Bio.motifs.matrix import PositionSpecificScoringMatrix
    out = {}
    for name, cnt in mats:
        lo = logodds(cnt)
        pssm = PositionSpecificScoringMatrix("ACGT", {b: lo[:, i].tolist() for i, b in enumerate("ACGT")})
        out[name] = sum(sum(1 for _ in pssm.search(Seq(s), threshold=THR, both=False))
                        for s in seqs)
    return out


def eng_moods(seqs, mats):
    import MOODS.scan
    import MOODS.tools
    bg = [0.25] * 4
    lo = [logodds(c).T.tolist() for _, c in mats]          # (4, W) rows A C G T
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(lo, bg, [THR] * len(lo))
    tot = np.zeros(len(mats), np.int64)
    for s in seqs:
        for k, hits in enumerate(scanner.scan(s)):
            tot[k] += len(hits)
    return {name: int(tot[k]) for k, (name, _) in enumerate(mats)}


def eng_numpy(seqs, mats):
    code = np.full(256, 0, np.uint8)
    for i, b in enumerate(b"ACGT"):
        code[b] = i
    X = np.stack([code[np.frombuffer(s.encode(), np.uint8)] for s in seqs])  # (N, L)
    out = {}
    for name, cnt in mats:
        lo = logodds(cnt)                                   # (W, 4)
        W = lo.shape[0]
        win = np.lib.stride_tricks.sliding_window_view(X, W, axis=1)   # (N, L-W+1, W)
        score = lo[np.arange(W), win].sum(axis=-1)
        out[name] = int((score >= THR).sum())
    return out


FIMO = os.environ.get("FIMO", "fimo")


def write_meme(mats, path):
    """MEME-format motifs holding exactly genomeblocks' probabilities
    (pseudocount already applied), so FIMO scores the same PSSMs."""
    with open(path, "w") as f:
        f.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: +\n\n"
                "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for name, cnt in mats:
            p = (cnt + PSEUDO) / (cnt + PSEUDO).sum(axis=1, keepdims=True)
            f.write(f"MOTIF {name}\nletter-probability matrix: alength= 4 w= {len(p)} nsites= 20 E= 0\n")
            f.write("\n".join(" ".join(f"{x:.8f}" for x in row) for row in p) + "\n\n")


def write_fasta(seqs, path):
    with open(path, "w") as f:
        f.write("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))


def eng_fimo(fa, meme, mats):
    """MEME-suite FIMO, fastest mode (--text), forward strand, same PSSMs.
    p-value cut-off 1e-3 is loose enough to keep every score >= THR hit;
    hits are then counted at score >= THR like the other engines."""
    out = subprocess.run([FIMO, "--text", "--norc", "--thresh", "1e-3", "--bfile", "--uniform--",
                          "--motif-pseudo", "0", "--skip-matched-sequence", "--verbosity", "1",
                          meme, fa], capture_output=True, text=True, check=True).stdout
    tot = {name: 0 for name, _ in mats}
    for line in out.splitlines():
        if line.startswith(("motif_id", "#")) or not line:
            continue
        f = line.split("\t")
        if float(f[6]) >= THR - 1e-6:
            tot[f[0]] += 1
    return tot


def agreement(ref: dict, other: dict) -> dict:
    keys = sorted(ref)
    a = np.array([ref[k] for k in keys], float)
    b = np.array([other.get(k, 0) for k in keys], float)
    return {"total_hits": int(b.sum()), "ref_hits": int(a.sum()),
            "pearson_r": float(np.corrcoef(a, b)[0, 1]) if a.std() and b.std() else 1.0}


def subset_file(limit):
    """A JASPAR file holding the first ``limit`` motifs."""
    if limit is None:
        return JASPAR
    p = DATA / f"jaspar_{limit}.txt"
    blocks = open(JASPAR).read().split(">")[1:limit + 1]
    p.write_text("".join(">" + b for b in blocks))
    return str(p)


# ── parts ────────────────────────────────────────────────────────────────────

def part_engines(rec, genome):
    print("\n== engines (N=1000 x 100 motifs) ==")
    N, M = 1000, 100
    L = windows(genome, N)
    seqs = L.sequences(genome, r=R, upper=True)
    mats = motif_counts(M)
    path = subset_file(M)
    td = tempfile.mkdtemp()
    fa, meme = os.path.join(td, "seqs.fa"), os.path.join(td, "motifs.meme")
    write_fasta(seqs, fa); write_meme(mats, meme)
    res = {}
    engines = [
        ("genomeblocks scan_motifs_matrix (lightmotif)", lambda: eng_genomeblocks(L, genome, M, path), 3),
        ("MEME FIMO --text (CLI)", lambda: eng_fimo(fa, meme, mats), 3),
        ("lightmotif, re-striped per motif", lambda: eng_lightmotif_restripe(seqs, path), 3),
        ("MOODS (C++, all motifs per pass)", lambda: eng_moods(seqs, mats), 3),
        ("numpy sliding window", lambda: eng_numpy(seqs, mats), 1),
        ("Biopython PSSM.search", lambda: eng_biopython(seqs, mats), 1),
    ]
    ref = None
    for name, fn, rep in engines:
        t = timeit(lambda: res.__setitem__("o", fn()), repeat=rep, warmup=1 if rep > 1 else 0)
        ref = res["o"] if ref is None else ref
        bp = N * 2 * R * M
        rec.add(part="engines", engine=name, n_seqs=N, n_motifs=M, seconds=t["median"],
                runs=t["runs"], gbp_motif_per_s=bp / t["median"] / 1e9,
                **agreement(ref, res["o"]))


def part_library(rec, genome):
    print("\n== full JASPAR library (1019 motifs) ==")
    mats = motif_counts()
    for N in (1000, 5000):
        L = windows(genome, N, seed=1)
        seqs = L.sequences(genome, r=R, upper=True)
        bp = N * 2 * R * len(mats)
        t = timeit(lambda: eng_genomeblocks(L, genome, None, JASPAR), repeat=3)
        rec.add(part="library", engine="genomeblocks scan_motifs_matrix (lightmotif)",
                n_seqs=N, n_motifs=len(mats), seconds=t["median"], runs=t["runs"],
                gbp_motif_per_s=bp / t["median"] / 1e9)
        t = timeit(lambda: eng_moods(seqs, mats), repeat=3 if N == 1000 else 1)
        rec.add(part="library", engine="MOODS (C++, all motifs per pass)",
                n_seqs=N, n_motifs=len(mats), seconds=t["median"], runs=t["runs"],
                gbp_motif_per_s=bp / t["median"] / 1e9)
        if N == 1000:
            td = tempfile.mkdtemp()
            fa, meme = os.path.join(td, "seqs.fa"), os.path.join(td, "motifs.meme")
            write_fasta(seqs, fa); write_meme(mats, meme)
            t = timeit(lambda: eng_fimo(fa, meme, mats), repeat=1, warmup=0)
            rec.add(part="library", engine="MEME FIMO --text (CLI)", n_seqs=N,
                    n_motifs=len(mats), seconds=t["median"], runs=t["runs"],
                    gbp_motif_per_s=bp / t["median"] / 1e9)


def part_workers(rec, genome):
    print("\n== workers ==")
    N = 5000
    L = windows(genome, N, seed=2)
    for w in (1, 2, 4, 8):
        t = timeit(lambda: scan_motifs_matrix(L, genome, JASPAR, format="jaspar16",
                                              r=R, threshold=THR, workers=w, verbose=False),
                   repeat=3)
        rec.add(part="workers", workers=w, n_seqs=N, n_motifs=1019, seconds=t["median"],
                runs=t["runs"])


def part_stripe(rec, genome):
    import lightmotif
    print("\n== stripe vs scan cost ==")
    L = windows(genome, 2000, seed=3)
    seqs = L.sequences(genome, r=R, upper=True)
    m = next(iter(lightmotif.load(JASPAR, format="jaspar16")))
    pssm = m.counts.normalize(PSEUDO).log_odds()
    striped = [lightmotif.stripe(s) for s in seqs]
    t = timeit(lambda: [lightmotif.stripe(s) for s in seqs], repeat=5)
    rec.add(part="stripe", step="stripe 2000 seqs", seconds=t["median"], runs=t["runs"])
    t = timeit(lambda: [sum(1 for _ in lightmotif.scan(pssm, s, threshold=THR)) for s in striped],
               repeat=5)
    rec.add(part="stripe", step="scan 2000 striped seqs, 1 motif", seconds=t["median"], runs=t["runs"])


if __name__ == "__main__":
    parts = sys.argv[1:] or ["engines", "library", "workers", "stripe"]
    genome = read_fasta(str(DATA / "genome.fa"))     # {chrom: sequence}, the memory source
    rec = Recorder("motifs")
    for p in parts:
        globals()[f"part_{p}"](rec, genome)
    rec.save(r=R, threshold=THR, pseudocount=PSEUDO)
