#!/usr/bin/env python3
"""Write cre_hub_view.html: a genomeblocks view laid out like a two-model
gene_tracks() figure (anchor O/E, node strength, ATAC, SE, prime hubs, copy
number + SV calls, CREs, loops, genes), on the synthetic benchmark data.

    python examples/columnar_prototype/make_view_demo.py [--fragment]

Model B is model A's Hi-C with multiplicative noise and an amplified region,
plus the trans contacts, so the two models differ visibly. All synthetic.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import genomeblocks.columnar as gbc
from genomeblocks.columnar.view import View

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "benchmarks" / "data"
CELLS = ["model A", "model B"]
COL = {"model A": "#46a8e4", "model B": "#ffa600"}
rng = np.random.default_rng(11)

cre = gbc.Loci.make(str(DATA / "peaks_A_100000.bed"))
genes = gbc.Genes.make(str(DATA / "genes.gtf"))
A = gbc.Architecture.make(cre, str(DATA / "loops_trans.bedpe"), verbose=False)
A.add_mcool(str(DATA / "hic_5kb.mcool"), resolution=5000, name="w_A", verbose=False)
A.add_mcool(str(DATA / "hic_trans_5kb.cool"), name="w_B", verbose=False)
A.ep.w_B = A.ep.w_B * np.exp(rng.normal(0, 0.6, A.n_links))
A.normalize(source="w_A", name="n_A", verbose=False).normalize(source="w_B", name="n_B", verbose=False)
A.annotate(genes, key="n_A", name="gene", verbose=False)
A.strength(key="n_A", name="strength_A", verbose=False).strength(key="n_B", name="strength_B", verbose=False)
hub = int(np.argmax(A.vp.strength_A))
hc, hpos = cre.chroms[hub], int(cre.centers[hub])
prime = {}
for s, k in zip(CELLS, ("n_A", "n_B")):
    uids = A.prime_hubs(key=k, gene="gene", verbose=False)["hub_uids"]
    prime[s] = np.isin(cre.uid, uids)

# super-enhancers, ROSE-style: stitch at 12.5 kb, rank by mean signal x width, slope-1 knee
def rose(bw):
    st = cre.slop(6250).merge().slop(-6250)
    sig = np.nan_to_num(st.signal([bw], span=True, n_bins=1, progress=False, verbose=False)[:, 0, 0])
    score = sig * st.lengths
    o = np.argsort(-score)
    y = score[o][::-1]
    yn = (y - y[0]) / max(y[-1] - y[0], 1e-9)
    slope = np.gradient(yn, np.linspace(0, 1, len(y)))
    cut = len(y) - int(np.argmax(slope >= 1))
    return st.take(np.sort(o[:cut]))
SE = {s: rose(str(DATA / f"signal_{i}.bw")) for i, s in enumerate(CELLS)}

# copy number (CNVkit-like 50 kb bins) and SV calls; model B carries an amplicon at the top hub
sizes = {l.split()[0]: int(l.split()[1]) for l in open(DATA / "hg38.chrom.sizes") if l.strip()}
def cnr(gain):
    rows = []
    for c, n in sizes.items():
        st = np.arange(0, n - 50_000, 50_000)
        seg = np.repeat(rng.normal(0, 0.25, len(st) // 60 + 1), 60)[:len(st)]
        l2 = seg + rng.normal(0, 0.18, len(st))
        if gain and c == hc:
            l2 += np.where(np.abs(st - hpos) < 1_500_000, 1.2, 0)
        rows.append(pd.DataFrame({"chromosome": c, "start": st, "end": st + 50_000, "log2": l2}))
    return pd.concat(rows, ignore_index=True)
CN = {"model A": cnr(False), "model B": cnr(True)}
SV = {"model A": pd.DataFrame({"chrom": [hc], "start": [hpos - 900_000], "end": [hpos - 300_000], "label": ["+-"]}),
      "model B": pd.DataFrame({"chrom": [hc, hc], "start": [hpos - 1_500_000, hpos + 200_000],
                               "end": [hpos + 1_500_000, hpos + 800_000], "label": ["-+", "++"]})}

# the gene to open on: the largest anchor budget (Σ O/E from its TSS anchor to its partners)
indptr, nbr, eid = A._adj()
budget = {}
for gname, rows in A.support(genes, r=5000, mode="center", rows=True).items():
    rs = set(rows.tolist())
    e = [eid[k] for v in rows for k in range(indptr[v], indptr[v + 1]) if nbr[k] not in rs]
    if e:
        budget[gname] = (A.ep.n_A[e].sum() + A.ep.n_B[e].sum(), e)
gene = max(budget, key=lambda g: budget[g][0])
e_best = max(budget[gene][1], key=lambda e: A.ep.n_A[e] + A.ep.n_B[e])
anchor_rows = set(A.support(genes, r=5000, mode="center", rows=True)[gene].tolist())
p_best = int(A.tgt[e_best] if A.src[e_best] in anchor_rows else A.src[e_best])
g_obj = genes[gene]
g_tss = g_obj.start if g_obj.strand == "+" else g_obj.end
t = A.trans
k = int(np.argmax(t.ep.n_B))
i, j = int(t.src[k]), int(t.tgt[k])
pair = " ".join(f"{cre.chroms[x]}:{int(cre.centers[x]) - 150_000}-{int(cre.centers[x]) + 150_000}" for x in (i, j))

v = View(A, genes=genes, samples=COL, title="CRE Hub View",
         subtitle="Synthetic benchmark data in a two-model layout. Search any gene, click a CRE to see its partners.",
         chrom_sizes=sizes, anchor_r=5000, anchor_mode="center", gene_half=1_000_000)
v.anchor_profile("graph O/E (anchor)")
v.cre_values("node strength", {"model A": A.vp.strength_A, "model B": A.vp.strength_B})
v.signal("ATAC", {s: str(DATA / f"signal_{i}.bw") for i, s in enumerate(CELLS)})
v.intervals("SE", SE)
v.intervals("prime hubs", prime)
v.points("copy number (log2)", CN, segments=SV)
v.cres("CREs (anchor yellow · partners black/grey)")
v.loops({"model A": "n_A", "model B": "n_B"})
v.genes()
v.region(gene, f"{g_obj.chrom}:{max(0, g_tss - 1_500_000)}-{g_tss + 1_500_000}", gene=gene,
         note=f"largest anchor budget · Σ O/E {budget[gene][0]:.0f} over {len(budget[gene][1])} contacts")
v.region("strongest hub", f"{hc}:{hpos - 2_000_000}-{hpos + 2_000_000}", anchor=hub,
         note=f"strength {A.vp.strength_A[hub]:.1f} (A) vs {A.vp.strength_B[hub]:.1f} (B) · amplified in B")
v.region("strongest trans loop", pair, anchor=i, note=f"O/E {t.ep.n_B[k]:.1f} in model B, side by side")
v.hubs(A.vp.strength_A, label="strength A")
v.mark("E-HUB", hc, hpos)
v.mark("E-1", g_obj.chrom, int(cre.centers[p_best]))
out = HERE / "cre_hub_view.html"
sizes_out = v.save(str(out), standalone="--fragment" not in sys.argv)
for k_, b in sizes_out.items():
    print(f"{k_:<28}{b / 1e6:6.2f} MB")
print("wrote", out)
