#!/usr/bin/env python3
"""Build and execute columnar_prototype.ipynb (outputs are stored in the file).

    python benchmarks/make_data.py peaks gtf pairs loops loops_trans hic_trans bw
    python examples/columnar_prototype/make_notebook.py
"""
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell

cells = [
md("""# Columnar genomeblocks: a notebook tour (prototype)

Branch `columnar-prototype` · package `genomeblocks.columnar`

**The idea in one picture.** Every table shares one `Genome`. The *row number* links the tables:

```
Genome  (chrom <-> code, shared by every table)
 ├─ Loci     chrom | start | end | strand          row i = CRE i, everywhere
 │    ├─ labels[i]    signal_cube[i]    A.vp[...][i]
 │    └─ Architecture  edges table: src | tgt | w | n | d   (src, tgt = Loci rows)
 │                     sorted:  | cis chr1 | cis chr2 | ... | trans |
 └─ Genes    genes ⇄ transcripts ⇄ exons   (linked by row numbers)
```

Data: the synthetic hg38-shaped benchmark set. It has 100k CREs, 20k genes and 52.5k loops (2.5k of them inter-chromosomal), plus 5 kb Hi-C. Make it with `benchmarks/make_data.py`."""),
code("""import time
import warnings
from contextlib import contextmanager
warnings.filterwarnings("ignore", message=".*(cairo|draw module).*")   # graph-tool without its drawing backend

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

@contextmanager
def timer(label):
    t = time.perf_counter()
    yield
    print(f"⏱  {label}: {1000 * (time.perf_counter() - t):,.0f} ms")

with timer("import genomeblocks.columnar"):
    import genomeblocks.columnar as gbc
    from genomeblocks.locus import Locus

DATA = "../../benchmarks/data\""""),
md("""## 1 · Load CREs: a table, not a list of objects

Four numpy columns. `make` sorts rows into genome order, so each chromosome becomes one block of rows."""),
code("""with timer("Loci.make, 100k peaks"):
    cre = gbc.Loci.make(f"{DATA}/peaks_A_100000.bed")
cre"""),
md("Looking at one row still gives you a normal `Locus`. It's a *view* onto the columns, so nothing is copied:"),
code("""v = cre[0]
print(v, "| isinstance(v, Locus):", isinstance(v, Locus), "| center:", v.center)
print("lookup by uid -> row", cre[v.uid].row)

chr21 = cre.by_chrom("chr21")
print(chr21, "| shares memory with cre:", np.shares_memory(chr21.starts, cre.starts))"""),
md("## 2 · Genes: three linked tables"),
code("""with timer("Genes.make, GTF with 877k lines"):
    genes = gbc.Genes.make(f"{DATA}/genes.gtf")
genes"""),
code("""g = genes["GENE0"]
print(g, "\\nTSS:", g.tss)
g.transcripts.to_pandas()"""),
md("""## 3 · Annotate every CRE at once

The annotation index (promoters, exons, UTRs, gene bodies) is built once and cached. After that, labelling is five vectorised overlap tests."""),
code("""from genomeblocks.columnar.genes import LABELS

with timer("annotation index (first use only)"):
    genes.annot
with timer("label 100k CREs"):
    labels = genes.labels(cre)          # one small int per CRE row
pd.Series(LABELS[labels]).value_counts().to_frame("CREs")"""),
md("""## 4 · Architecture: one graph stored as two tables

Vertices are the CRE rows. Edges are a table of `src | tgt | w | n | d`, sorted into per-chromosome **cis blocks** followed by one **trans block**. These loops include 2.5k inter-chromosomal ones."""),
code("""with timer("make → add_mcool → normalize → annotate → strength"):
    A = (gbc.Architecture.make(cre, f"{DATA}/loops_trans.bedpe")
           .add_mcool(f"{DATA}/hic_trans_5kb.cool")
           .normalize()
           .annotate(genes)
           .strength())
A"""),
md("""### The edge table drawn as a matrix

Rows and columns are CREs in genome order. Cis edges sit in the diagonal blocks (one block per chromosome). Trans edges are the red dots between blocks. It is one graph, and nothing was split apart."""),
code("""from matplotlib.patches import Rectangle

def draw_edges(ax, lo, hi, title):
    for c, (a, b) in cre.chrom_offsets.items():          # one shaded square per chromosome block
        if b > lo and a < hi:
            ax.add_patch(Rectangle((a, a), b - a, b - a, color="#1f6f8b", alpha=0.07, lw=0))
            ax.axvline(b, lw=0.3, c="0.75"); ax.axhline(b, lw=0.3, c="0.75")
    keep = (A.src >= lo) & (A.tgt < hi)
    cis = A.is_cis & keep
    trans = ~A.is_cis & keep
    for x, y in ((A.src, A.tgt), (A.tgt, A.src)):        # both triangles
        ax.scatter(x[cis], y[cis], s=0.6, c="#1f6f8b", rasterized=True)
        ax.scatter(x[trans], y[trans], s=0.6, c="#d1495b", rasterized=True)
    ax.set(xlim=(lo, hi), ylim=(hi, lo), xlabel="CRE row", ylabel="CRE row", title=title)

fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
draw_edges(ax[0], 0, len(cre), "all edges (shaded = chromosome blocks)")
z0, z1 = cre.chrom_offsets["chr1"][0], cre.chrom_offsets["chr3"][1]
draw_edges(ax[1], z0, z1, "zoom: chr1–chr3")
ax[1].scatter([], [], s=20, c="#1f6f8b", label=f"cis edges ({A.cis.n_links:,})")
ax[1].scatter([], [], s=20, c="#d1495b", label=f"trans edges ({A.trans.n_links:,})")
ax[1].legend(loc="upper right", frameon=True)

M = A.block_counts()
im = ax[2].imshow(np.log10(M.to_numpy() + 1), cmap="Blues")
ax[2].set_xticks(range(len(M)), M.columns, rotation=90, fontsize=7)
ax[2].set_yticks(range(len(M)), M.index, fontsize=7)
ax[2].set_title("edges per chromosome pair (log10)")
fig.colorbar(im, ax=ax[2], shrink=0.8)
plt.tight_layout()"""),
md("Cis edges hug the diagonal inside their chromosome's square, because most loops are short-range. Trans edges land in the off-diagonal squares. The heatmap counts the same thing per chromosome pair."),
md("## 5 · Views: cut the graph without copying it"),
code("""with timer("A.chrom('chr8')"):
    a8 = A.chrom("chr8")
with timer("A.cis and A.trans"):
    cis_part, trans_part = A.cis, A.trans
print(a8)
print(trans_part)
print("zero-copy:", np.shares_memory(a8.src, A.src), np.shares_memory(trans_part.ep.n, A.ep.n))"""),
code("""with timer("region chr8:20-40 Mb (first call also builds the interval index)"):
    r = A.region("chr8:20,000,000-40,000,000")
r"""),
md("## 6 · A CRE's partners include its trans partners"),
code("""t = A.trans
k = int(np.argmax(t.ep.w))           # the strongest trans edge
row = int(t.src[k])
A.neighbors(row)"""),
md("""## 7 · Per-chromosome work

Each cis block is independent, so you can loop over `A.chroms()` or hand the blocks to a process pool. Trans edges are one more small block at the end."""),
code("""with timer("per-chromosome summary"):
    per_chrom = pd.DataFrame([
        {"chrom": c, "edges": v.n_links, "loci": v.n_loci,
         "median distance (kb)": np.median(v.ep.d) / 1e3,
         "edges with Hi-C contacts": int((v.ep.w > 0).sum())}
        for c, v in A.chroms()])
per_chrom.head(8)"""),
md("## 8 · Graph algorithms see every edge, cis and trans"),
code("""with timer("graph-tool Graph from the tables (first call also imports graph-tool)"):
    g = A.graph()
with timer("connected components"):
    comp = A.components()
comp_cis = A.cis.components(name="component_cis_only")
n_comp = lambda c: len(np.unique(c[c >= 0]))
big = lambda c: np.bincount(c[c >= 0]).max()
print(f"components   with trans: {n_comp(comp):,}   cis only: {n_comp(comp_cis):,}")
print(f"largest one  with trans: {big(comp):,} CREs   cis only: {big(comp_cis):,} CREs")"""),
code("""import graph_tool.all as gt

with timer("pagerank"):
    pr = gt.pagerank(g)
A.vp.pagerank = np.asarray(pr.a)       # vertex i = Loci row i, so it just lines up
(A.to_frame()[["uid", "annot", "gene", "strength", "component", "pagerank"]]
   .sort_values("pagerank", ascending=False).head())"""),
md("""## 9 · Row alignment: signal, labels and graph columns line up by position

`signal()` here is the **classic** genomeblocks function, unchanged, running on the columnar Loci. Its cube rows are CRE rows, so selecting hubs is just an index."""),
code("""hubs = A.prime_hubs(verbose=False)
hub_rows = np.array([cre.row(u) for u in hubs["hub_uids"]])
rest = np.setdiff1d(np.flatnonzero(A.degree > 0), hub_rows)
pick = np.concatenate([hub_rows, np.random.default_rng(0).choice(rest, len(hub_rows), replace=False)])

with timer(f"signal() for {len(pick):,} CREs × 2 bigWigs"):
    cube = cre.take(pick).signal([f"{DATA}/signal_0.bw", f"{DATA}/signal_1.bw"],
                                 n_bins=60, progress=False, verbose=False)
print("cube:", cube.shape, "→ row j is CRE pick[j]")

x = np.linspace(-3, 3, 60)
fig, ax = plt.subplots(figsize=(6, 3.2))
ax.plot(x, np.nanmean(cube[:len(hub_rows), 0], 0), label=f"hub CREs ({len(hub_rows):,})")
ax.plot(x, np.nanmean(cube[len(hub_rows):, 0], 0), label="other linked CREs")
ax.set(xlabel="kb from CRE centre", ylabel="mean signal", title="signal_0 (synthetic)")
ax.legend(frameon=False); plt.tight_layout()"""),
md("""## 10 · A genome-wide test in one cell

Do edges that touch a promoter have higher O/E than other edges? Shuffle the promoter labels across linked CREs 1,000 times. Every permutation is a few array operations over all 120k edges.

*The data is synthetic, so expect no biological effect. The point is that the test is quick to run.*"""),
code("""is_prom = A.vp.annot == "Promoter-TSS"
oe = A.ep.n
touch = is_prom[A.src] | is_prom[A.tgt]
obs = oe[touch].mean() - oe[~touch].mean()

linked = np.flatnonzero(A.degree > 0)
k = int(is_prom[linked].sum())
rng = np.random.default_rng(1)
null = np.empty(1000)
with timer("1,000 permutations over 120k edges"):
    for i in range(1000):
        lab = np.zeros(len(cre), bool)
        lab[rng.choice(linked, k, replace=False)] = True
        t = lab[A.src] | lab[A.tgt]
        null[i] = oe[t].mean() - oe[~t].mean()
p = (1 + (null >= obs).sum()) / 1001
print(f"observed Δ O/E = {obs:.3f}   permutation p = {p:.3f}")"""),
md("## 11 · Save the session, reload in milliseconds"),
code("""import os, tempfile
SESSION = tempfile.mkdtemp(prefix="gb_session_")
with timer("save Architecture (+ its Loci) and Genes as parquet"):
    A.save(f"{SESSION}/architecture")
    genes.save(f"{SESSION}/genes")
size = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(SESSION) for f in fs)
print(f"on disk: {size / 1e6:.1f} MB")

with timer("load both back"):
    A2 = gbc.Architecture.load(f"{SESSION}/architecture")
    genes2 = gbc.Genes.load(f"{SESSION}/genes")
print(A2)
print("identical:", np.array_equal(A2.ep.n, A.ep.n), np.array_equal(A2.vp.gene, A.vp.gene))"""),
md("## 12 · Hand the tables to other tools"),
code("""A.edges_frame().head()"""),
code("""print(cre.to_arrow().schema)
cre.to_polars().head(3)"""),
code("""with timer("to_legacy(): the classic graph-tool Architecture (e.g. for its draw helpers)"):
    O = A.to_legacy()
O"""),
md("""## Cheat-sheet

| you want | write | cost |
|---|---|---|
| one CRE as an object | `cre[i]`, `cre["chr1:…"]` | µs |
| one chromosome of CREs | `cre.by_chrom("chr8")` | view, no copy |
| labels for all CREs | `genes.labels(cre)` | ~25 ms / 100k |
| graph for one chromosome | `A.chrom("chr8")` | view, no copy |
| cis only / trans only | `A.cis`, `A.trans` | view, no copy |
| a CRE's partners (cis + trans) | `A.neighbors(row)` | µs |
| graph algorithms | `A.graph()` → graph-tool | built once |
| per-CRE results | `A.vp[name][row]`, aligned with `cre` | — |
| save / reload | `A.save(path)` / `Architecture.load(path)` | parquet |
| classic API | `cre.signal(...)`, `A.to_legacy()` | unchanged |"""),
]

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
    "language_info": {"name": "python"}})
NotebookClient(nb, timeout=900, kernel_name="python3",
               resources={"metadata": {"path": str(HERE)}}).execute()
nbf.write(nb, HERE / "columnar_prototype.ipynb")
print("wrote", HERE / "columnar_prototype.ipynb")
