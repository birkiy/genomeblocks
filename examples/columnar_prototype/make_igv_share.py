#!/usr/bin/env python3
"""Write igv_share.html: a one-file genome browser to send to anyone.

    python examples/columnar_prototype/make_igv_share.py [path/to/igv.min.js] [--fragment]

With no argument, the page loads igv.js from jsDelivr (needs internet,
page ~5 MB). Pass a local igv.min.js to embed it, so the file also works
offline. --fragment writes page content only (for hosts that add their
own <html>/<head>).
"""
import sys
from pathlib import Path

import numpy as np

import genomeblocks.columnar as gbc
from genomeblocks.columnar.igv import igv_html

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "benchmarks" / "data"

cre = gbc.Loci.make(str(DATA / "peaks_A_100000.bed"))
genes = gbc.Genes.make(str(DATA / "genes.gtf"))
A = (gbc.Architecture.make(cre, str(DATA / "loops_trans.bedpe"), verbose=False)
       .add_mcool(str(DATA / "hic_trans_5kb.cool"), verbose=False)
       .normalize(verbose=False).annotate(genes, verbose=False).strength(verbose=False))

# regions: the five strongest hubs, plus the strongest trans loop shown side by side
s = A.vp.strength
top = np.argsort(-s)[:5]
regions, notes = [], {}
for r in top:
    c, a, b = cre.chroms[r], int(cre.starts[r]), int(cre.ends[r])
    key = f"{c}:{max(0, a - 150_000):,}-{b + 150_000:,}"
    regions.append(key)
    g = A.vp.gene[r]
    notes[key] = f"hub #{len(regions)} · strength {s[r]:.1f} · {A.vp.annot[r]}{' · ' + g if g else ''}"
t = A.trans
k = int(np.argmax(t.ep.w))
i, j = int(t.src[k]), int(t.tgt[k])
pair = " ".join(f"{cre.chroms[x]}:{max(0, int(cre.starts[x]) - 50_000):,}-{int(cre.ends[x]) + 50_000:,}"
                for x in (i, j))
regions.append(pair)
notes[pair] = f"strongest trans loop · {int(t.ep.w[k])} Hi-C reads · O/E {t.ep.n[k]:.1f}"

args = [a for a in sys.argv[1:] if not a.startswith("--")]
sizes = igv_html(str(HERE / "igv_share.html"), regions=regions, notes=notes,
                 loci={"CREs (ATAC peaks)": cre}, genes=genes,
                 signal={"signal 0": str(DATA / "signal_0.bw"), "signal 1": str(DATA / "signal_1.bw")},
                 architecture=A, title="CRE Hub Browser",
                 subtitle="Synthetic benchmark data: the five strongest CRE hubs and the strongest trans loop.",
                 chrom_sizes={l.split()[0]: int(l.split()[1]) for l in
                              open(DATA / "hg38.chrom.sizes") if l.strip()},
                 igv_js=args[0] if args else "cdn", standalone="--fragment" not in sys.argv)
for k, v in sizes.items():
    print(f"{k:<22}{v / 1e6:6.2f} MB")
