#!/usr/bin/env python3
"""Build the design diagrams for the docs site.

    python docs/_tools/build_diagrams.py      # writes docs/_includes/diagrams/*.svg

Pages pull them in with ``{% include diagrams/<name>.svg %}`` inside a
``<figure class="gb-fig">``. Each function below draws one figure on explicit
coordinates (one claim per figure; the caption lives in the page).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svg import Diagram  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "_includes" / "diagrams"
FIGS = {}


def fig(fn):
    FIGS[fn.__name__.replace("_", "-")] = fn
    return fn


def row_label(d, y, title, sub=None, x=88):
    d.text(x, y, title.upper(), "cap", "end")
    if sub:
        d.text(x, y + 14, sub, "s", "end")


# ════════════════════════════════════════════════════════════════════════
# Overview
# ════════════════════════════════════════════════════════════════════════

@fig
def overview():
    d = Diagram("overview", 880, 652,
                "Package map: plotting modules read the domain modules; every domain module builds on "
                "Genome, Loci and Locus; whole-set work runs through one backend per family; and the "
                "boundary layer turns any table, file or region into Loci and back.")
    xs = [100 + i * 128.4 for i in range(6)]
    W = 118
    # ── plot row
    row_label(d, 40, "Plot", "matplotlib, lazy")
    d.node(xs[0], 20, W, 46, "signal_draw", "heatmaps · profiles", mono=True)
    d.node(xs[1], 20, W * 2 + 10.4, 46, "browserview · view · igv", "region figures · one-file pages", mono=True)
    d.node(xs[3], 20, W, 46, "motifs_draw", "motif heatmaps · logos", mono=True)
    d.node(xs[5], 20, W, 46, "architecture_draw", "layouts, any graph backend", mono=True)
    for x in (xs[0], xs[1], xs[2], xs[3], xs[5]):
        d.arrow([(x + W / 2, 68), (x + W / 2, 118)])
    d.text(xs[0] + W / 2 + 8, 97, "reads", "lbl", "start")
    # ── domain row
    row_label(d, 152, "Domain", "one table per data type")
    dom = [("signal", "bigWig → cube", "bigwig backend", "green"),
           ("Pairs · bedpe", "BEDPE table · Hi-C pairs", "tables backend", "purple"),
           ("Genes", "GTF → 3 linked tables", "tables · intervals", "navy"),
           ("motifs", "Library · hit matrices", "motifs · fasta backends", "navy"),
           ("Atlas", "bins × tracks, CSR", "scipy.sparse", "navy"),
           ("Architecture", "vertex + edge tables", "intervals · graph backends", "purple")]
    for x, (t, s, dep, k) in zip(xs, dom):
        d.rect(x, 120, W, 76, f"bx {k}")
        d.text(x + W / 2, 146, t, "t")
        d.text(x + W / 2, 162, s, "s")
        d.text(x + W / 2, 184, dep, "ms")
    # domain -> core bus
    for x in xs:
        d.line([(x + W / 2, 197), (x + W / 2, 226)])
    d.line([(xs[0] + W / 2, 226), (xs[5] + W / 2, 226)])
    d.arrow([(578, 226), (578, 262)])
    d.text(586, 248, "every module builds on Genome · Loci · Locus", "lbl", "start")
    # ── core row
    row_label(d, 300, "Core", "numpy columns")
    d.rect(100, 264, 156, 96, "bx")
    d.text(178, 290, "Genome", "t")
    d.text(178, 306, "names ↔ int codes · sizes", "s")
    d.chip(112, 318, 132, 24, "codes[i] → names", "white")
    d.text(178, 354, "one per table, re-coded on contact", "s")
    d.rect(272, 264, 436, 96, "bx green")
    d.text(288, 290, "Loci", "t", "start")
    d.text(324, 290, "intervals as columns · row number = the join key", "s", "start")
    for x, w, s in ((288, 84, "codes int32"), (378, 84, "starts int64"), (468, 78, "ends int64"),
                    (552, 86, "strands int8"), (644, 50, "cols {}")):
        d.chip(x, 300, w, 22, s, "white")
    d.text(288, 347, "lookup index built on first use, one per interval backend, dropped when rows change",
           "s", "start")
    d.rect(724, 264, 136, 96, "bx green")
    d.text(792, 290, "Locus", "t")
    d.text(792, 306, "one interval", "s")
    d.chip(734, 318, 116, 24, "L[i] → LocusView", "white")
    d.text(792, 354, "reads the columns", "s")
    # ── backends band
    row_label(d, 456, "Backends", "one engine per family")
    d.rect(130, 400, 700, 124, "bx page")
    d.text(146, 422, "gb.backends()", "tm", "start")
    d.text(258, 422, "default first, pip-installable · backend= on a call · with gb.use_backend(...) · "
                     "a missing engine raises with its install command", "s", "start")
    fams = [("intervals", "genomeblocks", "cgranges · ncls · bioframe", "pyranges · bedtools", "green"),
            ("bigwig", "pybigtools", "pybigwig · python", "", "green"),
            ("motifs", "lightmotif", "moods · biopython", "", "navy"),
            ("fasta", "genomeblocks", "pysam · pyfaidx", "memory · biopython", "navy"),
            ("tables", "polars, else pandas", "", "", ""),
            ("graph", "graph-tool, else scipy", "igraph · networkx", "", "purple")]
    for j, (fam, default, others, others2, k) in enumerate(fams):
        x = 146 + j * 113.6
        d.rect(x, 434, 104, 80, f"bx {k}" if k else "bx")
        d.text(x + 52, 452, fam, "t")
        d.text(x + 52, 470, default, "m")
        d.text(x + 52, 488, others, "ms")
        d.text(x + 52, 502, others2, "ms")
    d.arrow([(490, 362), (490, 398)])
    d.text(498, 384, "whole-set work: overlap, signal, scanning, graphs", "lbl", "start")
    # ── boundary band
    row_label(d, 584, "Boundary", "interop, in and out")
    d.rect(100, 548, 760, 84, "bx page")
    d.text(116, 570, "as_loci(x)", "tm", "start")
    d.text(200, 570, "every public function calls it on its inputs · columns found by name, lenient on spelling",
           "s", "start")
    for x, w, s in ((116, 56, "pandas"), (178, 52, "polars"), (236, 60, "pyarrow"), (302, 64, "bioframe"),
                    (372, 64, "pyranges"), (442, 74, "pybedtools"), (522, 60, "AnnData"), (588, 112, "BED · CSV · parquet"),
                    (706, 72, "region strs"), (784, 60, "lists · dicts")):
        d.chip(x, 580, w, 20, s, "white", "ms")
    d.text(116, 622, "out: to_pandas · to_polars · to_arrow · to_bioframe · to_pyranges · to_bedtool · to_anndata · "
                     "to_bed · save   ·   protocols: Arrow C stream · __dataframe__ · narwhals", "s", "start")
    d.arrow([(114, 546), (114, 362)], both=True)
    d.arrow([(846, 546), (846, 362)], both=True)
    return d


# ════════════════════════════════════════════════════════════════════════
# Loci
# ════════════════════════════════════════════════════════════════════════

@fig
def loci_setops():
    d = Diagram("loci-setops", 880, 232,
                "Set algebra keeps whole rows: A & B keeps every A interval that touches B, A − B keeps the "
                "rest, and merge fuses overlapping or book-ended intervals.")
    X0, X1 = 150, 860
    sx = lambda p: X0 + p * (X1 - X0) / 1000
    A = [(40, 130), (190, 300), (360, 430), (520, 640), (700, 760), (860, 950)]
    B = [(110, 160), (270, 330), (560, 590), (600, 680), (880, 900)]
    hit = lambda a: any(s < a[1] and a[0] < e for s, e in B)
    rows = [("A", "Loci", A, "fn"), ("B", "Loci", B, "fp")]
    y = 30
    for name, sub, ivs, k in rows:
        d.text(20, y + 4, name, "tm", "start")
        d.line([(X0, y), (X1, y)], "ln faint")
        for s, e in ivs:
            d.interval(sx(s), sx(e), y, 12, k)
        y += 36
    # results
    res = [("A & B", [a for a in A if hit(a)], "rows of A that touch B, kept whole"),
           ("A − B", [a for a in A if not hit(a)], "rows of A that touch nothing in B"),
           ("(A + B).merge()", None, "overlapping or book-ended rows fuse")]
    y += 10
    d.line([(20, y - 18), (X1, y - 18)], "grid")
    for name, ivs, note in res:
        d.text(20, y + 4, name, "tm", "start")
        d.line([(X0, y), (X1, y)], "ln faint")
        if ivs is None:
            allv = sorted(A + B)
            merged = [list(allv[0])]
            for s, e in allv[1:]:
                if s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            ivs = merged
            k = "fg"
        else:
            k = "fn"
        for s, e in ivs:
            d.interval(sx(s), sx(e), y, 12, k)
        d.text(X1, y + 20, note, "s", "end")
        y += 44
    # ghost outlines of the dropped A rows on the A & B line
    return d


@fig
def loci_index():
    d = Diagram("loci-index", 880, 372,
                "A single-window lookup on the numpy point index: rows are sorted by start with a running "
                "maximum of their ends; lo is the first row whose running max passes the query start, hi the "
                "first start at or past the query end, and only rows lo:hi are tested. Other interval "
                "backends answer the same window with their own index and are normalised to the same rows.")
    X0, X1 = 60, 840
    sx = lambda p: X0 + p * (X1 - X0) / 1000
    ivs = [(40, 130), (90, 380), (150, 210), (230, 300), (300, 420), (360, 470), (430, 520), (440, 460),
           (500, 560), (590, 690), (650, 760), (720, 800), (820, 900)]
    qs, qe = 470, 540
    srt = sorted(ivs)
    run = []
    m = -1
    for s, e in srt:
        m = max(m, e)
        run.append(m)
    lo = next(i for i, r in enumerate(run) if r > qs)                  # searchsorted(run, qs, 'right')
    hi = next((i for i, (s, e) in enumerate(srt) if s >= qe), len(srt))  # searchsorted(starts, qe, 'left')
    cand = {srt[i] for i in range(lo, hi)}
    # query band
    d.rect(sx(qs), 18, sx(qe) - sx(qs), 112, "bx hl", r=4)
    d.text(sx(qs) + (sx(qe) - sx(qs)) / 2, 34, "query", "t")
    d.text(sx(qs) + (sx(qe) - sx(qs)) / 2, 48, "[qs, qe)", "ms")
    # pile-up rows (explicit, so nothing collides with the labels)
    row_of = {(40, 130): 0, (150, 210): 0, (230, 300): 0, (360, 470): 0, (500, 560): 0, (650, 760): 0,
              (820, 900): 0, (300, 420): 1, (430, 520): 1, (590, 690): 1, (720, 800): 1, (90, 380): 2,
              (440, 460): 2}
    for s, e in ivs:
        y = 66 + row_of[(s, e)] * 22
        is_c = (s, e) in cand
        hitq = is_c and e > qs
        k = "fg" if hitq else ("fm" if is_c else "fl")
        d.interval(sx(s), sx(e), y, 12, k, title=f"[{s}, {e})")
    # legend (top right, clear of the pile-up)
    for x, k, lab in ((560, "fg", "overlaps the query"), (690, "fm", "tested, no overlap"),
                      (818, "fl", "skipped")):
        d.interval(x - 24, x - 6, 30, 10, k)
        d.text(x, 34, lab, "s", "start")
    # axis 1: rows sorted by start
    ay = 158
    d.line([(X0, ay), (X1, ay)], "ln faint")
    d.text(X0, ay - 10, "rows sorted by start", "s", "start")
    for s, e in srt:
        d.circle(sx(s), ay, 4.5, "dot " + ("navy" if (s, e) in cand else "muted"))
    # axis 2: running maximum of the ends
    by = 206
    d.line([(X0, by), (X1, by)], "ln faint")
    d.text(X0, by - 10, "running max of the ends, in that order", "s", "start")
    pts = [(sx(srt[i][0]), by) for i in range(len(srt))]
    for i, ((s, e), r) in enumerate(zip(srt, run)):
        d.circle(sx(s), by, 4.5, "dot " + ("navy" if i == lo else "muted"))
        if i == lo:
            d.text(sx(s), by + 18, f"run_max = {r} > qs", "ms")
    # pointers
    d.arrow([(sx(srt[lo][0]), by - 10), (sx(srt[lo][0]), ay + 10)], "navy")
    d.text(sx(srt[lo][0]) - 8, by - 20, "lo = searchsorted(run_max, qs, 'right')", "m", "end")
    d.arrow([(sx(qe), ay + 40), (sx(qe), ay + 10)], "navy")
    d.text(sx(qe) + 8, ay + 36, "hi = searchsorted(starts, qe, 'left')", "m", "start")
    # bracket
    a, b = sx(srt[lo][0]), sx(srt[hi][0]) if hi < len(srt) else X1
    d.line([(a, by + 28), (a, by + 36), (b, by + 36), (b, by + 28)], "ln navy")
    d.text((a + b) / 2, by + 52, "rows lo:hi are tested: keep end > qs, return their row numbers, sorted", "lbl")
    # backend seam
    yb = 290
    d.text(20, yb, "SAME WINDOW, OTHER ENGINES", "cap", "start")
    engines = [("genomeblocks", "PointIndex, numpy (default)", "green"),
               ("cgranges · ncls", "one tree per Loci, cached on it", ""),
               ("bioframe · pyranges · bedtools", "overlap_pairs of a 1-row Loci", "")]
    for j, (t, s, k) in enumerate(engines):
        x = 20 + j * 226
        d.node(x, yb + 10, 212, 44, t, s, kind=k, mono=True)
    d.arrow([(700, yb + 32), (730, yb + 32)])
    d.node(732, yb + 10, 128, 44, "same rows", "half-open · sorted", kind="green")
    return d


# ════════════════════════════════════════════════════════════════════════
# Genes
# ════════════════════════════════════════════════════════════════════════

@fig
def genes_model():
    d = Diagram("genes-model", 880, 318,
                "Genes.make reads GTF lines into three row-aligned tables: genes, transcripts (gene = the "
                "gene's row) and features (kind, transcript = the transcript's row); every table is a Loci "
                "with 0-based starts, so start − 1 is applied once while parsing.")
    d.rect(20, 20, 300, 278, "bx dark")
    d.text(36, 42, "GTF LINES · 1-BASED, CLOSED", "cap", "start")
    lines = [("gene", "chr8  127735434  +", "MYC"), ("transcript", "chr8  127735434  +", "T1"),
             ("exon", "chr8  127735434  +", "T1 #1"), ("five_prime_UTR", "chr8  127735434  +", "T1"),
             ("CDS", "chr8  127736231  +", "T1"), ("exon", "chr8  127740396  +", "T1 #2"),
             ("transcript", "chr8  127736623  +", "T2"), ("exon", "chr8  127736623  +", "T2 #1"), ("…", "", "")]
    for i, (f, c, k) in enumerate(lines):
        y = 68 + i * 20
        d.text(36, y, f, "m w", "start")
        d.text(146, y, c, "ms w2", "start")
        d.text(304, y, k, "m w", "end")
    d.text(36, 268, "polars (else pandas) parses the file;", "s w2", "start")
    d.text(36, 284, "joins on gene_id / transcript_id give the links", "s w2", "start")
    d.arrow([(322, 150), (356, 150)])
    d.text(339, 138, "make", "ms")
    # three tables
    def table(x, y, title, cols, rows, hl_col=None, k="navy"):
        w = 16 + sum(c[1] for c in cols)
        d.rect(x, y, w, 24 + 20 * len(rows) + 8, f"bx {k}")
        d.text(x + 10, y + 17, title, "t", "start")
        cx = x + 10
        for name, cw in cols:
            d.text(cx + 2, y + 34, name, "s", "start")
            cx += cw
        for i, r in enumerate(rows):
            yy = y + 52 + i * 20
            cx = x + 10
            for (name, cw), v in zip(cols, r):
                d.text(cx + 2, yy, str(v), "m", "start")
                cx += cw
        return w
    gx = 360
    base = (("row", 34), ("chrom", 46), ("start", 82), ("end", 82), ("±", 24))
    table(gx, 20, "genes", base + (("gene_id", 70), ("gene_name", 76), ("gene_type", 90)),
          [(0, "chr8", "127735433", "127742951", "+", "ENSG…", "MYC", "protein_coding")])
    d.text(gx + 10, 86, "0-based: 127735434 − 1 · TSS = start on '+', end − 1 on '−'", "s", "start")
    table(gx, 104, "transcripts", base + (("transcript_id", 96), ("gene →", 56)),
          [(0, "chr8", "127735433", "127742951", "+", "T1", 0), (1, "chr8", "127736622", "127742951", "+", "T2", 0)])
    table(gx, 196, "features", base + (("kind", 44), ("transcript →", 88), ("exon_number", 84)),
          [(0, "chr8", "127735433", "127736230", "+", "exon", 0, 1), (1, "chr8", "127735433", "127735582", "+", "5UTR", 0, 1),
           (2, "chr8", "127736230", "127736623", "+", "CDS", 0, 1), (3, "chr8", "127736622", "127737000", "+", "exon", 1, 1)])
    # links: row numbers, not ids
    tx_gene_x = gx + 10 + 34 + 46 + 82 + 82 + 24 + 96 + 12
    d.path(f"M{tx_gene_x + 14},{150} C{tx_gene_x + 60},{150} {tx_gene_x + 60},{66} {gx + 44},{66}", "ln navy")
    d.head(gx + 44, 66, math.pi, "hd navy", 6)
    d.path(f"M{tx_gene_x + 14},{170} C{tx_gene_x + 60},{170} {tx_gene_x + 60},{66} {gx + 44},{66}", "ln navy thin")
    ft_tx_x = gx + 10 + 34 + 46 + 82 + 82 + 24 + 44 + 12
    d.path(f"M{ft_tx_x + 12},{246} C{ft_tx_x + 50},{246} {ft_tx_x + 50},{156} {gx + 44},{156}", "ln navy thin")
    d.path(f"M{ft_tx_x + 12},{306} C{ft_tx_x + 56},{306} {ft_tx_x + 56},{176} {gx + 44},{176}", "ln navy thin")
    d.head(gx + 44, 156, math.pi, "hd navy", 6)
    d.head(gx + 44, 176, math.pi, "hd navy", 6)
    d.text(gx + 10, 316, "links are row numbers: transcripts['gene'][k] indexes the genes table, "
                         "features['transcript'][j] the transcripts table", "s", "start")
    return d


@fig
def genes_annot():
    d = Diagram("genes-annot", 880, 292,
                "Each CRE is checked against five merged interval sets in a fixed order; the first set it "
                "touches names it, and a CRE that touches none is Intergenic.")
    X0, X1 = 180, 860
    sx = lambda p: X0 + p * (X1 - X0) / 1000
    # gene model
    d.text(20, 36, "gene model", "s", "start")
    d.line([(sx(120), 32), (sx(820), 32)], "ln navy thin")
    for a, b, k, h in ((120, 215, "fnl", 8), (215, 230, "fn", 14), (330, 390, "fn", 14), (520, 580, "fn", 14),
                       (700, 735, "fn", 14), (735, 820, "fnl", 8)):
        d.interval(sx(a), sx(b), 32, h, k)
    d.path(f"M{sx(120)},26 L{sx(120)},14 L{sx(160)},14", "ln navy")
    d.head(sx(160) + 1, 14, 0, "hd navy", 6)
    d.text(sx(120), 8, "TSS", "ms", "middle")
    # CRE row
    cres = [("Promoter-TSS", 80, 100), ("5UTR", 196, 210), ("Exonic", 340, 360), ("Intronic", 440, 460),
            ("3UTR", 760, 780), ("Intergenic", 900, 920)]
    d.text(20, 66, "CREs", "s", "start")
    for _, a, b in cres:
        d.interval(sx(a), sx(b), 62, 10, "fg")
    tracks = [("prom", "TSS ± r", [(60, 180)]), ("utr5", "5′ UTRs", [(120, 215)]), ("utr3", "3′ UTRs", [(735, 820)]),
              ("exon", "exons", [(120, 230), (330, 390), (520, 580), (700, 820)]), ("body", "gene bodies", [(120, 820)])]
    ys = [96 + i * 30 for i in range(5)]
    for (name, sub, ivs), y in zip(tracks, ys):
        d.text(20, y + 4, f"{ys.index(y) + 1}", "cap", "start")
        d.text(36, y + 4, name, "m", "start")
        d.text(80, y + 4, sub, "s", "start")
        d.line([(X0, y), (X1, y)], "ln faint")
        for a, b in ivs:
            d.interval(sx(a), sx(b), y, 10, "fnl")
    # guides + first-hit dots
    for lab, a, b in cres:
        cx = sx((a + b) / 2)
        d.line([(cx, 70), (cx, ys[-1] + 14)], "ln faint dash")
        first = True
        for (name, _, ivs), y in zip(tracks, ys):
            if any(s0 < b and a < e0 for s0, e0 in ivs):
                d.circle(cx, y, 5.5 if first else 4, "dot green" if first else "dot muted")
                first = False
        w = max(52, 12 + 6.3 * len(lab))
        d.chip(cx - w / 2, ys[-1] + 26, w, 22, lab, "green" if lab != "Intergenic" else "", "s")
    d.text(X1, ys[-1] + 66, "filled dot = first set hit (wins) · grey = later sets, ignored", "s", "end")
    return d


@fig
def genes_isoforms():
    d = Diagram("genes-isoforms", 880, 262,
                "select_isoforms keeps isoforms whose TSS window overlaps a peak, picks the longest of those, "
                "and moves the gene's body and TSS onto it.")
    X0, X1 = 110, 690
    sx = lambda p: X0 + p * (X1 - X0) / 1000
    r = 35
    iso = [("T_long", 100, False, "✗  no peak in its TSS window"), ("T_mid", 380, True, "✓  peak · longest supported"),
           ("T_short", 560, True, "✓  peak")]
    ys = [80, 118, 156]
    # TSS windows as columns from the peak track down to the isoform
    for (name, tss, ok, _), y in zip(iso, ys):
        d.rect(sx(tss - r), 18, sx(tss + r) - sx(tss - r), y - 18 + 14, "bx hl" if ok else "bx ghost", r=4)
    d.text(20, 34, "ATAC peaks", "s", "start")
    d.line([(X0, 30), (X1, 30)], "ln faint")
    for a, b in ((358, 398), (545, 575), (760, 790)):
        d.interval(sx(a), sx(b), 30, 12, "fg")
    for (name, tss, ok, note), y in zip(iso, ys):
        d.text(20, y + 4, name, "m", "start")
        d.line([(sx(tss), y), (sx(950), y)], "ln navy thin")
        d.interval(sx(tss), sx(tss + 40), y, 12, "fn")
        d.interval(sx(880), sx(950), y, 12, "fn")
        d.text(712, y + 4, note, "s", "start")
    d.text(sx(100), 186, "TSS ± r", "s")
    # result
    y = 222
    d.text(20, y + 4, "Gene", "t", "start")
    d.rect(sx(100), y - 7, sx(950) - sx(100), 14, "bx ghost", r=3)
    d.interval(sx(380), sx(950), y, 14, "fn")
    d.text(712, y - 2, "canonical = T_mid: body and", "s", "start")
    d.text(712, y + 12, "TSS move to it; nothing dropped", "s", "start")
    d.text(sx(100) + 4, y + 26, "span before", "s", "start")
    return d


# ════════════════════════════════════════════════════════════════════════
# signal
# ════════════════════════════════════════════════════════════════════════

@fig
def signal_cube():
    d = Diagram("signal-cube", 880, 262,
                "signal() turns each region into a fixed window around its centre, asks the bigWig for "
                "n_bins summary values in one call, and writes them into one row of the cube.")
    # loci table
    d.text(20, 26, "LOCI", "cap", "start")
    cols = (("chrom", 20), ("start", 76), ("end", 146))
    for name, x in cols:
        d.text(x + 6, 46, name, "s", "start")
    rows = [("chr1", "10,400", "10,900"), ("chr1", "88,050", "88,700"), ("chr2", "3,100", "3,350"),
            ("chr2", "41,000", "41,600"), ("chrX", "8,800", "9,420")]
    for i, row in enumerate(rows):
        y = 56 + i * 26
        d.rect(20, y, 200, 24, "bx hl" if i == 1 else "bx", r=4)
        for (name, x), v in zip(cols, row):
            d.text(x + 6, y + 16, v, "m", "start")
    # signal panel
    X0, X1 = 300, 590
    d.text(X0, 26, "BIGWIG TRACK t", "cap", "start")
    base = 158
    pts = []
    import math
    for k in range(0, 101):
        x = X0 + k * (X1 - X0) / 100
        v = 18 + 70 * math.exp(-((k - 48) / 13) ** 2) + 22 * math.exp(-((k - 80) / 6) ** 2) + 6 * math.sin(k / 3)
        pts.append((x, base - v))
    d.path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + f" L{X1},{base} L{X0},{base} Z", "iv fl")
    d.line(pts, "ln thin")
    d.line([(X0, base), (X1, base)], "ln faint")
    nb = 10
    for b in range(nb):
        x0 = X0 + b * (X1 - X0) / nb
        x1 = X0 + (b + 1) * (X1 - X0) / nb
        seg = [y for x, y in pts if x0 <= x <= x1]
        mean = sum(base - y for y in seg) / len(seg)
        d.rect(x0 + 3, base - mean, x1 - x0 - 6, mean, "iv fn", r=2)
        if b:
            d.line([(x0, 50), (x0, base)], "grid")
    d.line([(X0, base + 12), (X0, base + 18), (X1, base + 18), (X1, base + 12)], "ln")
    d.text((X0 + X1) / 2, base + 34, "centre ± flank, split into n_bins", "s")
    d.text((X0 + X1) / 2, base + 52, "stats_array(chrom, L, R, n_bins): one call to the bigwig backend", "ms")
    d.arrow([(222, 81), (X0 - 10, 81)])
    # cube
    gx, gy, cw, ch = 652, 62, 18, 22
    for t, off in ((2, 24), (1, 12), (0, 0)):
        for i in range(5):
            for b in range(nb):
                x, y = gx + b * cw + off, gy + i * ch - off
                cls = "iv fn" if (t == 0 and i == 1) else ("iv fl" if t == 0 else "iv fl")
                d.rect(x, y, cw - 2, ch - 2, cls if t == 0 else "chip", r=2)
    d.text(gx - 8, gy + 2.5 * ch, "loci", "s", "end")
    d.text(gx + nb * cw / 2, gy + 5 * ch + 18, "bins", "s")
    d.text(gx + 24, gy - 32, "one grid per track", "s", "start")
    d.arrow([(X1 + 10, 81), (gx - 6, 88)])
    d.text(gx + nb * cw / 2, gy + 5 * ch + 40, "cube[i, t, :]", "m")
    d.text(gx + nb * cw / 2, gy + 5 * ch + 56, "(n_loci × n_tracks × n_bins)", "ms")
    return d


@fig
def signal_parallel():
    d = Diagram("signal-parallel", 880, 262,
                "With workers > 1 the cube lives in shared memory; each worker process opens its own "
                "bigWig handles and fills one slab in place, so nothing is copied back through pipes.")
    d.rect(20, 10, 330, 242, "bx page")
    d.text(36, 32, "main process", "t", "start")
    d.text(36, 48, "SharedMemory cube · copied out once, then unlinked", "s", "start")
    slabs = ["0–3", "4–7", "8–11", "12–15"]
    for k, lab in enumerate(slabs):
        y = 62 + k * 46
        on = k == 1
        d.rect(36, y, 298, 38, "bx green" if on else "bx", r=4)
        d.text(52, y + 24, f"tracks {lab}", "m", "start")
        for b in range(10):
            d.rect(160 + b * 16, y + 9, 13, 20, "iv fg" if on else "iv fl", r=2)
        d.rect(470, y, 390, 38, "bx green" if on else "bx")
        d.text(486, y + 17, f"worker {k + 1}", "t", "start")
        d.text(486, y + 31, f"own handles for bigWigs {lab} · all loci", "s", "start")
        d.arrow([(468, y + 19), (338, y + 19)], "green" if on else "faint")
    d.text(403, 54, "writes its slab in place", "lbl")
    return d


# ════════════════════════════════════════════════════════════════════════
# motifs
# ════════════════════════════════════════════════════════════════════════

@fig
def motifs_block():
    d = Diagram("motifs-block", 880, 290,
                "The windows are joined and striped once; each motif is scanned over the whole block in one "
                "call, and hits are mapped back to their window, dropping any that cross a boundary.")
    X0 = 20
    W = 120
    # windows
    d.text(X0, 24, "WINDOWS · centre ± r", "cap", "start")
    for i in range(5):
        x = X0 + i * (W + 10)
        d.rect(x, 34, W, 18, "iv fnl", r=3)
        d.text(x + W / 2, 47, f"w{i}", "ms")
    d.arrow([(320, 58), (320, 86)])
    d.text(328, 76, "join, stripe once", "lbl", "start")
    # block
    y = 128
    d.text(X0, 100, "ONE STRIPED BLOCK", "cap", "start")
    d.text(152, 100, "· hits of motif j from one engine call (lightmotif · MOODS · Biopython)", "s", "start")
    BW = 5 * W
    for i in range(5):
        d.rect(X0 + i * W, y, W, 22, "iv fnl" if i % 2 == 0 else "iv fl", r=0)
        d.text(X0 + i * W, y + 38, f"{2 * i}r" if i else "0", "ms", "middle")
    d.text(X0 + BW, y + 38, "10r", "ms", "middle")
    # windows are W=120 wide: w0 2 · w1 0 · w2 3 · w3 1 · w4 1 (the 476 hit straddles w3|w4)
    hits = [(46, "ok"), (88, "ok"), (262, "ok"), (300, "ok"), (330, "ok"), (410, "ok"), (476, "x"), (560, "ok")]
    for hx, kind in hits:
        x = X0 + hx
        if kind == "x":
            d.rect(x - 2, y - 14, 18, 8, "iv fw", r=2)
            d.text(x - 6, y - 20, "crosses w3 | w4: dropped", "s", "start")
        else:
            d.rect(x - 2, y - 14, 12, 8, "iv fg", r=2)
    # counts
    yc = 210
    d.arrow([(320, 176), (320, yc - 8)])
    d.text(328, 194, "window = searchsorted(offsets, pos) · bincount", "lbl", "start")
    for i, c in enumerate((2, 0, 3, 1, 1)):
        x = X0 + i * W
        d.rect(x + 30, yc, 60, 30, "bx green" if c else "bx", r=6)
        d.text(x + 60, yc + 20, str(c), "t")
    d.text(X0, yc + 52, "= column j of the hit matrix (÷ motif width if norm=True)", "s", "start")
    # matrix
    mx, my, cw, ch = 690, 60, 22, 26
    d.text(mx, 40, "LOCI × MOTIFS", "cap", "start")
    for i in range(5):
        for j in range(7):
            d.rect(mx + j * cw, my + i * ch, cw - 3, ch - 3, "iv fg" if j == 3 else "chip", r=3)
    d.text(mx + 3 * cw + 9, my + 5 * ch + 16, "j", "m")
    d.text(mx, my + 5 * ch + 48, "motifs are split across", "s", "start")
    d.text(mx, my + 5 * ch + 62, "worker processes; each", "s", "start")
    d.text(mx, my + 5 * ch + 76, "holds one copy of the block", "s", "start")
    d.arrow([(X0 + BW + 12, yc + 15), (mx + 3 * cw + 2, my + 5 * ch + 22)], "green")
    return d


# ════════════════════════════════════════════════════════════════════════
# BEDPE & pairs
# ════════════════════════════════════════════════════════════════════════

@fig
def bedpe_reader():
    d = Diagram("bedpe-reader", 880, 232,
                "Pairs.make is the one BEDPE parser: it reads the file through the tables backend into two "
                "row-aligned Loci, anchor a and anchor b, plus name / score columns; Architecture.make, the "
                "browser and Pairs.overlapping all consume that table.")
    d.rect(20, 30, 214, 172, "bx dark")
    d.text(34, 50, "LOOPS.BEDPE", "cap", "start")
    for i, l in enumerate(("chr8 127.73M  chr8 127.80M  l1 5", "chr8 127.73M  chr8 128.21M  l2 3",
                           "chr8 127.74M  chr2 41.05M   l3 1", "…")):
        d.text(34, 74 + i * 20, l, "m w", "start")
    d.text(34, 166, "any frame with these columns:", "s w2", "start")
    d.text(34, 182, "Pairs.from_frame · as_pairs", "s w2", "start")
    d.arrow([(236, 116), (292, 116)])
    d.text(264, 104, "Pairs.make", "ms")
    d.text(264, 132, "tables backend", "ms")
    d.rect(294, 30, 346, 172, "bx purple")
    d.text(310, 52, "Pairs", "t", "start")
    d.text(354, 52, "two Loci aligned by row + columns", "s", "start")
    # anchor a / anchor b / cols as three column groups over 3 rows
    hdr_y = 76
    groups = (("a", 310, (("codes", 46), ("starts", 54), ("ends", 54), ("±", 18)), "fp"),
              ("b", 494, (("codes", 46), ("starts", 54), ("ends", 54), ("±", 18)), "fp"),)
    for name, x, cols, k in groups:
        d.text(x, hdr_y - 10, f"P.{name}", "tm", "start")
        d.text(x + 26, hdr_y - 10, "Loci", "s", "start")
        cx = x
        for c, w in cols:
            d.text(cx, hdr_y + 6, c, "s", "start")
            cx += w
    vals_a = (("0", "127,73…", "127,73…", "+"), ("0", "127,73…", "127,73…", "+"), ("0", "127,74…", "127,74…", "−"))
    vals_b = (("0", "127,80…", "127,80…", "−"), ("0", "128,21…", "128,21…", "+"), ("1", "41,05…", "41,05…", "+"))
    for i in range(3):
        y = hdr_y + 28 + i * 20
        d.rect(304, y - 13, 326, 18, "bx hl" if i == 2 else "chip white", r=3)
        for (name, x, cols, k), vals in zip(groups, (vals_a, vals_b)):
            cx = x
            for (c, w), v in zip(cols, vals[i]):
                d.text(cx, y, v, "ms", "start")
                cx += w
    d.text(310, 164, "P.cols: name · score · …", "m", "start")
    d.text(310, 184, "is_cis = a.codes == b.codes · distance = |mid2 − mid1|, ∞ across chromosomes", "s", "start")
    for i, (t, s) in enumerate((("Architecture.make", "anchors (mid ± r) → overlap_pairs → edges"),
                                ("browser · View · igv_html", "arcs between anchor midpoints"),
                                ("Pairs.overlapping(loci)", "pairtobed: an anchor touches the set"))):
        y = 22 + i * 62
        d.node(690, y, 170, 50, t, s, mono=True)
        d.arrow([(642, 116), (666, 116), (666, y + 25), (688, y + 25)])
    d.text(690, 212, "every consumer calls as_pairs()", "s", "start")
    return d


@fig
def bedpe_count():
    d = Diagram("bedpe-count", 880, 300,
                "count_pairs places every window on one sorted key (chromosome code, start); each read end is "
                "located with one searchsorted, then counts are a bincount.")
    d.text(20, 24, "PAIRS CHUNK (pandas)", "cap", "start")
    hdr = ("chrom1", "pos1", "chrom2", "pos2")
    xs = (20, 82, 150, 212)
    for x, h in zip(xs, hdr):
        d.text(x + 4, 44, h, "s", "start")
    rows = (("chr1", "12,040", "chr1", "61,500"), ("chr1", "60,380", "chr2", "5,120"),
            ("chr2", "4,990", "chr2", "18,200"), ("chr1", "33,700", "chrX", "9,410"))
    for i, r in enumerate(rows):
        y = 52 + i * 24
        d.rect(20, y, 260, 22, "bx hl" if i == 1 else "bx", r=4)
        for x, v in zip(xs, r):
            d.text(x + 4, y + 15, v, "m", "start")
    d.text(20, 172, "chrom columns are categorical:", "s", "start")
    d.text(20, 188, "names → codes once per category", "s", "start")
    for k, (c, n) in enumerate((("chr1", 0), ("chr2", 1), ("chrX", 2))):
        d.chip(20 + k * 84, 198, 76, 22, f"{c} → {n}", "white")
    # key axis
    X0, X1 = 330, 640
    d.text(X0, 24, "WINDOW KEYS · code · 2⁴⁰ + start (sorted)", "cap", "start")
    ay = 92
    wins = [(0, 0), (0, 1), (0, 2), (0, 4), (1, 0), (1, 1), (2, 0)]
    pos = {}
    x = X0
    prev = 0
    for c, w in wins:
        if c != prev:
            x += 14
        if (c, w) == (0, 4):
            d.text(x + 12, ay + 4, "·", "s")
            x += 24
        prev = c
        d.rect(x, ay - 14, 34, 28, "bx navy" if (c, w) == (0, 1) else "bx", r=3)
        d.text(x + 17, ay + 4, f"w{w}", "ms")
        pos[(c, w)] = x
        x += 37
    for c, lab in ((0, "chr1"), (1, "chr2"), (2, "chrX")):
        xs_ = [pos[k] for k in pos if k[0] == c]
        d.text((min(xs_) + max(xs_) + 34) / 2, ay - 22, lab, "ms")
    d.text(pos[(0, 4)] - 12, ay + 32, "gap", "s")
    d.arrow([(pos[(0, 1)] + 18, ay + 64), (pos[(0, 1)] + 18, ay + 16)], "navy")
    d.text(pos[(0, 1)] + 18, ay + 78, "60,380 on chr1 → window 1", "s")
    d.text(X0, ay + 104, "i = searchsorted(keys, code·2⁴⁰ + pos) − 1", "ms", "start")
    d.text(X0, ay + 120, "keep i if same code and pos < end, else −1", "ms", "start")
    # outputs
    ox = 680
    d.text(ox, 24, "COUNTS", "cap", "start")
    d.text(ox, 44, "count_pairs: windows × partner", "s", "start")
    for i in range(4):
        for j in range(3):
            d.rect(ox + j * 30, 52 + i * 18, 27, 15, "iv fg" if (i, j) == (1, 1) else "chip", r=2)
    for j, c in enumerate(("chr1", "chr2", "chrX")):
        d.text(ox + j * 30 + 13, 136, c, "ms")
    d.text(ox, 172, "count_pairs_2d: windows × windows", "s", "start")
    for i in range(5):
        for j in range(5):
            on = (i, j) in {(0, 0), (1, 1), (0, 2), (2, 0), (3, 3), (1, 4), (4, 1), (2, 2), (4, 4)}
            d.rect(ox + j * 18, 182 + i * 18, 15, 15, "iv fp" if on else "chip", r=2)
    d.text(ox + 100, 230, "sparse CSR", "ms", "start")
    d.text(ox, 290, "bincount / unique per chunk", "s", "start")
    return d


# ════════════════════════════════════════════════════════════════════════
# Atlas
# ════════════════════════════════════════════════════════════════════════

@fig
def atlas_index():
    d = Diagram("atlas-index", 880, 318,
                "Atlas stores one bit per (genome bin, track); a query is the column sum over its bins, "
                "which fills the 2×2 table of a Fisher test for every track at once.")
    mx, my, cw, ch = 110, 44, 30, 18
    nt, nb = 9, 12
    rng = [((i * 7 + j * 13) % 11) < 3 for i in range(nb) for j in range(nt)]
    q = {3, 4, 9}
    d.text(mx, 24, "M: GENOME BINS × TRACKS (CSR, 1 bit each)", "cap", "start")
    for i in range(nb):
        if i in q:
            d.rect(mx - 6, my + i * ch - 1, nt * cw + 8, ch, "bx hl", r=3)
        for j in range(nt):
            on = rng[i * nt + j]
            d.rect(mx + j * cw + 4, my + i * ch + 2, cw - 8, ch - 6, ("iv fn" if i in q else "iv fm") if on else "chip", r=2)
    for i, lab in ((0, "chr1"), (7, "chr2")):
        d.text(mx - 12, my + i * ch + 12, lab, "ms", "end")
    d.line([(mx - 8, my + 7 * ch - 1), (mx + nt * cw, my + 7 * ch - 1)], "ln faint")
    d.text(mx - 12, my + 4 * ch + 12, "query", "s", "end")
    d.text(mx - 12, my + 4 * ch + 26, "bins", "s", "end")
    for j in range(nt):
        d.text(mx + j * cw + cw / 2, my + nb * ch + 14, f"t{j + 1}", "ms")
    sums = [sum(rng[i * nt + j] for i in q) for j in range(nt)]
    yb = my + nb * ch + 30
    d.text(mx - 12, yb + 14, "a =", "m", "end")
    for j, a in enumerate(sums):
        d.rect(mx + j * cw + 3, yb, cw - 6, 22, "bx green" if a else "bx", r=4)
        d.text(mx + j * cw + cw / 2, yb + 15, str(a), "t")
    d.text(mx, yb + 44, "M[query_bins].sum(axis=0): one sparse product for all tracks", "ms", "start")
    # 2x2
    tx, ty = 520, 70
    d.text(tx, 50, "PER TRACK, VECTORISED", "cap", "start")
    cells = (("a", "query ∩ track"), ("b", "query only"), ("c", "track only"), ("d", "neither"))
    for k, (v, lab) in enumerate(cells):
        x = tx + (k % 2) * 120
        y = ty + (k // 2) * 56
        d.rect(x, y, 112, 48, "bx green" if v == "a" else "bx", r=6)
        d.text(x + 56, y + 21, v, "tm")
        d.text(x + 56, y + 38, lab, "s")
    d.text(tx + 248, ty + 28, "c = track bins − a", "ms", "start")
    d.text(tx + 248, ty + 46, "b = query bins − a", "ms", "start")
    d.text(tx + 248, ty + 84, "d = all bins − …", "ms", "start")
    d.arrow([(tx + 116, ty + 118), (tx + 116, ty + 150)])
    d.rect(tx, ty + 154, 340, 44, "bx navy")
    d.text(tx + 170, ty + 174, "Fisher's exact test", "t")
    d.text(tx + 170, ty + 189, "log2 odds · p · GIGGLE-style score per track", "s")
    d.text(tx, ty + 222, "ref = a second set replaces the genome in b, c, d", "s", "start")
    d.text(tx, ty + 238, "bootstrap() reshuffles the query per chromosome", "s", "start")
    return d


# ════════════════════════════════════════════════════════════════════════
# Architecture
# ════════════════════════════════════════════════════════════════════════

@fig
def arch_pipeline():
    d = Diagram("arch-pipeline", 880, 198,
                "The Architecture pipeline: six chained calls, each reading the previous step's output "
                "and writing an edge column (ep) or a vertex column (vp).")
    steps = [("make", "Loci + BEDPE / Pairs", "vertices · edge table"),
             ("add_mcool", ".mcool at 5 kb", "ep.w (Hi-C count)"),
             ("normalize", "CRE positions", "ep.d · ep.n (O/E)"),
             ("annotate", "Genes", "vp.annot · vp.gene"),
             ("strength", "ep.n", "vp.strength"),
             ("prime_hubs", "vp.strength", "hub CREs → genes")]
    w, g, x0 = 116, 12.4, 110
    row_label(d, 52, "reads")
    row_label(d, 104, "step")
    row_label(d, 160, "writes")
    for i, (t, r, wr) in enumerate(steps):
        x = x0 + i * (w + g)
        d.text(x + w / 2, 56, r, "s")
        d.rect(x, 76, w, 50, "bx purple")
        d.text(x + w / 2, 106, "." + t + "()" if i else "make()", "tm")
        d.chip(x + 4, 148, w - 8, 24, wr, "", "s")
        d.line([(x + w / 2, 64), (x + w / 2, 74)], "ln faint")
        d.line([(x + w / 2, 128), (x + w / 2, 146)], "ln faint")
        if i:
            d.arrow([(x - g + 1, 101), (x - 1, 101)], "purple")
    return d


@fig
def arch_weights():
    d = Diagram("arch-weights", 880, 300,
                "make turns each loop into edges between the CREs under its two anchors; add_mcool gives every "
                "edge its Hi-C pixel count, split between edges sharing the pixel; normalize divides by the "
                "distance-decay expectation (cis) or the mean trans weight (trans).")
    # panel 1: make
    d.text(20, 22, "1 · MAKE", "cap", "start")
    gy = 112
    d.line([(20, gy), (280, gy)], "ln faint")
    for a, b in ((40, 66), (84, 104), (222, 250)):
        d.interval(a, b, gy, 12, "fg")
    d.rect(30, gy - 22, 88, 44, "bx ghost", r=6)
    d.rect(206, gy - 22, 60, 44, "bx ghost", r=6)
    d.text(74, gy + 38, "anchor 1 ± r", "s")
    d.text(236, gy + 38, "anchor 2 ± r", "s")
    d.path(f"M74,{gy - 24} Q155,{gy - 90} 236,{gy - 24}", "ln purple")
    d.text(155, gy - 64, "loop", "lbl")
    for (x, lab) in ((53, "a"), (94, "b"), (236, "c")):
        d.text(x, gy - 10, lab, "m")
    V = {"a": (70, 206), "b": (110, 246), "c": (230, 222)}
    d.line([V["a"], V["c"]], "ln purple thick")
    d.line([V["b"], V["c"]], "ln purple thick")
    for k, (x, y) in V.items():
        d.circle(x, y, 9, "dot green")
        d.text(x - 16, y + 4, k, "m", "end") if k != "c" else d.text(x + 16, y + 4, k, "m", "start")
    d.text(150, 280, "edges a–c, b–c (never within one anchor)", "s")
    # panel 2: add_mcool
    px, py, c = 330, 44, 26
    d.text(px, 22, "2 · ADD_MCOOL", "cap", "start")
    for i in range(6):
        for j in range(i, 6):
            on = (i, j) == (1, 4)
            d.rect(px + j * c, py + i * c, c - 3, c - 3, "iv fp" if on else ("iv fpl" if j - i < 2 else "chip"), r=2)
    d.text(px + 4 * c + 11, py + c + 17, "12", "s w")
    d.text(px + 1 * c + 11, py + 6 * c + 14, "bin i", "ms")
    d.text(px + 4 * c + 11, py + 6 * c + 14, "bin j", "ms")
    d.text(px, py + 6 * c + 40, "a, b → bin i · c → bin j", "ms", "start")
    d.text(px, py + 6 * c + 58, "edges a–c, b–c share pixel (i, j)", "s", "start")
    d.text(px, py + 6 * c + 76, "ep.w = 12 / 2 = 6 each", "m", "start")
    # panel 3: normalize
    X0, Y0, X1, Y1 = 610, 52, 860, 210
    d.text(X0, 22, "3 · NORMALIZE", "cap", "start")
    d.line([(X0, Y0), (X0, Y1), (X1, Y1)], "ln faint")
    d.text(X0 + 4, Y0 + 4, "w (log)", "s", "start")
    d.text(X1, Y1 + 16, "distance d (log)", "s", "end")
    import math
    for k in range(34):
        t = k / 33
        x = X0 + 10 + t * 220
        yfit = Y0 + 20 + t * 120
        y = yfit + 22 * math.sin(k * 2.3) * (0.6 + 0.4 * ((k * 7) % 5) / 5)
        d.circle(x, y, 2.6, "dot muted")
    d.line([(X0 + 6, Y0 + 16), (X0 + 236, Y0 + 144)], "ln navy thick")
    d.text(X0 + 128, Y0 + 66, "fit: E(d) = C / d^α", "m", "start")
    d.text(X0, Y1 + 40, "cis: ep.n = w / E(d)", "m", "start")
    d.text(X0, Y1 + 58, "trans: d = ∞, E = mean trans w", "m", "start")
    return d


@fig
def arch_hubs():
    d = Diagram("arch-hubs", 880, 290,
                "annotate gives promoter CREs their nearest gene and every other CRE the gene of its "
                "highest-O/E promoter neighbour; prime_hubs cuts the ranked node strengths at the "
                "slope-1 knee and returns the hubs' genes.")
    d.text(20, 22, "ANNOTATE · GENE ASSIGNMENT", "cap", "start")
    P = {"MYC": (110, 90), "PVT1": (330, 92)}
    E = {"e1": (220, 60), "e2": (210, 196), "e3": (60, 210), "e4": (380, 214)}
    edges = [("e1", "MYC", 3.1), ("e1", "PVT1", 1.2), ("e2", "MYC", 0.9), ("e2", "PVT1", 2.4),
             ("e3", "MYC", 1.7), ("e4", "PVT1", 1.1)]
    for a, b, w in edges:
        (x1, y1), (x2, y2) = E[a], P[b]
        top = (a == "e1" and b == "MYC") or (a == "e2" and b == "PVT1") or a in ("e3", "e4")
        d.line([(x1, y1), (x2, y2)], "ln purple thick" if top else "ln faint")
        mx_, my_ = (x1 + x2) / 2, (y1 + y2) / 2
        d.rect(mx_ - 13, my_ - 9, 26, 16, "chip white", r=8)
        d.text(mx_, my_ + 3, f"{w}", "ms")
    for g, (x, y) in P.items():
        d.rect(x - 34, y - 14, 68, 28, "bx navy", r=6)
        d.text(x, y + 4, g, "t")
    assigned = {"e1": "MYC", "e2": "PVT1", "e3": "MYC", "e4": "PVT1"}
    for e, (x, y) in E.items():
        d.circle(x, y, 10, "dot green")
        d.text(x, y + 26, f"{e} → {assigned[e]}", "ms")
    d.text(230, 268, "each enhancer takes the gene of its highest-O/E promoter neighbour", "s")
    # strength curve
    X0, Y0, X1, Y1 = 500, 40, 860, 230
    d.text(X0, 22, "PRIME_HUBS · SLOPE-1 KNEE", "cap", "start")
    d.line([(X0, Y0), (X0, Y1), (X1, Y1)], "ln faint")
    import math
    pts = []
    for k in range(60):
        t = k / 59
        y = Y1 - 6 - (Y1 - Y0 - 20) * math.exp(-6.5 * t)
        pts.append((X0 + 6 + t * (X1 - X0 - 12), y))
    cut = 11
    d.rect(X0 + 2, Y0, pts[cut][0] - X0, Y1 - Y0 - 2, "bx hl", r=4)
    d.line(pts, "ln navy thick")
    d.circle(*pts[cut], 5, "dot navy")
    d.text(pts[cut][0] + 10, pts[cut][1] - 4, "knee: tangent slope = 1", "s", "start")
    d.text(pts[cut][0] + 10, pts[cut][1] + 10, "on the [0, 1]-scaled curve", "s", "start")
    d.text(X0 + 10, Y0 + 16, "hubs", "t", "start")
    d.text(X1, Y1 + 16, "CREs ranked by node strength (Σ incident O/E)", "s", "end")
    d.text(X0, Y1 + 44, "hub CREs → their vp.gene → prime genes", "m", "start")
    return d


# ════════════════════════════════════════════════════════════════════════
# browser
# ════════════════════════════════════════════════════════════════════════

@fig
def browser_tracks():
    d = Diagram("browser-tracks", 880, 336,
                "browser() gives every track its own axis on a shared x range and picks the drawer from "
                "the track's type: intervals through as_loci, loops through as_pairs, bigWigs as binned "
                "summaries of the region only, through the bigwig backend.")
    X0, X1 = 150, 600
    rows = [("ruler", None, "region → chrom, start, end"),
            ("ATAC", ".bw", "binned means, one native call"),
            ("H3K27ac ×2", "[.bw, .bw]", "replicates averaged per bin"),
            ("reads", ".bam", "per-base pileup, mismatches coloured"),
            ("peaks", "Loci · frame · .bed", "rectangles, via as_loci"),
            ("loops", "Pairs · .bedpe", "half-sine arcs between anchors"),
            ("genes", "Genes", "stacked models from the 3 tables")]
    import math
    for i, (name, src, note) in enumerate(rows):
        y = 26 + i * 44
        d.text(20, y + 14, name, "t", "start")
        if src:
            d.text(20, y + 29, src, "ms", "start")
        d.text(632, y + 14, note, "s", "start")
        base = y + 30
        if i == 0:
            d.line([(X0, y + 18), (X1, y + 18)], "ln")
            for k in range(6):
                x = X0 + k * (X1 - X0) / 5
                d.line([(x, y + 12), (x, y + 18)], "ln")
                d.text(x, y + 8, f"{127.70 + k * 0.03:.2f} Mb", "ms", "end" if k == 5 else ("start" if k == 0 else "middle"))
        elif i in (1, 2):
            pts = [(X0 + k * (X1 - X0) / 60, base - (4 + 22 * math.exp(-((k - (24 if i == 1 else 36)) / 4) ** 2)
                                                      + 10 * math.exp(-((k - 48) / 3) ** 2)))
                   for k in range(61)]
            d.path("M" + " L".join(f"{x:.1f},{yy:.1f}" for x, yy in pts) + f" L{X1},{base} L{X0},{base} Z",
                   "iv fg" if i == 1 else "iv fn")
        elif i == 3:
            for k in range(90):
                x = X0 + k * (X1 - X0) / 90
                h = 6 + 16 * math.exp(-((k - 50) / 14) ** 2)
                cls = "iv fw" if k in (44, 61) else "iv fm"
                d.rect(x, base - h, (X1 - X0) / 90 - 0.6, h, cls, r=0)
        elif i == 4:
            for a, b in ((60, 110), (180, 210), (260, 330), (420, 440)):
                d.rect(X0 + a, y + 12, b - a, 12, "iv fg", r=2)
        elif i == 5:
            for a, b in ((85, 290), (195, 430), (300, 340)):
                d.path(f"M{X0 + a},{base} Q{X0 + (a + b) / 2},{base - (b - a) * 0.18 - 8} {X0 + b},{base}", "ln purple")
        else:
            d.line([(X0 + 40, y + 12), (X0 + 330, y + 12)], "ln navy thin")
            for a, b, h in ((40, 60, 6), (60, 80, 12), (150, 175, 12), (300, 330, 6)):
                d.rect(X0 + a, y + 12 - h / 2, b - a, h, "iv fn", r=1)
            d.line([(X0 + 200, y + 28), (X0 + 420, y + 28)], "ln navy thin")
            for a, b in ((200, 230), (380, 420)):
                d.rect(X0 + a, y + 22, b - a, 12, "iv fn", r=1)
    d.line([(X0, 20), (X0, 26 + 7 * 44 - 6)], "ln faint")
    d.line([(X1, 20), (X1, 26 + 7 * 44 - 6)], "ln faint")
    return d


# ════════════════════════════════════════════════════════════════════════
# columnar
# ════════════════════════════════════════════════════════════════════════

@fig
def columnar_rows():
    d = Diagram("columnar-rows", 880, 330,
                "Every table describing the CREs is aligned by row, and the Architecture's edges store "
                "row numbers, so joins are array indexing.")
    y0, rh = 62, 26
    # edges (left)
    d.text(20, 24, "ARCHITECTURE EDGES", "cap", "start")
    for h, x in (("src", 28), ("tgt", 64), ("w", 100), ("n", 136)):
        d.text(x, 48, h, "s", "start")
    E = [(0, 1, 12, "1.4"), (0, 4, 9, "2.2"), (3, 4, 30, "1.1"), (1, 5, 7, "0.8")]
    for k, (a, b, w, n) in enumerate(E):
        y = y0 + k * rh
        d.rect(20, y, 160, rh - 4, "bx hl" if k == 1 else "bx", r=3)
        for x, v in ((28, a), (64, b), (100, w), (136, n)):
            d.text(x, y + 15, str(v), "m", "start")
    d.text(20, y0 + 4 * rh + 14, "src and tgt are row numbers:", "s", "start")
    d.text(20, y0 + 4 * rh + 28, "no uid → vertex lookup", "s", "start")
    # loci (middle)
    lx = 262
    d.text(lx, 24, "LOCI (CREs)", "cap", "start")
    hdr = (("row", 0), ("code", 36), ("start", 82), ("end", 146), ("±", 208))
    for h, x in hdr:
        d.text(lx + x + 6, 48, h, "s", "start")
    vals = [(0, "9,800", "10,100", "."), (0, "15,200", "15,640", "+"), (0, "88,050", "88,700", "."),
            (1, "3,100", "3,350", "."), (1, "8,800", "9,420", "−"), (2, "1,010", "1,400", ".")]
    for i, v in enumerate(vals):
        y = y0 + i * rh
        d.rect(lx, y, 228, rh - 4, "bx hl" if i in (0, 4) else "bx", r=3)
        d.text(lx + 6, y + 15, str(i), "ms", "start")
        for (h, x), val in zip(hdr[1:], v):
            d.text(lx + x + 6, y + 15, str(val), "m", "start")
    for a, b, c in ((0, 3, "chr1"), (3, 5, "chr2"), (5, 6, "chrX")):
        d.line([(lx + 234, y0 + a * rh + 2), (lx + 240, y0 + a * rh + 2), (lx + 240, y0 + b * rh - 6),
                (lx + 234, y0 + b * rh - 6)], "ln faint")
        d.text(lx + 246, y0 + (a + b) / 2 * rh + 2, c, "ms", "start")
    # the highlighted edge points at rows 0 and 4
    d.path(f"M182,{y0 + rh + 11} C222,{y0 + rh + 11} 222,{y0 + 11} {lx - 6},{y0 + 11}", "ln navy")
    d.head(lx - 2, y0 + 11, 0, "hd navy", 6)
    d.path(f"M182,{y0 + rh + 11} C222,{y0 + rh + 11} 222,{y0 + 4 * rh + 11} {lx - 6},{y0 + 4 * rh + 11}", "ln navy")
    d.head(lx - 2, y0 + 4 * rh + 11, 0, "hd navy", 6)
    # same rows (right)
    ax = 570
    d.text(ax, 24, "SAME ROWS: ANNOTATIONS · SIGNAL", "cap", "start")
    for h, x in (("label", 0), ("gene", 78), ("strength", 130)):
        d.text(ax + x + 4, 48, h, "s", "start")
    ann = [("Promoter", "MYC", "3.1"), ("Intronic", "MYC", "0.8"), ("Intergenic", "—", "0"),
           ("Promoter", "SOX2", "5.2"), ("Exonic", "SOX2", "1.4"), ("Intronic", "—", "0.2")]
    for i, (a, g, st) in enumerate(ann):
        y = y0 + i * rh
        d.rect(ax, y, 188, rh - 4, "bx hl" if i in (0, 4) else "bx", r=3)
        d.text(ax + 4, y + 15, a, "m", "start")
        d.text(ax + 82, y + 15, g, "m", "start")
        d.text(ax + 144, y + 15, st, "m", "start")
        for b in range(7):
            d.rect(ax + 196 + b * 10, y + 3, 8, rh - 10, "iv fn" if (b in (3, 4) and i in (0, 3)) else "iv fl", r=1)
    d.text(ax + 230, 48, "cube[i]", "ms")
    # genes chain
    gy = 252
    d.text(20, gy - 12, "GENES · THREE TABLES LINKED BY ROW NUMBER", "cap", "start")
    for k, (t, st) in enumerate((("genes", "gene_id · gene_name"), ("transcripts", "transcript_id · gene → g"),
                                 ("features", "kind · transcript → t"))):
        x = 20 + k * 290
        d.node(x, gy, 250, 50, t, st, kind="navy", mono=True)
        if k:
            d.arrow([(x - 2, gy + 25), (x - 38, gy + 25)], "navy")
    return d


@fig
def columnar_edges():
    d = Diagram("columnar-edges", 880, 300,
                "The edge table is sorted so each chromosome's cis edges are one block and all trans edges "
                "are the last block; chrom(), cis and trans are slices of the same arrays.")
    X0, X1 = 40, 840
    blocks = [("chr1", 150), ("chr2", 135), ("chr3", 110), ("…", 160), ("chrX", 110), ("trans", 70)]
    tot = sum(w for _, w in blocks)
    x = X0
    pos = {}
    y = 70
    d.text(X0, 30, "EDGES (src < tgt), SORTED: | cis chr1 | cis chr2 | … | cis chrX | trans |", "cap", "start")
    for name, w in blocks:
        ww = w * (X1 - X0) / tot
        cls = "iv fp" if name == "trans" else ("iv fpl" if name != "chr2" else "iv fh")
        d.rect(x, y, ww - 3, 34, cls, r=3)
        d.text(x + ww / 2, y + 22, name, "m w" if name in ("trans", "chr2") else "m")
        pos[name] = (x, x + ww - 3)
        x += ww
    # brackets
    a, b = pos["chr1"][0], pos["chrX"][1]
    d.line([(a, y - 6), (a, y - 14), (b, y - 14), (b, y - 6)], "ln")
    d.text((a + b) / 2, y - 20, "A.cis = edges[: trans_start]", "ms")
    for name, lab in (("chr2", "A.chrom('chr2') = edges[lo:hi]"), ("trans", "A.trans")):
        a, b = pos[name]
        d.line([(a, y + 40), (a, y + 48), (b, y + 48), (b, y + 40)], "ln navy")
        d.text((a + b) / 2, y + 64, lab, "ms")
    d.text(X0, y + 64, "views: two offsets, no copy", "s", "start")
    # adjacency
    ax, ay, n, c = 60, 182, 10, 10
    d.text(ax, ay - 8, "AS A MATRIX", "cap", "start")
    for i in range(n):
        for j in range(n):
            blk = (i < 4 and j < 4) or (4 <= i < 7 and 4 <= j < 7) or (i >= 7 and j >= 7)
            tr = (i, j) in {(1, 8), (8, 1), (2, 5), (5, 2)}
            cls = "iv fpl" if blk else ("iv fp" if tr else "chip")
            d.rect(ax + j * c, ay + i * c, c - 1.5, c - 1.5, cls, r=1)
    d.text(ax + n * c + 14, ay + 34, "cis blocks on the diagonal,", "s", "start")
    d.text(ax + n * c + 14, ay + 48, "a few trans dots off it", "s", "start")
    for k, (t, s_) in enumerate((("A.neighbors(row)", "one CSR adjacency over all edges: trans partners included"),
                                 ("A.graph(backend=)", "the graph backend's object, built from the arrays on demand"),
                                 ("A.save(path)", "parquet: vertices + edges tables, readable from R / polars"))):
        yy = 176 + k * 38
        d.text(420, yy, t, "tm", "start")
        d.text(420, yy + 16, s_, "s", "start")
    return d


@fig
def columnar_se():
    d = Diagram("columnar-se", 880, 200,
                "HiChIP short-range read ends behave like ChIP reads: they give a coverage track and MACS3 "
                "peaks, and stitched peaks scored by signal × width are cut at the same slope-1 knee as hubs.")
    nodes = [(20, 76, 120, "allValidPairs", "HiC-Pro pairs"),
             (170, 76, 140, "shortrange_ends", "cis, |Δpos| ≤ 1 kb"),
             (360, 20, 150, "write_bed → macs3", "narrowPeak"),
             (360, 132, 150, "fragments → bigWig", "147 bp, coverage"),
             (560, 76, 130, "call_se", "stitch 12.5 kb"),
             (740, 76, 120, "SEs", "slope-1 knee")]
    for x, y, w, t, s_ in nodes:
        d.node(x, y, w, 48, t, s_, kind="purple" if t in ("call_se", "SEs") else "", mono=True)
    d.arrow([(142, 100), (168, 100)])
    d.arrow([(312, 92), (336, 92), (336, 44), (358, 44)])
    d.arrow([(312, 108), (336, 108), (336, 156), (358, 156)])
    d.arrow([(512, 44), (536, 44), (536, 92), (558, 92)])
    d.arrow([(512, 156), (536, 156), (536, 108), (558, 108)])
    d.arrow([(692, 100), (738, 100)])
    d.text(240, 150, "5′ ends, stranded", "s")
    d.text(436, 82, "peaks", "lbl")
    d.text(436, 196, "signal", "lbl")
    return d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[1:])
    for name, fn in FIGS.items():
        if only and name not in only:
            continue
        (OUT / f"{name}.svg").write_text(fn().svg())
        print(f"[diagram] {name}")


if __name__ == "__main__":
    main()
