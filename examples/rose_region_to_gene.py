"""
Replicate ROSE's *_REGION_TO_GENE.txt byte-identically using Genomeblocks-style
primitives (pyranges + pandas). The .ucsc file is parsed directly to a row-per-
occurrence table because a refseq ID can appear multiple times — on alt contigs
and even on the primary chromosome (MHC paralogs etc.) — and `Genes.transcripts`
is a dict keyed by transcript_id and would collapse them.

ROSE behavior matched here (`ROSE/bin/ROSE_geneMapper.py::mapEnhancerToGene`,
called with `-r`, byRefseq=True):

  body  = transcript span (txStart..txEnd) for every refseq row
          deduped by (chrom,start,end,strand) — see Locus.__eq__/__hash__ which
          ignore ID, so LocusCollection.__addLocus silently drops coord-dups.
  tss   = single-point TSS per refseq row, also coord-deduped.
  start_dict = TSS of the FIRST occurrence per refseq ID (ROSE's makeStartDict /
               getTSSs / refseqFromKey only ever index refseqDict[name][0]).

  overlap   = body rows intersecting the enhancer.
  proximal  = tss rows within +/-50 kb of the enhancer, minus overlap members.
  distal    = tss rows within +/-50 Mb of the enhancer (≈ same chrom), minus
              proximal-pruned members.
  candidates = overlap + proximal_pruned + distal_pruned, in ROSE's
               bin-iteration order ((strand, first_hit_bin, file_order)).
  closest   = candidates[argmin(|enh_center - start_dict[id].tss|)] — first wins.

ROSE uses two bin sizes: 500 for transcribedCollection, 50 for tssCollection.
We replicate both so the iteration order matches.
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
import pyranges as pr


ANNOT  = "/home/ualtintas/apps/ROSE/annotations/hg38_refseq.ucsc"
STITCH = ("/groups/lackgrp/projects/inv-berkay-dominance/nf-chipseq-hg38/"
          "results/ROSE/SE_K562_H3K27ac/"
          "K562_H3K27ac_REP1_peaks_SuperStitched.table.txt")
OUT    = "K562_H3K27ac_REP1_peaks_SuperStitched_REGION_TO_GENE.gb.txt"

PROX_WIN  = 50_000
DIST_WIN  = 50_000_000
WIN_BODY  = 500   # ROSE transcribedCollection winSize
WIN_TSS   = 50    # ROSE tssCollection winSize

REFSEQ_COLS = ["bin", "name", "chrom", "strand", "txStart", "txEnd",
               "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
               "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]


def read_stitched(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    return df.rename(columns={"CHROM": "Chromosome", "START": "Start", "STOP": "End"})


def load_refseq(path):
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=REFSEQ_COLS,
                     dtype={"txStart": int, "txEnd": int})
    df["order"] = np.arange(len(df))   # file order; ROSE relies on it for ties
    return df


def build_tables(refseq: pd.DataFrame):
    # body: every refseq row, deduped by (chrom,start,end,strand)
    body = refseq[["chrom", "txStart", "txEnd", "strand", "name", "order"]].copy()
    body.columns = ["Chromosome", "Start", "End", "Strand", "Name", "order"]
    body = body.drop_duplicates(["Chromosome", "Start", "End", "Strand"], keep="first")

    # ROSE's tssCollection holds ONE tssLocus per refseq ID (the first-occurrence
    # TSS via makeStartDict/getTSSs/refseqFromKey, all of which index [name][0]),
    # then dedups by coord+strand. Replicate that here.
    first = refseq.drop_duplicates("name", keep="first").copy()
    first["tss"] = np.where(first["strand"] == "+", first["txStart"], first["txEnd"])
    tss = pd.DataFrame({
        "Chromosome": first["chrom"],
        "Start":      first["tss"],
        "End":        first["tss"] + 1,
        "Strand":     first["strand"],
        "Name":       first["name"],
        "order":      first["order"],
    })
    tss = tss.drop_duplicates(["Chromosome", "Start", "End", "Strand"], keep="first")

    # start_dict: TSS of first occurrence per refseq ID (ROSE's makeStartDict)
    start_dict = first.set_index("name")[["chrom", "tss"]]

    return body.reset_index(drop=True), tss.reset_index(drop=True), start_dict


def _ordered_match_ids(matches: pd.DataFrame, query_start: int, win: int,
                        start_col: str = "Start_b") -> list[str]:
    """Return matching IDs in ROSE's LocusCollection iteration order:
    grouped by strand ('+' before '-' — __subsetHelper iterates senses=['+','-']),
    then by first bin shared with the query (max(locus.start, query.start)//win),
    then by file order. uniquified, preserving first occurrence.
    `matches` is the dataframe from pyranges join: body coords live under
    Start_b/End_b, body strand under Strand, body file-order under `order`.
    """
    if matches.empty: return []
    m = matches.copy()
    m["_strand_ord"] = (m["Strand"] == "-").astype(int)        # '+' -> 0, '-' -> 1
    m["_fhb"] = np.maximum(m[start_col], query_start) // win
    m = m.sort_values(["_strand_ord", "_fhb", "order"], kind="stable")
    seen, out = set(), []
    for n in m["Name"]:
        if n in seen: continue
        seen.add(n); out.append(n)
    return out


def map_regions(enh: pd.DataFrame, body: pd.DataFrame, tss: pd.DataFrame,
                start_dict: pd.DataFrame) -> pd.DataFrame:
    enh = enh.reset_index(drop=True).copy()
    enh["rid"] = enh.index

    enh_pr  = pr.PyRanges(enh[["Chromosome", "Start", "End", "rid"]])
    body_pr = pr.PyRanges(body)
    tss_pr  = pr.PyRanges(tss)

    # All overlap hits (every refseq row, no dedup yet — order built per region)
    ov_all = enh_pr.join(body_pr, suffix="_b").df

    # Slopped enh for proximal / distal
    def _slop(df, n):
        s = df.copy()
        s["Start"] = (s["Start"] - n).clip(lower=0)
        s["End"]   =  s["End"]   + n
        return s
    px_all = pr.PyRanges(_slop(enh, PROX_WIN)).join(tss_pr, suffix="_t").df
    di_all = pr.PyRanges(_slop(enh, DIST_WIN)).join(tss_pr, suffix="_t").df

    # Pre-group by region for cheap lookup
    def grp(df):
        return {rid: sub for rid, sub in df.groupby("rid")} if len(df) else {}
    ov_g, px_g, di_g = grp(ov_all), grp(px_all), grp(di_all)

    out_rows = []
    for row in enh.itertuples(index=False):
        rid = row.rid
        empty = ov_all.iloc[0:0]
        empty_t = px_all.iloc[0:0]
        ov_ids = _ordered_match_ids(ov_g.get(rid, empty),   row.Start,                       WIN_BODY, "Start_b")
        px_ids = _ordered_match_ids(px_g.get(rid, empty_t), max(0, row.Start - PROX_WIN),    WIN_TSS,  "Start_t")
        di_ids = _ordered_match_ids(di_g.get(rid, empty_t), max(0, row.Start - DIST_WIN),    WIN_TSS,  "Start_t")

        # ROSE pruning: proximal -= overlap, distal -= proximal_pruned
        ov_set = set(ov_ids)
        px_pruned = [x for x in px_ids if x not in ov_set]
        px_pruned_set = set(px_pruned)
        di_pruned = [x for x in di_ids if x not in px_pruned_set]

        # ROSE order: overlap + proximal_pruned + distal_pruned (may have dup IDs
        # if an overlap-id reappears in distal; argmin/index returns first).
        all_ids = ov_ids + px_pruned + di_pruned

        if not all_ids:
            closest = ""
        else:
            center = (row.Start + row.End) / 2
            dists = [abs(center - start_dict.at[i, "tss"]) for i in all_ids]
            closest = all_ids[dists.index(min(dists))]

        out_rows.append((",".join(ov_ids), ",".join(px_pruned), closest))

    enh[["OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE"]] = pd.DataFrame(
        out_rows, index=enh.index
    )
    return enh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot",  default=ANNOT)
    ap.add_argument("--stitch", default=STITCH)
    ap.add_argument("--out",    default=OUT)
    args = ap.parse_args()

    print(f"[INFO] reading refseq: {args.annot}")
    refseq = load_refseq(args.annot)
    body, tss, start_dict = build_tables(refseq)
    print(f"[INFO] refseq rows={len(refseq)}  body(dedup)={len(body)}  tss(dedup)={len(tss)}")

    print(f"[INFO] reading stitched regions: {args.stitch}")
    enh = read_stitched(args.stitch)
    print(f"[INFO] {len(enh)} regions")

    mapped = map_regions(enh, body, tss, start_dict)

    out = mapped.rename(columns={"Chromosome": "CHROM", "Start": "START", "End": "STOP"})
    out = out[["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE",
               "OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE",
               "enhancerRank", "isSuper"]]
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[INFO] wrote {args.out}")


if __name__ == "__main__":
    main()
