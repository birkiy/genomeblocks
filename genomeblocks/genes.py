"""Genes: three linked tables instead of a tree of objects.

    genes        one row per gene        (gene_id, gene_name, gene_type)
    transcripts  one row per transcript  (transcript_id, gene -> genes row)
    features     one row per exon / CDS / UTR (kind, transcript -> transcripts row,
                                               exon_number)

Each table is a :class:`~genomeblocks.Loci` on the session Genome, in the
same 0-based half-open coordinates as every other table (a GTF's 1-based
start becomes ``start - 1``), so genes and CREs compare directly. The links
are row numbers: ``transcripts['gene'][k]`` is the gene row of transcript
``k``. A gene's TSS is the 1-bp interval at its 5' end (``start`` on '+',
``end - 1`` on '-').

The annotation index (promoter / exon / UTR / gene-body intervals) is built
once, cached, and reused by every :meth:`Genes.labels` call.
"""
from __future__ import annotations

import inspect
from typing import Optional

import numpy as np

from ._table import TableMixin
from .genome import Genome
from .loci import SCODE, STRANDS, Loci
from .locus import Locus

KINDS = ("exon", "CDS", "5UTR", "3UTR")
_FEATURE = {"exon": 0, "CDS": 1, "five_prime_UTR": 2, "three_prime_UTR": 3, "UTR": 4,
            "five_prime_utr": 2, "three_prime_utr": 3, "utr": 4}      # GENCODE and Ensembl spellings
_GTF_FEATURE = {0: "exon", 1: "CDS", 2: "five_prime_UTR", 3: "three_prime_UTR"}
LABELS = np.array(["Intergenic", "Intronic", "Exonic", "3UTR", "5UTR", "Promoter-TSS"], dtype=object)


def tss_base(starts, ends, strands) -> np.ndarray:
    """0-based TSS position: start on '+' (and unstranded), end - 1 on '-'."""
    return np.where(np.asarray(strands) == 2, np.asarray(ends) - 1, np.asarray(starts))


# ── GTF / GFF3 parsing (polars, with a pandas fallback) ────────────────────

def _gtf_frame(path, keys):
    """GTF -> polars frame: chrom, feature, start, end, strand + one column per attribute."""
    import polars as pl
    pats = {k: rf'(?:^|;\s*){k} "?([^";]*)"?' for k in keys}
    df = pl.read_csv(path, separator="\t", comment_prefix="#", has_header=False,
                     columns=[0, 2, 3, 4, 6, 8], quote_char=None, infer_schema_length=0,
                     new_columns=["chrom", "feature", "start", "end", "strand", "attr"])
    return df.with_columns([pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64)]
                           + [pl.col("attr").str.extract(p, 1).alias(k) for k, p in pats.items()]
                           ).drop("attr").with_row_index("line")


def _is_gff3(path) -> bool:
    """GFF3 by extension, else by sniffing the attribute column (key=value vs key "value")."""
    name = str(path).lower().removesuffix(".gz")
    if name.endswith((".gff3", ".gff")):
        return True
    if name.endswith(".gtf"):
        return False
    import gzip
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            attr = line.rstrip("\n").split("\t")[-1]
            return "=" in attr and '"' not in attr
    return False


def _gff3_frame(path, name_k, type_k):
    """GFF3 -> the same frame _gtf_frame returns. Links follow ID/Parent; gene_id,
    transcript_id, gene_name and gene_type attributes are used when present
    (GENCODE), otherwise ID / Name / biotype (Ensembl, RefSeq)."""
    import polars as pl

    def attr(k):
        return pl.col("attr").str.extract(rf"(?:^|;){k}=([^;]*)", 1)

    def strip(c):
        return c.str.replace(r"^(gene|transcript):", "")
    df = pl.read_csv(path, separator="\t", comment_prefix="#", has_header=False,
                     columns=[0, 2, 3, 4, 6, 8], quote_char=None, infer_schema_length=0,
                     new_columns=["chrom", "feature", "start", "end", "strand", "attr"])
    df = df.with_columns(pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64),
                         attr("ID").alias("ID"), attr("Parent").alias("Parent"),
                         attr("gene_id").alias("a_gene_id"), attr("transcript_id").alias("a_tx_id"),
                         pl.coalesce(attr("exon_number"), attr("rank")).alias("exon_number"),
                         pl.coalesce(attr(name_k), attr("Name")).alias("a_name"),
                         pl.coalesce(attr(type_k), attr("biotype"), attr("gene_biotype")).alias("a_type")
                         ).with_row_index("line")
    feats = list(_FEATURE)
    genes = (df.filter(pl.col("Parent").is_null() & pl.col("ID").is_not_null()
                       & ~pl.col("feature").is_in(feats + ["chromosome", "region", "scaffold"]))
             .with_columns(pl.coalesce("a_gene_id", strip(pl.col("ID"))).alias("gene_id")))
    gid = genes.select(pl.col("ID").alias("Parent"), "gene_id", pl.col("a_name").alias("g_name"),
                       pl.col("a_type").alias("g_type"))
    tx = (df.filter(pl.col("Parent").is_not_null() & ~pl.col("feature").is_in(feats))
          .join(gid, on="Parent", how="inner", maintain_order="left")
          .with_columns(pl.coalesce("a_tx_id", strip(pl.col("ID"))).alias("transcript_id")))
    txid = tx.select(pl.col("ID").alias("Parent"), pl.col("transcript_id").alias("t_id"),
                     pl.col("gene_id").alias("t_gene"), "g_name", "g_type")
    ft = df.filter(pl.col("feature").is_in(feats)).with_columns(pl.col("Parent").str.split(","))
    kw = {"empty_as_null": True} if "empty_as_null" in inspect.signature(ft.explode).parameters else {}
    ft = ft.explode("Parent", **kw).join(txid, on="Parent", how="inner", maintain_order="left")
    out = [
        genes.select("line", "chrom", pl.lit("gene").alias("feature"), "start", "end", "strand", "gene_id",
                     pl.lit(None, pl.Utf8).alias("transcript_id"), pl.col("a_name").alias(name_k),
                     pl.col("a_type").alias(type_k), "exon_number"),
        tx.select("line", "chrom", pl.lit("transcript").alias("feature"), "start", "end", "strand", "gene_id",
                  "transcript_id", pl.col("g_name").alias(name_k), pl.col("g_type").alias(type_k), "exon_number"),
        ft.select("line", "chrom", "feature", "start", "end", "strand", pl.col("t_gene").alias("gene_id"),
                  pl.col("t_id").alias("transcript_id"), pl.col("g_name").alias(name_k),
                  pl.col("g_type").alias(type_k), "exon_number"),
    ]
    return pl.concat(out).sort("line")


def _tables_polars(path, g, name_k, type_k, chr_map):
    """Build the three tables with polars joins; only numbers cross into numpy."""
    import polars as pl
    if _is_gff3(path):
        df = _gff3_frame(path, name_k, type_k)
    else:
        df = _gtf_frame(path, ("gene_id", "transcript_id", name_k, type_k, "exon_number"))
    if chr_map:
        df = df.with_columns(pl.col("chrom").replace(chr_map))
    df = df.with_columns(pl.col("start") - 1)                  # GTF is 1-based, closed
    chroms = df["chrom"].unique(maintain_order=True).to_list()
    df = df.with_columns(
        pl.col("chrom").replace_strict(chroms, [g._add(c) for c in chroms],
                                       return_dtype=pl.Int32).alias("code"),
        pl.col("strand").replace_strict([".", "+", "-"], [0, 1, 2], default=0,
                                        return_dtype=pl.Int8).alias("s"),
        pl.col("feature").replace_strict(list(_FEATURE), list(_FEATURE.values()), default=-1,
                                         return_dtype=pl.Int8).alias("kind"))
    genes = (df.filter(pl.col("feature") == "gene").unique("gene_id", keep="first").sort("line")
             .with_columns(pl.col(name_k).fill_null(pl.col("gene_id"))))
    tx = df.filter(pl.col("feature") == "transcript")
    orphan = (tx.join(genes.select("gene_id"), on="gene_id", how="anti")
              .unique("gene_id", keep="first").sort("line")
              .with_columns(pl.col("gene_id").alias(name_k), pl.lit("").alias(type_k)))
    genes = pl.concat([genes, orphan.select(genes.columns)]).with_row_index("gene_row")
    tx = (tx.join(genes.select("gene_id", "gene_row"), on="gene_id", how="left",
                  maintain_order="left").with_row_index("tx_row"))
    ft = (df.filter(pl.col("kind") >= 0)
          .join(tx.select("gene_id", "transcript_id", "tx_row"), on=["gene_id", "transcript_id"],
                how="left", maintain_order="left")
          .with_columns(pl.col("exon_number").cast(pl.Int32, strict=False).fill_null(0)))

    def num(f, c, dt):
        return f[c].to_numpy().astype(dt, copy=False)

    G = Loci(num(genes, "code", np.int32), num(genes, "start", np.int64), num(genes, "end", np.int64),
             num(genes, "s", np.int8), genome=g,
             cols={"gene_id": genes["gene_id"].to_numpy().astype(object),
                   "gene_name": genes[name_k].to_numpy().astype(object),
                   "gene_type": genes[type_k].fill_null("").to_numpy().astype(object)})
    T = Loci(num(tx, "code", np.int32), num(tx, "start", np.int64), num(tx, "end", np.int64),
             num(tx, "s", np.int8), genome=g,
             cols={"transcript_id": tx["transcript_id"].to_numpy().astype(object),
                   "gene": num(tx, "gene_row", np.int32)})
    tx_row = ft["tx_row"].fill_null(-1).to_numpy().astype(np.int32)
    F = Loci(num(ft, "code", np.int32), num(ft, "start", np.int64), num(ft, "end", np.int64),
             num(ft, "s", np.int8), genome=g,
             cols={"kind": num(ft, "kind", np.int8), "transcript": tx_row,
                   "exon_number": num(ft, "exon_number", np.int32)})
    return G, T, F.take(tx_row >= 0)


def _read_gtf_pandas(path, keys):
    """GTF -> dict of numpy columns with pandas."""
    import pandas as pd
    cols = ["chrom", "feature", "start", "end", "strand", "attr"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                     names=cols, dtype={"chrom": str, "attr": str}, quoting=3)
    out = {c: df[c].to_numpy() for c in cols if c != "attr"}
    for k in keys:
        out[k] = df["attr"].str.extract(rf'(?:^|;\s*){k} "?([^";]*)"?', expand=False).to_numpy(object)
    return out


def _read_gff3_pandas(path, name_k, type_k):
    """GFF3 -> the columns _read_gtf_pandas returns (pandas twin of _gff3_frame)."""
    import pandas as pd
    df = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                     names=["chrom", "feature", "start", "end", "strand", "attr"],
                     dtype={"chrom": str, "attr": str}, quoting=3)
    df["line"] = np.arange(len(df))

    def attr(k):
        return df["attr"].str.extract(rf"(?:^|;){k}=([^;]*)", expand=False)

    def strip(c):
        return c.str.replace(r"^(gene|transcript):", "", regex=True)
    df["ID"], df["Parent"] = attr("ID"), attr("Parent")
    df["exon_number"] = attr("exon_number").fillna(attr("rank"))      # Ensembl GFF3 uses rank=
    df["a_name"] = attr(name_k).fillna(attr("Name"))
    df["a_type"] = attr(type_k).fillna(attr("biotype")).fillna(attr("gene_biotype"))
    a_gene, a_tx = attr("gene_id"), attr("transcript_id")
    feats = list(_FEATURE)
    is_feat = df["feature"].isin(feats)
    gm = df["Parent"].isna() & df["ID"].notna() & ~df["feature"].isin(feats + ["chromosome", "region", "scaffold"])
    genes = df[gm].assign(gene_id=a_gene[gm].fillna(strip(df["ID"][gm])))
    gid = genes[["ID", "gene_id", "a_name", "a_type"]].rename(
        columns={"ID": "Parent", "a_name": "g_name", "a_type": "g_type"})
    tm = df["Parent"].notna() & ~is_feat
    tx = df[tm].assign(transcript_id=a_tx[tm].fillna(strip(df["ID"][tm])))
    tx = tx.drop(columns=["a_name", "a_type"]).merge(gid, on="Parent", how="inner")
    txid = tx[["ID", "transcript_id", "gene_id", "g_name", "g_type"]].rename(
        columns={"ID": "Parent", "transcript_id": "t_id", "gene_id": "t_gene"})
    ft = df[is_feat].assign(Parent=df["Parent"][is_feat].str.split(",")).explode("Parent")
    ft = ft.merge(txid, on="Parent", how="inner")
    keep = ["line", "chrom", "feature", "start", "end", "strand", "gene_id", "transcript_id",
            name_k, type_k, "exon_number"]
    out = pd.concat([
        genes.assign(feature="gene", transcript_id=None).rename(columns={"a_name": name_k, "a_type": type_k}),
        tx.assign(feature="transcript").rename(columns={"g_name": name_k, "g_type": type_k}),
        ft.rename(columns={"t_gene": "gene_id", "t_id": "transcript_id", "g_name": name_k, "g_type": type_k}),
    ])[keep].sort_values("line", kind="stable")
    return {c: out[c].to_numpy() if c in ("start", "end") else out[c].to_numpy(object)
            for c in keep if c != "line"}


def _tables_numpy(d, g, gene_name_key, gene_type_key):
    """The three tables from numpy columns of GTF records (0-based starts)."""
    import pandas as pd
    feat, gid, tid = d["feature"], d["gene_id"].astype(object), d["transcript_id"].astype(object)
    codes = g.encode(d["chrom"])
    strands = pd.Series(d["strand"]).map(SCODE).fillna(0).to_numpy(np.int8)
    starts, ends = d["start"].astype(np.int64), d["end"].astype(np.int64)
    gi = np.flatnonzero(feat == "gene")
    _, keep = np.unique(gid[gi].astype(str), return_index=True)      # first record per gene_id
    gi = np.sort(gi[keep])
    ti = np.flatnonzero(feat == "transcript")
    known = pd.Index(gid[gi])
    orphan = ti[known.get_indexer(gid[ti]) < 0]
    if len(orphan):
        _, first = np.unique(gid[orphan].astype(str), return_index=True)
        gi_all = np.concatenate([gi, orphan[np.sort(first)]])
    else:
        gi_all = gi
    n_named = len(gi)
    name = np.asarray(d[gene_name_key], dtype=object)[gi_all].copy()
    gtype = np.asarray(d[gene_type_key], dtype=object)[gi_all].copy()
    nan_name = pd.isna(name)
    name[nan_name] = gid[gi_all][nan_name]
    name[n_named:] = gid[gi_all][n_named:]              # genes known only from transcripts
    gtype[pd.isna(gtype)] = ""
    gtype[n_named:] = ""
    G = Loci(codes[gi_all], starts[gi_all], ends[gi_all], strands[gi_all], genome=g,
             cols={"gene_id": gid[gi_all], "gene_name": name, "gene_type": gtype})
    gene_index = pd.Index(gid[gi_all])
    T = Loci(codes[ti], starts[ti], ends[ti], strands[ti], genome=g,
             cols={"transcript_id": tid[ti],
                   "gene": gene_index.get_indexer(gid[ti]).astype(np.int32)})
    kind = pd.Series(feat).map(_FEATURE).fillna(-1).to_numpy(np.int8)
    fi = np.flatnonzero(kind >= 0)
    tx_row = pd.Index(tid[ti]).get_indexer(tid[fi]).astype(np.int32)
    en = pd.to_numeric(pd.Series(d["exon_number"][fi]), errors="coerce").fillna(0).to_numpy(np.int32)
    F = Loci(codes[fi], starts[fi], ends[fi], strands[fi], genome=g,
             cols={"kind": kind[fi].copy(), "transcript": tx_row, "exon_number": en})
    return G, T, F.take(tx_row >= 0)


def _resolve_generic_utr(F, T):
    """Ensembl-style 'UTR' rows: 5' if before the transcript's first CDS base."""
    k = F.cols["kind"]
    generic = np.flatnonzero(k == 4)
    if not len(generic):
        return
    tx_row = F.cols["transcript"]
    is_cds = k == 1
    big = np.iinfo(np.int64).max
    cds_start = np.full(len(T), big)
    np.minimum.at(cds_start, tx_row[is_cds], F.starts[is_cds])
    cs = cds_start[tx_row[generic]]
    plus = F.strands[generic] != 2
    before = F.ends[generic] <= cs
    k[generic] = np.where(cs == big, 2, np.where(plus, np.where(before, 2, 3), np.where(before, 3, 2)))


# ── UCSC genePred tables (refGene / ncbiRefSeq / knownGene) ────────────────

def _tables_ucsc(path, g, chr_map, keep_alt_contigs):
    import pandas as pd
    raw = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, quoting=3)
    layout = ("a UCSC genePred table: [bin] name chrom strand txStart txEnd cdsStart cdsEnd exonCount "
              "exonStarts exonEnds [...] (refGene, ncbiRefSeq, knownGene, genePredExt) or refFlat "
              "(geneName name chrom strand ...); for GTF / GFF3 use Genes.make()")
    ncol = raw.shape[1]
    # the strand column fixes the layout: it is column 2 (name chrom strand ...), 3 with a
    # leading bin column (bin name chrom strand ...) or 3 for refFlat (geneName name chrom strand ...)
    sc = next((i for i in range(min(4, ncol)) if raw[i].isin(["+", "-", "."]).all()), None)
    if sc is None or sc < 2 or ncol < sc + 8:
        raise ValueError(f"{path}: expected {layout}")
    try:
        tx_s, tx_e = raw[sc + 1].astype(np.int64).to_numpy(), raw[sc + 2].astype(np.int64).to_numpy()
        cds_s, cds_e = raw[sc + 3].astype(np.int64).to_numpy(), raw[sc + 4].astype(np.int64).to_numpy()
        n_ex = raw[sc + 5].astype(np.int64).to_numpy()
    except ValueError as e:
        raise ValueError(f"{path}: expected {layout}; {e}") from None
    with_bin = sc == 3 and raw[0].str.fullmatch(r"\d+").all()
    name = raw[sc - 2].to_numpy(object)
    chrom = raw[sc - 1].to_numpy(object)
    if chr_map:
        chrom = np.array([chr_map.get(c, c) for c in chrom], dtype=object)
    strand = raw[sc].to_numpy(object)
    ex_s = raw[sc + 6].to_numpy(object)
    ex_e = raw[sc + 7].to_numpy(object)
    if sc == 3 and not with_bin:                           # refFlat: the symbol comes first
        sym = raw[0].fillna("").to_numpy(object)
    elif ncol >= sc + 13:                                   # genePredExt: score name2 cdsStartStat ...
        sym = raw[sc + 9].fillna("").to_numpy(object)
    else:                                                   # knownGene and plain genePred: no symbol
        sym = name.copy()
    sym = np.where(sym == "", name, sym)
    # one gene per (symbol, chromosome); the first chromosome seen is the primary one
    first_chrom = {}
    gene_key = []
    keep = np.ones(len(raw), bool)
    for i, (s, c) in enumerate(zip(sym.tolist(), chrom.tolist())):
        p = first_chrom.setdefault(s, c)
        if c == p:
            gene_key.append(s)
        elif keep_alt_contigs:
            gene_key.append(f"{s}__{c}")
        else:
            gene_key.append(None)
            keep[i] = False
    gene_key = np.array(gene_key, dtype=object)
    rows = np.flatnonzero(keep)
    gk = gene_key[rows]
    uniq, inv = np.unique(gk.astype(str), return_inverse=True)
    first_row = np.full(len(uniq), len(rows))
    np.minimum.at(first_row, inv, np.arange(len(rows)))
    order = np.argsort(first_row)                        # genes in file order
    rank = np.empty(len(order), np.int64)
    rank[order] = np.arange(len(order))
    gene_of_tx = rank[inv]
    ng = len(uniq)
    gs = np.full(ng, np.iinfo(np.int64).max)
    ge = np.full(ng, np.iinfo(np.int64).min)
    np.minimum.at(gs, gene_of_tx, tx_s[rows])
    np.maximum.at(ge, gene_of_tx, tx_e[rows])
    grow = rows[first_row[order]]
    G = Loci(g.encode(chrom[grow]), gs, ge, np.array([SCODE.get(x, 0) for x in strand[grow]], np.int8), genome=g,
             cols={"gene_id": uniq[order].astype(object), "gene_name": sym[grow].astype(object),
                   "gene_type": np.full(ng, "", dtype=object)})
    # transcript ids: disambiguate repeats within a gene with __2, __3, ...
    tids = []
    seen = {}
    for k, t in zip(gene_of_tx.tolist(), name[rows].tolist()):
        key = (k, t)
        n = seen.get(key, 0) + 1
        seen[key] = n
        tids.append(t if n == 1 else f"{t}__{n}")
    T = Loci(g.encode(chrom[rows]), tx_s[rows], tx_e[rows],
             np.array([SCODE.get(x, 0) for x in strand[rows]], np.int8), genome=g,
             cols={"transcript_id": np.array(tids, dtype=object), "gene": gene_of_tx.astype(np.int32)})
    # exons, then CDS / UTR pieces cut from them
    fs, fe, ft, fk, fn = [], [], [], [], []
    for k, i in enumerate(rows.tolist()):
        es = np.array([int(x) for x in str(ex_s[i]).rstrip(",").split(",") if x], np.int64)
        ee = np.array([int(x) for x in str(ex_e[i]).rstrip(",").split(",") if x], np.int64)
        m = len(es)
        num = np.arange(1, m + 1) if strand[i] != "-" else np.arange(m, 0, -1)
        fs.append(es); fe.append(ee); ft.append(np.full(m, k)); fk.append(np.zeros(m, np.int8)); fn.append(num)
        cs, ce = cds_s[i], cds_e[i]
        if ce <= cs:
            continue
        a, b = np.maximum(es, cs), np.minimum(ee, ce)
        ok = a < b
        fs.append(a[ok]); fe.append(b[ok]); ft.append(np.full(ok.sum(), k)); fk.append(np.ones(ok.sum(), np.int8))
        fn.append(num[ok])
        left = es < cs
        right = ee > ce
        lk, rk = (2, 3) if strand[i] != "-" else (3, 2)
        fs.append(es[left]); fe.append(np.minimum(ee[left], cs)); ft.append(np.full(left.sum(), k))
        fk.append(np.full(left.sum(), lk, np.int8)); fn.append(num[left])
        fs.append(np.maximum(es[right], ce)); fe.append(ee[right]); ft.append(np.full(right.sum(), k))
        fk.append(np.full(right.sum(), rk, np.int8)); fn.append(num[right])
    ftx = np.concatenate(ft).astype(np.int32) if ft else np.zeros(0, np.int32)
    F = Loci(T.codes[ftx], np.concatenate(fs) if fs else np.zeros(0, np.int64),
             np.concatenate(fe) if fe else np.zeros(0, np.int64), T.strands[ftx], genome=g,
             cols={"kind": np.concatenate(fk) if fk else np.zeros(0, np.int8), "transcript": ftx,
                   "exon_number": np.concatenate(fn).astype(np.int32) if fn else np.zeros(0, np.int32)})
    return G, T, F


# ── a gene as one object ───────────────────────────────────────────────────

class GeneView(Locus):
    """One gene as a Locus, with its transcripts and exons one attribute away."""

    def __init__(self, genes: "Genes", i: int):
        self.__dict__.update(_G=genes, _i=i)

    _L = property(lambda s: s._G.genes)
    chrom = property(lambda s: s._L.genome.names[s._L.codes[s._i]])
    start = property(lambda s: int(s._L.starts[s._i]))
    end = property(lambda s: int(s._L.ends[s._i]))
    strand = property(lambda s: STRANDS[s._L.strands[s._i]])
    gene_id = property(lambda s: s._L.cols["gene_id"][s._i])
    gene_name = property(lambda s: s._L.cols["gene_name"][s._i])
    gene_type = property(lambda s: s._L.cols["gene_type"][s._i])

    @property
    def row(self) -> int:
        return self.__dict__["_i"]

    @property
    def tss(self) -> Locus:
        """The 1-bp TSS (on the gene's strand)."""
        t = int(tss_base(self.start, self.end, 2 if self.strand == "-" else 1))
        return Locus(self.chrom, t, t + 1, self.strand)

    @property
    def transcripts(self) -> Loci:
        T = self._G.transcripts
        return T.take(np.flatnonzero(T.cols["gene"] == self._i))

    @property
    def canonical(self) -> Optional[str]:
        """transcript_id chosen by :meth:`Genes.select_isoforms` (None if not run)."""
        c = self._G.genes.cols.get("canonical")
        if c is None or c[self._i] < 0:
            return None
        return self._G.transcripts.cols["transcript_id"][c[self._i]]

    @property
    def exons(self) -> Loci:
        F, T = self._G.features, self._G.transcripts
        tx = np.flatnonzero(T.cols["gene"] == self._i)
        return F.take(np.flatnonzero(np.isin(F.cols["transcript"], tx) & (F.cols["kind"] == 0)))

    def __repr__(self):
        return f"Gene({self.gene_name}, {self.chrom}:{self.start:,}-{self.end:,}({self.strand}), {self.gene_type})"


class Genes(TableMixin):
    """Genes, transcripts and features as linked tables.

    The object itself behaves like its genes table: ``len``, ``shape``,
    ``columns``, ``head()``, ``describe()`` and the Arrow / dataframe
    protocols read the genes; ``to_pandas('transcripts')`` and friends give
    the others."""

    def __init__(self, genes: Loci, transcripts: Loci, features: Loci, *, promoter_r: int = 1000,
                 filename: Optional[str] = None):
        need = {"genes": ("gene_id", "gene_name", "gene_type"), "transcripts": ("transcript_id", "gene"),
                "features": ("kind", "transcript", "exon_number")}
        for what, L in (("genes", genes), ("transcripts", transcripts), ("features", features)):
            if not isinstance(L, Loci):
                raise TypeError(f"Genes: {what} must be a Loci, got {type(L).__name__}")
            missing = [k for k in need[what] if k not in L.cols]
            if missing:
                raise ValueError(f"Genes: the {what} table lacks the column(s) {missing}; "
                                 f"build with Genes.make / make_ucsc / from_frame")
        self.genes, self.transcripts, self.features = genes, transcripts, features
        self.promoter_r = promoter_r
        self.filename = filename
        self._annot = None
        self._by_key = None

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, filename: str, *, gene_name_key: str = "gene_name", gene_type_key: str = "gene_type",
             promoter_r: int = 1000, genome: Optional[Genome] = None, chr_map: Optional[dict] = None,
             cre=None, bw=None, r=None, kw=None, backend: Optional[str] = None) -> "Genes":
        """Parse a GTF or GFF3 (GENCODE / Ensembl / RefSeq style, ``.gz`` too).

        ``chr_map`` renames chromosomes while parsing (``{'1': 'chr1'}``).
        ``cre`` / ``bw`` / ``r`` / ``kw`` run :meth:`select_isoforms` right after
        parsing. ``backend`` picks the parser (polars or pandas)."""
        from .backends import resolve
        g = genome if genome is not None else Genome()
        if resolve("tables", backend) == "polars":
            G, T, F = _tables_polars(filename, g, gene_name_key, gene_type_key, chr_map)
        else:
            if _is_gff3(filename):
                d = _read_gff3_pandas(filename, gene_name_key, gene_type_key)
            else:
                d = _read_gtf_pandas(filename, ("gene_id", "transcript_id", gene_name_key,
                                                gene_type_key, "exon_number"))
            if chr_map:
                d["chrom"] = np.array([chr_map.get(c, c) for c in d["chrom"]], dtype=object)
            d["start"] = d["start"].astype(np.int64) - 1          # GTF is 1-based, closed
            G, T, F = _tables_numpy(d, g, gene_name_key, gene_type_key)
        _resolve_generic_utr(F, T)
        out = cls(G, T, F, promoter_r=promoter_r, filename=str(filename))
        if cre is not None or bw is not None:
            out.select_isoforms(cre, bw, r=r, **(kw or {}))
        return out

    @classmethod
    def make_ucsc(cls, filename: str, *, promoter_r: int = 1000, genome: Optional[Genome] = None,
                  chr_map: Optional[dict] = None, keep_alt_contigs: bool = False, cre=None, bw=None,
                  r=None, kw=None) -> "Genes":
        """Parse a UCSC genePred table (refGene, ncbiRefSeq, knownGene; with or
        without the leading ``bin`` column). One gene per symbol and chromosome;
        copies on other chromosomes (alt haplotypes) are skipped unless
        ``keep_alt_contigs`` (then keyed ``SYMBOL__chrom``)."""
        g = genome if genome is not None else Genome()
        G, T, F = _tables_ucsc(filename, g, chr_map, keep_alt_contigs)
        out = cls(G, T, F, promoter_r=promoter_r, filename=str(filename))
        if cre is not None or bw is not None:
            out.select_isoforms(cre, bw, r=r, **(kw or {}))
        return out

    @classmethod
    def from_frame(cls, df, *, one_based: Optional[bool] = None, promoter_r: int = 1000,
                   genome: Optional[Genome] = None) -> "Genes":
        """From a GTF-like table — one row per gene / transcript / exon / CDS / UTR
        record — e.g. ``pyranges.read_gtf(...)``, gtfparse, or a polars frame.

        Columns are found by name (Chromosome / seqname / chrom, Feature /
        feature, Start / start, End / end, Strand / strand, gene_id,
        transcript_id, gene_name, gene_type, exon_number). ``one_based`` says
        whether starts are GTF-style 1-based; by default PyRanges objects (and
        capitalised ``Start``) are taken as 0-based and anything else as 1-based.
        """
        from .interop import _as_pandas, _find
        is_pr = (type(df).__module__ or "").startswith("pyranges")
        pdf = _as_pandas(df)
        cols = list(pdf.columns)
        c = _find(cols, ("chromosome", "chrom", "seqname", "seqnames", "seqid", "chr"), None)
        f = _find(cols, ("feature", "type", "feature_type"), None)
        s = _find(cols, ("start",), None)
        e = _find(cols, ("end",), None)
        d = _find(cols, ("strand",), None)
        missing = [w for w, k in (("chrom / seqname / Chromosome", c), ("feature / type", f),
                                  ("start", s), ("end", e)) if k is None]
        if missing:
            raise ValueError(f"Genes.from_frame needs the columns {', '.join(missing)} (one GTF record "
                             f"per row); got columns {cols}")
        for req in ("gene_id", "transcript_id"):
            k = _find(cols, (req,), None)
            if k is None or pdf[k].isna().all():
                raise ValueError(f"Genes.from_frame needs a {req!r} column (gene / transcript records "
                                 f"are linked by it); got columns {cols}")
        if one_based is None:
            one_based = not (is_pr or s == "Start")

        def col(name):
            k = _find(cols, (name,), None)
            return pdf[k].to_numpy(object) if k is not None else np.full(len(pdf), None, dtype=object)
        data = {"chrom": pdf[c].astype(str).to_numpy(object), "feature": pdf[f].astype(str).to_numpy(object),
                "start": pdf[s].to_numpy(np.int64) - (1 if one_based else 0), "end": pdf[e].to_numpy(np.int64),
                "strand": pdf[d].astype(str).to_numpy(object) if d else np.full(len(pdf), ".", dtype=object),
                "gene_id": col("gene_id"), "transcript_id": col("transcript_id"), "gene_name": col("gene_name"),
                "gene_type": col("gene_type") if _find(cols, ("gene_type",), None) else col("gene_biotype"),
                "exon_number": col("exon_number")}
        g = genome if genome is not None else Genome()
        G, T, F = _tables_numpy(data, g, "gene_name", "gene_type")
        _resolve_generic_utr(F, T)
        return cls(G, T, F, promoter_r=promoter_r)

    # ── ATAC-supported isoform selection ─────────────────────────────────
    def select_isoforms(self, cre=None, bw=None, *, r=None, agg="max", min_signal=0.0, min_frac=0.5,
                        rank="longest", collapse=True, verbose=True) -> "Genes":
        """Point each gene at the isoform its cells actually use.

        A transcript is *supported* when its TSS window (TSS ± ``r``) overlaps a
        peak in ``cre`` and, with ``bw``, its window score is > ``min_signal``
        and >= ``min_frac`` x the best candidate score of its gene (max over the
        bigWigs). Per gene the winner is the longest supported isoform
        (``rank='longest'``) or the strongest (``'signal'``), the other breaking
        ties, then the transcript id; genes with no supported isoform fall back
        to all of theirs. ``collapse`` moves the gene body (and so its TSS)
        onto the winner. ``cre`` is anything :func:`~genomeblocks.as_loci` takes
        (a list of peak sets is merged).

        Writes ``transcripts['tss_score']`` (NaN = not scored), ``['tss_support']``
        and ``genes['canonical']`` (transcript row, -1 = none).
        """
        if rank not in ("longest", "signal"):
            raise ValueError(f"rank must be 'longest' or 'signal', got {rank!r}")
        T, G = self.transcripts, self.genes
        r = self.promoter_r if r is None else r
        peaks = _peaks(cre, T.genome)
        import os
        if bw is None:
            bigwigs = []
        elif isinstance(bw, (str, os.PathLike)) or not hasattr(bw, "__iter__"):
            bigwigs = [bw]                                   # one path or one open handle
        elif isinstance(bw, dict):
            bigwigs = list(bw.values())
        else:
            bigwigs = list(bw)
        bigwigs = [str(b) if isinstance(b, os.PathLike) else b for b in bigwigs]
        if peaks is None and not bigwigs:
            raise ValueError("select_isoforms() needs `cre` (peaks) and/or `bw` (bigWig).")
        n = len(T)
        pos = tss_base(T.starts, T.ends, T.strands)
        ws, we = np.maximum(0, pos - r), pos + r + 1
        cand = (Loci(T.codes, ws, we, genome=T.genome).overlap_any(peaks) if peaks is not None
                else np.ones(n, bool))
        score = np.full(n, np.nan)
        if bigwigs and cand.any():
            idx = np.flatnonzero(cand)
            key = np.stack([T.codes[idx].astype(np.int64), ws[idx], we[idx]], 1)
            uk, inv = np.unique(key, axis=0, return_inverse=True)
            wins = Loci(uk[:, 0], uk[:, 1], uk[:, 2], genome=T.genome)
            cube = wins.signal(bigwigs, n_bins=1, span=True, agg=agg, progress=False, verbose=False)
            score[idx] = np.asarray(cube[:, :, 0], float).max(axis=1)[inv.ravel()]
        support = cand.copy()
        gene = T.cols["gene"].astype(np.int64)
        ng = len(G)
        if bigwigs:
            sc = np.where(support, score, -np.inf)
            best = np.full(ng, -np.inf)
            np.maximum.at(best, gene[support], sc[support])
            best[~np.isfinite(best)] = 0.0
            with np.errstate(invalid="ignore"):
                support &= (score > min_signal) & (score >= min_frac * best[gene])
        has_sup = np.bincount(gene[support], minlength=ng) > 0
        pool = support | ~has_sup[gene]
        length = T.ends - T.starts
        s0 = np.nan_to_num(score, nan=0.0)
        tid_rank = np.unique(T.cols["transcript_id"].astype(str), return_inverse=True)[1].ravel()
        first_key, second_key = (length, s0) if rank == "longest" else (s0, length)
        cand_rows = np.flatnonzero(pool)
        o = np.lexsort((tid_rank[cand_rows], -second_key[cand_rows], -first_key[cand_rows], gene[cand_rows]))
        rows = cand_rows[o]
        g_sorted = gene[rows]
        winners = rows[np.r_[True, g_sorted[1:] != g_sorted[:-1]]] if len(rows) else rows
        canonical = np.full(ng, -1, np.int64)
        canonical[gene[winners]] = winners
        T.cols["tss_score"], T.cols["tss_support"] = score, support
        G.cols["canonical"] = canonical
        if collapse and len(winners):
            G.starts, G.ends = G.starts.copy(), G.ends.copy()
            G.starts[gene[winners]] = T.starts[winners]
            G.ends[gene[winners]] = T.ends[winners]
            G._dirty()
        self._annot = None
        if verbose:
            n_genes = int((np.bincount(gene, minlength=ng) > 0).sum())
            n_sup = int(has_sup.sum())
            print(f"[INFO] Isoform support: {n_sup}/{n_genes} genes with an open TSS "
                  f"({int(support.sum())}/{n} isoforms), {n_genes - n_sup} fell back to the longest isoform.")
            if n_sup == 0:
                print("[WARN] No isoform was supported — check that the peaks/bigwig use the "
                      "same chromosome names as the annotation.")
        return self

    # ── lookups ───────────────────────────────────────────────────────────
    def __len__(self):
        return len(self.genes)

    def __iter__(self):
        for i in range(len(self.genes)):
            yield GeneView(self, i)

    def _key_index(self):
        if self._by_key is None:
            ids = self.genes.cols["gene_id"].tolist()
            names = self.genes.cols["gene_name"].tolist()
            by = {n: i for i, n in enumerate(names)}
            by.update({g: i for i, g in enumerate(ids)})
            self._by_key = by
        return self._by_key

    def __getitem__(self, key) -> GeneView:
        """``genes['TP53']``, ``genes['ENSG...']`` or ``genes[row]``."""
        if isinstance(key, (int, np.integer)):
            return GeneView(self, int(key))
        return GeneView(self, self._key_index()[key])

    def __contains__(self, key):
        return key in self._key_index()

    def find(self, pattern: str):
        """Genes whose name contains ``pattern`` (case-insensitive)."""
        names = self.genes.cols["gene_name"].astype(str)
        hit = np.char.find(np.char.lower(names), pattern.lower()) >= 0
        return [GeneView(self, int(i)) for i in np.flatnonzero(hit)]

    def rows(self, names) -> np.ndarray:
        """Gene rows for gene names / ids (KeyError on an unknown one)."""
        ix = self._key_index()
        return np.array([ix[n] for n in ([names] if isinstance(names, str) else names)], np.int64)

    # ── TSS and the annotation index ─────────────────────────────────────
    def get_tss(self, gene_type=None) -> Loci:
        """1-bp TSS of every gene (optionally only ``gene_type`` ones), on the
        gene's strand, with ``gene_name``, ``gene_id`` and ``gene`` (row) columns."""
        G = self.genes
        rows = np.arange(len(G))
        if gene_type is not None:
            types = [gene_type] if isinstance(gene_type, str) else list(gene_type)
            rows = rows[np.isin(G.cols["gene_type"].astype(str), types)]
        t = tss_base(G.starts[rows], G.ends[rows], G.strands[rows])
        return Loci(G.codes[rows], t, t + 1, G.strands[rows], genome=G.genome,
                    cols={"gene_name": G.cols["gene_name"][rows], "gene_id": G.cols["gene_id"][rows],
                          "gene": rows})

    @property
    def annot(self) -> dict:
        """The annotation index: body, prom (TSS ± promoter_r), exon, utr5, utr3 (merged)."""
        if self._annot is None:
            F = self.features
            k = F.cols["kind"]
            self._annot = {
                "body": self.genes,
                "prom": self.get_tss().slop(self.promoter_r).merge(),
                "exon": F.take(k == 0).merge(),
                "utr5": F.take(k == 2).merge(),
                "utr3": F.take(k == 3).merge(),
            }
        return self._annot

    def labels(self, L, *, backend: Optional[str] = None) -> np.ndarray:
        """Region label code per row of ``L`` (index into ``LABELS``); the
        highest-priority class wins: Promoter-TSS > 5UTR > 3UTR > Exonic >
        Intronic > Intergenic."""
        from .interop import as_loci
        L = as_loci(L, genome=self.genes.genome)
        a = self.annot
        out = np.zeros(len(L), np.int8)
        for code, key in ((1, "body"), (2, "exon"), (3, "utr3"), (4, "utr5"), (5, "prom")):
            out[L.overlap_any(a[key], backend=backend)] = code
        return out

    def annotations(self, L, *, backend: Optional[str] = None):
        """Region label per row of ``L`` as a DataFrame (uid, annotation)."""
        import pandas as pd
        from .interop import as_loci
        L = as_loci(L, genome=self.genes.genome)
        return pd.DataFrame({"uid": L.uid, "annotation": LABELS[self.labels(L, backend=backend)]})

    def nearest_tss(self, L, *, backend: Optional[str] = None):
        """(gene_name, distance) of the nearest promoter window (TSS ± promoter_r)
        per row of ``L``; ('', -1) where the chromosome has no gene."""
        from .interop import as_loci
        L = as_loci(L, genome=self.genes.genome)
        win = self.get_tss().slop(self.promoter_r)
        j, dist = L.nearest(win, backend=backend)
        names = np.where(j >= 0, win.cols["gene_name"][np.maximum(j, 0)], "")
        return names.astype(object), dist

    def nearest_genes(self, L, *, backend: Optional[str] = None):
        """DataFrame(uid, gene_name, distance) for the rows of ``L`` that have a gene
        on their chromosome."""
        import pandas as pd
        from .interop import as_loci
        L = as_loci(L, genome=self.genes.genome)
        names, dist = self.nearest_tss(L, backend=backend)
        ok = dist >= 0
        return pd.DataFrame({"uid": L.uid[ok], "gene_name": names[ok], "distance": dist[ok]})

    # ── export ────────────────────────────────────────────────────────────
    def to_pandas(self, table: str = "genes"):
        """One table as a DataFrame with its links spelled out: 'genes'
        (+ canonical transcript), 'transcripts' (+ gene_id / gene_name),
        'features' (+ feature name, transcript_id, gene_id)."""
        G, T, F = self.genes, self.transcripts, self.features
        if table == "genes":
            df = G.to_pandas()
            if "canonical" in G.cols:
                c = G.cols["canonical"]
                df["canonical"] = np.where(c >= 0, T.cols["transcript_id"][np.maximum(c, 0)], None)
            return df
        if table == "transcripts":
            df = T.to_pandas()
            gr = T.cols["gene"]
            df["gene_id"] = G.cols["gene_id"][gr]
            df["gene_name"] = G.cols["gene_name"][gr]
            return df
        if table == "features":
            df = F.to_pandas()
            tr = F.cols["transcript"]
            df["feature"] = np.array(["exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"],
                                     dtype=object)[F.cols["kind"]]
            df["transcript_id"] = T.cols["transcript_id"][tr]
            df["gene_id"] = G.cols["gene_id"][T.cols["gene"][tr]]
            return df
        raise ValueError("table must be 'genes', 'transcripts' or 'features'")

    def to_polars(self, table: str = "genes"):
        import polars as pl
        return pl.from_pandas(self.to_pandas(table))

    def to_arrow(self, table: str = "genes"):
        import pyarrow as pa
        return pa.Table.from_pandas(self.to_pandas(table), preserve_index=False)

    def to_gtf(self, path: str) -> str:
        """Write the tables back as a GTF (1-based), genes then transcripts then features."""
        G, T, F = self.genes, self.transcripts, self.features
        gid, gname, gtype = G.cols["gene_id"], G.cols["gene_name"], G.cols["gene_type"]
        tid, tgene = T.cols["transcript_id"], T.cols["gene"]
        kinds = np.array([_GTF_FEATURE.get(k, "UTR") for k in range(5)], dtype=object)
        with open(path, "w") as f:
            for i in range(len(G)):
                f.write(f"{G.genome.names[G.codes[i]]}\tgenomeblocks\tgene\t{G.starts[i] + 1}\t{G.ends[i]}\t.\t"
                        f"{STRANDS[G.strands[i]]}\t.\tgene_id \"{gid[i]}\"; gene_name \"{gname[i]}\"; "
                        f"gene_type \"{gtype[i]}\";\n")
            for k in range(len(T)):
                gi = tgene[k]
                f.write(f"{T.genome.names[T.codes[k]]}\tgenomeblocks\ttranscript\t{T.starts[k] + 1}\t{T.ends[k]}\t.\t"
                        f"{STRANDS[T.strands[k]]}\t.\tgene_id \"{gid[gi]}\"; transcript_id \"{tid[k]}\"; "
                        f"gene_name \"{gname[gi]}\"; gene_type \"{gtype[gi]}\";\n")
            ftr, fkind, fnum = F.cols["transcript"], F.cols["kind"], F.cols["exon_number"]
            for j in range(len(F)):
                k = ftr[j]
                gi = tgene[k]
                f.write(f"{F.genome.names[F.codes[j]]}\tgenomeblocks\t{kinds[fkind[j]]}\t{F.starts[j] + 1}\t"
                        f"{F.ends[j]}\t.\t{STRANDS[F.strands[j]]}\t.\tgene_id \"{gid[gi]}\"; "
                        f"transcript_id \"{tid[k]}\"; gene_name \"{gname[gi]}\"; exon_number \"{fnum[j]}\";\n")
        return path

    def to_bed12(self, path: Optional[str] = None, *, canonical: bool = True):
        """One BED12 line per gene — its canonical isoform when selected, else its
        longest — with exons as blocks (what IGV / UCSC draw)."""
        T, F, G = self.transcripts, self.features, self.genes
        tx = self.representative(canonical=canonical)
        ex = F.cols["kind"] == 0
        f_tx, f_s, f_e = F.cols["transcript"][ex], F.starts[ex], F.ends[ex]
        o = np.lexsort((f_s, f_tx))
        f_tx, f_s, f_e = f_tx[o], f_s[o], f_e[o]
        lo, hi = np.searchsorted(f_tx, tx), np.searchsorted(f_tx, tx, side="right")
        # thickStart / thickEnd: the CDS span of the transcript (thickStart == thickEnd when non-coding)
        cds = F.cols["kind"] == 1
        cs = np.full(len(T), np.iinfo(np.int64).max)
        ce = np.full(len(T), -1, np.int64)
        if cds.any():
            np.minimum.at(cs, F.cols["transcript"][cds], F.starts[cds])
            np.maximum.at(ce, F.cols["transcript"][cds], F.ends[cds])
        names = G.cols["gene_name"]
        lines = []
        for t, a, b in zip(tx.tolist(), lo.tolist(), hi.tolist()):
            s, e = int(T.starts[t]), int(T.ends[t])
            bs, be = (f_s[a:b], f_e[a:b]) if b > a else (np.array([s]), np.array([e]))
            ts, te = (int(cs[t]), int(ce[t])) if ce[t] > cs[t] else (e, e)
            lines.append(f"{T.genome.names[T.codes[t]]}\t{s}\t{e}\t{names[T.cols['gene'][t]]}\t0\t"
                         f"{STRANDS[T.strands[t]]}\t{ts}\t{te}\t0\t{len(bs)}\t"
                         f"{','.join(str(int(x)) for x in be - bs)}\t{','.join(str(int(x)) for x in bs - s)}")
        text = "\n".join(lines)
        if path is None:
            return text
        with open(path, "w") as f:
            f.write(text + "\n")
        return path

    def representative(self, *, canonical: bool = True) -> np.ndarray:
        """One transcript row per gene: the canonical one (when selected), else the longest."""
        T = self.transcripts
        if not len(T):
            return np.zeros(0, np.int64)
        tlen = T.ends - T.starts
        order = np.lexsort((-tlen, T.cols["gene"]))
        g_sorted = T.cols["gene"][order]
        tx = order[np.r_[True, g_sorted[1:] != g_sorted[:-1]]]
        c = self.genes.cols.get("canonical")
        if canonical and c is not None:
            have = c[T.cols["gene"][tx]]
            tx = np.where(have >= 0, have, tx)
        return tx.astype(np.int64)

    def save(self, path: str):
        """``path/`` gets genes / transcripts / features parquet tables + meta.json."""
        import json
        import os
        os.makedirs(path, exist_ok=True)
        for name in ("genes", "transcripts", "features"):
            getattr(self, name).save(os.path.join(path, f"{name}.parquet"))
        with open(os.path.join(path, "meta.json"), "w") as f:
            json.dump({"promoter_r": self.promoter_r, "filename": self.filename,
                       "format": "genomeblocks.Genes/2", "coordinates": "0-based half-open"}, f)

    @classmethod
    def load(cls, path: str, *, genome: Optional[Genome] = None) -> "Genes":
        import json
        import os
        g = genome if genome is not None else Genome()
        t = [Loci.load(os.path.join(path, f"{n}.parquet"), genome=g)
             for n in ("genes", "transcripts", "features")]
        with open(os.path.join(path, "meta.json")) as f:
            meta = json.load(f)
        return cls(*t, promoter_r=meta["promoter_r"], filename=meta.get("filename"))

    # ── like a DataFrame (the genes table) ────────────────────────────────
    @property
    def columns(self) -> list:
        return self.genes.columns

    def head(self, n: int = 5):
        return self.to_pandas("genes").head(n)

    def tail(self, n: int = 5):
        return self.to_pandas("genes").tail(n)

    def describe(self):
        """One-table summary of the three tables."""
        import pandas as pd
        c = self.counts()
        G = self.genes
        rows = [("genes", c["genes"]), ("transcripts", c["transcripts"]), ("exons", c["exons"]),
                ("CDS", c["CDS"]), ("UTR", c["UTR"]),
                ("chromosomes", int(len(np.unique(G.codes))) if len(G) else 0),
                ("gene types", int(len(np.unique(G.cols["gene_type"].astype(str)))) if len(G) else 0),
                ("canonical isoforms", "selected" if "canonical" in G.cols else "not selected"),
                ("promoter_r", self.promoter_r), ("coordinates", "0-based, half-open")]
        return pd.DataFrame(rows, columns=["", "value"]).set_index("")

    summary = describe

    # ── display ───────────────────────────────────────────────────────────
    def counts(self) -> dict:
        k = self.features.cols["kind"]
        return {"genes": len(self.genes), "transcripts": len(self.transcripts),
                "exons": int((k == 0).sum()), "CDS": int((k == 1).sum()),
                "UTR": int((k >= 2).sum())}

    def __repr__(self):
        c = self.counts()
        return f"Genes({c['genes']:,} genes, {c['transcripts']:,} transcripts, {c['exons']:,} exons)"

    def _repr_html_(self):
        from ._display import kv_html
        c = self.counts()
        mb = sum(a.nbytes for L in (self.genes, self.transcripts, self.features)
                 for a in (L.codes, L.starts, L.ends, L.strands)) / 1e6
        return kv_html("Genes · 3 linked tables", [
            ("genes", f"{c['genes']:,} rows  (gene_id, gene_name, gene_type)"),
            ("transcripts", f"{c['transcripts']:,} rows  (transcript_id, gene → genes row)"),
            ("features", f"{c['exons']:,} exons · {c['CDS']:,} CDS · {c['UTR']:,} UTR "
                         f"(kind, transcript → transcripts row)"),
            ("coordinates", f"0-based half-open · {mb:.1f} MB of arrays"),
            ("annotation index", "built" if self._annot is not None else "built on first use"),
        ])


def _peaks(cre, genome):
    """Peaks for isoform selection: anything as_loci takes; several sets are merged."""
    if cre is None:
        return None
    from .interop import as_loci
    if isinstance(cre, (list, tuple)) and len(cre) > 1 and not (len(cre) in (3, 4) and isinstance(cre[1], int)):
        parts = [as_loci(p, genome=genome) for p in cre]
        out = parts[0]
        for p in parts[1:]:
            out = out + p
        return out.sort().merge()
    return as_loci(cre, genome=genome)
