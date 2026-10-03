"""Columnar Genes: three linked tables instead of a tree of objects.

    genes        one row per gene        (gene_id, gene_name, gene_type)
    transcripts  one row per transcript  (transcript_id, gene -> genes row)
    features     one row per exon/CDS/UTR (kind, transcript -> transcripts row,
                                           exon_number)

Each table is a columnar :class:`Loci` sharing the session Genome, so gene
coordinates and CRE coordinates are directly comparable. The "links" are
integer row numbers (``transcripts.cols['gene'][k]`` is the gene row of
transcript ``k``), not nested objects.

The annotation index (promoter / exon / UTR / gene-body intervals) is built
once from the tables, cached, and reused by every ``annotations`` call.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..locus import Locus
from .genome import Genome, default_genome
from .loci import STRANDS, Loci, SCODE

KINDS = ("exon", "CDS", "5UTR", "3UTR")
_FEATURE = {"exon": 0, "CDS": 1, "five_prime_UTR": 2, "three_prime_UTR": 3, "UTR": 4}
LABELS = np.array(["Intergenic", "Intronic", "Exonic", "3UTR", "5UTR", "Promoter-TSS"], dtype=object)


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


def _tables_polars(path, g, name_k, type_k):
    """Build the three tables with polars joins; only numbers cross into numpy
    (turning 900k strings into Python objects would cost more than the parse)."""
    import polars as pl
    df = _gtf_frame(path, ("gene_id", "transcript_id", name_k, type_k, "exon_number"))
    chroms = df["chrom"].unique(maintain_order=True).to_list()
    df = df.with_columns(
        pl.col("chrom").replace_strict(chroms, [g._add(c) for c in chroms],
                                       return_dtype=pl.Int32).alias("code"),
        pl.col("strand").replace_strict([".", "+", "-"], [0, 1, 2], default=0,
                                        return_dtype=pl.Int8).alias("s"),
        pl.col("feature").replace_strict(list(_FEATURE), list(_FEATURE.values()), default=-1,
                                         return_dtype=pl.Int8).alias("kind"))
    genes = (df.filter(pl.col("feature") == "gene").unique("gene_id", keep="last").sort("line")
             .with_columns(pl.col(name_k).fill_null(pl.col("gene_id"))))
    tx = df.filter(pl.col("feature") == "transcript")
    orphan = (tx.join(genes.select("gene_id"), on="gene_id", how="anti")
              .unique("gene_id", keep="first").sort("line")
              .with_columns(pl.lit("").alias(name_k), pl.lit("").alias(type_k)))
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
                   "gene_type": genes[type_k].to_numpy().astype(object)})
    T = Loci(num(tx, "code", np.int32), num(tx, "start", np.int64), num(tx, "end", np.int64),
             num(tx, "s", np.int8), genome=g,
             cols={"transcript_id": tx["transcript_id"].to_numpy().astype(object),
                   "gene": num(tx, "gene_row", np.int32)})
    F = Loci(num(ft, "code", np.int32), num(ft, "start", np.int64), num(ft, "end", np.int64),
             num(ft, "s", np.int8), genome=g,
             cols={"kind": num(ft, "kind", np.int8), "transcript": num(ft, "tx_row", np.int32),
                   "exon_number": num(ft, "exon_number", np.int32)})
    return G, T, F


def _read_gtf_pandas(path, keys):
    """GTF -> dict of numpy columns with pandas (fallback when polars is missing)."""
    import pandas as pd
    cols = ["chrom", "feature", "start", "end", "strand", "attr"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                     names=cols, dtype={"chrom": str})
    out = {c: df[c].to_numpy() for c in cols if c != "attr"}
    for k in keys:
        out[k] = df["attr"].str.extract(rf'(?:^|;\s*){k} "?([^";]*)"?', expand=False).to_numpy(object)
    return out


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


class GeneView(Locus):
    """One gene as a Locus, with its transcripts/exons one attribute away."""

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
    def tss(self) -> Locus:
        return Locus(self.chrom, self.start, self.start + 1) if self.strand == "+" else \
            Locus(self.chrom, self.end, self.end - 1)

    @property
    def transcripts(self) -> Loci:
        T = self._G.transcripts
        return T.take(np.flatnonzero(T.cols["gene"] == self._i))

    @property
    def exons(self) -> Loci:
        F, T = self._G.features, self._G.transcripts
        tx = np.flatnonzero(T.cols["gene"] == self._i)
        return F.take(np.flatnonzero(np.isin(F.cols["transcript"], tx) & (F.cols["kind"] == 0)))

    def __repr__(self):
        return f"Gene({self.gene_name}, {self.chrom}:{self.start:,}-{self.end:,}({self.strand}), {self.gene_type})"


class Genes:
    """Genes, transcripts and features as linked columnar tables."""

    def __init__(self, genes: Loci, transcripts: Loci, features: Loci, *, promoter_r: int = 1000,
                 filename: Optional[str] = None):
        self.genes, self.transcripts, self.features = genes, transcripts, features
        self.promoter_r = promoter_r
        self.filename = filename
        self._annot = None
        self._by_key = None

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, filename: str, *, gene_name_key: str = "gene_name", gene_type_key: str = "gene_type",
             promoter_r: int = 1000, genome: Optional[Genome] = None) -> "Genes":
        """Parse a GTF/GFF into the three tables (polars if installed, else pandas)."""
        g = genome or default_genome()
        try:
            G, T, F = _tables_polars(filename, g, gene_name_key, gene_type_key)
        except ImportError:
            G, T, F = _tables_numpy(_read_gtf_pandas(
                filename, ("gene_id", "transcript_id", gene_name_key, gene_type_key, "exon_number")),
                g, gene_name_key, gene_type_key)
        _resolve_generic_utr(F, T)
        return cls(G, T, F, promoter_r=promoter_r, filename=filename)

    # ── lookups ───────────────────────────────────────────────────────────
    def __len__(self):
        return len(self.genes)

    def _key_index(self):
        if self._by_key is None:
            ids = self.genes.cols["gene_id"].tolist()
            names = self.genes.cols["gene_name"].tolist()
            by = {n: i for i, n in enumerate(names)}
            by.update({g: i for i, g in enumerate(ids)})
            self._by_key = by
        return self._by_key

    def __getitem__(self, key) -> GeneView:
        """``genes['TP53']`` or ``genes['ENSG...']`` or ``genes[row]``."""
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

    # ── TSS and the annotation index ─────────────────────────────────────
    def _tss_rows(self, gene_type=None):
        """One gene row per gene name (the last, like the classic dict)."""
        G = self.genes
        rows = np.arange(len(G))
        if gene_type is not None:
            types = [gene_type] if isinstance(gene_type, str) else list(gene_type)
            rows = rows[np.isin(G.cols["gene_type"].astype(str), types)]
        names = G.cols["gene_name"][rows]
        _, last = np.unique(names[::-1], return_index=True)
        return np.sort(rows[::-1][last])

    def get_tss(self, gene_type=None) -> Loci:
        """1-bp TSS per gene name ('-' keeps the classic ``(end, end-1)`` orientation)."""
        G = self.genes
        r = self._tss_rows(gene_type)
        plus = G.strands[r] != 2
        s = np.where(plus, G.starts[r], G.ends[r])
        e = np.where(plus, G.starts[r] + 1, G.ends[r] - 1)
        return Loci(G.codes[r], s, e, np.zeros(len(r), np.int8), genome=G.genome,
                    cols={"gene_name": G.cols["gene_name"][r]})

    @property
    def annot(self) -> dict:
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

    def labels(self, L: Loci) -> np.ndarray:
        """Region label code per row of ``L`` (index into ``LABELS``)."""
        a = self.annot
        out = np.zeros(len(L), np.int8)
        for code, key in ((1, "body"), (2, "exon"), (3, "utr3"), (4, "utr5"), (5, "prom")):
            out[L.overlap_any(a[key])] = code               # later = higher priority
        return out

    def annotations(self, L: Loci):
        """Classic layout: DataFrame(uid, annotation)."""
        import pandas as pd
        return pd.DataFrame({"uid": L.uid, "annotation": LABELS[self.labels(L)]})

    def nearest_tss(self, L: Loci):
        """(gene_name, distance) of the nearest TSS window (±promoter_r) per row."""
        win = self.get_tss().slop(self.promoter_r)
        j, dist = L.nearest(win)
        names = np.where(j >= 0, win.cols["gene_name"][np.maximum(j, 0)], "")
        return names.astype(object), dist

    def nearest_genes(self, L: Loci):
        """Classic layout: DataFrame(Name=uid, Name_b=gene, Distance)."""
        import pandas as pd
        names, dist = self.nearest_tss(L)
        ok = dist >= 0
        return pd.DataFrame({"Name": L.uid[ok], "Name_b": names[ok], "Distance": dist[ok]})

    # ── display ───────────────────────────────────────────────────────────
    def counts(self) -> dict:
        k = self.features.cols["kind"]
        return {"genes": len(self.genes), "transcripts": len(self.transcripts),
                "exons": int((k == 0).sum()), "CDS": int((k == 1).sum()),
                "UTR": int((k >= 2).sum())}

    def table(self) -> str:
        c = self.counts()
        lines = [f"{'Name':<20} {'Count':<10}", "-" * 30]
        lines += [f"{'Transcripts':<20} {c['transcripts']:<10}", f"{'Exons':<20} {c['exons']:<10}",
                  f"{'CDS':<20} {c['CDS']:<10}", f"{'UTR':<20} {c['UTR']:<10}"]
        return "\n".join(lines)

    def __repr__(self):
        c = self.counts()
        return (f"Genes({c['genes']:,} genes, {c['transcripts']:,} transcripts, "
                f"{c['exons']:,} exons, columnar)")

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
            ("coordinates", f"{mb:.1f} MB of arrays"),
            ("annotation index", "built" if self._annot is not None else "built on first use"),
        ])


def _tables_numpy(d, g, gene_name_key, gene_type_key):
    """Same tables from numpy/pandas columns (no polars)."""
    import pandas as pd
    feat, gid, tid = d["feature"], d["gene_id"].astype(object), d["transcript_id"].astype(object)
    codes = g.encode(d["chrom"])
    strands = pd.Series(d["strand"]).map(SCODE).fillna(0).to_numpy(np.int8)
    starts, ends = d["start"].astype(np.int64), d["end"].astype(np.int64)

    # genes: the 'gene' lines (last line wins for a repeated gene_id, like a dict)
    gi = np.flatnonzero(feat == "gene")
    _, keep = np.unique(gid[gi][::-1], return_index=True)
    gi = np.sort(gi[::-1][keep])
    ti = np.flatnonzero(feat == "transcript")
    # genes that only appear through their transcripts get a row from the first one
    known = pd.Index(gid[gi])
    orphan = ti[known.get_indexer(gid[ti]) < 0]
    if len(orphan):
        _, first = np.unique(gid[orphan], return_index=True)
        gi_all = np.concatenate([gi, orphan[np.sort(first)]])
    else:
        gi_all = gi
    n_named = len(gi)
    name = np.asarray(d[gene_name_key], dtype=object)[gi_all]
    gtype = np.asarray(d[gene_type_key], dtype=object)[gi_all]
    nan_name = pd.isna(name)
    name[nan_name] = gid[gi_all][nan_name]
    name[n_named:] = ""                        # classic Gene() default for orphans
    gtype[pd.isna(gtype)] = None
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
    k = kind[fi].copy()
    F = Loci(codes[fi], starts[fi], ends[fi], strands[fi], genome=g,
             cols={"kind": k, "transcript": tx_row, "exon_number": en})
    return G, T, F
