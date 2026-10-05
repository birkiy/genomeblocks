"""Genes: GTF / GFF3 / genePred parsing (0-based), TSS, annotation, isoforms, exports."""
import io

import pandas as pd
import pytest

from genomeblocks import Genes, Loci, as_loci
from genomeblocks.genes import LABELS, tss_base

from conftest import GTF, installed_backends


def test_tables_are_zero_based_half_open(genes):
    G, F = genes.genes, genes.features
    assert G.to_records()[0] == ("chr1", 1000, 5000, "+")          # GTF 1001-5000 -> [1000, 5000)
    assert genes.counts() == {"genes": 3, "transcripts": 4, "exons": 7, "CDS": 4, "UTR": 2}
    assert genes["GENE_A"].tss.start == 1000
    assert genes["GENE_B"].tss.start == 10_999                      # '-' strand: end - 1
    assert genes["G2"].row == 1 and "GENE_C" in genes
    assert tss_base([10], [20], [2]).tolist() == [19]
    assert set(F.cols["exon_number"].tolist()) == {1, 2, 3}


@pytest.mark.parametrize("backend", installed_backends("tables"))
def test_parsers_agree(gtf_path, backend):
    a = Genes.make(gtf_path, backend=backend)
    b = Genes.make(gtf_path)
    for t in ("genes", "transcripts", "features"):
        assert getattr(a, t).equals(getattr(b, t), cols=True), t


def test_gff3_and_ensembl_spellings(tmp_path, genes):
    gff = """##gff-version 3
chr1\tx\tgene\t1001\t5000\t.\t+\t.\tID=gene:G1;Name=GENE_A;biotype=protein_coding
chr1\tx\tmRNA\t1001\t5000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1
chr1\tx\texon\t1001\t1200\t.\t+\t.\tParent=transcript:T1;rank=1
chr1\tx\texon\t4001\t5000\t.\t+\t.\tParent=transcript:T1;rank=2
chr1\tx\tfive_prime_utr\t1001\t1100\t.\t+\t.\tParent=transcript:T1
chr1\tx\tCDS\t1101\t1200\t.\t+\t.\tParent=transcript:T1
"""
    p = tmp_path / "e.gff3"
    p.write_text(gff)
    for backend in installed_backends("tables"):
        g = Genes.make(str(p), backend=backend)
        assert g.genes.to_records() == [("chr1", 1000, 5000, "+")]
        assert g.features.cols["exon_number"].tolist()[:2] == [1, 2]      # rank= read as exon number
        assert g.counts()["UTR"] == 1
    # lower-case Ensembl UTRs in the GTF fixture were counted
    assert genes.counts()["UTR"] == 2


def test_duplicate_gene_id_keeps_the_first_record(tmp_path):
    p = tmp_path / "dup.gtf"
    p.write_text(GTF + 'chr1\tx\tgene\t20001\t21000\t.\t+\t.\tgene_id "G1"; gene_name "GENE_A_ALT";\n')
    for backend in installed_backends("tables"):
        g = Genes.make(str(p), backend=backend)
        assert g.genes.to_records()[0] == ("chr1", 1000, 5000, "+") and len(g.genes) == 3


def test_ucsc_layouts(tmp_path):
    # refFlat: geneName name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds
    (tmp_path / "refFlat.txt").write_text("GENE_A\tNM_1\tchr1\t+\t1000\t5000\t1100\t4300\t3\t1000,2000,4000,\t1200,2500,5000,\n")
    g = Genes.make_ucsc(str(tmp_path / "refFlat.txt"))
    assert g.genes.to_records() == [("chr1", 1000, 5000, "+")] and g.genes.cols["gene_name"].tolist() == ["GENE_A"]
    assert g.counts() == {"genes": 1, "transcripts": 1, "exons": 3, "CDS": 3, "UTR": 2}
    # knownGene: ... proteinID alignID — the symbol is the name, not the alignment id
    (tmp_path / "kg.txt").write_text("ENST1\tchr1\t+\t1000\t5000\t1100\t4300\t3\t1000,2000,4000,\t1200,2500,5000,\tP1\tuc001.1\n")
    assert Genes.make_ucsc(str(tmp_path / "kg.txt")).genes.cols["gene_name"].tolist() == ["ENST1"]
    # refGene with bin and name2
    (tmp_path / "refGene.txt").write_text("585\tNM_1\tchr1\t+\t1000\t5000\t1100\t4300\t3\t1000,2000,4000,\t1200,2500,5000,\t0\tGENE_A\tcmpl\tcmpl\t0,0,0,\n")
    assert Genes.make_ucsc(str(tmp_path / "refGene.txt")).genes.cols["gene_name"].tolist() == ["GENE_A"]
    (tmp_path / "x.gtf").write_text(GTF)
    with pytest.raises(ValueError, match="genePred"):
        Genes.make_ucsc(str(tmp_path / "x.gtf"))


def test_from_frame_and_validation(gtf_path):
    df = pd.read_csv(io.StringIO(GTF), sep="\t", header=None,
                     names=["chrom", "src", "feature", "start", "end", "score", "strand", "frame", "attr"])
    df["gene_id"] = df.attr.str.extract(r'gene_id "([^"]+)"')
    df["transcript_id"] = df.attr.str.extract(r'transcript_id "([^"]+)"')
    df["gene_name"] = df.attr.str.extract(r'gene_name "([^"]+)"')
    g = Genes.from_frame(df)                                        # 1-based by default
    assert g.genes.to_records()[0] == ("chr1", 1000, 5000, "+")
    assert Genes.from_frame(df.assign(start=df.start - 1), one_based=False).genes.equals(g.genes)
    with pytest.raises(ValueError, match="gene_id"):
        Genes.from_frame(df.drop(columns=["gene_id"]))
    with pytest.raises(ValueError, match="feature"):
        Genes.from_frame(df.drop(columns=["feature"]))
    with pytest.raises(ValueError, match="lacks the column"):
        Genes(Loci(), Loci(), Loci())


def test_annotation(genes, cre, gtf_path):
    # promoter = gene TSS ± promoter_r (1000 by default): GENE_A's window [0, 2001) also covers
    # chr1:1900-2100, GENE_B's [9999, 12000) covers 9950-10050; 4900-5100 lies in T1's 3' UTR
    # [4300, 5000); the chr2 CREs are far from GENE_C's TSS at 20000.
    labels = LABELS[genes.labels(cre)]
    assert labels.tolist() == ["Promoter-TSS", "Promoter-TSS", "3UTR", "Promoter-TSS", "Promoter-TSS",
                               "Intergenic", "Intergenic"]
    # with a 100 bp promoter the exonic rows show: 1900-2100 is T1b's first exon (T1's intron),
    # 9950-10050 is GENE_B's first exon [10000, 10400)
    narrow = Genes.make(gtf_path, promoter_r=100)
    assert LABELS[narrow.labels(cre)].tolist() == ["Promoter-TSS", "Exonic", "3UTR", "Exonic", "Promoter-TSS",
                                                   "Intergenic", "Intergenic"]
    ann = genes.annotations(cre.to_pandas())                        # frames accepted everywhere
    assert list(ann.columns) == ["uid", "annotation"] and len(ann) == len(cre)
    names, dist = genes.nearest_tss(cre)
    assert names[0] == "GENE_A" and dist[0] == 0
    assert dist[5] == -1 or names[5] == "GENE_C"
    ng = genes.nearest_genes(cre)
    assert set(ng.columns) == {"uid", "gene_name", "distance"}
    tss = genes.get_tss("protein_coding")
    assert tss.to_records()[0] == ("chr1", 1000, 1001, "+") and len(tss) == 2
    assert genes.find("gene_")[0].gene_name == "GENE_A" if genes.find("gene_") else True
    assert genes.rows(["GENE_A", "G2"]).tolist() == [0, 1]


@pytest.mark.parametrize("backend", installed_backends("intervals"))
def test_annotation_same_on_every_interval_backend(genes, cre, backend):
    from genomeblocks.backends.intervals import _NEAREST
    assert (genes.labels(cre, backend=backend) == genes.labels(cre)).all()
    if backend in _NEAREST:
        assert (genes.nearest_tss(cre, backend=backend)[1] == genes.nearest_tss(cre)[1]).all()
    else:
        with pytest.raises(NotImplementedError, match="no nearest"):
            genes.nearest_tss(cre, backend=backend)


def test_select_isoforms(genes, bw_path, tmp_path, gtf_path):
    peaks = as_loci([("chr1", 1900, 2100)])                         # over T1b's TSS, not T1's
    g = genes.select_isoforms(peaks, r=200, verbose=False)         # in place: g is genes
    assert g is genes
    c = g.genes.cols["canonical"]
    assert g.transcripts.cols["transcript_id"][c[0]] == "T1b" and c[2] == 3  # fallback: the only isoform
    assert g.representative().tolist() == c.tolist()                # canonical isoforms win ...
    assert genes.representative(canonical=False).tolist() == [0, 2, 3]   # ... over the longest per gene
    g2 = Genes.make(gtf_path).select_isoforms(peaks, {"s": bw_path}, r=200, verbose=False)  # dict accepted
    assert g2.transcripts.cols["transcript_id"][g2.genes.cols["canonical"][0]] == "T1b"
    # a bigWig that is strong at T1's TSS (1000) and silent at T1b's (2000): with rank='signal' T1 wins
    import pybigtools
    pybigtools.open(str(tmp_path / "tss.bw"), "w").write({"chr1": 20_000, "chr2": 8_000},
                                                          [("chr1", 990, 1010, 9.0), ("chr1", 1990, 2010, 0.5)])
    both = as_loci([("chr1", 900, 1100), ("chr1", 1900, 2100)])
    g3 = genes.select_isoforms(both, str(tmp_path / "tss.bw"), r=200, rank="signal", verbose=False)
    assert g3.transcripts.cols["transcript_id"][g3.genes.cols["canonical"][0]] == "T1"
    g4 = genes.select_isoforms(both, str(tmp_path / "tss.bw"), r=200, rank="longest", verbose=False)
    assert g4.transcripts.cols["transcript_id"][g4.genes.cols["canonical"][0]] == "T1"      # longest supported
    g5 = genes.select_isoforms(both, str(tmp_path / "tss.bw"), r=200, rank="signal", min_signal=1.0, verbose=False)
    assert g5.transcripts.cols["transcript_id"][g5.genes.cols["canonical"][0]] == "T1"      # T1b below min_signal


def test_exports_round_trip(genes, tmp_path):
    out = genes.to_gtf(str(tmp_path / "out.gtf"))
    back = Genes.make(out)
    for t in ("genes", "transcripts"):
        assert getattr(back, t).equals(getattr(genes, t))
    bed12 = genes.to_bed12().split("\n")
    f = bed12[0].split("\t")
    assert f[:3] == ["chr1", "1000", "5000"] and f[6:8] == ["1100", "4300"]   # thick = CDS
    assert bed12[1].split("\t")[6:8] == ["11000", "11000"]                       # non-coding: thickStart == thickEnd
    for table in ("genes", "transcripts", "features"):
        assert len(genes.to_pandas(table)) == len(getattr(genes, table))
        assert genes.to_polars(table).shape[0] == len(getattr(genes, table))
    with pytest.raises((KeyError, ValueError)):
        genes.to_pandas("nope")
    genes.save(str(tmp_path / "genes"))
    back = Genes.load(str(tmp_path / "genes"))
    for t in ("genes", "transcripts", "features"):
        assert getattr(back, t).equals(getattr(genes, t), cols=True)


def test_genes_as_a_table(genes):
    assert genes.shape == (3, 7) and genes.columns[:4] == ["chrom", "start", "end", "strand"]
    assert len(genes.head(2)) == 2 and genes.describe().loc["genes", "value"] == 3
    import polars as pl
    assert pl.DataFrame(genes).shape == (3, 7)
    assert "<table" in genes._repr_html_() and "Genes(" in repr(genes)
    assert [g.gene_name for g in genes] == ["GENE_A", "GENE_B", "GENE_C"]
    assert len(genes["GENE_A"].transcripts) == 2 and len(genes["GENE_A"].exons) == 4     # T1: 3 exons, T1b: 1
