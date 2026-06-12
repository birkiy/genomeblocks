from genomeblocks import Genes, Loci
from genomeblocks.locus import Locus


def test_make_parses_gene_and_transcript(gtf_path):
    genes = Genes.make(gtf_path, promoter_r=1000)
    assert len(genes) == 1
    g = genes["G1"]
    assert g.gene_name == "GENE1"
    assert g.gene_type == "protein_coding"
    assert "T1" in g.transcripts
    assert len(g.transcripts["T1"].exons) == 2


def test_tss_on_plus_strand(gtf_path):
    genes = Genes.make(gtf_path, promoter_r=1000)
    tss = genes.get_tss()
    assert tss["GENE1"].start == 1000           # + strand TSS = gene start


def test_annot_building_blocks(gtf_path):
    genes = Genes.make(gtf_path, promoter_r=1000)
    annot = genes.annot
    for key in ("body", "prom", "exon"):
        assert key in annot and len(annot[key]) >= 1


def test_annotations_region_classes(gtf_path):
    genes = Genes.make(gtf_path, promoter_r=1000)
    loci = Loci([
        Locus("chr1", 1050, 1150),    # in promoter window (TSS +-1000) -> Promoter-TSS
        Locus("chr1", 2500, 2600),    # in gene body, not exon/prom      -> Intronic
        Locus("chr1", 4850, 4950),    # exon 2                            -> Exonic
        Locus("chr1", 20000, 20100),  # far away                          -> Intergenic
    ])
    df = genes.annotations(loci)
    by_uid = dict(zip(df["uid"], df["annotation"]))
    assert by_uid[loci[0].uid] == "Promoter-TSS"
    assert by_uid[loci[1].uid] == "Intronic"
    assert by_uid[loci[2].uid] == "Exonic"
    assert by_uid[loci[3].uid] == "Intergenic"


def test_nearest_genes(gtf_path):
    genes = Genes.make(gtf_path, promoter_r=1000)
    loci = Loci([Locus("chr1", 20000, 20100)])
    df = genes.nearest_genes(loci)
    assert df.iloc[0]["Name_b"] == "GENE1"
    assert df.iloc[0]["Distance"] > 0
