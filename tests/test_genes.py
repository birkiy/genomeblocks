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


# ── ATAC-supported isoform selection ─────────────────────────────────────────

def test_select_isoforms_peaks_pick_longest_supported(gtf_isoforms_path, atac_bed_path):
    genes = Genes.make(gtf_isoforms_path, cre=atac_bed_path,
                       kw={"verbose": False})
    g1 = genes["G1"]
    assert g1.canonical == "T_mid"                 # longest isoform with an open TSS
    assert (g1.start, g1.end) == (4000, 9000)      # body collapsed onto it
    assert g1.tss.start == 4000
    assert [t.tss_support for t in g1.transcripts.values()] == [False, True, True]
    assert len(g1.transcripts) == 3                # nothing dropped


def test_select_isoforms_is_strand_aware(gtf_isoforms_path, atac_bed_path):
    genes = Genes.make(gtf_isoforms_path, cre=atac_bed_path,
                       kw={"verbose": False})
    g2 = genes["G2"]                               # '-' strand: TSS is the end
    assert g2.canonical == "T_near"
    assert (g2.start, g2.end) == (20000, 26000)
    assert g2.tss.start == 26000


def test_select_isoforms_fallback_when_nothing_open(gtf_isoforms_path, tmp_path):
    empty = tmp_path / "none.bed"
    empty.write_text("chr9\t100\t200\tx\t0\t.\n")
    genes = Genes.make(gtf_isoforms_path, cre=str(empty),
                       kw={"verbose": False})
    g1 = genes["G1"]
    assert g1.canonical == "T_long"                # longest annotated isoform
    assert (g1.start, g1.end) == (1000, 9000)
    assert not any(t.tss_support for t in g1.transcripts.values())


def test_select_isoforms_collapse_false_keeps_gene_span(gtf_isoforms_path, atac_bed_path):
    genes = Genes.make(gtf_isoforms_path)
    genes.select_isoforms(atac_bed_path, collapse=False, verbose=False)
    g1 = genes["G1"]
    assert g1.canonical == "T_mid"
    assert (g1.start, g1.end) == (1000, 9000)


def test_select_isoforms_signal_relative_cut(gtf_isoforms_path, atac_bw_path):
    # T_short's TSS is 4x T_mid's: at min_frac=0.5 only T_short survives.
    genes = Genes.make(gtf_isoforms_path, bw=atac_bw_path,
                       kw={"verbose": False})
    g1 = genes["G1"]
    assert g1.canonical == "T_short"
    assert g1.transcripts["T_short"].tss_score == 4.0
    assert g1.transcripts["T_mid"].tss_score == 1.0
    assert g1.transcripts["T_long"].tss_score == 0.0

    # a looser cut lets T_mid back in, and 'longest' then prefers it
    genes = Genes.make(gtf_isoforms_path, bw=atac_bw_path,
                       kw={"min_frac": 0.2, "verbose": False})
    assert genes["G1"].canonical == "T_mid"

    # ...unless the ranking asks for the strongest TSS instead
    genes = Genes.make(gtf_isoforms_path, bw=atac_bw_path,
                       kw={"min_frac": 0.2, "rank": "signal", "verbose": False})
    assert genes["G1"].canonical == "T_short"


def test_select_isoforms_bed_and_signal_intersect(gtf_isoforms_path, atac_bed_path,
                                                  atac_bw_path):
    # peaks gate first (T_long is out), then the signal cut (T_mid is out)
    genes = Genes.make(gtf_isoforms_path, cre=atac_bed_path,
                       bw=atac_bw_path, kw={"verbose": False})
    g1 = genes["G1"]
    assert g1.canonical == "T_short"
    assert g1.transcripts["T_long"].tss_score is None      # never queried
    assert g1.transcripts["T_mid"].tss_support is False


def test_select_isoforms_needs_evidence(gtf_isoforms_path):
    import pytest
    genes = Genes.make(gtf_isoforms_path)
    with pytest.raises(ValueError):
        genes.select_isoforms()


def test_select_isoforms_rebuilds_annot(gtf_isoforms_path, atac_bed_path):
    genes = Genes.make(gtf_isoforms_path)
    assert any(l.start == 1000 for l in genes.annot["body"])
    genes.select_isoforms(atac_bed_path, verbose=False)
    assert all(l.start != 1000 for l in genes.annot["body"])   # body follows T_mid


def test_browser_shows_canonical_first(gtf_isoforms_path, atac_bed_path):
    from genomeblocks.browserview import _select_transcripts
    genes = Genes.make(gtf_isoforms_path, cre=atac_bed_path,
                       kw={"verbose": False})
    picked = _select_transcripts(genes["G1"], 1)
    assert [t.transcript_id for t in picked] == ["T_mid"]


def test_select_isoforms_takes_loci_and_promoter_r_window(gtf_isoforms_path):
    # a Loci works wherever a BED path does; the TSS window defaults to promoter_r
    peak = Loci([Locus("chr1", 6800, 6850)])       # 800 bp from T_short's TSS
    genes = Genes.make(gtf_isoforms_path, promoter_r=1000, cre=peak,
                       kw={"verbose": False})
    assert genes["G1"].canonical == "T_short"

    genes = Genes.make(gtf_isoforms_path, promoter_r=500, cre=peak,
                       kw={"verbose": False})
    assert genes["G1"].canonical == "T_long"       # out of reach: fell back

    genes = Genes.make(gtf_isoforms_path, promoter_r=500, cre=peak, r=1000,
                       kw={"verbose": False})
    assert genes["G1"].canonical == "T_short"      # explicit r wins
