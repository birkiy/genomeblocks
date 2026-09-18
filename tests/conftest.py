"""Shared fixtures for the genomeblocks test suite.

Everything here is synthetic — small in-memory objects and tiny temp files — so
the suite runs fast and needs no external genome / bigwig / ChIP-Atlas data.
"""
import textwrap

import pytest

from genomeblocks import Loci
from genomeblocks.locus import Locus


# ── tiny on-disk fixtures ────────────────────────────────────────────────────

@pytest.fixture
def bed_path(tmp_path):
    """A 4-interval BED (tab-separated, with a comment line)."""
    p = tmp_path / "peaks.bed"
    p.write_text(
        "# a comment\n"
        "chr1\t1000\t1200\tp1\t0\t+\n"
        "chr1\t1100\t1300\tp2\t0\t+\n"      # overlaps p1 (for merge tests)
        "chr1\t5000\t5200\tp3\t0\t-\n"
        "chr2\t2000\t2200\tp4\t0\t+\n"
    )
    return str(p)


@pytest.fixture
def bedpe_path(tmp_path):
    """A small BEDPE connecting anchors near the `cre` fixture loci."""
    p = tmp_path / "loops.bedpe"
    p.write_text(
        "chr1\t1000\t1100\tchr1\t5000\t5100\tl1\t10\n"
        "chr1\t1000\t1100\tchr1\t9000\t9100\tl2\t5\n"
        "chr1\t5000\t5100\tchr1\t9000\t9100\tl3\t3\n"
    )
    return str(p)


@pytest.fixture
def gtf_path(tmp_path):
    """A minimal but valid GTF: one + strand protein-coding gene with two exons."""
    p = tmp_path / "genes.gtf"
    p.write_text(textwrap.dedent('''\
        chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";
        chr1\tsrc\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1"; gene_type "protein_coding";
        chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
        chr1\tsrc\texon\t4800\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";
        chr1\tsrc\tCDS\t1050\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
    '''))
    return str(p)


@pytest.fixture
def gtf_isoforms_path(tmp_path):
    """A GTF with two multi-isoform genes, one per strand.

    G1 (+, TSS = start): T_long 1000-9000, T_mid 4000-9000, T_short 6000-9000.
    G2 (-, TSS = end):   T_far  20000-30000, T_near 20000-26000.
    Only the non-longest isoforms sit under the `atac_bed_path` peaks.
    """
    p = tmp_path / "isoforms.gtf"
    rows = [
        ("gene",       "G1", None,      1000,  9000, "+"),
        ("transcript", "G1", "T_long",  1000,  9000, "+"),
        ("exon",       "G1", "T_long",  1000,  9000, "+"),
        ("transcript", "G1", "T_mid",   4000,  9000, "+"),
        ("exon",       "G1", "T_mid",   4000,  9000, "+"),
        ("transcript", "G1", "T_short", 6000,  9000, "+"),
        ("exon",       "G1", "T_short", 6000,  9000, "+"),
        ("gene",       "G2", None,     20000, 30000, "-"),
        ("transcript", "G2", "T_far",  20000, 30000, "-"),
        ("exon",       "G2", "T_far",  20000, 30000, "-"),
        ("transcript", "G2", "T_near", 20000, 26000, "-"),
        ("exon",       "G2", "T_near", 20000, 26000, "-"),
    ]
    lines = []
    for feat, gid, tid, start, end, strand in rows:
        attrs = f'gene_id "{gid}"; gene_name "{gid}"; gene_type "protein_coding";'
        if tid: attrs += f' transcript_id "{tid}"; exon_number "1";'
        lines.append(f"chr1\tsrc\t{feat}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


@pytest.fixture
def atac_bed_path(tmp_path):
    """Peaks over the TSSs of G1's T_mid / T_short and G2's T_near."""
    p = tmp_path / "atac.bed"
    p.write_text(
        "chr1\t3900\t4100\ta1\t0\t.\n"
        "chr1\t5950\t6050\ta2\t0\t.\n"
        "chr1\t25900\t26100\ta3\t0\t.\n"
    )
    return str(p)


@pytest.fixture
def atac_bw_path(tmp_path):
    """A tiny bigwig: weak at T_mid's TSS, strong at T_short's and T_near's."""
    pybigtools = pytest.importorskip("pybigtools")
    p = tmp_path / "atac.bw"
    pybigtools.open(str(p), "w").write(
        {"chr1": 40000},
        [("chr1", 3990, 4010, 1.0),
         ("chr1", 5990, 6010, 4.0),
         ("chr1", 25990, 26010, 4.0)],
    )
    return str(p)


# ── in-memory objects ────────────────────────────────────────────────────────

@pytest.fixture
def cre():
    """Five non-overlapping loci on chr1 at 1k, 5k, 9k, 13k, 17k."""
    return Loci([Locus("chr1", b, b + 100) for b in (1000, 5000, 9000, 13000, 17000)])


@pytest.fixture
def arch(cre):
    """A small Architecture over `cre` with ep.w / ep.n set.

    Vertex 0 (chr1:1000-1100) is a hub: connected to 1, 2, 3 with strong weights;
    a weak 1-2 edge is the background.
    """
    from genomeblocks import Architecture
    g = Architecture(name="test")
    uids = [l.uid for l in cre]
    edges = [(0, 1, 5.0), (0, 2, 5.0), (0, 3, 5.0), (1, 2, 0.1)]
    for a, b, w in edges:
        e = g.add_edge(g._add_vertex(uids[a]), g._add_vertex(uids[b]))
        g.ep.w[e] = w
        g.ep.n[e] = w
    return g
