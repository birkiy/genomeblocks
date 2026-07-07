"""Synthetic-fixture tests for the BAM coverage / mismatch viewer.

Builds a tiny reference FASTA and a coordinate-sorted+indexed BAM on the fly
(no external data), then checks pileup extraction and the browser drawer.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from matplotlib.colors import to_rgba

from genomeblocks.browserview import _detect_track_type, browser
from genomeblocks import bam as bammod
from genomeblocks import coverage


# ── fixtures ─────────────────────────────────────────────────────────────────

REF = ("ACGT" * 10)          # 40 bp on chr1
N_MATCH, N_VAR = 6, 4        # reads at [10,30); variants carry 'T' at pos 20


@pytest.fixture
def fasta_path(tmp_path):
    p = tmp_path / "ref.fa"
    p.write_text(f">chr1\n{REF}\n")
    pysam.faidx(str(p))
    return str(p)


@pytest.fixture
def bam_path(tmp_path):
    unsorted = tmp_path / "reads.unsorted.bam"
    header = {"HD": {"VN": "1.0"},
              "SQ": [{"SN": "chr1", "LN": len(REF)}]}
    read = REF[10:30]                       # perfect-match 20-mer
    var = read[:10] + "T" + read[11:]       # mismatch at genomic pos 20 (ref 'A')
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for i in range(N_MATCH + N_VAR):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = read if i < N_MATCH else var
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 10
            a.mapping_quality = 60
            a.cigartuples = [(0, 20)]        # 20M
            a.query_qualities = pysam.qualitystring_to_array("I" * 20)
            out.write(a)
    sorted_bam = tmp_path / "reads.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


# ── detection ────────────────────────────────────────────────────────────────

def test_detect_bam():
    assert _detect_track_type("x.bam") == "bam"


# ── extraction ───────────────────────────────────────────────────────────────

def test_pileup_counts_shape_and_total(bam_path):
    counts = bammod.pileup_counts(bam_path, "chr1", 10, 30)
    assert counts.shape == (4, 20)
    total = counts.sum(axis=0)
    assert np.all(total == N_MATCH + N_VAR)          # every position fully covered


def test_pileup_variant_allele(bam_path):
    counts = bammod.pileup_counts(bam_path, "chr1", 10, 30)
    pos = 20 - 10                                     # column for genomic 20
    a_idx, t_idx = 0, 3                               # A, T rows
    assert counts[a_idx, pos] == N_MATCH             # reference 'A'
    assert counts[t_idx, pos] == N_VAR               # variant 'T'


def test_coverage_helper_matches(bam_path):
    cov = coverage(bam_path, ("chr1", 10, 30))
    assert cov.tolist() == [N_MATCH + N_VAR] * 20


def test_pileup_pads_contig_edge(bam_path):
    # Window runs past the 40 bp contig; the tail must be zero-padded, not error.
    counts = bammod.pileup_counts(bam_path, "chr1", 30, 50)
    assert counts.shape == (4, 20)
    assert counts[:, 10:].sum() == 0                 # positions 40-50 don't exist


def test_reference_seq_uppercase_and_pad(fasta_path):
    assert bammod.reference_seq(fasta_path, "chr1", 10, 30) == REF[10:30]
    padded = bammod.reference_seq(fasta_path, "chr1", 35, 45)   # runs off contig
    assert padded == REF[35:40] + "N" * 5


def test_contig_name_resolution(bam_path):
    # A '1' request against a 'chr1' BAM should still resolve.
    counts = bammod.pileup_counts(bam_path, "1", 10, 30)
    assert counts.sum(axis=0).tolist() == [N_MATCH + N_VAR] * 20


# ── browser drawing ──────────────────────────────────────────────────────────

def test_browser_requires_reference_for_bam(bam_path):
    with pytest.raises(ValueError, match="reference"):
        browser(("chr1", 10, 30), {"aln": bam_path})


def test_browser_bam_draws_mismatch(fasta_path, bam_path):
    fig, ax = browser(("chr1", 10, 30), {"aln": bam_path}, reference=fasta_path)
    # The mismatch bars are a red (T) PatchCollection layered over the gray fill.
    red = to_rgba("#ff0000")
    reds = [c for c in ax["aln"].collections
            if len(c.get_facecolor()) and np.allclose(c.get_facecolor()[0], red)]
    assert reds, "expected a red mismatch collection at the variant position"


def test_browser_low_freq_variant_stays_gray(fasta_path, bam_path):
    # allele_freq above the 0.4 observed fraction => no mismatch coloring.
    fig, ax = browser(("chr1", 10, 30), {"aln": bam_path}, reference=fasta_path,
                      bam_allele_freq=0.9)
    red = to_rgba("#ff0000")
    reds = [c for c in ax["aln"].collections
            if len(c.get_facecolor()) and np.allclose(c.get_facecolor()[0], red)]
    assert not reds


def test_browser_sequence_track(fasta_path, bam_path):
    fig, ax = browser(("chr1", 10, 30), {"aln": bam_path}, reference=fasta_path)
    assert "_sequence" in ax
    assert len(ax["_sequence"].texts) == 20          # one colored letter per base


def test_browser_bam_ymax(fasta_path, bam_path):
    fig, ax = browser(("chr1", 10, 30), {"aln": bam_path}, reference=fasta_path,
                      bam_ymax=25.0)
    assert ax["aln"].get_ylim()[1] == pytest.approx(25.0)
