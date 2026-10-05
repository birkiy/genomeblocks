"""BAM pileups, super-enhancers and HiChIP short-range tracks."""
import numpy as np
import pytest

import genomeblocks as gb
from genomeblocks import Genome, as_loci, hichip, se

from conftest import CHROM_SIZES


def test_bam_pileup_needs_pysam_or_counts(tmp_path):
    pysam = pytest.importorskip("pysam")
    from genomeblocks.bam import pileup_counts, reference_seq
    hdr = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 2000, "SN": "1"}]}       # no 'chr' prefix
    with pysam.AlignmentFile(str(tmp_path / "t.bam"), "wb", header=hdr) as f:
        for i in range(5):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag, a.reference_id = f"r{i}", "A" * 20, 0, 0
            a.reference_start, a.mapping_quality, a.cigar = 100 + i, 60, [(0, 20)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 20)
            f.write(a)
    pysam.index(str(tmp_path / "t.bam"))
    c = pileup_counts(str(tmp_path / "t.bam"), "chr1", 90, 130)       # chr1 -> 1 resolved
    assert c.shape == (4, 40) and c[0].max() == 5 and c[1:].sum() == 0
    assert gb.coverage(str(tmp_path / "t.bam"), ("chr1", 90, 130)).sum() == 100
    (tmp_path / "r.fa").write_text(">1\n" + "C" * 2000 + "\n")
    pysam.faidx(str(tmp_path / "r.fa"))
    assert reference_seq(str(tmp_path / "r.fa"), "chr1", -5, 10) == "N" * 5 + "C" * 10


def test_stitch_peaks_from_peak_boundaries():
    peaks = as_loci([("chr1", 1000, 1500), ("chr1", 1600, 2000), ("chr1", 2500, 3000), ("chr1", 10_000, 11_000),
                     ("chr1", 40_000, 40_500), ("chr2", 100, 200)])
    st = se.stitch_peaks(peaks)
    assert st.to_records() == [("chr1", 1000, 11_000, "."), ("chr1", 40_000, 40_500, "."), ("chr2", 100, 200, ".")]
    assert st.cols["n_peaks"].tolist() == [4, 1, 1]
    assert se.stitch_peaks(peaks.to_pandas(), stitch=1000).to_records()[0] == ("chr1", 1000, 3000, ".")
    assert len(se.stitch_peaks(as_loci([]))) == 0


def test_call_se(bw_path):
    peaks = as_loci([("chr1", 1000, 1500), ("chr1", 1600, 2000), ("chr1", 2500, 3000), ("chr1", 10_000, 11_000),
                     ("chr1", 15_000, 15_500), ("chr2", 100, 200)])
    ses, every = se.call_se(peaks.to_pandas(), bw_path, stitch=1000, return_all=True)
    assert list(every.cols) == ["n_peaks", "score", "rank"] and len(every) == 4
    pos = every.cols["score"] > 0
    assert every.cols["rank"][pos].min() == 1 and (every.cols["rank"][~pos] == 0).all()    # zero scores are unranked
    assert every.cols["score"][every.cols["rank"] == 1][0] == every.cols["score"].max()
    assert 0 < len(ses) <= len(every) and ses.is_sorted
    assert list(peaks.call_se([bw_path], stitch=1000, backend="python").cols) == ["n_peaks", "score", "rank"]
    cut, order = se.knee([10, 9, 8, 1, 0.5, 0.2, 0])
    assert cut <= 6 and len(order) == 6


def test_nearest_gene_within_recodes_genomes():
    regions = as_loci([("chr3", 5, 10), ("chr2", 0, 20)])
    tss = as_loci([("chr2", 5, 6)])
    assert se.nearest_gene_within(regions, tss, 100).tolist() == [-1, 0]


def test_hichip_short_range_track(tmp_path, sizes_path):
    lines = "".join(f"r{i}\tchr1\t{100 + i * 10}\t+\tchr1\t{600 + i * 10}\t-\t500\n" for i in range(20))
    lines += "r99\tchr1\t100\t+\tchr2\t5000\t-\t0\nr100\tchr1\t100\t+\tchr1\t9000\t-\t8900\n"
    (tmp_path / "t.allValidPairs").write_text(lines)
    ends = hichip.shortrange_ends(str(tmp_path / "t.allValidPairs"), 1000)
    assert len(ends) == 40 and ends.is_sorted and set(ends.strand) == {"+", "-"}
    ends_pd = hichip._shortrange_ends_pandas(str(tmp_path / "t.allValidPairs"), 1000, hichip._AVP, Genome())
    assert ends_pd.equals(ends)
    frags = hichip.fragments(ends.to_pandas(), 147, Genome.from_sizes({"chr1": 300}))
    assert frags.ends.max() == 300 and frags.lengths.max() <= 147
    bw = hichip.to_bigwig(hichip.fragments(ends), str(tmp_path / "cov.bw"), sizes_path)
    from genomeblocks.backends.bigwig import open_bigwig
    h = open_bigwig(bw)
    assert h.stats_array("chr1", 100, 200, n_bins=1, stat="max")[0] > 0
    h.close()
    runs = list(hichip.coverage(hichip.fragments(ends)))
    assert runs[0][0] == "chr1" and (runs[0][3] > 0).all()
    out = hichip.shortrange_track(str(tmp_path / "t.allValidPairs"), str(tmp_path / "sr"), CHROM_SIZES)
    assert set(out) == {"ends", "bed", "bigwig"}
    with pytest.raises(RuntimeError, match="not found"):
        hichip.macs3("x.bed", "n", str(tmp_path / "out"), exe="macs3-nope")
