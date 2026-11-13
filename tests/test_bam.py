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
    total = every.signal(bw_path, span=True, n_bins=1, agg="sum", dtype=np.float64, progress=False, verbose=False)[:, 0, 0]
    assert np.allclose(every.cols["score"], np.nan_to_num(total))                         # score = summed signal, no width factor
    assert list(peaks.call_se([bw_path], stitch=1000, backend="python").cols) == ["n_peaks", "score", "rank"]
    cut, order = se.knee([10, 9, 8, 1, 0.5, 0.2, 0])
    assert cut <= 6 and len(order) == 6


def test_call_se_forwards_signal_keywords(bw_path):
    peaks = as_loci([("chr1", 1000, 1500), ("chr1", 1600, 2000), ("chr1", 2500, 3000), ("chr1", 10_000, 11_000),
                     ("chr1", 15_000, 15_500), ("chr2", 100, 200)])
    plain = se.call_se(peaks, bw_path, stitch=1000)
    shown = se.call_se(peaks, bw_path, stitch=1000, progress=True)
    assert shown.to_records() == plain.to_records()
    assert shown.cols["score"].tolist() == plain.cols["score"].tolist()
    assert se.call_se(peaks, bw_path, stitch=1000, verbose=True).to_records() == plain.to_records()


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
    frags = hichip.fragments(ends)
    runs = list(hichip.coverage(frags))
    assert [r[0] for r in runs] == ["chr1"]
    chrom, a, z, d = runs[0]
    depth = np.zeros(int(frags.ends.max()) + 1, int)
    for s, e in zip(frags.starts, frags.ends):
        depth[s:e] += 1
    got = np.zeros_like(depth)
    for s, e, v in zip(a, z, d):
        got[s:e] = v
    assert np.array_equal(got, depth)                                      # the runs are the per-base depth
    assert (a[1:] != z[:-1]).any() or (d[1:] != d[:-1]).all()               # adjacent equal runs merged
    out = hichip.shortrange_track(str(tmp_path / "t.allValidPairs"), str(tmp_path / "sr"), CHROM_SIZES)
    assert set(out) == {"ends", "bed", "bigwig"}
    with pytest.raises(RuntimeError, match="not found"):
        hichip.macs3("x.bed", "n", str(tmp_path / "out"), exe="macs3-nope")


# ── bam_to_bigwig ────────────────────────────────────────────────────────────

def _write_test_bam(path, pysam):
    """chr1 (2,000 bp) and chr2 (500 bp): forward / reverse reads, a spliced read, a deletion, a duplicate, a
    secondary, a low-MAPQ read, a proper pair and an unmapped read. Returns the reads as dicts."""
    hdr = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 2000, "SN": "chr1"}, {"LN": 500, "SN": "chr2"}]}
    M, D, N = 0, 2, 3
    reads = [dict(tid=0, pos=100, cigar=[(M, 30)], flag=0), dict(tid=0, pos=110, cigar=[(M, 30)], flag=16),
             dict(tid=0, pos=125, cigar=[(M, 30)], flag=0),
             dict(tid=0, pos=300, cigar=[(M, 20), (N, 100), (M, 30)], flag=0),          # spliced: 300-320, 420-450
             dict(tid=0, pos=600, cigar=[(M, 10), (D, 5), (M, 10)], flag=16),           # deletion: 600-610, 615-625
             dict(tid=0, pos=700, cigar=[(M, 30)], flag=1024),                          # duplicate
             dict(tid=0, pos=800, cigar=[(M, 30)], flag=256),                           # secondary
             dict(tid=0, pos=900, cigar=[(M, 30)], flag=0, mapq=5),                     # low MAPQ
             dict(tid=0, pos=1000, cigar=[(M, 30)], flag=1 | 2 | 32 | 64, mpos=1170, tlen=200),    # proper pair, read 1
             dict(tid=0, pos=1170, cigar=[(M, 30)], flag=1 | 2 | 16 | 128, mpos=1000, tlen=-200),  # its mate
             dict(tid=1, pos=50, cigar=[(M, 30)], flag=0), dict(tid=1, pos=480, cigar=[(M, 20)], flag=0),
             dict(tid=1, pos=60, cigar=[(M, 30)], flag=4)]                              # unmapped (placed)
    reads.sort(key=lambda r: (r["tid"], r["pos"]))
    with pysam.AlignmentFile(str(path), "wb", header=hdr) as f:
        for i, r in enumerate(reads):
            a = pysam.AlignedSegment()
            qlen = sum(n for op, n in r["cigar"] if op == M)
            a.query_name, a.query_sequence, a.flag = f"r{i}", "A" * qlen, r["flag"]
            a.reference_id, a.reference_start, a.cigartuples = r["tid"], r["pos"], r["cigar"]
            a.mapping_quality = r.get("mapq", 60)
            a.next_reference_id = r["tid"] if "mpos" in r else -1
            a.next_reference_start, a.template_length = r.get("mpos", -1), r.get("tlen", 0)
            a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
            f.write(a)
    pysam.index(str(path))
    return reads


def _blocks(r):
    out, x = [], r["pos"]
    for op, n in r["cigar"]:
        if op == 0:
            out.append((x, x + n))
        x += n
    return out


def _counted(reads, *, dup=False, mapq=0):
    bad = 0x4 | 0x100 | 0x200 | 0x800 | (0x400 if dup else 0)
    return [r for r in reads if not r["flag"] & bad and r.get("mapq", 60) >= mapq]


def _depth(intervals, size):
    d = np.zeros(size)
    for a, b in intervals:
        d[max(a, 0):min(b, size)] += 1
    return d


def _read_back(path, chrom, size):
    from genomeblocks.backends.bigwig import open_bigwig
    return np.nan_to_num(open_bigwig(str(path)).values(chrom, 0, size))


def test_bam_to_bigwig_per_base_depth_follows_the_aligned_blocks(tmp_path):
    pysam = pytest.importorskip("pysam")
    reads = _write_test_bam(tmp_path / "t.bam", pysam)
    out = gb.bam_to_bigwig(str(tmp_path / "t.bam"), str(tmp_path / "t.bw"))
    for tid, chrom, size in ((0, "chr1", 2000), (1, "chr2", 500)):
        want = _depth([b for r in _counted(reads) if r["tid"] == tid for b in _blocks(r)], size)
        assert np.array_equal(_read_back(out, chrom, size), want)
    got = _read_back(out, "chr1", 2000)
    assert got[320:420].sum() == 0 and got[610:615].sum() == 0          # the intron and the deletion hold nothing
    assert got[700:730].sum() == 30 and got[800:830].sum() == 0         # duplicates count unless ignored; secondary never
    dup = gb.bam_to_bigwig(str(tmp_path / "t.bam"), str(tmp_path / "d.bw"), ignore_duplicates=True, min_mapq=10)
    want = _depth([b for r in _counted(reads, dup=True, mapq=10) if r["tid"] == 0 for b in _blocks(r)], 2000)
    assert np.array_equal(_read_back(dup, "chr1", 2000), want)
    assert _read_back(dup, "chr1", 2000)[900:930].sum() == 0


def test_bam_to_bigwig_bins_count_each_read_once(tmp_path):
    pysam = pytest.importorskip("pysam")
    reads = _write_test_bam(tmp_path / "t.bam", pysam)
    out = gb.bam_to_bigwig(str(tmp_path / "t.bam"), str(tmp_path / "b.bw"), bin_size=50)
    for tid, chrom, size in ((0, "chr1", 2000), (1, "chr2", 500)):
        nb = (size + 49) // 50
        want = np.zeros(nb)
        for r in _counted(reads):
            if r["tid"] == tid:
                bins = {x // 50 for a, b in _blocks(r) for x in range(a, b)}
                for k in bins:
                    want[k] += 1
        got = _read_back(out, chrom, size)
        assert np.array_equal(got, np.repeat(want, 50)[:size])        # every base of a bin holds the bin's count


def test_bam_to_bigwig_extension_pairs_normalisation_and_workers(tmp_path):
    pysam = pytest.importorskip("pysam")
    reads = _write_test_bam(tmp_path / "t.bam", pysam)
    bam = str(tmp_path / "t.bam")
    ok = _counted(reads)
    # extend=60: every read to 60 bases along its strand (longer reads kept)
    frag = [(r["pos"], max(_blocks(r)[-1][1], r["pos"] + 60)) if not r["flag"] & 16
            else (min(r["pos"], _blocks(r)[-1][1] - 60), _blocks(r)[-1][1]) for r in ok if r["tid"] == 0]
    got = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "e.bw"), extend=60), "chr1", 2000)
    assert np.array_equal(got, _depth(frag, 2000))
    # pairs: the proper pair once over 1000-1200; other reads over their span
    pairs = [(1000, 1200)] + [(r["pos"], _blocks(r)[-1][1]) for r in ok if r["tid"] == 0 and not r["flag"] & 1]
    got = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "p.bw"), extend="pairs"), "chr1", 2000)
    assert np.array_equal(got, _depth(pairs, 2000))
    # CPM / RPKM / BPM / RPGC scale the raw counts by deepTools' factors
    raw = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "r.bw"), bin_size=50), "chr1", 2000)
    n = len(ok)
    for norm, factor in (("CPM", 1e6 / n), ("RPKM", 1e9 / (n * 50))):
        got = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / f"{norm}.bw"), bin_size=50, normalize=norm), "chr1", 2000)
        assert np.allclose(got, raw * factor, rtol=1e-6)
    bins_sum = sum(len({x // 50 for a, b in _blocks(r) for x in range(a, b)}) for r in ok)
    got = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "bpm.bw"), bin_size=50, normalize="BPM"), "chr1", 2000)
    assert np.allclose(got, raw * 1e6 / bins_sum, rtol=1e-6)
    got = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "g.bw"), normalize="RPGC", effective_genome_size=2500, extend=60),
                     "chr1", 2000)
    assert np.allclose(got, _depth(frag, 2000) * 2500 / (n * 60), rtol=1e-6)
    # processes give the same file; chroms= writes only those
    one = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "w1.bw")), "chr1", 2000)
    two = _read_back(gb.bam_to_bigwig(bam, str(tmp_path / "w2.bw"), workers=2), "chr1", 2000)
    assert np.array_equal(one, two)
    only = gb.bam_to_bigwig(bam, str(tmp_path / "c.bw"), chroms="chr2")
    assert _read_back(only, "chr1", 2000).sum() == 0 and _read_back(only, "chr2", 500).sum() > 0
    with pytest.raises(ValueError, match="effective_genome_size"):
        gb.bam_to_bigwig(bam, str(tmp_path / "x.bw"), normalize="RPGC")
    with pytest.raises(ValueError, match="normalize"):
        gb.bam_to_bigwig(bam, str(tmp_path / "x.bw"), normalize="TPM")
    with pytest.raises(ValueError, match="not in the BAM"):
        gb.bam_to_bigwig(bam, str(tmp_path / "x.bw"), chroms="chrX")
    import shutil
    shutil.copy(bam, tmp_path / "noindex.bam")
    with pytest.raises(ValueError, match="no index"):
        gb.bam_to_bigwig(str(tmp_path / "noindex.bam"), str(tmp_path / "x.bw"))
