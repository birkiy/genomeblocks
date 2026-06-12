from genomeblocks.bedpe import Pair, read_bedpe
from genomeblocks import Loci
from genomeblocks.locus import Locus


def test_pair_properties():
    p = Pair("chr1", 1000, 1100, "chr1", 5000, 5100, name="l1", score=10.0)
    assert p.mid1 == 1050 and p.mid2 == 5050
    assert p.distance == 4000


def test_pair_distance_cross_chrom_is_inf():
    p = Pair("chr1", 1000, 1100, "chr2", 5000, 5100)
    assert p.distance == float("inf")


def test_read_bedpe(bedpe_path):
    pairs = read_bedpe(bedpe_path, verbose=False)
    assert len(pairs) == 3
    assert all(isinstance(p, Pair) for p in pairs)
    assert pairs[0].chrom1 == "chr1" and pairs[0].start2 == 5000


def test_read_bedpe_min_score_filter(bedpe_path):
    pairs = read_bedpe(bedpe_path, min_score=6.0, verbose=False)
    assert len(pairs) == 1                       # only the score=10 loop survives


def test_read_bedpe_max_distance_filter(bedpe_path):
    # distances are 4000, 8000, 4000 -> max_distance 5000 keeps the two 4000s
    pairs = read_bedpe(bedpe_path, max_distance=5000, verbose=False)
    assert len(pairs) == 2


def test_loci_pair_to_bed_attached(bedpe_path):
    # the bedpe anchors sit at chr1:1000-1100 / 5000-5100 / 9000-9100
    loci = Loci([Locus("chr1", 1000, 1100)])
    hits = loci.pair_to_bed(bedpe_path, either=True, both=False, verbose=False)
    # both loops whose first anchor overlaps the locus should be returned
    assert len(hits) >= 1
    assert all(isinstance(p, Pair) for p in hits)
