"""Locus and region parsing."""
import pytest

from genomeblocks import Locus
from genomeblocks.locus import parse_region, parse_regions


def test_locus_basics():
    a, b = Locus("chr1", 100, 200, "+"), Locus("chr1", 150, 250)
    assert a.uid == "chr1:100-200(+)" and Locus.from_uid(a.uid) == a
    assert a.overlaps(b) and not a.overlaps(Locus("chr1", 200, 300))      # half-open
    assert a.length == 100 and a.center == 150
    assert a.distance_to(Locus("chr1", 300, 400)) == 200          # centre to centre
    assert hash(a) == hash(Locus("chr1", 100, 200, "+")) and a != b
    assert a.copy() == a and a.copy() is not a


@pytest.mark.parametrize("text, want", [
    ("chr1:1,000-2,000", ("chr1", 1000, 2000)),
    ("chr8:127.7-128.1 Mb", ("chr8", 127_700_000, 128_100_000)),
    ("chr2:5kb-12kb", ("chr2", 5000, 12000)),
    ("chr1:1_000–2_000", ("chr1", 1000, 2000)),
    (("chr3", 5, 9), ("chr3", 5, 9)),
    (Locus("chrX", 1, 2), ("chrX", 1, 2)),
])
def test_parse_region(text, want):
    assert parse_region(text) == want


def test_parse_region_rejects_garbage():
    with pytest.raises(ValueError, match="cannot parse region"):
        parse_region("kb")
    with pytest.raises(ValueError):
        parse_region("chr1")


def test_parse_regions_splits_a_two_locus_string():
    assert parse_regions("chr8:127.7-128.0 Mb chr1:1,000-2,000") == [("chr8", 127_700_000, 128_000_000),
                                                                      ("chr1", 1000, 2000)]
    assert parse_regions(("chr1", 1, 2)) == [("chr1", 1, 2)]
    assert parse_regions(["chr1:1-2", ("chr2", 3, 4), Locus("chr3", 5, 6)]) == [("chr1", 1, 2), ("chr2", 3, 4),
                                                                              ("chr3", 5, 6)]
