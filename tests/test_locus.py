from genomeblocks.locus import Locus, Exon, CDS, UTR


def test_uid_format():
    assert Locus("chr1", 100, 200, "+").uid == "chr1:100-200(+)"
    assert Locus("chr1", 100, 200).uid == "chr1:100-200(.)"


def test_length_and_center():
    l = Locus("chr1", 100, 200)
    assert l.length == 100
    assert l.center == 150


def test_distance_to_same_chrom():
    a = Locus("chr1", 100, 200)   # center 150
    b = Locus("chr1", 300, 500)   # center 400
    assert a.distance_to(b) == 250
    assert b.distance_to(a) == 250


def test_distance_to_cross_chrom_is_notimplemented():
    a = Locus("chr1", 100, 200)
    b = Locus("chr2", 100, 200)
    assert a.distance_to(b) is NotImplemented


def test_overlaps():
    a = Locus("chr1", 100, 200)
    assert a.overlaps(Locus("chr1", 150, 250))
    assert not a.overlaps(Locus("chr1", 200, 300))   # half-open: touching is not overlap
    assert not a.overlaps(Locus("chr2", 150, 250))


def test_equality_by_uid():
    assert Locus("chr1", 100, 200, "+") == Locus("chr1", 100, 200, "+")
    assert Locus("chr1", 100, 200, "+") != Locus("chr1", 100, 200, "-")


def test_ordering_within_chrom():
    assert Locus("chr1", 100, 200) < Locus("chr1", 300, 400)


def test_cross_chrom_ordering_is_consistent():
    # Regression: __gt__/__ge__ used to return `chrom < chrom` across
    # chromosomes, so `a > b` and `b > a` could both be True/False.
    a = Locus("chr1", 100, 200)
    b = Locus("chr2", 100, 200)
    assert a < b and b > a
    assert not (a > b) and not (b < a)
    assert a <= b and b >= a
    assert not (a >= b) and not (b <= a)


def test_locus_is_hashable_by_uid():
    # Regression: Locus was unhashable (dataclass eq without __hash__), so it
    # could not go in a set / dict key despite being value-like.
    a = Locus("chr1", 100, 200, "+")
    b = Locus("chr1", 100, 200, "+")
    c = Locus("chr1", 100, 200, "-")
    assert hash(a) == hash(b)
    assert len({a, b, c}) == 2          # a and b collapse; strand differs for c


def test_ne_against_non_locus():
    a = Locus("chr1", 100, 200)
    assert a != "chr1:100-200(.)"       # different type → not equal
    assert not (a == "chr1:100-200(.)")


def test_copy_is_independent():
    a = Locus("chr1", 100, 200, "+")
    b = a.copy()
    b.start = 150
    assert a.start == 100 and b.start == 150


def test_subclasses_carry_extra_fields():
    assert Exon("chr1", 1, 2, "+", exon_number=3).exon_number == 3
    assert isinstance(CDS("chr1", 1, 2, "+"), Exon)
    assert UTR("chr1", 1, 2, "+", type="5'").type == "5'"
