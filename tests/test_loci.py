import genomeblocks.loci as loci_mod
from genomeblocks import Loci
from genomeblocks.locus import Locus


def test_pure_python_interval_index(monkeypatch):
    # cgranges isn't on PyPI, so pip-only installs use the pure-Python index.
    # Force it and check overlap-based ops are still correct.
    monkeypatch.setattr(loci_mod, "_get_cgranges", lambda: None)

    a = Loci([Locus("chr1", s, s + 100) for s in (0, 100, 250, 400)])
    b = Loci([Locus("chr1", s, s + 100) for s in (50, 300)])  # 50-150, 300-400

    assert isinstance(a.cgr, loci_mod._PyIntervalIndex)   # fallback engaged
    # 0-100, 100-200, 250-350 overlap b; 400-500 does not (half-open at 400)
    assert sorted(l.start for l in a.intersect(b)) == [0, 100, 250]
    assert sorted(l.start for l in a.difference(b)) == [400]
    # a far-away interval overlaps nothing in b
    assert len(Loci([Locus("chr1", 5000, 5100)]).intersect(b)) == 0
    # different chromosome never overlaps
    assert len(Loci([Locus("chr2", 60, 140)]).intersect(b)) == 0


def test_make_from_bed(bed_path):
    loci = Loci.make(bed_path)
    assert len(loci) == 4                       # comment line skipped
    assert loci[0].chrom == "chr1" and loci[0].start == 1000
    assert loci[2].strand == "-"


def test_uid_lookup_and_index(cre):
    uid = cre[0].uid
    assert cre[uid].start == 1000
    assert cre.uids[uid] == 0


def test_intersect_keeps_lhs_overlapping(cre):
    other = Loci([Locus("chr1", 1050, 1060), Locus("chr1", 9020, 9030)])
    hit = cre & other                            # cre loci overlapping `other`
    starts = sorted(l.start for l in hit)
    assert starts == [1000, 9000]


def test_difference(cre):
    other = Loci([Locus("chr1", 1050, 1060)])
    diff = cre - other
    assert len(diff) == 4
    assert all(l.start != 1000 for l in diff)


def test_concat_operators(cre):
    a = Loci(cre[:2])
    b = Loci(cre[2:])
    assert len(a + b) == len(cre)
    assert len(a | b) == len(cre)                # | concatenates (no dedup)


def test_symmetric_difference(cre):
    a = Loci(cre[:3])
    b = Loci(cre[2:])
    sym = a ^ b
    starts = sorted(l.start for l in sym)
    assert starts == [1000, 5000, 13000, 17000]  # the 9000 shared one drops out


def test_sort_and_merge_fuses_overlaps():
    raw = Loci([
        Locus("chr1", 1100, 1300),
        Locus("chr1", 1000, 1200),               # overlaps the first
        Locus("chr1", 5000, 5200),
    ])
    merged = raw.sort().merge()
    assert len(merged) == 2
    assert merged[0].start == 1000 and merged[0].end == 1300


def test_slop_expands_and_clamps():
    loci = Loci([Locus("chr1", 50, 100)])
    s = loci.slop(80)
    assert s[0].start == 0 and s[0].end == 180   # start clamped at 0


def test_overlaps_query_forms(cre):
    by_coords = cre.overlaps("chr1", 990, 1010)
    by_locus = cre.overlaps(Locus("chr1", 990, 1010))
    assert len(by_coords) == 1 and len(by_locus) == 1
    assert by_coords[0].start == 1000


def test_map_returns_uid_lists(cre):
    other = Loci([Locus("chr1", 1050, 1060)])
    m = cre.map(other)
    assert m[cre[0].uid] == [other[0].uid]
    assert m[cre[1].uid] == []


def test_nearest_dataframe(cre):
    other = Loci([Locus("chr1", 1300, 1400)])
    df = cre.nearest(other)
    assert "Distance" in df.columns
    assert len(df) == len(cre)


def test_to_frame_columns(cre):
    df = cre.to_frame()
    assert list(df.columns) == ["Chr", "Start", "End", "Strand", "Name"]
    assert len(df) == len(cre)
