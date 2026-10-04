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


def _write_pairs(path, rows):
    with open(path, "w") as f:
        f.write("## pairs format v1.0\n")
        for i, (c1, p1, c2, p2) in enumerate(rows):
            f.write(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\n")


def _random_pairs(n=3000, seed=0):
    import random
    rng = random.Random(seed)
    sizes = {"chr1": 50_000, "chr2": 30_000, "chrM": 1_000}
    rows = []
    for _ in range(n):
        c1 = rng.choice(list(sizes))
        c2 = c1 if rng.random() < 0.7 else rng.choice(list(sizes))
        rows.append((c1, rng.randrange(1, sizes[c1]), c2, rng.randrange(1, sizes[c2])))
    return rows


# windows with gaps, out of order, none on chrM
_WINDOWS = [("chr2", 10_000, 20_000), ("chr1", 0, 5_000), ("chr1", 40_000, 50_000),
            ("chr1", 5_000, 9_000)]


def _window_of(windows, c, p):
    hit = [i for i, (wc, s, e) in enumerate(windows) if wc == c and s <= p < e]
    return hit[0] if hit else None


def test_count_pairs_matches_brute_force(tmp_path):
    from genomeblocks.bedpe import count_pairs
    rows = _random_pairs()
    path = tmp_path / "x.pairs"
    _write_pairs(path, rows)
    loci = Loci([Locus(*w) for w in _WINDOWS])

    want = {}
    for c1, p1, c2, p2 in rows:
        for (ca, pa), cb in (((c1, p1), c2), ((c2, p2), c1)):
            w = _window_of(_WINDOWS, ca, pa)
            if w is not None:
                want[(w, cb)] = want.get((w, cb), 0) + 1

    df = count_pairs(loci, str(path), chunksize=700, verbose=False)
    partners = sorted({cb for _, cb in want})
    assert list(df.columns) == ["chrom", "start", "end", "uid"] + partners
    for i in range(len(loci)):
        for cb in partners:
            assert df[cb].iloc[i] == want.get((i, cb), 0)

    one = count_pairs(loci, str(path), target_chrom="chr2", chunksize=700, verbose=False)
    assert one["count"].tolist() == [want.get((i, "chr2"), 0) for i in range(len(loci))]


def test_count_pairs_2d_matches_brute_force(tmp_path):
    import numpy as np
    from genomeblocks.bedpe import count_pairs_2d
    rows = _random_pairs(seed=1)
    path = tmp_path / "x.pairs"
    _write_pairs(path, rows)
    a = Loci([Locus(*w) for w in _WINDOWS])
    b_windows = [("chr1", 0, 25_000), ("chr2", 0, 30_000)]
    b = Loci([Locus(*w) for w in b_windows])

    for loci_b, wb in ((None, _WINDOWS), (b, b_windows)):
        want = np.zeros((len(a), len(wb)), dtype=np.int64)
        for c1, p1, c2, p2 in rows:
            for (ca, pa), (cb, pb) in (((c1, p1), (c2, p2)), ((c2, p2), (c1, p1))):
                i, j = _window_of(_WINDOWS, ca, pa), _window_of(wb, cb, pb)
                if i is not None and j is not None:
                    want[i, j] += 1
        got = count_pairs_2d(a, str(path), loci_b=loci_b, chunksize=700, verbose=False)
        assert got.dtype == np.int64
        assert np.array_equal(got.toarray(), want)
