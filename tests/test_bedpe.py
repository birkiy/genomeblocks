"""Pairs (BEDPE) and pairs-file counting."""
import numpy as np
import pandas as pd
import polars as pl
import pytest

from genomeblocks import Pairs, as_loci
from genomeblocks.bedpe import _detect_pairs_format, as_pairs, count_pairs, count_pairs_2d, read_pairs_chunks

from conftest import installed_backends


def test_make_and_columns(pairs):
    assert len(pairs) == 4 and pairs.is_cis.tolist() == [True, True, False, True]
    assert pairs.cols["name"].tolist() == ["l1", "l2", "l3", "l4"] and pairs.cols["score"].tolist() == [5, 3, 1, 2]
    assert pairs.a.strand.tolist() == ["+", "+", "-", "+"]
    assert pairs.distance.tolist()[2] == np.inf
    assert pairs.columns == ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
    assert pairs.shape == (4, 10)


@pytest.mark.parametrize("backend", installed_backends("tables"))
def test_header_lines_and_backends(tmp_path, backend):
    p = tmp_path / "h.bedpe"
    p.write_text("#c1\ts1\te1\tc2\ts2\te2\nchr1\t100\t200\tchr1\t1000\t1100\nchr1\t300\t400\tchr2\t500\t600\n")
    P = Pairs.make(str(p), backend=backend)
    assert len(P) == 2 and P.cols == {}
    with pytest.raises(ValueError, match="at least 6"):
        (tmp_path / "x.bedpe").write_text("chr1\t1\t2\n")
        Pairs.make(str(tmp_path / "x.bedpe"))


def test_from_frame_round_trips_and_validation(pairs, tmp_path):
    for df in (pairs.to_pandas(), pairs.to_polars(), pairs.to_arrow()):
        back = Pairs.from_frame(df)
        assert back.a.equals(pairs.a) and back.b.equals(pairs.b) and back.cols["score"].tolist() == [5, 3, 1, 2]
    out = pairs.to_bedpe(str(tmp_path / "o.bedpe"))
    assert Pairs.make(out).a.equals(pairs.a)
    pairs.save(str(tmp_path / "p.parquet"))
    assert Pairs.load(str(tmp_path / "p.parquet")).b.equals(pairs.b)
    assert as_pairs(str(tmp_path / "p.parquet")).a.equals(pairs.a)
    with pytest.raises(ValueError, match="BEDPE frame"):
        Pairs.from_frame(pd.DataFrame({"a": [1]}))


def test_filters_and_overlaps(pairs, cre):
    assert len(pairs.filter(min_score=3)) == 2
    assert len(pairs.filter(max_distance=5000)) == 2                        # trans (inf) dropped
    m1, m2 = pairs.anchors_overlap(as_loci([("chr1", 900, 1100)]))
    assert m1.tolist() == [True, False, True, False] and not m2.any()
    assert len(pairs.overlapping(cre.to_pandas(), both=True)) == 4
    for backend in installed_backends("intervals"):
        assert (pairs.anchors_overlap(cre, backend=backend)[0] == m1 | True).all()
    assert pairs[0][0].uid == "chr1:900-1100(+)" and pairs["score"][0] == 5
    assert len(pairs.head(2)) == 2 and pairs.describe().loc["trans", "value"] == 1
    assert "<table" in pairs._repr_html_()
    assert pl.DataFrame(pairs).shape == (4, 10)


def test_pairs_file_formats(tmp_path):
    (tmp_path / "t.pairs").write_text("## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
                                      "r1\tchr1\t100\tchr2\t200\t+\t-\n")
    (tmp_path / "t.avp").write_text("r1\tchr1\t100\t+\tchr2\t200\t-\t300\tHIC_1\tHIC_2\t42\t42\n")
    (tmp_path / "t.juicer").write_text("r0\t0\tchr1\t123\t1\t16\tchr2\t456\t2\t60\t60\n")
    assert _detect_pairs_format(str(tmp_path / "t.pairs")) == "pairs"
    assert _detect_pairs_format(str(tmp_path / "t.avp")) == "allvalidpairs"
    assert _detect_pairs_format(str(tmp_path / "t.juicer")) == "juicer"
    for f in ("t.pairs", "t.avp", "t.juicer"):
        ch = next(read_pairs_chunks(str(tmp_path / f)))
        assert list(ch.columns) == ["chrom1", "pos1", "chrom2", "pos2"] and ch["chrom2"].tolist() == ["chr2"]
    with pytest.raises(ValueError, match="Unknown pairs format"):
        next(read_pairs_chunks(str(tmp_path / "t.pairs"), format="nope"))


def test_count_pairs_matches_brute_force(tmp_path):
    rng = np.random.default_rng(0)
    chroms = rng.choice(["chr1", "chr2"], 400)
    p1, p2 = rng.integers(0, 10_000, 400), rng.integers(0, 10_000, 400)
    c2 = np.where(rng.random(400) < 0.3, np.where(chroms == "chr1", "chr2", "chr1"), chroms)
    lines = "".join(f"r{i}\t{a}\t{x}\t+\t{b}\t{y}\t-\t100\n" for i, (a, x, b, y) in enumerate(zip(chroms, p1, c2, p2)))
    (tmp_path / "x.avp").write_text(lines)
    L = as_loci([("chr1", s, s + 1000) for s in range(0, 10_000, 2000)] + [("chr2", 0, 5000)])
    df = count_pairs(L, str(tmp_path / "x.avp"), verbose=False)
    want = np.zeros((len(L), 2), int)
    for a, x, b, y in zip(chroms, p1, c2, p2):
        for (ca, pa, cb) in ((a, x, b), (b, y, a)):
            for i in range(len(L)):
                if L.chroms[i] == ca and L.starts[i] <= pa < L.ends[i]:
                    want[i, 0 if cb == "chr1" else 1] += 1
    assert df[["chr1", "chr2"]].to_numpy().tolist() == want.tolist()
    M = count_pairs_2d(L, str(tmp_path / "x.avp"), verbose=False)
    assert M.shape == (len(L), len(L)) and (M.toarray() == M.toarray().T).all()
    only = L.count_pairs(str(tmp_path / "x.avp"), target_chrom="chr2", verbose=False)
    assert len(only) == len(L) and "count" in only.columns and only["count"].tolist() == want[:, 1].tolist()
