"""Loci: construction, set algebra, merge / nearest, protocols, round trips."""
import numpy as np
import pandas as pd
import pytest

import genomeblocks as gb
from genomeblocks import Genome, Loci, as_loci

from conftest import installed_backends, random_loci, brute_pairs


# ── construction ──────────────────────────────────────────────────────────────

def test_make_bed_reads_header_lines_once(bed_path):
    L = Loci.make(bed_path, keep=True)
    assert len(L) == 7                                   # the '#' and track lines cost no data rows
    assert L.is_sorted
    assert L.strand.tolist() == ["+", "-", "+", ".", "-", "+", "+"]
    assert L.cols["name"].tolist() == ["p1", "p2", "p3", "p4", "p5", "p6", "p7"]
    assert L.cols["score"].tolist() == [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]


@pytest.mark.parametrize("backend", installed_backends("tables"))
def test_make_same_on_every_table_backend(bed_path, backend):
    a = Loci.make(bed_path, keep=True, backend=backend)
    b = Loci.make(bed_path, keep=True)
    assert a.equals(b, cols=True)


def test_make_keep_list_and_hash_in_field(tmp_path):
    p = tmp_path / "x.bed"
    p.write_text("chr1\t1\t5\tpeak#1\t10\t+\nchr1\t2\t6\tpeak#2\t20\t-\n")
    L = Loci.make(str(p), keep=["name"], sort=False)
    assert L.cols["name"].tolist() == ["peak#1", "peak#2"]
    assert L.strand.tolist() == ["+", "-"]
    assert list(L.cols) == ["name"]


def test_make_empty_and_bad_files(tmp_path):
    (tmp_path / "e.bed").write_text("")
    (tmp_path / "c.bed").write_text("# only a comment\n")
    (tmp_path / "s.bed").write_text("chr1 1 5\n")
    (tmp_path / "h.bed").write_text("chrom\tstart\tend\nchr1\t1\t5\n")
    assert len(Loci.make(str(tmp_path / "e.bed"))) == 0
    assert len(Loci.make(str(tmp_path / "c.bed"))) == 0
    with pytest.raises(ValueError, match="at least 3 columns"):
        Loci.make(str(tmp_path / "s.bed"))
    with pytest.raises(ValueError, match="Header lines"):
        Loci.make(str(tmp_path / "h.bed"))
    with pytest.raises(ValueError, match="Genes.make"):
        Loci.make("genes.gtf")


def test_from_records_uids_tile():
    L = Loci.from_records([("chr1", 5, 15, "+"), "chr2:100-200", gb.Locus("chr1", 50, 60)])
    assert L.uid.tolist() == ["chr1:5-15(+)", "chr2:100-200(.)", "chr1:50-60(.)"]
    assert Loci.from_uids(L.uid).equals(L)
    T = Loci.tile_genome({"chr1": 1000, "chr2": 250}, 400)
    assert T.to_records() == [("chr1", 0, 400, "."), ("chr1", 400, 800, "."), ("chr1", 800, 1000, "."),
                              ("chr2", 0, 250, ".")]


def test_null_coordinates_raise():
    with pytest.raises(ValueError, match="missing"):
        Loci.from_frame(pd.DataFrame({"chrom": ["chr1", None], "start": [1, 2], "end": [5, 6]}))
    with pytest.raises(ValueError, match="missing"):
        Loci.from_frame(pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [1.0, np.nan], "end": [5, 6]}))
    with pytest.raises(ValueError, match="missing"):
        Genome().encode([None])


# ── set algebra against brute force, every backend ────────────────────────────

@pytest.mark.parametrize("backend", installed_backends("intervals"))
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_overlaps_match_brute_force(backend, seed):
    rng = np.random.default_rng(seed)
    q, r = random_loci(rng, 60), random_loci(rng, 50)          # different Genomes, zero-length rows
    qi, ri = q.overlap_pairs(r, backend=backend)
    assert list(zip(qi.tolist(), ri.tolist())) == brute_pairs(q, r)
    ri2, qi2 = r.overlap_pairs(q, backend=backend)
    assert sorted(zip(qi2.tolist(), ri2.tolist())) == brute_pairs(q, r)      # symmetric
    want = np.zeros(len(q), bool)
    want[[i for i, _ in brute_pairs(q, r)]] = True
    assert (q.overlap_any(r, backend=backend) == want).all()
    assert q.intersect(r, backend=backend).equals(q.take(want))
    assert (q - r).equals(q.take(~want))


@pytest.mark.parametrize("backend", installed_backends("intervals"))
def test_zero_length_rule_is_uniform(backend):
    r = as_loci([("chr1", 5, 15)])
    inside, at_start = as_loci([("chr1", 10, 10)]), as_loci([("chr1", 5, 5)])
    assert inside.overlap_pairs(r, backend=backend)[0].tolist() == [0]
    assert r.overlap_pairs(inside, backend=backend)[0].tolist() == [0]
    assert at_start.overlap_pairs(r, backend=backend)[0].tolist() == []
    assert r.overlap_rows("chr1", 10, 10, backend=backend).tolist() == [0]
    assert r.overlap_rows("chr1", 5, 5, backend=backend).tolist() == []


@pytest.mark.parametrize("backend", [b for b in installed_backends("intervals") if b in ("genomeblocks", "bioframe", "pyranges", "bedtools")])
def test_merge_and_nearest_match_default(backend):
    rng = np.random.default_rng(3)
    q, r = random_loci(rng, 80), random_loci(rng, 70)
    m0, m1 = q.merge(), q.merge(backend=backend)
    assert m0.equals(m1)
    _, d0 = q.nearest(r)
    _, d1 = q.nearest(r, backend=backend)
    assert (d0 == d1).all()


def test_merge_fuses_book_ended_and_keeps_first_strand():
    L = as_loci([("chr1", 10, 20, "+"), ("chr1", 20, 30, "-"), ("chr1", 40, 50, "-"), ("chr2", 0, 5)])
    assert L.merge().to_records() == [("chr1", 10, 30, "+"), ("chr1", 40, 50, "-"), ("chr2", 0, 5, ".")]


def test_nearest_distances_and_misses():
    q = as_loci([("chr1", 100, 110), ("chr1", 200, 210), ("chr3", 0, 10)])
    r = as_loci([("chr1", 110, 120), ("chr1", 150, 160)])
    rows, dist = q.nearest(r)
    assert dist.tolist() == [0, 40, -1]                      # book-ended = 0, gap in bases, no chromosome
    assert rows.tolist() == [0, 1, -1]


def test_operators_and_region_lookups(cre):
    other = as_loci([("chr1", 1000, 1050), ("chr2", 0, 10_000)])
    assert len(cre & other) == 3 and len(cre - other) == 4
    assert len(cre | other) == len(cre) + len(other)
    assert cre.overlaps("chr1:1,000-2,000").to_records() == [("chr1", 900, 1100, "."), ("chr1", 1900, 2100, ".")]
    assert cre.overlap_rows("chr9", 0, 100).tolist() == []
    assert cre.slop(100)[0].start == 800 and cre.slop(2000)[0].start == 0


def test_sort_is_natural_order():
    L = as_loci([("chr10", 5, 6), ("chr2", 5, 6), ("chr1", 9, 10), ("chr1", 1, 2), ("chrX", 0, 1)])
    assert L.sort().chroms.tolist() == ["chr1", "chr1", "chr2", "chr10", "chrX"]
    assert L.sort().chrom_offsets["chr1"] == (0, 2)


# ── like a table ──────────────────────────────────────────────────────────────

def test_table_ergonomics(cre):
    assert cre.shape == (7, 4) and cre.columns == ["chrom", "start", "end", "strand"]
    assert len(cre.head(2)) == 2 and len(cre.tail(3)) == 3
    assert cre.describe().loc["rows", "value"] == 7
    assert "<table" in cre._repr_html_()
    cre["score"] = 1.5                                       # scalar broadcast
    assert cre["score"].tolist() == [1.5] * 7 and "score" in cre.columns
    with pytest.raises(KeyError):
        cre["start"] = np.zeros(7)
    assert cre[0].uid == "chr1:900-1100(.)" and cre["chr1:900-1100(.)"].row == 0
    assert cre[["chr2:500-600(.)"]].to_records() == [("chr2", 500, 600, ".")]
    assert len(cre[np.array([True] * 2 + [False] * 5)]) == 2


def test_protocols_hand_the_table_to_other_libraries(cre):
    cre["name"] = [f"p{i}" for i in range(7)]
    import polars as pl
    import pyarrow as pa
    assert pl.DataFrame(cre).shape == (7, 5)
    assert pa.table(cre).num_rows == 7
    assert pd.api.interchange.from_dataframe(cre).shape == (7, 5)
    import duckdb
    assert duckdb.sql("select count(*) from cre").fetchall() == [(7,)]
    import narwhals as nw
    assert nw.from_native(cre, eager_only=True).shape == (7, 5)
    assert np.asarray(cre)["start"].tolist() == cre.starts.tolist()


@pytest.mark.parametrize("kind", ["pandas", "polars", "arrow", "bioframe", "pyranges", "records", "dict"])
def test_round_trips(cre, kind):
    cre["name"] = [f"p{i}" for i in range(7)]
    cre["score"] = np.arange(7, dtype=float)
    if kind == "pandas":
        back = Loci.from_frame(cre.to_pandas())
    elif kind == "polars":
        back = Loci.from_frame(cre.to_polars().lazy())
    elif kind == "arrow":
        back = Loci.from_arrow(cre.to_arrow())
    elif kind == "bioframe":
        back = Loci.from_bioframe(cre.to_bioframe())
    elif kind == "pyranges":
        pytest.importorskip("pyranges")
        back = Loci.from_pyranges(cre.to_pyranges())          # sorted input: pyranges keeps the order
    elif kind == "records":
        back = Loci.from_records(cre.to_records())
        back.cols = {}
        cre.cols = {}
    else:
        back = as_loci({"chrom": cre.chroms, "start": cre.starts, "end": cre.ends, "strand": cre.strand,
                        "name": cre["name"], "score": cre["score"]})
    assert back.equals(cre, cols=True)


def test_pyranges_round_trip_keeps_strand_with_unstranded_rows():
    pytest.importorskip("pyranges")
    L = as_loci([("chr1", 5, 15, "-"), ("chr1", 50, 60, "."), ("chr2", 100, 200, "+")])
    assert Loci.from_pyranges(L.to_pyranges()).strand.tolist() == ["-", ".", "+"]


def test_parquet_and_bed_round_trip(cre, tmp_path):
    cre["name"] = [f"p{i}" for i in range(7)]
    cre["note"] = ["a", None, "c", "d", "e", "f", "g"]        # a null in an object column
    cre.save(str(tmp_path / "x.parquet"))
    back = Loci.load(str(tmp_path / "x.parquet"))
    assert back.equals(cre, cols=True)
    cre.to_bed(str(tmp_path / "x.bed"), name="name")
    assert Loci.make(str(tmp_path / "x.bed"), keep=True).cols["name"].tolist() == cre["name"].tolist()


def test_equals_with_nan_columns():
    L = Loci.from_frame(pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2], "score": [np.nan]}))
    assert L.equals(L.copy(), cols=True)


def test_anndata_round_trip(cre):
    ad = pytest.importorskip("anndata")
    X = np.random.default_rng(0).random((3, len(cre)))
    a = cre.to_anndata(X, obs=["s1", "s2", "s3"])
    assert a.shape == (3, 7) and list(a.var.index)[0] == "chr1:900-1100"
    assert Loci.from_anndata(a).equals(cre)
    assert Loci.from_anndata(ad.AnnData(X, var=pd.DataFrame(index=["chr1:1-2", "chr1_3_4"] + list(range(5))[:0] + [f"chr2:{i}-{i+1}" for i in range(5)]))).to_records()[1] == ("chr1", 3, 4, ".")


def test_sequences_and_fasta(cre, fasta_path, tmp_path):
    seqs = cre.sequences(fasta_path)
    assert [len(s) for s in seqs] == cre.lengths.tolist()
    s_minus = as_loci([("chr1", 1000, 1012, "-")]).sequences(fasta_path, strand=True)[0]
    assert s_minus == "ACGT" * 3                              # the planted site is its own reverse complement
    out = cre.to_fasta(str(tmp_path / "s.fa"), fasta_path, r=10)
    assert open(out).read().count(">") == 7
    recs = cre.to_seqrecords(fasta_path)
    assert recs[0].id == cre.uid[0]
