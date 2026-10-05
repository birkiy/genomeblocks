"""The boundary: as_loci on every input, file readers, informative errors."""
import numpy as np
import pandas as pd
import polars as pl
import pytest

from genomeblocks import Loci, as_loci


def test_as_loci_accepts_every_table_kind(cre):
    df = cre.to_pandas()
    for x in (cre, df, pl.from_pandas(df), pl.from_pandas(df).lazy(), cre.to_arrow(), cre.to_bioframe(),
              {"chrom": cre.chroms, "start": cre.starts, "end": cre.ends}, cre.to_numpy(),
              (cre.chroms, cre.starts, cre.ends), cre.to_records(), cre.uid.tolist()):
        assert as_loci(x).equals(cre), type(x)
    assert as_loci("chr1:1,000-2,000").to_records() == [("chr1", 1000, 2000, ".")]
    assert as_loci(cre[0]).equals(cre.head(1))
    assert len(as_loci([])) == 0
    assert as_loci([cre.head(2), cre.tail(2)]).equals(cre.head(2) + cre.tail(2))


def test_as_loci_lenient_column_names():
    for cols in (("Chromosome", "Start", "End"), ("seqnames", "start", "end"), ("#chrom", "chromStart", "chromEnd"),
                 ("chr", "begin", "stop")):
        df = pd.DataFrame({cols[0]: ["chr1"], cols[1]: [5], cols[2]: [9]})
        assert as_loci(df).to_records() == [("chr1", 5, 9, ".")]
    headerless = pd.DataFrame([["chr1", 5, 9, "x"]])                 # text + int + int first: by position
    assert as_loci(headerless).to_records() == [("chr1", 5, 9, ".")]
    with pytest.raises(ValueError, match="expected columns"):
        as_loci(pd.DataFrame({"a": ["chr1"], "b": ["x"], "c": [1]}))


def test_as_loci_errors_name_the_input():
    with pytest.raises(TypeError, match="dict of columns"):
        as_loci({"a": 1})
    with pytest.raises(ValueError, match="different lengths"):
        as_loci({"chrom": ["chr1", "chr1"], "start": [1], "end": [2, 3]})
    with pytest.raises(TypeError, match="Series"):
        as_loci(pd.Series([1, 2]))
    import pyarrow as pa
    with pytest.raises(TypeError, match="ChunkedArray"):        # streams Arrow data, but one column, not a table
        as_loci(pa.chunked_array([[1, 2]]))                     # (what a pandas >= 3 Series does too)
    with pytest.raises(FileNotFoundError, match="not a region"):
        as_loci("no_such_file.bed")
    with pytest.raises(TypeError, match="cannot read intervals"):
        as_loci(object())


def test_file_inputs_keep_file_order(tmp_path, bed_path):
    (tmp_path / "o.bed").write_text("chr2\t1\t2\nchr1\t1\t2\n")
    assert as_loci(str(tmp_path / "o.bed")).chroms.tolist() == ["chr2", "chr1"]
    assert Loci.make(str(tmp_path / "o.bed")).chroms.tolist() == ["chr1", "chr2"]    # make() sorts by default
    assert len(as_loci(bed_path)) == 7


def test_csv_tsv_readers(tmp_path):
    (tmp_path / "u.tsv").write_text("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\nchr1\t5\t15\tb\t0\t-\n")
    L = as_loci(str(tmp_path / "u.tsv"))
    assert L.strand.tolist() == ["-"] and list(L.cols) == ["name", "score"]
    (tmp_path / "h.tsv").write_text("chr1\t5\t15\tb\t0\t-\n")
    L = as_loci(str(tmp_path / "h.tsv"))
    assert L.strand.tolist() == ["-"] and list(L.cols) == ["name", "score"]
    (tmp_path / "c.csv").write_text("# a comment\nchrom,start,end\nchr2,100,200\nchr1,1,2\n")
    assert as_loci(str(tmp_path / "c.csv")).chroms.tolist() == ["chr2", "chr1"]
    (tmp_path / "x.csv").write_text("a,b\n1,2\n")
    with pytest.raises(ValueError, match="at least 3"):
        as_loci(str(tmp_path / "x.csv"))


def test_parquet_path_and_genome_argument(cre, tmp_path):
    cre.save(str(tmp_path / "c.parquet"))
    g = cre.genome
    L = as_loci(str(tmp_path / "c.parquet"), genome=g)
    assert L.genome is g and L.equals(cre)


def test_arrow_export_with_nulls_in_object_columns():
    df = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [1, 7], "end": [5, 9], "name": ["a", np.nan]})
    L = Loci.from_frame(df)
    assert L.to_polars()["name"].to_list() == ["a", None]
    assert L.to_arrow().column("name").null_count == 1


def test_cube_converters(cre):
    S = np.random.default_rng(0).random((len(cre), 2, 5))
    from genomeblocks import interop
    pytest.importorskip("xarray")
    da = interop.cube_to_xarray(S, cre, ["a", "b"], flank=250)
    assert da.dims == ("region", "track", "bin") and da.coords["bin"].values[0] == -200
    df = interop.cube_to_pandas(S, cre, ["a", "b"])
    assert df.shape[0] == len(cre)
    pytest.importorskip("anndata")
    ad = interop.cube_to_anndata(S, cre, ["a", "b"])
    assert ad.shape == (2, len(cre))


def test_liftover_needs_pyliftover_or_works(cre, tmp_path):
    try:
        import pyliftover  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="pip install pyliftover"):
            cre.liftover(str(tmp_path / "x.chain"))
        return
    chain = tmp_path / "t.chain"
    chain.write_text("chain 1000 chr1 20000 + 0 20000 chrA 30000 + 10000 30000 1\n20000\n\n")
    out = cre.liftover(str(chain), verbose=False)
    assert out.chroms.tolist()[:1] == ["chrA"] and out.starts[0] == cre.starts[0] + 10000
