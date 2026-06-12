import pytest

from genomeblocks import Atlas, Loci
from genomeblocks.locus import Locus


@pytest.fixture
def atlas_dir(tmp_path):
    """Two synthetic ChIP-Atlas-style BED tracks."""
    (tmp_path / "trackA.bed").write_text(
        "chr1\t1000\t2000\n"
        "chr1\t3000\t4000\n"
    )
    (tmp_path / "trackB.bed").write_text(
        "chr1\t50000\t51000\n"
    )
    return str(tmp_path)


def test_make_and_search(atlas_dir):
    atlas = Atlas.make(atlas_dir, chromsizes={"chr1": 100000},
                       bin_size=1000, workers=1, verbose=False)
    assert set(atlas.track_names) == {"trackA", "trackB"}

    query = Loci([Locus("chr1", 1000, 2000)])     # overlaps trackA only
    df = atlas.search(query)
    for col in ["name", "overlaps", "log2_odds", "giggle_score"]:
        assert col in df.columns
    by_name = df.set_index("name")
    assert by_name.loc["trackA", "overlaps"] >= 1
    # trackA is more enriched in the query than trackB
    assert by_name.loc["trackA", "giggle_score"] >= by_name.loc["trackB", "giggle_score"]


def test_attach_meta_headerless(tmp_path, atlas_dir):
    meta = tmp_path / "meta.tsv"
    meta.write_text("trackA\tFOX\ntrackB\tAR\n")    # header-less id\tantigen
    atlas = Atlas.make(atlas_dir, chromsizes={"chr1": 100000},
                       bin_size=1000, workers=1, verbose=False,
                       meta=str(meta), meta_columns=["id", "antigen"],
                       meta_id_col="id")
    df = atlas.search(Loci([Locus("chr1", 1000, 2000)]))
    assert "antigen" in df.columns
    assert df.set_index("name").loc["trackA", "antigen"] == "FOX"
