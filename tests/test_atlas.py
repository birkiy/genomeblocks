"""Atlas: build, search, bootstrap, metadata, persistence, table protocols."""
import numpy as np
import pandas as pd
import pytest

import genomeblocks as gb
from genomeblocks import Atlas, as_loci

from conftest import CHROM_SIZES


@pytest.fixture
def atlas(tmp_path):
    for k in range(3):
        (tmp_path / f"t{k}.bed").write_text("".join(f"chr1\t{s}\t{s + 200}\n" for s in range(k * 1000, 18_000, 3000))
                                            + "chr2\t100\t300\n")
    return Atlas.make([str(tmp_path / f"t{k}.bed") for k in range(3)], chromsizes=CHROM_SIZES, bin_size=500,
                      verbose=False, workers=1)


def test_make_inputs(tmp_path, atlas, sizes_path):
    assert repr(atlas).startswith("Atlas(tracks=3") and atlas.track_names == ["t0", "t1", "t2"]
    a = Atlas.make(str(tmp_path), chromsizes=sizes_path, bin_size=500, verbose=False)
    assert a.track_names == ["t0", "t1", "t2"]
    assert Atlas.make(str(tmp_path / "t0.bed"), chromsizes=gb.Genome.from_sizes(CHROM_SIZES), bin_size=500,
                      verbose=False).n_bins == atlas.n_bins
    assert Atlas.make(str(tmp_path / "t0.bed"), chromsizes=pd.Series(CHROM_SIZES), bin_size=500, verbose=False).shape == (1, 3)
    with pytest.raises(ValueError, match="No BED files"):
        Atlas.make(str(tmp_path / "none" / "*.bed"), chromsizes=CHROM_SIZES)


def test_bin_ranges_are_columnar(atlas, cre):
    r = atlas._intervals_to_bin_ranges(cre)
    assert r.tolist() == [[1, 3], [3, 5], [9, 11], [19, 21], [21, 23], [41, 42], [50, 51]]
    assert atlas._intervals_to_bin_ranges(cre.to_pandas()).tolist() == r.tolist()
    assert atlas._intervals_to_bin_ranges(as_loci([("chrZ", 1, 2)])).shape == (0, 2)


def test_search_accepts_any_interval_input(atlas, cre, tmp_path):
    for q in (cre, cre.to_pandas(), cre.to_polars(), str(tmp_path / "t0.bed")):
        df = atlas.search(q)
        assert set(df.columns) >= {"name", "overlaps", "log2_odds", "p", "giggle_score"} and len(df) == 3
    df = atlas.search(cre, ref=str(tmp_path / "t1.bed"))
    assert "n_ref_bins" in df.columns
    assert cre.enrich(atlas).shape == df.shape
    with pytest.raises(ValueError, match="no bins"):
        atlas.search(as_loci([("chrZ", 1, 2)]))


def test_bootstrap_modes(atlas, cre, tmp_path):
    single = atlas.bootstrap(cre, n=3, seed=0, verbose=False)
    assert set(single.columns) >= {"name", "observed", "expected", "z", "p_emp"} and len(single) == 3
    groups = atlas.bootstrap({"a": cre.to_pandas(), "b": str(tmp_path / "t1.bed")}, n=3, seed=0, verbose=False)
    assert "group" in groups.columns and len(groups) == 6
    cols = atlas.bootstrap({"chrom": ["chr1"], "start": [100], "end": [900]}, n=2, verbose=False)   # a dict of columns is one query
    assert "group" not in cols.columns
    pooled = atlas.bootstrap(cre, pool=cre.to_polars(), n=2, sample=3, seed=1, verbose=False)
    assert len(pooled) == 3
    assert cre.enrich_mc(atlas, n=2, verbose=False).shape[0] == 3


def test_meta_and_persistence(atlas, tmp_path):
    meta = pd.DataFrame({"id": ["t0", "t1", "t9"], "factor": ["AR", "FOXA1", "X"]})
    atlas.attach_meta(meta, id_col="id")
    assert atlas.meta.loc["t0", "factor"] == "AR" and pd.isna(atlas.meta.loc["t2", "factor"])
    assert "factor" in atlas.search(as_loci([("chr1", 0, 3000)])).columns
    (tmp_path / "meta.tsv").write_text("t0\tAR\nt1\tFOXA1\n")
    atlas.attach_meta(str(tmp_path / "meta.tsv"), columns=["id", "factor"])
    assert atlas.meta.loc["t1", "factor"] == "FOXA1"
    atlas.save(str(tmp_path / "at.npz"))
    back = Atlas.load(str(tmp_path / "at.npz"))
    assert back.track_names == atlas.track_names and back.meta.loc["t0", "factor"] == "AR"
    assert (back.M != atlas.M).nnz == 0


def test_atlas_as_a_table(atlas):
    assert atlas.shape == (3, 3) and atlas.columns == ["name", "n_peaks", "n_bins"]
    assert list(atlas) == ["t0", "t1", "t2"] and atlas["t1"]["n_peaks"] == 7
    assert len(atlas.head(2)) == 2 and atlas.describe().loc["tracks", "value"] == 3
    import polars as pl
    assert pl.DataFrame(atlas).shape == (3, 3) and "<table" in atlas._repr_html_()
