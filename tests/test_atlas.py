"""Atlas: build, search, bootstrap, metadata, persistence, table protocols."""
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
    # brute force: query bins (500 bp) shared with each track's bins
    def bins(L):
        return {(c, b) for c, s, e in zip(L.chroms, L.starts, L.ends) for b in range(s // 500, (e - 1) // 500 + 1)}
    qb = bins(cre)
    want = {f"t{k}": len(qb & bins(as_loci(str(tmp_path / f"t{k}.bed")))) for k in range(3)}
    for q in (cre, cre.to_pandas(), cre.to_polars(), str(tmp_path / "t0.bed")):
        df = atlas.search(q)
        assert set(df.columns) >= {"name", "overlaps", "log2_odds", "p", "giggle_score"} and len(df) == 3
    df = atlas.search(cre).set_index("name")
    assert df["overlaps"].to_dict() == want and (df["n_query_bins"] == len(qb)).all()
    assert (df["p"] > 0).all() and (df["p"] <= 1).all()
    assert df["giggle_score"].is_monotonic_decreasing
    df = atlas.search(cre, ref=str(tmp_path / "t1.bed"))
    assert "n_ref_bins" in df.columns
    assert cre.enrich(atlas, ref=str(tmp_path / "t1.bed")).equals(df)
    with pytest.raises(ValueError, match="no bins"):
        atlas.search(as_loci([("chrZ", 1, 2)]))


def test_bootstrap_modes(atlas, cre, tmp_path):
    single = atlas.bootstrap(cre, n=3, seed=0, verbose=False)
    assert set(single.columns) >= {"name", "observed", "expected", "z", "p_emp"} and len(single) == 3

    def bins(c, s, e):
        return {(c, b) for b in range(s // 500, (e - 1) // 500 + 1)}

    def all_bins(L):
        return set().union(*(bins(c, s, e) for c, s, e in zip(L.chroms, L.starts, L.ends)))
    want = {}
    for k in range(3):
        tb = all_bins(as_loci(str(tmp_path / f"t{k}.bed")))
        want[f"t{k}"] = sum(1 for c, s, e in zip(cre.chroms, cre.starts, cre.ends) if bins(c, s, e) & tb)
    assert single.set_index("name")["observed"].to_dict() == want     # observed = query regions sharing a bin, each once
    assert (single["expected"] >= 0).all() and ((single["p_emp"] > 0) & (single["p_emp"] <= 1)).all()
    groups = atlas.bootstrap({"a": cre.to_pandas(), "b": str(tmp_path / "t1.bed")}, n=3, seed=0, verbose=False)
    assert "group" in groups.columns and len(groups) == 6
    cols = atlas.bootstrap({"chrom": ["chr1"], "start": [100], "end": [900]}, n=2, verbose=False)   # a dict of columns is one query
    assert "group" not in cols.columns
    pooled = atlas.bootstrap(cre, pool=cre.to_polars(), n=2, sample=3, seed=1, verbose=False)
    assert len(pooled) == 3
    assert cre.enrich_mc(atlas, n=2, verbose=False).shape[0] == 3


def test_bootstrap_counts_regions_so_a_clustered_query_is_not_read_as_depleted(tmp_path):
    """A query whose regions share bins (240 regions on the track in 120 bins, 160 off it) drawn from a pool
    where 500 of 1,200 regions are on the track: 1.44x enriched. Counting distinct bins (up to 2.0.1) read it
    as depleted (120 bins against ~144 for random draws); counting regions gives 240 against 400 x 500/1200."""
    (tmp_path / "track.bed").write_text("chr1\t0\t300000\n")              # bins 0-299 of 1 kb
    A = Atlas.make([str(tmp_path / "track.bed")], chromsizes={"chr1": 2_000_000}, bin_size=1000, verbose=False, workers=1)
    single = [("chr1", i * 1000 + 100, i * 1000 + 300) for i in range(1000)]
    twins = [("chr1", i * 1000 + 500, i * 1000 + 700) for i in range(200)]     # a second region in bins 0-199
    pool = as_loci(single + twins)
    query = as_loci(single[:120] + twins[:120] + single[300:460])
    r = A.bootstrap(query, pool=pool, n=200, seed=0, verbose=False).iloc[0]
    assert r.observed == 240 and abs(r.expected - 400 * 500 / 1200) < 4 and r.z > 5 and r.log2fc > 0.4
    sub = A.bootstrap(query, pool=pool, n=200, sample=200, seed=0, verbose=False).iloc[0]
    assert abs(sub.observed / sub.expected - 240 / (400 * 500 / 1200)) < 0.05 and sub.z > 2
    shuf = A.bootstrap(query, n=50, seed=0, verbose=False).iloc[0]                # 15 % of the chromosome is the track
    assert shuf.observed == 240 and abs(shuf.expected - 0.15 * 400) < 6
    # a region spanning several bins counts once however many of them the track covers
    long_q = as_loci([("chr1", 100, 5_100), ("chr1", 299_500, 300_500), ("chr1", 400_000, 403_000)])
    assert A._region_counts(A._intervals_to_bin_ranges(long_q)).tolist() == [2]
    assert A._region_counts(A._intervals_to_bin_ranges(as_loci(single[:3] * 2))).tolist() == [6]   # repeats count


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


def test_search_columns_are_the_documented_ones(atlas, cre, tmp_path):
    df = atlas.search(cre)
    assert df.columns.tolist() == ["name", "n_query_bins", "track_n_bins", "track_n_peaks", "overlaps",
                                   "log2_odds", "p", "giggle_score"]
    assert df.set_index("name")["track_n_bins"].to_dict() == dict(zip(atlas.track_names, atlas.track_n_bins.tolist()))
    assert cre.enrich(atlas).columns.tolist() == df.columns.tolist()
    ref = atlas.search(cre, ref=str(tmp_path / "t1.bed"))
    assert ref.columns.tolist() == ["name", "n_query_bins", "n_ref_bins", "track_n_bins", "track_n_peaks", "overlaps",
                                    "log2_odds", "p", "giggle_score"]
    assert ref.set_index("name")["track_n_bins"].to_dict() == df.set_index("name")["track_n_bins"].to_dict()


def test_load_restores_missing_metadata_exactly(atlas, cre, tmp_path):
    atlas.attach_meta(pd.DataFrame({"id": ["t0", "t1"], "factor": ["AR", "FOXA1"]}), id_col="id")
    assert pd.isna(atlas.meta.loc["t2", "factor"])
    atlas.save(str(tmp_path / "at.npz"))
    back = Atlas.load(str(tmp_path / "at.npz"))
    assert back.meta.equals(atlas.meta)
    assert back.search(cre).equals(atlas.search(cre))
    assert back.to_pandas().equals(atlas.to_pandas())
    assert back.bootstrap(cre, n=2, seed=0, verbose=False).equals(atlas.bootstrap(cre, n=2, seed=0, verbose=False))


def test_atlas_as_a_table(atlas):
    assert atlas.shape == (3, 3) and atlas.columns == ["name", "n_peaks", "n_bins"]
    assert list(atlas) == ["t0", "t1", "t2"] and atlas["t1"]["n_peaks"] == 7
    assert len(atlas.head(2)) == 2 and atlas.describe().loc["tracks", "value"] == 3
    import polars as pl
    assert pl.DataFrame(atlas).shape == (3, 3) and "<table" in atlas._repr_html_()
