"""The table contract every container keeps (the acceptance checklist)."""
import numpy as np
import polars as pl
import pyarrow as pa
import pytest

import genomeblocks as gb


@pytest.fixture
def containers(cre, genes, pairs, arch, jaspar_path, tmp_path):
    (tmp_path / "t0.bed").write_text("chr1\t0\t200\nchr2\t100\t300\n")
    atlas = gb.Atlas.make(str(tmp_path / "t0.bed"), chromsizes={"chr1": 20_000, "chr2": 8_000}, bin_size=500,
                          verbose=False)
    return {"Loci": cre, "Genes": genes, "Pairs": pairs, "Architecture": arch, "Atlas": atlas,
            "Library": gb.load_motifs(jaspar_path)}


def test_every_container_is_a_table(containers):
    import duckdb
    for name, obj in containers.items():
        n_rows, n_cols = obj.shape
        assert n_cols == len(obj.columns), name
        assert obj.head(1) is not None and obj.describe().shape[1] == 1, name
        assert obj.to_pandas().shape[1] == n_cols, name
        assert obj.to_polars().shape == (n_rows, n_cols), name
        assert obj.to_arrow().num_rows == n_rows, name
        assert pl.DataFrame(obj).shape == (n_rows, n_cols), name                 # __arrow_c_stream__
        assert pa.table(obj).num_rows == n_rows, name
        assert duckdb.sql("select count(*) from obj").fetchall() == [(n_rows,)], name
        assert obj.__dataframe__() is not None, name
        html = obj._repr_html_()
        assert "<" in html and repr(obj), name
        assert len(obj) >= 0


def test_plotting_libraries_take_a_loci(cre):
    cre["score"] = np.arange(len(cre), dtype=float)
    sns = pytest.importorskip("seaborn")
    ax = sns.scatterplot(data=cre, x="start", y="score")
    assert ax is not None
    px = pytest.importorskip("plotly.express")
    assert px.scatter(cre, x="start", y="score") is not None
    alt = pytest.importorskip("altair")
    spec = alt.Chart(cre).mark_point().encode(x="start:Q", y="score:Q").to_dict()
    assert spec["mark"] in ("point", {"type": "point"})


def test_no_hidden_global_state(cre):
    a = gb.as_loci([("chrZ", 1, 2)])
    assert "chrZ" not in cre.genome                                   # separate tables, separate Genomes
    b = cre.intersect(a)
    assert len(b) == 0
    with gb.use_backend(intervals="genomeblocks"):
        pass
    from genomeblocks.backends import _scoped
    assert not _scoped


def test_backend_differences_do_not_leak(cre):
    other = gb.as_loci([("chr1", 1000, 1050), ("chr2", 0, 10_000), ("chr1", 2000, 2000)])
    results = {b: cre.overlap_pairs(other, backend=b) for b in gb.backends.families()["intervals"]
               if gb.backends.installed("intervals", b)}
    ref = results["genomeblocks"]
    for b, (qi, ri) in results.items():
        assert qi.tolist() == ref[0].tolist() and ri.tolist() == ref[1].tolist(), b
