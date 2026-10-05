"""Architecture: build, weights, normalisation, annotation, graph backends, drawing."""
import numpy as np
import pandas as pd
import pytest

import genomeblocks as gb
from genomeblocks import Architecture, Genes, Pairs, as_loci

from conftest import installed_backends


def test_make_links_the_right_cres(cre, pairs, bedpe_path):
    A = Architecture.make(cre, pairs, r=100, verbose=False)
    assert list(A) == [(0, 2), (1, 4), (5, 6), (0, 5)]              # cis blocks in genome order, then trans
    assert A.blocks == {"chr1": (0, 2), "chr2": (2, 3), "trans": (3, 4)}
    assert A.n_links == 4 and A.n_trans == 1 and A.degree.tolist() == [2, 1, 1, 0, 1, 2, 1]
    assert list(A.ep) == ["w"]                                      # n and d only after normalize()
    for src in (bedpe_path, pairs.to_pandas(), pairs.to_polars()):
        assert list(Architecture.make(cre, src, r=100, verbose=False)) == list(A)
    assert Architecture.make(cre, pairs, r=100, trans=False, verbose=False).n_trans == 0
    assert Architecture.make(cre, pairs, r=100, dmax=5000, verbose=False).n_links == 3


@pytest.mark.parametrize("backend", installed_backends("intervals"))
def test_make_same_on_every_interval_backend(cre, pairs, backend):
    assert list(Architecture.make(cre, pairs, r=100, verbose=False, backend=backend)) == [(0, 2), (1, 4), (5, 6), (0, 5)]


def test_prune_and_annotate_need_normalize(cre, pairs, genes):
    A = Architecture.make(cre, pairs, r=100, verbose=False)
    with pytest.raises(ValueError, match="normalize"):
        A.prune(verbose=False)
    with pytest.raises(ValueError, match="normalize"):
        A.annotate(genes, verbose=False)
    A.ep["w"][:] = [5, 3, 2, 1]
    A.normalize(verbose=False)
    assert A.ep["d"].tolist()[3] == np.inf and A.ep["n"][3] == 1.0   # one trans edge: w / mean trans w
    assert A.prune(verbose=False).n_links == 4                          # nothing co-located


def test_normalize_with_few_edges():
    L = as_loci([("chr1", 0, 100), ("chr1", 1000, 1100), ("chr2", 0, 100), ("chr2", 1000, 1100)])
    A = Architecture.from_edges(L, [0, 0, 1], [1, 2, 3], w=[5.0, 2.0, 6.0])
    A.normalize(verbose=False)                                           # one cis edge: no fit, no crash
    assert np.isnan(A.fit["alpha"]) and np.isfinite(A.ep["n"]).all()


def test_annotate_strength_hubs_support(arch, genes):
    arch.annotate(genes, verbose=False)
    assert arch.vp["annot"][0] == "Promoter-TSS" and arch.vp["gene"][0] == "GENE_A"
    arch.strength(verbose=False)
    assert arch.vp["strength"][0] == pytest.approx(arch.ep["n"][0] + arch.ep["n"][3])
    hubs = arch.prime_hubs(verbose=False)
    assert set(hubs) >= {"hub_uids", "prime_genes"}
    sup = arch.support(genes, r=5000)
    assert sup["GENE_A"][0] == "chr1:900-1100(.)"
    assert arch.support(genes, r=5000, mode="center")["GENE_A"] == sup["GENE_A"]
    assert arch.support(genes, r=5000, rows=True)["GENE_A"].tolist() == [0, 1, 2]


def test_views_and_lookups(arch, cre):
    assert arch.chrom("chr1").n_links == 2 and arch.cis.n_links == 3 and arch.trans.n_links == 1
    assert arch.region("chr1:0-6,000").n_links == 1
    assert arch.neighbor_rows(0)[0].tolist() == [2, 5] if hasattr(arch.neighbor_rows(0), "__getitem__") else True
    assert arch[cre.uid[0]].keys() >= {cre.uid[2], cre.uid[5]}
    assert cre.uid[0] in arch and arch.near("chr1:0-3,000").to_records()[0] == ("chr1", 900, 1100, ".")
    sub = arch.subgraph(rows=[0, 2])
    assert sub.n_links == 1


@pytest.mark.parametrize("backend", installed_backends("graph"))
def test_graph_backends_agree(arch, backend):
    assert arch.components(backend=backend).tolist() == arch.components().tolist()
    assert np.allclose(arch.pagerank("w", backend=backend), arch.pagerank("w"), atol=1e-6)
    g = arch.graph(backend=backend)
    assert g is not None


def test_networkx_merges_parallel_edges_like_scipy():
    pytest.importorskip("networkx")
    L = as_loci([("chr1", 0, 100), ("chr1", 1000, 1100), ("chr1", 2000, 2100)])
    A = Architecture.from_edges(L, [0, 0, 1], [1, 1, 2], w=[1.0, 3.0, 2.0])
    assert np.allclose(A.pagerank("w", backend="networkx"), A.pagerank("w", backend="scipy"))
    assert A.to_networkx()[0][1]["w"] == 4.0 and A.to_scipy()[0, 1] == 4.0


def test_missing_graph_tool_is_an_error_not_a_fallback(arch):
    if gb.backends.installed("graph", "graph-tool"):
        pytest.skip("graph-tool present")
    with pytest.raises(ImportError, match="conda install"):
        arch.components(backend="graph-tool")


def test_exports_and_persistence(arch, tmp_path):
    df = arch.edges_frame()
    assert list(df.columns)[:7] == ["src", "tgt", "uid1", "uid2", "chrom1", "chrom2", "cis"]
    back = Architecture.from_frame(arch.loci, df)
    assert list(back) == list(arch) and back.ep["w"].tolist() == arch.ep["w"].tolist()
    assert Architecture.from_frame(arch.loci, df[["uid1", "uid2", "w"]]).n_links == 4
    with pytest.raises(ValueError, match="edge table"):
        Architecture.from_frame(arch.loci, pd.DataFrame({"x": [1]}))
    M = arch.to_scipy("w")
    assert Architecture.from_scipy(arch.loci, M).n_links == 4
    arch.save(str(tmp_path / "A"))
    back = Architecture.load(str(tmp_path / "A"))
    assert list(back) == list(arch) and back.loci.equals(arch.loci)
    assert arch.shape == (4, len(arch.columns)) and len(arch.head(2)) == 2
    assert arch.describe().loc["edges", "value"] == 4 and "<table" in arch._repr_html_()
    import polars as pl
    assert pl.DataFrame(arch).shape == (4, len(arch.columns))
    assert arch.block_counts().loc["chr1", "chr2"] == 1
    pytest.importorskip("igraph")
    assert arch.to_igraph().ecount() == 4
    pytest.importorskip("anndata")
    assert arch.to_anndata().n_obs == len(arch.loci)


def test_add_mcool_weights(cre, pairs, tmp_path):
    cooler = pytest.importorskip("cooler")
    bins = pd.DataFrame({"chrom": ["chr1"] * 12 + ["chr2"] * 8, "start": list(range(0, 12_000, 1000)) + list(range(0, 8_000, 1000)),
                         "end": list(range(1000, 13_000, 1000)) + list(range(1000, 9_000, 1000))})
    pix = pd.DataFrame({"bin1_id": [0, 1, 0], "bin2_id": [4, 10, 12], "count": [6, 3, 9]})
    cooler.create_cooler(str(tmp_path / "t.cool"), bins, pix)
    A = Architecture.make(cre, pairs, r=100, verbose=False).add_mcool(str(tmp_path / "t.cool"), verbose=False)
    assert A.ep["w"].tolist() == [6.0, 3.0, 0.0, 9.0]                  # cis, cis, chr2 (no pixel), trans
    with pytest.raises(ValueError, match="single-resolution"):
        A.add_mcool(str(tmp_path / "t.cool"), resolution=1000)


def test_draw_layouts(arch, genes):
    arch.annotate(genes, verbose=False).strength(verbose=False)
    for layout in ("spring", "circular", "genomic"):
        ax = arch.draw("chr1:0-12 kb", layout=layout)
        assert ax.get_title().startswith(arch.name)
    ax = arch.draw(("chr1", 0, 12_000), merge_distance=1500, vertex_size_by="strength", vertex_color="annot")
    assert "merged" in ax.get_title()
    with pytest.raises(ValueError, match="no CREs"):
        arch.draw("chr2:7,000-8,000")
    for backend in installed_backends("graph"):
        arch.draw("chr1:0-12kb", backend=backend)
