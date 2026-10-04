import pickle

import pytest

from genomeblocks import Architecture, Genes, Loci
from genomeblocks.locus import Locus


def test_make_from_bedpe(cre, bedpe_path):
    g = Architecture.make(cre, bedpe_path, r=200, verbose=False)
    # anchors at 1000/5000/9000 map to cre[0],cre[1],cre[2]; 3 loops -> 3 edges
    assert g.n_loci == 3
    assert g.n_links == 3


def test_basic_accessors(arch, cre):
    assert arch.n_loci == 4
    assert arch.n_links == 4
    assert cre[0].uid in arch
    neighbors = arch[cre[0].uid]                 # {neighbor_uid: weight}
    assert len(neighbors) == 3
    edge_props = arch[(cre[0].uid, cre[1].uid)]
    assert edge_props["w"] == 5.0


def test_strength(arch, cre):
    arch.strength(key="n", name="strength", verbose=False)
    s = {arch.vp.uid[v]: arch.vp.strength[v] for v in arch.vertices()}
    assert s[cre[0].uid] == pytest.approx(15.0)  # hub: 5+5+5
    assert s[cre[3].uid] == pytest.approx(5.0)
    # no normalized companion property is created
    assert "nstrength" not in arch.vp


def test_elbow_sorts_by_strength(arch, cre):
    # On a 4-node graph the knee detector may return 0; the guarantee that
    # matters is the descending sort (hub first). A real knee is exercised in
    # test_prime_hubs on a larger graph.
    arch.strength(key="n", name="strength", verbose=False)
    cutoff, sorted_uids = arch.elbow("strength", verbose=False)
    assert 0 <= cutoff <= arch.n_loci
    assert sorted_uids[0] == cre[0].uid          # hub ranks first
    strengths = {arch.vp.uid[v]: arch.vp.strength[v] for v in arch.vertices()}
    assert sorted_uids == sorted(sorted_uids, key=lambda u: -strengths[u])


def _star(hub, leaves, n_strong=5.0, n_weak=0.2):
    """Hub connected to every leaf; leaves weakly chained for a strength spread."""
    g = Architecture("star")
    h = g._add_vertex(hub.uid)
    vs = [g._add_vertex(l.uid) for l in leaves]
    for v in vs:
        e = g.add_edge(h, v); g.ep.n[e] = n_strong; g.ep.w[e] = n_strong
    for i in range(0, len(vs) - 1, 2):
        e = g.add_edge(vs[i], vs[i + 1]); g.ep.n[e] = n_weak; g.ep.w[e] = n_weak
    return g


def test_normalize_sets_distance_and_oe(arch, cre):
    arch.normalize(cre, source="w", name="n", verbose=False)
    dists = [arch.ep.d[e] for e in arch.edges()]
    assert all(d > 0 for d in dists)             # all distances populated
    import math
    assert all(math.isfinite(arch.ep.n[e]) for e in arch.edges())


def test_annotate(arch, cre, gtf_path):
    genes = Genes.make(gtf_path, promoter_r=1000)
    arch.annotate(cre, genes, key="n", verbose=False)
    # cre[0] sits in GENE1's promoter window
    v0 = arch.index[cre[0].uid]
    assert "Promoter" in arch.vp.annot[v0]
    assert arch.vp.gene[v0] == "GENE1"
    # a non-promoter neighbour inherits the gene of its top promoter contact
    v1 = arch.index[cre[1].uid]
    assert arch.vp.gene[v1] == "GENE1"


def test_prime_hubs(gtf_path):
    hub = Locus("chr1", 1000, 1100)               # in GENE1's promoter window
    leaves = [Locus("chr1", 50000 + i * 500, 50000 + i * 500 + 100) for i in range(60)]
    cre = Loci([hub] + leaves)
    g = _star(hub, leaves)

    genes = Genes.make(gtf_path, promoter_r=1000)
    g.annotate(cre, genes, key="n", verbose=False)
    res = g.prime_hubs(key="n", verbose=False)

    assert {"prime_genes", "promoter_genes", "enhancer_genes",
            "hub_uids", "cutoff"}.issubset(res)
    assert res["cutoff"] >= 1
    assert hub.uid in res["hub_uids"]
    assert "GENE1" in res["prime_genes"]


def test_subgraph_by_uids(arch, cre):
    sub = arch.subgraph(uids=[cre[0].uid, cre[1].uid])
    assert sub.n_loci == 2
    assert sub.n_links == 1                       # only the 0-1 edge survives


def test_copy_is_deep(arch):
    c = arch.copy()
    assert c.n_loci == arch.n_loci and c.n_links == arch.n_links
    # mutating the copy's weights doesn't touch the original
    for e in c.edges():
        c.ep.w[e] = 99.0
    assert any(arch.ep.w[e] != 99.0 for e in arch.edges())


def test_copy_preserves_property_value_types(arch):
    # Regression: type was inferred from a sampled Python value, but graph-tool
    # returns bools as ints — so bool props were silently copied as "double".
    arch.vp["is_hub"] = arch.new_vertex_property("bool")
    arch.vp["rank"] = arch.new_vertex_property("int32_t")
    for v in arch.vertices():
        arch.vp["is_hub"][v] = True
        arch.vp["rank"][v] = 3
    c = arch.copy()
    assert c.vp["is_hub"].value_type() == "bool"
    assert c.vp["rank"].value_type() == "int32_t"
    assert bool(c.vp["is_hub"][c.vertex(0)]) is True
    assert int(c.vp["rank"][c.vertex(0)]) == 3


def test_set_operations(arch):
    assert (arch | arch).n_loci == arch.n_loci
    assert (arch & arch).n_loci == arch.n_loci


def test_pickle_roundtrip(arch, cre):
    data = pickle.dumps(arch)
    g2 = pickle.loads(data)
    assert g2.n_loci == arch.n_loci
    assert g2.n_links == arch.n_links
    assert cre[0].uid in g2
    e = g2.edge(g2.index[cre[0].uid], g2.index[cre[1].uid])
    assert g2.ep.w[e] == 5.0


def test_pickle_preserves_bool_property_type(arch):
    # Regression: __getstate__ read bool maps as ints, so __setstate__ rebuilt
    # them as float. Value types are now recorded and restored.
    arch.vp["is_hub"] = arch.new_vertex_property("bool")
    for v in arch.vertices():
        arch.vp["is_hub"][v] = True
    g2 = pickle.loads(pickle.dumps(arch))
    assert g2.vp["is_hub"].value_type() == "bool"
    assert bool(g2.vp["is_hub"][g2.vertex(0)]) is True


def test_normalize_handles_trans_edges():
    # Regression: a single inter-chromosomal edge made normalize raise
    # TypeError (no distance between chromosomes).
    import math
    cis = [Locus("chr1", 1000 + i * 3000, 1100 + i * 3000) for i in range(12)]
    far = [Locus("chr2", 5000, 5100), Locus("chr3", 7000, 7100)]
    loci = Loci(cis + far)
    g = Architecture("cis+trans")
    for i in range(len(cis) - 1):
        for j in range(i + 1, min(i + 4, len(cis))):
            e = g.add_edge(g._add_vertex(cis[i].uid), g._add_vertex(cis[j].uid))
            g.ep.w[e] = 30.0 / (j - i)
    trans_w = {(0, far[0]): 2.0, (3, far[1]): 6.0}
    for (i, other), w in trans_w.items():
        e = g.add_edge(g._add_vertex(cis[i].uid), g._add_vertex(other.uid))
        g.ep.w[e] = w

    g.normalize(loci, source="w", name="n", verbose=False)
    for e in g.edges():
        a, b = g.vp.uid[e.source()], g.vp.uid[e.target()]
        if a.split(":")[0] == b.split(":")[0]:
            assert math.isfinite(g.ep.d[e]) and g.ep.d[e] > 0
            assert math.isfinite(g.ep.n[e])
        else:                                       # trans: flat expectation
            assert math.isinf(g.ep.d[e])
            assert g.ep.n[e] == pytest.approx(g.ep.w[e] / 4.0)   # mean trans w = 4
    g.prune(verbose=False)                          # trans edges survive prune
    assert g.n_links == 30 + 2


def test_add_mcool_splits_pixel_counts(tmp_path):
    cooler = pytest.importorskip("cooler")
    import numpy as np
    import pandas as pd
    bins = cooler.binnify(pd.Series({"chr1": 10_000, "chr2": 5_000}), 1000)
    pixels = pd.DataFrame({"bin1_id": [0, 0, 2, 3], "bin2_id": [2, 5, 2, 12],
                           "count": [12, 7, 4, 9]})
    uri = str(tmp_path / "t.cool")
    cooler.create_cooler(uri, bins, pixels)

    loci = Loci([Locus("chr1", 100, 200), Locus("chr1", 600, 700),      # bin 0
                 Locus("chr1", 2100, 2200), Locus("chr1", 2500, 2600),  # bin 2
                 Locus("chr1", 5100, 5200), Locus("chr1", 8100, 8200), # bins 5, 8
                 Locus("chr2", 2100, 2200)])                           # bin 12
    g = Architecture("mcool")
    def link(i, j):
        return g.add_edge(g._add_vertex(loci[i].uid), g._add_vertex(loci[j].uid))
    e02, e12, e03 = link(0, 2), link(1, 2), link(0, 3)   # three edges in pixel (0, 2)
    e23 = link(2, 3)                                    # pixel (2, 2)
    e04 = link(0, 4)                                    # pixel (0, 5)
    e05 = link(0, 5)                                    # bin pair (0, 8): no pixel
    e36 = link(3, 6)                                    # trans pixel (2, 12)
    g.add_mcool(loci, uri, name="w", verbose=False)
    assert [g.ep.w[e] for e in (e02, e12, e03)] == pytest.approx([4.0, 4.0, 4.0])
    assert g.ep.w[e23] == pytest.approx(4.0)
    assert g.ep.w[e04] == pytest.approx(7.0)
    assert g.ep.w[e05] == 0.0
    assert g.ep.w[e36] == 0.0                           # pixel (3, 12) is bin 3, not 2
