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
