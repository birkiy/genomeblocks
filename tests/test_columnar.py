"""genomeblocks.columnar (prototype): behaviour + parity with the classic objects.

A small random "world" (3 chromosomes, 1.5k CREs, 600 loops incl. trans, 120
genes, a 5 kb cooler) is generated once; every columnar result is checked
against the classic implementation on the same files.
"""
import contextlib
import io

import numpy as np
import pytest

import genomeblocks.columnar as gbc
from genomeblocks import Genes as ClassicGenes
from genomeblocks import Loci as ClassicLoci
from genomeblocks.locus import Locus


def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return f()


SIZES = {"chr1": 2_000_000, "chr2": 1_500_000, "chr10": 1_000_000}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    rng = np.random.default_rng(0)
    d = tmp_path_factory.mktemp("columnar")
    cre = []
    for c, n in SIZES.items():
        st = np.sort(rng.choice(n // 1000 - 2, 500, replace=False)) * 1000 + rng.integers(0, 400, 500)
        cre += [(c, int(s), int(s + rng.integers(150, 600))) for s in st]
    order = rng.permutation(len(cre))                      # unsorted on disk
    (d / "peaks.bed").write_text("".join(f"{c}\t{s}\t{e}\tp{i}\t0\t.\n"
                                         for i, (c, s, e) in enumerate(cre[k] for k in order)))
    # loops: cis between nearby CREs, plus some trans
    mids = {c: [(s + e) // 2 for cc, s, e in cre if cc == c] for c in SIZES}
    loops = []
    for _ in range(570):
        c = rng.choice(list(SIZES))
        i = int(rng.integers(0, len(mids[c]) - 6))
        j = i + int(rng.integers(1, 6))
        loops.append((c, mids[c][i], c, mids[c][j]))
    for _ in range(30):
        c1, c2 = rng.choice(list(SIZES), 2, replace=False)
        loops.append((c1, rng.choice(mids[c1]), c2, rng.choice(mids[c2])))
    (d / "loops.bedpe").write_text("".join(
        f"{a}\t{x - 500}\t{x + 500}\t{b}\t{y - 500}\t{y + 500}\tL{k}\t1\n"
        for k, (a, x, b, y) in enumerate(loops)))
    (d / "loops_cis.bedpe").write_text("".join(
        f"{a}\t{x - 500}\t{x + 500}\t{b}\t{y - 500}\t{y + 500}\tL{k}\t1\n"
        for k, (a, x, b, y) in enumerate(loops) if a == b))
    # GTF: genes with 1-3 transcripts, exons, CDS and UTRs
    lines = []
    for c, n in SIZES.items():
        for g in range(40):
            gs = int(rng.integers(10_000, n - 60_000))
            ge = gs + int(rng.integers(5_000, 50_000))
            sd = "+" if rng.random() < 0.5 else "-"
            gid, name = f"{c}_G{g}", f"GENE_{c}_{g}"
            a = f'gene_id "{gid}"; gene_name "{name}"; gene_type "protein_coding";'
            lines.append(f"{c}\tS\tgene\t{gs}\t{ge}\t.\t{sd}\t.\t{a}")
            for t in range(int(rng.integers(1, 4))):
                ts = gs + int(rng.integers(0, 2000)) * t
                tid = f"{gid}.T{t}"
                ta = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{name}";'
                lines.append(f"{c}\tS\ttranscript\t{ts}\t{ge}\t.\t{sd}\t.\t{ta}")
                ex = np.sort(rng.choice(np.arange(ts, ge - 300, 300), 4, replace=False))
                for k, x in enumerate(ex):
                    ea = ta + f' exon_number "{k + 1}";'
                    lines.append(f"{c}\tS\texon\t{x}\t{x + 200}\t.\t{sd}\t.\t{ea}")
                    kind = "five_prime_UTR" if k == 0 else "three_prime_UTR" if k == 3 else "CDS"
                    lines.append(f"{c}\tS\t{kind}\t{x}\t{x + 200}\t.\t{sd}\t.\t{ea}")
    (d / "genes.gtf").write_text("\n".join(lines) + "\n")
    # cooler: background near the diagonal + counts at every loop pixel
    cooler = pytest.importorskip("cooler")
    import pandas as pd
    bins = pd.concat([pd.DataFrame({"chrom": c, "start": np.arange(0, n, 5000),
                                    "end": np.minimum(np.arange(0, n, 5000) + 5000, n)})
                      for c, n in SIZES.items()], ignore_index=True)
    off = {c: int(np.flatnonzero(bins.chrom == c)[0]) for c in SIZES}
    px = {}
    for a, x, b, y in loops:
        i, j = sorted((off[a] + x // 5000, off[b] + y // 5000))
        px[(i, j)] = px.get((i, j), 0) + int(rng.integers(1, 30))
    for i in rng.integers(0, len(bins) - 5, 3000):
        px[(int(i), int(i) + int(rng.integers(0, 4)))] = int(rng.integers(1, 9))
    k = sorted(px)
    pixels = pd.DataFrame({"bin1_id": [a for a, _ in k], "bin2_id": [b for _, b in k],
                           "count": [px[x] for x in k]})
    cooler.create_cooler(str(d / "hic.cool"), bins, pixels, ordered=True)
    return d


@pytest.fixture
def L(world):
    return gbc.Loci.make(str(world / "peaks.bed"))


@pytest.fixture
def classic(world):
    return ClassicLoci.make(str(world / "peaks.bed"))


def as_set(loci):
    return sorted((l.chrom, l.start, l.end) for l in loci)


# ── Genome / Loci ─────────────────────────────────────────────────────────────

def test_shared_genome_codes(world, L):
    G = gbc.Genes.make(str(world / "genes.gtf"))
    assert L.genome is G.genes.genome
    c = L.genome.code["chr2"]
    assert set(L.chroms[L.codes == c]) == {"chr2"}
    assert set(G.genes.chroms[G.genes.codes == c]) == {"chr2"}


def test_make_sorts_into_blocks(L, classic):
    assert len(L) == len(classic) and L.is_sorted
    assert list(L.chrom_offsets) == ["chr1", "chr2", "chr10"]       # natural order
    a, b = L.chrom_offsets["chr2"]
    assert set(L.chroms[a:b]) == {"chr2"}
    assert np.shares_memory(L.by_chrom("chr2").starts, L.starts)     # a view, not a copy
    assert as_set(L) == as_set(classic)


def test_row_access_and_write_through(L):
    v = L[3]
    assert isinstance(v, Locus) and v.row == 3
    assert L[v.uid].row == 3
    v.end = v.end + 7
    assert L.ends[3] == v.end
    assert L[-1].row == len(L) - 1


@pytest.mark.parametrize("op", ["intersect", "difference"])
def test_set_ops_match_classic(world, L, classic, op):
    other = ClassicLoci(classic.slop(2000)[::3])        # classic slicing returns a list
    mine = getattr(L, op)(gbc.Loci.from_loci(other))
    assert as_set(mine) == as_set(getattr(classic, op)(other))


def test_merge_and_slop_match_classic(L, classic):
    assert as_set(L.slop(3000).merge()) == as_set(classic.slop(3000).merge())


def test_overlaps_region(L, classic):
    got = L.overlaps("chr1", 100_000, 400_000)
    want = classic.overlaps("chr1", 100_000, 400_000)
    assert as_set(got) == as_set(want)


def test_exports_roundtrip(tmp_path, L):
    df = L.to_pandas()
    assert list(df.columns[:4]) == ["chrom", "start", "end", "strand"]
    back = gbc.Loci.from_frame(df, "chrom", "start", "end", "strand")
    assert np.array_equal(back.starts, L.starts) and np.array_equal(back.codes, L.codes)
    p = tmp_path / "l.parquet"
    L.save(str(p))
    R = gbc.Loci.load(str(p))
    assert np.array_equal(R.uid, L.uid)


def test_classic_methods_still_work(L):
    uids = L.uid[:5].tolist()
    sub = L.subloci(uids)                           # a classic Loci method, via the bridge
    assert [l.uid for l in sub] == uids


# ── Genes ─────────────────────────────────────────────────────────────────────

@pytest.fixture
def genes_pair(world):
    path = str(world / "genes.gtf")
    return gbc.Genes.make(path), quiet(lambda: ClassicGenes.make(path))


def test_genes_tables(genes_pair):
    G, C = genes_pair
    assert len(G.genes) == len(C)
    assert len(G.transcripts) == sum(len(g.transcripts) for g in C.values())
    g = G["GENE_chr2_5"]
    c = next(x for x in C.values() if x.gene_name == "GENE_chr2_5")
    assert (g.chrom, g.start, g.end, g.strand) == (c.chrom, c.start, c.end, c.strand)
    assert len(g.transcripts) == len(c.transcripts)
    assert (g.tss.start, g.tss.end) == (c.tss.start, c.tss.end)


@pytest.mark.parametrize("key", ["prom", "exon", "utr5", "utr3"])
def test_annotation_index_matches_classic(genes_pair, key):
    G, C = genes_pair
    assert as_set(G.annot[key]) == as_set(C.annot[key])


def test_annotations_and_nearest_match_classic(genes_pair, L, classic):
    G, C = genes_pair
    a = dict(zip(*G.annotations(L).T.values))
    b = dict(zip(*C.annotations(classic).T.values))
    assert a == b
    n = G.nearest_genes(L)
    m = C.nearest_genes(classic)
    prom = {u for u, x in a.items() if x == "Promoter-TSS"}
    assert prom
    mine, theirs = dict(zip(n.Name, n.Name_b)), dict(zip(m.Name, m.Name_b))
    assert all(mine[u] == theirs[u] for u in prom)


# ── Architecture ─────────────────────────────────────────────────────────────

@pytest.fixture
def arch_pair(world, L, classic):
    pytest.importorskip("graph_tool")
    from genomeblocks import Architecture as Classic
    bedpe, cool = str(world / "loops_cis.bedpe"), str(world / "hic.cool")
    G, C = gbc.Genes.make(str(world / "genes.gtf")), quiet(lambda: ClassicGenes.make(str(world / "genes.gtf")))
    A = gbc.Architecture.make(L, bedpe, verbose=False)
    A.add_mcool(cool, verbose=False).normalize(verbose=False).annotate(G, verbose=False)
    A.strength(verbose=False)
    O = quiet(lambda: Classic.make(classic, bedpe, verbose=False))
    quiet(lambda: O.add_mcool(classic, cool, verbose=False))
    quiet(lambda: O.normalize(classic, verbose=False))
    quiet(lambda: O.annotate(classic, C, verbose=False))
    quiet(lambda: O.strength(verbose=False))
    return A, O


def _edge_props(O, prop):
    uid = np.array(list(O.vp.uid), dtype=object)
    return {tuple(sorted((uid[int(a)], uid[int(b)]))): v
            for a, b, v in O.get_edges(eprops=[O.ep[prop]])}


def test_architecture_matches_classic(arch_pair):
    A, O = arch_pair
    L = A.loci
    keys = [tuple(sorted(k)) for k in zip(L.uid[A.src], L.uid[A.tgt])]
    assert A.n_links == O.n_links and A.n_loci == O.n_loci
    for prop in ("w", "n", "d"):
        ref = _edge_props(O, prop)
        assert set(ref) == set(keys)
        assert np.allclose([ref[k] for k in keys], A.ep[prop], rtol=1e-6), prop
    assert (A.ep.w > 0).sum() > 100
    rows = np.array([L.uids[u] for u in O.vp.uid])
    assert np.array_equal(np.array(list(O.vp.annot), object), A.vp.annot[rows])
    assert np.allclose(np.asarray(O.vp.strength.a), A.vp.strength[rows])
    gene = np.array(list(O.vp.gene), object)
    assert (gene == A.vp.gene[rows]).mean() > 0.98          # rest: tied top weights


def test_prime_hubs_match_classic(arch_pair):
    A, O = arch_pair
    a, b = A.prime_hubs(verbose=False), quiet(lambda: O.prime_hubs(verbose=False))
    assert a["cutoff"] == b["cutoff"] and set(a["hub_uids"]) == set(b["hub_uids"])


@pytest.fixture
def mixed(world, L):
    G = gbc.Genes.make(str(world / "genes.gtf"))
    A = gbc.Architecture.make(L, str(world / "loops.bedpe"), verbose=False)
    return A.add_mcool(str(world / "hic.cool"), verbose=False).normalize(verbose=False) \
            .annotate(G, verbose=False).strength(verbose=False)


def test_blocks_are_contiguous_and_trans_last(mixed):
    A = mixed
    b = A.blocks
    assert list(b)[:3] == ["chr1", "chr2", "chr10"] and list(b)[-1] == "trans"
    for c, (lo, hi) in b.items():
        cs, ct = A.loci.chroms[A.src[lo:hi]], A.loci.chroms[A.tgt[lo:hi]]
        if c == "trans":
            assert hi - lo > 0 and np.all(cs != ct)
        else:
            assert set(cs) == {c} and set(ct) == {c}
    assert np.all(A.src < A.tgt)


def test_views_share_memory(mixed):
    A = mixed
    v = A.chrom("chr2")
    assert np.shares_memory(v.src, A.src) and np.shares_memory(v.ep.n, A.ep.n)
    assert v.vp is A.vp
    assert A.cis.n_links + A.trans.n_links == A.n_links
    assert sum(v.n_links for _, v in A.chroms()) == A.cis.n_links


def test_trans_edges_normalize_and_stay_in_neighbours(mixed):
    A = mixed
    t = A.trans
    assert np.all(np.isinf(t.ep.d)) and np.isfinite(t.ep.n).all() and (t.ep.n > 0).any()
    r = int(t.src[0])
    nb = A.neighbors(r)
    assert (~nb.cis).any() and nb.cis.any() | (~nb.cis).all()
    assert int(t.tgt[0]) in set(nb.row)


def test_graph_sees_trans_edges(mixed):
    pytest.importorskip("graph_tool")
    A = mixed
    g = A.graph()
    assert g.num_edges() == A.n_links and g.num_vertices() == len(A.loci)
    both = A.components()
    cis_only = A.cis.components(name="cc")
    n = lambda c: len(np.unique(c[c >= 0]))
    assert n(both) < n(cis_only)                    # trans edges join chromosomes


def test_region_subgraph_and_set_ops(mixed):
    A = mixed
    r = A.region("chr1", 0, 1_000_000)
    assert r.n_links > 0 and set(A.loci.chroms[r.src]) == {"chr1"}
    assert np.all(A.loci.starts[r.tgt] < 1_000_000)
    out = A.region("chr1", 0, 1_000_000, both=False)
    assert out.n_links >= r.n_links
    prom = A.subgraph(vp="annot", values="Promoter-TSS")
    assert np.all(A.vp.annot[prom.src] == "Promoter-TSS")
    u = A.cis | A.trans
    assert u.n_links == A.n_links
    assert (A & A.trans).n_links == A.trans.n_links


def test_save_load_roundtrip(tmp_path, mixed):
    A = mixed
    A.save(str(tmp_path / "arch"))
    B = gbc.Architecture.load(str(tmp_path / "arch"))
    assert np.array_equal(B.src, A.src) and np.array_equal(B.tgt, A.tgt)
    assert np.array_equal(B.loci.uid, A.loci.uid)
    for k in A.ep:
        assert np.array_equal(B.ep[k], A.ep[k])
    assert np.array_equal(B.vp.gene, A.vp.gene) and np.allclose(B.vp.strength, A.vp.strength)
    assert B.blocks == A.blocks


def test_classic_lookups_and_legacy_export(mixed):
    pytest.importorskip("graph_tool")
    A = mixed
    i, j = int(A.src[0]), int(A.tgt[0])
    ui, uj = A.loci.uid[i], A.loci.uid[j]
    assert uj in A[ui]
    assert A[ui, uj]["w"] == A.ep.w[0]
    O = A.to_legacy()
    assert O.n_links == A.n_links and O.n_loci == A.n_loci
    assert set(O.index) == set(A.loci.uid[A.degree > 0])


def test_genes_save_load(tmp_path, genes_pair, L):
    G, _ = genes_pair
    G.save(str(tmp_path / "genes"))
    H = gbc.Genes.load(str(tmp_path / "genes"))
    assert np.array_equal(H.labels(L), G.labels(L))
    assert H["GENE_chr1_3"].gene_id == G["GENE_chr1_3"].gene_id


def test_igv_html_embeds_tracks(tmp_path, mixed, world):
    import base64
    import gzip
    import json
    import re
    from genomeblocks.columnar.genes import Genes
    from genomeblocks.columnar.igv import igv_html
    A = mixed
    G = Genes.make(str(world / "genes.gtf"))
    hub = int(np.argmax(A.vp.strength))
    region = f"{A.loci.chroms[hub]}:{A.loci.starts[hub]}-{A.loci.ends[hub]}"
    out = tmp_path / "share.html"
    sizes = igv_html(str(out), regions=[region], loci={"CREs": A.loci}, genes=G,
                     architecture=A, title="test")
    page = out.read_text()
    cfg = json.loads(re.search(r"const CONFIG = (\{.*?\});\n", page, re.S).group(1))
    names = [t["name"] for t in cfg["tracks"]]
    assert names == ["CREs", "loops (n)", "genes"] and sizes["total"] == len(page)
    text = lambda t: gzip.decompress(base64.b64decode(t["url"].split(",", 1)[1])).decode()
    assert len(text(cfg["tracks"][0]).splitlines()) == len(A.loci)
    loops = text(cfg["tracks"][1]).splitlines()
    assert loops and all(len(l.split("\t")) == 8 for l in loops)
    genes = text(cfg["tracks"][2]).splitlines()
    assert len(genes) == len(G.genes) and all(len(l.split("\t")) == 12 for l in genes)
    assert cfg["reference"]["format"] == "chromsizes" and "const SIZES = \"chr1" in page
