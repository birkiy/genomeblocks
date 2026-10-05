"""Motif libraries, scanning engines, the FASTA backends and motif statistics."""
import numpy as np
import pytest

import genomeblocks as gb
from genomeblocks import as_loci, motifs as gm
from genomeblocks.backends.fasta import open_fasta, read_fasta
from genomeblocks.backends.motifs import Library, load_motifs

from conftest import installed_backends

MOTIF_BACKENDS = installed_backends("motifs")
FASTA_BACKENDS = installed_backends("fasta")


# ── libraries ─────────────────────────────────────────────────────────────────

def test_load_formats(jaspar_path, tmp_path):
    lib = load_motifs(jaspar_path)
    assert lib.names == ["M1", "M2"] and lib.descriptions == ["TFA", "TFB"] and lib.widths.tolist() == [4, 4]
    assert lib.consensus(0) == "ACGT" and lib["M2"].shape == (4, 4) and lib[0] is lib.counts[0]
    assert lib.select("tfa").names == ["M1"] and lib.select(["M2"], match="exact").names == ["M2"]
    lib.to_jaspar(str(tmp_path / "j16.txt"))
    assert load_motifs(str(tmp_path / "j16.txt"), format="jaspar16").names == ["M1", "M2"]
    lib.to_meme(str(tmp_path / "m.meme"))
    meme = load_motifs(str(tmp_path / "m.meme"), format="meme")
    assert meme.names == ["M1", "M2"] and np.allclose(meme.pfm(0), lib.pfm(0), atol=0.01)
    (tmp_path / "t.transfac").write_text("ID M1\nBF M1_TF\nP0  A  C  G  T\n01 1 1 17 1 G\n02 1 1 17 1 G\n03 1 17 1 1 C\nXX\n//\n")
    assert load_motifs(str(tmp_path / "t.transfac"), format="transfac").names == ["M1"]
    (tmp_path / "u.txt").write_text("U1\nA:\t0.1\t0.7\t0.1\nC:\t0.1\t0.1\t0.7\nG:\t0.7\t0.1\t0.1\nT:\t0.1\t0.1\t0.1\n\n")
    assert load_motifs(str(tmp_path / "u.txt"), format="uniprobe").consensus(0) == "GAC"
    assert lib.consensus(1) == "GGAA"
    assert load_motifs({"X": np.eye(4)[[0, 1, 2, 3]]}).counts[0].sum() == pytest.approx(400)   # probabilities x 100
    assert load_motifs([np.ones((5, 4)), np.ones((6, 4))]).names == ["motif1", "motif2"]
    assert load_motifs(np.ones((5, 4))).widths.tolist() == [5]
    bio = pytest.importorskip("Bio.motifs")
    assert load_motifs(lib.to_biopython()).names == ["M1", "M2"]
    with pytest.raises(ValueError, match="unknown motif format"):
        load_motifs(jaspar_path, format="homer")
    with pytest.raises(FileNotFoundError):
        load_motifs("nope.jaspar")
    (tmp_path / "freq.jaspar").write_text(">M1 F\n0.5 0.5\n0.5 0\n0 0.5\n0 0\n")   # non-integer counts load
    assert load_motifs(str(tmp_path / "freq.jaspar")).counts[0].tolist() == [[0.5, 0.5, 0, 0], [0.5, 0, 0.5, 0]]
    (tmp_path / "bad.jaspar").write_text(">M1\n1 2 3\n4 5 6\n")
    with pytest.raises(ValueError, match="rows"):
        load_motifs(str(tmp_path / "bad.jaspar"))


def test_library_as_a_table(jaspar_path):
    lib = load_motifs(jaspar_path)
    assert lib.shape == (2, 4) and lib.columns == ["name", "description", "width", "consensus"]
    assert list(lib) == ["M1", "M2"] and len(lib.head(1)) == 1 and lib.describe().loc["motifs", "value"] == 2
    import polars as pl
    assert pl.DataFrame(lib).shape == (2, 4) and "<table" in lib._repr_html_()
    with pytest.raises(KeyError):
        lib["nope"]


# ── scanning ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("backend", MOTIF_BACKENDS)
def test_engines_find_the_planted_sites(fasta_path, jaspar_path, backend):
    L = as_loci([("chr1", 950, 1050), ("chr1", 4950, 5050), ("chr1", 7000, 7100)])
    M = gm.scan_motifs_matrix(L, fasta_path, jaspar_path, r=50, threshold=7.0, norm=False, backend=backend,
                              verbose=False)
    assert list(M.columns) == ["M1", "M2"] and list(M.index) == L.uid.tolist()
    assert M.iloc[0, 0] >= 3 and M.iloc[1, 1] >= 3                   # ACGT x3 at 1000, GGAA x3 at 5000
    ref = gm.scan_motifs_matrix(L, fasta_path, jaspar_path, r=50, threshold=7.0, norm=False, verbose=False)
    assert M.equals(ref)                                             # every engine, same hits
    M2 = gm.scan_motifs_matrix(L, fasta_path, jaspar_path, r=50, pvalue=0.01, norm=False, backend=backend,
                               both_strands=True, verbose=False)              # 4-mers: p >= 1/256 per site
    assert M2.shape == (3, 2) and (M2.to_numpy() >= M.to_numpy()).all()


def test_matrix_masked_profile_and_totals(fasta_path, jaspar_path, cre):
    L = as_loci([("chr1", 950, 1050), ("chr1", 4950, 5050)])
    masked = gm.scan_motifs_matrix_masked(L, fasta_path, jaspar_path, ["M1"], r=50, threshold=[7.0, 7.0], verbose=False)
    assert list(masked.columns) == ["M2"]
    masked2 = gm.scan_motifs_matrix_masked(L, fasta_path, jaspar_path, ["M1"], r=50, pvalue=1e-3, skip_anchors=False,
                                           verbose=False)
    assert list(masked2.columns) == ["M1", "M2"]
    prof, names = gm.scan_motifs_profile(L, fasta_path, jaspar_path, ["M1"], r=60, n_bins=12, threshold=7.0,
                                         smooth=0, verbose=False)
    assert prof.shape == (2, 1, 12) and names == ["M1"] and prof[0].sum() >= 3
    tot = gm.scan_motifs(L, fasta_path, jaspar_path, r=50, threshold=7.0, norm=False, verbose=False)
    assert tot["M1"] >= 3
    assert cre.scan_motifs_matrix(fasta_path, jaspar_path, r=50, threshold=7.0, verbose=False).shape == (7, 2)


def test_bootstrap_and_archetypes(fasta_path, jaspar_path, cre):
    L = as_loci([("chr1", s, s + 100) for s in range(500, 6500, 200)])
    M = gm.scan_motifs_matrix(L, fasta_path, jaspar_path, r=50, threshold=7.0, norm=False, verbose=False)
    res = gm.bootstrap_enrichment({"a": M.iloc[:10], "b": M.iloc[10:]}, M, boot=20, sample=5, seed=0, verbose=False)
    assert res is not None
    D, names, pwms = gm.pwm_distance_matrix(jaspar_path, min_overlap=3, verbose=False)
    assert D.shape == (2, 2) and D[0, 0] == 0 and names == ["M1", "M2"]
    Z = gm.cluster_motifs(D, cutoff=0.5)
    assert Z is not None
    pfm = gm.archetype([load_motifs(jaspar_path).pfm(0), load_motifs(jaspar_path).pfm(1)], min_overlap=3)
    assert pfm.shape[0] == 4


# ── FASTA backends ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("backend", FASTA_BACKENDS)
def test_fasta_backends_agree(fasta_path, backend):
    seqs = read_fasta(fasta_path)
    src = open_fasta(fasta_path, backend=backend) if backend != "memory" else open_fasta(seqs)
    assert src.sizes() == {"chr1": 20_000, "chr2": 8_000}
    want = [seqs["chr1"][1000:1012], seqs["chr2"][100:160], seqs["chr1"][19_990:20_000]]
    got = src.fetch_many(np.array(["chr1", "chr2", "chr1"], dtype=object), np.array([1000, 100, 19_990]),
                         np.array([1012, 160, 20_000]))
    assert list(got) == want
    assert gb.Genome.from_fasta(fasta_path).sizes == {"chr1": 20_000, "chr2": 8_000}


def test_uneven_fasta_lines(tmp_path):
    p = tmp_path / "u.fa"
    p.write_text(">c1\nACGTACGT\nACG\nACGTACGT\n")
    with pytest.raises(ValueError, match="uneven"):
        open_fasta(str(p))
    assert open_fasta(str(p), backend="memory").sizes() == {"c1": 19}


# ── no lightmotif needed: parsers, log-odds and p-value cutoffs of our own ───

_FORMATS = {
    "jaspar": ">MA0001.1 AGL3\n0 3 79 40\n94 75 4 3\n1 0 3 4\n2 19 11 50\n"
              ">MA0002.2\tRUNX1 extra words\n10 20\n30 40\n50 60\n70 80\n",
    "jaspar16": ">MA0001.1 AGL3\nA  [ 0  3 79 40 ]\nC  [94 75  4  3 ]\nG  [ 1  0  3  4 ]\nT  [ 2 19 11 50 ]\n",
    "transfac": "AC  M00001\nXX\nID  V$MYOD_01\nXX\nNA  MyoD\nXX\nDE  myoblast determination\nXX\n"
                "P0      A      C      G      T\n01      1      2      2      0      S\n"
                "02      2      1      2      0      R\n03      3      0      1      1      A\nXX\n//\n"
                "AC  M00002\nXX\nID  V$E47_01\nXX\nP0 A C G T\n01 4 4 3 1 N\n02 2 5 4 1 N\nXX\n//\n",
    "uniprobe": "Foxa2 primary\nA:\t0.25\t0.1\t0.7\nC:\t0.25\t0.2\t0.1\nG:\t0.25\t0.3\t0.1\nT:\t0.25\t0.4\t0.1\n",
}


@pytest.mark.parametrize("fmt", sorted(_FORMATS))
def test_parsers_match_lightmotif(tmp_path, fmt):
    """Our pure-Python readers give lightmotif's names, descriptions, counts and
    (bit for bit) log-odds, so dropping lightmotif changes no result."""
    lightmotif = pytest.importorskip("lightmotif")
    path = tmp_path / f"m.{fmt}"
    path.write_text(_FORMATS[fmt])
    ours = load_motifs(str(path), format=fmt)
    ref = list(lightmotif.load(str(path), format=fmt))
    assert len(ours) == len(ref)
    for i, m in enumerate(ref):
        name = next(str(getattr(m, a)) for a in ("name", "id", "accession") if getattr(m, a, None))
        assert ours.names[i] == name and ours.descriptions[i] == (getattr(m, "description", None) or "")
        if getattr(m, "counts", None) is not None:
            acgt = [0, 1, 3, 2]
            assert np.array_equal(ours.counts[i], np.asarray(m.counts)[:, acgt])
            p = m.counts.normalize(0.1).log_odds()
            lm = np.array([[p[r][j] for j in acgt] for r in range(len(p))], np.float64)
            assert np.array_equal(ours.logodds(i), lm)                 # bit for bit


def test_pvalue_threshold_is_exact():
    """Against brute force over every 6-mer: P(score >= t) <= p, and one grid
    step lower it is > p."""
    from itertools import product
    from genomeblocks.backends.motifs import logodds_matrix, threshold_from_pvalue
    rng = np.random.default_rng(1)
    for _ in range(5):
        lo = logodds_matrix(rng.integers(0, 50, (6, 4)))
        scores = np.array([sum(lo[j, b] for j, b in enumerate(seq)) for seq in product(range(4), repeat=6)])
        for p in (1e-1, 1e-2, 1e-3):
            t = threshold_from_pvalue(lo, p)
            assert np.mean(scores >= t - 6 * 5e-4) > p or t >= scores.max() - 1e-9
            assert np.mean(scores >= t + 6 * 5e-4) <= p


def test_pvalue_threshold_matches_moods():
    MOODS_tools = pytest.importorskip("MOODS.tools")
    from genomeblocks.backends.motifs import logodds_matrix, threshold_from_pvalue
    rng = np.random.default_rng(2)
    for w in (8, 12, 20):
        lo = logodds_matrix(rng.integers(0, 40, (w, 4)))
        for p in (1e-3, 1e-4, 1e-5):
            ours = threshold_from_pvalue(lo, p)
            if p < 4.0 ** -w:                       # unreachable: we keep only perfect matches
                assert ours == pytest.approx(lo.max(1).sum())
                continue
            theirs = MOODS_tools.threshold_from_p(lo.T.tolist(), [0.25] * 4, p)
            assert abs(ours - theirs) < 0.01


def test_moods_is_the_default_engine():
    pytest.importorskip("MOODS.scan")
    assert gb.backends.resolve("motifs") == "moods"


def test_motifs_work_without_lightmotif(monkeypatch, fasta_path, jaspar_path):
    """The conda package has no lightmotif: loading, p-value cutoffs and scanning
    must not import it."""
    pytest.importorskip("MOODS.scan")
    import sys
    from genomeblocks import backends
    monkeypatch.setitem(sys.modules, "lightmotif", None)             # `import lightmotif` now fails
    monkeypatch.setattr(backends, "_installed_cache", {})
    assert not backends.installed("motifs", "lightmotif")
    assert backends.resolve("motifs") == "moods"
    L = as_loci([("chr1", 950, 1050), ("chr1", 4950, 5050)])
    M = gm.scan_motifs_matrix(L, fasta_path, jaspar_path, r=50, pvalue=0.01, norm=False, verbose=False)
    assert M.iloc[0, 0] >= 3 and M.iloc[1, 1] >= 3
    with pytest.raises(ImportError, match="lightmotif"):
        backends.resolve("motifs", "lightmotif")
