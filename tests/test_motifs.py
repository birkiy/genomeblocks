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
    (tmp_path / "bad.jaspar").write_text(">M1\n0.5 0.5 0 0\n")
    with pytest.raises(ValueError, match="integers"):
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
