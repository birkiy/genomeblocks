"""End-to-end functional tests — does the package actually *work*?

Unlike the unit tests (which poke individual helpers), these drive whole
pipelines across module boundaries on synthetic-but-real inputs:

* a real bigWig written with pybigtools → ``loci.signal()`` extraction on both
  the Rust and pure-Python backends, checked for agreement;
* a FASTA + JASPAR PWM → ``scan_motifs`` finds a planted motif;
* BED files → ``Atlas.make`` / ``Atlas.search`` enrichment;
* a Loci → Genes → Architecture pipeline end to end;
* a ``browser()`` render (headless Agg).

Plus regression tests for bugs fixed in the v1 hardening pass.

Optional compiled deps (pybigtools, lightmotif, pysam) are ``importorskip``-ed
so a pip-only install still collects a green subset.
"""
import textwrap

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")  # headless: browser tests must not need a display

from genomeblocks import Loci, Genes, Architecture, Atlas
from genomeblocks.locus import Locus


CHROM_LEN = 20_000


# ── fixtures: synthetic-but-real files ───────────────────────────────────────

@pytest.fixture
def bigwig_path(tmp_path):
    """A real bigWig: baseline 1.0, a plateau of 8.0 over chr1:9500-10500."""
    pbt = pytest.importorskip("pybigtools")
    p = tmp_path / "signal.bw"
    h = pbt.open(str(p), "w")
    h.write(
        {"chr1": CHROM_LEN},
        iter([
            ("chr1", 0, 9500, 1.0),
            ("chr1", 9500, 10500, 8.0),
            ("chr1", 10500, CHROM_LEN, 1.0),
        ]),
    )
    h.close()
    return str(p)


@pytest.fixture
def genome_and_motif(tmp_path):
    """A 1-chrom FASTA with a TGACTCA (AP-1) core planted at the window center,
    plus a matching single-PWM JASPAR file (raw 4-line count format)."""
    core = "TGACTCA"
    seq = "A" * 3000 + core + "A" * (3000 - len(core))  # center of a 6 kb window
    fa = tmp_path / "genome.fa"
    fa.write_text(f">chr1\n{seq}\n")

    counts = {b: [] for b in "ACGT"}
    for ch in core:
        for b in counts:
            counts[b].append("100" if b == ch else "0")
    jaspar = tmp_path / "motif.jaspar"
    jaspar.write_text(
        ">M001 APONE\n" + "\n".join(" ".join(counts[b]) for b in "ACGT") + "\n"
    )
    return str(fa), str(jaspar)


@pytest.fixture
def atlas_beds(tmp_path):
    a = tmp_path / "trackA.bed"
    b = tmp_path / "trackB.bed"
    a.write_text("chr1\t1000\t2000\nchr1\t9000\t11000\n")
    b.write_text("chr1\t9000\t11000\nchr1\t15000\t16000\n")
    return str(a), str(b)


# ── signal extraction: the hot path, on both backends ────────────────────────

def test_signal_backends_agree_on_plateau(bigwig_path):
    loci = Loci([Locus("chr1", 9500, 10500)])  # centered on the plateau
    kw = dict(n_bins=10, flank=500, agg="mean", progress=False, verbose=False)

    cube_rust = loci.signal([bigwig_path], backend="pybigtools", **kw)
    cube_py = loci.signal([bigwig_path], backend="bigwig", **kw)

    assert cube_rust.shape == (1, 1, 10)
    # the whole ±500 window sits inside the 8.0 plateau
    assert np.allclose(cube_rust, 8.0, atol=1e-3)
    # the pure-Python reader must agree with the Rust reader
    assert np.allclose(cube_rust, cube_py, atol=1e-3)


def test_signal_span_mode(bigwig_path):
    # span=True uses the full locus, not center±flank
    loci = Loci([Locus("chr1", 0, CHROM_LEN)])
    cube = loci.signal([bigwig_path], n_bins=20, agg="mean", span=True,
                       progress=False, verbose=False)
    assert cube.shape == (1, 1, 20)
    # the plateau (9500-10500) straddles the 10000 bin boundary, so the two
    # central bins average to ~4.5 — clearly above the 1.0 baseline.
    assert cube.max() > 4.0
    assert cube.min() == pytest.approx(1.0, abs=1e-3)


def test_pure_python_reader_matches_pybigtools(bigwig_path):
    pbt = pytest.importorskip("pybigtools")
    from genomeblocks import bigwig as gb_bigwig

    r = gb_bigwig.open(bigwig_path)
    h = pbt.open(bigwig_path, "r")
    try:
        assert r.chroms() == {"chr1": CHROM_LEN}
        py = r.stats("chr1", 9000, 11000, n_bins=4, stat="mean")
        rust = h.values("chr1", 9000, 11000, bins=4, summary="mean",
                        exact=True, missing=0.0)
        py = [0.0 if v is None else v for v in py]
        assert np.allclose(py, rust, atol=1e-3)
    finally:
        r.close()
        h.close()


# ── motif scanning ───────────────────────────────────────────────────────────

def test_scan_motifs_finds_planted_motif(genome_and_motif):
    pytest.importorskip("lightmotif")
    from genomeblocks.motifs import make_genome, scan_motifs

    fa, jaspar = genome_and_motif
    genome = make_genome(fa)
    assert list(genome) == ["chr1"]

    loci = Loci([Locus("chr1", 0, 6000)])  # center 3000, r=3000 → full window
    hits = scan_motifs(loci, genome, jaspar, r=3000, threshold=10.0,
                       norm=False, verbose=False)
    assert hits["M001"] >= 1  # the planted TGACTCA is found


# ── Atlas enrichment ─────────────────────────────────────────────────────────

def test_atlas_make_and_search(atlas_beds):
    a, b = atlas_beds
    atlas = Atlas.make([a, b], chromsizes={"chr1": CHROM_LEN},
                       bin_size=1000, workers=1, verbose=False)
    assert len(atlas) == 2

    query = Loci([Locus("chr1", 9000, 11000)])  # overlaps both tracks
    res = atlas.search(query)
    assert list(res.columns[:3]) == ["name", "n_query_bins", "n_track_bins"]
    assert set(res["name"]) == {"trackA", "trackB"}
    assert (res["overlaps"] > 0).all()


def test_atlas_bootstrap_runs(atlas_beds):
    a, b = atlas_beds
    atlas = Atlas.make([a, b], chromsizes={"chr1": CHROM_LEN},
                       bin_size=1000, workers=1, verbose=False)
    res = atlas.bootstrap(Loci([Locus("chr1", 9000, 11000)]),
                          n=20, seed=0, verbose=False)
    assert {"observed", "expected", "z", "p_emp"} <= set(res.columns)
    assert len(res) == 2


# ── cross-module pipeline: Loci → Genes → Architecture ───────────────────────

def test_full_pipeline(bed_path, gtf_path, bedpe_path):
    # 1. peaks → set algebra
    cre = Loci.make(bed_path).sort().merge()
    assert len(cre) >= 1
    slopped = cre.slop(50)
    assert all(l2.length >= l1.length for l1, l2 in zip(cre, slopped))

    # 2. genes → region annotation
    genes = Genes.make(gtf_path, promoter_r=500)
    annot = genes.annotations(cre)
    assert "annotation" in annot.columns
    assert len(annot) == len(cre)

    # 3. architecture from loops → normalize → strength → hubs
    #    (use the 5-locus cre grid that the bedpe fixture anchors onto)
    grid = Loci([Locus("chr1", b, b + 100) for b in (1000, 5000, 9000)])
    arch = Architecture.make(grid, bedpe_path, r=2500, verbose=False)
    assert arch.n_loci >= 2 and arch.n_links >= 1

    arch.normalize(grid, source="w", name="n", verbose=False)
    arch.strength(key="n", name="strength", verbose=False)
    assert "strength" in arch.vp
    cutoff, uids = arch.elbow("strength", verbose=False)
    assert cutoff <= len(uids)


# ── browser render (headless) ────────────────────────────────────────────────

def test_browser_renders_mixed_tracks(bigwig_path, bed_path):
    # `from genomeblocks import browser` (the function) now resolves cleanly
    # since the drawing module was renamed to `browserview`.
    from genomeblocks import browser
    fig, axes = browser(("chr1", 8000, 12000),
                        {"atac": bigwig_path, "peaks": bed_path})
    assert fig is not None
    assert "atac" in axes and "peaks" in axes and "_axis" in axes
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_browser_averages_bigwig_list(bigwig_path):
    from genomeblocks.browserview import browser
    # a list of bigwigs is averaged into one track — must not raise
    fig, axes = browser(("chr1", 8000, 12000),
                        {"mean_atac": [bigwig_path, bigwig_path]})
    assert "mean_atac" in axes
    import matplotlib.pyplot as plt
    plt.close(fig)


# ── regressions for the v1 hardening pass ────────────────────────────────────

def test_signal_ram_guard_raises(monkeypatch):
    # The guard must refuse a cube larger than half of available RAM. Pin the
    # RAM estimate low so the check is deterministic regardless of the host.
    import genomeblocks.signal as sig
    monkeypatch.setattr(sig, "_available_ram_bytes", lambda: 1000)
    with pytest.raises(MemoryError):
        sig.signal(Loci([Locus("chr1", 0, 100)]), ["nonexistent.bw"],
                   n_bins=200, verbose=False, progress=False)


def test_generic_utr_without_cds_does_not_crash(tmp_path):
    # A non-coding transcript with a generic 'UTR' feature (no CDS) used to
    # crash Genes.make with IndexError on cds[0].
    gtf = tmp_path / "noncoding.gtf"
    gtf.write_text(textwrap.dedent('''\
        chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tgene_id "G1"; gene_name "NC1"; gene_type "lncRNA";
        chr1\tsrc\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
        chr1\tsrc\texon\t1000\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
        chr1\tsrc\tUTR\t1000\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
    '''))
    genes = Genes.make(str(gtf), promoter_r=500)
    assert len(genes) == 1
