"""Public-surface guards: what should be importable, what must stay gone."""
import subprocess
import sys

import genomeblocks as gb


EXPECTED_ALL = [
    "Architecture", "Atlas", "CDS", "Exon", "Gene", "Genes", "Loci", "Locus",
    "Transcript", "UTR", "browser", "compare_heatmap", "coverage",
    "make_genome", "scan_motifs", "tmm",
]


def test_public_all():
    assert gb.__all__ == EXPECTED_ALL


def test_kept_names_importable():
    for name in EXPECTED_ALL:
        assert getattr(gb, name) is not None


def test_removed_top_level_names_are_gone():
    for name in ["Tags", "draw"]:
        assert not hasattr(gb, name)


def test_removed_architecture_methods_gone():
    from genomeblocks import Architecture
    for m in ["make_spread", "make_clique", "surrounds", "add_patches",
              "weighted_activity", "cluster", "focus", "focus_genes",
              "supernode_metrics", "supernode_abc", "mutual", "aggregate", "draw"]:
        assert not hasattr(Architecture, m), m


def test_removed_genes_methods_gone():
    from genomeblocks import Genes
    for m in ["get_tss_transcripts", "nearest_transcripts", "enhancer_to_genes",
              "cre_supernodes", "nearby"]:
        assert not hasattr(Genes, m), m


def test_tags_module_removed():
    import importlib.util
    assert importlib.util.find_spec("genomeblocks.tags") is None


def test_loci_lost_tag_method():
    from genomeblocks import Loci
    assert not hasattr(Loci, "tag")


def test_signal_draw_separated():
    from genomeblocks import signal as sig
    from genomeblocks import signal_draw as sd
    assert hasattr(sd, "plot_heatmap") and hasattr(sd, "compare_heatmap")
    assert not hasattr(sig, "plot_heatmap")       # processing module is plot-free


def test_signal_first_import_order():
    # importing genomeblocks.signal first must not deadlock on signal_draw
    r = subprocess.run([sys.executable, "-c", "import genomeblocks.signal"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_browser_export_is_the_function_not_the_module():
    # Regression: the drawing module was renamed to `browserview` so the public
    # `browser` name is unambiguously the callable, regardless of import order.
    import types
    assert callable(gb.browser)
    assert not isinstance(gb.browser, types.ModuleType)


def test_browser_submodule_name_is_gone():
    import importlib.util
    assert importlib.util.find_spec("genomeblocks.browser") is None
    assert importlib.util.find_spec("genomeblocks.browserview") is not None
