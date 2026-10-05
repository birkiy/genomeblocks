"""genomeblocks: building blocks for regulatory genomics, as columnar tables.

    import genomeblocks as gb

    cre   = gb.Loci.make("atac.narrowPeak")            # intervals: numpy columns
    genes = gb.Genes.make("gencode.gtf")               # genes / transcripts / features tables
    A     = (gb.Architecture.make(cre, "loops.bedpe")  # vertices = cre rows, edges = a table
               .add_mcool("hic.mcool", resolution=5000)
               .normalize().annotate(genes))

Tables built from one another share a :class:`Genome` and join on the row
number. Heavy work runs through swappable backends (``gb.backends()`` lists
them; ``backend=`` on a call or :func:`use_backend` picks one), and every
table converts to and from pandas, polars, Arrow, bioframe, pyranges,
pybedtools, AnnData, Biopython and friends (:mod:`genomeblocks.interop`).

Names are imported on first use, so ``import genomeblocks`` stays light.
"""
from importlib import import_module
from typing import Dict

__version__ = "2.0.0"

_EXPORTS: Dict[str, tuple] = {
    # tables
    "Genome": (".genome", "Genome"),
    "Locus": (".locus", "Locus"),
    "Loci": (".loci", "Loci"),
    "Genes": (".genes", "Genes"),
    "Pairs": (".bedpe", "Pairs"),
    "Architecture": (".architecture", "Architecture"),
    "Atlas": (".atlas", "Atlas"),
    # backends
    "backends": (".backends", None),          # the module; calling it lists the backends
    "use_backend": (".backends", "use_backend"),
    # interop
    "as_loci": (".interop", "as_loci"),
    "read_fasta": (".backends.fasta", "read_fasta"),
    "load_motifs": (".backends.motifs", "load_motifs"),
    # analysis helpers
    "tmm": (".signal", "tmm"),
    "compare_heatmap": (".signal_draw", "compare_heatmap"),
    "plot_motif_heatmap": (".motifs_draw", "plot_motif_heatmap"),
    "coverage": (".bam", "coverage"),
    # viewers
    "browser": (".browserview", "browser"),
    "View": (".view", "View"),
    "igv_html": (".igv", "igv_html"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    if name in _EXPORTS:
        module, attr = _EXPORTS[name]
        mod = import_module(module, __name__)
        value = mod if attr is None else getattr(mod, attr)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + __all__)
