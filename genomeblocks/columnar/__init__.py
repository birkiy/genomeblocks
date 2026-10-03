"""genomeblocks.columnar — prototype of a table-first genomeblocks.

Everything is a table of numpy columns that shares one :class:`Genome`, and
the *row number* is the join key between tables:

    Genome ── Loci (CREs) ──┬── annotations / signal / motif hits   (row i = CRE i)
                            └── Architecture (vertices = Loci rows, edges table)
           └─ Genes (genes ⇄ transcripts ⇄ features, linked by row numbers)

Objects (Locus, Gene) are still there when you look at a single row, and the
classic ``genomeblocks`` functions keep working on these tables.

    >>> import genomeblocks.columnar as gbc
    >>> cre = gbc.Loci.make("peaks.bed")
    >>> genes = gbc.Genes.make("genes.gtf")
    >>> A = gbc.Architecture.make(cre, "loops.bedpe")
    >>> A.add_mcool("hic.mcool", resolution=5000).normalize().annotate(genes)

Heavy dependencies (pandas, graph-tool, cooler, scipy) load on first use.
"""
from importlib import import_module

_EXPORTS = {
    "Genome": (".genome", "Genome"),
    "default_genome": (".genome", "default_genome"),
    "set_default_genome": (".genome", "set_default_genome"),
    "Loci": (".loci", "Loci"),
    "LocusView": (".loci", "LocusView"),
    "Genes": (".genes", "Genes"),
    "Architecture": (".architecture", "Architecture"),
}
__all__ = sorted(_EXPORTS)


def __getattr__(name):
    if name in _EXPORTS:
        mod, attr = _EXPORTS[name]
        value = getattr(import_module(mod, __name__), attr)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + __all__)
