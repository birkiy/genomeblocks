"""End to end: peaks → genes → architecture → signal → motifs → enrichment → views."""
import numpy as np

import genomeblocks as gb
from genomeblocks import Architecture, Genes, Loci


def test_pipeline(bed_path, gtf_path, bedpe_path, bw_path, bw2_path, fasta_path, jaspar_path, tmp_path):
    cre = Loci.make(bed_path).slop(50).sort().merge()
    assert len(cre) == 7
    genes = Genes.make(gtf_path).select_isoforms(cre, bw_path, r=500, verbose=False)
    labels = genes.annotations(cre)
    assert "Promoter-TSS" in set(labels["annotation"])

    A = (Architecture.make(cre, bedpe_path, r=200, verbose=False))
    A.ep["w"][:] = np.arange(1, A.n_links + 1, dtype=float)
    A.normalize(verbose=False).annotate(genes, verbose=False).strength(verbose=False)
    hubs = A.prime_hubs(verbose=False)
    assert "prime_genes" in hubs and A.components().max() >= 0
    assert A.support(genes)["GENE_A"]

    S = cre.signal([bw_path, bw2_path], n_bins=20, flank=500, verbose=False, progress=False)
    assert S.shape == (7, 2, 20) and S.max() > 0
    norm = gb.tmm(S)
    assert norm.shape == S.shape
    fig = cre.plot_heatmap(S, groups={"promoter": labels["annotation"].to_numpy() == "Promoter-TSS"})
    assert fig.axes

    M = cre.scan_motifs_matrix(fasta_path, jaspar_path, r=100, threshold=7.0, verbose=False)
    assert M.shape == (7, 2)
    prof, names = cre.scan_motifs_profile(fasta_path, jaspar_path, r=100, n_bins=10, threshold=7.0, verbose=False)
    assert prof.shape == (7, 2, 10) and names == ["M1", "M2"]

    atlas = gb.Atlas.make(bed_path, chromsizes={"chr1": 20_000, "chr2": 8_000}, bin_size=500, verbose=False)
    assert len(cre.enrich(atlas)) == 1

    se_regions = cre.call_se(bw_path, stitch=2000)
    assert "rank" in se_regions.cols

    fig, axes = gb.browser("chr1:0-12 kb", {"ATAC": bw_path, "CRE": cre, "loops": bedpe_path, "genes": genes})
    assert "genes" in axes
    gb.igv_html(str(tmp_path / "igv.html"), regions=["chr1:0-12 kb"], loci={"CRE": cre}, genes=genes,
                signal={"ATAC": bw_path}, architecture=A)
    v = gb.View(A, genes=genes).signal("ATAC", bw_path).cres().loops().genes()
    assert v.save(str(tmp_path / "view.html"))["total"] > 0

    # everything hands itself to pandas / polars and back
    assert Loci.from_frame(cre.to_polars()).equals(cre)
    assert Architecture.from_frame(cre, A.to_polars()).n_links == A.n_links
    A.save(str(tmp_path / "A"))
    assert Architecture.load(str(tmp_path / "A")).n_links == A.n_links
