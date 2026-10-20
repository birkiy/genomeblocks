"""Expression tables: every quantifier format reads to the same numbers, the
id rules line them up with a Genes / Loci table, transcripts summarise to
genes exactly as tximport, normalisation matches edgeR / DESeq2, and the table
behaves like the others (protocols, save / load, AnnData)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import genomeblocks as gb
from conftest import (RNASEQ_GENES, RNASEQ_SAMPLES, RNASEQ_TX, installed_backends, requires, rnaseq_gene_counts,
                      rnaseq_matrices)
from genomeblocks import Expression
from genomeblocks.expression import sample_name, sniff_format

# ── reference numbers computed in R on the fixture files ───────────────────────
# tximport 1.34.0: tximport(files, type="salmon", tx2gene, countsFromAbundance=..., ignoreTxVersion=TRUE,
# ignoreAfterBar=TRUE); columns WT_1, WT_2, WT_3, KO_1, KO_2, KO_3
TXI_LENGTH = [[1872.40793201, 1807.39795918, 1900.23360965, 2234.000832640, 2293.407968390, 2229.88888889],
              [620.0, 625.0, 630.0, 629.931016637, 629.931016637, 645.0],       # KO_1 / KO_2: geometric mean of the rest
              [9632.5] * 6]                                                     # nothing anywhere: plain mean
TXI_ABUNDANCE = [[136548.565197, 156753.156283, 149512.845307, 845421.215147, 871233.755131, 309956.721807],
                 [824756.827678, 793282.937333, 822351.946362, 0.0, 0.0, 623464.682166],
                 [0.0] * 6]
TXI_SCALED_TPM = [[63.920222225, 72.5987035945, 73.8437757624, 300.0, 310.0, 144.44834178],
                  [386.079777775, 367.4012964055, 406.1562242376, 0.0, 0.0, 290.55165822],
                  [0.0] * 6]
TXI_LENGTH_SCALED_TPM = [[157.866049144, 172.516506308, 178.762211502, 300.0, 310.0, 269.139551437],
                         [292.133950856, 267.483493692, 301.237788498, 0.0, 0.0, 165.860448563],
                         [0.0] * 6]

# a 20 x 6 count matrix (negative binomial, seeded; one all-zero row, zeros in some samples)
MATRIX20 = [
    [124, 219, 67, 132, 163, 176], [26, 41, 8, 18, 22, 31], [93, 360, 220, 235, 145, 583], [0, 0, 0, 0, 0, 0],
    [2, 13, 2, 3, 5, 15], [7, 14, 7, 14, 22, 18], [102, 52, 49, 141, 71, 184], [0, 0, 0, 52, 126, 343],
    [153, 104, 68, 113, 72, 114], [58, 69, 17, 59, 28, 108], [283, 682, 271, 497, 220, 752],
    [344, 190, 263, 286, 139, 565], [154, 0, 95, 84, 70, 169], [403, 224, 277, 606, 430, 1055],
    [150, 85, 100, 179, 131, 203], [22, 61, 15, 67, 23, 52], [154, 82, 75, 196, 162, 273], [25, 30, 28, 20, 10, 57],
    [254, 261, 159, 579, 279, 809], [55, 105, 73, 238, 42, 118],
]
# edgeR 4.4.0: calcNormFactors(DGEList(x))$samples$norm.factors
EDGER_TMM = [0.992391291459, 1.234910059062, 0.963966821255, 0.985796821890, 1.049536654095, 0.818153168846]
# edgeR: cpm(calcNormFactors(y), log=TRUE, prior.count=2) rows g1, g4, g8
EDGER_LOGCPM_TMM = {0: [15.67938415772, 16.07439852647, 15.26476396245, 15.23922837512, 16.14534005573, 15.24637043622],
                    3: [9.40645206712] * 6,
                    7: [9.40645206712, 9.40645206712, 9.40645206712, 13.93369987355, 15.77785279884, 16.19668789611]}
EDGER_CPM_G1 = [51473.640514736, 84490.74074074, 37346.71125975, 37510.656436488, 75462.96296296, 31288.88888889]
# edgeR: filterByExpr(y, group=WT/KO) and filterByExpr(x) (no group)
EDGER_FILTER_GROUP = [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
EDGER_FILTER_NOGROUP = [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1]
# DESeq2 1.46.0: estimateSizeFactorsForMatrix(x)
DESEQ2_SIZE_FACTORS = [0.924202879550, 1.179665237662, 0.560282189486, 1.230718449121, 0.802741783002, 1.922092948553]


def in_order(E):
    """Samples in WT_1..KO_3 order whatever order the files were read in."""
    return E.take_samples(RNASEQ_SAMPLES)


@pytest.fixture
def matrix20():
    df = pd.DataFrame(np.array(MATRIX20), columns=RNASEQ_SAMPLES)
    df.insert(0, "gene", [f"g{i + 1}" for i in range(20)])
    return Expression.from_frame(df, verbose=False)


# ── formats ──────────────────────────────────────────────────────────────────

def test_sniff_and_sample_names(rnaseq_dir):
    d = rnaseq_dir
    assert sniff_format(d / "salmon/WT_1/quant.sf")["format"] == "salmon"
    assert sniff_format(d / "kallisto/WT_1/abundance.tsv") == {"format": "kallisto", "header": True, "skip": 0,
                                                              "sep": "\t", "level": "transcript"}
    assert sniff_format(d / "rsem/WT_1.isoforms.results")["level"] == "transcript"
    assert sniff_format(d / "rsem/WT_1.genes.results")["level"] == "gene"
    assert sniff_format(d / "counts.featureCounts.txt") == {"format": "featurecounts", "header": True, "skip": 1,
                                                           "sep": "\t", "level": None}
    assert sniff_format(d / "star/WT_1.ReadsPerGene.out.tab")["format"] == "star"
    assert sniff_format(d / "htseq/WT_1.counts.txt")["format"] == "htseq"
    assert sniff_format(d / "nfcore/salmon.merged.gene_counts.tsv") == {"format": "matrix", "header": True, "skip": 0,
                                                                      "sep": "\t", "level": "gene"}
    assert sniff_format(d / "by_name.csv")["sep"] == ","
    assert sample_name(d / "salmon/WT_1/quant.sf") == "WT_1"
    assert sample_name("x/KO_2/abundance.tsv") == "KO_2"
    assert sample_name("star/KO_3.ReadsPerGene.out.tab") == "KO_3"
    assert sample_name("rsem/WT_2.isoforms.results") == "WT_2"
    assert sample_name("htseq/WT_1.counts.txt.gz") == "WT_1"


def test_every_format_holds_the_same_numbers(rnaseq_dir, genes):
    d = rnaseq_dir
    gc = rnaseq_gene_counts()
    tx_counts, eff, tpm = rnaseq_matrices()
    sources = {
        "salmon": (str(d / "salmon/*/quant.sf"), ["counts", "tpm", "length"]),
        "kallisto": (str(d / "kallisto/*/abundance.tsv"), ["counts", "tpm", "length"]),
        "rsem isoforms": (str(d / "rsem/*.isoforms.results"), ["counts", "tpm", "length"]),
        "rsem genes": (str(d / "rsem/*.genes.results"), ["counts", "tpm", "length"]),
        "featureCounts": (str(d / "counts.featureCounts.txt"), ["counts", "length"]),
        "STAR": (str(d / "star/*.ReadsPerGene.out.tab"), ["counts"]),
        "HTSeq": (str(d / "htseq/*.counts.txt"), ["counts"]),
        "nf-core": ([str(d / "nfcore/salmon.merged.gene_counts.tsv"), str(d / "nfcore/salmon.merged.gene_tpm.tsv")],
                    ["counts", "tpm"]),
    }
    for name, (src, layers) in sources.items():
        E = in_order(Expression.make(src, genes=genes, how="counts", verbose=False))
        assert E.level == "gene" and len(E) == 3 and E.n_samples == 6, name
        assert list(E.layers) == layers, name
        assert E.ids.tolist() == list(RNASEQ_GENES) and E.names.tolist() == list(RNASEQ_GENES.values()), name
        np.testing.assert_allclose(E.counts, gc, rtol=1e-6, err_msg=name)
        assert not E.missing.any(), name
        assert E.sample_names == RNASEQ_SAMPLES
    # transcript level: the four GTF transcripts, the spike-in dropped
    E = in_order(Expression.make(str(d / "salmon/*/quant.sf"), genes=genes, level="transcript", verbose=False))
    assert E.level == "transcript" and E.ids.tolist() == ["T1", "T1b", "T2", "T3"]
    assert E.names.tolist() == ["GENE_A", "GENE_A", "GENE_B", "GENE_C"]
    np.testing.assert_allclose(E.counts, tx_counts[:4])
    np.testing.assert_allclose(E.tpm, tpm[:4], rtol=1e-6)
    np.testing.assert_allclose(E.length, eff[:4])
    assert E.unmatched.tolist() == ["ERCC-00002"] and E.meta["key_rule"].startswith("transcript_id, pipe header cut")
    # a dict of files names the samples; a list keeps its order
    E = Expression.make({f"s{j}": str(d / f"salmon/{s}/quant.sf") for j, s in enumerate(RNASEQ_SAMPLES)}, verbose=False)
    assert E.sample_names == [f"s{j}" for j in range(6)] and E.rows is None and E.level == "transcript"
    assert len(E) == 5 and "ERCC-00002" in E
    E = Expression.make([str(d / f"salmon/{s}/quant.sf") for s in RNASEQ_SAMPLES[::-1]], verbose=False)
    assert E.sample_names == RNASEQ_SAMPLES[::-1]


@pytest.mark.parametrize("backend", installed_backends("tables"))
def test_parsers_agree(rnaseq_dir, genes, backend):
    ref = Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, verbose=False)
    E = Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, verbose=False, backend=backend)
    for k in ref.layers:
        np.testing.assert_array_equal(E.layers[k], ref.layers[k])
    for src in ("counts.featureCounts.txt", "star/*.ReadsPerGene.out.tab", "htseq/*.counts.txt", "by_name.csv"):
        a = Expression.make(str(rnaseq_dir / src), verbose=False)
        b = Expression.make(str(rnaseq_dir / src), verbose=False, backend=backend)
        np.testing.assert_array_equal(a.counts, b.counts)
        assert a.ids.tolist() == b.ids.tolist() and a.sample_names == b.sample_names


def test_featurecounts_star_htseq_details(rnaseq_dir, genes, capsys):
    d = rnaseq_dir
    E = Expression.make(str(d / "counts.featureCounts.txt"), verbose=False)
    assert E.sample_names == RNASEQ_SAMPLES                      # bam/WT_1.sorted.bam -> WT_1
    assert E.layers["length"][:, 0].tolist() == [1800, 800, 9800] and E.rows is None and E.level is None
    E = Expression.make(str(d / "counts.featureCounts.txt"), sample_names=list("abcdef"), verbose=False)
    assert E.sample_names == list("abcdef")
    with pytest.raises(ValueError, match="sample_names"):
        Expression.make(str(d / "counts.featureCounts.txt"), sample_names=["a"], verbose=False)
    # STAR: the strand column
    gc = rnaseq_gene_counts()
    star = str(d / "star/*.ReadsPerGene.out.tab")
    rev = in_order(Expression.make(star, genes=genes, verbose=False))
    np.testing.assert_array_equal(rev.counts, gc)
    np.testing.assert_array_equal(in_order(Expression.make(star, genes=genes, strand="reverse", verbose=False)).counts, gc)
    np.testing.assert_array_equal(in_order(Expression.make(star, genes=genes, strand="none", verbose=False)).counts, gc + 7)
    np.testing.assert_array_equal(in_order(Expression.make(star, genes=genes, strand="forward", verbose=False)).counts,
                                  gc // 20)
    with pytest.raises(ValueError, match="strand"):
        Expression.make(star, strand="sideways", verbose=False)
    Expression.make(star, verbose=True)
    out = capsys.readouterr().out
    assert "reverse-stranded" in out and "6 files" in out
    # HTSeq: the summary rows are gone, once per file
    Expression.make(str(d / "htseq/*.counts.txt"), verbose=True)
    out = capsys.readouterr().out
    assert "HTSeq summary rows" in out and "6 files" in out


def test_union_of_files_with_different_features(rnaseq_dir, tmp_path):
    a = rnaseq_dir / "salmon/WT_1/quant.sf"
    lines = a.read_text().splitlines()
    (tmp_path / "b").mkdir()
    (tmp_path / "b" / "quant.sf").write_text("\n".join(lines[:3] + ["NEW.1\t500\t400\t10\t7"]) + "\n")
    E = Expression.make([str(a), str(tmp_path / "b" / "quant.sf")], verbose=False)
    assert len(E) == 6 and E.sample_names == ["WT_1", "b"]
    assert np.isnan(E.counts[5, 0]) and E.counts[5, 1] == 7 and np.isnan(E.counts[3, 1])


def test_nfcore_matrices_and_errors(rnaseq_dir, tmp_path):
    d = rnaseq_dir
    E = Expression.make(str(d / "nfcore/salmon.merged.gene_tpm.tsv"), verbose=False)
    assert list(E.layers) == ["tpm"] and E.names.tolist() == list(RNASEQ_GENES.values())
    assert E.default_layer == "tpm" and E.X is E.tpm and E.shape == (3, 6) and not E.missing.any()
    assert E.to_pandas().columns[0] == "id" and E["GENE_A"].tolist() == E.tpm[0].tolist()
    with pytest.raises(AttributeError, match="no counts"):
        E.counts
    with pytest.raises(AttributeError, match="no counts"):
        E.tmm()
    both = Expression.make([str(d / "nfcore/salmon.merged.gene_counts.tsv"), str(d / "nfcore/salmon.merged.gene_tpm.tsv")],
                           verbose=False)
    assert list(both.layers) == ["counts", "tpm"]
    with pytest.raises(ValueError, match="mix"):
        Expression.make([str(d / "nfcore/salmon.merged.gene_counts.tsv"), str(d / "salmon/WT_1/quant.sf")], verbose=False)
    (tmp_path / "other.tsv").write_text("gene_id\tgene_name\ta\tb\nG1\tx\t1\t2\nG2\ty\t3\t4\nG3\tz\t5\t6\n")
    with pytest.raises(ValueError, match="samples"):
        Expression.make([str(d / "nfcore/salmon.merged.gene_counts.tsv"), str(tmp_path / "other.tsv")], verbose=False)
    (tmp_path / "text.tsv").write_text("gene_id\tnote\nG1\tx\n")
    with pytest.raises(ValueError, match="numeric"):
        Expression.make(str(tmp_path / "text.tsv"), verbose=False)
    with pytest.raises(FileNotFoundError):
        Expression.make(str(tmp_path / "nothing/*/quant.sf"))


# ── alignment ────────────────────────────────────────────────────────────────

def test_alignment_rules(rnaseq_dir, genes):
    d = rnaseq_dir
    # a gene-name keyed CSV: GENE_C missing, the text column ignored
    E = in_order(Expression.make(str(d / "by_name.csv"), genes=genes, verbose=False))
    assert E.meta["key_rule"] == "gene_name" and E.missing.tolist() == [False, False, True]
    assert np.isnan(E.counts[2]).all() and E.counts[0].tolist() == rnaseq_gene_counts()[0].tolist()
    assert E.lib_size.tolist() == (rnaseq_gene_counts()[:2].sum(axis=0)).tolist()
    # explicit keys
    assert Expression.make(str(d / "by_name.csv"), genes=genes, key="gene_name", verbose=False).missing.sum() == 1
    with pytest.raises(ValueError, match="none of the"):
        Expression.make(str(d / "by_name.csv"), genes=genes, key="gene_id", verbose=False)
    with pytest.raises(ValueError, match="key"):
        Expression.make(str(d / "by_name.csv"), genes=genes, key="uid", verbose=False)
    # strict
    with pytest.raises(ValueError, match="ERCC-00002"):
        Expression.make(str(d / "salmon/*/quant.sf"), genes=genes, strict=True, verbose=False)
    # align later, and re-align
    U = Expression.make(str(d / "salmon/*/quant.sf"), verbose=False)
    assert U.rows is None and U.table is None
    A = U.align(genes=genes, verbose=False)
    assert A.level == "gene" and A.table is genes.genes and A.rows is genes
    A2 = A.align(genes=genes, verbose=False)
    np.testing.assert_allclose(A2.counts, A.counts)
    T = U.align(genes=genes, level="transcript", verbose=False)
    assert T.table is genes.transcripts
    with pytest.raises(ValueError, match="gene-level"):
        A.align(genes=genes, level="transcript")
    with pytest.raises(ValueError, match="exactly one"):
        U.align()
    with pytest.raises(ValueError, match="needs genes="):
        Expression.make(str(d / "salmon/*/quant.sf"), level="gene", verbose=False)
    with pytest.raises(TypeError):
        U.align(genes=genes.genes)
    # duplicated input ids on one row are summed; versions stripped on both sides
    df = pd.DataFrame({"id": ["G1.5", "G1.5", "G2", "G9"], "a": [1.0, 2.0, 3.0, 4.0], "b": [5.0, 6.0, 7.0, 8.0]})
    E = Expression.from_frame(df, genes=genes, verbose=False)
    assert E.counts[0].tolist() == [3.0, 11.0] and E.counts[1].tolist() == [3.0, 7.0] and np.isnan(E.counts[2]).all()
    assert E.unmatched.tolist() == ["G9"] and E.meta["key_rule"] == "gene_id, version stripped"
    assert repr(E).startswith("Expression(3 genes × 2 samples")


def test_loci_alignment(cre):
    df = pd.DataFrame({"peak": ["chr1:900-1100", "chr1_1900_2100", "chr9:1-2", "chr2:500-600(.)"],
                       "x.bam": [1, 2, 3, 4], "y.bam": [5, 6, 7, 8]})
    P = Expression.from_frame(df, loci=cre, verbose=False)
    assert P.level == "loci" and P.rows is cre and len(P) == len(cre) and P.sample_names == ["x", "y"]
    assert P.ids[0] == cre.uid[0] and P.names is None
    assert P.counts[0].tolist() == [1, 5] and P.counts[1].tolist() == [2, 6] and P.counts[5].tolist() == [4, 8]
    assert P.missing.sum() == 4 and P.unmatched.tolist() == ["chr9:1-2"]
    assert P.meta["key_rule"] == "region (chrom:start-end)"
    named = cre.copy()
    named["name"] = [f"p{i}" for i in range(len(cre))]
    df = pd.DataFrame({"GeneID": ["p1", "p0"], "a": [1, 2]})
    Q = Expression.from_frame(df, loci=named, verbose=False)
    assert Q.counts[:2, 0].tolist() == [2, 1] and Q.meta["key_rule"] == "name" and Q.ids[0] == "p0"
    assert Q.get("p1").tolist() == [1] and Q.row(cre.uid[0]) == 0
    with pytest.raises(ValueError, match="level"):
        Expression.from_frame(df, loci=named, level="gene", verbose=False)
    # anything as_loci takes
    R = Expression.from_frame(df, loci=named.to_pandas(), key="name", verbose=False)
    assert R.table.equals(named)


def test_genes_rows_missing(genes):
    assert genes.rows(["GENE_A", "G3", "nope"], missing=-1).tolist() == [0, 2, -1]
    with pytest.raises(KeyError):
        genes.rows(["nope"])


# ── tximport ─────────────────────────────────────────────────────────────────

def test_to_genes_matches_tximport(rnaseq_dir, genes):
    tx = in_order(Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, level="transcript", verbose=False))
    E = tx.to_genes(how="counts", verbose=False)
    assert E.level == "gene" and E.table is genes.genes and E.meta["counts_from"] == "counts"
    np.testing.assert_allclose(E.counts, rnaseq_gene_counts())
    np.testing.assert_allclose(E.length, TXI_LENGTH, rtol=1e-8)
    np.testing.assert_allclose(E.tpm, TXI_ABUNDANCE, rtol=1e-8)
    np.testing.assert_allclose(tx.to_genes(how="scaled_tpm", verbose=False).counts, TXI_SCALED_TPM, rtol=1e-8)
    E = tx.to_genes(how="length_scaled_tpm", verbose=False)
    np.testing.assert_allclose(E.counts, TXI_LENGTH_SCALED_TPM, rtol=1e-8)
    np.testing.assert_allclose(E.length, TXI_LENGTH, rtol=1e-8)
    # the default in make
    D = in_order(Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, verbose=False))
    np.testing.assert_allclose(D.counts, TXI_LENGTH_SCALED_TPM, rtol=1e-8)
    assert D.meta["counts_from"] == "length_scaled_tpm" and D.unmatched.tolist() == ["ERCC-00002"]
    with pytest.raises(ValueError, match="how"):
        tx.to_genes(how="magic")
    with pytest.raises(ValueError, match="transcript-level"):
        E.to_genes()
    # a transcript with no TPM layer cannot be TPM-scaled
    fc = Expression.make(str(rnaseq_dir / "counts.featureCounts.txt"), verbose=False)
    tx_like = Expression(dict(counts=np.ones((4, 2))), ["a", "b"], ids=["T1", "T1b", "T2", "T3"]).align(
        genes=genes, level="transcript", verbose=False)
    np.testing.assert_array_equal(tx_like.to_genes(how="counts", verbose=False).counts, [[2, 2], [1, 1], [1, 1]])
    with pytest.raises(ValueError, match="tpm"):
        tx_like.to_genes()
    assert fc.level is None


# ── normalisation ────────────────────────────────────────────────────────────

def test_normalisation_matches_edger_and_deseq2(matrix20):
    E = matrix20
    np.testing.assert_allclose(E.tmm(), EDGER_TMM, rtol=1e-9)
    np.testing.assert_allclose(E.size_factors(), DESEQ2_SIZE_FACTORS, rtol=1e-9)
    np.testing.assert_allclose(E.cpm()[0], EDGER_CPM_G1, rtol=1e-9)
    lc = E.log_cpm(norm="tmm")
    for i, row in EDGER_LOGCPM_TMM.items():
        np.testing.assert_allclose(lc[i], row, rtol=1e-9)
    assert E.filter_by_expr(group=["WT"] * 3 + ["KO"] * 3).astype(int).tolist() == EDGER_FILTER_GROUP
    E.samples["condition"] = ["WT"] * 3 + ["KO"] * 3
    assert E.filter_by_expr("condition").astype(int).tolist() == EDGER_FILTER_GROUP
    assert E.filter_by_expr().astype(int).tolist() == EDGER_FILTER_NOGROUP
    np.testing.assert_allclose(E.cpm(norm="size")[0], E.counts[0] / (E.lib_size * E.size_factors()) * 1e6)
    np.testing.assert_array_equal(E.log_cpm(factors=E.tmm()), lc)
    with pytest.raises(ValueError, match="norm"):
        E.cpm(norm="quantile")
    with pytest.raises(ValueError, match="factors"):
        E.cpm(norm="tmm", factors=E.tmm())
    # NaN rows take no part and stay NaN
    X = E.counts.copy()
    X[0] = np.nan
    N = Expression(dict(counts=X), E.samples, ids=E.ids)
    assert np.isnan(N.cpm()[0]).all() and not N.filter_by_expr()[0]
    assert np.isfinite(N.tmm()).all() and np.isfinite(N.size_factors()).all()
    assert N.lib_size.tolist() == (E.lib_size - E.counts[0]).tolist()
    # tpm from counts and lengths when the quantifier gave none
    L = Expression(dict(counts=np.array([[10.0, 20.0], [30.0, 40.0]]), length=np.array([[1000.0, 1000.0], [500.0, 500.0]])),
                   ["a", "b"], ids=["x", "y"])
    np.testing.assert_allclose(L.tpm.sum(axis=0), [1e6, 1e6])
    np.testing.assert_allclose(L.tpm[1] / L.tpm[0], [6.0, 4.0])
    with pytest.raises(AttributeError):
        Expression(dict(counts=np.ones((2, 2))), ["a", "b"], ids=["x", "y"]).tpm


def test_signal_tmm_unchanged_rule(bw_path, cre):
    from genomeblocks._norm import tmm_factors
    S = cre.signal(bw_path, n_bins=4, flank=200)
    cube = np.concatenate([S, S * 2.0], axis=1)
    out = gb.tmm(cube)
    np.testing.assert_allclose(out[:, 0], out[:, 1])
    assert tmm_factors(np.array([[10.0, 20.0], [5.0, 10.0], [1.0, 2.0]])).tolist() == pytest.approx([1.0, 1.0])


# ── samples, lookups, subsets ────────────────────────────────────────────────

def test_samples_metadata(rnaseq_dir, genes):
    src = str(rnaseq_dir / "counts.featureCounts.txt")
    E = Expression.make(src, genes=genes, samples=str(rnaseq_dir / "samples.tsv"), verbose=False)
    assert E.samples.index.tolist() == RNASEQ_SAMPLES and E.samples["condition"].tolist() == ["WT"] * 3 + ["KO"] * 3
    assert E.samples.index.name == "sample"
    meta = pd.DataFrame({"condition": ["KO", "WT"]}, index=["KO_1", "WT_1"])
    with pytest.raises(KeyError, match="lacks"):
        Expression.make(src, samples=meta, verbose=False)
    E = Expression.make(src, verbose=False)
    df = pd.DataFrame({"g": range(6)}, index=RNASEQ_SAMPLES[::-1])
    E.samples = df                                              # a permutation is put back in column order
    assert E.samples["g"].tolist() == [5, 4, 3, 2, 1, 0]
    with pytest.raises(ValueError):
        E.samples = df.iloc[:3]
    with pytest.raises(TypeError):
        E.samples = {"a": 1}
    with pytest.raises(ValueError, match="duplicated"):
        Expression(dict(counts=np.ones((1, 2))), ["a", "a"], ids=["x"])
    with pytest.raises(ValueError, match="layer"):
        Expression({}, ["a", "b"], ids=["x"])
    with pytest.raises(ValueError, match="shape"):
        Expression(dict(counts=np.ones((1, 2)), tpm=np.ones((2, 2))), ["a", "b"], ids=["x"])


def test_lookups_and_subsets(rnaseq_dir, genes):
    E = in_order(Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, how="counts",
                                 samples=str(rnaseq_dir / "samples.tsv"), verbose=False))
    gc = rnaseq_gene_counts()
    assert E["GENE_A"].tolist() == gc[0].tolist() and E["G3"].name == "G3" and E[1].tolist() == gc[1].tolist()
    assert E.get("GENE_B", "tpm").index.tolist() == RNASEQ_SAMPLES
    assert E[["GENE_A", "G3"]].shape == (2, 6) and E[["GENE_A", "G3"]].index.tolist() == ["GENE_A", "G3"]
    assert "GENE_A" in E and "G2" in E and "nope" not in E and E.row("GENE_C") == 2
    with pytest.raises(KeyError):
        E["nope"]
    with pytest.raises(IndexError):
        E[3]
    got = E.lookup(["GENE_C", "", "nope", "GENE_A"])
    np.testing.assert_allclose(got, [gc[2].mean(), np.nan, np.nan, gc[0].mean()])
    np.testing.assert_allclose(E.lookup(["G1"], agg="max"), [gc[0].max()])
    np.testing.assert_allclose(E.lookup(np.array(["GENE_B"]), pd.Series(E.tpm.mean(axis=1))), [E.tpm[1].mean()])
    assert E.lookup(["GENE_B", "x"], np.array(["b1", "b2", "b3"], dtype=object), missing="").tolist() == ["b2", ""]
    assert E.lookup("GENE_A", missing=-1.0).tolist() == [gc[0].mean()]
    with pytest.raises(ValueError, match="values"):
        E.lookup(["GENE_A"], [1.0])
    # sample subsets keep the alignment, row subsets drop it
    K = E[:, E.samples.condition == "KO"]
    assert K.sample_names == RNASEQ_SAMPLES[3:] and K.rows is genes and K.samples["batch"].tolist() == ["b1", "b2", "b3"]
    assert E.take_samples(["KO_1", "WT_1"]).sample_names == ["KO_1", "WT_1"]
    assert E.take_samples("WT_2").counts.shape == (3, 1) and E[:, 1:3].n_samples == 2
    S = E[[True, False, True]]
    assert S.rows is None and S.ids.tolist() == ["G1", "G3"] and S.names.tolist() == ["GENE_A", "GENE_C"]
    assert E[1:, :2].shape == (2, 2) and E.take([2]).ids.tolist() == ["G3"]
    with pytest.raises(IndexError):
        E[0, 1, 2]
    assert list(E) == ["G1", "G2", "G3"]


# ── table behaviour, exits, persistence ──────────────────────────────────────

def test_as_a_table(rnaseq_dir, genes):
    E = in_order(Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, how="counts",
                                 samples=str(rnaseq_dir / "samples.tsv"), verbose=False))
    assert E.shape == (3, 6) and E.columns == RNASEQ_SAMPLES and len(E) == 3
    wide = E.to_pandas()
    assert wide.columns.tolist() == ["gene_id", "gene_name"] + RNASEQ_SAMPLES and len(wide) == 3
    assert E.to_pandas("tpm").iloc[0, 2] == pytest.approx(E.tpm[0, 0])
    long = E.to_pandas(long=True)
    assert long.columns.tolist() == ["gene_id", "gene_name", "sample", "counts", "condition", "batch"]
    assert len(long) == 18 and long["counts"].iloc[1] == E.counts[0, 1] and long["condition"].iloc[4] == "KO"
    assert len(E.head(2)) == 2 and len(E.tail(1)) == 1
    d = E.describe()
    assert d.loc["rows", "value"] == 3 and d.loc["level", "value"] == "gene" and d.loc["input ids dropped", "value"] == 1
    assert E.to_arrow().schema.names == wide.columns.tolist()
    import polars as pl
    assert pl.DataFrame(E).shape == (3, 8) and E.to_polars("tpm").shape == (3, 8)
    assert "Expression ·" in E._repr_html_() and "genes.genes" in repr(E)
    U = Expression.make(str(rnaseq_dir / "salmon/WT_1/quant.sf"), verbose=False)
    assert U.to_pandas().columns[0] == "id" and U.names is None and U.describe().loc["aligned to", "value"] == "nothing"


def test_from_frame_inputs(genes):
    import polars as pl
    df = pl.DataFrame({"gene_id": ["G1", "G2", "G3"], "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
                       "s1": [1.0, 2.0, 3.0], "s2": [4.0, 5.0, 6.0]})
    E = Expression.from_frame(df, genes=genes, verbose=False)
    assert E.counts.tolist() == [[1, 4], [2, 5], [3, 6]] and E.names.tolist() == ["GENE_A", "GENE_B", "GENE_C"]
    E = Expression.from_frame({"symbol": ["GENE_A", "GENE_B"], "a": [1, 2]}, ids="symbol", layer="tpm", genes=genes,
                              verbose=False)
    assert list(E.layers) == ["tpm"] and E.meta["key_rule"] == "gene_name" and E.tpm[:, 0].tolist()[:2] == [1, 2]
    assert np.isnan(E.tpm[2, 0]) and E.missing.tolist() == [False, False, True]
    with pytest.raises(KeyError, match="id column"):
        Expression.from_frame(df, ids="nope")
    with pytest.raises(TypeError):
        Expression.from_frame(42)
    with pytest.raises(ValueError, match="level"):
        Expression(dict(counts=np.ones((1, 1))), ["a"], ids=["x"], level="exon")
    with pytest.raises(ValueError, match="rows"):
        Expression(dict(counts=np.ones((2, 1))), ["a"], ids=["x", "y"], rows=genes, level="gene")


@requires("anndata")
def test_anndata_round_trip(rnaseq_dir, genes):
    E = in_order(Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes,
                                 samples=str(rnaseq_dir / "samples.tsv"), verbose=False))
    a = E.to_anndata()
    assert a.shape == (6, 3) and a.obs["condition"].tolist()[:3] == ["WT"] * 3
    assert {k for k in a.layers.keys() if k is not None} == {"tpm", "length"}     # anndata 0.13 lists a None key
    assert a.var.index.tolist() == ["G1", "G2", "G3"] and a.var["gene_name"].tolist() == ["GENE_A", "GENE_B", "GENE_C"]
    np.testing.assert_allclose(a.X, E.counts.T)
    B = Expression.from_anndata(a, genes=genes, verbose=False)
    for k in E.layers:
        np.testing.assert_allclose(B.layers[k], E.layers[k])
    assert B.samples["batch"].tolist() == E.samples["batch"].tolist() and B.rows is genes
    C = Expression.from_anndata(E.to_anndata("tpm"), verbose=False)         # X was tpm: it comes back as tpm
    np.testing.assert_allclose(C.tpm, E.tpm)
    np.testing.assert_allclose(C.counts, E.counts)
    assert C.rows is None and C.names.tolist() == ["GENE_A", "GENE_B", "GENE_C"] and C.level == "gene"
    import anndata as ad
    foreign = ad.AnnData(X=np.ones((6, 3)), layers={"counts": np.full((6, 3), 2.0)})
    foreign.var_names = ["G1", "G2", "G3"]
    F = Expression.from_anndata(foreign, verbose=False)                       # a foreign 'counts' layer keeps its name
    assert list(F.layers) == ["counts", "X"] and F.counts[0, 0] == 2.0 and F.layers["X"][0, 0] == 1.0
    assert Expression.from_anndata(foreign, layer="counts", verbose=False).counts[0, 0] == 2.0


def test_save_load(rnaseq_dir, genes, tmp_path):
    E = in_order(Expression.make(str(rnaseq_dir / "by_name.csv"), genes=genes, samples=str(rnaseq_dir / "samples.tsv"),
                                 verbose=False))
    E.save(str(tmp_path / "expr"))
    import os
    assert sorted(os.listdir(tmp_path / "expr")) == ["counts.parquet", "features.parquet", "meta.json", "samples.parquet"]
    B = Expression.load(str(tmp_path / "expr"))
    assert B.rows is None and B.level == "gene" and B.ids.tolist() == ["G1", "G2", "G3"]
    np.testing.assert_allclose(B.counts, E.counts, equal_nan=True)
    assert B.samples.equals(E.samples) and B.meta["key_rule"] == "gene_name"
    C = Expression.load(str(tmp_path / "expr"), genes=genes)
    assert C.rows is genes and C.missing.tolist() == [False, False, True]
    np.testing.assert_allclose(C.counts, E.counts, equal_nan=True)
    T = Expression.make(str(rnaseq_dir / "salmon/*/quant.sf"), genes=genes, level="transcript", verbose=False)
    T.save(str(tmp_path / "tx"))
    T2 = Expression.load(str(tmp_path / "tx"), genes=genes)
    assert T2.level == "transcript" and T2.table is genes.transcripts
    np.testing.assert_allclose(T2.tpm, T.tpm)
    np.testing.assert_allclose(T2.to_genes(verbose=False).counts, T.to_genes(verbose=False).counts)
