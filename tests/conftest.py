"""Shared fixtures for the genomeblocks 2.0 test suite.

Everything is synthetic — tiny files written into ``tmp_path`` and small
in-memory tables — so the suite runs in seconds and needs no genome, bigWig
or ChIP-Atlas data. Backend parity tests run over the engines installed in
the environment (``installed_backends``); the conda-only ones (graph-tool,
cgranges, bedtools) are skipped when absent.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

import genomeblocks as gb  # noqa: E402
from genomeblocks import Loci, as_loci  # noqa: E402


def installed_backends(family):
    """The backends of ``family`` that can run here."""
    return [b for b in gb.backends.families()[family] if gb.backends.installed(family, b)]


def requires(module):
    return pytest.mark.skipif(__import__("importlib").util.find_spec(module) is None,
                              reason=f"{module} not installed")


# ── tiny on-disk fixtures ────────────────────────────────────────────────────

GTF = """\
chr1\tx\tgene\t1001\t5000\t.\t+\t.\tgene_id "G1"; gene_name "GENE_A"; gene_type "protein_coding";
chr1\tx\ttranscript\t1001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE_A"; gene_type "protein_coding";
chr1\tx\texon\t1001\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
chr1\tx\texon\t2001\t2500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";
chr1\tx\texon\t4001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "3";
chr1\tx\tCDS\t1101\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
chr1\tx\tCDS\t2001\t2500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";
chr1\tx\tCDS\t4001\t4300\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "3";
chr1\tx\tfive_prime_utr\t1001\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
chr1\tx\tthree_prime_utr\t4301\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "3";
chr1\tx\ttranscript\t2001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1b"; gene_name "GENE_A"; gene_type "protein_coding";
chr1\tx\texon\t2001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1b"; exon_number "1";
chr1\tx\tgene\t10001\t11000\t.\t-\t.\tgene_id "G2"; gene_name "GENE_B"; gene_type "lncRNA";
chr1\tx\ttranscript\t10001\t11000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE_B"; gene_type "lncRNA";
chr1\tx\texon\t10001\t10400\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; exon_number "2";
chr1\tx\texon\t10601\t11000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; exon_number "1";
chr2\tx\tgene\t20001\t30000\t.\t+\t.\tgene_id "G3"; gene_name "GENE_C"; gene_type "protein_coding";
chr2\tx\ttranscript\t20001\t30000\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; gene_name "GENE_C"; gene_type "protein_coding";
chr2\tx\texon\t20001\t30000\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; exon_number "1";
chr2\tx\tCDS\t20101\t29900\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; exon_number "1";
"""

BED6 = """\
# a comment
track name=peaks
chr1\t900\t1100\tp1\t10\t+
chr1\t1900\t2100\tp2\t20\t-
chr1\t4900\t5100\tp3\t30\t+
chr1\t9950\t10050\tp4\t40\t.
chr1\t10900\t11100\tp5\t50\t-
chr2\t500\t600\tp6\t60\t+
chr2\t5000\t5100\tp7\t70\t+
"""

BEDPE = """\
chr1\t900\t1100\tchr1\t4900\t5100\tl1\t5\t+\t-
chr1\t1900\t2100\tchr1\t10900\t11100\tl2\t3\t+\t+
chr1\t900\t1100\tchr2\t500\t600\tl3\t1\t-\t+
chr2\t500\t600\tchr2\t5000\t5100\tl4\t2\t+\t+
"""

JASPAR = """\
>M1 TFA
10 0 0 0
0 10 0 0
0 0 10 0
0 0 0 10
>M2 TFB
0 0 10 10
0 0 0 0
10 10 0 0
0 0 0 0
"""                                                   # rows A C G T: M1 = ACGT, M2 = GGAA

CHROM_SIZES = {"chr1": 20_000, "chr2": 8_000}


def _write_bigwig(pybigtools, path, recs):
    """Write and read back: a truncated file (disk full, an interrupted writer)
    would otherwise surface far away as a BBIReadError inside signal()."""
    pybigtools.open(str(path), "w").write(CHROM_SIZES, recs)
    h = pybigtools.open(str(path), "r")
    try:
        for chrom in h.chroms():
            n = sum(1 for _ in h.records(chrom))
            if n != sum(r[0] == chrom for r in recs):
                raise RuntimeError(f"{path} is incomplete after writing ({n} records on {chrom})")
    except Exception as exc:
        raise RuntimeError(f"bigWig fixture {path} did not write completely (disk full?): {exc}") from exc


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    """All the synthetic files, written once per session."""
    d = tmp_path_factory.mktemp("data")
    (d / "genes.gtf").write_text(GTF)
    (d / "peaks.bed").write_text(BED6)
    (d / "loops.bedpe").write_text(BEDPE)
    (d / "motifs.jaspar").write_text(JASPAR)
    (d / "genome.chrom.sizes").write_text("".join(f"{c}\t{n}\n" for c, n in CHROM_SIZES.items()))
    # a bigWig with runs and gaps on chr1, nothing on chr2
    import pybigtools
    rng = np.random.default_rng(0)
    recs = [("chr1", s, s + 50, float(v)) for s, v in zip(range(0, 19_900, 100), rng.uniform(1, 10, 199))]
    _write_bigwig(pybigtools, d / "signal.bw", recs)
    recs2 = [("chr1", s, s + 50, float(v)) for s, v in zip(range(0, 19_900, 100), rng.uniform(1, 10, 199))]
    recs2 += [("chr2", s, s + 20, 2.0) for s in range(0, 7_900, 100)]
    _write_bigwig(pybigtools, d / "signal2.bw", recs2)
    # a FASTA with planted motif sites (M1 = ACGT, M2 = GGAA), 60-column lines
    seq = "".join(rng.choice(list("ACGT"), 20_000))
    seq = seq[:1000] + "ACGT" * 3 + seq[1012:5000] + "GGAA" * 3 + seq[5012:]
    chr2 = "".join(rng.choice(list("ACGT"), 8_000))
    with open(d / "genome.fa", "w") as f:
        for name, s in (("chr1", seq), ("chr2", chr2)):
            f.write(f">{name}\n")
            for k in range(0, len(s), 60):
                f.write(s[k:k + 60] + "\n")
    return d


@pytest.fixture
def gtf_path(data_dir):
    return str(data_dir / "genes.gtf")


@pytest.fixture
def bed_path(data_dir):
    return str(data_dir / "peaks.bed")


@pytest.fixture
def bedpe_path(data_dir):
    return str(data_dir / "loops.bedpe")


@pytest.fixture
def bw_path(data_dir):
    return str(data_dir / "signal.bw")


@pytest.fixture
def bw2_path(data_dir):
    return str(data_dir / "signal2.bw")


@pytest.fixture
def fasta_path(data_dir):
    return str(data_dir / "genome.fa")


@pytest.fixture
def jaspar_path(data_dir):
    return str(data_dir / "motifs.jaspar")


@pytest.fixture
def sizes_path(data_dir):
    return str(data_dir / "genome.chrom.sizes")


# ── in-memory tables ──────────────────────────────────────────────────────────

@pytest.fixture
def cre():
    """Seven CREs: five on chr1 (around the fixture genes), two on chr2."""
    return as_loci([("chr1", 900, 1100), ("chr1", 1900, 2100), ("chr1", 4900, 5100), ("chr1", 9950, 10050),
                    ("chr1", 10900, 11100), ("chr2", 500, 600), ("chr2", 5000, 5100)])


@pytest.fixture
def genes(gtf_path):
    from genomeblocks import Genes
    return Genes.make(gtf_path)


@pytest.fixture
def pairs(bedpe_path):
    from genomeblocks import Pairs
    return Pairs.make(bedpe_path)


@pytest.fixture
def arch(cre, pairs):
    """An Architecture over `cre` with weights set and normalised."""
    from genomeblocks import Architecture
    A = Architecture.make(cre, pairs, r=100, verbose=False)
    A.ep["w"][:] = [5.0, 3.0, 2.0, 1.0][:A.n_links]
    return A.normalize(verbose=False)


def random_loci(rng, n, chroms=("chr1", "chr2", "chr10"), span=10_000, zero_length=0.1, genome=None):
    """Random, unsorted intervals (some zero-length, some duplicated) for parity tests."""
    c = rng.choice(chroms, n)
    s = rng.integers(0, span, n)
    ln = rng.integers(1, 300, n)
    ln[rng.random(n) < zero_length] = 0
    if n > 4:
        s[-2:], ln[-2:] = s[:2], ln[:2]                     # duplicates
    return Loci.from_frame({"chrom": c, "start": s, "end": s + ln}, genome=genome)


def brute_pairs(q, r):
    """Overlapping (query row, reference row) pairs by the half-open rule, sorted."""
    qc, rc = q.chroms, r.chroms
    out = [(i, j) for i in range(len(q)) for j in range(len(r))
           if qc[i] == rc[j] and q.starts[i] < r.ends[j] and r.starts[j] < q.ends[i]]
    return sorted(out)
