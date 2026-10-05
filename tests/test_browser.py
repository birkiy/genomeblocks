"""The three viewers: the matplotlib browser, the IGV page and the one-file View."""
import base64
import gzip

import numpy as np
import pytest

import genomeblocks as gb
from genomeblocks.browserview import _detect_track_type
from genomeblocks.view import View


def test_detect_track_types(genes, pairs, cre, bw_path, bedpe_path, bed_path):
    assert _detect_track_type(genes) == "genes" and _detect_track_type(pairs) == "bedpe"
    assert _detect_track_type(bedpe_path) == "bedpe" and _detect_track_type(bw_path) == "bw"
    assert _detect_track_type([bw_path, bw_path]) == "bw" and _detect_track_type("x.bam") == "bam"
    assert _detect_track_type(cre) == "bed" and _detect_track_type(cre.to_pandas()) == "bed"
    assert _detect_track_type(bed_path) == "bed" and _detect_track_type("chr1:1-2") == "bed"
    import pyBigWig
    assert _detect_track_type(pyBigWig.open(bw_path)) == "bw"


def test_browser_draws_every_track_kind(genes, pairs, cre, bw_path, bw2_path, bedpe_path, bed_path):
    fig, axes = gb.browser("chr1:0-12 kb", {"ATAC": bw_path, "ATAC x2": [bw_path, bw2_path], "CRE": cre,
                                            "frame": cre.to_pandas(), "bed": bed_path, "loops": pairs,
                                            "loops file": bedpe_path, "genes": genes,
                                            "regions": ["chr1:500-700", "chr1:3000-3500"]},
                           bw_share=[["ATAC", "ATAC x2"]], backend="python")
    assert set(axes) >= {"ATAC", "genes", "_axis"} and len(fig.axes) == len(axes)
    assert axes["ATAC"].get_ylim() == axes["ATAC x2"].get_ylim()        # shared scale
    from matplotlib.collections import LineCollection, PatchCollection
    n_rects = sum(len(c.get_paths()) for c in axes["CRE"].collections if isinstance(c, PatchCollection))
    assert n_rects == 5                                                  # the five chr1 CREs in view
    arcs = [c for c in axes["loops"].collections if isinstance(c, LineCollection)]
    assert sum(len(c.get_paths()) for c in arcs) == 2                     # the two cis chr1 loops
    assert sum(len(c.get_paths()) for c in axes["regions"].collections) == 2
    gene_boxes = sum(len(c.get_paths()) for c in axes["genes"].collections if isinstance(c, PatchCollection))
    assert gene_boxes == 6 + 3                                           # exons (T1 3 + T1b 1 + T2 2) + CDS (3)
    assert len(axes["genes"].texts) == 2                                 # one label per gene on chr1
    fig, axes = gb.browser(("chr1", 0, 12_000), {"genes": genes}, genes_max_transcripts=1, bw_ymax=5)
    assert "genes" in axes
    with pytest.raises(ValueError, match="end"):
        gb.browser("chr1:100-50", {"CRE": cre})


def test_browser_bam_track(tmp_path, cre):
    pysam = pytest.importorskip("pysam")
    hdr = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 20_000, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(tmp_path / "t.bam"), "wb", header=hdr) as f:
        for i in range(20):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag, a.reference_id = f"r{i}", "ACGT" * 10, 0, 0
            a.reference_start, a.mapping_quality, a.cigar = 1000 + i * 3, 60, [(0, 40)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 40)
            f.write(a)
    pysam.index(str(tmp_path / "t.bam"))
    (tmp_path / "ref.fa").write_text(">chr1\n" + "ACGT" * 5000 + "\n")
    pysam.faidx(str(tmp_path / "ref.fa"))
    assert gb.coverage(str(tmp_path / "t.bam"), "chr1:1,000-1,100").max() == 14
    with pytest.raises(ValueError, match="reference"):
        gb.browser("chr1:900-1200", {"bam": str(tmp_path / "t.bam")})
    fig, axes = gb.browser("chr1:950-1150", {"bam": str(tmp_path / "t.bam"), "CRE": cre}, reference=str(tmp_path / "ref.fa"))
    assert "_sequence" in axes


def test_igv_html(tmp_path, genes, cre, arch, bw_path):
    out = tmp_path / "share.html"
    sizes = gb.igv_html(str(out), regions=["chr1:0.0-0.012 Mb chr2:1-6 kb", ("chr1", 0, 5000)],
                        loci={"CRE": cre.to_pandas()}, genes=genes, signal={"ATAC": bw_path}, architecture=arch)
    assert set(sizes) == {"ATAC", "CRE", "genes", "loops (n)", "total"}
    page = out.read_text()
    assert "chr1:1-12000 chr2:1001-6000" in page                       # split view, units parsed, 1-based for igv.js
    assert page.count("data:application/gzip") == 4
    gz = page.split("genes")[-1]
    assert gz                                                           # the genes track is to_bed12 (0-based starts)


def _unpack(p):
    return np.frombuffer(gzip.decompress(base64.b64decode(p["b"])), dtype=p["t"])


def test_view_payload_and_html(tmp_path, genes, cre, arch, bw_path):
    v = (View(arch, genes=genes, samples=["S1"]).signal("ATAC", {"S1": bw_path})
         .intervals("peaks", {"S1": cre.to_pandas()}).cres().loops().genes())
    v.region("GENE_B", gene="GENE_B").region("split", "chr1:0-5 kb chr2:1-6 kb").mark("m", "chr1", 1000)
    D, parts = v._payload()
    assert np.cumsum(_unpack(D["genes"]["start"])).tolist() == [1000, 10_000, 20_000]     # 0-based gene starts
    assert _unpack(D["genes"]["tss"]).tolist() == [1000, 10_999, 20_000]                  # TSS = end - 1 on '-'
    assert D["regions"][0]["loci"] == [["chr1", 10_999 - 500_000 if 10_999 > 500_000 else 0, 10_999 + 500_000]]
    assert D["regions"][1]["loci"] == [["chr1", 0, 5000], ["chr2", 1000, 6000]]
    assert len(D["tracks"]) == 5 and D["samples"][0]["name"] == "S1"
    sizes = v.save(str(tmp_path / "v.html"))
    assert sizes["total"] > 1000 and (tmp_path / "v.html").read_text().startswith("<!doctype html>")
    assert "<iframe" in v._repr_html_()
    assert len(View(cre=cre.to_pandas()).to_html()[0]) > 1000
    with pytest.raises(ValueError, match="Architecture or cre"):
        View()
