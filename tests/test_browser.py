import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

from genomeblocks.browserview import _detect_track_type, _parse_region, browser
from genomeblocks import Loci
from genomeblocks.locus import Locus


def test_parse_region_forms():
    assert _parse_region("chr1:1,000-2,000") == ("chr1", 1000, 2000)
    assert _parse_region(("chr2", 5, 9)) == ("chr2", 5, 9)
    assert _parse_region(Locus("chr3", 10, 20)) == ("chr3", 10, 20)


def test_detect_track_type():
    assert _detect_track_type("a.bw") == "bw"
    assert _detect_track_type("a.narrowPeak") == "narrowPeak"
    assert _detect_track_type("a.bed") == "bed"
    assert _detect_track_type("a.bedpe") == "bedpe"
    assert _detect_track_type(["a.bw", "b.bigwig"]) == "bw"      # rep list -> bw
    assert _detect_track_type(Loci([Locus("chr1", 1, 2)])) == "bed"


class _StubBW:
    """Minimal bigwig handle returning a constant value across bins."""
    def __init__(self, v): self.v = v
    def stats_array(self, chrom, start, end, *, n_bins, stat, missing):
        return np.full(n_bins, self.v, dtype=float)
    def close(self): pass


@pytest.fixture
def stub_bw(monkeypatch):
    maxes = {"a.bw": 2.0, "b.bw": 5.0, "c.bw": 9.0}
    import genomeblocks.signal as sig
    monkeypatch.setattr(sig, "_bw_open", lambda p: _StubBW(maxes[p]))
    return maxes


def test_bigwig_list_is_averaged(stub_bw):
    fig, ax = browser(("chr1", 0, 100), {"reps": ["a.bw", "b.bw"]})
    line = ax["reps"].get_lines()[-1]
    assert np.allclose(line.get_ydata(), 3.5)    # mean(2, 5)


def test_bw_share_uses_group_max(stub_bw):
    tracks = {"AR 0h": "a.bw", "AR 4h": "b.bw", "other": "c.bw"}
    fig, ax = browser(("chr1", 0, 100), tracks, bw_share=[["AR 0h", "AR 4h"]])
    y0 = ax["AR 0h"].get_ylim()[1]
    y1 = ax["AR 4h"].get_ylim()[1]
    assert y0 == pytest.approx(y1)               # shared scale
    assert y0 == pytest.approx(5.0 * 1.05)       # group max * headroom
    assert ax["other"].get_ylim()[1] == pytest.approx(9.0 * 1.05)


def test_bw_ymax_scalar_applies_to_all(stub_bw):
    tracks = {"AR 0h": "a.bw", "AR 4h": "b.bw"}
    fig, ax = browser(("chr1", 0, 100), tracks, bw_ymax=12.0)
    assert ax["AR 0h"].get_ylim()[1] == pytest.approx(12.0)
    assert ax["AR 4h"].get_ylim()[1] == pytest.approx(12.0)
