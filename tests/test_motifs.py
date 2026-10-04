import numpy as np
import pandas as pd
import pytest

from genomeblocks.motifs import bootstrap_enrichment


def _matrix(n, m1, m2, seed):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "M1": np.abs(rng.normal(m1, 0.05, n)),
        "M2": np.abs(rng.normal(m2, 0.05, n)),
    })


def test_bootstrap_enrichment_columns_and_sign():
    ref = _matrix(200, m1=0.1, m2=1.0, seed=0)
    a   = _matrix(80,  m1=2.0, m2=1.0, seed=1)   # M1 enriched in A
    b   = _matrix(80,  m1=0.1, m2=1.0, seed=2)

    enr = bootstrap_enrichment({"A": a, "B": b}, ref=ref,
                               boot=50, sample=40, seed=0, verbose=False)

    for col in ["Factor", "mean_ref", "mean_A", "LFC_A", "mean_B", "LFC_B", "LFC"]:
        assert col in enr.columns

    row = enr.set_index("Factor").loc["M1"]
    assert row["LFC_A"] > 0                       # M1 up in A vs pool
    assert row["LFC"] > 0                         # and up in A vs B (LFC_A - LFC_B)


def test_bootstrap_enrichment_single_group_no_diff_column():
    ref = _matrix(100, 0.5, 0.5, seed=3)
    a   = _matrix(50, 0.5, 0.5, seed=4)
    enr = bootstrap_enrichment({"A": a}, ref=ref, boot=20, sample=20,
                               seed=0, verbose=False)
    assert "LFC_A" in enr.columns
    assert "LFC" not in enr.columns               # only emitted for exactly 2 groups


def test_block_scan_matches_per_window_scan():
    # Windows are scanned as one concatenated block; per-window counts must
    # equal scanning each window alone (no hit may straddle two windows), and
    # windows lightmotif can't stripe count 0.
    lightmotif = pytest.importorskip("lightmotif")
    import random
    from genomeblocks.motifs import _Block
    rng = random.Random(3)
    core = "TGACTCA"
    seqs = []
    for i in range(60):
        s = "".join(rng.choice("ACGT") for _ in range(100))
        if i % 3 == 0:
            s = s[:40] + core + s[47:]               # a hit inside the window
        if i % 5 == 0:
            s = s[:-4] + core[:4]                    # half a hit at the right edge...
        if i % 5 == 1:
            s = core[4:] + s[3:]                     # ...completed by the next window
        seqs.append(s)
    seqs[7] = seqs[7][:50] + "X" + seqs[7][51:]      # invalid character
    motif = lightmotif.create([core] * 50)          # simple consensus PWM
    pssm = motif.counts.normalize(0.1).log_odds()

    got = _Block(seqs).counts(pssm, len(core), threshold=5.0)
    want = []
    for s in seqs:
        try:
            st = lightmotif.stripe(s)
        except ValueError:
            want.append(0)
            continue
        want.append(sum(1 for _ in lightmotif.scan(pssm, st, threshold=5.0)))
    assert got.tolist() == want
    assert got[7] == 0 and sum(want) > 0
