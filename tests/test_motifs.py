import numpy as np
import pandas as pd

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
