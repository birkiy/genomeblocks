"""The public surface: lazy exports, import cost, the backends registry."""
import subprocess
import sys

import pytest

import genomeblocks as gb


def test_every_export_resolves():
    for name in gb.__all__:
        assert getattr(gb, name) is not None, name
    assert sorted(gb.__all__) == gb.__all__


def test_removed_names_stay_gone():
    for name in ("set_backend", "make_genome", "Pair", "read_bedpe", "columnar", "to_legacy"):
        assert not hasattr(gb, name), name
    with pytest.raises(ImportError):
        __import__("genomeblocks.columnar")


def test_import_is_light():
    code = ("import sys, time; t = time.perf_counter(); import genomeblocks; from genomeblocks import Loci; "
            "dt = time.perf_counter() - t; heavy = [m for m in ('pandas', 'scipy', 'matplotlib', 'polars', "
            "'pyarrow', 'narwhals', 'cooler', 'pybigtools') if m in sys.modules]; print(dt, heavy)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.split()
    assert float(out[0]) < 1.0
    assert out[1:] == ["[]"], out


def test_backends_table_and_resolve():
    df = gb.backends()
    assert set(df["family"]) == {"intervals", "bigwig", "motifs", "fasta", "tables", "graph"}
    assert df.loc[df["default"], "installed"].all()
    assert gb.backends.resolve("intervals") == "genomeblocks"
    assert gb.backends.resolve("intervals", "numpy") == "genomeblocks"        # alias
    with pytest.raises(ValueError, match="unknown intervals backend"):
        gb.backends.resolve("intervals", "nope")
    with pytest.raises(ValueError, match="unknown intervals backend"):
        gb.backends.installed("intervals", "nope")
    with pytest.raises(ValueError, match="unknown backend family"):
        gb.backends.resolve("nope")


def test_missing_backend_names_the_install():
    missing = [b for b in ("cgranges", "bedtools") if not gb.backends.installed("intervals", b)]
    for b in missing:
        with pytest.raises(ImportError, match="install"):
            gb.backends.resolve("intervals", b)
        with pytest.raises(ImportError):
            with gb.use_backend(intervals=b):
                pass


def test_use_backend_scopes_and_restores():
    from genomeblocks.backends import _scoped, resolve
    assert resolve("tables") in ("polars", "pandas")
    with gb.use_backend(tables="pandas"):
        assert resolve("tables") == "pandas"
        with gb.use_backend(tables="auto"):
            assert resolve("tables") == ("polars" if gb.backends.installed("tables", "polars") else "pandas")
        assert resolve("tables") == "pandas"
        assert resolve("tables", "polars" if gb.backends.installed("tables", "polars") else "pandas") != "nope"
    assert "tables" not in _scoped
