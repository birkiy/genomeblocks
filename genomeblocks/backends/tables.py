"""Table backends: parsing tab-separated text (BED, BEDPE, narrowPeak, ...)
with polars (the default when installed) or pandas.

Both return the same thing: a dict of numpy columns, positions as int64 and
everything else as strings (object arrays), so callers never touch the frame
library. Leading ``#`` / ``track`` / ``browser`` lines are skipped and
``.gz`` files are read directly.
"""
from __future__ import annotations

import gzip
from typing import Dict, Optional, Sequence

import numpy as np

from . import resolve


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _first_data_line(path, max_lines: int = 50):
    """(header lines before it, tab-separated fields of the first data line; [] when none)."""
    skip = 0
    with _open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                skip += 1
                if skip > max_lines:
                    break
                continue
            return skip, line.rstrip("\n").split("\t")
    return skip, []


def sniff(path, max_lines: int = 50):
    """(header lines to skip, number of columns of the first data line)."""
    skip, fields = _first_data_line(path, max_lines)
    return skip, len(fields)


def first_fields(path, max_lines: int = 50) -> list:
    """The tab-separated fields of the first data line ([] when there is none)."""
    return _first_data_line(path, max_lines)[1]


def _parse_error(path, exc, column=None) -> ValueError:
    where = f" (column {column + 1})" if column is not None else ""
    return ValueError(
        f"{path}: cannot read as a tab-separated BED-like table{where}: {exc}. Header lines must start "
        f"with '#', 'track' or 'browser'; a file with a column-name header is read with "
        f"Loci.from_frame(pandas.read_csv(path, sep='\\t')) or as_loci(path) for .tsv / .csv")


def read_columns(path, columns: Sequence[int], *, ints: Sequence[int] = (), floats: Sequence[int] = (),
                 backend: Optional[str] = None) -> Dict[int, np.ndarray]:
    """Columns ``columns`` (0-based) of a headerless TSV as numpy arrays.

    Columns in ``ints`` become int64, in ``floats`` float64, the rest strings.
    The header lines counted by :func:`sniff` are skipped; a ``#`` inside a
    data field is data (as for bedtools / UCSC).
    """
    skip, ncol = sniff(path)
    columns = [c for c in columns if c < max(ncol, 1)]
    b = resolve("tables", backend)
    if ncol == 0 or not columns:
        return {c: np.zeros(0, np.int64 if c in ints else (np.float64 if c in floats else object))
                for c in columns}
    if b == "polars":
        import polars as pl
        order = sorted(set(columns))                 # polars returns columns in file order
        names = {c: f"c{c}" for c in order}
        try:
            df = pl.read_csv(path, separator="\t", has_header=False, skip_rows=skip, columns=order,
                             new_columns=[names[c] for c in order], infer_schema_length=0, quote_char=None,
                             truncate_ragged_lines=True)
        except Exception as e:                       # noqa: BLE001 — polars has many error classes
            raise _parse_error(path, e) from None
        out = {}
        for c in columns:
            s = df[names[c]]
            if c in ints:
                try:
                    out[c] = s.cast(pl.Int64).to_numpy()
                except Exception as e:               # noqa: BLE001
                    raise _parse_error(path, e, column=c) from None
            elif c in floats:
                out[c] = s.cast(pl.Float64, strict=False).to_numpy()
            else:
                out[c] = s.to_numpy().astype(object)
        return out
    import pandas as pd
    dtypes = {c: (np.int64 if c in ints else str) for c in columns}
    try:
        df = pd.read_csv(path, sep="\t", header=None, skiprows=skip, usecols=list(columns), dtype=dtypes,
                         quoting=3, engine="c")
    except (ValueError, TypeError) as e:
        raise _parse_error(path, e) from None
    out = {}
    for c in columns:
        if c in ints:
            out[c] = df[c].to_numpy(np.int64)
        elif c in floats:
            out[c] = pd.to_numeric(df[c], errors="coerce").to_numpy(np.float64)
        else:
            out[c] = df[c].to_numpy(object)
    return out
