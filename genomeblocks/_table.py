"""What every genomeblocks table offers the rest of the ecosystem.

A class that can build an Arrow table of itself (``to_arrow``) gets, from
this mixin, the protocols other libraries look for — so polars, pandas,
DuckDB, seaborn, plotly, altair, narwhals and friends take the object as is:

    __arrow_c_stream__   Arrow PyCapsule stream (polars, pyarrow, duckdb, pandas >= 2.2)
    __dataframe__        dataframe interchange protocol (seaborn, pandas)
    __narwhals_dataframe__  a full narwhals frame backed by our Arrow table, which
                         narwhals-based tools (altair, plotly, marimo, ...) prefer
                         over the interchange protocol
    to_polars()          via Arrow (zero-copy for the numeric columns)
    shape / columns      like a DataFrame
"""
from __future__ import annotations


class TableMixin:
    """Protocols and conversions derived from ``to_arrow()``."""

    def __arrow_c_stream__(self, requested_schema=None):
        return self.to_arrow().__arrow_c_stream__(requested_schema)

    def __dataframe__(self, nan_as_null: bool = False, allow_copy: bool = True):
        return self.to_arrow().__dataframe__(nan_as_null=nan_as_null, allow_copy=allow_copy)

    def __narwhals_dataframe__(self):
        import narwhals as nw
        return nw.from_native(self.to_arrow(), eager_only=True)._compliant_frame

    def to_polars(self):
        """polars DataFrame (through Arrow: numeric columns are not copied)."""
        import polars as pl
        return pl.from_arrow(self.to_arrow())

    @property
    def shape(self):
        return (len(self), len(self.columns))
