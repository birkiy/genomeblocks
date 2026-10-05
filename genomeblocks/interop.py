"""Interoperability: genomeblocks tables in and out of the Python ecosystem.

Anything interval-like goes in through :func:`as_loci` — which every
genomeblocks function uses on its inputs — and every table converts back out:

    ===================  =====================================  ===================================
    library / format     in                                     out
    ===================  =====================================  ===================================
    BED / narrowPeak     ``Loci.make(path)``, ``as_loci(path)``  ``L.to_bed(path)``
    parquet              ``Loci.load(path)``                     ``L.save(path)``
    pandas               ``Loci.from_frame(df)``                 ``L.to_pandas()``
    polars               ``Loci.from_frame(df)``                 ``L.to_polars()``
    pyarrow              ``Loci.from_frame(table)``              ``L.to_arrow()``
    bioframe             ``Loci.from_frame(df)``                 ``L.to_bioframe()``
    pyranges             ``Loci.from_pyranges(gr)``              ``L.to_pyranges()``
    pybedtools           ``Loci.from_bedtool(bt)``               ``L.to_bedtool()``
    cgranges             —                                       ``L.to_cgranges()``
    AnnData (scATAC)     ``Loci.from_anndata(adata)``            ``L.to_anndata(X)``
    Biopython            FASTA handles everywhere                ``L.to_seqrecords(fasta)``
    xarray               —                                       ``cube_to_xarray(S, L)``
    ===================  =====================================  ===================================

Frames are read as 0-based, half-open intervals (BED / bioframe / pyranges
convention). Row order is always kept, so a frame's rows, an AnnData's
``var`` and the resulting Loci line up one to one.
"""
from __future__ import annotations

import os
import re
from typing import Optional

import numpy as np

from .genome import Genome
from .locus import Locus

_CHROM = ("chrom", "chromosome", "chr", "seqnames", "seqname", "seqid", "contig", "#chrom",
          "chrom_name", "chromosome_name", "ref", "reference")
_START = ("start", "chromstart", "begin", "start_position", "txstart")
_END = ("end", "chromend", "stop", "end_position", "txend")
_STRAND = ("strand",)
_NAME_RE = re.compile(r"^(?P<chrom>.+?)[:_-](?P<start>\d+)[-_:](?P<end>\d+)$")
_SCODE = {".": 0, "+": 1, "-": 2, "*": 0, "1": 1, "-1": 2, "": 0}
_EXPECTED = ("expected columns for the chromosome, start and end, e.g. chrom/start/end (bioframe), "
             "Chromosome/Start/End (pyranges), seqnames/start/end or chr/chromStart/chromEnd — "
             "or pass chrom=, start=, end= to name them")


def _mod(x) -> str:
    return (type(x).__module__ or "").split(".")[0]


def _find(columns, wanted, given):
    if given is not None:
        if given not in columns:
            raise KeyError(f"column {given!r} not found; columns are {list(columns)}")
        return given
    low = {str(c).lower(): c for c in columns}
    for w in wanted:
        if w in low:
            return low[w]
    return None


def _as_pandas(df):
    """pandas DataFrame of any table genomeblocks accepts (used where pandas is the
    natural tool, e.g. GTF-like frames); goes through :func:`frame`."""
    import pandas as pd
    if isinstance(df, pd.DataFrame) and type(df) is pd.DataFrame:
        return df
    f = frame(df)
    if f is None:
        raise TypeError(f"expected a table (pandas, polars, pyarrow, ...), got {type(df).__name__}")
    return f.to_pandas()


# ── the boundary: any table in, one representation out ────────────────────

def frame(obj):
    """An eager narwhals DataFrame for anything table-like, else None.

    Takes pandas, polars (eager or lazy), pyarrow, modin, cuDF, DuckDB and the
    other frames narwhals knows; PyRanges and pybedtools objects; dicts of
    equal-length columns; numpy structured arrays; and any object that speaks
    the Arrow C stream (``__arrow_c_stream__``) or dataframe interchange
    (``__dataframe__``) protocol.
    """
    import narwhals as nw
    m = _mod(obj)
    if m == "pyranges":
        import pandas as pd
        obj = obj.df if hasattr(obj, "df") else pd.DataFrame(obj)
    elif m == "pybedtools":
        obj = obj.to_dataframe(disable_auto_names=False)
    elif isinstance(obj, dict):
        import pandas as pd
        scalars = [k for k, v in obj.items() if isinstance(v, str) or not hasattr(v, "__len__")]
        if scalars:
            raise TypeError(f"a dict of columns needs an array per key, but {scalars[:3]} hold scalars; "
                            f"expected e.g. {{'chrom': [...], 'start': [...], 'end': [...]}}")
        lens = {k: len(v) for k, v in obj.items()}
        if len(set(lens.values())) > 1:
            raise ValueError(f"the columns have different lengths: {lens}")
        obj = pd.DataFrame({k: np.asarray(v) for k, v in obj.items()})
    elif isinstance(obj, np.ndarray) and obj.dtype.names:
        import pandas as pd
        obj = pd.DataFrame({k: obj[k] for k in obj.dtype.names})
    try:
        df = nw.from_native(obj)
    except TypeError:
        if hasattr(obj, "__arrow_c_stream__"):
            import pyarrow as pa
            df = nw.from_native(pa.table(obj))
        elif hasattr(obj, "__dataframe__"):
            import pyarrow.interchange
            df = nw.from_native(pyarrow.interchange.from_dataframe(obj))
        else:
            return None
    if isinstance(df, nw.LazyFrame):
        df = df.collect()
    if not isinstance(df, nw.DataFrame):
        return None
    return df


def _col(df, name) -> np.ndarray:
    """One column as numpy (categoricals and dictionaries come back as values)."""
    import narwhals as nw
    s = df.get_column(name)
    if s.dtype == nw.Categorical or s.dtype == nw.Enum:
        s = s.cast(nw.String)
    return s.to_numpy()


def _strand_codes(values) -> np.ndarray:
    v = np.asarray(values).astype(str)
    out = np.zeros(len(v), np.int8)
    out[(v == "+") | (v == "1")] = 1
    out[(v == "-") | (v == "-1")] = 2
    return out


def loci_from_frame(df, chrom=None, start=None, end=None, strand=None, *, keep=True,
                    genome: Optional[Genome] = None):
    """Loci from any table :func:`frame` accepts; row order kept.

    Columns are found by name (case-insensitive: chrom / chromosome / chr /
    seqnames / Chromosome ..., start / chromStart / Start, end / chromEnd /
    End, strand / Strand), or given with ``chrom=``, ``start=``, ``end=``,
    ``strand=``. A frame without such names but with a text column and two
    integer columns first (``read_csv(header=None)`` on a BED) is read by
    position. ``keep=True`` keeps every other column (a list keeps those)."""
    import narwhals as nw
    from .loci import Loci
    f = df if isinstance(df, nw.DataFrame) else frame(df)
    if f is None:
        raise TypeError(f"expected a table (pandas / polars / pyarrow / dict of columns), "
                        f"got {type(df).__name__}")
    g = genome if genome is not None else Genome()
    cols = list(f.columns)
    c, s, e = _find(cols, _CHROM, chrom), _find(cols, _START, start), _find(cols, _END, end)
    if c is None or s is None or e is None:
        sch = f.schema
        first = cols[:3]
        if (len(first) == 3 and (sch[first[0]] == nw.String or sch[first[0]] == nw.Categorical
                                 or sch[first[0]] == nw.Object)
                and sch[first[1]].is_integer() and sch[first[2]].is_integer()):
            c, s, e = first
        else:
            raise ValueError(f"{_EXPECTED}; got columns {cols}")
    for k in (c, s, e):
        col = f.get_column(k)
        n = col.null_count()
        if not n and col.dtype.is_float():
            n = int(col.is_nan().sum())
        if n:
            raise ValueError(f"column {k!r} has {n} missing value(s) in {len(f)} rows; chrom, start "
                             f"and end must be complete — drop or fill those rows first")
    d = _find(cols, _STRAND, strand)
    used = {c, s, e} | ({d} if d is not None else set())
    if keep is True:
        extra = [k for k in cols if k not in used]
    elif keep:
        extra = [k for k in keep if k in cols]
    else:
        extra = []
    codes = g.encode(_col(f, c).astype(object)) if len(f) else np.zeros(0, np.int32)
    return Loci(codes, _col(f, s).astype(np.int64), _col(f, e).astype(np.int64),
                _strand_codes(_col(f, d)) if d is not None else None,
                genome=g, cols={str(k): _col(f, k) for k in extra})


def loci_from_pyranges(gr, **kw):
    """Loci from a PyRanges (0.x or 1.x)."""
    return loci_from_frame(gr, **kw)


def loci_from_bedtool(bt, **kw):
    """Loci from a pybedtools BedTool (any BED-like file it wraps)."""
    return loci_from_frame(bt, **kw)


def loci_from_anndata(adata, axis: str = "var", *, genome: Optional[Genome] = None, keep=True):
    """The regions of an AnnData axis ('var' for scATAC peaks), in its order:
    chrom/start/end columns when present, else parsed from the names."""
    import pandas as pd
    fr = adata.var if axis == "var" else adata.obs
    cols = list(fr.columns)
    if all(_find(cols, w, None) is not None for w in (_CHROM, _START, _END)):
        return loci_from_frame(fr.reset_index(drop=True), keep=keep, genome=genome)
    names = pd.Series(fr.index.astype(str))
    parts = names.str.extract(_NAME_RE)
    if parts["chrom"].isna().any():
        bad = names[parts["chrom"].isna()].iloc[0]
        raise ValueError(f"cannot read a region from the {axis} name {bad!r}: expected chr1:100-200, "
                         f"chr1-100-200 or chr1_100_200 names, or chrom/start/end columns in .{axis}")
    data = {"chrom": parts["chrom"].to_numpy(object), "start": parts["start"].astype(np.int64).to_numpy(),
            "end": parts["end"].astype(np.int64).to_numpy()}
    for k in (cols if keep is True else (keep or [])):
        data[k] = fr[k].to_numpy()
    return loci_from_frame(data, keep=True, genome=genome)


def _read_table_file(path: str, *, genome=None, keep=True):
    """CSV / TSV with a header row (a UCSC ``#chrom  chromStart ...`` header
    counts), or a header-less BED-like .txt / .tsv whose columns then get the
    BED names (name, score, strand, ...). Leading ``#`` / track / browser /
    blank lines are skipped. Row order is kept."""
    import pandas as pd
    from .backends.tables import _open_text
    from .loci import _bed_names
    low = path.lower().removesuffix(".gz")
    sep = "," if low.endswith(".csv") else "\t"
    skip, commented, first = 0, None, ""
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "browser")):
                skip += 1
            elif line.startswith("#"):
                skip += 1
                commented = line                       # a '#chrom ...' header, maybe
            else:
                first = line
                break
    fields = first.rstrip("\n").split(sep)
    if not fields or len(fields) < 3:
        if not first:
            return loci_from_frame({"chrom": np.zeros(0, object), "start": np.zeros(0, np.int64),
                                    "end": np.zeros(0, np.int64)}, genome=genome)
        raise ValueError(f"{path}: expected at least 3 {'comma' if sep == ',' else 'tab'}-separated "
                         f"columns (chrom, start, end), found {len(fields)}")
    headerless = fields[1].strip().lstrip("-").isdigit() and fields[2].strip().lstrip("-").isdigit()
    if not headerless:
        df = pd.read_csv(path, sep=sep, header=0, skiprows=skip)
    else:
        names = None
        if commented is not None:
            hf = [h.strip() for h in commented.lstrip("#").rstrip("\n").split(sep)]
            if len(hf) == len(fields) and not hf[1].isdigit():
                names = hf
        if names is None:
            names = _bed_names(low, len(fields))
        df = pd.read_csv(path, sep=sep, header=None, skiprows=skip, names=names)
    return loci_from_frame(df, keep=keep, genome=genome)


def as_loci(x, *, genome: Optional[Genome] = None, keep=True):
    """Anything interval-like as a Loci — the one entry point every genomeblocks
    function uses on its inputs.

    Accepts a Loci; a path (BED / narrowPeak / broadPeak / bedGraph, CSV or TSV
    with a header, parquet; ``.gz`` too — rows stay in file order); a region
    string (``'chr1:1-2,000'``);
    a Locus; any table :func:`frame` takes (pandas, polars, pyarrow, bioframe,
    PyRanges, BedTool, dict of columns, structured array, Arrow / interchange
    protocol objects); an AnnData (its ``var``); a ``(chroms, starts, ends)``
    tuple of arrays; or a list of records (Locus, tuples, region strings) or of
    any of the above (concatenated in order).
    """
    from .loci import Loci
    if isinstance(x, Loci):
        if genome is None or x.genome is genome:
            return x
        return Loci(genome=genome)._check(x)
    if isinstance(x, Locus):
        return Loci.from_records([x], genome=genome)
    if isinstance(x, (str, os.PathLike)):
        p = str(x)
        if not os.path.exists(p):
            if ":" in p:
                return Loci.from_records([p], genome=genome)
            raise FileNotFoundError(f"no such file: {p!r} (and not a region like 'chr1:1,000-2,000')")
        low = p.lower().removesuffix(".gz")
        if low.endswith((".parquet", ".pq")):
            return Loci.load(p, genome=genome)
        if low.endswith((".csv", ".tsv", ".txt")):
            return _read_table_file(p, genome=genome, keep=keep)
        return Loci.make(p, genome=genome, keep=keep, sort=False)   # file order, like every input
    if _mod(x) == "anndata":
        return loci_from_anndata(x, genome=genome, keep=keep)
    if _mod(x) in ("pandas", "polars", "pyarrow") and frame(x) is None:
        raise TypeError(f"cannot read intervals from a {type(x).__name__}: pass a DataFrame / Table "
                        f"with chrom, start and end columns")
    if isinstance(x, tuple) and len(x) in (3, 4) and all(hasattr(v, "__len__") and not isinstance(v, str)
                                                         for v in x):
        data = dict(zip(("chrom", "start", "end", "strand"), x))
        return loci_from_frame(data, genome=genome)
    f = frame(x) if not isinstance(x, (list, tuple)) else None
    if f is not None:
        return loci_from_frame(f, keep=keep, genome=genome)
    if isinstance(x, tuple) and len(x) in (3, 4) and isinstance(x[0], str):
        return Loci.from_records([x], genome=genome)
    if hasattr(x, "__iter__"):
        items = list(x)
        if not items:
            return Loci(genome=genome)
        if all(isinstance(i, (Locus, tuple, list)) or (isinstance(i, str) and ":" in i and not os.path.exists(i))
               for i in items):
            return Loci.from_records(items, genome=genome)
        parts = [as_loci(i, genome=genome, keep=keep) for i in items]
        out = parts[0]
        for q in parts[1:]:
            out = out + q
        return out
    raise TypeError(f"cannot read intervals from {type(x).__name__}: pass a Loci, a BED / CSV / parquet "
                    f"path, a region string, a pandas / polars / pyarrow frame, or a list of regions")


# ── out of Loci ────────────────────────────────────────────────────────────

def _chrom_categorical(L):
    import pandas as pd
    used = np.unique(L.codes) if len(L) else np.zeros(0, np.int32)
    used = used[np.argsort(L.genome.rank[used])] if len(used) else used
    pos = np.full(len(L.genome), -1, np.int64)
    pos[used] = np.arange(len(used))
    return pd.Categorical.from_codes(pos[L.codes] if len(L) else np.zeros(0, np.int64),
                                     [L.genome.names[c] for c in used])


def loci_to_pandas(L, uid: bool = False):
    import pandas as pd
    d = {"chrom": _chrom_categorical(L), "start": L.starts, "end": L.ends,
         "strand": pd.Categorical.from_codes(L.strands, [".", "+", "-"])}
    for k, v in L.cols.items():
        d[k] = v
    if uid:
        d["uid"] = L.uid
    return pd.DataFrame(d)


def loci_to_bioframe(L):
    """bioframe expects chrom / start / end; strand and the extra columns follow."""
    import pandas as pd
    d = {"chrom": L.chroms.astype(str), "start": L.starts, "end": L.ends, "strand": L.strand.astype(str)}
    d.update(L.cols)
    return pd.DataFrame(d)


def loci_to_arrow(L):
    import pyarrow as pa
    cat = _chrom_categorical(L)
    chrom = pa.DictionaryArray.from_arrays(pa.array(np.asarray(cat.codes, np.int32)),
                                           pa.array([str(c) for c in cat.categories], pa.string()))
    strand = pa.DictionaryArray.from_arrays(pa.array(L.strands.astype(np.int8)), pa.array([".", "+", "-"]))
    arrs = {"chrom": chrom, "start": pa.array(L.starts), "end": pa.array(L.ends), "strand": strand}
    for k, v in L.cols.items():
        v = np.asarray(v)
        if v.dtype == object:
            import pandas as pd
            lst = v.tolist()
            na = pd.isna(v)
            if na.any():                                # NaN / None / pd.NA -> Arrow null
                lst = [None if m else x for x, m in zip(lst, na.tolist())]
            arrs[k] = pa.array(lst)
        else:
            arrs[k] = pa.array(v)
    return pa.table(arrs)


def loci_to_polars(L):
    import polars as pl
    return pl.from_arrow(loci_to_arrow(L))


def loci_to_pyranges(L):
    """PyRanges (Chromosome / Start / End / Strand and the extra columns).

    The Strand column keeps '.' rows (pyranges treats such an object as
    unstranded but keeps the values). pyranges 0.x stores rows per chromosome,
    so it returns them grouped by chromosome: sort the Loci first when the
    row order matters."""
    import pandas as pd
    import pyranges as pr
    d = {"Chromosome": _chrom_categorical(L), "Start": L.starts, "End": L.ends,
         "Strand": L.strand.astype(str)}
    df = pd.DataFrame(d)
    for k, v in L.cols.items():
        df[k] = v
    return pr.PyRanges(df)


def loci_to_bedtool(L):
    """pybedtools BedTool (BED6: name = cols['name'] or the uid; score = cols['score'] or 0)."""
    import pandas as pd
    import pybedtools
    df = pd.DataFrame({"chrom": L.chroms.astype(str), "start": L.starts, "end": L.ends,
                       "name": L.cols["name"] if "name" in L.cols else L.uid,
                       "score": L.cols["score"] if "score" in L.cols else np.zeros(len(L), int),
                       "strand": L.strand.astype(str)})
    return pybedtools.BedTool.from_dataframe(df)


def loci_to_cgranges(L):
    import cgranges
    ix = cgranges.cgranges()
    names = L.genome.names
    for i, (c, s, e) in enumerate(zip(L.codes.tolist(), L.starts.tolist(), L.ends.tolist())):
        ix.add(names[c], s, e, i)
    ix.index()
    return ix


def loci_to_anndata(L, X=None, *, obs=None, layers=None, **kw):
    """AnnData with the loci as ``var`` (index ``chrom:start-end``).

    ``X`` is (n_obs x n_loci) — e.g. samples x regions; ``obs`` a DataFrame or a
    list of observation names."""
    import anndata as ad
    import pandas as pd
    var = loci_to_pandas(L)
    var.index = pd.Index(L.names, dtype=str)
    if X is not None:
        X = np.asarray(X) if not hasattr(X, "tocsr") else X
        if X.shape[1] != len(L):
            raise ValueError(f"X has {X.shape[1]} columns for {len(L)} loci (X is obs x loci)")
    if obs is not None and not isinstance(obs, pd.DataFrame):
        obs = pd.DataFrame(index=pd.Index([str(o) for o in obs], dtype=str))
    if X is None and obs is None:
        X = np.zeros((0, len(L)), np.float32)
    return ad.AnnData(X=X, obs=obs, var=var, layers=layers, **kw)


# ── signal cubes ───────────────────────────────────────────────────────────

def cube_to_xarray(S, L, tracks=None, *, flank: Optional[int] = None, name: str = "signal"):
    """A (rows x tracks x bins) cube as an ``xarray.DataArray`` with the loci as
    coordinates (region names plus chrom / start / end) and, when ``flank`` is
    given, the bin centres in bp relative to the locus centre."""
    import xarray as xr
    S = np.asarray(S)
    n, t, b = S.shape
    tracks = [f"track_{i}" for i in range(t)] if tracks is None else [str(x) for x in tracks]
    coords = {"region": L.names, "chrom": ("region", L.chroms.astype(str)),
              "start": ("region", L.starts), "end": ("region", L.ends), "track": tracks}
    if flank is not None:
        edges = np.linspace(-flank, flank, b + 1)
        coords["bin"] = (edges[:-1] + edges[1:]) / 2
    else:
        coords["bin"] = np.arange(b)
    return xr.DataArray(S, dims=("region", "track", "bin"), coords=coords, name=name)


def cube_to_anndata(S, L, tracks=None, *, agg: str = "mean"):
    """A (rows x tracks x bins) cube as AnnData: tracks are ``obs``, loci are ``var``,
    ``X`` is the per-locus ``agg`` over bins; with more than one bin, each track's
    profile is kept in ``varm['bins:<track>']`` (loci x bins)."""
    S = np.asarray(S)
    n, t, b = S.shape
    tracks = [f"track_{i}" for i in range(t)] if tracks is None else [str(x) for x in tracks]
    red = {"mean": np.nanmean, "sum": np.nansum, "max": np.nanmax}[agg]
    a = loci_to_anndata(L, red(S, axis=2).T, obs=tracks)
    if b > 1:
        for i, tr in enumerate(tracks):
            a.varm[f"bins:{tr}"] = S[:, i, :]
    return a


def cube_to_pandas(S, L, tracks=None):
    """A cube as a DataFrame indexed by region name: one column per track for a
    single bin, else (track, bin) columns."""
    import pandas as pd
    S = np.asarray(S)
    n, t, b = S.shape
    tracks = [f"track_{i}" for i in range(t)] if tracks is None else [str(x) for x in tracks]
    if b == 1:
        return pd.DataFrame(S[:, :, 0], index=pd.Index(L.names, name="region"), columns=tracks)
    cols = pd.MultiIndex.from_product([tracks, range(b)], names=["track", "bin"])
    return pd.DataFrame(S.reshape(n, t * b), index=pd.Index(L.names, name="region"), columns=cols)


# ── sequences ──────────────────────────────────────────────────────────────

_RC = str.maketrans("ACGTNacgtnRYKMSWrykmsw", "TGCANtgcanYRMKSWyrmksw")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def windows(L, r: Optional[int] = None):
    """(start, end) of each row, or of ``center ± r``."""
    if r is None:
        return L.starts, L.ends
    c = L.centers
    return c - r, c + r


def sequences(L, fasta, r: Optional[int] = None, *, strand: bool = False, upper: bool = False,
              backend: Optional[str] = None) -> list:
    from .backends.fasta import open_fasta
    src = open_fasta(fasta, backend=backend)
    a, b = windows(L, r)
    sizes = src.sizes()
    names = L.genome.names
    lim = np.array([sizes.get(n, -1) for n in names], np.int64)
    size = lim[L.codes] if len(L) else np.zeros(0, np.int64)
    present = size >= 0
    a2 = np.clip(a, 0, None)
    b2 = np.minimum(b, np.where(present, size, 0))
    ok = present & (b2 > a2)
    out = [""] * len(L)
    rows = np.flatnonzero(ok)
    got = src.fetch_many(L.chroms[rows], a2[rows], b2[rows])
    for i, s in zip(rows.tolist(), got):
        out[i] = s
    if upper:
        out = [s.upper() for s in out]
    if strand:
        neg = np.flatnonzero(L.strands == 2)
        for i in neg.tolist():
            out[i] = revcomp(out[i])
    return out


def to_seqrecords(L, fasta, r: Optional[int] = None, *, strand: bool = False, **kw):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    seqs = sequences(L, fasta, r=r, strand=strand, **kw)
    return [SeqRecord(Seq(s), id=u, description="") for s, u in zip(seqs, L.uid)]


def write_fasta(L, path: str, fasta, r: Optional[int] = None, *, strand: bool = False,
                width: int = 60, **kw) -> str:
    seqs = sequences(L, fasta, r=r, strand=strand, **kw)
    with open(path, "w") as f:
        for u, s in zip(L.uid, seqs):
            f.write(f">{u}\n")
            for k in range(0, len(s), width):
                f.write(s[k:k + width] + "\n")
    return path


# ── liftover ───────────────────────────────────────────────────────────────

_LIFTOVER_CACHE: dict = {}


def liftover(L, chain_file: str, *, min_match: float = 0.95, verbose: bool = True):
    """UCSC-style liftover through pyliftover (see ``Loci.liftover``)."""
    try:
        from pyliftover import LiftOver
    except ImportError:
        raise ImportError("Loci.liftover needs pyliftover: pip install pyliftover") from None
    from .loci import SCODE, STRANDS, Loci
    key = (os.path.abspath(chain_file), os.path.getmtime(chain_file), os.path.getsize(chain_file))
    lo = _LIFTOVER_CACHE.get(key)
    if lo is None:
        _LIFTOVER_CACHE.clear()                         # one chain at a time is the common case
        lo = _LIFTOVER_CACHE[key] = LiftOver(chain_file)

    def walk(chrom, pos, q_strand, step, budget):
        for off in range(budget + 1):
            hits = lo.convert_coordinate(chrom, pos + step * off, q_strand)
            if hits:
                return hits[0], off
        return None, None

    ch, st, en, sd, src = [], [], [], [], []
    for i, (c, s, e, d) in enumerate(zip(L.chroms.tolist(), L.starts.tolist(), L.ends.tolist(),
                                         STRANDS[L.strands].tolist())):
        length = e - s
        if length <= 0:
            continue
        q = d if d in ("+", "-") else "+"
        budget = max(0, int((1.0 - min_match) * length))
        s_hit, s_off = walk(c, s, q, +1, budget)
        e_hit, e_off = walk(c, e - 1, q, -1, budget)
        if s_hit is None or e_hit is None or s_hit[0] != e_hit[0] or s_hit[2] != e_hit[2] \
                or s_off + e_off > budget:
            continue
        sp, ep = s_hit[1], e_hit[1]
        a, b = (sp, ep + 1) if sp <= ep else (ep, sp + 1)
        ch.append(s_hit[0])
        st.append(a)
        en.append(b)
        sd.append(SCODE[s_hit[2]] if d in ("+", "-") else 0)
        src.append(i)
    g = L.genome
    out = Loci(g.encode(ch) if ch else np.zeros(0, np.int32), np.asarray(st, np.int64), np.asarray(en, np.int64),
               np.asarray(sd, np.int8), genome=g, cols={"source_row": np.asarray(src, np.int64)})
    if verbose:
        print(f"[INFO] liftover: {len(L)} → {len(out)} ({len(L) - len(out)} dropped, min_match={min_match})")
    return out
