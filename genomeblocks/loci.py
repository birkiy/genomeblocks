"""Loci: genomic intervals as a table of numpy columns.

    codes    int32   chromosome code into the shared Genome
    starts   int64   0-based start (BED)
    ends     int64   end, exclusive
    strands  int8    0 '.', 1 '+', 2 '-'
    cols     dict    any other columns (name, score, gene_name, ...), aligned to the rows

The row number is the join key: a signal cube, a motif matrix, an annotation
vector or an Architecture vertex column computed from these loci has one row
per locus, in this order. ``L[i]`` hands out a :class:`LocusView` (a real
:class:`~genomeblocks.Locus` reading the columns); whole-set work (overlap,
merge, nearest, ...) runs on the arrays through the interval backend.

Every Loci converts to and from the rest of the ecosystem — pandas, polars,
Arrow / parquet, bioframe, pyranges, pybedtools, cgranges, AnnData, Biopython —
see :mod:`genomeblocks.interop`.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from . import _intervals as K
from ._table import TableMixin
from .genome import Genome, read_sizes
from .locus import Locus, parse_region

STRANDS = np.array([".", "+", "-"], dtype=object)
SCODE = {".": 0, "+": 1, "-": 2}


class LocusView(Locus):
    """A Locus whose fields live in a Loci (read / write-through)."""

    def __init__(self, owner: "Loci", i: int):          # nothing is copied
        self.__dict__["_o"] = owner
        self.__dict__["_i"] = i

    def _set(self, name, value):
        o, i = self.__dict__["_o"], self.__dict__["_i"]
        if name == "chrom":
            o.codes[i] = o.genome._add(value)
        elif name == "strand":
            o.strands[i] = SCODE[value]
        else:
            getattr(o, name + "s")[i] = value
        o._dirty()

    chrom = property(lambda s: s._o.genome.names[s._o.codes[s._i]], lambda s, v: s._set("chrom", v))
    start = property(lambda s: int(s._o.starts[s._i]), lambda s, v: s._set("start", v))
    end = property(lambda s: int(s._o.ends[s._i]), lambda s, v: s._set("end", v))
    strand = property(lambda s: STRANDS[s._o.strands[s._i]], lambda s, v: s._set("strand", v))

    @property
    def row(self) -> int:
        return self.__dict__["_i"]

    def __getattr__(self, name):                       # extra columns: L[i].score
        cols = self.__dict__["_o"].cols
        if name in cols:
            v = cols[name][self.__dict__["_i"]]
            return v.item() if hasattr(v, "item") else v
        raise AttributeError(name)

    def __setattr__(self, name, value):
        if name in ("chrom", "start", "end", "strand"):
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value

    def __repr__(self):
        return f"Locus[{self.row}]({self.uid})"


class Loci(TableMixin):
    """Genomic intervals as columns sharing one :class:`~genomeblocks.Genome`.

    Build one with :meth:`make` (BED-like files), :meth:`from_frame` (pandas,
    polars, pyarrow, bioframe, PyRanges, dicts of columns) or
    :func:`genomeblocks.as_loci` (anything). It behaves like a table:
    ``len``, ``shape``, ``columns``, ``head()``, ``describe()``, ``L['start']``,
    and it hands itself to pandas / polars / Arrow / seaborn / plotly / altair
    through the standard protocols."""

    def __init__(self, codes=(), starts=(), ends=(), strands=None, *, genome: Optional[Genome] = None,
                 cols: Optional[Dict[str, np.ndarray]] = None, filename: Optional[str] = None,
                 is_sorted: Optional[bool] = None):
        self.genome = genome if genome is not None else Genome()
        self.codes = np.asarray(codes, np.int32)
        self.starts = np.asarray(starts, np.int64)
        self.ends = np.asarray(ends, np.int64)
        self.strands = (np.zeros(len(self.starts), np.int8) if strands is None
                        else np.asarray(strands, np.int8))
        self.cols = {k: np.asarray(v) for k, v in (cols or {}).items()}
        self.filename = filename
        self._sorted = is_sorted
        n = len(self.starts)
        if not (len(self.codes) == len(self.ends) == len(self.strands) == n):
            raise ValueError("codes, starts, ends and strands must have the same length")
        for k, v in self.cols.items():
            if len(v) != n:
                raise ValueError(f"column {k!r} has {len(v)} values for {n} rows")
        self._dirty()

    def _dirty(self):
        self._uid = None
        self._uids = None
        self._offsets = None
        self._indexes = {}

    def _index_cache(self, name, build):
        """A lookup index (built once per backend, dropped when the rows change)."""
        ix = self._indexes.get(name)
        if ix is None:
            ix = self._indexes[name] = build()
        return ix

    def __getstate__(self):
        d = self.__dict__.copy()
        for k in ("_uid", "_uids", "_offsets"):
            d[k] = None
        d["_indexes"] = {}
        return d

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, filename: str, *, genome: Optional[Genome] = None, sort: bool = True, keep=False,
             backend: Optional[str] = None) -> "Loci":
        """Read a BED / narrowPeak / broadPeak file (``.gz`` too).

        Coordinates and strand are always read. ``keep=True`` also keeps the
        other standard columns under their usual names (name, score,
        signalValue, pValue, qValue, peak for narrowPeak); a list keeps just
        those. Rows are sorted into genome order unless ``sort=False``.
        ``backend`` picks the table parser (polars or pandas).
        """
        from .backends.tables import read_columns, sniff
        low = str(filename).lower().removesuffix(".gz")
        if low.endswith((".gtf", ".gff", ".gff3")):
            raise ValueError("this is a gene annotation: use Genes.make()")
        if low.endswith(".bedpe"):
            raise ValueError("this is a BEDPE file: use Pairs.make()")
        _, ncol = sniff(filename)
        names = _bed_names(low, ncol)
        wanted = []
        if keep is True:
            wanted = [i for i in range(3, ncol) if i != 5]
        elif keep:
            wanted = [names.index(k) for k in keep if k in names and names.index(k) >= 3]
        use = [0, 1, 2] + ([5] if ncol > 5 else []) + wanted
        ints = [1, 2] + [i for i in wanted if names[i] in ("peak", "thickStart", "thickEnd", "blockCount")]
        floats = [i for i in wanted if names[i] in ("score", "signalValue", "pValue", "qValue")]
        d = read_columns(filename, use, ints=ints, floats=floats, backend=backend)
        g = genome if genome is not None else Genome()
        strands = (np.array([SCODE.get(x, 0) for x in d[5]], np.int8) if 5 in d else None)
        L = cls(g.encode(d[0]), d[1], d[2], strands, genome=g,
                cols={names[i]: d[i] for i in wanted}, filename=str(filename))
        return L.sort() if sort else L

    @classmethod
    def from_frame(cls, df, chrom=None, start=None, end=None, strand=None, *, keep=True,
                   genome: Optional[Genome] = None) -> "Loci":
        """From any table: pandas, polars (eager or lazy), pyarrow, or bioframe /
        pyranges frames. Columns are found by name (chrom / chromosome / Chromosome
        / seqnames / chr, start / Start / chromStart, end / End / chromEnd,
        strand / Strand) or else taken as the first three. ``keep=True`` keeps the
        other columns (a list keeps those). Row order is kept."""
        from . import interop
        return interop.loci_from_frame(df, chrom, start, end, strand, keep=keep, genome=genome)

    from_pandas = from_frame

    @classmethod
    def from_polars(cls, df, **kw) -> "Loci":
        return cls.from_frame(df, **kw)

    @classmethod
    def from_arrow(cls, table, **kw) -> "Loci":
        return cls.from_frame(table, **kw)

    @classmethod
    def from_bioframe(cls, df, **kw) -> "Loci":
        return cls.from_frame(df, **kw)

    @classmethod
    def from_pyranges(cls, gr, **kw) -> "Loci":
        from . import interop
        return interop.loci_from_pyranges(gr, **kw)

    @classmethod
    def from_bedtool(cls, bt, **kw) -> "Loci":
        from . import interop
        return interop.loci_from_bedtool(bt, **kw)

    @classmethod
    def from_anndata(cls, adata, axis: str = "var", **kw) -> "Loci":
        """Regions of an AnnData (scATAC peaks as ``var`` by default), in its order:
        chrom/start/end columns if present, else parsed from the names
        (``chr1:100-200``, ``chr1-100-200``, ``chr1_100_200``)."""
        from . import interop
        return interop.loci_from_anndata(adata, axis=axis, **kw)

    @classmethod
    def from_records(cls, items, *, genome: Optional[Genome] = None) -> "Loci":
        """From Locus objects, ``(chrom, start, end[, strand])`` tuples or region
        strings, in order."""
        g = genome if genome is not None else Genome()
        ch, st, en, sd = [], [], [], []
        for x in items:
            if isinstance(x, Locus):
                c, s, e, d = x.chrom, x.start, x.end, x.strand
            elif isinstance(x, str):
                if x.rstrip().endswith(")"):                       # a uid: chrom:start-end(strand)
                    l = Locus.from_uid(x.strip())
                    c, s, e, d = l.chrom, l.start, l.end, l.strand
                else:
                    (c, s, e), d = parse_region(x), "."
            else:
                c, s, e = x[0], x[1], x[2]
                d = x[3] if len(x) > 3 else "."
            ch.append(c)
            st.append(s)
            en.append(e)
            sd.append(SCODE.get(d, 0))
        return cls(g.encode(ch) if ch else np.zeros(0, np.int32), np.asarray(st, np.int64),
                   np.asarray(en, np.int64), np.asarray(sd, np.int8), genome=g)

    @classmethod
    def from_uids(cls, uids, *, genome: Optional[Genome] = None) -> "Loci":
        """From ``chrom:start-end(strand)`` strings (what :attr:`uid` makes)."""
        return cls.from_records([Locus.from_uid(u) for u in uids], genome=genome)

    @classmethod
    def tile(cls, chrom: str, size: int, chromsizes, *, genome: Optional[Genome] = None) -> "Loci":
        """Uniform ``size``-bp tiles across one chromosome (the last one clipped).
        ``chromsizes`` is a length, a dict, a .chrom.sizes path, a Series, a
        Genome or a cooler."""
        n = int(chromsizes) if isinstance(chromsizes, (int, np.integer)) else read_sizes(chromsizes)[chrom]
        return cls.tile_genome({chrom: n}, size, genome=genome)

    @classmethod
    def tile_genome(cls, chromsizes, size: int, chroms=None, *, genome: Optional[Genome] = None) -> "Loci":
        """Uniform ``size``-bp tiles over every chromosome of ``chromsizes``
        (or just ``chroms``, in that order)."""
        sizes = read_sizes(chromsizes)
        g = genome if genome is not None else Genome()
        chroms = list(sizes) if chroms is None else list(chroms)
        parts_c, parts_s, parts_e = [], [], []
        for c in chroms:
            n = sizes[c]
            s = np.arange(0, n, size, dtype=np.int64)
            parts_s.append(s)
            parts_e.append(np.minimum(s + size, n))
            parts_c.append(np.full(len(s), g._add(c), np.int32))
        if not parts_s:
            return cls(genome=g)
        return cls(np.concatenate(parts_c), np.concatenate(parts_s), np.concatenate(parts_e), genome=g)

    def take(self, idx) -> "Loci":
        """Rows ``idx`` (int array, bool mask or slice) as a new Loci."""
        if isinstance(idx, list):
            idx = np.asarray(idx, dtype=bool if idx and isinstance(idx[0], (bool, np.bool_)) else np.int64)
        out = Loci(self.codes[idx], self.starts[idx], self.ends[idx], self.strands[idx],
                   genome=self.genome, cols={k: v[idx] for k, v in self.cols.items()})
        if isinstance(idx, slice) or (np.asarray(idx).dtype == bool):
            out._sorted = self._sorted                    # order-preserving selections
        return out

    def copy(self) -> "Loci":
        return Loci(self.codes.copy(), self.starts.copy(), self.ends.copy(), self.strands.copy(),
                    genome=self.genome, cols={k: v.copy() for k, v in self.cols.items()},
                    filename=self.filename, is_sorted=self._sorted)

    def equals(self, other: "Loci", cols: bool = False) -> bool:
        """Same rows (chromosome, start, end, strand) in the same order."""
        if len(self) != len(other):
            return False
        o = self._check(other)
        same = (np.array_equal(self.codes, o.codes) and np.array_equal(self.starts, o.starts)
                and np.array_equal(self.ends, o.ends) and np.array_equal(self.strands, o.strands))
        if same and cols:
            same = set(self.cols) == set(other.cols) and all(
                np.array_equal(np.asarray(self.cols[k]), np.asarray(other.cols[k])) for k in self.cols)
        return same

    # ── sequence protocol ─────────────────────────────────────────────────
    def __len__(self):
        return len(self.starts)

    def __bool__(self):
        return len(self.starts) > 0

    def __iter__(self):
        for i in range(len(self.starts)):
            yield LocusView(self, i)

    def __getitem__(self, key):
        """``L[i]`` a row (Locus), ``L['start']`` / ``L['score']`` a column (numpy),
        ``L['chr1:5-9(.)']`` the row with that uid, ``L[mask]`` / ``L[rows]`` /
        ``L[a:b]`` a new Loci."""
        if isinstance(key, (int, np.integer)):
            n = len(self.starts)
            if not -n <= key < n:
                raise IndexError(key)
            return LocusView(self, int(key) % n)
        if isinstance(key, str):
            core = _CORE.get(key)
            if core is not None:
                return core(self)
            if key in self.cols:
                return self.cols[key]
            try:
                return LocusView(self, self.uids[key])
            except KeyError:
                raise KeyError(f"No column or locus {key!r}.") from None
        if isinstance(key, (list, tuple, np.ndarray)) and len(key) and isinstance(
                np.asarray(key).ravel()[0], str):
            return self.take(np.array([self.uids[k] for k in key], np.int64))
        return self.take(key)

    def __setitem__(self, key: str, values):
        """``L['score'] = values`` adds or replaces a column (one value per row)."""
        if not isinstance(key, str):
            raise TypeError("only columns can be set: L['name'] = values")
        if key in ("chrom", "start", "end", "strand"):
            raise KeyError(f"{key!r} is a coordinate; build a new Loci instead")
        v = np.asarray(values)
        if v.ndim == 0:
            v = np.full(len(self), v.item())
        if len(v) != len(self):
            raise ValueError(f"{len(v)} values for {len(self)} rows")
        self.cols[key] = v

    def __contains__(self, key) -> bool:
        if isinstance(key, Locus):
            key = key.uid
        return key in self.uids

    def row(self, uid: str) -> int:
        return self.uids[uid]

    # ── like a DataFrame ─────────────────────────────────────────────────
    @property
    def columns(self) -> list:
        return ["chrom", "start", "end", "strand"] + list(self.cols)

    def head(self, n: int = 5) -> "Loci":
        return self.take(slice(0, n))

    def tail(self, n: int = 5) -> "Loci":
        return self.take(slice(max(len(self) - n, 0), len(self)))

    def describe(self):
        """A one-table summary: rows, chromosomes, bases covered, length
        distribution, strands."""
        import pandas as pd
        n = len(self)
        lens = self.lengths
        covered = int(self.merge().lengths.sum()) if n else 0
        q = np.percentile(lens, [0, 25, 50, 75, 100]) if n else [0] * 5
        rows = [("rows", n), ("chromosomes", int(len(np.unique(self.codes))) if n else 0),
                ("bases covered (merged)", covered), ("length min", int(q[0])), ("length 25%", float(q[1])),
                ("length median", float(q[2])), ("length 75%", float(q[3])), ("length max", int(q[4])),
                ("length mean", float(lens.mean()) if n else 0.0),
                ("strand + / - / .", "{} / {} / {}".format(*np.bincount(self.strands, minlength=3)[[1, 2, 0]])),
                ("genome-sorted", self.is_sorted)]
        return pd.DataFrame(rows, columns=["", "value"]).set_index("")

    summary = describe

    def to_numpy(self) -> np.ndarray:
        """A structured array (chrom, start, end, strand, then the extra columns)."""
        fields = [("chrom", object), ("start", np.int64), ("end", np.int64), ("strand", object)]
        fields += [(k, np.asarray(v).dtype) for k, v in self.cols.items()]
        out = np.empty(len(self), dtype=fields)
        out["chrom"], out["start"], out["end"], out["strand"] = self.chroms, self.starts, self.ends, self.strand
        for k, v in self.cols.items():
            out[k] = v
        return out

    def __array__(self, dtype=None, copy=None):
        a = self.to_numpy()
        return a if dtype is None else a.astype(dtype)

    # ── derived columns ──────────────────────────────────────────────────
    @property
    def chroms(self) -> np.ndarray:
        return self.genome.decode(self.codes)

    @property
    def centers(self) -> np.ndarray:
        return (self.starts + self.ends) // 2

    @property
    def lengths(self) -> np.ndarray:
        return self.ends - self.starts

    @property
    def strand(self) -> np.ndarray:
        """Strand symbols ('.', '+', '-') per row."""
        return STRANDS[self.strands]

    @property
    def uid(self) -> np.ndarray:
        """``chrom:start-end(strand)`` per row (built once)."""
        if self._uid is None:
            self._uid = _join_uid(self.chroms, self.starts, self.ends, STRANDS[self.strands])
        return self._uid

    @property
    def names(self) -> np.ndarray:
        """``chrom:start-end`` per row (the usual region name in AnnData / scATAC)."""
        return _join_uid(self.chroms, self.starts, self.ends, None)

    @property
    def uids(self) -> dict:
        """uid -> row (built once)."""
        if self._uids is None:
            self._uids = dict(zip(self.uid.tolist(), range(len(self))))
        return self._uids

    # ── chromosome blocks ─────────────────────────────────────────────────
    @property
    def is_sorted(self) -> bool:
        if self._sorted is None:
            k = self.genome.rank[self.codes] * K.STRIDE + self.starts
            self._sorted = bool(np.all(k[1:] >= k[:-1])) if len(k) > 1 else True
        return self._sorted

    @property
    def chrom_offsets(self) -> Dict[str, tuple]:
        """{chrom: (first_row, last_row + 1)} — needs genome-sorted rows."""
        if self._offsets is None:
            if not self.is_sorted:
                raise ValueError("chrom_offsets needs sorted rows: use L.sort() first.")
            cut = np.flatnonzero(np.diff(self.codes)) + 1
            lo = np.concatenate([[0], cut]).astype(np.int64)
            hi = np.concatenate([cut, [len(self)]]).astype(np.int64)
            names = self.genome.names
            self._offsets = {names[self.codes[a]]: (int(a), int(b)) for a, b in zip(lo, hi) if b > a}
        return self._offsets

    def chrom_slice(self, chrom: str) -> slice:
        a, b = self.chrom_offsets.get(chrom, (0, 0))
        return slice(a, b)

    def by_chrom(self, chrom: str) -> "Loci":
        """The rows on ``chrom`` (a view on sorted Loci, a selection otherwise)."""
        if self.is_sorted:
            return self.take(self.chrom_slice(chrom))
        return self.take(self.codes == self.genome.code.get(chrom, -1))

    # ── point lookups ─────────────────────────────────────────────────────
    def overlap_rows(self, key, start=None, end=None, *, backend: Optional[str] = None) -> np.ndarray:
        """Row numbers overlapping one region: ``(chrom, start, end)``, a Locus or
        ``'chr1:1,000-2,000'``. The lookup index is built once and cached."""
        from .backends.intervals import point_rows
        if start is None:
            key, start, end = parse_region(key)
        elif isinstance(key, Locus):
            key, start, end = key.chrom, key.start, key.end
        code = self.genome.code.get(key)
        if code is None:
            return np.zeros(0, np.int64)
        return point_rows(self, code, int(start), int(end), backend=backend)

    def overlaps(self, key, start=None, end=None, *, backend: Optional[str] = None) -> "Loci":
        """Rows overlapping one region, as a Loci."""
        return self.take(self.overlap_rows(key, start, end, backend=backend))

    # ── whole-set operations ─────────────────────────────────────────────
    def _check(self, o) -> "Loci":
        """Bring another interval set onto this Genome (any input :func:`as_loci` takes)."""
        if not isinstance(o, Loci):
            from .interop import as_loci
            o = as_loci(o, genome=self.genome)
        if o.genome is not self.genome:
            lut = np.array([self.genome._add(n) for n in o.genome.names], np.int32)
            o = Loci(lut[o.codes] if len(o) else o.codes, o.starts, o.ends, o.strands, genome=self.genome,
                     cols=o.cols)
        return o

    def overlap_any(self, o, *, backend: Optional[str] = None) -> np.ndarray:
        """Boolean mask: row overlaps anything in ``o``."""
        from .backends.intervals import overlap_any
        return overlap_any(self, self._check(o), backend=backend)

    def overlap_pairs(self, o, *, backend: Optional[str] = None):
        """(rows of self, rows of o) for every overlapping pair, sorted."""
        from .backends.intervals import overlap_pairs
        return overlap_pairs(self, self._check(o), backend=backend)

    def intersect(self, o, *, backend: Optional[str] = None) -> "Loci":
        """Rows that overlap ``o`` (``A & B``, bedtools intersect -u)."""
        return self.take(self.overlap_any(o, backend=backend))

    def difference(self, o, *, backend: Optional[str] = None) -> "Loci":
        """Rows that do not overlap ``o`` (``A - B``, bedtools intersect -v)."""
        return self.take(~self.overlap_any(o, backend=backend))

    __and__ = intersect
    __sub__ = difference
    __truediv__ = difference

    def __add__(self, o) -> "Loci":
        """Concatenation (columns both sides have are kept)."""
        o = self._check(o)
        keys = [k for k in self.cols if k in o.cols]
        return Loci(np.concatenate([self.codes, o.codes]), np.concatenate([self.starts, o.starts]),
                    np.concatenate([self.ends, o.ends]), np.concatenate([self.strands, o.strands]),
                    genome=self.genome,
                    cols={k: np.concatenate([self.cols[k], o.cols[k]]) for k in keys})

    __or__ = __add__
    union = __add__

    def __xor__(self, o) -> "Loci":
        """Rows of either set that do not overlap the other."""
        o = self._check(o)
        return (self - o) + (o - self)

    def slop(self, n: int) -> "Loci":
        """Every interval widened by ``n`` bp on each side (clipped at 0)."""
        return Loci(self.codes, np.maximum(self.starts - n, 0), self.ends + n, self.strands,
                    genome=self.genome, cols=self.cols)

    def sort(self) -> "Loci":
        """Genome order: natural chromosome order (chr1, chr2, ..., chr10), start, end."""
        r = self.genome.rank[self.codes] if len(self) else np.zeros(0, np.int64)
        out = self.take(np.lexsort((self.ends, self.starts, r)))
        out._sorted, out.filename = True, self.filename
        return out

    def merge(self, *, backend: Optional[str] = None) -> "Loci":
        """bedtools merge: overlapping or book-ended rows fuse (strand of the first)."""
        from .backends.intervals import merge
        c, s, e, first = merge(self, backend=backend)
        return Loci(c, s, e, self.strands[first] if len(first) else None, genome=self.genome).sort()

    def nearest(self, o, *, backend: Optional[str] = None):
        """(row in ``o``, distance) of the nearest interval of ``o`` for every row.

        Overlaps are 0 apart; otherwise the gap in bases (book-ended = 0). Rows
        with nothing on their chromosome get (-1, -1)."""
        from .backends.intervals import nearest
        return nearest(self, self._check(o), backend=backend)

    def liftover(self, chain_file: str, *, min_match: float = 0.95, verbose: bool = True) -> "Loci":
        """Lift to another assembly with a UCSC chain file (pyliftover).

        Each end that falls in a chain gap walks inward by up to ``1 - min_match``
        of the length before giving up (UCSC's base-fraction rule). Rows that do
        not lift are dropped; ``cols['source_row']`` says where each came from."""
        from .interop import liftover
        return liftover(self, chain_file, min_match=min_match, verbose=verbose)

    # ── sequences ─────────────────────────────────────────────────────────
    def sequences(self, fasta, r: Optional[int] = None, *, strand: bool = False, upper: bool = False,
                  backend: Optional[str] = None) -> list:
        """Sequence of every row (or of ``center ± r``) from a FASTA path, dict or
        open handle. Windows are clipped to the chromosome (so they can come back
        shorter); ``strand=True`` reverse-complements '-' rows."""
        from .interop import sequences
        return sequences(self, fasta, r=r, strand=strand, upper=upper, backend=backend)

    def to_seqrecords(self, fasta, r: Optional[int] = None, *, strand: bool = False, **kw):
        """Biopython ``SeqRecord`` per row (id = uid)."""
        from .interop import to_seqrecords
        return to_seqrecords(self, fasta, r=r, strand=strand, **kw)

    def to_fasta(self, path: str, fasta, r: Optional[int] = None, *, strand: bool = False, **kw) -> str:
        """Write the rows' sequences as FASTA (headers are uids)."""
        from .interop import write_fasta
        return write_fasta(self, path, fasta, r=r, strand=strand, **kw)

    # ── export ────────────────────────────────────────────────────────────
    def to_pandas(self, uid: bool = False):
        """pandas DataFrame: chrom (categorical, genome order), start, end, strand,
        then the extra columns — bioframe's column names."""
        from .interop import loci_to_pandas
        return loci_to_pandas(self, uid=uid)

    def to_bioframe(self):
        """A bioframe-ready DataFrame (chrom / start / end as plain columns)."""
        from .interop import loci_to_bioframe
        return loci_to_bioframe(self)

    def to_arrow(self):
        from .interop import loci_to_arrow
        return loci_to_arrow(self)

    def to_pyranges(self):
        from .interop import loci_to_pyranges
        return loci_to_pyranges(self)

    def to_bedtool(self):
        from .interop import loci_to_bedtool
        return loci_to_bedtool(self)

    def to_cgranges(self):
        """A built cgranges index (label = row)."""
        from .interop import loci_to_cgranges
        return loci_to_cgranges(self)

    def to_anndata(self, X=None, *, obs=None, layers=None, **kw):
        """AnnData with these loci as ``var`` (names ``chrom:start-end``) and
        ``X`` (n_obs x n_loci) as the data, e.g. a sample x region matrix."""
        from .interop import loci_to_anndata
        return loci_to_anndata(self, X, obs=obs, layers=layers, **kw)

    def to_records(self) -> list:
        """``(chrom, start, end, strand)`` tuples."""
        return list(zip(self.chroms.tolist(), self.starts.tolist(), self.ends.tolist(),
                        STRANDS[self.strands].tolist()))

    def to_bed(self, path: Optional[str] = None, *, name: Optional[str] = None,
               score: Optional[str] = None):
        """BED6 text (written to ``path`` when given): name = ``cols[name]`` or the uid."""
        nm = self.cols[name] if name else self.uid
        sc = self.cols[score] if score else np.zeros(len(self), int)
        lines = "\n".join(f"{c}\t{s}\t{e}\t{u}\t{v}\t{d}" for c, s, e, u, v, d in
                          zip(self.chroms, self.starts.tolist(), self.ends.tolist(), nm, sc,
                              STRANDS[self.strands]))
        if path is None:
            return lines
        with open(path, "w") as f:
            f.write(lines + ("\n" if lines else ""))
        return path

    def save(self, path: str):
        """Parquet (needs pyarrow); :meth:`load` reads it back on any Genome."""
        import pyarrow.parquet as pq
        pq.write_table(self.to_arrow(), path)

    @classmethod
    def load(cls, path: str, *, genome: Optional[Genome] = None) -> "Loci":
        import pyarrow.parquet as pq
        L = cls.from_frame(pq.read_table(path), genome=genome)
        L.filename = path
        return L

    # ── signal, motifs, enrichment, contacts (in their own modules) ─────
    def signal(self, bigwigs, **kw):
        """bigWig signal cube (rows x tracks x bins); see :func:`genomeblocks.signal.signal`."""
        from .signal import signal
        return signal(self, bigwigs, **kw)

    def scan_motifs(self, fasta, motifs, **kw):
        """Total hits per motif; see :func:`genomeblocks.motifs.scan_motifs`."""
        from .motifs import scan_motifs
        return scan_motifs(self, fasta, motifs, **kw)

    def scan_motifs_matrix(self, fasta, motifs, **kw):
        """(rows x motifs) hit counts; see :func:`genomeblocks.motifs.scan_motifs_matrix`."""
        from .motifs import scan_motifs_matrix
        return scan_motifs_matrix(self, fasta, motifs, **kw)

    def scan_motifs_matrix_masked(self, fasta, motifs, anchors, **kw):
        from .motifs import scan_motifs_matrix_masked
        return scan_motifs_matrix_masked(self, fasta, motifs, anchors, **kw)

    def scan_motifs_profile(self, fasta, motifs, select=None, **kw):
        """(rows x motifs x bins) positional hits; see :func:`genomeblocks.motifs.scan_motifs_profile`."""
        from .motifs import scan_motifs_profile
        return scan_motifs_profile(self, fasta, motifs, select, **kw)

    def enrich(self, atlas, **kw):
        """Atlas enrichment of these loci; see :meth:`genomeblocks.Atlas.search`."""
        return atlas.search(self, **kw)

    def enrich_mc(self, atlas, *, n: int = 10, **kw):
        return atlas.bootstrap(self, n=n, **kw)

    def count_pairs(self, pairs_file: str, **kw):
        """Contacts per window by partner chromosome; see :func:`genomeblocks.bedpe.count_pairs`."""
        from .bedpe import count_pairs
        return count_pairs(self, pairs_file, **kw)

    def count_pairs_2d(self, pairs_file: str, **kw):
        from .bedpe import count_pairs_2d
        return count_pairs_2d(self, pairs_file, **kw)

    def call_se(self, bigwigs, **kw):
        """ROSE-style super-enhancers from these peaks; see :func:`genomeblocks.se.call_se`."""
        from .se import call_se
        return call_se(self, bigwigs, **kw)

    def plot_heatmap(self, S, **kw):
        from .signal_draw import plot_heatmap
        return plot_heatmap(self, S, **kw)

    def plot_profiles(self, S, **kw):
        from .signal_draw import plot_profiles
        return plot_profiles(self, S, **kw)

    # ── display ───────────────────────────────────────────────────────────
    def __repr__(self):
        nch = len(np.unique(self.codes)) if len(self) else 0
        extra = f", cols=[{', '.join(self.cols)}]" if self.cols else ""
        return f"Loci(n={len(self):,}, chroms={nch}{', sorted' if self._sorted else ''}{extra})"

    def _repr_html_(self):
        from ._display import table_html
        rows = list(range(min(5, len(self))))
        if len(self) > 10:
            rows += list(range(len(self) - 3, len(self)))
        body = [[str(i), self.genome.names[self.codes[i]], f"{self.starts[i]:,}", f"{self.ends[i]:,}",
                 STRANDS[self.strands[i]]] + [str(v[i]) for v in self.cols.values()] for i in rows]
        head = ["row", "chrom", "start", "end", "strand"] + list(self.cols)
        nch = len(np.unique(self.codes)) if len(self) else 0
        mb = sum(a.nbytes for a in (self.codes, self.starts, self.ends, self.strands)) / 1e6
        return table_html(f"Loci · {len(self):,} rows · {nch} chromosomes · {mb:.1f} MB",
                          head, body, gap_after=4 if len(self) > 10 else None)


_CORE = {"chrom": lambda L: L.chroms, "start": lambda L: L.starts, "end": lambda L: L.ends,
         "strand": lambda L: L.strand}


def _join_uid(chroms, starts, ends, strands):
    if len(starts) == 0:
        return np.zeros(0, object)
    import pandas as pd
    s = (pd.Series(chroms, dtype=object) + ":" + pd.Series(starts).astype(str) + "-"
         + pd.Series(ends).astype(str))
    if strands is not None:
        s = s + "(" + pd.Series(strands, dtype=object) + ")"
    return s.to_numpy(object)


_NARROWPEAK = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
_BED12 = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb",
          "blockCount", "blockSizes", "blockStarts"]


def _bed_names(low: str, ncol: int):
    base = _NARROWPEAK if low.endswith(("narrowpeak", "broadpeak")) else _BED12
    names = list(base[:ncol])
    names += [f"col{i + 1}" for i in range(len(names), ncol)]
    return names
