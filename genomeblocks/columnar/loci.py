"""Columnar Loci: four numpy columns, objects only when you look at one row.

    codes   int32   chromosome code (into the shared Genome)
    starts  int64
    ends    int64
    strands int8    0 '.', 1 '+', 2 '-'
    cols    dict    optional extra columns (name, score, gene_id, ...)

``L[i]`` and iteration hand out :class:`LocusView` objects — real ``Locus``
subclasses whose fields read from (and write to) the columns — so code that
expects Locus objects keeps working. Whole-set operations (intersect, merge,
nearest, ...) run on the arrays.

Row order is meaningful: by default ``make`` sorts into genome order, so
each chromosome is one contiguous block of rows (``L.chrom_slice('chr2')``).
Every other table that describes these loci (annotations, signal cubes,
Architecture vertices) is aligned to these rows: row ``i`` is the same
locus everywhere.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from ..locus import Locus
from . import _intervals as K
from .genome import Genome, default_genome

STRANDS = np.array([".", "+", "-"], dtype=object)
SCODE = {".": 0, "+": 1, "-": 2}


class LocusView(Locus):
    """A Locus whose fields live in a columnar Loci (read/write-through)."""

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


class Loci:
    """Columnar interval set aligned to a shared Genome."""

    def __init__(self, codes=(), starts=(), ends=(), strands=None, *, genome: Optional[Genome] = None,
                 cols: Optional[Dict[str, np.ndarray]] = None, filename: Optional[str] = None,
                 is_sorted: Optional[bool] = None):
        self.genome = genome if genome is not None else default_genome()
        self.codes = np.asarray(codes, np.int32)
        self.starts = np.asarray(starts, np.int64)
        self.ends = np.asarray(ends, np.int64)
        self.strands = (np.zeros(len(self.starts), np.int8) if strands is None
                        else np.asarray(strands, np.int8))
        self.cols = dict(cols or {})
        self.filename = filename
        self._sorted = is_sorted
        self._dirty()

    def _dirty(self):
        self._uid = None
        self._uids = None
        self._cgr = None
        self._offsets = None
        self._order = None

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, filename: str, *, genome: Optional[Genome] = None, sort: bool = True,
             names: bool = False) -> "Loci":
        """Read a BED / narrowPeak file into columns (genome order by default)."""
        import pandas as pd
        ncol = 3
        with open(filename) as f:
            for line in f:
                if line.strip() and not line.startswith(("#", "track", "browser")):
                    ncol = len(line.rstrip("\n").split("\t"))
                    break
        use = [0, 1, 2] + ([3] if names and ncol > 3 else []) + ([5] if ncol > 5 else [])
        df = pd.read_csv(filename, sep="\t", header=None, comment="#", usecols=use,
                         dtype={0: "category", 5: "category", 3: str})
        g = genome if genome is not None else default_genome()
        cat = df[0]
        lut = np.fromiter((g._add(str(c)) for c in cat.cat.categories), np.int32,
                          len(cat.cat.categories))
        codes = lut[cat.cat.codes.to_numpy()]
        strands = (df[5].astype(str).map(SCODE).fillna(0).to_numpy(np.int8)
                   if 5 in df else None)
        cols = {"name": df[3].to_numpy(object)} if 3 in df else None
        L = cls(codes, df[1].to_numpy(np.int64), df[2].to_numpy(np.int64), strands,
                genome=g, cols=cols, filename=filename)
        return L.sort() if sort else L

    @classmethod
    def from_frame(cls, df, chrom=None, start=None, end=None, strand=None, *,
                   genome: Optional[Genome] = None, keep=()) -> "Loci":
        """From a pandas/polars frame (first three columns unless named)."""
        if hasattr(df, "to_pandas") and not hasattr(df, "iloc"):
            df = df.to_pandas()
        c = df[chrom] if chrom else df.iloc[:, 0]
        s = df[start] if start else df.iloc[:, 1]
        e = df[end] if end else df.iloc[:, 2]
        g = genome if genome is not None else default_genome()
        st = df[strand].map(SCODE).fillna(0).to_numpy(np.int8) if strand else None
        return cls(g.encode(c.to_numpy()), s.to_numpy(np.int64), e.to_numpy(np.int64), st,
                   genome=g, cols={k: df[k].to_numpy() for k in keep})

    @classmethod
    def from_loci(cls, items, *, genome: Optional[Genome] = None) -> "Loci":
        """From any iterable of Locus objects (e.g. a classic ``genomeblocks.Loci``)."""
        items = list(items)
        g = genome if genome is not None else default_genome()
        n = len(items)
        return cls(g.encode([l.chrom for l in items]),
                   np.fromiter((l.start for l in items), np.int64, n),
                   np.fromiter((l.end for l in items), np.int64, n),
                   np.fromiter((SCODE.get(l.strand, 0) for l in items), np.int8, n), genome=g)

    def take(self, idx) -> "Loci":
        """Rows ``idx`` (int array, bool mask or slice) as a new Loci."""
        out = Loci(self.codes[idx], self.starts[idx], self.ends[idx], self.strands[idx],
                   genome=self.genome, cols={k: v[idx] for k, v in self.cols.items()})
        if isinstance(idx, slice) or (np.asarray(idx).dtype == bool):
            out._sorted = self._sorted                    # order-preserving selections
        return out

    def copy(self) -> "Loci":
        return Loci(self.codes.copy(), self.starts.copy(), self.ends.copy(), self.strands.copy(),
                    genome=self.genome, cols={k: v.copy() for k, v in self.cols.items()},
                    filename=self.filename, is_sorted=self._sorted)

    # ── sequence protocol ─────────────────────────────────────────────────
    def __len__(self):
        return len(self.starts)

    def __bool__(self):
        return len(self.starts) > 0

    def __iter__(self):
        for i in range(len(self.starts)):
            yield LocusView(self, i)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            n = len(self.starts)
            if not -n <= key < n:
                raise IndexError(key)
            return LocusView(self, int(key) % n)
        if isinstance(key, str):
            if key in self.cols:
                return self.cols[key]
            try:
                return LocusView(self, self.uids[key])
            except KeyError:
                raise KeyError(f"No Locus with UID '{key}' found.") from None
        return self.take(key)

    def row(self, uid: str) -> int:
        return self.uids[uid]

    # ── derived columns (cheap, cached where it matters) ─────────────────
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
    def uid(self) -> np.ndarray:
        """``chrom:start-end(strand)`` per row, built once on first use."""
        if self._uid is None:
            if len(self) == 0:
                self._uid = np.zeros(0, object)
            else:
                import pandas as pd
                c = pd.Series(self.chroms)
                self._uid = (c + ":" + pd.Series(self.starts).astype(str) + "-"
                             + pd.Series(self.ends).astype(str) + "("
                             + pd.Series(STRANDS[self.strands]) + ")").to_numpy(object)
        return self._uid

    @property
    def uids(self) -> dict:
        """uid -> row (classic-API compatibility; built once)."""
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
        """The rows on ``chrom`` (a zero-copy view on sorted Loci)."""
        return self.take(self.chrom_slice(chrom))

    # ── point lookups (cgranges, built lazily, like the classic Loci) ────
    @property
    def cgr(self):
        if self._cgr is None:
            from ..loci import _PyIntervalIndex, _get_cgranges
            cg = _get_cgranges()
            idx = cg.cgranges() if cg is not None else _PyIntervalIndex()
            names = self.genome.names
            for i, (c, s, e) in enumerate(zip(self.codes.tolist(), self.starts.tolist(),
                                              self.ends.tolist())):
                idx.add(names[c], s, e, i)
            idx.index()
            self._cgr = idx
        return self._cgr

    def overlap_rows(self, key, start=None, end=None) -> np.ndarray:
        """Row numbers overlapping one region (a Locus, 'chr:a-b', or chrom, a, b)."""
        if isinstance(key, Locus):
            key, start, end = key.chrom, key.start, key.end
        elif start is None:
            key, start, end = _parse_region(key)
        return np.array(sorted(j for *_, j in self.cgr.overlap(key, int(start), int(end))), np.int64)

    def overlaps(self, key, start=None, end=None) -> "Loci":
        """Rows overlapping one region: ``L.overlaps('chr1', a, b)`` or a Locus / 'chr:a-b'."""
        return self.take(self.overlap_rows(key, start, end))

    # ── whole-set operations (vectorised) ────────────────────────────────
    def _c(self):
        return self.codes, self.starts, self.ends

    def _check(self, o):
        if not isinstance(o, Loci):
            o = Loci.from_loci(o, genome=self.genome)
        if o.genome is not self.genome:
            lut = np.array([self.genome._add(n) for n in o.genome.names], np.int32)
            o = Loci(lut[o.codes], o.starts, o.ends, o.strands, genome=self.genome)
        return o

    def overlap_any(self, o) -> np.ndarray:
        """Boolean mask: row overlaps anything in ``o``."""
        o = self._check(o)
        return K.overlaps_any(*self._c(), *o._c())

    def overlap_pairs(self, o):
        """(rows of self, rows of o) for every overlapping pair."""
        o = self._check(o)
        return K.overlap_pairs(*self._c(), *o._c())

    def intersect(self, o) -> "Loci":
        return self.take(self.overlap_any(o))

    def difference(self, o) -> "Loci":
        return self.take(~self.overlap_any(o))

    __and__ = intersect
    __sub__ = difference
    __truediv__ = difference

    def __add__(self, o) -> "Loci":
        o = self._check(o)
        keys = set(self.cols) & set(o.cols)
        return Loci(np.concatenate([self.codes, o.codes]), np.concatenate([self.starts, o.starts]),
                    np.concatenate([self.ends, o.ends]), np.concatenate([self.strands, o.strands]),
                    genome=self.genome,
                    cols={k: np.concatenate([self.cols[k], o.cols[k]]) for k in keys})
    __or__ = __add__

    def slop(self, n: int) -> "Loci":
        return Loci(self.codes, np.maximum(self.starts - n, 0), self.ends + n, self.strands,
                    genome=self.genome, cols=self.cols)

    def sort(self) -> "Loci":
        """Genome order (natural chromosome order, then start, then end)."""
        r = self.genome.rank[self.codes] if len(self) else np.zeros(0, np.int64)
        out = self.take(np.lexsort((self.ends, self.starts, r)))
        out._sorted, out.filename = True, self.filename
        return out

    def merge(self) -> "Loci":
        """bedtools-merge: overlapping or book-ended rows fuse (strand of the first)."""
        c, s, e, first = K.merge(*self._c())
        return Loci(c, s, e, self.strands[first], genome=self.genome).sort()

    def nearest(self, o):
        """(row in o, distance) of the nearest interval of ``o`` for every row.

        Overlaps have distance 0 (the overlapping interval with the lowest
        start wins); otherwise the closer of the left / right neighbour.
        Rows with nothing on their chromosome get (-1, -1).
        """
        o = self._check(o)
        og = K.gpos(o.codes, o.starts)
        order = np.lexsort((o.starts, o.codes.astype(np.int64)))
        og, oe = og[order], K.gpos(o.codes, o.ends)[order]
        qs, qe = K.gpos(self.codes, self.starts), K.gpos(self.codes, self.ends)
        n = len(self)
        best = np.full(n, -1, np.int64)
        dist = np.full(n, np.iinfo(np.int64).max, np.int64)
        # left: the interval with the largest end among those starting before q.end
        k = np.searchsorted(og, qe, side="left")
        run_max = np.maximum.accumulate(oe) if len(oe) else oe
        arg_max = _running_argmax(oe)
        has = k > 0
        j = np.where(has, arg_max[np.maximum(k - 1, 0)], -1)
        same = has & (o.codes[order][np.maximum(j, 0)] == self.codes)
        d_left = np.where(same, np.maximum(qs - run_max[np.maximum(k - 1, 0)], 0), dist)
        # overlap: first interval (by start) overlapping q — scan forward from the
        # first start > q.start - maxlen is exact but slower; for an overlap any
        # overlapping interval has distance 0, pick the lowest start among them.
        pairs_q, pairs_r = K.overlap_pairs(self.codes, self.starts, self.ends,
                                           o.codes, o.starts, o.ends)
        ov = np.zeros(n, bool)
        if len(pairs_q):
            ostart = o.starts[pairs_r]
            srt = np.lexsort((pairs_r, ostart, pairs_q))
            pq, pr = pairs_q[srt], pairs_r[srt]
            firstq = np.r_[True, pq[1:] != pq[:-1]]
            best[pq[firstq]] = pr[firstq]
            dist[pq[firstq]] = 0
            ov[pq[firstq]] = True
        # right: first interval starting at or after q.end
        k2 = np.searchsorted(og, qe, side="left")
        ok2 = k2 < len(og)
        j2 = np.where(ok2, np.minimum(k2, len(og) - 1), 0)
        same2 = ok2 & (o.codes[order][j2] == self.codes)
        d_right = np.where(same2, og[j2] - qe, dist) if len(og) else dist
        use_left = ~ov & same & (d_left <= d_right)
        use_right = ~ov & same2 & ~use_left
        best[use_left] = order[j[use_left]]
        dist[use_left] = d_left[use_left]
        best[use_right] = order[j2[use_right]]
        dist[use_right] = d_right[use_right]
        dist[best < 0] = -1
        return best, dist

    # ── export ────────────────────────────────────────────────────────────
    def to_pandas(self, uid: bool = False):
        import pandas as pd
        d = {"chrom": pd.Categorical.from_codes(self.codes, self.genome.names),
             "start": self.starts, "end": self.ends,
             "strand": pd.Categorical.from_codes(self.strands, [".", "+", "-"])}
        d.update(self.cols)
        if uid:
            d["uid"] = self.uid
        return pd.DataFrame(d)

    def to_frame(self, names=None):
        """Classic ``Loci.to_frame`` layout (Chr/Start/End/Strand/Name)."""
        import pandas as pd
        return pd.DataFrame({"Chr": self.chroms, "Start": self.starts, "End": self.ends,
                             "Strand": STRANDS[self.strands],
                             "Name": names if names is not None else self.uid})

    def to_polars(self):
        import polars as pl
        return pl.from_arrow(self.to_arrow())

    def to_arrow(self):
        import pyarrow as pa
        chrom = pa.DictionaryArray.from_arrays(pa.array(self.codes), pa.array(self.genome.names))
        strand = pa.DictionaryArray.from_arrays(pa.array(self.strands), pa.array([".", "+", "-"]))
        arrs = {"chrom": chrom, "start": pa.array(self.starts), "end": pa.array(self.ends),
                "strand": strand}
        arrs.update({k: pa.array(v) for k, v in self.cols.items()})
        return pa.table(arrs)

    def to_legacy(self):
        """A classic ``genomeblocks.Loci`` of plain Locus objects."""
        from ..loci import Loci as Classic
        names = self.genome.names
        return Classic(Locus(names[c], s, e, STRANDS[d]) for c, s, e, d in
                       zip(self.codes.tolist(), self.starts.tolist(), self.ends.tolist(),
                           self.strands.tolist()))

    def to_bed(self, path=None):
        lines = "\n".join(f"{c}\t{s}\t{e}\t{u}\t0\t{d}" for c, s, e, u, d in
                          zip(self.chroms, self.starts, self.ends, self.uid, STRANDS[self.strands]))
        if path is None:
            return lines
        with open(path, "w") as f:
            f.write(lines + "\n")

    # ── classic methods (signal, motifs, atlas, ...) keep working ────────
    def __getattr__(self, name):
        """Anything the classic Loci has (``signal``, ``scan_motifs``, ...) runs
        unchanged: those functions only iterate rows / use ``cgr`` / ``uids``.
        The classic module (and its plotting deps) is imported on first use."""
        if name.startswith("_"):
            raise AttributeError(name)
        from ..loci import Loci as Classic
        fn = getattr(Classic, name, None)
        if fn is None or not callable(fn):
            raise AttributeError(f"'Loci' has no attribute {name!r}")
        return fn.__get__(self)

    # ── display ───────────────────────────────────────────────────────────
    def __repr__(self):
        nch = len(np.unique(self.codes)) if len(self) else 0
        return f"Loci(n={len(self):,}, chroms={nch}, columnar{', sorted' if self._sorted else ''})"

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

    # ── persistence ───────────────────────────────────────────────────────
    def save(self, path: str):
        import pyarrow.parquet as pq
        pq.write_table(self.to_arrow(), path)

    @classmethod
    def load(cls, path: str, *, genome: Optional[Genome] = None) -> "Loci":
        import pyarrow.parquet as pq
        t = pq.read_table(path)
        g = genome if genome is not None else default_genome()
        ch = t.column("chrom").combine_chunks()
        lut = np.array([g._add(str(n)) for n in ch.dictionary.to_pylist()], np.int32)
        st = t.column("strand").combine_chunks()
        cols = {k: t.column(k).to_numpy() for k in t.column_names
                if k not in ("chrom", "start", "end", "strand")}
        L = cls(lut[ch.indices.to_numpy()], t.column("start").to_numpy(),
                t.column("end").to_numpy(), st.indices.to_numpy().astype(np.int8),
                genome=g, cols=cols, filename=path)
        return L


def _running_argmax(a: np.ndarray) -> np.ndarray:
    """idx[i] = argmax(a[:i+1]) (first index of the running maximum)."""
    if len(a) == 0:
        return np.zeros(0, np.int64)
    run = np.maximum.accumulate(a)
    is_new = np.r_[True, a[1:] > run[:-1]]
    idx = np.where(is_new, np.arange(len(a)), 0)
    return np.maximum.accumulate(idx)


def _parse_region(s: str):
    chrom, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return chrom, int(a), int(b)
