"""
Pure-Python BigWig reader — thread-safe, zero compiled deps beyond numpy.

Reads the BigWig binary format directly with struct + zlib + numpy + mmap.
All file I/O and decompression release the GIL, enabling true parallelism
via threading.Thread instead of multiprocessing.Process.

Usage::

    from genomeblocks import bigwig
    f = bigwig.open("signal.bw")
    print(f.chroms())                                       # {chrom: size}
    vals = f.stats("chr1", 1000, 2000, n_bins=100)          # binned mean
    bp   = f.values("chr1", 1000, 2000)                     # base-pair resolution
    f.close()
"""
from __future__ import annotations

import builtins as _builtins
import mmap
import struct
import zlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── constants ────────────────────────────────────────────────────────────────
_BW_MAGIC_LE = 0x888FFC26
_BW_MAGIC_BE = 0x26FC8F88
_BP_MAGIC    = 0x78CA8C91
_RT_MAGIC    = 0x2468ACE0

_BED_GRAPH   = 1
_VAR_STEP    = 2
_FIXED_STEP  = 3


# ── reader ───────────────────────────────────────────────────────────────────

class BigWigReader:
    """Thread-safe read-only BigWig file handle.

    Uses mmap for I/O — multiple threads can call stats() / values()
    concurrently since mmap reads and zlib.decompress both release the GIL.
    """

    __slots__ = (
        '_path', '_fh', '_mm', '_end',
        '_version', '_n_zooms', '_chrom_tree_off',
        '_full_data_off', '_full_index_off',
        '_total_summary_off', '_uncompress_buf',
        '_zooms', '_chroms', '_name_to_id', '_id_to_name',
        '_rtree_cache',
    )

    def __init__(self, path: str | Path):
        self._path = str(path)
        self._fh = _builtins.open(self._path, 'rb')
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        self._rtree_cache: Dict[int, object] = {}
        self._parse_header()
        self._parse_chrom_tree()

    # ── low-level helpers ────────────────────────────────────────────────

    def _u(self, fmt: str, off: int) -> tuple:
        """Unpack *fmt* at byte offset *off* using file endianness."""
        return struct.unpack_from(self._end + fmt, self._mm, off)

    # ── header ───────────────────────────────────────────────────────────

    def _parse_header(self):
        raw = struct.unpack_from('<I', self._mm, 0)[0]
        if raw == _BW_MAGIC_LE:
            self._end = '<'
        elif raw == _BW_MAGIC_BE:
            self._end = '>'
        else:
            raise ValueError(f"Not a BigWig file: {self._path} "
                             f"(magic={raw:#010x})")

        self._version, self._n_zooms = self._u('HH', 4)
        (self._chrom_tree_off,)  = self._u('Q', 8)
        (self._full_data_off,)   = self._u('Q', 16)
        (self._full_index_off,)  = self._u('Q', 24)
        (self._total_summary_off,) = self._u('Q', 44)
        (self._uncompress_buf,)  = self._u('I', 52)

        # zoom level headers — 24 bytes each starting at offset 64
        self._zooms: List[Tuple[int, int, int]] = []
        for i in range(self._n_zooms):
            off = 64 + i * 24
            rl, _reserved = self._u('II', off)
            d_off, i_off  = self._u('QQ', off + 8)
            self._zooms.append((rl, d_off, i_off))

    # ── chromosome B+ tree ───────────────────────────────────────────────

    def _parse_chrom_tree(self):
        off = self._chrom_tree_off
        magic = self._u('I', off)[0]
        if magic != _BP_MAGIC:
            raise ValueError(f"Bad B+ tree magic: {magic:#010x}")
        _block_size, key_size, val_size = self._u('III', off + 4)
        # item_count, reserved  — not needed for traversal
        self._chroms: Dict[str, int] = {}
        self._name_to_id: Dict[str, int] = {}
        self._id_to_name: Dict[int, str] = {}
        self._walk_bp(off + 32, key_size, val_size)

    def _walk_bp(self, node: int, ksz: int, vsz: int):
        is_leaf, _, count = self._u('BBH', node)
        off = node + 4
        if is_leaf:
            for _ in range(count):
                raw = bytes(self._mm[off:off + ksz])
                name = raw.split(b'\x00', 1)[0].decode('ascii')
                cid, csz = self._u('II', off + ksz)
                self._chroms[name] = csz
                self._name_to_id[name] = cid
                self._id_to_name[cid] = name
                off += ksz + vsz
        else:
            for _ in range(count):
                child, = self._u('Q', off + ksz)
                self._walk_bp(child, ksz, vsz)
                off += ksz + 8

    # ── R-tree spatial index ─────────────────────────────────────────────
    #
    # The tree is loaded once into Python tuples and then traversed
    # purely in Python — no struct unpacking on the hot path.
    # Typical tree size: a few hundred KB of tuples.

    def _load_rtree(self, idx_off: int):
        """Load entire R-tree rooted at *idx_off* into nested tuples."""
        if idx_off in self._rtree_cache:
            return self._rtree_cache[idx_off]
        magic = self._u('I', idx_off)[0]
        if magic != _RT_MAGIC:
            raise ValueError(f"Bad R-tree magic: {magic:#010x}")
        root = self._read_rtnode(idx_off + 48)
        self._rtree_cache[idx_off] = root
        return root

    def _read_rtnode(self, node: int):
        """Recursively read an R-tree node into a Python tuple.

        Returns
        -------
        For leaf:     (True,  [(sc,sb,ec,eb, d_off,d_sz), ...])
        For non-leaf: (False, [(sc,sb,ec,eb, child_node), ...])
        """
        is_leaf, _, count = self._u('BBH', node)
        off = node + 4
        if is_leaf:
            entries = []
            for _ in range(count):
                sc, sb, ec, eb, d_off, d_sz = self._u('IIIIQQ', off)
                entries.append((sc, sb, ec, eb, d_off, d_sz))
                off += 32
            return (True, entries)
        else:
            entries = []
            for _ in range(count):
                sc, sb, ec, eb, child_off = self._u('IIIIQ', off)
                child = self._read_rtnode(child_off)
                entries.append((sc, sb, ec, eb, child))
                off += 24
            return (False, entries)

    @staticmethod
    def _rtree_query(tree_node, qc: int, qs: int, qe: int,
                     out: List[Tuple[int, int]]):
        """Walk cached R-tree, appending (data_off, data_sz) to *out*."""
        is_leaf, entries = tree_node
        if is_leaf:
            for sc, sb, ec, eb, d_off, d_sz in entries:
                if qc > ec or (qc == ec and qs >= eb):
                    continue
                if qc < sc or (qc == sc and qe <= sb):
                    continue
                out.append((d_off, d_sz))
        else:
            for sc, sb, ec, eb, child in entries:
                if qc > ec or (qc == ec and qs >= eb):
                    continue
                if qc < sc or (qc == sc and qe <= sb):
                    continue
                BigWigReader._rtree_query(child, qc, qs, qe, out)

    def _rtree_find(self, idx_off: int, cid: int,
                    start: int, end: int) -> List[Tuple[int, int]]:
        """Return [(data_offset, data_size), ...] for overlapping blocks."""
        root = self._load_rtree(idx_off)
        blocks: List[Tuple[int, int]] = []
        self._rtree_query(root, cid, start, end, blocks)
        return blocks

    # ── block decompression & parsing ────────────────────────────────────

    def _decompress(self, offset: int, size: int) -> bytes:
        raw = bytes(self._mm[offset:offset + size])
        if self._uncompress_buf > 0:
            return zlib.decompress(raw)
        return raw

    def _parse_data_block(self, offset: int, size: int,
                          expect_cid: int | None = None,
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Parse a full-resolution data block.

        Returns
        -------
        starts, ends, values : np.ndarray (int64, int64, float64)
        """
        data = self._decompress(offset, size)
        e = self._end

        # block header — 24 bytes
        cid, blk_s, blk_e, step, span = struct.unpack_from(e + 'IIIII', data, 0)
        btype, _, n = struct.unpack_from(e + 'BBH', data, 20)

        if expect_cid is not None and cid != expect_cid:
            _empty = np.empty(0, dtype=np.int64)
            return _empty, _empty.copy(), np.empty(0, dtype=np.float64)

        body = data[24:]

        if btype == _BED_GRAPH:
            dt = np.dtype([('s', f'{e}u4'), ('e', f'{e}u4'), ('v', f'{e}f4')])
            rec = np.frombuffer(body, dtype=dt, count=n)
            return (rec['s'].astype(np.int64),
                    rec['e'].astype(np.int64),
                    rec['v'].astype(np.float64))

        if btype == _VAR_STEP:
            dt = np.dtype([('s', f'{e}u4'), ('v', f'{e}f4')])
            rec = np.frombuffer(body, dtype=dt, count=n)
            s = rec['s'].astype(np.int64)
            return s, s + int(span), rec['v'].astype(np.float64)

        if btype == _FIXED_STEP:
            v = np.frombuffer(body, dtype=f'{e}f4', count=n).astype(np.float64)
            s = np.arange(n, dtype=np.int64) * int(step) + int(blk_s)
            return s, s + int(span), v

        raise ValueError(f"Unknown block type {btype}")

    def _parse_zoom_block(self, offset: int, size: int) -> np.ndarray:
        """Parse zoom data → structured array of 32-byte summary records."""
        data = self._decompress(offset, size)
        e = self._end
        dt = np.dtype([
            ('cid', f'{e}u4'), ('s', f'{e}u4'), ('e', f'{e}u4'),
            ('valid', f'{e}u4'),
            ('mn', f'{e}f4'), ('mx', f'{e}f4'),
            ('sum', f'{e}f4'), ('ssq', f'{e}f4'),
        ])
        return np.frombuffer(data, dtype=dt)

    # ── zoom level selection ─────────────────────────────────────────────

    def _select_zoom(self, bin_size: float) -> Optional[int]:
        """Pick finest zoom with reductionLevel ≤ bin_size, or None."""
        best: Optional[int] = None
        for i, (rl, _, _) in enumerate(self._zooms):
            if rl <= bin_size:
                if best is None or rl > self._zooms[best][0]:
                    best = i
        return best

    # ── public API ───────────────────────────────────────────────────────

    def chroms(self) -> Dict[str, int]:
        """Return ``{chrom_name: size}`` dictionary."""
        return dict(self._chroms)

    # -- stats ---------------------------------------------------------

    def stats(self, chrom: str, start: int, end: int, *,
              n_bins: int = 1, nBins: int | None = None,
              stat: str = 'mean', type: str | None = None,
              exact: bool = False) -> list:
        """Binned summary statistics over a region.

        Parameters
        ----------
        chrom : str
        start, end : int  (0-based half-open)
        n_bins / nBins : int
        stat / type : mean | min | max | std | sum | coverage
        exact : bool  — force full-resolution (skip zoom levels)

        Returns
        -------
        list of float | None  (length *n_bins*)
        """
        if nBins is not None:
            n_bins = nBins
        if type is not None:
            stat = type

        cid = self._name_to_id.get(chrom)
        if cid is None:
            return [None] * n_bins

        if end <= start or n_bins <= 0:
            return [None] * n_bins

        bin_size = (end - start) / n_bins
        zoom_idx = None if exact else self._select_zoom(bin_size)

        if zoom_idx is not None:
            return self._stats_zoom(cid, start, end, n_bins, stat, zoom_idx)
        return self._stats_full(cid, start, end, n_bins, stat)

    def stats_array(self, chrom: str, start: int, end: int, *,
                    n_bins: int = 1, stat: str = 'mean',
                    exact: bool = True, missing: float = 0.0,
                    ) -> np.ndarray:
        """ndarray variant of :meth:`stats` — missing bins filled with *missing*.

        Matches the pybigtools adapter's fast-path signature so the
        worker in ``signal.py`` can call it uniformly. ``exact`` is
        accepted for API compatibility; the pure-Python reader always
        uses full-resolution data when it has it.
        """
        raw = self.stats(chrom, start, end, n_bins=n_bins,
                         stat=stat, exact=exact)
        return np.fromiter(
            (missing if v is None else v for v in raw),
            dtype=np.float64, count=len(raw),
        )

    def _stats_full(self, cid, start, end, n_bins, stat):
        bs = (end - start) / n_bins
        vc = np.zeros(n_bins)
        sd = np.zeros(n_bins)
        sq = np.zeros(n_bins)
        mn = np.full(n_bins, np.inf)
        mx = np.full(n_bins, -np.inf)

        for d_off, d_sz in self._rtree_find(self._full_index_off,
                                            cid, start, end):
            s, e, v = self._parse_data_block(d_off, d_sz, expect_cid=cid)
            if len(s) == 0:
                continue

            # clip to query range
            cs = np.maximum(s, start)
            ce = np.minimum(e, end)
            ok = cs < ce
            if not ok.any():
                continue
            cs, ce, v = cs[ok], ce[ok], v[ok]

            fb = np.clip(((cs - start) / bs).astype(np.int64), 0, n_bins - 1)
            lb = np.clip(((ce - 1 - start) / bs).astype(np.int64), 0, n_bins - 1)

            # ---- single-bin records (vectorised) ----
            single = fb == lb
            if single.any():
                bi = fb[single]
                bases = (ce[single] - cs[single]).astype(np.float64)
                vs = v[single]
                vc += np.bincount(bi, weights=bases, minlength=n_bins)
                sd += np.bincount(bi, weights=vs * bases, minlength=n_bins)
                sq += np.bincount(bi, weights=vs * vs * bases, minlength=n_bins)
                np.minimum.at(mn, bi, vs)
                np.maximum.at(mx, bi, vs)

            # ---- multi-bin records (loop, uncommon) ----
            for idx in np.where(~single)[0]:
                for b in range(int(fb[idx]), int(lb[idx]) + 1):
                    b_lo = start + b * bs
                    b_hi = start + (b + 1) * bs
                    ov = min(float(ce[idx]), b_hi) - max(float(cs[idx]), b_lo)
                    if ov > 0:
                        vc[b] += ov
                        sd[b] += v[idx] * ov
                        sq[b] += v[idx] ** 2 * ov
                        if v[idx] < mn[b]:
                            mn[b] = v[idx]
                        if v[idx] > mx[b]:
                            mx[b] = v[idx]

        return self._finalise(vc, sd, mn, mx, sq, stat, n_bins, bs)

    def _stats_zoom(self, cid, start, end, n_bins, stat, zoom_idx):
        _, _, idx_off = self._zooms[zoom_idx]
        bs = (end - start) / n_bins
        vc = np.zeros(n_bins)
        sd = np.zeros(n_bins)
        sq = np.zeros(n_bins)
        mn = np.full(n_bins, np.inf)
        mx = np.full(n_bins, -np.inf)

        for d_off, d_sz in self._rtree_find(idx_off, cid, start, end):
            rec = self._parse_zoom_block(d_off, d_sz)
            mask = ((rec['cid'] == cid) &
                    (rec['s'] < end) &
                    (rec['e'] > start))
            rec = rec[mask]
            if len(rec) == 0:
                continue

            rs = np.maximum(rec['s'].astype(np.int64), start)
            re = np.minimum(rec['e'].astype(np.int64), end)
            rec_span = (rec['e'] - rec['s']).astype(np.float64)
            frac = (re - rs).astype(np.float64) / rec_span

            fb = np.clip(((rs - start) / bs).astype(np.int64), 0, n_bins - 1)
            lb = np.clip(((re - 1 - start) / bs).astype(np.int64), 0, n_bins - 1)

            single = fb == lb
            if single.any():
                bi = fb[single]
                f = frac[single]
                vc += np.bincount(bi, weights=rec['valid'][single].astype(np.float64) * f, minlength=n_bins)
                sd += np.bincount(bi, weights=rec['sum'][single].astype(np.float64) * f, minlength=n_bins)
                sq += np.bincount(bi, weights=rec['ssq'][single].astype(np.float64) * f, minlength=n_bins)
                np.minimum.at(mn, bi, rec['mn'][single].astype(np.float64))
                np.maximum.at(mx, bi, rec['mx'][single].astype(np.float64))

            for idx in np.where(~single)[0]:
                r = rec[idx]
                total = float(r['e'] - r['s'])
                for b in range(int(fb[idx]), int(lb[idx]) + 1):
                    b_lo = start + b * bs
                    b_hi = start + (b + 1) * bs
                    ov_s = max(float(rs[idx]), b_lo)
                    ov_e = min(float(re[idx]), b_hi)
                    f = (ov_e - ov_s) / total
                    if f > 0:
                        vc[b] += float(r['valid']) * f
                        sd[b] += float(r['sum']) * f
                        sq[b] += float(r['ssq']) * f
                        mn[b] = min(mn[b], float(r['mn']))
                        mx[b] = max(mx[b], float(r['mx']))

        return self._finalise(vc, sd, mn, mx, sq, stat, n_bins, bs)

    @staticmethod
    def _finalise(vc, sd, mn, mx, sq, stat, n_bins, bin_size):
        out: list = []
        for i in range(n_bins):
            if vc[i] == 0:
                out.append(None)
                continue
            if stat == 'mean':
                out.append(sd[i] / vc[i])
            elif stat == 'max':
                out.append(float(mx[i]))
            elif stat == 'min':
                out.append(float(mn[i]))
            elif stat == 'std':
                m = sd[i] / vc[i]
                out.append(max(0.0, sq[i] / vc[i] - m * m) ** 0.5)
            elif stat == 'sum':
                out.append(float(sd[i]))
            elif stat in ('coverage', 'cov'):
                out.append(vc[i] / bin_size)
            else:
                raise ValueError(f"Unknown stat: {stat!r}")
        return out

    # -- base-pair resolution ------------------------------------------

    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        """Per-base float64 array.  NaN where no data exists."""
        cid = self._name_to_id.get(chrom)
        n = end - start
        out = np.full(n, np.nan, dtype=np.float64)
        if cid is None or n <= 0:
            return out

        for d_off, d_sz in self._rtree_find(self._full_index_off,
                                            cid, start, end):
            s, e, v = self._parse_data_block(d_off, d_sz, expect_cid=cid)
            if len(s) == 0:
                continue
            cs = np.maximum(s, start)
            ce = np.minimum(e, end)
            ok = cs < ce
            for lo, hi, val in zip(cs[ok] - start, ce[ok] - start, v[ok]):
                out[lo:hi] = val

        return out

    # -- raw intervals -------------------------------------------------

    def intervals(self, chrom: str, start: int, end: int,
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(starts, ends, values)`` for overlapping records."""
        cid = self._name_to_id.get(chrom)
        if cid is None:
            z = np.empty(0, dtype=np.int64)
            return z, z.copy(), np.empty(0, dtype=np.float64)

        all_s: List[np.ndarray] = []
        all_e: List[np.ndarray] = []
        all_v: List[np.ndarray] = []
        for d_off, d_sz in self._rtree_find(self._full_index_off,
                                            cid, start, end):
            s, e, v = self._parse_data_block(d_off, d_sz, expect_cid=cid)
            mask = (s < end) & (e > start)
            if mask.any():
                all_s.append(s[mask])
                all_e.append(e[mask])
                all_v.append(v[mask])

        if all_s:
            return (np.concatenate(all_s),
                    np.concatenate(all_e),
                    np.concatenate(all_v))
        z = np.empty(0, dtype=np.int64)
        return z, z.copy(), np.empty(0, dtype=np.float64)

    # -- lifecycle -----------------------------------------------------

    def close(self):
        if self._mm is not None:
            self._mm.close()
            self._mm = None
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self):
        return f"BigWigReader({self._path!r})"


# ── module-level convenience ─────────────────────────────────────────────────

def open(path: str | Path) -> BigWigReader:
    """Open a BigWig file for reading.  Thread-safe."""
    return BigWigReader(path)
