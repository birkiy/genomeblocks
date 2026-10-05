"""FASTA backends: where sequence comes from.

    src = open_fasta("hg38.fa")                 # path, dict, or an open handle
    src.sizes()                                 # {chrom: length}
    src.fetch("chr1", 1_000, 1_500)             # str, as stored (case kept)
    src.fetch_many(chroms, starts, ends)        # list[str]

Backends: ``genomeblocks`` (the default: a pure-Python reader on the
samtools ``.fai`` index — built next to the file when missing — so only the
requested bases are read), ``pysam`` (C, htslib), ``pyfaidx``, ``memory``
(the whole file read into a ``{chrom: str}`` dict) and ``biopython``
(``Bio.SeqIO.index``). All return the same bases. Besides a path,
``open_fasta`` takes a ``{chrom: str}`` dict, a Biopython ``SeqIO.index`` /
``to_dict`` mapping, a ``pyfaidx.Fasta`` or a ``pysam.FastaFile`` and uses it
as is. A gzip (not bgzip) FASTA cannot be indexed and is read into memory.
"""
from __future__ import annotations

from typing import Dict, List

from . import resolve


def read_fasta(path: str) -> Dict[str, str]:
    """Every record of a FASTA (plain or .gz) as ``{name: sequence}``."""
    import gzip
    out, name, chunks = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks)
                name, chunks = line[1:].split()[0] if line[1:].strip() else "", []
            else:
                chunks.append(line.rstrip("\n\r"))
    if name is not None:
        out[name] = "".join(chunks)
    return out


class _Source:
    backend = "?"

    def sizes(self) -> Dict[str, int]:
        raise NotImplementedError

    def fetch(self, chrom: str, start: int, end: int) -> str:
        raise NotImplementedError

    def fetch_many(self, chroms, starts, ends) -> List[str]:
        f = self.fetch
        return [f(c, int(s), int(e)) for c, s, e in zip(chroms, starts, ends)]

    def close(self):
        pass


class IndexedSource(_Source):
    """Random access through a samtools-style .fai (name, length, offset,
    line bases, line width) and a memory map of the file."""
    backend = "genomeblocks"

    def __init__(self, src):
        import mmap
        import os
        path = str(src)
        if path.endswith(".gz"):
            raise ValueError("a gzip FASTA cannot be indexed")
        fai = path + ".fai"
        if os.path.exists(fai) and os.path.getmtime(fai) >= os.path.getmtime(path):
            self.index = _read_fai(fai)
        else:
            self.index = _build_fai(path)
            try:
                with open(fai, "w") as f:
                    for k, (n, off, lb, lw) in self.index.items():
                        f.write(f"{k}\t{n}\t{off}\t{lb}\t{lw}\n")
            except OSError:
                pass                                      # read-only location: keep it in memory
        self._fh = open(path, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)

    def sizes(self):
        return {k: v[0] for k, v in self.index.items()}

    def fetch(self, chrom, start, end):
        n, off, lb, lw = self.index[chrom]
        start, end = max(0, int(start)), min(int(end), n)
        if end <= start:
            return ""
        a = off + (start // lb) * lw + start % lb
        b = off + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        return self._mm[a:b].replace(b"\n", b"").replace(b"\r", b"").decode("ascii")

    def close(self):
        self._mm.close()
        self._fh.close()


def _read_fai(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                out[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
    return out


def _build_fai(path):
    """Index a FASTA in one pass (lines of each record must have one length)."""
    out = {}
    name, length, off, lb, lw, last_short = None, 0, 0, 0, 0, False
    pos = 0
    with open(path, "rb") as f:
        for raw in f:
            if raw.startswith(b">"):
                if name is not None:
                    out[name] = (length, off, lb, lw)
                name = raw[1:].split()[0].decode() if raw[1:].strip() else ""
                length, off, lb, lw, last_short = 0, pos + len(raw), 0, 0, False
            elif name is not None:
                bases = len(raw.rstrip(b"\r\n"))
                if lb == 0:
                    lb, lw = bases, len(raw)
                elif last_short or bases > lb:
                    raise ValueError(f"{path}: uneven line lengths in {name}; cannot index")
                if bases < lb:
                    last_short = True
                length += bases
            pos += len(raw)
    if name is not None:
        out[name] = (length, off, lb, lw)
    return out


class MemorySource(_Source):
    backend = "memory"

    def __init__(self, src):
        self.d = read_fasta(str(src)) if isinstance(src, str) or hasattr(src, "__fspath__") else src
        self._sizes = None

    def sizes(self):
        if self._sizes is None:
            self._sizes = {k: len(v) for k, v in self.d.items()}
        return self._sizes

    def fetch(self, chrom, start, end):
        return self.d[chrom][start:end]

    def fetch_many(self, chroms, starts, ends):
        d = self.d
        return [d[c][s:e] for c, s, e in zip(chroms, starts.tolist() if hasattr(starts, "tolist") else starts,
                                             ends.tolist() if hasattr(ends, "tolist") else ends)]


class PysamSource(_Source):
    backend = "pysam"

    def __init__(self, src):
        import pysam
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self.h = pysam.FastaFile(str(src)) if self._own else src

    def sizes(self):
        return dict(zip(self.h.references, self.h.lengths))

    def fetch(self, chrom, start, end):
        return self.h.fetch(chrom, start, end)

    def close(self):
        if self._own:
            self.h.close()


class PyfaidxSource(_Source):
    backend = "pyfaidx"

    def __init__(self, src):
        import pyfaidx
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        self.h = pyfaidx.Fasta(str(src), as_raw=True, sequence_always_upper=False,
                               read_ahead=1 << 16) if self._own else src

    def sizes(self):
        return {k: len(v) for k, v in self.h.items()}

    def fetch(self, chrom, start, end):
        s = self.h[chrom][start:end]
        return s if isinstance(s, str) else str(s)

    def close(self):
        if self._own:
            self.h.close()


class BiopythonSource(_Source):
    """``Bio.SeqIO.index`` on a path, or any mapping of name -> SeqRecord."""
    backend = "biopython"

    def __init__(self, src):
        self._own = isinstance(src, str) or hasattr(src, "__fspath__")
        if self._own:
            from Bio import SeqIO
            src = SeqIO.index(str(src), "fasta")
        self.h = src
        self._cache = {}

    def _seq(self, chrom):
        s = self._cache.get(chrom)
        if s is None:
            s = self._cache[chrom] = str(self.h[chrom].seq)
        return s

    def sizes(self):
        return {k: len(self.h[k].seq) for k in self.h.keys()}

    def fetch(self, chrom, start, end):
        return self._seq(chrom)[start:end]

    def close(self):
        if self._own and hasattr(self.h, "close"):
            self.h.close()


_SOURCES = {"genomeblocks": IndexedSource, "memory": MemorySource, "pysam": PysamSource,
            "pyfaidx": PyfaidxSource, "biopython": BiopythonSource}


def _kind(src):
    """Backend of an already-open object, or None for a path."""
    if isinstance(src, _Source):
        return src.backend
    mod = type(src).__module__ or ""
    if mod.startswith("pysam"):
        return "pysam"
    if mod.startswith("pyfaidx"):
        return "pyfaidx"
    if hasattr(src, "keys") and hasattr(src, "__getitem__"):
        for k in src.keys():
            return "biopython" if hasattr(src[k], "seq") else "memory"
        return "memory"
    return None


def open_fasta(src, *, backend=None) -> _Source:
    """A sequence source for a FASTA path, a dict, or an open handle.

    Open handles and dicts are used as they are. A path goes to the requested
    backend; by default the genomeblocks reader, or memory for a gzip file
    (which cannot be indexed — every reader returns the same bases)."""
    if isinstance(src, _Source):
        return src
    kind = _kind(src)
    if kind is not None:
        return _SOURCES[kind](src)
    name = resolve("fasta", backend)
    if name == "genomeblocks" and backend is None:
        try:
            return IndexedSource(src)
        except ValueError:                                # gzip, or uneven lines
            return MemorySource(src)
    return _SOURCES[name](src)
