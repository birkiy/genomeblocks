"""Loci container and interval operations."""
from typing import Iterable, Optional, Union, List, Dict
from collections import defaultdict

# Local imports are delayed to avoid circular imports
from .locus import Locus
import cgranges as cr

class Loci(list):

    def __init__(self, iterable: Iterable[Locus] = (), *, filename: Optional[str] = None):
        super().__init__(iterable)
        self.filename = filename
        self._uids: Optional[dict] = None
        self._cgr: Optional[cr.cgranges] = None

    def _build_cgr(self) -> None:
        idx = cr.cgranges()
        for i, l in enumerate(self):
            idx.add(l.chrom, int(l.start), int(l.end), i)
        idx.index()
        self._cgr = idx

    @property
    def cgr(self) -> cr.cgranges:
        if self._cgr is None: self._build_cgr()
        return self._cgr


    @property
    def uids(self) -> dict:
        if self._uids is None: self._uids = {l.uid: i for i, l in enumerate(self)}
        return self._uids

    def to_pyranges(s, names=None):
        import pyranges as pr
        df = s.to_frame().rename(columns={'Chr': 'Chromosome'})
        if names is None: names = [l.uid for l in s]
        df['Name'] = names
        return pr.PyRanges(df)

    def to_frame(s, names=None):
        import pandas as pd
        if names is None: names = [l.uid for l in s]
        return pd.DataFrame([
            {"Chr": l.chrom, "Start": l.start, "End": l.end, 'Strand': l.strand, 'Name': names[i]}
            for i, l in enumerate(s)
        ])

    def to_bed(s, path=None):
        lines = []
        for l in s:
            lines.append(f"{l.chrom}\t{l.start}\t{l.end}\t{l.uid}\t0\t{l.strand}")
        content = "\n".join(lines)
        if path is not None:
            with open(path, 'w') as f:
                f.write(content)


    def __and__(s, o: "Loci") -> "Loci": return s.intersect(o)
    def __add__(s, o: "Loci") -> "Loci": return Loci([*s, *o])
    def __sub__(s, o: "Loci") -> "Loci": return s.difference(o)
    def __truediv__(s, o: "Loci") -> "Loci": return s.difference(o)
    def __or__(s, o: "Loci") -> "Loci": return Loci([*s, *o])
    def __xor__(s, o: "Loci") -> "Loci": return (s - o) + (o - s)

    def copy(self) -> "Loci":
        return Loci(list(self), filename=self.filename)

    def __str__(self):
        return f"Loci(n={len(self)})"
    __repr__ = __str__

    def __getstate__(self):
        s = self.__dict__.copy()
        s.pop("_cgr", None)
        return s

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._cgr = None

    def __getitem__(self, key: Union[int, str, slice]):
        if isinstance(key, int) or isinstance(key, slice):
            return super().__getitem__(key)
        elif isinstance(key, str):
            try:
                index = self.uids[key]
                return super().__getitem__(index)
            except KeyError:
                raise KeyError(f"No Locus with UID '{key}' found.")
        else:
            raise TypeError(f"Invalid key type: {type(key)}. Expected int, str, or slice.")


# helper functions attached to Loci (kept as in original)

def _auto_format(p: str) -> str:
    if p.endswith(".gtf"): return "Features"
    if p.endswith((".gff", ".gff3")): return "Features"
    if p.endswith(".narrowPeak"): return "Loci"
    return "Loci"

@classmethod
def make(cls, filename: str, filetype: Optional[str] = None) -> Loci:
    if filetype is None: filetype = _auto_format(filename)
    if filetype == "Features": raise NotImplementedError("GTF/GFF parsing not implemented. Please use Features.make() instead.")
    elif filetype != 'Loci': raise ValueError(f"Unsupported file type for Loci!")
    L = Loci(filename=filename)
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            strand = fields[5] if len(fields) > 5 else "."
            L.append(Locus(chrom, start, end, strand))
    return L


@classmethod
def tile(cls, chrom: str, size: int, chromsizes) -> "Loci":
    """Build uniform `size`-bp tiles across `chrom`.

    chromsizes can be:
      - int: chromosome length directly
      - dict: {chrom: length}
      - str: path to a UCSC .chrom.sizes file (chrom\\tlength)
      - any object with a .chromsizes mapping (e.g. cooler.Cooler)
    """
    if isinstance(chromsizes, int):
        clen = chromsizes
    elif isinstance(chromsizes, dict):
        clen = int(chromsizes[chrom])
    elif isinstance(chromsizes, str):
        sizes = {}
        with open(chromsizes) as f:
            for line in f:
                if not line.strip() or line.startswith("#"): continue
                parts = line.split()
                sizes[parts[0]] = int(parts[1])
        clen = sizes[chrom]
    else:
        clen = int(chromsizes.chromsizes[chrom])
    return cls(
        Locus(chrom, s, min(s + size, clen))
        for s in range(0, clen, size)
    )


@classmethod
def tile_genome(cls, chromsizes, size: int, chroms: Optional[List[str]] = None) -> "Loci":
    """Tile every chromosome in `chromsizes` at uniform `size` bp.

    chromsizes: dict, .chrom.sizes path, or any object with a .chromsizes mapping.
    chroms: optional whitelist (and ordering) of chroms to include.
    """
    if isinstance(chromsizes, dict):
        sizes = {k: int(v) for k, v in chromsizes.items()}
    elif isinstance(chromsizes, str):
        sizes = {}
        with open(chromsizes) as f:
            for line in f:
                if not line.strip() or line.startswith("#"): continue
                parts = line.split()
                sizes[parts[0]] = int(parts[1])
    else:
        sizes = {k: int(v) for k, v in dict(chromsizes.chromsizes).items()}
    if chroms is None:
        chroms = list(sizes.keys())
    out = cls()
    for c in chroms:
        clen = sizes[c]
        for s in range(0, clen, size):
            out.append(Locus(c, s, min(s + size, clen)))
    return out


def subloci(s, uids: List[str]):
    sub = Loci(s[s.uids[u]] for u in uids)
    return sub

@classmethod
def from_frame(cls, df) -> "Loci":
    """Build a Loci from the first 3 columns of a DataFrame (chrom, start, end).

    Column names are ignored — positional only. Strand defaults to '.'.
    """
    return cls(
        Locus(str(r[0]), int(r[1]), int(r[2]))
        for r in df.iloc[:, :3].itertuples(index=False)
    )


Loci.make = make
Loci.from_frame = from_frame
Loci.tile = tile
Loci.tile_genome = tile_genome
Loci.subloci = subloci


def overlaps(s, key:Union[Locus, str], start:Optional[int]=None, end:Optional[int]=None):
    if isinstance(key, Locus):
        start = key.start
        end = key.end
        key = key.chrom
    return Loci(s[i] for *_ , i in s.cgr.overlap(key, start, end))


def intersect(s, o: Loci) -> Loci:
    return Loci([l for l in s if any(True for *_ , _ in o.cgr.overlap(l.chrom, l.start, l.end))])


def difference(s, o: Loci) -> Loci:
    return Loci([l for l in s if not any(True for *_ , _ in o.cgr.overlap(l.chrom, l.start, l.end))])


def slop(s, n:int) -> Loci:
    return Loci([Locus(l.chrom, max(0, l.start - n), l.end + n, l.strand) for l in s])


def sort(s) -> Loci:
    out = Loci(sorted(s))
    out.uids
    return out


def merge(s) -> Loci:
    if not s: return Loci()
    sorted_loci = sort(s)
    merged = [sorted_loci[0].copy()]
    for l in sorted_loci[1:]:
        last = merged[-1]
        if last.chrom != l.chrom: merged.append(l.copy()); continue
        if last.overlaps(l) or last.end == l.start:
            last.end = max(last.end, l.end)
        else:
            merged.append(l.copy())
    return Loci(merged)


def nearest(s: Loci, o: Loci, s_names=None, o_names=None):
    pr_s = s.to_pyranges(s_names)
    pr_o = o.to_pyranges(o_names)
    return pr_s.nearest(pr_o).df.rename(columns={'Chromosome': 'Chr'})


def map(s: "Loci", o: "Loci") -> Dict[str, List[str]]:
    """Overlap-based mapper: returns {a.uid: [b.uid, ...]} for each a in self
    against loci in o. Multi-overlap preserved; empty list if no hit.
    """
    out: Dict[str, List[str]] = {}
    for a in s:
        out[a.uid] = [o[j].uid for *_, j in o.cgr.overlap(a.chrom, a.start, a.end)]
    return out


_LIFTOVER_CACHE: Dict[str, object] = {}


def liftover(s: Loci, chain_file: str, *, min_match: float = 0.95, verbose: bool = True) -> Loci:
    """Lift to another assembly via a UCSC chain file (pyliftover wrapper).

    Mirrors UCSC liftOver defaults: per-base tolerance of ``1 - min_match``.
    When an endpoint falls in a chain gap, it walks inward up to that budget
    before giving up — so small gaps at the edges (very common) don't kill
    the whole interval, matching UCSC's base-fraction semantics.
    """
    from pyliftover import LiftOver
    lo = _LIFTOVER_CACHE.get(chain_file)
    if lo is None:
        lo = LiftOver(chain_file)
        _LIFTOVER_CACHE[chain_file] = lo

    def _walk(chrom, pos, q_strand, step, budget):
        for off in range(budget + 1):
            hits = lo.convert_coordinate(chrom, pos + step * off, q_strand)
            if hits:
                return hits[0], off
        return None, None

    out = []
    n_dropped = 0
    for l in s:
        length = l.end - l.start
        if length <= 0:
            n_dropped += 1
            continue
        q_strand = l.strand if l.strand in ('+', '-') else '+'
        budget = max(0, int((1.0 - min_match) * length))
        s_hit, s_off = _walk(l.chrom, l.start,     q_strand, +1, budget)
        e_hit, e_off = _walk(l.chrom, l.end - 1,   q_strand, -1, budget)
        if s_hit is None or e_hit is None:
            n_dropped += 1
            continue
        sc, sp, sst, _ = s_hit
        ec, ep, est, _ = e_hit
        if sc != ec or sst != est:
            n_dropped += 1
            continue
        if s_off + e_off > budget:
            n_dropped += 1
            continue
        new_start, new_end = (sp, ep + 1) if sp <= ep else (ep, sp + 1)
        out_strand = sst if l.strand in ('+', '-') else '.'
        out.append(Locus(sc, new_start, new_end, out_strand))
    if verbose:
        print(f"[INFO] liftover: {len(s)} → {len(out)} ({n_dropped} dropped, min_match={min_match})")
    return Loci(out)


Loci.overlaps = overlaps
Loci.intersect = intersect
Loci.difference = difference
Loci.slop = slop
Loci.sort = sort
Loci.merge = merge
Loci.nearest = nearest
Loci.map = map
Loci.liftover = liftover

# Import domain modules to attach their methods to Loci. Each *_draw module is
# imported alongside its processing module so the plotting methods attach too.
from . import signal  # noqa: F401
from . import signal_draw  # noqa: F401
from . import motifs  # noqa: F401
from . import bedpe  # noqa: F401
from . import atlas  # noqa: F401
