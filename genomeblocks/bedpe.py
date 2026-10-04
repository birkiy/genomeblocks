"""BEDPE file operations and pair intersection utilities."""
from __future__ import annotations
from typing import Union, Optional, List, Tuple, Iterator, TYPE_CHECKING
from dataclasses import dataclass

if TYPE_CHECKING:
    import pandas as pd
    from .loci import Loci


@dataclass
class Pair:
    """Represents a paired-end genomic interval from BEDPE format."""
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str = ""
    score: float = 0.0
    strand1: str = "."
    strand2: str = "."
    extra_fields: Tuple = ()
    
    @property
    def mid1(self) -> int:
        """Midpoint of first interval."""
        return (self.start1 + self.end1) // 2
    
    @property
    def mid2(self) -> int:
        """Midpoint of second interval."""
        return (self.start2 + self.end2) // 2
    
    @property
    def distance(self) -> int:
        """Distance between midpoints (only if on same chromosome)."""
        if self.chrom1 == self.chrom2:
            return abs(self.mid2 - self.mid1)
        return float('inf')
    
    def __str__(self):
        return f"Pair({self.chrom1}:{self.start1}-{self.end1} <-> {self.chrom2}:{self.start2}-{self.end2})"
    
    __repr__ = __str__


def read_bedpe(
    filename: str, 
    *, 
    min_score: Optional[float] = None,
    max_distance: Optional[float] = None,
    verbose: bool = True
) -> List[Pair]:
    """Read a BEDPE file and return a list of Pair objects.
    
    Args:
        filename: Path to BEDPE file
        min_score: Optional minimum score filter
        max_distance: Optional maximum distance filter (for intra-chromosomal pairs)
        verbose: Print loading statistics
    
    Returns:
        List of Pair objects
    """
    pairs = []
    skipped = 0
    
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                skipped += 1
                continue
            
            try:
                chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
                chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
                
                # Optional fields
                name = fields[6] if len(fields) > 6 else ""
                score = float(fields[7]) if len(fields) > 7 else 0.0
                strand1 = fields[8] if len(fields) > 8 else "."
                strand2 = fields[9] if len(fields) > 9 else "."
                extra_fields = tuple(fields[10:]) if len(fields) > 10 else ()
                
                pair = Pair(
                    chrom1, start1, end1,
                    chrom2, start2, end2,
                    name, score, strand1, strand2, extra_fields
                )
                
                # Apply filters
                if min_score is not None and pair.score < min_score:
                    skipped += 1
                    continue
                
                if max_distance is not None and pair.distance > max_distance:
                    skipped += 1
                    continue
                
                pairs.append(pair)
                
            except (ValueError, IndexError):
                skipped += 1
                continue
    
    if verbose:
        print(f"[INFO] Loaded {len(pairs)} pairs from {filename}")
        if skipped > 0:
            print(f"[INFO] Skipped {skipped} lines (malformed or filtered)")
    
    return pairs


def pair_to_bed(
    loci,
    bedpe: Union[str, List[Pair]],
    *,
    r: int = 0,
    either: bool = True,
    both: bool = False,
    min_score: Optional[float] = None,
    max_distance: Optional[float] = None,
    verbose: bool = True
) -> List[Pair]:
    """Find pairs from BEDPE that intersect with given loci.
    
    Similar to bedtools pairtoBed, returns pairs where at least one anchor
    (or both, depending on parameters) overlaps with the provided loci.
    
    Args:
        loci: Loci object or single Locus to intersect with
        bedpe: Path to BEDPE file or list of Pair objects
        r: Slop/extension radius around pair anchors (default: 0)
        either: Return pairs where at least one anchor overlaps (default: True)
        both: Return pairs where both anchors overlap (default: False)
        min_score: Minimum score threshold for pairs
        max_distance: Maximum distance between pair anchors
        verbose: Print statistics
    
    Returns:
        List of Pair objects that satisfy the intersection criteria
    """
    from .loci import Loci
    from .locus import Locus
    
    # Ensure loci is a Loci object
    if isinstance(loci, Locus):
        loci = Loci([loci])
    elif not isinstance(loci, Loci):
        loci = Loci(loci)
    
    # Load pairs if needed
    if isinstance(bedpe, str):
        pairs = read_bedpe(
            bedpe, 
            min_score=min_score, 
            max_distance=max_distance, 
            verbose=verbose
        )
    else:
        pairs = bedpe
    
    # Find overlapping pairs
    overlapping = []
    
    for pair in pairs:
        # Check overlap for anchor 1
        anchor1_overlaps = any(
            True for *_, _ in loci.cgr.overlap(
                pair.chrom1, 
                max(0, pair.start1 - r), 
                pair.end1 + r
            )
        )
        
        # Check overlap for anchor 2
        anchor2_overlaps = any(
            True for *_, _ in loci.cgr.overlap(
                pair.chrom2, 
                max(0, pair.start2 - r), 
                pair.end2 + r
            )
        )
        
        # Apply logic
        if both and anchor1_overlaps and anchor2_overlaps:
            overlapping.append(pair)
        elif either and (anchor1_overlaps or anchor2_overlaps):
            overlapping.append(pair)
    
    if verbose:
        print(f"[INFO] Found {len(overlapping)}/{len(pairs)} pairs overlapping loci")
    
    return overlapping


def pairs_to_frame(pairs: List[Pair]) -> pd.DataFrame:
    """Convert list of Pair objects to a pandas DataFrame.
    
    Args:
        pairs: List of Pair objects
    
    Returns:
        DataFrame with standard BEDPE columns
    """
    import pandas as pd
    return pd.DataFrame([
        {
            'chrom1': p.chrom1, 'start1': p.start1, 'end1': p.end1,
            'chrom2': p.chrom2, 'start2': p.start2, 'end2': p.end2,
            'name': p.name, 'score': p.score,
            'strand1': p.strand1, 'strand2': p.strand2,
        }
        for p in pairs
    ])


def pairs_to_bedpe(pairs: List[Pair], filename: str) -> None:
    """Write list of Pair objects to a BEDPE file.
    
    Args:
        pairs: List of Pair objects
        filename: Output file path
    """
    with open(filename, 'w') as f:
        for p in pairs:
            fields = [
                p.chrom1, str(p.start1), str(p.end1),
                p.chrom2, str(p.start2), str(p.end2),
                p.name, str(p.score), p.strand1, p.strand2
            ]
            if p.extra_fields:
                fields.extend(p.extra_fields)
            f.write('\t'.join(fields) + '\n')


# (chrom1_col, pos1_col, chrom2_col, pos2_col) — 0-indexed positions in the row.
PAIRS_FORMAT_COLUMNS = {
    "allvalidpairs": (1, 2, 4, 5),  # HiC-Pro: readID chr1 pos1 strand1 chr2 pos2 strand2 ...
    "pairs":         (1, 2, 3, 4),  # pairtools / 4DN: readID chrom1 pos1 chrom2 pos2 strand1 strand2 ...
    "juicer":        (2, 3, 6, 7),  # Juicer medium: readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2
}


def _detect_pairs_format(filename: str) -> str:
    """Peek at the first non-comment line of a pairs file to identify the format."""
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8 and fields[1] in ("+", "-") and fields[5] in ("+", "-"):
                return "juicer"
            if len(fields) >= 7 and fields[3] in ("+", "-"):
                return "allvalidpairs"
            return "pairs"
    raise ValueError(f"No data lines in {filename}")


def read_pairs_chunks(
    filename: str,
    *,
    format: str = "auto",
    columns: Optional[Tuple[int, int, int, int]] = None,
    chunksize: int = 2_000_000,
) -> Iterator[pd.DataFrame]:
    """Stream a pairs / allValidPairs file in chunks.

    Yields DataFrames with columns chrom1, pos1, chrom2, pos2 only — the
    minimum needed for window-vs-chromosome contact counting. chrom columns
    are categorical for memory efficiency on large files.

    Args:
        filename: path to the pairs file (optionally gzipped — pandas auto-detects).
        format: one of 'allvalidpairs' (HiC-Pro), 'pairs' (pairtools/4DN),
            'juicer' (Juicer medium: readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2),
            or 'auto' to detect.
        columns: explicit (chrom1_col, pos1_col, chrom2_col, pos2_col) 0-indexed
            override for non-standard layouts. Wins over `format` if provided.
        chunksize: rows per pandas chunk.
    """
    import pandas as pd
    if columns is not None:
        usecols = list(columns)
    else:
        if format == "auto":
            format = _detect_pairs_format(filename)
        if format not in PAIRS_FORMAT_COLUMNS:
            raise ValueError(
                f"Unknown pairs format: {format!r}. "
                f"Known: {sorted(PAIRS_FORMAT_COLUMNS)}, or pass columns=(c1,p1,c2,p2)."
            )
        usecols = list(PAIRS_FORMAT_COLUMNS[format])

    reader = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        comment="#",
        usecols=usecols,
        names=["chrom1", "pos1", "chrom2", "pos2"],
        dtype={"chrom1": "category", "pos1": "int64",
               "chrom2": "category", "pos2": "int64"},
        chunksize=chunksize,
        engine="c",
    )
    yield from reader


def count_pairs(
    loci: "Loci",
    pairs_file: str,
    *,
    target_chrom: Optional[str] = None,
    format: str = "auto",
    columns: Optional[Tuple[int, int, int, int]] = None,
    chunksize: int = 2_000_000,
    verbose: bool = True,
) -> pd.DataFrame:
    """Count valid pairs landing in each window of `loci`, broken down by partner chromosome.

    Single streaming pass over the pairs/allValidPairs file. For each pair (cA, pA, cB, pB):
      - if pA is inside a window W on chromosome cA, increment cell (W, cB)
      - if pB is inside a window W on chromosome cB, increment cell (W, cA)

    Both anchors are scanned independently, so a cis pair (both ends on the
    same chrom and both inside windows) contributes to two cells — once per anchor.

    Args:
        loci: Loci of windows. Windows on a given chromosome must be non-overlapping.
        pairs_file: .allValidPairs (HiC-Pro) or .pairs (pairtools/4DN).
        target_chrom: if given, only contacts whose partner is on this chrom are
            counted, and the output has a single 'count' column. If None, every
            partner chrom seen in the file gets its own column (whole-genome scan
            in one pass).
        format: 'allvalidpairs', 'pairs', 'juicer', or 'auto'.
        columns: explicit (chrom1_col, pos1_col, chrom2_col, pos2_col) 0-indexed
            override for non-standard layouts.
        chunksize: rows per pandas chunk.

    Returns:
        DataFrame with columns chrom, start, end, uid, then either 'count' (when
        target_chrom is set) or one int column per partner chromosome.
    """
    import numpy as np
    import pandas as pd

    n = len(loci)
    index = _window_index(loci)
    names: dict = {}                       # partner chrom -> column of `mat`
    mat = np.zeros((n, 0), dtype=np.int64)
    target = np.zeros(n, dtype=np.int64)
    n_seen = 0

    for chunk in read_pairs_chunks(pairs_file, format=format, columns=columns, chunksize=chunksize):
        n_seen += len(chunk)

        for self_chrom_col, self_pos_col, other_chrom_col in (
            ("chrom1", "pos1", "chrom2"),
            ("chrom2", "pos2", "chrom1"),
        ):
            w = _locate_codes(index, chunk[self_chrom_col],
                              chunk[self_pos_col].to_numpy(dtype=np.int64, copy=False))
            other = chunk[other_chrom_col].astype("category")
            other_names = [str(c) for c in other.cat.categories] + ["nan"]  # code -1 = NaN
            if target_chrom is not None:
                hit = (w >= 0) & np.array([c == target_chrom for c in other_names])[other.cat.codes.to_numpy()]
                target += np.bincount(w[hit], minlength=n)
                continue
            p_lut = np.array([names.setdefault(c, len(names)) for c in other_names], dtype=np.int64)
            if len(names) > mat.shape[1]:
                mat = np.pad(mat, ((0, 0), (0, len(names) - mat.shape[1])))
            hit = w >= 0
            p = p_lut[other.cat.codes.to_numpy()][hit]
            P = mat.shape[1]
            mat += np.bincount(w[hit] * P + p, minlength=n * P).reshape(n, P)

        if verbose:
            print(f"[INFO] processed {n_seen:,} pairs", end="\r")

    if verbose:
        print(f"\n[INFO] processed {n_seen:,} pairs total")

    out = pd.DataFrame({
        "chrom": [l.chrom for l in loci],
        "start": [l.start for l in loci],
        "end":   [l.end   for l in loci],
        "uid":   [l.uid   for l in loci],
    })
    if target_chrom is not None:
        out["count"] = target
    else:
        for chrom in sorted(names):
            col = mat[:, names[chrom]]
            if col.any():                  # only partners with at least one hit
                out[chrom] = col
    return out


_KEY_SHIFT = 40     # window key = chrom code << 40 | start; positions < 1.1e12


def _window_index(loci):
    """Windows on one genome-wide sorted key for vectorised lookup.

    Returns ``(chrom -> code, keys, codes, ends, rows)`` sorted by
    (chrom, start); ``rows`` maps back to positions in ``loci``.
    """
    import numpy as np
    n = len(loci)
    cid: dict = {}
    code = np.fromiter((cid.setdefault(l.chrom, len(cid)) for l in loci), dtype=np.int64, count=n)
    start = np.fromiter((l.start for l in loci), dtype=np.int64, count=n)
    end = np.fromiter((l.end for l in loci), dtype=np.int64, count=n)
    order = np.lexsort((start, code))
    return cid, (code[order] << _KEY_SHIFT) + start[order], code[order], end[order], order


def _locate_codes(index, chroms, pos):
    """Row in ``loci`` of the window holding each (chrom, pos); -1 if none.

    ``chroms`` is a (categorical) Series, so chromosome names are matched
    once per category rather than once per row.
    """
    import numpy as np
    cid, keys, wcode, wend, rows = index
    chroms = chroms.astype("category")
    lut = np.array([cid.get(str(c), -1) for c in chroms.cat.categories] + [-1], dtype=np.int64)
    c = lut[chroms.cat.codes.to_numpy()]               # NaN code -1 -> trailing -1
    i = np.searchsorted(keys, (c << _KEY_SHIFT) + pos, side="right") - 1
    j = np.maximum(i, 0)
    ok = (c >= 0) & (i >= 0) & (wcode[j] == c) & (pos < wend[j])
    return np.where(ok, rows[j], -1)


def count_pairs_2d(
    loci_a: "Loci",
    pairs_file: str,
    *,
    loci_b: Optional["Loci"] = None,
    format: str = "auto",
    columns: Optional[Tuple[int, int, int, int]] = None,
    chunksize: int = 2_000_000,
    verbose: bool = True,
):
    """Window-to-window pair counts. Returns a scipy.sparse.csr_matrix of shape
    (len(loci_a), len(loci_b)).

    Single streaming pass. For each pair (cA, pA, cB, pB), both directions are tried:
      - if pA in loci_a window i and pB in loci_b window j → M[i, j] += 1
      - if pB in loci_a window i and pA in loci_b window j → M[i, j] += 1

    With loci_b=None (default), uses loci_a on both axes — the resulting matrix is
    symmetric (each pair contributes to both M[i,j] and M[j,i]).

    Caveat for cis self-cells: a pair whose two anchors fall in the same window W
    contributes M[W,W] += 2 in the symmetric (loci_b=None) case (once per direction).
    Negligible at 50kb resolution.

    Args:
        loci_a: Loci for the rows (must be non-overlapping per chrom).
        pairs_file: .allValidPairs / .pairs / juicer-medium file.
        loci_b: Loci for the columns; if None, uses loci_a (all-vs-all).
        format, columns, chunksize: forwarded to read_pairs_chunks.
        verbose: print progress.

    Returns:
        scipy.sparse.csr_matrix of shape (len(loci_a), len(loci_b)) with int64 counts.
        Use pair_2d_block / pair_2d_to_frame to slice or materialize results.
    """
    import numpy as np
    import scipy.sparse as sp

    same = loci_b is None
    if same:
        loci_b = loci_a

    index_a = _window_index(loci_a)
    index_b = index_a if same else _window_index(loci_b)

    n_a, n_b = len(loci_a), len(loci_b)

    rows_acc: list = []
    cols_acc: list = []
    data_acc: list = []
    n_seen = 0

    for chunk in read_pairs_chunks(pairs_file, format=format, columns=columns, chunksize=chunksize):
        n_seen += len(chunk)

        p1 = chunk["pos1"].to_numpy(dtype=np.int64, copy=False)
        p2 = chunk["pos2"].to_numpy(dtype=np.int64, copy=False)
        a1 = _locate_codes(index_a, chunk["chrom1"], p1)
        a2 = _locate_codes(index_a, chunk["chrom2"], p2)
        if same:
            b1, b2 = a1, a2
        else:
            b1 = _locate_codes(index_b, chunk["chrom1"], p1)
            b2 = _locate_codes(index_b, chunk["chrom2"], p2)

        # both directions: anchor 1 in a row window and anchor 2 in a column
        # window, then the reverse
        keys = []
        for ga, gb in ((a1, b2), (a2, b1)):
            ok = (ga >= 0) & (gb >= 0)
            keys.append(ga[ok] * n_b + gb[ok])
        keys = np.concatenate(keys)
        if keys.size:
            uk, uc = np.unique(keys, return_counts=True)
            rows_acc.append(uk // n_b)
            cols_acc.append(uk % n_b)
            data_acc.append(uc.astype(np.int64))

        if verbose:
            print(f"[INFO] processed {n_seen:,} pairs", end="\r")

    if verbose:
        print(f"\n[INFO] processed {n_seen:,} pairs total")

    if not data_acc:
        return sp.csr_matrix((n_a, n_b), dtype=np.int64)

    rows = np.concatenate(rows_acc)
    cols = np.concatenate(cols_acc)
    data = np.concatenate(data_acc)
    coo = sp.coo_matrix((data, (rows, cols)), shape=(n_a, n_b), dtype=np.int64)
    coo.sum_duplicates()
    return coo.tocsr()


def pair_2d_block(mat, loci_a, loci_b, chrom_a: str, chrom_b: str):
    """Extract the dense submatrix for one chromosome pair from a count_pairs_2d result.

    Returns (block, windows_a_on_chrom_a, windows_b_on_chrom_b).
    """
    import numpy as np
    from .loci import Loci
    rows = [i for i, l in enumerate(loci_a) if l.chrom == chrom_a]
    cols = [j for j, l in enumerate(loci_b) if l.chrom == chrom_b]
    if not rows or not cols:
        return np.zeros((len(rows), len(cols)), dtype=np.int64), Loci(), Loci()
    block = mat[rows, :][:, cols].toarray()
    return block, Loci(loci_a[i] for i in rows), Loci(loci_b[j] for j in cols)


def pair_2d_to_frame(mat, loci_a, loci_b) -> pd.DataFrame:
    """Convert a count_pairs_2d sparse matrix to a long-format DataFrame of non-zero cells:
    chrom1 start1 end1 chrom2 start2 end2 count.
    """
    import pandas as pd
    coo = mat.tocoo()
    rows = [int(i) for i in coo.row]
    cols = [int(j) for j in coo.col]
    return pd.DataFrame({
        "chrom1": [loci_a[i].chrom for i in rows],
        "start1": [loci_a[i].start for i in rows],
        "end1":   [loci_a[i].end   for i in rows],
        "chrom2": [loci_b[j].chrom for j in cols],
        "start2": [loci_b[j].start for j in cols],
        "end2":   [loci_b[j].end   for j in cols],
        "count":  coo.data,
    })


# Attach to Loci class
def _attach_to_loci():
    """Attach pair_to_bed, count_pairs, and count_pairs_2d to Loci class."""
    from .loci import Loci
    Loci.pair_to_bed = pair_to_bed
    Loci.count_pairs = count_pairs
    Loci.count_pairs_2d = count_pairs_2d


# Auto-attach when module is imported
try:
    _attach_to_loci()
except ImportError:
    # Module not yet fully initialized, will be attached later
    pass
