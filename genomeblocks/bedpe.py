"""BEDPE file operations and pair intersection utilities."""
from typing import Union, Optional, List, Tuple, Iterator, TYPE_CHECKING
from dataclasses import dataclass
import pandas as pd

if TYPE_CHECKING:
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


def _detect_pairs_format(filename: str) -> str:
    """Peek at the first non-comment line of a pairs file to identify the format.

    Returns 'allvalidpairs' (HiC-Pro: col 3 is a strand) or 'pairs' (pairtools 4DN: col 3 is chrom2).
    """
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7 and fields[3] in ("+", "-"):
                return "allvalidpairs"
            return "pairs"
    raise ValueError(f"No data lines in {filename}")


def read_pairs_chunks(
    filename: str,
    *,
    format: str = "auto",
    chunksize: int = 2_000_000,
) -> Iterator[pd.DataFrame]:
    """Stream a HiC-Pro .allValidPairs or pairtools .pairs file in chunks.

    Yields DataFrames with columns chrom1, pos1, chrom2, pos2 only — the
    minimum needed for window-vs-chromosome contact counting. chrom columns
    are categorical for memory efficiency on large files.

    Args:
        filename: path to .allValidPairs or .pairs (optionally gzipped — pandas auto-detects).
        format: 'allvalidpairs' (HiC-Pro), 'pairs' (pairtools/4DN), or 'auto' to detect.
        chunksize: rows per pandas chunk.
    """
    if format == "auto":
        format = _detect_pairs_format(filename)

    if format == "allvalidpairs":
        usecols = [1, 2, 4, 5]
    elif format == "pairs":
        usecols = [1, 2, 3, 4]
    else:
        raise ValueError(f"Unknown pairs format: {format!r}")

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
        format: 'allvalidpairs', 'pairs', or 'auto'.
        chunksize: rows per pandas chunk.

    Returns:
        DataFrame with columns chrom, start, end, uid, then either 'count' (when
        target_chrom is set) or one int column per partner chromosome.
    """
    import numpy as np
    from collections import defaultdict

    by_chrom: dict = {}
    for gi, locus in enumerate(loci):
        d = by_chrom.setdefault(locus.chrom, {"starts": [], "ends": [], "gidx": []})
        d["starts"].append(locus.start)
        d["ends"].append(locus.end)
        d["gidx"].append(gi)
    for d in by_chrom.values():
        order = np.argsort(d["starts"])
        d["starts"] = np.asarray(d["starts"], dtype=np.int64)[order]
        d["ends"] = np.asarray(d["ends"], dtype=np.int64)[order]
        d["gidx"] = np.asarray(d["gidx"], dtype=np.int64)[order]

    counts: dict = defaultdict(lambda: np.zeros(len(loci), dtype=np.int64))
    n_seen = 0

    for chunk in read_pairs_chunks(pairs_file, format=format, chunksize=chunksize):
        n_seen += len(chunk)

        for self_chrom_col, self_pos_col, other_chrom_col in (
            ("chrom1", "pos1", "chrom2"),
            ("chrom2", "pos2", "chrom1"),
        ):
            self_chroms = chunk[self_chrom_col].astype(str).to_numpy()
            self_pos = chunk[self_pos_col].to_numpy(dtype=np.int64, copy=False)
            other_chroms = chunk[other_chrom_col].astype(str).to_numpy()

            if target_chrom is not None:
                m = other_chroms == target_chrom
                if not m.any():
                    continue
                self_chroms = self_chroms[m]
                self_pos = self_pos[m]
                other_chroms = other_chroms[m]

            for chrom, d in by_chrom.items():
                m = self_chroms == chrom
                if not m.any():
                    continue
                positions = self_pos[m]
                partner = other_chroms[m]

                idx = np.searchsorted(d["starts"], positions, side="right") - 1
                ok = idx >= 0
                if ok.any():
                    safe = idx[ok]
                    ok2 = positions[ok] < d["ends"][safe]
                    final = np.zeros_like(ok)
                    final[ok] = ok2
                    if not final.any():
                        continue
                    gidx_hit = d["gidx"][idx[final]]
                    partner_hit = partner[final]

                    if target_chrom is not None:
                        np.add.at(counts[target_chrom], gidx_hit, 1)
                    else:
                        for up in np.unique(partner_hit):
                            pm = partner_hit == up
                            np.add.at(counts[str(up)], gidx_hit[pm], 1)

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
        out["count"] = counts[target_chrom]
    else:
        for chrom in sorted(counts.keys()):
            out[chrom] = counts[chrom]
    return out


# Attach to Loci class
def _attach_to_loci():
    """Attach pair_to_bed and count_pairs methods to Loci class."""
    from .loci import Loci
    Loci.pair_to_bed = pair_to_bed
    Loci.count_pairs = count_pairs


# Auto-attach when module is imported
try:
    _attach_to_loci()
except ImportError:
    # Module not yet fully initialized, will be attached later
    pass
