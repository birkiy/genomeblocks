"""BEDPE file operations and pair intersection utilities."""
from typing import Union, Optional, List, Tuple
from dataclasses import dataclass
import pandas as pd


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


# Attach to Loci class
def _attach_to_loci():
    """Attach pair_to_bed method to Loci class."""
    from .loci import Loci
    Loci.pair_to_bed = pair_to_bed


# Auto-attach when module is imported
try:
    _attach_to_loci()
except ImportError:
    # Module not yet fully initialized, will be attached later
    pass
