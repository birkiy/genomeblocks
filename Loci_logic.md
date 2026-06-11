# Loci

`Loci` is a list of `Locus` objects — a sorted, queryable container for genomic intervals. It subclasses Python's built-in `list`, so all standard list operations work, plus it adds genomic interval logic on top.

---

## Construction

```python
from genomeblocks import Loci, Locus

# Empty
loci = Loci()

# From a list
loci = Loci([Locus("chr1", 1000, 2000), Locus("chr1", 5000, 6000)])

# From a file (auto-detects .bed, .narrowPeak, etc.)
loci = Loci.make("peaks.narrowPeak")

# Tile a genome into fixed-size windows
tiles = Loci.tile_genome("hg38.chrom.sizes", size=10000)

# Tile a single chromosome
tiles = Loci.tile("chr1", size=5000, chromsizes={"chr1": 248956422})
```

---

## Locus

Each element is a `Locus`:

```python
@dataclass
class Locus:
    chrom: str
    start: int   # 0-based
    end: int
    strand: str = "."
```

Key properties:
- `uid` — unique string ID: `"chr1:1000-2000(+)"`
- `length` — `end - start`
- `center` — `(start + end) // 2`

---

## Interval Operations

### Overlap query

```python
hits = loci.overlaps("chr1", 1000, 5000)
hits = loci.overlaps(some_locus)
```

Uses a [cgranges](https://github.com/lh3/cgranges) index (built lazily on first access) for fast lookups.

### Set operations

```python
common   = loci_a & loci_b    # intersection
only_a   = loci_a - loci_b    # difference (not in b)
xor      = loci_a ^ loci_b    # symmetric difference
combined = loci_a + loci_b    # concatenation
```

### Sort and merge

```python
merged = loci.sort().merge()  # canonical pre-processing step
```

### Slop (expand intervals)

```python
expanded = loci.slop(2500)  # +2500 bp on each side
```

---

## Mapping

`map` returns a dict of which loci in `other` overlap each locus in `self`:

```python
# {promoter.uid: [overlapping_enhancer.uid, ...]}
mapping = promoters.map(enhancers)
```

`tag` does the same but returns a `Tags` object for downstream labeling.

---

## Nearest

```python
df = loci.nearest(other_loci)  # pandas DataFrame with distances
```

---

## Export

```python
df  = loci.to_frame()     # pandas DataFrame (Chr, Start, End, Strand, Name)
pr  = loci.to_pyranges()  # pyranges PyRanges
loci.to_bed("out.bed")    # write BED file
```

---

## Signal Extraction

```python
# Returns numpy array of shape (n_loci, n_tracks, n_bins)
cube = loci.signal(["track1.bw", "track2.bw"], n_bins=200, flank=3000, workers=4)
```

---

## Indexing

```python
loci[0]                    # by position → Locus
loci["chr1:1000-2000(+)"]  # by UID string → Locus
loci[0:10]                 # slice → Loci
```

---

## Design Notes

- **Lazy indices**: the cgranges overlap index (`cgr`) and UID map (`uids`) are built only when first accessed — cheap to construct, fast to query.
- **Pickling-safe**: `cgr` is excluded from `__getstate__` and rebuilt on demand after deserialization.
- **Dynamic methods**: `signal`, `scan_motifs`, `enrich`, `pair_to_bed`, and plotting helpers are attached at import time from submodules, keeping the core class lightweight.
