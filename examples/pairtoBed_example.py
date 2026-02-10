"""Example demonstrating pairToBed-like functionality in genomeblocks."""

from genomeblocks import Loci, Locus
from genomeblocks.bedpe import pair_to_bed, read_bedpe, pairs_to_frame, Pair

# Example 1: Create some test loci
loci = Loci([
    Locus("chr1", 1000, 2000, "+"),
    Locus("chr1", 5000, 6000, "+"),
    Locus("chr2", 10000, 11000, "-"),
])

print("=" * 60)
print("Example 1: Using pair_to_bed with Loci object")
print("=" * 60)

# Example 2: Create some test pairs (simulating a BEDPE file)
test_pairs = [
    Pair("chr1", 1500, 1700, "chr1", 5500, 5700, name="pair1", score=100),
    Pair("chr1", 3000, 3200, "chr1", 7000, 7200, name="pair2", score=50),
    Pair("chr1", 1800, 2100, "chr2", 10500, 10700, name="pair3", score=75),
    Pair("chr2", 10200, 10400, "chr2", 15000, 15200, name="pair4", score=25),
]

print("")
print("Test pairs:")
for p in test_pairs:
    print("  " + str(p))

# Find pairs where either anchor overlaps our loci
print("")
print("--- Finding pairs with 'either' mode (at least one anchor overlaps) ---")
overlapping_either = pair_to_bed(loci, test_pairs, either=True, both=False, verbose=True)
print("")
print("Results:")
for p in overlapping_either:
    print("  " + p.name + ": " + str(p))

# Find pairs where both anchors overlap our loci
print("")
print("--- Finding pairs with 'both' mode (both anchors must overlap) ---")
overlapping_both = pair_to_bed(loci, test_pairs, either=False, both=True, verbose=True)
print("")
print("Results:")
for p in overlapping_both:
    print("  " + p.name + ": " + str(p))

# Example 3: Using with slop/extension
print("")
print("=" * 60)
print("Example 2: Using with slop (r=500bp extension)")
print("=" * 60)
overlapping_slop = pair_to_bed(loci, test_pairs, r=500, either=True, both=False, verbose=True)
print("")
print("Results with 500bp extension:")
for p in overlapping_slop:
    print("  " + p.name + ": " + str(p))

# Example 4: Using the Loci method directly
print("")
print("=" * 60)
print("Example 3: Using as a Loci method")
print("=" * 60)
overlapping = loci.pair_to_bed(test_pairs, verbose=True)
print("")
print("Results:")
for p in overlapping:
    print("  " + p.name + ": " + str(p))

# Example 5: Convert results to DataFrame
print("")
print("=" * 60)
print("Example 4: Convert results to DataFrame")
print("=" * 60)
df = pairs_to_frame(overlapping)
print(df)

# Example 6: Using with a single Locus
print("")
print("=" * 60)
print("Example 5: Using with a single Locus")
print("=" * 60)
single_locus = Locus("chr1", 1000, 2000, "+")
overlapping_single = pair_to_bed(single_locus, test_pairs, verbose=True)
print("")
print("Results:")
for p in overlapping_single:
    print("  " + p.name + ": " + str(p))

# Example 7: Filtering by score and distance
print("")
print("=" * 60)
print("Example 6: Filtering by score >= 50")
print("=" * 60)
overlapping_filtered = pair_to_bed(
    loci, 
    test_pairs, 
    min_score=50,
    verbose=True
)
print("")
print("Results:")
for p in overlapping_filtered:
    print("  " + p.name + ": score=" + str(p.score))

print("")
print("=" * 60)
print("Summary")
print("=" * 60)
print("""
The pair_to_bed function provides bedtools pairToBed-like functionality:

1. Find pairs where at least one anchor overlaps loci (either=True)
2. Find pairs where both anchors overlap loci (both=True)
3. Use slop/extension with the 'r' parameter
4. Filter by score and distance
5. Works with Loci objects, single Locus, or iterables
6. Can read from BEDPE files or use Pair objects directly
7. Convert results to DataFrame or write to BEDPE files

Usage:
    # From file
    pairs = loci.pair_to_bed('loops.bedpe', r=2500)
    
    # From Pair objects
    pairs = loci.pair_to_bed(my_pairs, both=True, min_score=50)
    
    # Standalone function
    from genomeblocks.bedpe import pair_to_bed
    pairs = pair_to_bed(loci, 'loops.bedpe', r=1000)
""")
