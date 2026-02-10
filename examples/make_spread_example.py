"""
Example: Using Architecture.make_spread() to explore chromatin loops from promoters

This example demonstrates how to:
1. Define source loci (e.g., promoters) 
2. Build a spreading network from loops (NO pre-defined CRE set needed!)
3. Explore with different hop levels and directionality
4. The network discovers anchors as it spreads through the BEDPE
"""

from genomeblocks import Loci, Architecture, Genes

# ============================================================================
# Setup: Load your data
# ============================================================================

# Load gene annotations to get promoters
genes = Genes.from_gencode("path/to/gencode.gtf")

# Define source loci - for example, promoters around TSS
promoters = genes.promoters(upstream=2000, downstream=500)

# BEDPE file with chromatin loops (e.g., from HiChIP, ChIA-PET, Hi-C)
loops_file = "path/to/loops.bedpe"

# NOTE: You DON'T need a pre-defined set of all CREs!
# The network is built directly from the BEDPE anchors


# ============================================================================
# Example 1: Direct connections only (1 hop, undirected)
# ============================================================================

print("\n=== Example 1: Direct promoter connections (1 hop, undirected) ===")
G1 = Architecture.make_spread(
    source_loci=promoters,
    bedpe=loops_file,
    name="Promoter_Direct",
    r=2500,          # Search ±2.5kb around loop anchor centers
    dmax=1e6,        # Only consider loops < 1Mb
    hops=1,          # Only direct connections
    directed=False,  # Undirected graph
    verbose=True
)

print(f"Result: {G1.n_loci} anchors, {G1.n_links} links")
# Vertices are loop anchor centers in format "chr:pos"


# ============================================================================
# Example 2: 2-hop exploration (directed)
# ============================================================================

print("\n=== Example 2: 2-hop exploration (directed) ===")
# Promoter → Anchor1 → Anchor2
# This finds:
#   - Hop 1: All anchors directly connected to promoters via loops
#   - Hop 2: All anchors connected to those first-level anchors

G2 = Architecture.make_spread(
    source_loci=promoters,
    bedpe=loops_file,
    name="Promoter_2Hop_Directed",
    r=2500,
    dmax=1e6,
    hops=2,          # 2 hops
    directed=True,   # Directed graph (edges point away from promoters)
    verbose=True
)

print(f"Result: {G2.n_loci} anchors, {G2.n_links} links")


# ============================================================================
# Example 3: 3-hop exploration for long-range spreading
# ============================================================================

print("\n=== Example 3: 3-hop exploration (undirected) ===")
G3 = Architecture.make_spread(
    source_loci=promoters,
    bedpe=loops_file,
    name="Promoter_3Hop",
    r=2500,
    dmax=2e6,        # Allow loops up to 2Mb
    hops=3,          # 3 hops
    directed=False,
    verbose=True
)

print(f"Result: {G3.n_loci} anchors, {G3.n_links} links")


# ============================================================================
# Example 4: Using a specific gene's promoter
# ============================================================================

print("\n=== Example 4: Spread from MYC promoter only ===")

# Get just MYC gene
myc_promoters = genes.filter(lambda g: g.name == "MYC").promoters(upstream=2000, downstream=500)

G4 = Architecture.make_spread(
    source_loci=myc_promoters,
    bedpe=loops_file,
    name="MYC_Network",
    r=2500,
    hops=2,
    directed=True,
    verbose=True
)

print(f"MYC network: {G4.n_loci} anchors, {G4.n_links} links")


# ============================================================================
# Example 5: Building complete network then filtering
# ============================================================================

print("\n=== Example 5: Alternative approach - Build full network ===")
# If you want to work with a full loci set, you can still use the old approach:
# 1. Build a complete architecture from bedpe
# 2. Filter to promoter neighborhoods

all_cres = Loci.from_bed("path/to/all_cres.bed")  # IF you have this

# Option A: Use make() to build full network
G5a = Architecture.make(
    loci=all_cres,
    bedpe=loops_file,
    name="Full_Network",
    verbose=True
)

# Then filter or analyze promoter neighborhoods
print(f"Full network: {G5a.n_loci} loci, {G5a.n_links} links")

# Option B: Use make_spread() without pre-defined loci
G5b = Architecture.make_spread(
    source_loci=promoters,
    bedpe=loops_file,
    name="Spread_Network",
    hops=1,
    verbose=True
)

print(f"Spread network: {G5b.n_loci} anchors, {G5b.n_links} links")
print("Note: Anchors are loop centers (chr:pos format)")


# ============================================================================
# Tips for interpretation
# ============================================================================

"""
Understanding the spreading network:

NEW BEHAVIOR:
- Vertices are LOOP ANCHOR CENTERS (format: "chr:pos")
- No need for pre-defined CRE collection
- Network grows organically from BEDPE

1. **Hops = 1**: Only direct loop connections
   - Source A connects to all X where A-X forms a loop
   
2. **Hops = 2**: Second-order connections
   - Source A -> X1, X2, ... (first hop)
   - X1 -> Y1, Y2, ... (second hop)
   - X2 -> Z1, Z2, ... (second hop)
   - Result includes A, all X's, and all Y's and Z's
   
3. **Hops = 3**: Third-order connections
   - Continues spreading from second-order targets
   
4. **Directed vs Undirected**:
   - directed=True: Edges point away from source (A → X)
     - Useful for understanding signal flow from promoters
   - directed=False: Bidirectional edges (A ↔ X)
     - Useful for general network analysis
     
5. **The r parameter**:
   - Controls how precisely loops must match loci
   - r=2500 means ±2.5kb around loop anchor center
   - Larger r = more permissive matching
   
6. **The dmax parameter**:
   - Filters loops by genomic distance
   - dmax=1e6 = only loops < 1Mb
   - Helps focus on local vs long-range interactions
"""
