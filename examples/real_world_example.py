"""
Real-World Example: Finding Regulatory Networks of Specific Genes

This example shows how to use Architecture.make_spread() to:
1. Start from specific gene promoters (e.g., disease-associated genes)
2. Find all connected regulatory elements
3. Analyze the resulting network
4. Visualize and export results
"""

from genomeblocks import Loci, Architecture, Genes

# ============================================================================
# Scenario: Find the regulatory network of cancer-associated genes
# ============================================================================

print("="*70)
print("REAL-WORLD EXAMPLE: Cancer Gene Regulatory Networks")
print("="*70)

# Load data
all_cres = Loci.from_bed("data/all_enhancers_and_promoters.bed")
genes = Genes.from_gencode("data/gencode.v38.annotation.gtf")
loops_file = "data/K562_HiChIP_loops.bedpe"  # Cell type specific

# Define genes of interest (example: common cancer genes)
cancer_genes = ["MYC", "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA"]

print(f"\n1. Starting from {len(cancer_genes)} cancer-associated genes")
print(f"   Genes: {', '.join(cancer_genes)}")

# Get promoters for these genes
cancer_gene_objs = genes.filter(lambda g: g.name in cancer_genes)
cancer_promoters = cancer_gene_objs.promoters(upstream=2000, downstream=500)

print(f"   Found {len(cancer_promoters)} promoter regions")

# ============================================================================
# Build multi-hop spreading network
# ============================================================================

print("\n2. Building 2-hop regulatory network...")

G = Architecture.make_spread(
    loci=all_cres,
    source_loci=cancer_promoters,
    bedpe=loops_file,
    name="CancerGene_Network",
    r=2500,          # Standard radius
    dmax=1e6,        # Focus on loops < 1Mb
    hops=2,          # 2-hop: promoters → enhancers → enhancers
    directed=True,   # Directed: shows regulatory flow
    verbose=True
)

print(f"\n   Result: {G.n_loci} regulatory elements, {G.n_links} connections")

# ============================================================================
# Add contact frequencies from Hi-C
# ============================================================================

print("\n3. Adding Hi-C contact frequencies...")

G.add_mcool(
    loci=all_cres,
    mcool="data/K562_hic.mcool",
    resolution=5000,
    name="w",
    verbose=True
)

# Normalize by genomic distance to get O/E
G.normalize(
    loci=all_cres,
    source="w",
    name="oe",  # Observed/Expected
    verbose=True
)

# ============================================================================
# Annotate with gene information
# ============================================================================

print("\n4. Annotating with nearest genes and genomic features...")

G.annotate(loci=all_cres, genes=genes, verbose=True)

# ============================================================================
# Analyze the network
# ============================================================================

print("\n5. Network Analysis")
print("="*70)

# Count by hop distance from sources
print("\n   a) Analyzing network structure:")

# Get all edges and their properties
edges_data = []
for e in G.edges():
    v1, v2 = e.source(), e.target()
    uid1 = G.vp.uid[v1]
    uid2 = G.vp.uid[v2]
    
    # Check if source is a cancer promoter
    is_source = uid1 in {p.uid for p in cancer_promoters}
    
    edges_data.append({
        'source': uid1,
        'target': uid2,
        'is_direct': is_source,
        'weight': G.ep.w[e],
        'oe': G.ep.n[e] if 'n' in G.ep else None,
        'distance': G.ep.d[e] if 'd' in G.ep else None
    })

import pandas as pd
edges_df = pd.DataFrame(edges_data)

n_direct = edges_df['is_direct'].sum()
n_indirect = len(edges_df) - n_direct

print(f"      - Direct connections (hop 1): {n_direct}")
print(f"      - Indirect connections (hop 2): {n_indirect}")
print(f"      - Total connections: {len(edges_df)}")

# ============================================================================
# Find highly connected enhancers (hubs)
# ============================================================================

print("\n   b) Finding regulatory hubs:")

# Count connections per locus
from collections import Counter
target_counts = Counter(edges_df['target'])

# Get top 10 most connected
top_hubs = target_counts.most_common(10)

print("\n      Top 10 most connected enhancers:")
for i, (uid, count) in enumerate(top_hubs, 1):
    v = G.index[uid]
    gene = G.vp.gene[v] if 'gene' in G.vp else "N/A"
    annot = G.vp.annot[v] if 'annot' in G.vp else "N/A"
    print(f"      {i:2d}. {uid}: {count} connections | near {gene} | {annot}")

# ============================================================================
# Identify high-confidence connections (strong Hi-C signal)
# ============================================================================

print("\n   c) High-confidence regulatory interactions (O/E > 2.0):")

if 'oe' in edges_df.columns and edges_df['oe'].notna().any():
    high_conf = edges_df[edges_df['is_direct'] & (edges_df['oe'] > 2.0)]
    high_conf = high_conf.sort_values('oe', ascending=False).head(10)
    
    print(f"\n      Found {len(high_conf)} direct connections with O/E > 2.0")
    print("      Top 10:")
    for i, row in enumerate(high_conf.itertuples(), 1):
        src_v = G.index[row.source]
        tgt_v = G.index[row.target]
        src_gene = G.vp.gene[src_v] if 'gene' in G.vp else "N/A"
        tgt_gene = G.vp.gene[tgt_v] if 'gene' in G.vp else "N/A"
        print(f"      {i:2d}. {src_gene} → {tgt_gene} (O/E={row.oe:.2f}, dist={row.distance/1000:.1f}kb)")

# ============================================================================
# Export results
# ============================================================================

print("\n6. Exporting Results")
print("="*70)

# Export edge list
edges_df.to_csv("cancer_gene_network_edges.csv", index=False)
print("   ✓ Saved edge list to: cancer_gene_network_edges.csv")

# Export vertex annotations
vertices_data = []
for v in G.vertices():
    uid = G.vp.uid[v]
    vertices_data.append({
        'uid': uid,
        'gene': G.vp.gene[v] if 'gene' in G.vp else "",
        'annotation': G.vp.annot[v] if 'annot' in G.vp else "",
        'is_source': uid in {p.uid for p in cancer_promoters}
    })

vertices_df = pd.DataFrame(vertices_data)
vertices_df.to_csv("cancer_gene_network_vertices.csv", index=False)
print("   ✓ Saved vertex list to: cancer_gene_network_vertices.csv")

# ============================================================================
# Visualize a specific region
# ============================================================================

print("\n7. Visualization Example")
print("="*70)

# Let's visualize the network around one of our genes, e.g., MYC
myc_promoter = cancer_promoters.filter(lambda p: 'MYC' in str(p))

if len(myc_promoter) > 0:
    myc_loc = myc_promoter[0]
    print(f"   Visualizing MYC locus: {myc_loc}")
    
    # Expand region to show context
    region = (myc_loc.chrom, myc_loc.start - 500000, myc_loc.end + 500000)
    
    print(f"   Region: {region[0]}:{region[1]:,}-{region[2]:,}")
    
    # Draw the network
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    G.draw(
        loci=all_cres,
        region=region,
        vertex_size_by='Agg_H1',  # Size by signal (if available in loci)
        edge_width_by='oe',        # Width by O/E
        vertex_color='annot',      # Color by genomic annotation
        edge_color='#CCCCCC',
        layout='spring',
        show_labels=True,
        label_prop='gene',
        ax=ax
    )
    
    plt.savefig('myc_regulatory_network.pdf', bbox_inches='tight', dpi=300)
    print("   ✓ Saved figure to: myc_regulatory_network.pdf")

# ============================================================================
# Summary Statistics
# ============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
Input:
  - {len(cancer_genes)} cancer genes
  - {len(cancer_promoters)} promoter regions
  - {len(all_cres)} total regulatory elements
  
Network:
  - {G.n_loci} nodes (regulatory elements)
  - {G.n_links} edges (regulatory connections)
  - {n_direct} direct connections (hop 1)
  - {n_indirect} indirect connections (hop 2)
  
Characterization:
  - Graph type: {'Directed' if G.is_directed() else 'Undirected'}
  - Enriched for: Hi-C contact frequency (O/E)
  - Annotated with: Nearest genes, genomic features
  
Outputs:
  ✓ Edge list: cancer_gene_network_edges.csv
  ✓ Vertex list: cancer_gene_network_vertices.csv
  ✓ Visualization: myc_regulatory_network.pdf
""")

print("="*70)
print("Analysis Complete! 🎉")
print("="*70)

# ============================================================================
# Next steps and interpretation
# ============================================================================

print("""
NEXT STEPS:

1. Functional enrichment:
   - Analyze which pathways/functions are enriched in the network
   - Look for GO term enrichment in connected genes
   
2. Compare cell types:
   - Run the same analysis with loops from different cell types
   - Identify cell-type-specific regulatory connections
   
3. Prioritize variants:
   - Overlay GWAS variants on the network
   - Identify variants in regulatory hubs
   
4. Predict gene expression:
   - Use network structure to predict gene expression changes
   - Model enhancer-promoter cooperativity
   
5. Network modules:
   - Find communities/modules in the network
   - Identify co-regulated gene sets

For more information, see:
  - examples/MAKE_SPREAD_GUIDE.md
  - examples/make_spread_example.py
""")
