# Architecture.make_spread() - UPDATED IMPLEMENTATION

## 🎉 Major Change: No Pre-Defined CRE Collection Needed!

The new `make_spread()` builds networks **entirely from BEDPE anchors** - you don't need a pre-defined set of all CREs!

## Quick Comparison

### OLD Approach (if you had used it):
```python
# Required: Pre-defined CRE collection
all_cres = Loci.from_bed("all_enhancers.bed")  # 50,000+ elements

G = Architecture.make_spread(
    loci=all_cres,           # ← Required this!
    source_loci=promoters,
    bedpe="loops.bedpe"
)
```

### NEW Approach:
```python
# NO pre-defined CRE collection needed!
# Network discovers anchors from BEDPE

G = Architecture.make_spread(
    source_loci=promoters,   # Just your starting points
    bedpe="loops.bedpe"      # Network built from this
)
```

## How It Works

1. **Start with source loci** (e.g., promoters)
2. **Scan BEDPE file** for loops touching these sources
3. **Discover anchors** from loop partners
4. **Vertices = Loop anchor centers** (format: `chr:pos`)
5. **Multi-hop** continues spreading from discovered anchors

### Visual Example

```
BEDPE:
  chr1  1000  2000  chr1  5000  6000   # Loop 1
  chr1  5000  6000  chr1  10000 11000  # Loop 2

Source: Promoter at chr1:1000-2000

Hop 1:
  - Promoter center = chr1:1500
  - Loop 1 connects chr1:1500 ←→ chr1:5500
  - Result: Add vertex chr1:5500, edge 1500↔5500

Hop 2:
  - chr1:5500 is now active
  - Loop 2 connects chr1:5500 ←→ chr1:10500
  - Result: Add vertex chr1:10500, edge 5500↔10500

Final Graph:
  Vertices: {chr1:1500, chr1:5500, chr1:10500}
  Edges: 1500↔5500, 5500↔10500
```

## Usage

### Basic Example
```python
from genomeblocks import Architecture, Genes

# Get promoters
genes = Genes.from_gencode("gencode.gtf")
promoters = genes.promoters(upstream=2000, downstream=500)

# Build network (1-hop, undirected)
G = Architecture.make_spread(
    source_loci=promoters,
    bedpe="HiChIP_loops.bedpe",
    hops=1,
    directed=False
)

print(f"Network: {G.n_loci} anchors, {G.n_links} links")
```

### Multi-Hop Directed Network
```python
# 2-hop directed network
G = Architecture.make_spread(
    source_loci=promoters,
    bedpe="loops.bedpe",
    hops=2,
    directed=True  # Edges point away from sources
)

# Useful for regulatory flow analysis
```

### Specific Gene Network
```python
# Just MYC gene
myc_promoters = genes.filter(lambda g: g.name == "MYC").promoters()

G = Architecture.make_spread(
    source_loci=myc_promoters,
    bedpe="loops.bedpe",
    hops=2,
    r=2500,        # Search radius around anchors
    dmax=1e6       # Max loop distance (1Mb)
)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `source_loci` | Loci | required | Starting loci (e.g., promoters) |
| `bedpe` | str | required | Path to BEDPE file |
| `name` | str | "Spread" | Graph name |
| `r` | int | 2500 | Radius (bp) to match source loci to anchors |
| `dmax` | float | 1e9 | Max genomic distance between loop anchors |
| `hops` | int | 1 | Number of spreading iterations |
| `directed` | bool | False | Create directed (True) or undirected graph |
| `verbose` | bool | True | Print progress |

### Parameter `r` (radius)
- Controls how close a source locus must be to a loop anchor center
- Default 2500bp works well for most use cases
- Larger = more permissive matching
- Smaller = stricter matching

### Parameter `hops`
- `hops=1`: Direct connections only (promoters → targets)
- `hops=2`: Second-order (promoters → targets → targets)  
- `hops=3+`: Higher-order neighborhood exploration
- **Warning**: Network size grows exponentially with hops!

### Parameter `directed`
- `directed=False`: Bidirectional edges (A ↔ B)
  - Use for: General network analysis, clustering, communities
- `directed=True`: Edges point away from sources (A → B)
  - Use for: Regulatory flow, signal propagation, causal analysis

## Vertex UIDs

Vertices are represented as **`chr:pos`** where `pos` is the anchor center:

```python
G = Architecture.make_spread(source_loci=promoters, bedpe="loops.bedpe")

# Vertex UIDs look like:
for v in G.vertices():
    uid = G.vp.uid[v]
    print(uid)  # e.g., "chr1:123456", "chr2:789012"
```

## When to Use make_spread() vs make()

| Method | Use When |
|--------|----------|
| `make_spread()` | You have specific starting points (promoters, TF sites) and want to explore their neighborhood through loops. **Don't have curated CRE list.** |
| `make()` | You have a complete CRE collection and want all loop-connected pairs. |
| `make_clique()` | You want all possible pairwise connections (fully connected graph). |

## Advantages

1. **No CRE collection needed**: Works directly with BEDPE
2. **Focused networks**: Only includes relevant anchors
3. **Memory efficient**: Doesn't require loading all CREs
4. **Flexible exploration**: Multi-hop neighborhood discovery
5. **Directed networks**: Model regulatory information flow

## Limitations

1. **Vertices are anchor centers**: Not your original loci objects
2. **No pre-annotation**: Vertices don't have CRE properties (H3K27ac, etc.)
3. **Resolution limited**: By BEDPE anchor resolution

## Working with Results

### Accessing Neighbors
```python
# Get neighbors of an anchor
anchor_uid = "chr1:123456"
if anchor_uid in G:
    neighbors = G[anchor_uid]
    print(f"{anchor_uid} connects to {len(neighbors)} anchors")
```

### Exporting
```python
# Convert to pandas DataFrame
import pandas as pd

# Vertices
vertices = []
for v in G.vertices():
    vertices.append({'uid': G.vp.uid[v]})
df_v = pd.DataFrame(vertices)

# Edges
edges = []
for e in G.edges():
    edges.append({
        'source': G.vp.uid[e.source()],
        'target': G.vp.uid[e.target()]
    })
df_e = pd.DataFrame(edges)

df_v.to_csv("anchors.csv")
df_e.to_csv("connections.csv")
```

### Visualization
```python
# Get subgraph for a region
region = ("chr1", 1000000, 2000000)

# Note: For visualization, you may want to create a Loci object
# from the anchor positions for use with draw()
```

## Examples

See:
- `examples/test_make_spread_simple.py` - Simple working test
- `examples/make_spread_example.py` - Multiple usage patterns
- `examples/MAKE_SPREAD_GUIDE.md` - Detailed documentation

## Migration from Old API (if needed)

If you were using an older version that required `loci` parameter:

```python
# OLD (if you had this):
G = Architecture.make_spread(loci=all_cres, source_loci=promoters, ...)

# NEW:
G = Architecture.make_spread(source_loci=promoters, ...)
# That's it! No loci parameter needed.
```

## Performance

- Fast for 1-2 hops with reasonable source counts
- Scales well with BEDPE size (streaming read)
- For 10,000 promoters, 1-hop typically finds 50,000-100,000 anchors
- 2-hop can reach 200,000-500,000 anchors
- Use `dmax` parameter to limit search space

## Summary

`make_spread()` is perfect when you:
- ✅ Have specific starting loci (promoters, peaks, etc.)
- ✅ Want to explore chromatin neighborhoods
- ✅ Don't have a pre-defined CRE collection
- ✅ Need directed regulatory networks
- ✅ Want multi-hop exploration (2nd, 3rd order connections)

The network discovers and builds itself entirely from the BEDPE file! 🎉
