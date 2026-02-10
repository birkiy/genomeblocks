# Architecture.make_spread() - Spreading Network from Source Loci

## Overview

`Architecture.make_spread()` builds a chromatin architecture graph by spreading outward from a set of source loci (e.g., promoters) through chromatin loops. This is perfect for exploring regulatory networks where you want to find all enhancers and other regulatory elements connected to your genes of interest.

## Key Concept: Center of Second Anchor

For each loop in the BEDPE file:
- If anchor1 overlaps a source locus → use the **center of anchor2** to find target loci
- If anchor2 overlaps a source locus → use the **center of anchor1** to find target loci

This ensures precise mapping while being flexible about exact loop boundaries.

---

## Visual Example: Multi-Hop Spreading

```
Source Loci (Promoters):
    P1 (chr1:1000-2000)

BEDPE Loops:
    chr1  1000  2000  chr1   5000   6000   # Loop 1: P1 ↔ E1
    chr1  5000  6000  chr1  10000  11000   # Loop 2: E1 ↔ E2
    chr1 10000 11000  chr1  15000  16000   # Loop 3: E2 ↔ E3
    chr1  1000  2000  chr1  20000  21000   # Loop 4: P1 ↔ E4

All CREs (Loci):
    L1 = chr1:1000-2000   (promoter P1)
    L2 = chr1:5000-6000   (enhancer E1)
    L3 = chr1:10000-11000 (enhancer E2)
    L4 = chr1:15000-16000 (enhancer E3)
    L5 = chr1:20000-21000 (enhancer E4)
```

### Hop 1: Direct Connections

```
Start: P1 (L1)

Process Loop 1: P1 (anchor1) ↔ center(5500) → finds E1 (L2)
  → Add edge: P1 → E1

Process Loop 4: P1 (anchor1) ↔ center(20500) → finds E4 (L5)
  → Add edge: P1 → E4

Result:
    P1 → E1
    P1 → E4

Graph: {P1, E1, E4}, 2 edges
New loci for next hop: {E1, E4}
```

### Hop 2: Connections of Connections

```
Active loci from Hop 1: {E1, E4}

Process Loop 2: E1 (anchor1) ↔ center(10500) → finds E2 (L3)
  → Add edge: E1 → E2

Result:
    P1 → E1
    P1 → E4
    E1 → E2

Graph: {P1, E1, E4, E2}, 3 edges
New loci for next hop: {E2}
```

### Hop 3: Third-Order Connections

```
Active loci from Hop 2: {E2}

Process Loop 3: E2 (anchor1) ↔ center(15500) → finds E3 (L4)
  → Add edge: E2 → E3

Final Result:
    P1 → E1
    P1 → E4
    E1 → E2
    E2 → E3

Graph: {P1, E1, E4, E2, E3}, 4 edges
```

---

## Directed vs Undirected

### Directed Graph (`directed=True`)
Edges point **away** from source loci, showing flow of information/regulation:
```
P1 → E1 → E2 → E3
P1 → E4
```

Use cases:
- Understanding regulatory cascades from promoters
- Signal flow analysis
- Identifying downstream targets

### Undirected Graph (`directed=False`)
Edges are bidirectional, showing symmetric connectivity:
```
P1 ↔ E1 ↔ E2 ↔ E3
P1 ↔ E4
```

Use cases:
- General network analysis
- Community detection
- Clustering related elements

---

## Parameter Guide

### Required Parameters

- **`loci`**: Complete Loci collection (all possible CREs/regulatory elements)
- **`source_loci`**: Starting Loci (e.g., promoters, specific enhancers)
- **`bedpe`**: Path to BEDPE file with chromatin loops

### Key Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `hops` | 1 | Number of spreading iterations (1=direct, 2=second-order, etc.) |
| `directed` | False | Create directed (True) or undirected (False) graph |
| `r` | 2500 | Radius (bp) around loop anchor centers to search for loci |
| `dmax` | 1e9 | Maximum genomic distance between loop anchors |
| `name` | "Spread" | Name for the Architecture graph |
| `verbose` | True | Print progress information |

---

## Common Usage Patterns

### Pattern 1: Find Direct Enhancer Targets of Genes

```python
# Get promoters for your genes of interest
gene_promoters = genes.filter(lambda g: g.name in ["MYC", "TP53", "BRCA1"])
                     .promoters(upstream=2000, downstream=500)

# Build 1-hop network
G = Architecture.make_spread(
    loci=all_cres,
    source_loci=gene_promoters,
    bedpe="H1_HiChIP_loops.bedpe",
    hops=1,
    directed=True
)
```

### Pattern 2: Explore Regulatory Neighborhoods

```python
# Start from active enhancers
active_enhancers = all_cres.filter(lambda l: l.H3K27ac > 10)

# Explore 2-hop neighborhood
G = Architecture.make_spread(
    loci=all_cres,
    source_loci=active_enhancers,
    bedpe="loops.bedpe",
    hops=2,
    directed=False,
    dmax=1e6  # Focus on local interactions
)
```

### Pattern 3: TAD-Scale Exploration

```python
# From promoters, explore large-scale structure
G = Architecture.make_spread(
    loci=all_cres,
    source_loci=promoters,
    bedpe="HiC_loops.bedpe",
    hops=3,
    directed=False,
    dmax=2e6  # Allow long-range loops
)
```

---

## Combining with Other Methods

### Add Contact Frequencies

```python
G = Architecture.make_spread(
    loci=all_cres,
    source_loci=promoters,
    bedpe="loops.bedpe",
    hops=2
)

# Add Hi-C contact counts
G.add_mcool(loci=all_cres, mcool="hic.mcool", resolution=5000, name="w")

# Normalize by distance
G.normalize(loci=all_cres, source="w", name="n")
```

### Annotate with Genes

```python
G = Architecture.make_spread(
    loci=all_cres,
    source_loci=promoters,
    bedpe="loops.bedpe",
    hops=1
)

# Add gene annotations
G.annotate(loci=all_cres, genes=genes)

# Access annotations
for v in G.vertices():
    uid = G.vp.uid[v]
    gene = G.vp.gene[v]
    annot = G.vp.annot[v]
    print(f"{uid}: {gene} ({annot})")
```

---

## Performance Tips

1. **Start with fewer hops**: More hops = exponentially larger graphs
   - Hop 1: Typically hundreds of loci
   - Hop 2: Typically thousands of loci
   - Hop 3: Can be tens of thousands

2. **Use `dmax` to filter**: Restricting loop distance reduces computation
   ```python
   dmax=1e6  # Only loops < 1 Mb (local interactions)
   ```

3. **Adjust `r` for precision**: Smaller `r` = more precise but may miss some loops
   ```python
   r=1000   # Strict: ±1kb
   r=2500   # Default: ±2.5kb
   r=5000   # Permissive: ±5kb
   ```

4. **Filter source loci**: Start with a focused set
   ```python
   # Instead of all promoters (20,000)
   gene_promoters = genes.filter(lambda g: g.tpm > 5).promoters()
   ```

---

## Comparison with Other Methods

| Method | Use Case |
|--------|----------|
| `make()` | All loops between any loci (full network) |
| `make_clique()` | Fully connected graph (all possible edges) |
| `make_spread()` | **Spreading from specific sources through loops** |

`make_spread()` is best when you have:
- Specific starting points (promoters, super-enhancers, etc.)
- Want to explore their regulatory neighborhoods
- Need to understand multi-order connections
- Want directed graphs showing regulatory flow

---

## Example Output

```
[INFO] Starting with 1247 source loci
[INFO] Hop 1/2: Scanning loops: 100%|████████████| 125847/125847
[INFO] Hop 1: 3521 loops mapped, 2847 new loci discovered
[INFO] Hop 2/2: Scanning loops: 100%|████████████| 125847/125847
[INFO] Hop 2: 8934 loops mapped, 5124 new loci discovered
[INFO] Complete: 125847 loops scanned | 12455 mapped (9.9%) | loci=9218, links=12455
[INFO] Graph is directed (→)
```

This shows:
- Started with 1,247 promoters
- After 2 hops: 9,218 loci total, 12,455 edges
- Only 9.9% of all loops were relevant (good filtering!)
