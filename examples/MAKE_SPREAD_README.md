# Architecture.make_spread() - New Feature

## Quick Start

```python
from genomeblocks import Loci, Architecture, Genes

# Load your data
all_cres = Loci.from_bed("cres.bed")
genes = Genes.from_gencode("gencode.gtf")
promoters = genes.promoters(upstream=2000, downstream=500)

# Build spreading network
G = Architecture.make_spread(
    loci=all_cres,           # All possible loci
    source_loci=promoters,   # Starting points
    bedpe="loops.bedpe",     # Chromatin loops
    hops=2,                  # Number of spreading iterations
    directed=True            # Directed edges (source → target)
)

print(f"Network: {G.n_loci} loci, {G.n_links} links")
```

## What Does It Do?

Creates an architecture graph by **spreading** from source loci (like promoters) through chromatin loops:

1. **Hop 1**: Find all loci connected to sources via loops
2. **Hop 2**: Find all loci connected to Hop 1 targets
3. **Hop N**: Continue spreading...

Key feature: Uses the **center of the second anchor** for precise mapping.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `loci` | Loci | required | Complete loci collection |
| `source_loci` | Loci | required | Starting loci (e.g., promoters) |
| `bedpe` | str | required | Path to BEDPE loops file |
| `name` | str | "Spread" | Graph name |
| `r` | int | 2500 | Search radius around anchor centers (bp) |
| `dmax` | float | 1e9 | Max loop distance (bp) |
| `hops` | int | 1 | Number of spreading iterations |
| `directed` | bool | False | Create directed graph |
| `verbose` | bool | True | Print progress |

## Examples

See:
- `examples/make_spread_example.py` - Comprehensive usage examples
- `examples/test_make_spread.py` - Unit tests
- `examples/MAKE_SPREAD_GUIDE.md` - Detailed documentation

## Use Cases

- **Find enhancers** connected to your genes of interest
- **Explore regulatory neighborhoods** around transcription factors
- **Build directed regulatory networks** from promoters
- **Multi-hop exploration** of chromatin contact networks
- **TAD-scale structure** analysis

## Comparison

| Method | Description |
|--------|-------------|
| `Architecture.make()` | All loops between any loci (full network) |
| `Architecture.make_clique()` | Fully connected (all possible pairs) |
| `Architecture.make_spread()` | **Spread from sources through loops** ⭐ |

## Integration

Works seamlessly with existing Architecture methods:

```python
# Build spreading network
G = Architecture.make_spread(loci, promoters, "loops.bedpe", hops=2)

# Add Hi-C weights
G.add_mcool(loci, "hic.mcool", resolution=5000)

# Normalize by distance
G.normalize(loci, source="w", name="n")

# Annotate with genes
G.annotate(loci, genes)

# Draw subgraph
G.draw(loci, region=("chr1", 1000000, 2000000))
```

## Implementation Details

- Uses `set_directed()` to support both directed and undirected graphs
- Tracks active loci at each hop level
- Prevents backtracking to source loci in multi-hop exploration
- Efficiently handles large BEDPE files with progress bars
- Compatible with existing Architecture serialization (`__getstate__`/`__setstate__`)
