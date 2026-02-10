"""
Visual ASCII diagrams for Architecture.make_spread()

Run this file to see visual representations of how spreading works.
"""

def print_diagram_1hop():
    """Show 1-hop spreading"""
    print("\n" + "="*70)
    print("DIAGRAM 1: 1-HOP SPREADING (directed=True)")
    print("="*70)
    print("""
Source Loci:
    P1 = chr1:1000-2000 (promoter)

BEDPE Loops:
    chr1  1000  2000  chr1  5000  6000  # P1 ←→ E1
    chr1  1000  2000  chr1  8000  9000  # P1 ←→ E2

All Loci:
    P1 = chr1:1000-2000
    E1 = chr1:5000-6000
    E2 = chr1:8000-9000
    E3 = chr1:12000-13000

Step-by-step:

1. Start with source: P1
   
   Graph: [P1]

2. Process loop 1: P1 (anchor1) ←→ center(5500)
   - Find loci near center(5500) with radius ±2500 → finds E1
   - P1 is active → Add edge P1 → E1
   
   Graph: [P1] → [E1]

3. Process loop 2: P1 (anchor1) ←→ center(8500)
   - Find loci near center(8500) with radius ±2500 → finds E2
   - P1 is active → Add edge P1 → E2
   
   Graph: [P1] → [E1]
                 ↓
                [E2]

Final Network (1 hop):
    - Vertices: {P1, E1, E2}
    - Edges: P1→E1, P1→E2
    """)


def print_diagram_2hop():
    """Show 2-hop spreading"""
    print("\n" + "="*70)
    print("DIAGRAM 2: 2-HOP SPREADING (directed=True)")
    print("="*70)
    print("""
Source Loci:
    P1 = chr1:1000-2000 (promoter)

BEDPE Loops:
    chr1  1000  2000  chr1  5000  6000   # P1 ←→ E1
    chr1  5000  6000  chr1  10000  11000 # E1 ←→ E2
    chr1  1000  2000  chr1  15000  16000 # P1 ←→ E3

All Loci:
    P1 = chr1:1000-2000
    E1 = chr1:5000-6000
    E2 = chr1:10000-11000
    E3 = chr1:15000-16000
    E4 = chr1:20000-21000

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

HOP 1: Starting from P1

Active loci: {P1}

Process loops:
  ✓ Loop 1: P1 → center(5500) → E1
  ✓ Loop 3: P1 → center(15500) → E3

Discovered: {E1, E3}

Graph after Hop 1:
         
         [P1]
        /    \\
       ↓      ↓
      [E1]   [E3]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

HOP 2: Spreading from newly discovered loci

Active loci: {E1, E3}  (Note: P1 excluded to prevent backtracking)

Process loops:
  ✓ Loop 2: E1 → center(10500) → E2

Discovered: {E2}

Final Graph after Hop 2:

         [P1]
        /    \\
       ↓      ↓
      [E1]   [E3]
       ↓
      [E2]

Final Network:
  - Vertices: {P1, E1, E3, E2}
  - Edges: P1→E1, P1→E3, E1→E2
    """)


def print_diagram_undirected():
    """Show undirected spreading"""
    print("\n" + "="*70)
    print("DIAGRAM 3: UNDIRECTED vs DIRECTED")
    print("="*70)
    print("""
Same setup as Diagram 1:
    Source: P1
    Loops: P1←→E1, P1←→E2

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DIRECTED (directed=True):
    Edges point AWAY from source
    
    [P1] ──→ [E1]
      │
      └───→ [E2]
    
    Use case: Regulatory flow from promoters
    
    You can traverse: P1 → E1, P1 → E2
    You CANNOT traverse: E1 → P1 (wrong direction)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

UNDIRECTED (directed=False):
    Edges are bidirectional
    
    [P1] ←→ [E1]
      ↕
    [E2]
    
    Use case: General network analysis, clustering
    
    You can traverse in any direction:
      P1 ← E1 ← P1 → E2 → P1 (all valid)
    """)


def print_diagram_radius():
    """Show how radius parameter works"""
    print("\n" + "="*70)
    print("DIAGRAM 4: HOW PARAMETER 'r' (RADIUS) WORKS")
    print("="*70)
    print("""
BEDPE Loop:
    chr1  1000  2000  chr1  5000  6000
    
Anchor 1: chr1:1000-2000 → center = 1500
Anchor 2: chr1:5000-6000 → center = 5500

All Loci:
    L1 = chr1:1400-1600  (overlaps center ± r)
    L2 = chr1:1450-1550  (overlaps center ± r)
    L3 = chr1:4500-4600  (outside r, won't be found)
    L4 = chr1:5400-5600  (overlaps center ± r)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

WITH r=2500 (default):

Anchor 1 center = 1500
Search region: 1500 - 2500 = -1000 to 1500 + 2500 = 4000
              chr1:(-1000)-4000  (negative becomes 0)

    ←─────────── r=2500 ─────────→
              ↓ center
    |------------------------------|
   -1000                         4000
    
    Loci found: L1 ✓, L2 ✓, L3 ✗ (too far)

Anchor 2 center = 5500
Search region: 5500 - 2500 = 3000 to 5500 + 2500 = 8000
              chr1:3000-8000

    Loci found: L4 ✓

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

WITH r=500 (strict):

Anchor 1 center = 1500
Search region: 1500 - 500 = 1000 to 1500 + 500 = 2000
              chr1:1000-2000

    ←─ r=500 ─→
        ↓ center
    |-----------|
   1000       2000
    
    Loci found: L1 ✓, L2 ✓

Anchor 2 center = 5500
Search region: 5500 - 500 = 5000 to 5500 + 500 = 6000
              chr1:5000-6000

    Loci found: L4 ✓

Summary:
  - Larger r: More permissive, finds more loci
  - Smaller r: More strict, only precise matches
  - Default r=2500: Good balance for most use cases
    """)


def print_all():
    """Print all diagrams"""
    print_diagram_1hop()
    print_diagram_2hop()
    print_diagram_undirected()
    print_diagram_radius()
    
    print("\n" + "="*70)
    print("KEY CONCEPTS")
    print("="*70)
    print("""
1. CENTER-BASED MATCHING
   - Loop anchors are converted to their CENTER point
   - We search ±r around this center
   - This gives precise, consistent mapping

2. SPREADING vs BACKTRACKING
   - Each hop discovers new loci
   - Next hop uses ONLY newly discovered loci
   - Source loci excluded to prevent infinite loops
   
3. DIRECTED vs UNDIRECTED
   - Directed: Models information flow (promoter → enhancer)
   - Undirected: Models connectivity (promoter ↔ enhancer)

4. HOP LEVELS
   - Hop 1: Direct neighbors
   - Hop 2: Neighbors of neighbors
   - Hop N: N-th order connections
   
   More hops → Larger network → More computation
    """)


if __name__ == "__main__":
    print_all()
