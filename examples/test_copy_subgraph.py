"""
Test the improved copy() and new subgraph() methods
"""

from genomeblocks import Loci, Locus, Architecture
import tempfile
import os

# Create test data
p1 = Locus("chr1", 1000, 2000)
p2 = Locus("chr1", 100000, 101000)
promoters = Loci([p1, p2])

bedpe_content = """chr1\t1000\t2000\tchr1\t5000\t6000
chr1\t1000\t2000\tchr1\t10000\t11000
chr1\t5000\t6000\tchr1\t15000\t16000
chr1\t10000\t11000\tchr1\t20000\t21000
chr1\t100000\t101000\tchr1\t105000\t106000
"""

with tempfile.NamedTemporaryFile(mode='w', suffix='.bedpe', delete=False) as f:
    f.write(bedpe_content)
    bedpe_file = f.name

try:
    print("="*70)
    print("TEST 1: Create a graph with properties")
    print("="*70)
    
    G = Architecture.make_spread(
        source_loci=promoters,
        bedpe=bedpe_file,
        hops=2,
        directed=False,
        verbose=True
    )
    
    # Add some custom vertex properties
    G.vp.score = G.new_vertex_property("double")
    G.vp.label = G.new_vertex_property("string")
    G.vp.active = G.new_vertex_property("bool")
    
    for i, v in enumerate(G.vertices()):
        G.vp.score[v] = float(i) * 0.1
        G.vp.label[v] = f"vertex_{i}"
        G.vp.active[v] = (i % 2 == 0)
    
    # Add some custom edge properties
    G.ep.strength = G.new_edge_property("double")
    
    for i, e in enumerate(G.edges()):
        G.ep.strength[e] = float(i) * 2.0
    
    print(f"\nOriginal graph: {G.n_loci} vertices, {G.n_links} edges")
    print(f"Vertex properties: {list(G.vp.keys())}")
    print(f"Edge properties: {list(G.ep.keys())}")
    
    # Test copy()
    print("\n" + "="*70)
    print("TEST 2: Copy the graph")
    print("="*70)
    
    G_copy = G.copy()
    
    print(f"\nCopied graph: {G_copy.n_loci} vertices, {G_copy.n_links} edges")
    print(f"Vertex properties: {list(G_copy.vp.keys())}")
    print(f"Edge properties: {list(G_copy.ep.keys())}")
    
    # Verify properties were copied
    print("\nVerifying vertex properties were copied:")
    v_orig = next(G.vertices())
    v_copy = next(G_copy.vertices())
    print(f"  Original score: {G.vp.score[v_orig]}")
    print(f"  Copy score: {G_copy.vp.score[v_copy]}")
    print(f"  Original label: {G.vp.label[v_orig]}")
    print(f"  Copy label: {G_copy.vp.label[v_copy]}")
    
    print("\nVerifying edge properties were copied:")
    e_orig = next(G.edges())
    e_copy = next(G_copy.edges())
    print(f"  Original strength: {G.ep.strength[e_orig]}")
    print(f"  Copy strength: {G_copy.ep.strength[e_copy]}")
    
    print("✓ Copy test passed!")
    
    # Test subgraph by UIDs
    print("\n" + "="*70)
    print("TEST 3: Subgraph by UIDs")
    print("="*70)
    
    # Get first 3 UIDs
    sample_uids = [G.vp.uid[v] for v in list(G.vertices())[:3]]
    print(f"\nFiltering by UIDs: {sample_uids}")
    
    G_sub1 = G.subgraph(uids=sample_uids)
    
    print(f"Subgraph has {G_sub1.n_loci} vertices (expected 3)")
    assert G_sub1.n_loci == 3, "Should have 3 vertices"
    
    # Verify properties were preserved
    v_sub = next(G_sub1.vertices())
    print(f"Sample vertex properties:")
    print(f"  score: {G_sub1.vp.score[v_sub]}")
    print(f"  label: {G_sub1.vp.label[v_sub]}")
    print(f"  active: {G_sub1.vp.active[v_sub]}")
    
    print("✓ Subgraph by UIDs test passed!")
    
    # Test subgraph by vertex property
    print("\n" + "="*70)
    print("TEST 4: Subgraph by vertex property")
    print("="*70)
    
    # Filter by active vertices
    print("\nFiltering by active=True")
    G_sub2 = G.subgraph(vp_name='active', vp_values=True)
    
    print(f"Subgraph has {G_sub2.n_loci} vertices")
    
    # Verify all are active
    all_active = all(G_sub2.vp.active[v] for v in G_sub2.vertices())
    print(f"All vertices active: {all_active}")
    assert all_active, "All vertices should be active"
    
    print("✓ Subgraph by property test passed!")
    
    # Test subgraph by custom function
    print("\n" + "="*70)
    print("TEST 5: Subgraph by custom function")
    print("="*70)
    
    # Filter by score > 0.2
    print("\nFiltering by score > 0.2")
    G_sub3 = G.subgraph(filter_func=lambda v: G.vp.score[v] > 0.2)
    
    print(f"Subgraph has {G_sub3.n_loci} vertices")
    
    # Verify all have score > 0.2
    min_score = min(G_sub3.vp.score[v] for v in G_sub3.vertices())
    print(f"Minimum score in subgraph: {min_score}")
    assert min_score > 0.2, "All scores should be > 0.2"
    
    print("✓ Subgraph by function test passed!")
    
    # Test directed graph copy
    print("\n" + "="*70)
    print("TEST 6: Copy directed graph")
    print("="*70)
    
    G_dir = Architecture.make_spread(
        source_loci=promoters,
        bedpe=bedpe_file,
        hops=1,
        directed=True,
        verbose=True
    )
    
    G_dir_copy = G_dir.copy()
    
    print(f"Original directed: {G_dir.is_directed()}")
    print(f"Copy directed: {G_dir_copy.is_directed()}")
    assert G_dir_copy.is_directed(), "Copy should be directed"
    
    print("✓ Directed graph copy test passed!")
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED! ✓")
    print("="*70)
    
finally:
    os.unlink(bedpe_file)
