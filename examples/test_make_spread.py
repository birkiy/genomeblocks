"""
Quick test for Architecture.make_spread() functionality
"""

def test_make_spread_basic():
    """Test basic make_spread functionality"""
    from genomeblocks import Loci, Locus, Architecture
    import tempfile
    import os
    
    # Create test loci
    all_loci = Loci([
        Locus("chr1", 1000, 2000, name="L1"),
        Locus("chr1", 5000, 6000, name="L2"),
        Locus("chr1", 10000, 11000, name="L3"),
        Locus("chr1", 15000, 16000, name="L4"),
        Locus("chr1", 20000, 21000, name="L5"),
    ])
    
    # Source loci (promoters)
    source = Loci([all_loci[0]])  # Just L1
    
    # Create temporary BEDPE file with loops
    # L1 loops to L2, L2 loops to L3, L3 loops to L4
    bedpe_content = """chr1\t1000\t2000\tchr1\t5000\t6000
chr1\t5000\t6000\tchr1\t10000\t11000
chr1\t10000\t11000\tchr1\t15000\t16000
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bedpe', delete=False) as f:
        f.write(bedpe_content)
        bedpe_file = f.name
    
    try:
        # Test 1 hop
        print("\n=== Test 1: 1 hop, undirected ===")
        G1 = Architecture.make_spread(
            loci=all_loci,
            source_loci=source,
            bedpe=bedpe_file,
            name="Test1Hop",
            r=1000,
            hops=1,
            directed=False,
            verbose=True
        )
        
        assert G1.n_loci >= 2, f"Expected at least 2 loci, got {G1.n_loci}"
        assert "L1" in G1, "L1 should be in graph"
        assert "L2" in G1, "L2 should be in graph (1 hop from L1)"
        print(f"✓ 1 hop test passed: {G1.n_loci} loci, {G1.n_links} links")
        
        # Test 2 hops
        print("\n=== Test 2: 2 hops, undirected ===")
        G2 = Architecture.make_spread(
            loci=all_loci,
            source_loci=source,
            bedpe=bedpe_file,
            name="Test2Hop",
            r=1000,
            hops=2,
            directed=False,
            verbose=True
        )
        
        assert G2.n_loci >= 3, f"Expected at least 3 loci, got {G2.n_loci}"
        assert "L3" in G2, "L3 should be in graph (2 hops from L1)"
        print(f"✓ 2 hop test passed: {G2.n_loci} loci, {G2.n_links} links")
        
        # Test 3 hops
        print("\n=== Test 3: 3 hops, undirected ===")
        G3 = Architecture.make_spread(
            loci=all_loci,
            source_loci=source,
            bedpe=bedpe_file,
            name="Test3Hop",
            r=1000,
            hops=3,
            directed=False,
            verbose=True
        )
        
        assert G3.n_loci >= 4, f"Expected at least 4 loci, got {G3.n_loci}"
        assert "L4" in G3, "L4 should be in graph (3 hops from L1)"
        print(f"✓ 3 hop test passed: {G3.n_loci} loci, {G3.n_links} links")
        
        # Test directed graph
        print("\n=== Test 4: 1 hop, directed ===")
        G4 = Architecture.make_spread(
            loci=all_loci,
            source_loci=source,
            bedpe=bedpe_file,
            name="TestDirected",
            r=1000,
            hops=1,
            directed=True,
            verbose=True
        )
        
        assert G4.is_directed(), "Graph should be directed"
        print(f"✓ Directed test passed: {G4.n_loci} loci, {G4.n_links} links")
        
        print("\n=== All tests passed! ✓ ===")
        
    finally:
        # Cleanup
        os.unlink(bedpe_file)


if __name__ == "__main__":
    test_make_spread_basic()
