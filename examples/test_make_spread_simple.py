"""
Simple test for the NEW Architecture.make_spread() - no pre-defined CRE collection needed!
"""

def test_spread_without_all_cres():
    """Test make_spread where network is built entirely from BEDPE anchors"""
    from genomeblocks import Loci, Locus, Architecture
    import tempfile
    import os
    
    # Create source loci (promoters)
    p1 = Locus("chr1", 1000, 2000)
    p2 = Locus("chr1", 100000, 101000)
    promoters = Loci([p1, p2])
    
    print(f"\nSource promoters: {len(promoters)}")
    for p in promoters:
        print(f"  {p.uid}: {p.chrom}:{p.start}-{p.end}")
    
    # Create BEDPE file with loops
    # P1 (center=1500) connects to enhancers at 5500, 10500
    # Those enhancers connect to more distal sites
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
        # Test 1 hop - should find direct connections
        print("\n" + "="*70)
        print("TEST 1: 1-hop spread")
        print("="*70)
        G1 = Architecture.make_spread(
            source_loci=promoters,
            bedpe=bedpe_file,
            name="Test1Hop",
            r=1000,
            hops=1,
            directed=False,
            verbose=True
        )
        
        print(f"\nGraph vertices (anchor centers):")
        for v in G1.vertices():
            uid = G1.vp.uid[v]
            print(f"  {uid}")
        
        # Should have: chr1:1500 (P1), chr1:100500 (P2), chr1:5500, chr1:10500, chr1:105500
        assert G1.n_loci >= 3, f"Expected at least 3 vertices, got {G1.n_loci}"
        assert "chr1:1500" in G1 or any("1500" in uid for uid in [G1.vp.uid[v] for v in G1.vertices()]), "P1 center should be in graph"
        print(f"✓ 1-hop test passed: {G1.n_loci} anchors, {G1.n_links} links\n")
        
        # Test 2 hops
        print("="*70)
        print("TEST 2: 2-hop spread")
        print("="*70)
        G2 = Architecture.make_spread(
            source_loci=promoters,
            bedpe=bedpe_file,
            name="Test2Hop",
            r=1000,
            hops=2,
            directed=False,
            verbose=True
        )
        
        print(f"\nGraph vertices:")
        for v in G2.vertices():
            uid = G2.vp.uid[v]
            print(f"  {uid}")
        
        # Should have more vertices now (15500, 20500 from hop 2)
        assert G2.n_loci > G1.n_loci, f"2-hop should have more vertices than 1-hop"
        print(f"✓ 2-hop test passed: {G2.n_loci} anchors, {G2.n_links} links\n")
        
        # Test directed
        print("="*70)
        print("TEST 3: Directed graph")
        print("="*70)
        G3 = Architecture.make_spread(
            source_loci=promoters,
            bedpe=bedpe_file,
            name="TestDirected",
            r=1000,
            hops=1,
            directed=True,
            verbose=True
        )
        
        assert G3.is_directed(), "Graph should be directed"
        print(f"✓ Directed test passed: {G3.n_loci} anchors, {G3.n_links} links\n")
        
        print("="*70)
        print("ALL TESTS PASSED! ✓")
        print("="*70)
        print("\nKEY INSIGHT:")
        print("  - Vertices are loop anchor CENTERS (chr:pos format)")
        print("  - NO need for pre-defined CRE collection")
        print("  - Network grows organically from BEDPE")
        print("  - Perfect when you don't have a curated CRE list!")
        
    finally:
        os.unlink(bedpe_file)


if __name__ == "__main__":
    test_spread_without_all_cres()
