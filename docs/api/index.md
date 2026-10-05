---
title: API Reference
layout: default
nav_order: 6
has_children: true
permalink: /api/
---

# API Reference

Signatures and short descriptions for every public name, grouped by module. Use Ctrl-F / the search box to jump around.

| Module | Public names |
|---|---|
| [`genomeblocks.locus`]({{ '/api/locus/' | relative_url }}) | `Locus`, `Exon`, `CDS`, `UTR` |
| [`genomeblocks.loci`]({{ '/api/loci/' | relative_url }}) | `Loci` |
| [`genomeblocks.genes`]({{ '/api/genes/' | relative_url }}) | `Gene`, `Transcript`, `Genes` |
| [`genomeblocks.architecture`]({{ '/api/architecture/' | relative_url }}) | `Architecture` |
| [`genomeblocks.signal`]({{ '/api/signal/' | relative_url }}) | `signal`, `tmm` (plotting: `plot_heatmap`, `plot_profiles`, `compare_heatmap` in `signal_draw`) |
| [`genomeblocks.browserview`]({{ '/api/browser/' | relative_url }}) | `browser` (re-exported as `genomeblocks.browser`) |
| [`genomeblocks.bedpe`]({{ '/api/bedpe/' | relative_url }}) | `Pair`, `read_bedpe`, `pair_to_bed`, `pairs_to_frame`, `pairs_to_bedpe`, `count_pairs`, `count_pairs_2d` |
| [`genomeblocks.motifs`]({{ '/api/motifs/' | relative_url }}) | `make_genome`, `scan_motifs` (matrix / masked / differential / archetypes in the module) |
| [`genomeblocks.atlas`]({{ '/api/atlas/' | relative_url }}) | `Atlas` |
