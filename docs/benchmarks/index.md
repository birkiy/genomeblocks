---
title: Benchmarks
layout: default
nav_order: 9
permalink: /benchmarks/
---

# Benchmarks
{: .no_toc }

How fast each building block is against the tools people would otherwise use,
measured on genomeblocks 1.1. Every comparison first checks that the engines
return the same answer; timings are medians of 3–5 runs after a warm-up.
{: .fs-5 .fw-300 }

<p class="gb-env">Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz · 64 cores · 811 GB RAM ·
Python 3.12.12 · synthetic hg38-shaped data (seeded) · benchmarks ran one at a time and use at
most 8 worker processes</p>

<div class="gb-tiles">
<div class="gb-tile green"><span class="tile-label">One overlap lookup</span><span class="tile-value">2.5 µs</span><span class="tile-ctx">on 100k indexed peaks: 454× faster than pyranges, 2,200× than bioframe</span></div>
<div class="gb-tile green"><span class="tile-label">bigWig → heatmap matrix</span><span class="tile-value">98×</span><span class="tile-ctx">faster than deepTools computeMatrix on one core (12k vs 124 regions/s)</span></div>
<div class="gb-tile navy"><span class="tile-label">Enrichment vs 500 peak files</span><span class="tile-value">21 ms</span><span class="tile-ctx">per 20k-peak query including Fisher tests: 256× faster than looping over the tracks</span></div>
<div class="gb-tile navy"><span class="tile-label">Motif scanning</span><span class="tile-value">127×</span><span class="tile-ctx">faster than MEME FIMO, 171× faster than Biopython</span></div>
<div class="gb-tile purple"><span class="tile-label">Hi-C pair counting</span><span class="tile-value">1.6M/s</span><span class="tile-ctx">read pairs into 50 kb windows: 2.7× the rate of cooler cload</span></div>
<div class="gb-tile purple"><span class="tile-label">Whole pipeline, fresh process</span><span class="tile-value">26 s → 6.2 s</span><span class="tile-ctx">classic vs columnar: load CREs and genes, build the graph, add Hi-C, normalise, annotate, hubs, save</span></div>
</div>

1. TOC
{:toc}

## Loci
{: .sec-green #loci }

A single overlap query against an indexed set is where the object model pays off: **2.5 µs** per lookup with cgranges and 3.1 µs with the pure-Python index, against 1.1 ms for pyranges, which has to slice a dataframe per call. That is the cost of every interactive question ("what overlaps this peak?").

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 182" role="img" aria-label="Single overlap query latency" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="156"/>
<text class="tick" x="300.0" y="172" text-anchor="middle">1 µs</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="156"/>
<text class="tick" x="421.0" y="172" text-anchor="middle">10 µs</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="156"/>
<text class="tick" x="542.0" y="172" text-anchor="middle">100 µs</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="156"/>
<text class="tick" x="663.0" y="172" text-anchor="middle">1 ms</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="156"/>
<text class="tick" x="784.0" y="172" text-anchor="middle">10 ms</text>
<g class="row gb green"><title>genomeblocks (cgranges): 2.5 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks (cgranges)</text>
<circle class="dot" cx="348.3" cy="22.0" r="5.5"/>
<text class="val" x="359.3" y="26.0">2.5 µs</text>
</g>
<g class="row gb green"><title>genomeblocks (pure-Python fallback): 3.1 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks (pure-Python fallback)</text>
<circle class="dot" cx="360.1" cy="50.0" r="5.5"/>
<text class="val" x="371.1" y="54.0">3.1 µs</text>
</g>
<g class="row"><title>intervaltree: 8.2 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">intervaltree</text>
<circle class="dot" cx="410.9" cy="78.0" r="5.5"/>
<text class="val" x="421.9" y="82.0">8.2 µs</text>
</g>
<g class="row"><title>pyranges: 1.1 ms · per query, 100k-peak index</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="669.8" cy="106.0" r="5.5"/>
<text class="val" x="680.8" y="110.0">1.1 ms</text>
</g>
<g class="row"><title>bioframe: 5.5 ms · per query, 100k-peak index</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="752.7" cy="134.0" r="5.5"/>
<text class="val" x="763.7" y="138.0">5.5 ms</text>
</g>
</svg>
</div><figcaption><strong>One overlap lookup</strong> against 100,000 indexed peaks (median per query; log scale).</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 210" role="img" aria-label="A &amp; B at 1,000,000 intervals" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="184"/>
<text class="tick" x="300.0" y="200" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="184"/>
<text class="tick" x="461.3" y="200" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="184"/>
<text class="tick" x="622.7" y="200" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="184"/>
<text class="tick" x="784.0" y="200" text-anchor="middle">10 s</text>
<g class="row gb green"><title>genomeblocks.columnar: 74 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks.columnar</text>
<circle class="dot" cx="440.2" cy="22.0" r="5.5"/>
<text class="val" x="451.2" y="26.0">74 ms</text>
</g>
<g class="row"><title>pyranges: 345 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="548.0" cy="50.0" r="5.5"/>
<text class="val" x="559.0" y="54.0">345 ms</text>
</g>
<g class="row"><title>bioframe: 1.1 s · n = 1,000,000</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="628.6" cy="78.0" r="5.5"/>
<text class="val" x="639.6" y="82.0">1.1 s</text>
</g>
<g class="row gb green"><title>genomeblocks (cgranges): 1.5 s · n = 1,000,000</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks (cgranges)</text>
<circle class="dot" cx="650.8" cy="106.0" r="5.5"/>
<text class="val" x="661.8" y="110.0">1.5 s</text>
</g>
<g class="row"><title>bedtools (CLI): 1.8 s · n = 1,000,000</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">bedtools (CLI)</text>
<circle class="dot" cx="662.0" cy="134.0" r="5.5"/>
<text class="val" x="673.0" y="138.0">1.8 s</text>
</g>
<g class="row gb green"><title>genomeblocks (pure-Python fallback): 2.2 s · n = 1,000,000</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">genomeblocks (pure-Python fallback)</text>
<circle class="dot" cx="678.3" cy="162.0" r="5.5"/>
<text class="val" x="689.3" y="166.0">2.2 s</text>
</g>
</svg>
</div><figcaption><strong>A &amp; B</strong> on 1M intervals per set (rows of A overlapping B); every engine returns the same rows.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 182" role="img" aria-label="merge at 1,000,000 intervals" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="156"/>
<text class="tick" x="300.0" y="172" text-anchor="middle">100 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="156"/>
<text class="tick" x="461.3" y="172" text-anchor="middle">1 s</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="156"/>
<text class="tick" x="622.7" y="172" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="156"/>
<text class="tick" x="784.0" y="172" text-anchor="middle">100 s</text>
<g class="row"><title>pyranges: 174 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="338.7" cy="22.0" r="5.5"/>
<text class="val" x="349.7" y="26.0">174 ms</text>
</g>
<g class="row gb green"><title>genomeblocks.columnar: 408 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks.columnar</text>
<circle class="dot" cx="398.6" cy="50.0" r="5.5"/>
<text class="val" x="409.6" y="54.0">408 ms</text>
</g>
<g class="row"><title>bioframe: 1.2 s · n = 1,000,000</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="474.2" cy="78.0" r="5.5"/>
<text class="val" x="485.2" y="82.0">1.2 s</text>
</g>
<g class="row"><title>bedtools (CLI): 4.9 s · n = 1,000,000</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">bedtools (CLI)</text>
<circle class="dot" cx="572.6" cy="106.0" r="5.5"/>
<text class="val" x="583.6" y="110.0">4.9 s</text>
</g>
<g class="row gb green"><title>genomeblocks: 7.7 s · n = 1,000,000</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="604.7" cy="134.0" r="5.5"/>
<text class="val" x="615.7" y="138.0">7.7 s</text>
</g>
</svg>
</div><figcaption><strong>merge</strong> on 1M intervals per set (sort + merge of A ∪ B); every engine returns the same rows.</figcaption></figure>

On whole-set operations at 1M intervals the object-per-interval `Loci` is not the fastest: `A & B` takes 1.5 s and `merge` 7.7 s, against 345 ms and 174 ms for pyranges. The [columnar `Loci`]({{ '/design/columnar/' | relative_url }}) does the same work on numpy columns: 74 ms and 408 ms.

<details class="gb-table"><summary>All Loci timings</summary><table><thead><tr><th>op</th><th>engine</th><th>n = 1k</th><th>n = 10k</th><th>n = 100k</th><th>n = 1M</th></tr></thead><tbody><tr><td>intersect</td><td>genomeblocks (cgranges)</td><td>1.6 ms</td><td>12 ms</td><td>120 ms</td><td>1.5 s</td></tr><tr><td>intersect</td><td>genomeblocks (cgranges, warm index)</td><td>1.1 ms</td><td>9.2 ms</td><td>98 ms</td><td>1.3 s</td></tr><tr><td>intersect</td><td>genomeblocks (pure-Python fallback)</td><td>1.9 ms</td><td>16 ms</td><td>176 ms</td><td>2.2 s</td></tr><tr><td>intersect</td><td>pyranges</td><td>45 ms</td><td>41 ms</td><td>62 ms</td><td>345 ms</td></tr><tr><td>intersect</td><td>bioframe</td><td>7.3 ms</td><td>11 ms</td><td>77 ms</td><td>1.1 s</td></tr><tr><td>intersect</td><td>bedtools (CLI)</td><td>15 ms</td><td>37 ms</td><td>182 ms</td><td>1.8 s</td></tr><tr><td>intersect</td><td>intervaltree</td><td>19 ms</td><td>369 ms</td><td>20 s</td><td>–</td></tr><tr><td>intersect</td><td>naive double loop</td><td>1.7 ms</td><td>132 ms</td><td>–</td><td>–</td></tr><tr><td>intersect</td><td>genomeblocks.columnar</td><td>170 µs</td><td>700 µs</td><td>7.2 ms</td><td>74 ms</td></tr><tr><td>merge</td><td>pyranges</td><td>28 ms</td><td>32 ms</td><td>48 ms</td><td>174 ms</td></tr><tr><td>merge</td><td>bioframe</td><td>17 ms</td><td>25 ms</td><td>124 ms</td><td>1.2 s</td></tr><tr><td>merge</td><td>bedtools (CLI)</td><td>28 ms</td><td>93 ms</td><td>529 ms</td><td>4.9 s</td></tr><tr><td>merge</td><td>genomeblocks</td><td>1.7 ms</td><td>26 ms</td><td>629 ms</td><td>7.7 s</td></tr><tr><td>merge</td><td>genomeblocks.columnar</td><td>423 µs</td><td>3.0 ms</td><td>33 ms</td><td>408 ms</td></tr><tr><td>make</td><td>pyranges</td><td>5.4 ms</td><td>13 ms</td><td>92 ms</td><td>866 ms</td></tr><tr><td>make</td><td>genomeblocks</td><td>1.2 ms</td><td>11 ms</td><td>150 ms</td><td>2.0 s</td></tr><tr><td>make</td><td>pandas.read_csv</td><td>1.5 ms</td><td>4.8 ms</td><td>29 ms</td><td>254 ms</td></tr><tr><td>make</td><td>genomeblocks.columnar</td><td>2.6 ms</td><td>5.2 ms</td><td>57 ms</td><td>333 ms</td></tr></tbody></table></details>

## Signal
{: .sec-green #signal }

Filling a heatmap matrix (200 bins, ±3 kb) is one native call per region and track: **12k regions/s** on one core, 98× deepTools `computeMatrix` and 373× pyBigWig's binned `stats()`. Every engine returns the same matrix (deepTools differs only at bin edges, r ≥ 0.99).

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 266" role="img" aria-label="Regions per second, one core" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="240"/>
<text class="tick" x="300.0" y="256" text-anchor="middle">10</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="240"/>
<text class="tick" x="421.0" y="256" text-anchor="middle">100</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="240"/>
<text class="tick" x="542.0" y="256" text-anchor="middle">1k</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="240"/>
<text class="tick" x="663.0" y="256" text-anchor="middle">10k</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="240"/>
<text class="tick" x="784.0" y="256" text-anchor="middle">100k</text>
<g class="row gb green"><title>genomeblocks · pybigtools (zoom, exact=False): 14k · 5,000 regions</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks · pybigtools (zoom, exact=False)</text>
<circle class="dot" cx="681.0" cy="22.0" r="5.5"/>
<text class="val" x="692.0" y="26.0">14k</text>
</g>
<g class="row gb green"><title>genomeblocks · pybigtools (exact): 12k · 5,000 regions</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks · pybigtools (exact)</text>
<circle class="dot" cx="673.5" cy="50.0" r="5.5"/>
<text class="val" x="684.5" y="54.0">12k</text>
</g>
<g class="row"><title>pybigtools values() + numpy binning: 6.6k · 5,000 regions</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">pybigtools values() + numpy binning</text>
<circle class="dot" cx="641.4" cy="78.0" r="5.5"/>
<text class="val" x="652.4" y="82.0">6.6k</text>
</g>
<g class="row"><title>pyBigWig values() + numpy binning: 4.1k · 5,000 regions</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyBigWig values() + numpy binning</text>
<circle class="dot" cx="615.9" cy="106.0" r="5.5"/>
<text class="val" x="626.9" y="110.0">4.1k</text>
</g>
<g class="row gb green"><title>genomeblocks · pure-Python reader: 2.7k · 5,000 regions</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">genomeblocks · pure-Python reader</text>
<circle class="dot" cx="594.5" cy="134.0" r="5.5"/>
<text class="val" x="605.5" y="138.0">2.7k</text>
</g>
<g class="row"><title>deepTools computeMatrix (-p 1): 124 · 5,000 regions</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">deepTools computeMatrix (-p 1)</text>
<circle class="dot" cx="432.4" cy="162.0" r="5.5"/>
<text class="val" x="443.4" y="166.0">124</text>
</g>
<g class="row"><title>pyBigWig stats(nBins) zoom: 36 · 300 regions</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">pyBigWig stats(nBins) zoom</text>
<circle class="dot" cx="367.9" cy="190.0" r="5.5"/>
<text class="val" x="378.9" y="194.0">36</text>
</g>
<g class="row"><title>pyBigWig stats(nBins) exact: 33 · 300 regions</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">pyBigWig stats(nBins) exact</text>
<circle class="dot" cx="362.2" cy="218.0" r="5.5"/>
<text class="val" x="373.2" y="222.0">33</text>
</g>
</svg>
</div><figcaption><strong>Regions per second</strong> into a 200-bin matrix from one bigWig, one core (higher is better; log scale).</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 300" role="img" aria-label="Throughput vs number of regions" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="64" y1="256.0" x2="650" y2="256.0"/>
<text class="tick" x="56" y="260.0" text-anchor="end">1k</text>
<line class="grid" x1="64" y1="135.0" x2="650" y2="135.0"/>
<text class="tick" x="56" y="139.0" text-anchor="end">10k</text>
<line class="grid" x1="64" y1="14.0" x2="650" y2="14.0"/>
<text class="tick" x="56" y="18.0" text-anchor="end">100k</text>
<text class="tick" x="64.0" y="274" text-anchor="middle">1k</text>
<text class="tick" x="357.0" y="274" text-anchor="middle">10k</text>
<text class="tick" x="650.0" y="274" text-anchor="middle">100k</text>
<line class="axis" x1="64" y1="256" x2="650" y2="256"/>
<text class="axl" x="357.0" y="294" text-anchor="middle">regions</text>
<g class="line"><title>genomeblocks, exact=False</title><polyline points="64.0,134.9 268.8,117.3 445.2,81.9 561.8,58.2"/>
<circle cx="64.0" cy="134.9" r="4"><title>genomeblocks, exact=False · 1k: 10k</title></circle>
<circle cx="268.8" cy="117.3" r="4"><title>genomeblocks, exact=False · 5k: 14k</title></circle>
<circle cx="445.2" cy="81.9" r="4"><title>genomeblocks, exact=False · 20k: 27k</title></circle>
<circle cx="561.8" cy="58.2" r="4"><title>genomeblocks, exact=False · 50k: 43k</title></circle>
</g>
<g class="line"><title>genomeblocks, pure-Python</title><polyline points="64.0,218.7 268.8,205.0 445.2,199.9 561.8,191.8"/>
<circle cx="64.0" cy="218.7" r="4"><title>genomeblocks, pure-Python · 1k: 2k</title></circle>
<circle cx="268.8" cy="205.0" r="4"><title>genomeblocks, pure-Python · 5k: 2.6k</title></circle>
<circle cx="445.2" cy="199.9" r="4"><title>genomeblocks, pure-Python · 20k: 2.9k</title></circle>
<circle cx="561.8" cy="191.8" r="4"><title>genomeblocks, pure-Python · 50k: 3.4k</title></circle>
</g>
<g class="line"><title>pyBigWig values()</title><polyline points="64.0,181.8 268.8,188.2 445.2,182.8 561.8,187.4"/>
<circle cx="64.0" cy="181.8" r="4"><title>pyBigWig values() · 1k: 4.1k</title></circle>
<circle cx="268.8" cy="188.2" r="4"><title>pyBigWig values() · 5k: 3.6k</title></circle>
<circle cx="445.2" cy="182.8" r="4"><title>pyBigWig values() · 20k: 4k</title></circle>
<circle cx="561.8" cy="187.4" r="4"><title>pyBigWig values() · 50k: 3.7k</title></circle>
</g>
<g class="line gb green"><title>genomeblocks (default)</title><polyline points="64.0,133.9 268.8,118.1 445.2,80.9 561.8,64.9"/>
<circle cx="64.0" cy="133.9" r="4"><title>genomeblocks (default) · 1k: 10k</title></circle>
<circle cx="268.8" cy="118.1" r="4"><title>genomeblocks (default) · 5k: 14k</title></circle>
<circle cx="445.2" cy="80.9" r="4"><title>genomeblocks (default) · 20k: 28k</title></circle>
<circle cx="561.8" cy="64.9" r="4"><title>genomeblocks (default) · 50k: 38k</title></circle>
</g>
<text class="lab" x="660" y="62.2">genomeblocks, exact=False</text>
<line class="leader" x1="567.8" y1="58.2" x2="656" y2="58.2"/>
<text class="lab gbl" x="660" y="76.2">genomeblocks (default)</text>
<line class="leader" x1="567.8" y1="72.2" x2="656" y2="72.2"/>
<text class="lab" x="660" y="191.4">pyBigWig values()</text>
<line class="leader" x1="567.8" y1="187.4" x2="656" y2="187.4"/>
<text class="lab" x="660" y="205.4">genomeblocks, pure-Python</text>
<line class="leader" x1="567.8" y1="201.4" x2="656" y2="201.4"/>
</svg>
</div><figcaption><strong>Throughput vs number of regions</strong> (regions/s, log–log). The native path speeds up with batch size as fixed costs amortise; per-base readers stay flat.</figcaption></figure>

Many tracks scale with processes: 16 bigWigs × 5,000 regions take 5.0 s sequentially and 1.5 s with `workers=8`. Threads do not help (pybigtools serialises Python threads), which is why `signal()` uses processes.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 322" role="img" aria-label="16 tracks, parallel" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="296"/>
<text class="tick" x="300.0" y="312" text-anchor="middle">100 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="296"/>
<text class="tick" x="421.0" y="312" text-anchor="middle">1 s</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="296"/>
<text class="tick" x="542.0" y="312" text-anchor="middle">10 s</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="296"/>
<text class="tick" x="663.0" y="312" text-anchor="middle">100 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="296"/>
<text class="tick" x="784.0" y="312" text-anchor="middle">1,000 s</text>
<g class="row gb green"><title>genomeblocks · processes · 8 workers: 1.5 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks · processes · 8 workers</text>
<circle class="dot" cx="440.9" cy="22.0" r="5.5"/>
<text class="val" x="451.9" y="26.0">1.5 s</text>
</g>
<g class="row gb green"><title>genomeblocks · processes · 4 workers: 2.1 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks · processes · 4 workers</text>
<circle class="dot" cx="459.1" cy="50.0" r="5.5"/>
<text class="val" x="470.1" y="54.0">2.1 s</text>
</g>
<g class="row gb green"><title>genomeblocks · processes · 2 workers: 3.8 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">genomeblocks · processes · 2 workers</text>
<circle class="dot" cx="491.2" cy="78.0" r="5.5"/>
<text class="val" x="502.2" y="82.0">3.8 s</text>
</g>
<g class="row gb green"><title>genomeblocks · processes · 1 worker: 5.0 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks · processes · 1 worker</text>
<circle class="dot" cx="505.8" cy="106.0" r="5.5"/>
<text class="val" x="516.8" y="110.0">5.0 s</text>
</g>
<g class="row"><title>threads (same chunks) · 2 workers: 19 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">threads (same chunks) · 2 workers</text>
<circle class="dot" cx="575.6" cy="134.0" r="5.5"/>
<text class="val" x="586.6" y="138.0">19 s</text>
</g>
<g class="row"><title>threads (same chunks) · 4 workers: 20 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">threads (same chunks) · 4 workers</text>
<circle class="dot" cx="577.6" cy="162.0" r="5.5"/>
<text class="val" x="588.6" y="166.0">20 s</text>
</g>
<g class="row"><title>threads (same chunks) · 8 workers: 20 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">threads (same chunks) · 8 workers</text>
<circle class="dot" cx="577.7" cy="190.0" r="5.5"/>
<text class="val" x="588.7" y="194.0">20 s</text>
</g>
<g class="row"><title>deepTools computeMatrix · 8 workers: 105 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">deepTools computeMatrix · 8 workers</text>
<circle class="dot" cx="665.6" cy="218.0" r="5.5"/>
<text class="val" x="676.6" y="222.0">105 s</text>
</g>
<g class="row"><title>deepTools computeMatrix · 4 workers: 3.1 min · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="232.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="246.0" x2="784" y2="246.0"/>
<text class="lab" x="286" y="250.0" text-anchor="end">deepTools computeMatrix · 4 workers</text>
<circle class="dot" cx="695.6" cy="246.0" r="5.5"/>
<text class="val" x="706.6" y="250.0">3.1 min</text>
</g>
<g class="row"><title>deepTools computeMatrix · 1 worker: 10.0 min · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="260.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="274.0" x2="784" y2="274.0"/>
<text class="lab" x="286" y="278.0" text-anchor="end">deepTools computeMatrix · 1 worker</text>
<circle class="dot" cx="757.0" cy="274.0" r="5.5"/>
<text class="val" x="768.0" y="278.0">10.0 min</text>
</g>
</svg>
</div><figcaption><strong>16 bigWigs × 5,000 regions</strong>: worker processes vs threads vs deepTools <code>-p</code> (wall time; log scale).</figcaption></figure>

## Atlas
{: .sec-navy #atlas }

Testing a 20,000-peak query against 500 peak files is one sparse row-sum plus a vectorised Fisher test: **21 ms** per query, including building the result table, 256× faster than a per-track loop over prebuilt cgranges indexes.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="Query against 500 tracks" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">1 ms</text>
<line class="grid" x1="396.8" y1="8" x2="396.8" y2="128"/>
<text class="tick" x="396.8" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="493.6" y1="8" x2="493.6" y2="128"/>
<text class="tick" x="493.6" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="590.4" y1="8" x2="590.4" y2="128"/>
<text class="tick" x="590.4" y="144" text-anchor="middle">1 s</text>
<line class="grid" x1="687.2" y1="8" x2="687.2" y2="128"/>
<text class="tick" x="687.2" y="144" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">100 s</text>
<g class="row gb navy"><title>Atlas sparse row-sum only: 3.8 ms · 500 tracks</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">Atlas sparse row-sum only</text>
<circle class="dot" cx="355.7" cy="22.0" r="5.5"/>
<text class="val" x="366.7" y="26.0">3.8 ms</text>
</g>
<g class="row gb navy"><title>Atlas.search (with Fisher + table): 21 ms · 500 tracks</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">Atlas.search (with Fisher + table)</text>
<circle class="dot" cx="428.1" cy="50.0" r="5.5"/>
<text class="val" x="439.1" y="54.0">21 ms</text>
</g>
<g class="row"><title>per-track cgranges loop (prebuilt): 5.4 s · 500 tracks</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">per-track cgranges loop (prebuilt)</text>
<circle class="dot" cx="661.2" cy="78.0" r="5.5"/>
<text class="val" x="672.2" y="82.0">5.4 s</text>
</g>
<g class="row"><title>pyranges overlap(), per track: 21 s · 500 tracks</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyranges overlap(), per track</text>
<circle class="dot" cx="718.4" cy="106.0" r="5.5"/>
<text class="val" x="729.4" y="110.0">21 s</text>
</g>
</svg>
</div><figcaption><strong>One query against 500 peak files</strong> (wall time; log scale). GIGGLE was not re-run for 1.1 (it has to be built from source).</figcaption></figure>

<details class="gb-table"><summary>Index build, load, accuracy</summary><table><thead><tr><th>step</th><th>bin size</th><th>workers</th><th>time</th><th>index size</th></tr></thead><tbody><tr><td>build</td><td>200 bp</td><td>3</td><td>7.8 s</td><td>190 MB</td></tr><tr><td>build</td><td>1,000 bp</td><td>1</td><td>7.9 s</td><td>68 MB</td></tr><tr><td>build</td><td>1,000 bp</td><td>3</td><td>3.5 s</td><td>68 MB</td></tr><tr><td>build</td><td>5,000 bp</td><td>3</td><td>2.8 s</td><td>43 MB</td></tr><tr><td>load (Atlas.load (npz))</td><td></td><td></td><td>271 ms</td><td></td></tr><tr><td>load (fresh python: import + load + search)</td><td></td><td></td><td>1.6 s</td><td></td></tr><tr><td>accuracy vs exact overlaps</td><td>200 bp</td><td></td><td>Spearman 0.999</td><td>top-25 shared 21</td></tr><tr><td>accuracy vs exact overlaps</td><td>1,000 bp</td><td></td><td>Spearman 0.999</td><td>top-25 shared 20</td></tr><tr><td>accuracy vs exact overlaps</td><td>5,000 bp</td><td></td><td>Spearman 0.993</td><td>top-25 shared 18</td></tr><tr><td>bootstrap (10 shuffles)</td><td></td><td></td><td>12 ms per shuffle</td><td></td></tr><tr><td>bootstrap (100 shuffles)</td><td></td><td></td><td>8.6 ms per shuffle</td><td></td></tr></tbody></table></details>

## Motifs
{: .sec-navy #motifs }

Counting hits of 100 JASPAR motifs in 1,000 windows of 500 bp: **57 ms** with genomeblocks (one block scan per motif), 127× faster than FIMO. The 993 hits are identical to those of lightmotif per window, MOODS, numpy and Biopython; FIMO reports 996.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 210" role="img" aria-label="Motif scanning engines" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="184"/>
<text class="tick" x="300.0" y="200" text-anchor="middle">10 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="184"/>
<text class="tick" x="421.0" y="200" text-anchor="middle">100 ms</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="184"/>
<text class="tick" x="542.0" y="200" text-anchor="middle">1 s</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="184"/>
<text class="tick" x="663.0" y="200" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="184"/>
<text class="tick" x="784.0" y="200" text-anchor="middle">100 s</text>
<g class="row gb navy"><title>genomeblocks scan_motifs_matrix: 57 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks scan_motifs_matrix</text>
<circle class="dot" cx="391.3" cy="22.0" r="5.5"/>
<text class="val" x="402.3" y="26.0">57 ms</text>
</g>
<g class="row"><title>MOODS (C++): 66 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">MOODS (C++)</text>
<circle class="dot" cx="399.5" cy="50.0" r="5.5"/>
<text class="val" x="410.5" y="54.0">66 ms</text>
</g>
<g class="row"><title>lightmotif, re-striped per motif: 280 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">lightmotif, re-striped per motif</text>
<circle class="dot" cx="475.1" cy="78.0" r="5.5"/>
<text class="val" x="486.1" y="82.0">280 ms</text>
</g>
<g class="row"><title>numpy sliding window: 5.2 s · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">numpy sliding window</text>
<circle class="dot" cx="628.8" cy="106.0" r="5.5"/>
<text class="val" x="639.8" y="110.0">5.2 s</text>
</g>
<g class="row"><title>MEME FIMO (CLI): 7.2 s · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">MEME FIMO (CLI)</text>
<circle class="dot" cx="645.7" cy="134.0" r="5.5"/>
<text class="val" x="656.7" y="138.0">7.2 s</text>
</g>
<g class="row"><title>Biopython PSSM.search: 9.7 s · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">Biopython PSSM.search</text>
<circle class="dot" cx="661.5" cy="162.0" r="5.5"/>
<text class="val" x="672.5" y="166.0">9.7 s</text>
</g>
</svg>
</div><figcaption><strong>Motif scanning engines</strong> on identical PSSMs and windows (wall time; log scale).</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 182" role="img" aria-label="Whole library" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="156"/>
<text class="tick" x="300.0" y="172" text-anchor="middle">100 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="156"/>
<text class="tick" x="421.0" y="172" text-anchor="middle">1 s</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="156"/>
<text class="tick" x="542.0" y="172" text-anchor="middle">10 s</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="156"/>
<text class="tick" x="663.0" y="172" text-anchor="middle">100 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="156"/>
<text class="tick" x="784.0" y="172" text-anchor="middle">1,000 s</text>
<g class="row gb navy"><title>genomeblocks scan_motifs_matrix · 1k windows: 538 ms · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks scan_motifs_matrix · 1k windows</text>
<circle class="dot" cx="388.4" cy="22.0" r="5.5"/>
<text class="val" x="399.4" y="26.0">538 ms</text>
</g>
<g class="row"><title>MOODS (C++) · 1k windows: 738 ms · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">MOODS (C++) · 1k windows</text>
<circle class="dot" cx="405.1" cy="50.0" r="5.5"/>
<text class="val" x="416.1" y="54.0">738 ms</text>
</g>
<g class="row"><title>MOODS (C++) · 5k windows: 2.3 s · 5,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">MOODS (C++) · 5k windows</text>
<circle class="dot" cx="464.5" cy="78.0" r="5.5"/>
<text class="val" x="475.5" y="82.0">2.3 s</text>
</g>
<g class="row gb navy"><title>genomeblocks scan_motifs_matrix · 5k windows: 2.5 s · 5,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks scan_motifs_matrix · 5k windows</text>
<circle class="dot" cx="469.3" cy="106.0" r="5.5"/>
<text class="val" x="480.3" y="110.0">2.5 s</text>
</g>
<g class="row"><title>MEME FIMO (CLI) · 1k windows: 72 s · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">MEME FIMO (CLI) · 1k windows</text>
<circle class="dot" cx="646.0" cy="134.0" r="5.5"/>
<text class="val" x="657.0" y="138.0">72 s</text>
</g>
</svg>
</div><figcaption><strong>The whole JASPAR library</strong> (1,019 motifs) on one core (wall time; log scale). Against MOODS: genomeblocks is 1.4× faster at 1,000 windows, MOODS is 1.1× faster at 5,000 windows. FIMO was timed at 1,000 windows only.</figcaption></figure>

<details class="gb-table"><summary>Worker processes</summary><table><thead><tr><th>workers</th><th>time</th><th>speed-up</th></tr></thead><tbody><tr><td>1</td><td>2.1 s</td><td>1.0×</td></tr><tr><td>2</td><td>1.7 s</td><td>1.2×</td></tr><tr><td>4</td><td>1.3 s</td><td>1.6×</td></tr><tr><td>8</td><td>1.2 s</td><td>1.7×</td></tr></tbody></table></details>

## Hi-C pairs
{: .sec-purple #pairs }

Streaming 5M read pairs into 50 kb windows × partner chromosome runs at **1.6M pairs/s**, close to the speed of just parsing the file, and 2.7× the rate of `cooler cload`.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 210" role="img" aria-label="Pairs per second" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="184"/>
<text class="tick" x="300.0" y="200" text-anchor="middle">100k</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="184"/>
<text class="tick" x="542.0" y="200" text-anchor="middle">1M</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="184"/>
<text class="tick" x="784.0" y="200" text-anchor="middle">10M</text>
<g class="row"><title>parse only (I/O floor): 2.8M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">parse only (I/O floor)</text>
<circle class="dot" cx="650.2" cy="22.0" r="5.5"/>
<text class="val" x="661.2" y="26.0">2.8M</text>
</g>
<g class="row gb purple"><title>count_pairs · one partner: 1.6M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">count_pairs · one partner</text>
<circle class="dot" cx="593.6" cy="50.0" r="5.5"/>
<text class="val" x="604.6" y="54.0">1.6M</text>
</g>
<g class="row gb purple"><title>count_pairs · 50 kb × partner: 1.6M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">count_pairs · 50 kb × partner</text>
<circle class="dot" cx="588.6" cy="78.0" r="5.5"/>
<text class="val" x="599.6" y="82.0">1.6M</text>
</g>
<g class="row gb purple"><title>count_pairs_2d · 500 kb: 1.5M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">count_pairs_2d · 500 kb</text>
<circle class="dot" cx="587.0" cy="106.0" r="5.5"/>
<text class="val" x="598.0" y="110.0">1.5M</text>
</g>
<g class="row"><title>cooler cload · 500 kb (CLI): 586k · 5,000,000 pairs</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">cooler cload · 500 kb (CLI)</text>
<circle class="dot" cx="485.8" cy="134.0" r="5.5"/>
<text class="val" x="496.8" y="138.0">586k</text>
</g>
<g class="row"><title>per-pair Python loop: 345k · 200,000 pairs</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">per-pair Python loop</text>
<circle class="dot" cx="430.2" cy="162.0" r="5.5"/>
<text class="val" x="441.2" y="166.0">345k</text>
</g>
</svg>
</div><figcaption><strong>Read pairs per second</strong> (higher is better; log scale). The naive loop was timed on the first 200k pairs.</figcaption></figure>

## Genes
{: .sec-navy #genes }

A GENCODE-shaped GTF with 20,000 genes, 70,028 transcripts and 876,947 lines.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 350" role="img" aria-label="Gene-model steps" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="324"/>
<text class="tick" x="300.0" y="340" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="324"/>
<text class="tick" x="461.3" y="340" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="324"/>
<text class="tick" x="622.7" y="340" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="324"/>
<text class="tick" x="784.0" y="340" text-anchor="middle">10 s</text>
<g class="row gb navy"><title>annotations (cgranges) · 10,000: 88 ms</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">annotations (cgranges) · 10,000</text>
<circle class="dot" cx="452.6" cy="22.0" r="5.5"/>
<text class="val" x="463.6" y="26.0">88 ms</text>
</g>
<g class="row gb navy"><title>annotations (pure-Python fallback) · 10,000: 128 ms</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">annotations (pure-Python fallback) · 10,000</text>
<circle class="dot" cx="478.6" cy="50.0" r="5.5"/>
<text class="val" x="489.6" y="54.0">128 ms</text>
</g>
<g class="row gb navy"><title>select_isoforms (peaks): 202 ms</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">select_isoforms (peaks)</text>
<circle class="dot" cx="510.6" cy="78.0" r="5.5"/>
<text class="val" x="521.6" y="82.0">202 ms</text>
</g>
<g class="row gb navy"><title>nearest_genes (pyranges) · 10,000: 293 ms</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">nearest_genes (pyranges) · 10,000</text>
<circle class="dot" cx="536.7" cy="106.0" r="5.5"/>
<text class="val" x="547.7" y="110.0">293 ms</text>
</g>
<g class="row gb navy"><title>select_isoforms (peaks + bigWig): 514 ms</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">select_isoforms (peaks + bigWig)</text>
<circle class="dot" cx="576.0" cy="134.0" r="5.5"/>
<text class="val" x="587.0" y="138.0">514 ms</text>
</g>
<g class="row gb navy"><title>nearest_genes (pyranges) · 100,000: 604 ms</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">nearest_genes (pyranges) · 100,000</text>
<circle class="dot" cx="587.3" cy="162.0" r="5.5"/>
<text class="val" x="598.3" y="166.0">604 ms</text>
</g>
<g class="row gb navy"><title>annotations (cgranges) · 100,000: 853 ms</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">annotations (cgranges) · 100,000</text>
<circle class="dot" cx="611.5" cy="190.0" r="5.5"/>
<text class="val" x="622.5" y="194.0">853 ms</text>
</g>
<g class="row gb navy"><title>annotations (pure-Python fallback) · 100,000: 1.1 s</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">annotations (pure-Python fallback) · 100,000</text>
<circle class="dot" cx="631.2" cy="218.0" r="5.5"/>
<text class="val" x="642.2" y="222.0">1.1 s</text>
</g>
<g class="row gb navy"><title>select_isoforms (bigWig only): 1.5 s</title>
<rect class="hit" x="0" y="232.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="246.0" x2="784" y2="246.0"/>
<text class="lab" x="286" y="250.0" text-anchor="end">select_isoforms (bigWig only)</text>
<circle class="dot" cx="650.7" cy="246.0" r="5.5"/>
<text class="val" x="661.7" y="250.0">1.5 s</text>
</g>
<g class="row gb navy"><title>annot index build (prom/exon/UTR merge): 2.2 s</title>
<rect class="hit" x="0" y="260.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="274.0" x2="784" y2="274.0"/>
<text class="lab" x="286" y="278.0" text-anchor="end">annot index build (prom/exon/UTR merge)</text>
<circle class="dot" cx="676.3" cy="274.0" r="5.5"/>
<text class="val" x="687.3" y="278.0">2.2 s</text>
</g>
<g class="row gb navy"><title>Genes.make (GTF parse): 4.9 s</title>
<rect class="hit" x="0" y="288.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="302.0" x2="784" y2="302.0"/>
<text class="lab" x="286" y="306.0" text-anchor="end">Genes.make (GTF parse)</text>
<circle class="dot" cx="733.5" cy="302.0" r="5.5"/>
<text class="val" x="744.5" y="306.0">4.9 s</text>
</g>
</svg>
</div><figcaption><strong>Gene-model steps</strong> (wall time; log scale). Labelling uses the cached annotation index; the counts are CREs labelled.</figcaption></figure>

## Architecture
{: .sec-purple #architecture }

The classic pipeline on 100k CREs and 50k loops (114,367 edges), step by step:

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 266" role="img" aria-label="Architecture steps" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="240"/>
<text class="tick" x="300.0" y="256" text-anchor="middle">1 ms</text>
<line class="grid" x1="396.8" y1="8" x2="396.8" y2="240"/>
<text class="tick" x="396.8" y="256" text-anchor="middle">10 ms</text>
<line class="grid" x1="493.6" y1="8" x2="493.6" y2="240"/>
<text class="tick" x="493.6" y="256" text-anchor="middle">100 ms</text>
<line class="grid" x1="590.4" y1="8" x2="590.4" y2="240"/>
<text class="tick" x="590.4" y="256" text-anchor="middle">1 s</text>
<line class="grid" x1="687.2" y1="8" x2="687.2" y2="240"/>
<text class="tick" x="687.2" y="256" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="240"/>
<text class="tick" x="784.0" y="256" text-anchor="middle">100 s</text>
<g class="row"><title>strength: graph-tool incident_edges_op: 3.0 ms</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">strength: graph-tool incident_edges_op</text>
<circle class="dot" cx="346.7" cy="22.0" r="5.5"/>
<text class="val" x="357.7" y="26.0">3.0 ms</text>
</g>
<g class="row gb purple"><title>strength (genomeblocks, numpy bincount): 9.3 ms</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">strength (genomeblocks, numpy bincount)</text>
<circle class="dot" cx="393.8" cy="50.0" r="5.5"/>
<text class="val" x="404.8" y="54.0">9.3 ms</text>
</g>
<g class="row gb purple"><title>annotate stage 2: vectorised (genomeblocks): 70 ms</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">annotate stage 2: vectorised (genomeblocks)</text>
<circle class="dot" cx="478.7" cy="78.0" r="5.5"/>
<text class="val" x="489.7" y="82.0">70 ms</text>
</g>
<g class="row gb purple"><title>normalize (power-law O/E): 219 ms</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">normalize (power-law O/E)</text>
<circle class="dot" cx="526.6" cy="106.0" r="5.5"/>
<text class="val" x="537.6" y="110.0">219 ms</text>
</g>
<g class="row gb purple"><title>make (50k loops → graph): 1.8 s</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">make (50k loops → graph)</text>
<circle class="dot" cx="614.7" cy="134.0" r="5.5"/>
<text class="val" x="625.7" y="138.0">1.8 s</text>
</g>
<g class="row gb purple"><title>annotate (total): 1.8 s</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">annotate (total)</text>
<circle class="dot" cx="616.0" cy="162.0" r="5.5"/>
<text class="val" x="627.0" y="166.0">1.8 s</text>
</g>
<g class="row gb purple"><title>add_mcool (5 kb): 9.5 s</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">add_mcool (5 kb)</text>
<circle class="dot" cx="685.1" cy="190.0" r="5.5"/>
<text class="val" x="696.1" y="194.0">9.5 s</text>
</g>
<g class="row"><title>annotate stage 2: per-vertex loop: 23 s</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">annotate stage 2: per-vertex loop</text>
<circle class="dot" cx="721.9" cy="218.0" r="5.5"/>
<text class="val" x="732.9" y="222.0">23 s</text>
</g>
</svg>
</div><figcaption><strong>Architecture pipeline steps</strong>, classic module (wall time; log scale). The grey rows are the per-vertex loop that annotate's vectorised stage replaced and graph-tool's own reduction for comparison.</figcaption></figure>

## Columnar vs classic
{: .sec-purple #columnar }

The same pipeline (load CREs and genes, build, add Hi-C, normalise, annotate, strength, hubs, save) in a fresh process: **26 s** with the classic modules, **6.2 s** with `genomeblocks.columnar`, peak memory 1,055 MB → 893 MB. Both give the same edges, weights, O/E, labels and hubs.

Reloading the saved graph in a new process takes 9.6 s from the classic pickle and 1.3 s from the columnar parquet tables (16 MB vs 3.2 MB on disk).

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 400" role="img" aria-label="Classic vs columnar per step" xmlns="http://www.w3.org/2000/svg">
<circle class="dot a" cx="256" cy="12" r="5"/><text class="lab" x="266" y="16">classic</text>
<g class="gb purple"><circle class="dot" cx="360" cy="12" r="5"/></g><text class="lab" x="370" y="16">columnar</text>
<line class="grid" x1="250.0" y1="30" x2="250.0" y2="374"/>
<text class="tick" x="250.0" y="390" text-anchor="middle">100 µs</text>
<line class="grid" x1="362.0" y1="30" x2="362.0" y2="374"/>
<text class="tick" x="362.0" y="390" text-anchor="middle">1 ms</text>
<line class="grid" x1="474.0" y1="30" x2="474.0" y2="374"/>
<text class="tick" x="474.0" y="390" text-anchor="middle">10 ms</text>
<line class="grid" x1="586.0" y1="30" x2="586.0" y2="374"/>
<text class="tick" x="586.0" y="390" text-anchor="middle">100 ms</text>
<line class="grid" x1="698.0" y1="30" x2="698.0" y2="374"/>
<text class="tick" x="698.0" y="390" text-anchor="middle">1 s</text>
<line class="grid" x1="810.0" y1="30" x2="810.0" y2="374"/>
<text class="tick" x="810.0" y="390" text-anchor="middle">10 s</text>
<g class="row gb purple"><title>Loci.make (100k CREs): classic 173 ms → columnar 42 ms (4.2× faster)</title>
<rect class="hit" x="0" y="30.0" width="880" height="28"/>
<text class="lab" x="236" y="48.0" text-anchor="end">Loci.make (100k CREs)</text>
<line class="span" x1="543.2" y1="44.0" x2="612.5" y2="44.0"/>
<circle class="dot a" cx="612.5" cy="44.0" r="5"/>
<circle class="dot" cx="543.2" cy="44.0" r="5.5"/>
<text class="val" x="874" y="48.0" text-anchor="end">4.2×</text>
</g>
<g class="row gb purple"><title>Genes.make (GTF, 877k lines): classic 4.9 s → columnar 536 ms (9.2× faster)</title>
<rect class="hit" x="0" y="58.0" width="880" height="28"/>
<text class="lab" x="236" y="76.0" text-anchor="end">Genes.make (GTF, 877k lines)</text>
<line class="span" x1="667.7" y1="72.0" x2="775.4" y2="72.0"/>
<circle class="dot a" cx="775.4" cy="72.0" r="5"/>
<circle class="dot" cx="667.7" cy="72.0" r="5.5"/>
<text class="val" x="874" y="76.0" text-anchor="end">9.2×</text>
</g>
<g class="row gb purple"><title>gene annotation index: classic 2.2 s → columnar 109 ms (20× faster)</title>
<rect class="hit" x="0" y="86.0" width="880" height="28"/>
<text class="lab" x="236" y="104.0" text-anchor="end">gene annotation index</text>
<line class="span" x1="590.1" y1="100.0" x2="735.7" y2="100.0"/>
<circle class="dot a" cx="735.7" cy="100.0" r="5"/>
<circle class="dot" cx="590.1" cy="100.0" r="5.5"/>
<text class="val" x="874" y="104.0" text-anchor="end">20×</text>
</g>
<g class="row gb purple"><title>Genes.annotations (100k CREs): classic 833 ms → columnar 25 ms (34× faster)</title>
<rect class="hit" x="0" y="114.0" width="880" height="28"/>
<text class="lab" x="236" y="132.0" text-anchor="end">Genes.annotations (100k CREs)</text>
<line class="span" x1="517.8" y1="128.0" x2="689.1" y2="128.0"/>
<circle class="dot a" cx="689.1" cy="128.0" r="5"/>
<circle class="dot" cx="517.8" cy="128.0" r="5.5"/>
<text class="val" x="874" y="132.0" text-anchor="end">34×</text>
</g>
<g class="row gb purple"><title>Architecture.make (50k loops): classic 1.6 s → columnar 99 ms (16× faster)</title>
<rect class="hit" x="0" y="142.0" width="880" height="28"/>
<text class="lab" x="236" y="160.0" text-anchor="end">Architecture.make (50k loops)</text>
<line class="span" x1="585.6" y1="156.0" x2="719.5" y2="156.0"/>
<circle class="dot a" cx="719.5" cy="156.0" r="5"/>
<circle class="dot" cx="585.6" cy="156.0" r="5.5"/>
<text class="val" x="874" y="160.0" text-anchor="end">16×</text>
</g>
<g class="row gb purple"><title>add_mcool (5 kb): classic 5.7 s → columnar 244 ms (23× faster)</title>
<rect class="hit" x="0" y="170.0" width="880" height="28"/>
<text class="lab" x="236" y="188.0" text-anchor="end">add_mcool (5 kb)</text>
<line class="span" x1="629.4" y1="184.0" x2="782.8" y2="184.0"/>
<circle class="dot a" cx="782.8" cy="184.0" r="5"/>
<circle class="dot" cx="629.4" cy="184.0" r="5.5"/>
<text class="val" x="874" y="188.0" text-anchor="end">23×</text>
</g>
<g class="row gb purple"><title>normalize (power-law O/E): classic 221 ms → columnar 8.5 ms (26× faster)</title>
<rect class="hit" x="0" y="198.0" width="880" height="28"/>
<text class="lab" x="236" y="216.0" text-anchor="end">normalize (power-law O/E)</text>
<line class="span" x1="465.9" y1="212.0" x2="624.5" y2="212.0"/>
<circle class="dot a" cx="624.5" cy="212.0" r="5"/>
<circle class="dot" cx="465.9" cy="212.0" r="5.5"/>
<text class="val" x="874" y="216.0" text-anchor="end">26×</text>
</g>
<g class="row gb purple"><title>annotate (labels + genes): classic 1.8 s → columnar 35 ms (51× faster)</title>
<rect class="hit" x="0" y="226.0" width="880" height="28"/>
<text class="lab" x="236" y="244.0" text-anchor="end">annotate (labels + genes)</text>
<line class="span" x1="535.5" y1="240.0" x2="726.3" y2="240.0"/>
<circle class="dot a" cx="726.3" cy="240.0" r="5"/>
<circle class="dot" cx="535.5" cy="240.0" r="5.5"/>
<text class="val" x="874" y="244.0" text-anchor="end">51×</text>
</g>
<g class="row gb purple"><title>strength: classic 9.0 ms → columnar 1.0 ms (8.9× faster)</title>
<rect class="hit" x="0" y="254.0" width="880" height="28"/>
<text class="lab" x="236" y="272.0" text-anchor="end">strength</text>
<line class="span" x1="362.7" y1="268.0" x2="469.0" y2="268.0"/>
<circle class="dot a" cx="469.0" cy="268.0" r="5"/>
<circle class="dot" cx="362.7" cy="268.0" r="5.5"/>
<text class="val" x="874" y="272.0" text-anchor="end">8.9×</text>
</g>
<g class="row gb purple"><title>prime_hubs: classic 109 ms → columnar 8.5 ms (13× faster)</title>
<rect class="hit" x="0" y="282.0" width="880" height="28"/>
<text class="lab" x="236" y="300.0" text-anchor="end">prime_hubs</text>
<line class="span" x1="466.0" y1="296.0" x2="590.0" y2="296.0"/>
<circle class="dot a" cx="590.0" cy="296.0" r="5"/>
<circle class="dot" cx="466.0" cy="296.0" r="5.5"/>
<text class="val" x="874" y="300.0" text-anchor="end">13×</text>
</g>
<g class="row gb purple"><title>save: classic 1.3 s → columnar 56 ms (23× faster)</title>
<rect class="hit" x="0" y="310.0" width="880" height="28"/>
<text class="lab" x="236" y="328.0" text-anchor="end">save</text>
<line class="span" x1="557.9" y1="324.0" x2="710.5" y2="324.0"/>
<circle class="dot a" cx="710.5" cy="324.0" r="5"/>
<circle class="dot" cx="557.9" cy="324.0" r="5.5"/>
<text class="val" x="874" y="328.0" text-anchor="end">23×</text>
</g>
<g class="row gb purple"><title>load: classic 1.8 s → columnar 58 ms (31× faster)</title>
<rect class="hit" x="0" y="338.0" width="880" height="28"/>
<text class="lab" x="236" y="356.0" text-anchor="end">load</text>
<line class="span" x1="559.7" y1="352.0" x2="727.1" y2="352.0"/>
<circle class="dot a" cx="727.1" cy="352.0" r="5"/>
<circle class="dot" cx="559.7" cy="352.0" r="5.5"/>
<text class="val" x="874" y="356.0" text-anchor="end">31×</text>
</g>
</svg>
</div><figcaption><strong>Per step, classic (grey) vs columnar</strong> on 100k CREs and 50k loops (in-process medians; log scale; right column = speed-up).</figcaption></figure>

<details class="gb-table" open><summary>Interactive operations</summary><table><thead><tr><th>operation</th><th>classic</th><th>columnar</th><th>speed-up</th></tr></thead><tbody><tr><td>one chromosome (chr8)</td><td>303 ms</td><td>32 µs</td><td>9,374×</td></tr><tr><td>all cis / all trans split</td><td>–</td><td>39 µs</td><td>–</td></tr><tr><td>region chr8:20-40 Mb</td><td>192 ms</td><td>998 µs</td><td>193×</td></tr><tr><td>neighbours of 1,000 CREs</td><td>321 ms</td><td>4.9 ms</td><td>66×</td></tr><tr><td>copy</td><td>2.8 s</td><td>2.2 ms</td><td>1,312×</td></tr><tr><td>graph-tool Graph from tables</td><td>–</td><td>14 ms</td><td>–</td></tr><tr><td>connected components (graph-tool)</td><td>10 ms</td><td>9.0 ms</td><td>1.1×</td></tr><tr><td>pagerank (graph-tool)</td><td>5.1 ms</td><td>3.9 ms</td><td>1.3×</td></tr></tbody></table></details>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 300" role="img" aria-label="Pipeline time vs loop count" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="64" y1="256.0" x2="650" y2="256.0"/>
<text class="tick" x="56" y="260.0" text-anchor="end">100 ms</text>
<line class="grid" x1="64" y1="175.3" x2="650" y2="175.3"/>
<text class="tick" x="56" y="179.3" text-anchor="end">1 s</text>
<line class="grid" x1="64" y1="94.7" x2="650" y2="94.7"/>
<text class="tick" x="56" y="98.7" text-anchor="end">10 s</text>
<line class="grid" x1="64" y1="14.0" x2="650" y2="14.0"/>
<text class="tick" x="56" y="18.0" text-anchor="end">100 s</text>
<text class="tick" x="64.0" y="274" text-anchor="middle">10k</text>
<text class="tick" x="357.0" y="274" text-anchor="middle">100k</text>
<text class="tick" x="650.0" y="274" text-anchor="middle">1M</text>
<line class="axis" x1="64" y1="256" x2="650" y2="256"/>
<text class="axl" x="357.0" y="294" text-anchor="middle">loops</text>
<g class="line"><title>classic</title><polyline points="92.4,111.6 268.8,103.5 445.2,86.2"/>
<circle cx="92.4" cy="111.6" r="4"><title>classic · 12k: 6.2 s</title></circle>
<circle cx="268.8" cy="103.5" r="4"><title>classic · 50k: 7.8 s</title></circle>
<circle cx="445.2" cy="86.2" r="4"><title>classic · 200k: 13 s</title></circle>
</g>
<g class="line gb purple"><title>columnar</title><polyline points="92.4,217.4 268.8,208.4 445.2,182.3 650.0,140.6"/>
<circle cx="92.4" cy="217.4" r="4"><title>columnar · 12k: 301 ms</title></circle>
<circle cx="268.8" cy="208.4" r="4"><title>columnar · 50k: 390 ms</title></circle>
<circle cx="445.2" cy="182.3" r="4"><title>columnar · 200k: 820 ms</title></circle>
<circle cx="650.0" cy="140.6" r="4"><title>columnar · 1M: 2.7 s</title></circle>
</g>
<text class="lab" x="660" y="90.2">classic</text>
<line class="leader" x1="451.2" y1="86.2" x2="656" y2="86.2"/>
<text class="lab gbl" x="660" y="144.6">columnar</text>
<line class="leader" x1="656.0" y1="140.6" x2="656" y2="140.6"/>
</svg>
</div><figcaption><strong>Pipeline time vs loop count</strong> (make → add_mcool → normalize → annotate → strength; log–log).</figcaption></figure>

## HiChIP short-range track
{: .sec-purple #shortrange }

From allValidPairs to a coverage bigWig: **17 s** and 1.0 GB peak memory with `columnar.hichip`, against 101 s and 1.9 GB for the `awk | sort | bedtools` recipe, with the same ends and a byte-identical bedGraph.

<details class="gb-table"><summary>Steps</summary><table><thead><tr><th>implementation</th><th>step</th><th>time</th><th>peak memory</th></tr></thead><tbody><tr><td>shell</td><td>extract ends (awk)</td><td>5.5 s</td><td>20 MB</td></tr><tr><td>shell</td><td>fragments + sort + genomecov</td><td>88 s</td><td>1,903 MB</td></tr><tr><td>shell</td><td>bedGraph -&gt; bigWig</td><td>7.8 s</td><td>–</td></tr><tr><td>genomeblocks</td><td>ends</td><td>5.3 s</td><td>–</td></tr><tr><td>genomeblocks</td><td>write BED</td><td>665 ms</td><td>–</td></tr><tr><td>genomeblocks</td><td>coverage (bedGraph)</td><td>3.6 s</td><td>–</td></tr><tr><td>genomeblocks</td><td>bigWig</td><td>7.6 s</td><td>–</td></tr></tbody></table></details>

## Import cost
{: .sec-navy #import }

`import genomeblocks` is lazy, and `from genomeblocks import Loci` costs **214 ms** on top of starting Python (35 ms): matplotlib, pandas, scipy and graph-tool load only when a function needs them.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 238" role="img" aria-label="Import cost" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="212"/>
<text class="tick" x="300.0" y="228" text-anchor="middle">1 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="212"/>
<text class="tick" x="421.0" y="228" text-anchor="middle">10 ms</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="212"/>
<text class="tick" x="542.0" y="228" text-anchor="middle">100 ms</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="212"/>
<text class="tick" x="663.0" y="228" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="212"/>
<text class="tick" x="784.0" y="228" text-anchor="middle">10 s</text>
<g class="row gb navy"><title>import genomeblocks: 15 ms</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">import genomeblocks</text>
<circle class="dot" cx="443.7" cy="22.0" r="5.5"/>
<text class="val" x="454.7" y="26.0">15 ms</text>
</g>
<g class="row gb navy"><title>+ Locus: 26 ms</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">+ Locus</text>
<circle class="dot" cx="471.6" cy="50.0" r="5.5"/>
<text class="val" x="482.6" y="54.0">26 ms</text>
</g>
<g class="row gb navy"><title>+ Genes: 34 ms</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">+ Genes</text>
<circle class="dot" cx="485.1" cy="78.0" r="5.5"/>
<text class="val" x="496.1" y="82.0">34 ms</text>
</g>
<g class="row gb navy"><title>+ Loci (pulls signal, motifs, atlas, bedpe): 214 ms</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">+ Loci (pulls signal, motifs, atlas, bedpe)</text>
<circle class="dot" cx="582.0" cy="106.0" r="5.5"/>
<text class="val" x="593.0" y="110.0">214 ms</text>
</g>
<g class="row gb navy"><title>+ browser (matplotlib): 583 ms</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">+ browser (matplotlib)</text>
<circle class="dot" cx="634.6" cy="134.0" r="5.5"/>
<text class="val" x="645.6" y="138.0">583 ms</text>
</g>
<g class="row gb navy"><title>+ Architecture (graph-tool): 1.7 s</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">+ Architecture (graph-tool)</text>
<circle class="dot" cx="691.6" cy="162.0" r="5.5"/>
<text class="val" x="702.6" y="166.0">1.7 s</text>
</g>
<g class="row gb navy"><title>everything: 1.7 s</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">everything</text>
<circle class="dot" cx="692.4" cy="190.0" r="5.5"/>
<text class="val" x="703.4" y="194.0">1.7 s</text>
</g>
</svg>
</div><figcaption><strong>Import cost</strong> above <code>python -c pass</code>, fresh process each (cumulative statements; log scale).</figcaption></figure>

## Reproduce
{: .sec-navy #reproduce }

```bash
cd benchmarks
python make_data.py            # synthetic hg38-shaped data, ~1.5 GB
PY=python ./run_all.sh         # every bench, one at a time
python build_site.py           # this page: docs/benchmarks/index.md
```

External baselines are found on `PATH` (bedtools, deepTools `computeMatrix`,
cooler) or via `FIMO=`. The data are synthetic and seeded, laid out on real hg38
chromosome sizes; motifs are the JASPAR CORE vertebrate collection.

<p class="gb-env">Versions: numpy 2.3.4, scipy 1.16.3, pandas 2.3.3, pybigtools 0.2.5, pyBigWig 0.3.24, cgranges installed, pyranges 0.1.4, bioframe 0.8.0, lightmotif 0.10.0, Bio 1.86, MOODS installed, graph_tool 2.98 (commit c96a6bf3, ), intervaltree installed, bedtools v2.31.1</p>
