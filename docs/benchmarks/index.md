---
title: Benchmarks
layout: default
nav_order: 9
permalink: /benchmarks/
---

# Benchmarks
{: .no_toc }

How fast each building block is against the tools people would otherwise use,
and what each swappable backend costs, measured on genomeblocks 2.0.0.
Every comparison first checks that the engines return the same answer;
timings are medians of 3–5 runs after a warm-up.
{: .fs-5 .fw-300 }

<p class="gb-env">Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz · 64 cores · 811 GB RAM ·
Python 3.12.12 · synthetic hg38-shaped data (seeded) · benchmarks ran one at a time and use at
most 8 worker processes</p>

<div class="gb-tiles">
<div class="gb-tile green"><span class="tile-label">intervals: fastest engine</span><span class="tile-value">genomeblocks</span><span class="tile-ctx">merge</span></div>
<div class="gb-tile green"><span class="tile-label">bigwig: fastest engine</span><span class="tile-value">pybigtools</span><span class="tile-ctx">signal 10,000 x 1 track x 200 bins</span></div>
<div class="gb-tile navy"><span class="tile-label">motifs: fastest engine</span><span class="tile-value">lightmotif</span><span class="tile-ctx">scan 1,000 windows x 100 motifs</span></div>
<div class="gb-tile navy"><span class="tile-label">fasta: fastest engine</span><span class="tile-value">memory</span><span class="tile-ctx">sequences 10,000 x 500 bp</span></div>
<div class="gb-tile green"><span class="tile-label">tables: fastest engine</span><span class="tile-value">polars</span><span class="tile-ctx">Loci.make 1,000,000 peaks (keep=True)</span></div>
<div class="gb-tile purple"><span class="tile-label">graph: fastest engine</span><span class="tile-value">graph-tool</span><span class="tile-ctx">build the engine's graph</span></div>
<div class="gb-tile green"><span class="tile-label">One overlap lookup</span><span class="tile-value">8.8 µs</span><span class="tile-ctx">on 100k indexed peaks: 167× faster than pyranges, 607× than bioframe</span></div>
<div class="gb-tile green"><span class="tile-label">bigWig → heatmap matrix</span><span class="tile-value">98×</span><span class="tile-ctx">faster than deepTools computeMatrix on one core (13k vs 137 regions/s)</span></div>
<div class="gb-tile navy"><span class="tile-label">Enrichment vs 500 peak files</span><span class="tile-value">11 ms</span><span class="tile-ctx">per 20k-peak query including Fisher tests</span></div>
<div class="gb-tile navy"><span class="tile-label">Motif scanning</span><span class="tile-value">131×</span><span class="tile-ctx">faster than MEME FIMO, 187× faster than Biopython</span></div>
<div class="gb-tile purple"><span class="tile-label">Hi-C pair counting</span><span class="tile-value">1.6M/s</span><span class="tile-ctx">read pairs into 50 kb windows: 2.9× the rate of cooler cload</span></div>
</div>

1. TOC
{:toc}

## Backends
{: .sec-navy #backends }

The same call through every installed engine (`backend=`); every row below returns the same answer as the default, so the choice is only about speed and what is installed.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 210" role="img" aria-label="intervals: A &amp; B" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="184"/>
<text class="tick" x="300.0" y="200" text-anchor="middle">100 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="184"/>
<text class="tick" x="461.3" y="200" text-anchor="middle">1 s</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="184"/>
<text class="tick" x="622.7" y="200" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="184"/>
<text class="tick" x="784.0" y="200" text-anchor="middle">100 s</text>
<g class="row"><title>ncls: 288 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">ncls</text>
<circle class="dot" cx="374.0" cy="22.0" r="5.5"/>
<text class="val" x="385.0" y="26.0">288 ms</text>
</g>
<g class="row gb green"><title>genomeblocks: 291 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="374.7" cy="50.0" r="5.5"/>
<text class="val" x="385.7" y="54.0">291 ms</text>
</g>
<g class="row"><title>bioframe: 966 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="458.9" cy="78.0" r="5.5"/>
<text class="val" x="469.9" y="82.0">966 ms</text>
</g>
<g class="row"><title>pyranges: 1.3 s · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="479.3" cy="106.0" r="5.5"/>
<text class="val" x="490.3" y="110.0">1.3 s</text>
</g>
<g class="row"><title>cgranges: 1.4 s · agrees</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">cgranges</text>
<circle class="dot" cx="482.7" cy="134.0" r="5.5"/>
<text class="val" x="493.7" y="138.0">1.4 s</text>
</g>
<g class="row"><title>bedtools: 9.3 s · agrees</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">bedtools</text>
<circle class="dot" cx="617.8" cy="162.0" r="5.5"/>
<text class="val" x="628.8" y="166.0">9.3 s</text>
</g>
</svg>
</div><figcaption><strong>intervals · A &amp; B</strong> at n = 1,000,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="intervals: merge" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="128"/>
<text class="tick" x="461.3" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="128"/>
<text class="tick" x="622.7" y="144" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">10 s</text>
<g class="row gb green"><title>genomeblocks: 86 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="450.8" cy="22.0" r="5.5"/>
<text class="val" x="461.8" y="26.0">86 ms</text>
</g>
<g class="row"><title>pyranges: 574 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="583.7" cy="50.0" r="5.5"/>
<text class="val" x="594.7" y="54.0">574 ms</text>
</g>
<g class="row"><title>bioframe: 661 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="593.7" cy="78.0" r="5.5"/>
<text class="val" x="604.7" y="82.0">661 ms</text>
</g>
<g class="row"><title>bedtools: 4.9 s · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">bedtools</text>
<circle class="dot" cx="733.3" cy="106.0" r="5.5"/>
<text class="val" x="744.3" y="110.0">4.9 s</text>
</g>
</svg>
</div><figcaption><strong>intervals · merge</strong> at n = 1,000,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="intervals: nearest (distances)" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="128"/>
<text class="tick" x="461.3" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="128"/>
<text class="tick" x="622.7" y="144" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">10 s</text>
<g class="row gb green"><title>genomeblocks: 30 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="377.1" cy="22.0" r="5.5"/>
<text class="val" x="388.1" y="26.0">30 ms</text>
</g>
<g class="row"><title>bioframe: 228 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="519.2" cy="50.0" r="5.5"/>
<text class="val" x="530.2" y="54.0">228 ms</text>
</g>
<g class="row"><title>pyranges: 392 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="557.0" cy="78.0" r="5.5"/>
<text class="val" x="568.0" y="82.0">392 ms</text>
</g>
<g class="row"><title>bedtools: 1.8 s · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">bedtools</text>
<circle class="dot" cx="663.3" cy="106.0" r="5.5"/>
<text class="val" x="674.3" y="110.0">1.8 s</text>
</g>
</svg>
</div><figcaption><strong>intervals · nearest (distances)</strong> at n = 100,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 126" role="img" aria-label="bigwig: signal 10,000 x 1 track x 200 bins" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="100"/>
<text class="tick" x="300.0" y="116" text-anchor="middle">100 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="100"/>
<text class="tick" x="421.0" y="116" text-anchor="middle">1 s</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="100"/>
<text class="tick" x="542.0" y="116" text-anchor="middle">10 s</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="100"/>
<text class="tick" x="663.0" y="116" text-anchor="middle">100 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="100"/>
<text class="tick" x="784.0" y="116" text-anchor="middle">1,000 s</text>
<g class="row gb green"><title>pybigtools: 455 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">pybigtools</text>
<circle class="dot" cx="379.6" cy="22.0" r="5.5"/>
<text class="val" x="390.6" y="26.0">455 ms</text>
</g>
<g class="row"><title>python: 3.1 s · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">python</text>
<circle class="dot" cx="480.2" cy="50.0" r="5.5"/>
<text class="val" x="491.2" y="54.0">3.1 s</text>
</g>
<g class="row"><title>pybigwig: 4.9 min · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">pybigwig</text>
<circle class="dot" cx="720.0" cy="78.0" r="5.5"/>
<text class="val" x="731.0" y="82.0">4.9 min</text>
</g>
</svg>
</div><figcaption><strong>bigwig · signal 10,000 x 1 track x 200 bins</strong> at n = 10,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 126" role="img" aria-label="motifs: scan 1,000 windows x 100 motifs" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="100"/>
<text class="tick" x="300.0" y="116" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="100"/>
<text class="tick" x="461.3" y="116" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="100"/>
<text class="tick" x="622.7" y="116" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="100"/>
<text class="tick" x="784.0" y="116" text-anchor="middle">10 s</text>
<g class="row"><title>lightmotif: 37 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">lightmotif</text>
<circle class="dot" cx="391.9" cy="22.0" r="5.5"/>
<text class="val" x="402.9" y="26.0">37 ms</text>
</g>
<g class="row gb navy"><title>moods: 53 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">moods</text>
<circle class="dot" cx="416.6" cy="50.0" r="5.5"/>
<text class="val" x="427.6" y="54.0">53 ms</text>
</g>
<g class="row"><title>biopython: 3.3 s · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">biopython</text>
<circle class="dot" cx="706.3" cy="78.0" r="5.5"/>
<text class="val" x="717.3" y="82.0">3.3 s</text>
</g>
</svg>
</div><figcaption><strong>motifs · scan 1,000 windows x 100 motifs</strong> at n = 1,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 182" role="img" aria-label="fasta: sequences 10,000 x 500 bp" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="156"/>
<text class="tick" x="300.0" y="172" text-anchor="middle">1 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="156"/>
<text class="tick" x="421.0" y="172" text-anchor="middle">10 ms</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="156"/>
<text class="tick" x="542.0" y="172" text-anchor="middle">100 ms</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="156"/>
<text class="tick" x="663.0" y="172" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="156"/>
<text class="tick" x="784.0" y="172" text-anchor="middle">10 s</text>
<g class="row"><title>memory: 5.4 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">memory</text>
<circle class="dot" cx="388.7" cy="22.0" r="5.5"/>
<text class="val" x="399.7" y="26.0">5.4 ms</text>
</g>
<g class="row gb navy"><title>genomeblocks: 37 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="489.9" cy="50.0" r="5.5"/>
<text class="val" x="500.9" y="54.0">37 ms</text>
</g>
<g class="row"><title>pysam: 40 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">pysam</text>
<circle class="dot" cx="494.0" cy="78.0" r="5.5"/>
<text class="val" x="505.0" y="82.0">40 ms</text>
</g>
<g class="row"><title>pyfaidx: 135 ms · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyfaidx</text>
<circle class="dot" cx="557.6" cy="106.0" r="5.5"/>
<text class="val" x="568.6" y="110.0">135 ms</text>
</g>
<g class="row"><title>biopython: 5.0 s · agrees</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">biopython</text>
<circle class="dot" cx="747.4" cy="134.0" r="5.5"/>
<text class="val" x="758.4" y="138.0">5.0 s</text>
</g>
</svg>
</div><figcaption><strong>fasta · sequences 10,000 x 500 bp</strong> at n = 10,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 98" role="img" aria-label="tables: Loci.make 1,000,000 peaks (keep=True)" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="72"/>
<text class="tick" x="300.0" y="88" text-anchor="middle">100 ms</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="72"/>
<text class="tick" x="542.0" y="88" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="72"/>
<text class="tick" x="784.0" y="88" text-anchor="middle">10 s</text>
<g class="row gb green"><title>polars: 462 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">polars</text>
<circle class="dot" cx="460.8" cy="22.0" r="5.5"/>
<text class="val" x="471.8" y="26.0">462 ms</text>
</g>
<g class="row"><title>pandas: 1.2 s · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pandas</text>
<circle class="dot" cx="558.4" cy="50.0" r="5.5"/>
<text class="val" x="569.4" y="54.0">1.2 s</text>
</g>
</svg>
</div><figcaption><strong>tables · Loci.make 1,000,000 peaks (keep=True)</strong> at n = 1,000,000; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 98" role="img" aria-label="tables: Genes.make (GTF)" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="72"/>
<text class="tick" x="300.0" y="88" text-anchor="middle">100 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="72"/>
<text class="tick" x="461.3" y="88" text-anchor="middle">1 s</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="72"/>
<text class="tick" x="622.7" y="88" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="72"/>
<text class="tick" x="784.0" y="88" text-anchor="middle">100 s</text>
<g class="row gb green"><title>polars: 550 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">polars</text>
<circle class="dot" cx="419.4" cy="22.0" r="5.5"/>
<text class="val" x="430.4" y="26.0">550 ms</text>
</g>
<g class="row"><title>pandas: 9.4 s · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pandas</text>
<circle class="dot" cx="618.4" cy="50.0" r="5.5"/>
<text class="val" x="629.4" y="54.0">9.4 s</text>
</g>
</svg>
</div><figcaption><strong>tables · Genes.make (GTF)</strong> at n = 70,028; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="graph: components (114,367 edges)" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">1 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="128"/>
<text class="tick" x="461.3" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="128"/>
<text class="tick" x="622.7" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">1 s</text>
<g class="row gb purple"><title>graph-tool: 14 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">graph-tool</text>
<circle class="dot" cx="484.4" cy="22.0" r="5.5"/>
<text class="val" x="495.4" y="26.0">14 ms</text>
</g>
<g class="row"><title>scipy: 18 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">scipy</text>
<circle class="dot" cx="501.4" cy="50.0" r="5.5"/>
<text class="val" x="512.4" y="54.0">18 ms</text>
</g>
<g class="row"><title>igraph: 59 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">igraph</text>
<circle class="dot" cx="586.0" cy="78.0" r="5.5"/>
<text class="val" x="597.0" y="82.0">59 ms</text>
</g>
<g class="row"><title>networkx: 447 ms · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">networkx</text>
<circle class="dot" cx="727.6" cy="106.0" r="5.5"/>
<text class="val" x="738.6" y="110.0">447 ms</text>
</g>
</svg>
</div><figcaption><strong>graph · components (114,367 edges)</strong> at n = 114,367; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="graph: pagerank (114,367 edges)" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="128"/>
<text class="tick" x="461.3" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="128"/>
<text class="tick" x="622.7" y="144" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">10 s</text>
<g class="row"><title>scipy: 106 ms · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">scipy</text>
<circle class="dot" cx="465.6" cy="22.0" r="5.5"/>
<text class="val" x="476.6" y="26.0">106 ms</text>
</g>
<g class="row gb purple"><title>graph-tool: 112 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">graph-tool</text>
<circle class="dot" cx="469.3" cy="50.0" r="5.5"/>
<text class="val" x="480.3" y="54.0">112 ms</text>
</g>
<g class="row"><title>igraph: 148 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">igraph</text>
<circle class="dot" cx="488.9" cy="78.0" r="5.5"/>
<text class="val" x="499.9" y="82.0">148 ms</text>
</g>
<g class="row"><title>networkx: 686 ms · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">networkx</text>
<circle class="dot" cx="596.3" cy="106.0" r="5.5"/>
<text class="val" x="607.3" y="110.0">686 ms</text>
</g>
</svg>
</div><figcaption><strong>graph · pagerank (114,367 edges)</strong> at n = 114,367; the default engine in colour.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="graph: build the engine&#x27;s graph" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">100 µs</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="128"/>
<text class="tick" x="421.0" y="144" text-anchor="middle">1 ms</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="128"/>
<text class="tick" x="542.0" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="128"/>
<text class="tick" x="663.0" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">1 s</text>
<g class="row gb purple"><title>graph-tool: 275 µs · agrees</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">graph-tool</text>
<circle class="dot" cx="353.1" cy="22.0" r="5.5"/>
<text class="val" x="364.1" y="26.0">275 µs</text>
</g>
<g class="row"><title>scipy: 6.6 ms · agrees</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">scipy</text>
<circle class="dot" cx="520.5" cy="50.0" r="5.5"/>
<text class="val" x="531.5" y="54.0">6.6 ms</text>
</g>
<g class="row"><title>igraph: 36 ms · agrees</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">igraph</text>
<circle class="dot" cx="609.6" cy="78.0" r="5.5"/>
<text class="val" x="620.6" y="82.0">36 ms</text>
</g>
<g class="row"><title>networkx: 338 ms · agrees</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">networkx</text>
<circle class="dot" cx="727.0" cy="106.0" r="5.5"/>
<text class="val" x="738.0" y="110.0">338 ms</text>
</g>
</svg>
</div><figcaption><strong>graph · build the engine&#x27;s graph</strong> at n = 114,367; the default engine in colour.</figcaption></figure>

## Loci
{: .sec-green #loci }

A single overlap query against an indexed set is the cost of every interactive question ("what overlaps this peak?"): **8.8 µs** per lookup through the cached point index (8.8 µs with genomeblocks (cgranges), 9.6 µs with genomeblocks (ncls), 14 µs with genomeblocks (numpy point index)), against 1.5 ms for pyranges, which has to slice a dataframe per call.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 210" role="img" aria-label="Single overlap query latency" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="184"/>
<text class="tick" x="300.0" y="200" text-anchor="middle">1 µs</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="184"/>
<text class="tick" x="421.0" y="200" text-anchor="middle">10 µs</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="184"/>
<text class="tick" x="542.0" y="200" text-anchor="middle">100 µs</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="184"/>
<text class="tick" x="663.0" y="200" text-anchor="middle">1 ms</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="184"/>
<text class="tick" x="784.0" y="200" text-anchor="middle">10 ms</text>
<g class="row gb green"><title>genomeblocks (cgranges): 8.8 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks (cgranges)</text>
<circle class="dot" cx="414.6" cy="22.0" r="5.5"/>
<text class="val" x="425.6" y="26.0">8.8 µs</text>
</g>
<g class="row gb green"><title>genomeblocks (ncls): 9.6 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks (ncls)</text>
<circle class="dot" cx="418.8" cy="50.0" r="5.5"/>
<text class="val" x="429.8" y="54.0">9.6 µs</text>
</g>
<g class="row"><title>intervaltree: 10 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">intervaltree</text>
<circle class="dot" cx="422.4" cy="78.0" r="5.5"/>
<text class="val" x="433.4" y="82.0">10 µs</text>
</g>
<g class="row gb green"><title>genomeblocks (numpy point index): 14 µs · per query, 100k-peak index</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks (numpy point index)</text>
<circle class="dot" cx="438.8" cy="106.0" r="5.5"/>
<text class="val" x="449.8" y="110.0">14 µs</text>
</g>
<g class="row"><title>pyranges: 1.5 ms · per query, 100k-peak index</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="683.5" cy="134.0" r="5.5"/>
<text class="val" x="694.5" y="138.0">1.5 ms</text>
</g>
<g class="row"><title>bioframe: 5.4 ms · per query, 100k-peak index</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="751.3" cy="162.0" r="5.5"/>
<text class="val" x="762.3" y="166.0">5.4 ms</text>
</g>
</svg>
</div><figcaption><strong>One overlap lookup</strong> against 100,000 indexed peaks (median per query; log scale).</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 182" role="img" aria-label="A &amp; B at 1,000,000 intervals" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="156"/>
<text class="tick" x="300.0" y="172" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="156"/>
<text class="tick" x="461.3" y="172" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="156"/>
<text class="tick" x="622.7" y="172" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="156"/>
<text class="tick" x="784.0" y="172" text-anchor="middle">10 s</text>
<g class="row gb green"><title>genomeblocks: 87 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="451.6" cy="22.0" r="5.5"/>
<text class="val" x="462.6" y="26.0">87 ms</text>
</g>
<g class="row"><title>pyranges: 336 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="546.2" cy="50.0" r="5.5"/>
<text class="val" x="557.2" y="54.0">336 ms</text>
</g>
<g class="row"><title>bioframe: 1.0 s · n = 1,000,000</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="623.7" cy="78.0" r="5.5"/>
<text class="val" x="634.7" y="82.0">1.0 s</text>
</g>
<g class="row gb green"><title>genomeblocks (cgranges): 1.6 s · n = 1,000,000</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks (cgranges)</text>
<circle class="dot" cx="655.6" cy="106.0" r="5.5"/>
<text class="val" x="666.6" y="110.0">1.6 s</text>
</g>
<g class="row"><title>bedtools (CLI): 1.8 s · n = 1,000,000</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">bedtools (CLI)</text>
<circle class="dot" cx="665.4" cy="134.0" r="5.5"/>
<text class="val" x="676.4" y="138.0">1.8 s</text>
</g>
</svg>
</div><figcaption><strong>A &amp; B</strong> on 1M intervals per set (rows of A overlapping B); every engine returns the same rows.</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 154" role="img" aria-label="merge at 1,000,000 intervals" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="128"/>
<text class="tick" x="300.0" y="144" text-anchor="middle">10 ms</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="128"/>
<text class="tick" x="461.3" y="144" text-anchor="middle">100 ms</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="128"/>
<text class="tick" x="622.7" y="144" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="128"/>
<text class="tick" x="784.0" y="144" text-anchor="middle">10 s</text>
<g class="row gb green"><title>genomeblocks: 135 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks</text>
<circle class="dot" cx="482.5" cy="22.0" r="5.5"/>
<text class="val" x="493.5" y="26.0">135 ms</text>
</g>
<g class="row"><title>pyranges: 164 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pyranges</text>
<circle class="dot" cx="495.8" cy="50.0" r="5.5"/>
<text class="val" x="506.8" y="54.0">164 ms</text>
</g>
<g class="row"><title>bioframe: 699 ms · n = 1,000,000</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">bioframe</text>
<circle class="dot" cx="597.6" cy="78.0" r="5.5"/>
<text class="val" x="608.6" y="82.0">699 ms</text>
</g>
<g class="row"><title>bedtools (CLI): 5.1 s · n = 1,000,000</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">bedtools (CLI)</text>
<circle class="dot" cx="736.4" cy="106.0" r="5.5"/>
<text class="val" x="747.4" y="110.0">5.1 s</text>
</g>
</svg>
</div><figcaption><strong>merge</strong> on 1M intervals per set (sort + merge of A ∪ B); every engine returns the same rows.</figcaption></figure>

On whole-set operations at 1M intervals the numpy kernel does `A & B` in 87 ms and `merge` in 135 ms, against 336 ms and 164 ms for pyranges; the [Backends]({{ '/benchmarks/#backends' | relative_url }}) section times the other interval engines behind the same call.

<details class="gb-table"><summary>All Loci timings</summary><table><thead><tr><th>op</th><th>engine</th><th>n = 1k</th><th>n = 10k</th><th>n = 100k</th><th>n = 1M</th></tr></thead><tbody><tr><td>intersect</td><td>genomeblocks</td><td>268 µs</td><td>901 µs</td><td>7.4 ms</td><td>87 ms</td></tr><tr><td>intersect</td><td>genomeblocks (cgranges)</td><td>1.1 ms</td><td>9.8 ms</td><td>99 ms</td><td>1.6 s</td></tr><tr><td>intersect</td><td>genomeblocks (cgranges, warm index)</td><td>1.1 ms</td><td>10 ms</td><td>99 ms</td><td>1.4 s</td></tr><tr><td>intersect</td><td>pyranges</td><td>49 ms</td><td>44 ms</td><td>62 ms</td><td>336 ms</td></tr><tr><td>intersect</td><td>bioframe</td><td>7.3 ms</td><td>12 ms</td><td>61 ms</td><td>1.0 s</td></tr><tr><td>intersect</td><td>bedtools (CLI)</td><td>16 ms</td><td>36 ms</td><td>169 ms</td><td>1.8 s</td></tr><tr><td>intersect</td><td>intervaltree</td><td>24 ms</td><td>464 ms</td><td>22 s</td><td>–</td></tr><tr><td>intersect</td><td>naive double loop</td><td>36 ms</td><td>2.2 s</td><td>–</td><td>–</td></tr><tr><td>merge</td><td>genomeblocks</td><td>366 µs</td><td>2.0 ms</td><td>18 ms</td><td>135 ms</td></tr><tr><td>merge</td><td>pyranges</td><td>28 ms</td><td>30 ms</td><td>52 ms</td><td>164 ms</td></tr><tr><td>merge</td><td>bioframe</td><td>16 ms</td><td>23 ms</td><td>102 ms</td><td>699 ms</td></tr><tr><td>merge</td><td>bedtools (CLI)</td><td>27 ms</td><td>92 ms</td><td>543 ms</td><td>5.1 s</td></tr><tr><td>make</td><td>genomeblocks</td><td>2.2 ms</td><td>6.6 ms</td><td>56 ms</td><td>394 ms</td></tr><tr><td>make</td><td>pyranges</td><td>5.7 ms</td><td>13 ms</td><td>91 ms</td><td>928 ms</td></tr><tr><td>make</td><td>pandas.read_csv</td><td>1.5 ms</td><td>4.2 ms</td><td>29 ms</td><td>253 ms</td></tr></tbody></table></details>

## Signal
{: .sec-green #signal }

Filling a heatmap matrix (200 bins, ±3 kb) is one native call per region and track: **13k regions/s** on one core, 98× deepTools `computeMatrix` and 386× pyBigWig's binned `stats()`. Every engine returns the same matrix (deepTools differs only at bin edges, r ≥ 0.99).

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
<g class="row gb green"><title>genomeblocks · pybigtools (exact): 13k · 5,000 regions</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks · pybigtools (exact)</text>
<circle class="dot" cx="678.5" cy="22.0" r="5.5"/>
<text class="val" x="689.5" y="26.0">13k</text>
</g>
<g class="row gb green"><title>genomeblocks · pybigtools (zoom, exact=False): 13k · 5,000 regions</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks · pybigtools (zoom, exact=False)</text>
<circle class="dot" cx="676.6" cy="50.0" r="5.5"/>
<text class="val" x="687.6" y="54.0">13k</text>
</g>
<g class="row"><title>pybigtools values() + numpy binning: 5.3k · 5,000 regions</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">pybigtools values() + numpy binning</text>
<circle class="dot" cx="629.2" cy="78.0" r="5.5"/>
<text class="val" x="640.2" y="82.0">5.3k</text>
</g>
<g class="row"><title>pyBigWig values() + numpy binning: 3.7k · 5,000 regions</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyBigWig values() + numpy binning</text>
<circle class="dot" cx="611.4" cy="106.0" r="5.5"/>
<text class="val" x="622.4" y="110.0">3.7k</text>
</g>
<g class="row gb green"><title>genomeblocks · pure-Python reader: 2.3k · 5,000 regions</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">genomeblocks · pure-Python reader</text>
<circle class="dot" cx="584.9" cy="134.0" r="5.5"/>
<text class="val" x="595.9" y="138.0">2.3k</text>
</g>
<g class="row"><title>deepTools computeMatrix (-p 1): 137 · 5,000 regions</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">deepTools computeMatrix (-p 1)</text>
<circle class="dot" cx="437.4" cy="162.0" r="5.5"/>
<text class="val" x="448.4" y="166.0">137</text>
</g>
<g class="row"><title>pyBigWig stats(nBins) zoom: 35 · 300 regions</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">pyBigWig stats(nBins) zoom</text>
<circle class="dot" cx="365.5" cy="190.0" r="5.5"/>
<text class="val" x="376.5" y="194.0">35</text>
</g>
<g class="row"><title>pyBigWig stats(nBins) exact: 35 · 300 regions</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">pyBigWig stats(nBins) exact</text>
<circle class="dot" cx="365.5" cy="218.0" r="5.5"/>
<text class="val" x="376.5" y="222.0">35</text>
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
<g class="line"><title>genomeblocks, exact=False</title><polyline points="64.0,129.8 268.8,124.6 445.2,80.8 561.8,55.4"/>
<circle cx="64.0" cy="129.8" r="4"><title>genomeblocks, exact=False · 1k: 11k</title></circle>
<circle cx="268.8" cy="124.6" r="4"><title>genomeblocks, exact=False · 5k: 12k</title></circle>
<circle cx="445.2" cy="80.8" r="4"><title>genomeblocks, exact=False · 20k: 28k</title></circle>
<circle cx="561.8" cy="55.4" r="4"><title>genomeblocks, exact=False · 50k: 45k</title></circle>
</g>
<g class="line"><title>genomeblocks, pure-Python</title><polyline points="64.0,221.0 268.8,202.9 445.2,193.8 561.8,188.0"/>
<circle cx="64.0" cy="221.0" r="4"><title>genomeblocks, pure-Python · 1k: 1.9k</title></circle>
<circle cx="268.8" cy="202.9" r="4"><title>genomeblocks, pure-Python · 5k: 2.7k</title></circle>
<circle cx="445.2" cy="193.8" r="4"><title>genomeblocks, pure-Python · 20k: 3.3k</title></circle>
<circle cx="561.8" cy="188.0" r="4"><title>genomeblocks, pure-Python · 50k: 3.6k</title></circle>
</g>
<g class="line"><title>pyBigWig values()</title><polyline points="64.0,192.4 268.8,184.5 445.2,180.9 561.8,180.9"/>
<circle cx="64.0" cy="192.4" r="4"><title>pyBigWig values() · 1k: 3.4k</title></circle>
<circle cx="268.8" cy="184.5" r="4"><title>pyBigWig values() · 5k: 3.9k</title></circle>
<circle cx="445.2" cy="180.9" r="4"><title>pyBigWig values() · 20k: 4.2k</title></circle>
<circle cx="561.8" cy="180.9" r="4"><title>pyBigWig values() · 50k: 4.2k</title></circle>
</g>
<g class="line gb green"><title>genomeblocks (default)</title><polyline points="64.0,128.0 268.8,124.5 445.2,80.7 561.8,54.9"/>
<circle cx="64.0" cy="128.0" r="4"><title>genomeblocks (default) · 1k: 11k</title></circle>
<circle cx="268.8" cy="124.5" r="4"><title>genomeblocks (default) · 5k: 12k</title></circle>
<circle cx="445.2" cy="80.7" r="4"><title>genomeblocks (default) · 20k: 28k</title></circle>
<circle cx="561.8" cy="54.9" r="4"><title>genomeblocks (default) · 50k: 46k</title></circle>
</g>
<text class="lab gbl" x="660" y="58.9">genomeblocks (default)</text>
<line class="leader" x1="567.8" y1="54.9" x2="656" y2="54.9"/>
<text class="lab" x="660" y="72.9">genomeblocks, exact=False</text>
<line class="leader" x1="567.8" y1="68.9" x2="656" y2="68.9"/>
<text class="lab" x="660" y="184.9">pyBigWig values()</text>
<line class="leader" x1="567.8" y1="180.9" x2="656" y2="180.9"/>
<text class="lab" x="660" y="198.9">genomeblocks, pure-Python</text>
<line class="leader" x1="567.8" y1="194.9" x2="656" y2="194.9"/>
</svg>
</div><figcaption><strong>Throughput vs number of regions</strong> (regions/s, log–log). The native path speeds up with batch size as fixed costs amortise; per-base readers stay flat.</figcaption></figure>

Many tracks scale with processes: 16 bigWigs × 5,000 regions take 5.3 s sequentially and 1.4 s with `workers=8`. Threads do not help (pybigtools serialises Python threads), which is why `signal()` uses processes.

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
<g class="row gb green"><title>genomeblocks · processes · 8 workers: 1.4 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks · processes · 8 workers</text>
<circle class="dot" cx="440.1" cy="22.0" r="5.5"/>
<text class="val" x="451.1" y="26.0">1.4 s</text>
</g>
<g class="row gb green"><title>genomeblocks · processes · 4 workers: 2.1 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks · processes · 4 workers</text>
<circle class="dot" cx="459.8" cy="50.0" r="5.5"/>
<text class="val" x="470.8" y="54.0">2.1 s</text>
</g>
<g class="row gb green"><title>genomeblocks · processes · 2 workers: 3.6 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">genomeblocks · processes · 2 workers</text>
<circle class="dot" cx="488.6" cy="78.0" r="5.5"/>
<text class="val" x="499.6" y="82.0">3.6 s</text>
</g>
<g class="row gb green"><title>genomeblocks · processes · 1 worker: 5.3 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks · processes · 1 worker</text>
<circle class="dot" cx="508.5" cy="106.0" r="5.5"/>
<text class="val" x="519.5" y="110.0">5.3 s</text>
</g>
<g class="row"><title>threads (same chunks) · 2 workers: 19 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">threads (same chunks) · 2 workers</text>
<circle class="dot" cx="575.5" cy="134.0" r="5.5"/>
<text class="val" x="586.5" y="138.0">19 s</text>
</g>
<g class="row"><title>threads (same chunks) · 4 workers: 19 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">threads (same chunks) · 4 workers</text>
<circle class="dot" cx="576.7" cy="162.0" r="5.5"/>
<text class="val" x="587.7" y="166.0">19 s</text>
</g>
<g class="row"><title>threads (same chunks) · 8 workers: 19 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">threads (same chunks) · 8 workers</text>
<circle class="dot" cx="576.9" cy="190.0" r="5.5"/>
<text class="val" x="587.9" y="194.0">19 s</text>
</g>
<g class="row"><title>deepTools computeMatrix · 8 workers: 96 s · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">deepTools computeMatrix · 8 workers</text>
<circle class="dot" cx="660.8" cy="218.0" r="5.5"/>
<text class="val" x="671.8" y="222.0">96 s</text>
</g>
<g class="row"><title>deepTools computeMatrix · 4 workers: 2.8 min · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="232.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="246.0" x2="784" y2="246.0"/>
<text class="lab" x="286" y="250.0" text-anchor="end">deepTools computeMatrix · 4 workers</text>
<circle class="dot" cx="689.4" cy="246.0" r="5.5"/>
<text class="val" x="700.4" y="250.0">2.8 min</text>
</g>
<g class="row"><title>deepTools computeMatrix · 1 worker: 9.8 min · 16 tracks × 5,000 regions</title>
<rect class="hit" x="0" y="260.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="274.0" x2="784" y2="274.0"/>
<text class="lab" x="286" y="278.0" text-anchor="end">deepTools computeMatrix · 1 worker</text>
<circle class="dot" cx="756.3" cy="274.0" r="5.5"/>
<text class="val" x="767.3" y="278.0">9.8 min</text>
</g>
</svg>
</div><figcaption><strong>16 bigWigs × 5,000 regions</strong>: worker processes vs threads vs deepTools <code>-p</code> (wall time; log scale).</figcaption></figure>

## Atlas
{: .sec-navy #atlas }

Testing a 20,000-peak query against 500 peak files is one sparse row-sum plus a vectorised Fisher test: **11 ms** per query, including building the result table.

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
<g class="row gb navy"><title>Atlas sparse row-sum only: 3.7 ms · 500 tracks</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">Atlas sparse row-sum only</text>
<circle class="dot" cx="355.5" cy="22.0" r="5.5"/>
<text class="val" x="366.5" y="26.0">3.7 ms</text>
</g>
<g class="row gb navy"><title>Atlas.search (with Fisher + table): 11 ms · 500 tracks</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">Atlas.search (with Fisher + table)</text>
<circle class="dot" cx="401.5" cy="50.0" r="5.5"/>
<text class="val" x="412.5" y="54.0">11 ms</text>
</g>
<g class="row"><title>per-track overlap_any loop: 442 ms · 500 tracks</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">per-track overlap_any loop</text>
<circle class="dot" cx="556.1" cy="78.0" r="5.5"/>
<text class="val" x="567.1" y="82.0">442 ms</text>
</g>
<g class="row"><title>pyranges overlap(), per track: 21 s · 500 tracks</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">pyranges overlap(), per track</text>
<circle class="dot" cx="718.9" cy="106.0" r="5.5"/>
<text class="val" x="729.9" y="110.0">21 s</text>
</g>
</svg>
</div><figcaption><strong>One query against 500 peak files</strong> (wall time; log scale). GIGGLE was not run (it is built from source; set <code>GIGGLE=</code> to include it).</figcaption></figure>

<details class="gb-table"><summary>Index build, load, accuracy</summary><table><thead><tr><th>step</th><th>bin size</th><th>workers</th><th>time</th><th>index size</th></tr></thead><tbody><tr><td>build</td><td>200 bp</td><td>3</td><td>10 s</td><td>190 MB</td></tr><tr><td>build</td><td>1,000 bp</td><td>1</td><td>7.8 s</td><td>68 MB</td></tr><tr><td>build</td><td>1,000 bp</td><td>3</td><td>3.5 s</td><td>68 MB</td></tr><tr><td>build</td><td>5,000 bp</td><td>3</td><td>2.8 s</td><td>43 MB</td></tr><tr><td>load (Atlas.load (npz))</td><td></td><td></td><td>266 ms</td><td></td></tr><tr><td>load (fresh python: import + load + search)</td><td></td><td></td><td>4.2 s</td><td></td></tr><tr><td>accuracy vs exact overlaps</td><td>200 bp</td><td></td><td>Spearman 0.999</td><td>top-25 shared 21</td></tr><tr><td>accuracy vs exact overlaps</td><td>1,000 bp</td><td></td><td>Spearman 0.999</td><td>top-25 shared 20</td></tr><tr><td>accuracy vs exact overlaps</td><td>5,000 bp</td><td></td><td>Spearman 0.993</td><td>top-25 shared 18</td></tr><tr><td>bootstrap (10 shuffles)</td><td></td><td></td><td>8.7 ms per shuffle</td><td></td></tr><tr><td>bootstrap (100 shuffles)</td><td></td><td></td><td>7.5 ms per shuffle</td><td></td></tr></tbody></table></details>

## Motifs
{: .sec-navy #motifs }

Counting hits of 100 JASPAR motifs in 1,000 windows of 500 bp: **58 ms** with genomeblocks on MOODS (the windows joined into one block, the library scanned in one pass), 131× faster than FIMO. The 993 hits are identical to those of lightmotif per window, MOODS, numpy and Biopython; FIMO reports 996.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 238" role="img" aria-label="Motif scanning engines" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="212"/>
<text class="tick" x="300.0" y="228" text-anchor="middle">10 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="212"/>
<text class="tick" x="421.0" y="228" text-anchor="middle">100 ms</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="212"/>
<text class="tick" x="542.0" y="228" text-anchor="middle">1 s</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="212"/>
<text class="tick" x="663.0" y="228" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="212"/>
<text class="tick" x="784.0" y="228" text-anchor="middle">100 s</text>
<g class="row gb navy"><title>genomeblocks · lightmotif: 42 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks · lightmotif</text>
<circle class="dot" cx="374.9" cy="22.0" r="5.5"/>
<text class="val" x="385.9" y="26.0">42 ms</text>
</g>
<g class="row gb navy"><title>genomeblocks · MOODS (default): 58 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks · MOODS (default)</text>
<circle class="dot" cx="392.2" cy="50.0" r="5.5"/>
<text class="val" x="403.2" y="54.0">58 ms</text>
</g>
<g class="row"><title>MOODS (C++): 69 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">MOODS (C++)</text>
<circle class="dot" cx="401.1" cy="78.0" r="5.5"/>
<text class="val" x="412.1" y="82.0">69 ms</text>
</g>
<g class="row"><title>lightmotif, re-striped per motif: 289 ms · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">lightmotif, re-striped per motif</text>
<circle class="dot" cx="476.7" cy="106.0" r="5.5"/>
<text class="val" x="487.7" y="110.0">289 ms</text>
</g>
<g class="row"><title>numpy sliding window: 5.4 s · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">numpy sliding window</text>
<circle class="dot" cx="631.0" cy="134.0" r="5.5"/>
<text class="val" x="642.0" y="138.0">5.4 s</text>
</g>
<g class="row"><title>MEME FIMO (CLI): 7.6 s · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">MEME FIMO (CLI)</text>
<circle class="dot" cx="648.3" cy="162.0" r="5.5"/>
<text class="val" x="659.3" y="166.0">7.6 s</text>
</g>
<g class="row"><title>Biopython PSSM.search: 11 s · 1,000 windows × 100 motifs</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">Biopython PSSM.search</text>
<circle class="dot" cx="667.0" cy="190.0" r="5.5"/>
<text class="val" x="678.0" y="194.0">11 s</text>
</g>
</svg>
</div><figcaption><strong>Motif scanning engines</strong> on identical PSSMs and windows (wall time; log scale).</figcaption></figure>

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 238" role="img" aria-label="Whole library" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="212"/>
<text class="tick" x="300.0" y="228" text-anchor="middle">100 ms</text>
<line class="grid" x1="421.0" y1="8" x2="421.0" y2="212"/>
<text class="tick" x="421.0" y="228" text-anchor="middle">1 s</text>
<line class="grid" x1="542.0" y1="8" x2="542.0" y2="212"/>
<text class="tick" x="542.0" y="228" text-anchor="middle">10 s</text>
<line class="grid" x1="663.0" y1="8" x2="663.0" y2="212"/>
<text class="tick" x="663.0" y="228" text-anchor="middle">100 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="212"/>
<text class="tick" x="784.0" y="228" text-anchor="middle">1,000 s</text>
<g class="row gb navy"><title>genomeblocks · lightmotif · 1k windows: 368 ms · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">genomeblocks · lightmotif · 1k windows</text>
<circle class="dot" cx="368.5" cy="22.0" r="5.5"/>
<text class="val" x="379.5" y="26.0">368 ms</text>
</g>
<g class="row gb navy"><title>genomeblocks · MOODS (default) · 1k windows: 565 ms · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">genomeblocks · MOODS (default) · 1k windows</text>
<circle class="dot" cx="391.0" cy="50.0" r="5.5"/>
<text class="val" x="402.0" y="54.0">565 ms</text>
</g>
<g class="row"><title>MOODS (C++) · 1k windows: 721 ms · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">MOODS (C++) · 1k windows</text>
<circle class="dot" cx="403.8" cy="78.0" r="5.5"/>
<text class="val" x="414.8" y="82.0">721 ms</text>
</g>
<g class="row gb navy"><title>genomeblocks · lightmotif · 5k windows: 1.4 s · 5,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">genomeblocks · lightmotif · 5k windows</text>
<circle class="dot" cx="439.0" cy="106.0" r="5.5"/>
<text class="val" x="450.0" y="110.0">1.4 s</text>
</g>
<g class="row gb navy"><title>genomeblocks · MOODS (default) · 5k windows: 1.4 s · 5,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">genomeblocks · MOODS (default) · 5k windows</text>
<circle class="dot" cx="439.7" cy="134.0" r="5.5"/>
<text class="val" x="450.7" y="138.0">1.4 s</text>
</g>
<g class="row"><title>MOODS (C++) · 5k windows: 2.2 s · 5,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">MOODS (C++) · 5k windows</text>
<circle class="dot" cx="462.9" cy="162.0" r="5.5"/>
<text class="val" x="473.9" y="166.0">2.2 s</text>
</g>
<g class="row"><title>MEME FIMO (CLI) · 1k windows: 78 s · 1,000 windows × 1019 motifs</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">MEME FIMO (CLI) · 1k windows</text>
<circle class="dot" cx="650.1" cy="190.0" r="5.5"/>
<text class="val" x="661.1" y="194.0">78 s</text>
</g>
</svg>
</div><figcaption><strong>The whole JASPAR library</strong> (1,019 motifs) on one core (wall time; log scale). Against MOODS: genomeblocks is 1.3× faster at 1,000 windows, genomeblocks is 1.6× faster at 5,000 windows. FIMO was timed at 1,000 windows only.</figcaption></figure>

<details class="gb-table"><summary>Worker processes</summary><table><thead><tr><th>workers</th><th>time</th><th>speed-up</th></tr></thead><tbody><tr><td>1</td><td>1.4 s</td><td>1.0×</td></tr><tr><td>2</td><td>872 ms</td><td>1.6×</td></tr><tr><td>4</td><td>641 ms</td><td>2.2×</td></tr><tr><td>8</td><td>530 ms</td><td>2.7×</td></tr></tbody></table></details>

## Hi-C pairs
{: .sec-purple #pairs }

Streaming 5M read pairs into 50 kb windows × partner chromosome runs at **1.6M pairs/s**, close to the speed of just parsing the file, and 2.9× the rate of `cooler cload`.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 210" role="img" aria-label="Pairs per second" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="184"/>
<text class="tick" x="300.0" y="200" text-anchor="middle">10k</text>
<line class="grid" x1="461.3" y1="8" x2="461.3" y2="184"/>
<text class="tick" x="461.3" y="200" text-anchor="middle">100k</text>
<line class="grid" x1="622.7" y1="8" x2="622.7" y2="184"/>
<text class="tick" x="622.7" y="200" text-anchor="middle">1M</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="184"/>
<text class="tick" x="784.0" y="200" text-anchor="middle">10M</text>
<g class="row"><title>parse only (I/O floor): 2.7M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">parse only (I/O floor)</text>
<circle class="dot" cx="692.7" cy="22.0" r="5.5"/>
<text class="val" x="703.7" y="26.0">2.7M</text>
</g>
<g class="row gb purple"><title>count_pairs · one partner: 1.6M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">count_pairs · one partner</text>
<circle class="dot" cx="657.6" cy="50.0" r="5.5"/>
<text class="val" x="668.6" y="54.0">1.6M</text>
</g>
<g class="row gb purple"><title>count_pairs · 50 kb × partner: 1.6M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">count_pairs · 50 kb × partner</text>
<circle class="dot" cx="654.2" cy="78.0" r="5.5"/>
<text class="val" x="665.2" y="82.0">1.6M</text>
</g>
<g class="row gb purple"><title>count_pairs_2d · 500 kb: 1.5M · 5,000,000 pairs</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">count_pairs_2d · 500 kb</text>
<circle class="dot" cx="651.4" cy="106.0" r="5.5"/>
<text class="val" x="662.4" y="110.0">1.5M</text>
</g>
<g class="row"><title>cooler cload · 500 kb (CLI): 542k · 5,000,000 pairs</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">cooler cload · 500 kb (CLI)</text>
<circle class="dot" cx="579.8" cy="134.0" r="5.5"/>
<text class="val" x="590.8" y="138.0">542k</text>
</g>
<g class="row"><title>per-pair Python loop: 37k · 200,000 pairs</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">per-pair Python loop</text>
<circle class="dot" cx="392.4" cy="162.0" r="5.5"/>
<text class="val" x="403.4" y="166.0">37k</text>
</g>
</svg>
</div><figcaption><strong>Read pairs per second</strong> (higher is better; log scale). The naive loop was timed on the first 200k pairs.</figcaption></figure>

## Genes
{: .sec-navy #genes }

A GENCODE-shaped GTF with 20,000 genes, 70,028 transcripts and 876,947 lines.

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 630" role="img" aria-label="Gene-model steps" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="604"/>
<text class="tick" x="300.0" y="620" text-anchor="middle">1 ms</text>
<line class="grid" x1="396.8" y1="8" x2="396.8" y2="604"/>
<text class="tick" x="396.8" y="620" text-anchor="middle">10 ms</text>
<line class="grid" x1="493.6" y1="8" x2="493.6" y2="604"/>
<text class="tick" x="493.6" y="620" text-anchor="middle">100 ms</text>
<line class="grid" x1="590.4" y1="8" x2="590.4" y2="604"/>
<text class="tick" x="590.4" y="620" text-anchor="middle">1 s</text>
<line class="grid" x1="687.2" y1="8" x2="687.2" y2="604"/>
<text class="tick" x="687.2" y="620" text-anchor="middle">10 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="604"/>
<text class="tick" x="784.0" y="620" text-anchor="middle">100 s</text>
<g class="row gb navy"><title>nearest_genes · 10,000: 16 ms</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">nearest_genes · 10,000</text>
<circle class="dot" cx="416.0" cy="22.0" r="5.5"/>
<text class="val" x="427.0" y="26.0">16 ms</text>
</g>
<g class="row gb navy"><title>annotations (genomeblocks) · 10,000: 17 ms</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">annotations (genomeblocks) · 10,000</text>
<circle class="dot" cx="419.8" cy="50.0" r="5.5"/>
<text class="val" x="430.8" y="54.0">17 ms</text>
</g>
<g class="row gb navy"><title>select_isoforms (peaks): 40 ms</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">select_isoforms (peaks)</text>
<circle class="dot" cx="454.8" cy="78.0" r="5.5"/>
<text class="val" x="465.8" y="82.0">40 ms</text>
</g>
<g class="row gb navy"><title>annotations (ncls) · 10,000: 52 ms</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">annotations (ncls) · 10,000</text>
<circle class="dot" cx="466.4" cy="106.0" r="5.5"/>
<text class="val" x="477.4" y="110.0">52 ms</text>
</g>
<g class="row gb navy"><title>nearest_genes · 100,000: 92 ms</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">nearest_genes · 100,000</text>
<circle class="dot" cx="490.0" cy="134.0" r="5.5"/>
<text class="val" x="501.0" y="138.0">92 ms</text>
</g>
<g class="row gb navy"><title>annotations (genomeblocks) · 100,000: 98 ms</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">annotations (genomeblocks) · 100,000</text>
<circle class="dot" cx="492.8" cy="162.0" r="5.5"/>
<text class="val" x="503.8" y="166.0">98 ms</text>
</g>
<g class="row gb navy"><title>annot index build (prom/exon/UTR merge): 103 ms</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">annot index build (prom/exon/UTR merge)</text>
<circle class="dot" cx="495.0" cy="190.0" r="5.5"/>
<text class="val" x="506.0" y="194.0">103 ms</text>
</g>
<g class="row gb navy"><title>annotations (ncls) · 100,000: 149 ms</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">annotations (ncls) · 100,000</text>
<circle class="dot" cx="510.4" cy="218.0" r="5.5"/>
<text class="val" x="521.4" y="222.0">149 ms</text>
</g>
<g class="row gb navy"><title>annotations (cgranges) · 10,000: 168 ms</title>
<rect class="hit" x="0" y="232.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="246.0" x2="784" y2="246.0"/>
<text class="lab" x="286" y="250.0" text-anchor="end">annotations (cgranges) · 10,000</text>
<circle class="dot" cx="515.4" cy="246.0" r="5.5"/>
<text class="val" x="526.4" y="250.0">168 ms</text>
</g>
<g class="row gb navy"><title>annotations (bioframe) · 10,000: 194 ms</title>
<rect class="hit" x="0" y="260.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="274.0" x2="784" y2="274.0"/>
<text class="lab" x="286" y="278.0" text-anchor="end">annotations (bioframe) · 10,000</text>
<circle class="dot" cx="521.5" cy="274.0" r="5.5"/>
<text class="val" x="532.5" y="278.0">194 ms</text>
</g>
<g class="row gb navy"><title>select_isoforms (peaks + bigWig): 292 ms</title>
<rect class="hit" x="0" y="288.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="302.0" x2="784" y2="302.0"/>
<text class="lab" x="286" y="306.0" text-anchor="end">select_isoforms (peaks + bigWig)</text>
<circle class="dot" cx="538.6" cy="302.0" r="5.5"/>
<text class="val" x="549.6" y="306.0">292 ms</text>
</g>
<g class="row gb navy"><title>annotations (bioframe) · 100,000: 416 ms</title>
<rect class="hit" x="0" y="316.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="330.0" x2="784" y2="330.0"/>
<text class="lab" x="286" y="334.0" text-anchor="end">annotations (bioframe) · 100,000</text>
<circle class="dot" cx="553.5" cy="330.0" r="5.5"/>
<text class="val" x="564.5" y="334.0">416 ms</text>
</g>
<g class="row gb navy"><title>annotations (pyranges) · 10,000: 426 ms</title>
<rect class="hit" x="0" y="344.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="358.0" x2="784" y2="358.0"/>
<text class="lab" x="286" y="362.0" text-anchor="end">annotations (pyranges) · 10,000</text>
<circle class="dot" cx="554.5" cy="358.0" r="5.5"/>
<text class="val" x="565.5" y="362.0">426 ms</text>
</g>
<g class="row gb navy"><title>annotations (cgranges) · 100,000: 461 ms</title>
<rect class="hit" x="0" y="372.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="386.0" x2="784" y2="386.0"/>
<text class="lab" x="286" y="390.0" text-anchor="end">annotations (cgranges) · 100,000</text>
<circle class="dot" cx="557.9" cy="386.0" r="5.5"/>
<text class="val" x="568.9" y="390.0">461 ms</text>
</g>
<g class="row gb navy"><title>Genes.make (polars parser): 598 ms</title>
<rect class="hit" x="0" y="400.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="414.0" x2="784" y2="414.0"/>
<text class="lab" x="286" y="418.0" text-anchor="end">Genes.make (polars parser)</text>
<circle class="dot" cx="568.8" cy="414.0" r="5.5"/>
<text class="val" x="579.8" y="418.0">598 ms</text>
</g>
<g class="row gb navy"><title>Genes.make (GTF parse): 609 ms</title>
<rect class="hit" x="0" y="428.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="442.0" x2="784" y2="442.0"/>
<text class="lab" x="286" y="446.0" text-anchor="end">Genes.make (GTF parse)</text>
<circle class="dot" cx="569.5" cy="442.0" r="5.5"/>
<text class="val" x="580.5" y="446.0">609 ms</text>
</g>
<g class="row gb navy"><title>annotations (pyranges) · 100,000: 649 ms</title>
<rect class="hit" x="0" y="456.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="470.0" x2="784" y2="470.0"/>
<text class="lab" x="286" y="474.0" text-anchor="end">annotations (pyranges) · 100,000</text>
<circle class="dot" cx="572.2" cy="470.0" r="5.5"/>
<text class="val" x="583.2" y="474.0">649 ms</text>
</g>
<g class="row gb navy"><title>select_isoforms (bigWig only): 1.3 s</title>
<rect class="hit" x="0" y="484.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="498.0" x2="784" y2="498.0"/>
<text class="lab" x="286" y="502.0" text-anchor="end">select_isoforms (bigWig only)</text>
<circle class="dot" cx="601.4" cy="498.0" r="5.5"/>
<text class="val" x="612.4" y="502.0">1.3 s</text>
</g>
<g class="row gb navy"><title>annotations (bedtools) · 10,000: 1.4 s</title>
<rect class="hit" x="0" y="512.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="526.0" x2="784" y2="526.0"/>
<text class="lab" x="286" y="530.0" text-anchor="end">annotations (bedtools) · 10,000</text>
<circle class="dot" cx="605.6" cy="526.0" r="5.5"/>
<text class="val" x="616.6" y="530.0">1.4 s</text>
</g>
<g class="row gb navy"><title>annotations (bedtools) · 100,000: 2.6 s</title>
<rect class="hit" x="0" y="540.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="554.0" x2="784" y2="554.0"/>
<text class="lab" x="286" y="558.0" text-anchor="end">annotations (bedtools) · 100,000</text>
<circle class="dot" cx="630.6" cy="554.0" r="5.5"/>
<text class="val" x="641.6" y="558.0">2.6 s</text>
</g>
<g class="row gb navy"><title>Genes.make (pandas parser): 9.3 s</title>
<rect class="hit" x="0" y="568.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="582.0" x2="784" y2="582.0"/>
<text class="lab" x="286" y="586.0" text-anchor="end">Genes.make (pandas parser)</text>
<circle class="dot" cx="684.3" cy="582.0" r="5.5"/>
<text class="val" x="695.3" y="586.0">9.3 s</text>
</g>
</svg>
</div><figcaption><strong>Gene-model steps</strong> (wall time; log scale). Labelling uses the cached annotation index; the counts are CREs labelled.</figcaption></figure>

## Architecture
{: .sec-purple #architecture }

The pipeline on 100k CREs and 50k loops (114,367 edges), step by step, then the graph algorithms on every installed graph engine and the parquet round trip:

<figure class="gb-fig gb-chart"><div class="gb-fig-body">
<svg class="gbc" viewBox="0 0 880 518" role="img" aria-label="Architecture steps" xmlns="http://www.w3.org/2000/svg">
<line class="grid" x1="300.0" y1="8" x2="300.0" y2="492"/>
<text class="tick" x="300.0" y="508" text-anchor="middle">100 µs</text>
<line class="grid" x1="396.8" y1="8" x2="396.8" y2="492"/>
<text class="tick" x="396.8" y="508" text-anchor="middle">1 ms</text>
<line class="grid" x1="493.6" y1="8" x2="493.6" y2="492"/>
<text class="tick" x="493.6" y="508" text-anchor="middle">10 ms</text>
<line class="grid" x1="590.4" y1="8" x2="590.4" y2="492"/>
<text class="tick" x="590.4" y="508" text-anchor="middle">100 ms</text>
<line class="grid" x1="687.2" y1="8" x2="687.2" y2="492"/>
<text class="tick" x="687.2" y="508" text-anchor="middle">1 s</text>
<line class="grid" x1="784.0" y1="8" x2="784.0" y2="492"/>
<text class="tick" x="784.0" y="508" text-anchor="middle">10 s</text>
<g class="row gb purple"><title>strength (numpy bincount): 994 µs</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">strength (numpy bincount)</text>
<circle class="dot" cx="396.6" cy="22.0" r="5.5"/>
<text class="val" x="407.6" y="26.0">994 µs</text>
</g>
<g class="row gb purple"><title>pagerank (O/E weights) · graph-tool: 5.2 ms</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">pagerank (O/E weights) · graph-tool</text>
<circle class="dot" cx="465.8" cy="50.0" r="5.5"/>
<text class="val" x="476.8" y="54.0">5.2 ms</text>
</g>
<g class="row gb purple"><title>prime_hubs (elbow on strength): 8.2 ms</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">prime_hubs (elbow on strength)</text>
<circle class="dot" cx="485.0" cy="78.0" r="5.5"/>
<text class="val" x="496.0" y="82.0">8.2 ms</text>
</g>
<g class="row gb purple"><title>normalize (power-law O/E): 9.6 ms</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">normalize (power-law O/E)</text>
<circle class="dot" cx="491.8" cy="106.0" r="5.5"/>
<text class="val" x="502.8" y="110.0">9.6 ms</text>
</g>
<g class="row gb purple"><title>components · graph-tool: 14 ms</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">components · graph-tool</text>
<circle class="dot" cx="507.2" cy="134.0" r="5.5"/>
<text class="val" x="518.2" y="138.0">14 ms</text>
</g>
<g class="row gb purple"><title>support (CREs per TSS ± 5 kb): 17 ms</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">support (CREs per TSS ± 5 kb)</text>
<circle class="dot" cx="517.1" cy="162.0" r="5.5"/>
<text class="val" x="528.1" y="166.0">17 ms</text>
</g>
<g class="row"><title>components · scipy: 18 ms</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">components · scipy</text>
<circle class="dot" cx="518.4" cy="190.0" r="5.5"/>
<text class="val" x="529.4" y="194.0">18 ms</text>
</g>
<g class="row gb purple"><title>annotate (labels + genes): 35 ms</title>
<rect class="hit" x="0" y="204.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="218.0" x2="784" y2="218.0"/>
<text class="lab" x="286" y="222.0" text-anchor="end">annotate (labels + genes)</text>
<circle class="dot" cx="546.5" cy="218.0" r="5.5"/>
<text class="val" x="557.5" y="222.0">35 ms</text>
</g>
<g class="row gb purple"><title>save (parquet): 70 ms</title>
<rect class="hit" x="0" y="232.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="246.0" x2="784" y2="246.0"/>
<text class="lab" x="286" y="250.0" text-anchor="end">save (parquet)</text>
<circle class="dot" cx="575.6" cy="246.0" r="5.5"/>
<text class="val" x="586.6" y="250.0">70 ms</text>
</g>
<g class="row"><title>pagerank (O/E weights) · igraph: 73 ms</title>
<rect class="hit" x="0" y="260.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="274.0" x2="784" y2="274.0"/>
<text class="lab" x="286" y="278.0" text-anchor="end">pagerank (O/E weights) · igraph</text>
<circle class="dot" cx="577.0" cy="274.0" r="5.5"/>
<text class="val" x="588.0" y="278.0">73 ms</text>
</g>
<g class="row"><title>components · igraph: 74 ms</title>
<rect class="hit" x="0" y="288.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="302.0" x2="784" y2="302.0"/>
<text class="lab" x="286" y="306.0" text-anchor="end">components · igraph</text>
<circle class="dot" cx="577.8" cy="302.0" r="5.5"/>
<text class="val" x="588.8" y="306.0">74 ms</text>
</g>
<g class="row"><title>pagerank (O/E weights) · scipy: 81 ms</title>
<rect class="hit" x="0" y="316.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="330.0" x2="784" y2="330.0"/>
<text class="lab" x="286" y="334.0" text-anchor="end">pagerank (O/E weights) · scipy</text>
<circle class="dot" cx="581.7" cy="330.0" r="5.5"/>
<text class="val" x="592.7" y="334.0">81 ms</text>
</g>
<g class="row gb purple"><title>load (parquet): 88 ms</title>
<rect class="hit" x="0" y="344.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="358.0" x2="784" y2="358.0"/>
<text class="lab" x="286" y="362.0" text-anchor="end">load (parquet)</text>
<circle class="dot" cx="585.2" cy="358.0" r="5.5"/>
<text class="val" x="596.2" y="362.0">88 ms</text>
</g>
<g class="row gb purple"><title>make (50k loops → graph): 155 ms</title>
<rect class="hit" x="0" y="372.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="386.0" x2="784" y2="386.0"/>
<text class="lab" x="286" y="390.0" text-anchor="end">make (50k loops → graph)</text>
<circle class="dot" cx="608.7" cy="386.0" r="5.5"/>
<text class="val" x="619.7" y="390.0">155 ms</text>
</g>
<g class="row"><title>components · networkx: 488 ms</title>
<rect class="hit" x="0" y="400.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="414.0" x2="784" y2="414.0"/>
<text class="lab" x="286" y="418.0" text-anchor="end">components · networkx</text>
<circle class="dot" cx="657.1" cy="414.0" r="5.5"/>
<text class="val" x="668.1" y="418.0">488 ms</text>
</g>
<g class="row"><title>pagerank (O/E weights) · networkx: 859 ms</title>
<rect class="hit" x="0" y="428.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="442.0" x2="784" y2="442.0"/>
<text class="lab" x="286" y="446.0" text-anchor="end">pagerank (O/E weights) · networkx</text>
<circle class="dot" cx="680.8" cy="442.0" r="5.5"/>
<text class="val" x="691.8" y="446.0">859 ms</text>
</g>
<g class="row gb purple"><title>add_mcool (5 kb): 4.1 s</title>
<rect class="hit" x="0" y="456.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="470.0" x2="784" y2="470.0"/>
<text class="lab" x="286" y="474.0" text-anchor="end">add_mcool (5 kb)</text>
<circle class="dot" cx="746.8" cy="470.0" r="5.5"/>
<text class="val" x="757.8" y="474.0">4.1 s</text>
</g>
</svg>
</div><figcaption><strong>Architecture pipeline steps</strong> (wall time; log scale). Graph algorithms run once per installed engine, with identical results; the default engine in colour.</figcaption></figure>

## HiChIP short-range track
{: .sec-purple #shortrange }

From allValidPairs to a coverage bigWig: **17 s** and 1.0 GB peak memory with `genomeblocks.hichip`, against 2.3 min and 1.9 GB for the `awk | sort | bedtools` recipe, with the same ends and a byte-identical bedGraph.

<details class="gb-table"><summary>Steps</summary><table><thead><tr><th>implementation</th><th>step</th><th>time</th><th>peak memory</th></tr></thead><tbody><tr><td>shell</td><td>extract ends (awk)</td><td>6.1 s</td><td>20 MB</td></tr><tr><td>shell</td><td>fragments + sort + genomecov</td><td>2.0 min</td><td>1,903 MB</td></tr><tr><td>shell</td><td>bedGraph -&gt; bigWig</td><td>8.9 s</td><td>–</td></tr><tr><td>genomeblocks</td><td>ends</td><td>4.5 s</td><td>–</td></tr><tr><td>genomeblocks</td><td>write BED</td><td>897 ms</td><td>–</td></tr><tr><td>genomeblocks</td><td>coverage (bedGraph)</td><td>3.6 s</td><td>–</td></tr><tr><td>genomeblocks</td><td>bigWig</td><td>7.8 s</td><td>–</td></tr></tbody></table></details>

## Import cost
{: .sec-navy #import }

`import genomeblocks` is lazy, and `from genomeblocks import Loci` costs **201 ms** on top of starting Python (37 ms): matplotlib, pandas, scipy and graph-tool load only when a function needs them.

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
<g class="row gb navy"><title>import genomeblocks: 13 ms</title>
<rect class="hit" x="0" y="8.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="22.0" x2="784" y2="22.0"/>
<text class="lab" x="286" y="26.0" text-anchor="end">import genomeblocks</text>
<circle class="dot" cx="436.3" cy="22.0" r="5.5"/>
<text class="val" x="447.3" y="26.0">13 ms</text>
</g>
<g class="row gb navy"><title>+ Locus: 28 ms</title>
<rect class="hit" x="0" y="36.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="50.0" x2="784" y2="50.0"/>
<text class="lab" x="286" y="54.0" text-anchor="end">+ Locus</text>
<circle class="dot" cx="474.8" cy="50.0" r="5.5"/>
<text class="val" x="485.8" y="54.0">28 ms</text>
</g>
<g class="row gb navy"><title>+ Architecture (graph-tool): 155 ms</title>
<rect class="hit" x="0" y="64.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="78.0" x2="784" y2="78.0"/>
<text class="lab" x="286" y="82.0" text-anchor="end">+ Architecture (graph-tool)</text>
<circle class="dot" cx="564.9" cy="78.0" r="5.5"/>
<text class="val" x="575.9" y="82.0">155 ms</text>
</g>
<g class="row gb navy"><title>+ Genes: 178 ms</title>
<rect class="hit" x="0" y="92.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="106.0" x2="784" y2="106.0"/>
<text class="lab" x="286" y="110.0" text-anchor="end">+ Genes</text>
<circle class="dot" cx="572.3" cy="106.0" r="5.5"/>
<text class="val" x="583.3" y="110.0">178 ms</text>
</g>
<g class="row gb navy"><title>+ Loci (pulls signal, motifs, atlas, bedpe): 201 ms</title>
<rect class="hit" x="0" y="120.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="134.0" x2="784" y2="134.0"/>
<text class="lab" x="286" y="138.0" text-anchor="end">+ Loci (pulls signal, motifs, atlas, bedpe)</text>
<circle class="dot" cx="578.8" cy="134.0" r="5.5"/>
<text class="val" x="589.8" y="138.0">201 ms</text>
</g>
<g class="row gb navy"><title>+ browser (matplotlib): 544 ms</title>
<rect class="hit" x="0" y="148.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="162.0" x2="784" y2="162.0"/>
<text class="lab" x="286" y="166.0" text-anchor="end">+ browser (matplotlib)</text>
<circle class="dot" cx="631.0" cy="162.0" r="5.5"/>
<text class="val" x="642.0" y="166.0">544 ms</text>
</g>
<g class="row gb navy"><title>everything: 627 ms</title>
<rect class="hit" x="0" y="176.0" width="880" height="28"/>
<line class="rowline" x1="300" y1="190.0" x2="784" y2="190.0"/>
<text class="lab" x="286" y="194.0" text-anchor="end">everything</text>
<circle class="dot" cx="638.5" cy="190.0" r="5.5"/>
<text class="val" x="649.5" y="194.0">627 ms</text>
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

<p class="gb-env">Versions: genomeblocks 2.0.0, numpy 2.3.4, scipy 1.16.3, pandas 2.3.3, polars 1.44.2, pyarrow 21.0.0, narwhals 2.26.0, pybigtools 0.2.5, pyBigWig 0.3.24, cgranges installed, ncls installed, pyranges 0.1.4, bioframe 0.8.0, pybedtools 0.12.1, lightmotif 0.10.0, Bio 1.86, MOODS installed, graph_tool 2.98 (commit c96a6bf3, ), igraph 1.0.0, networkx 3.6.1, bedtools v2.31.1</p>
