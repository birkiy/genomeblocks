# Replicating ROSE's `*_REGION_TO_GENE.txt` with Genomeblocks

A walkthrough of how ROSE maps enhancers to genes, what's quirky about it,
and how a Genomeblocks-based pipeline reproduces the output **byte-identically**
in <10 s vs ROSE's ~5 min.

- Reference output: `…/SE_K562_H3K27ac/K562_H3K27ac_REP1_peaks_SuperStitched_REGION_TO_GENE.txt`
- Replication script: [rose_region_to_gene.py](rose_region_to_gene.py)
- Annotation: `/home/ualtintas/apps/ROSE/annotations/hg38_refseq.ucsc`

```bash
micromamba run -n notebook_cn1 python examples/rose_region_to_gene.py
diff -q K562_H3K27ac_REP1_peaks_SuperStitched_REGION_TO_GENE.gb.txt \
        …/K562_H3K27ac_REP1_peaks_SuperStitched_REGION_TO_GENE.txt
# → identical ✅
```

| Column | Match |
|---|---|
| OVERLAP_GENES | 1049/1049 byte-identical |
| PROXIMAL_GENES | 1049/1049 byte-identical |
| CLOSEST_GENE | 1049/1049 byte-identical (incl. all 4 ties) |

---

## How ROSE does it

[`ROSE/bin/ROSE_geneMapper.py::mapEnhancerToGene`](../../ROSE/bin/ROSE_geneMapper.py)
called with `-r` (`byRefseq=True`):

1. **Parse refseq table** → flat list, one row per RefSeq isoform.
2. **`makeStartDict`** — dict `refseq_id → first-occurrence's chrom/strand/TSS`.
   Only `refseqDict[name][0]` is ever indexed.
3. **`makeTranscriptCollection`** (`winSize=500`) — every refseq row's
   `txStart..txEnd` goes into a chr+strand-keyed bin index.
4. **`tssCollection`** (`winSize=50`) — one zero-length point per refseq id,
   using its `startDict` TSS.
5. **Per enhancer**:
   - `getOverlap` against transcript collection → **OVERLAP**
   - `getOverlap` against TSS-collection slopped ±50 kb → **PROXIMAL**
   - same against ±50 Mb → **distal** (effectively whole chrom)
   - prune `proximal -= overlap`, `distal -= proximal_pruned`
6. **CLOSEST** = `argmin |enh_center - startDict[id].tss|` over
   `overlap + proximal_pruned + distal_pruned`, first-min wins.

## Where ROSE is right / wrong

| Behavior | Verdict | Why |
|---|---|---|
| Coord-collision dedup (`Locus.__eq__` ignores ID) | **Wrong / silent** | `__hash__ = start+end`, `__eq__` compares chr/start/end/sense — distinct isoforms with identical coords silently collapse to one. |
| Per-refseq TSS = first occurrence only | Right-ish | RefSeq IDs are biologically unique; the alt-contig duplicates are mapping artifacts. Sensible default but loses the alt-contig copy. |
| Closest = enhancer-center → TSS | Defensible | Fine for narrow peaks, biased for wide enhancers. Modern tools use interval-to-interval distance. |
| Tie-breaking by bin-iteration order | **Effectively random** | Depends on `winSize`, file order, dict insertion. Not biologically meaningful — but reproducible because everyone uses the same `winSize`. |
| `geneList.count(line[1])` in inner loop | **Wrong/slow** | O(n²) over 58k refseqs (~5 min wallclock). Should be a `set`. |
| TSS bin=50, transcript bin=500 | Right-ish | TSS is a point; finer bins reduce wasted overlap checks. Different bin sizes → different tie-break order between OVERLAP and PROXIMAL — surprising. |

## Quirks that mattered for parity

To reach byte-identity we had to mirror three undocumented ROSE behaviors:

1. **Coord-collision dedup** of body and TSS tables by
   `(chrom, start, end, strand)`, keep-first.
2. **Per-id-first-occurrence** for the TSS collection (NOT per-row), matching
   the `makeStartDict → tssCollection` chain.
3. **Bin-iteration tie-breaking**: sort matches by
   `(strand_priority, max(locus.start, query.start) // winSize, file_order)` —
   `'+'` strand iterated before `'-'`, then by first bin shared with the query,
   then by file order. `winSize` is **500 for OVERLAP**, **50 for PROXIMAL/distal**.

## How Genomeblocks fixes / could improve

What we already have after this session:

- **`Genes.make_ucsc`** in [genomeblocks/genes.py](../genomeblocks/genes.py) —
  clean GTF-style parser for UCSC RefSeq, skipping alt-contig pollution.
- **[examples/rose_region_to_gene.py](rose_region_to_gene.py)** —
  byte-identical ROSE replication using pyranges + pandas. <10 s end-to-end.

Real improvements over ROSE the GB stack already enables (no new code):

1. **Strand-aware overlap as opt-in**, not always-on default — `Loci.overlaps`
   is strand-agnostic by default; ROSE conflates strands via `'both'`.
2. **Real interval-distance closest** via `pyranges.nearest`
   (`Genes.nearest_genes`) instead of center→TSS — handles wide enhancers
   correctly.
3. **No silent dedup.** `Loci` is a list-backed cgranges; nothing collapses
   behind your back.
4. **Promoter / 5'UTR / 3'UTR / exon / intron labels** out of
   `Genes._build_annot`, with priority semantics — ROSE has none of this.
5. **Reproducible tie-breaking** by TSS coordinate or refseq ID — deterministic,
   documented, not a bin-iteration accident.

Added to GB after this session ([genomeblocks/genes.py](../genomeblocks/genes.py)):

- `Genes.make_ucsc(..., keep_alt_contigs=True)` — keep alt-contig copies as
  separate `Gene` entries keyed `f"{gene_name}__{chrom}"` (default False skips
  them). Same-chrom transcript-id collisions (MHC paralogs) get `__N` suffix
  so nothing is silently dropped.
- `Genes.get_tss_transcripts()` / `Genes.nearest_transcripts(loci)` —
  alt-promoter-aware TSS map + nearest transcript per query (relevant for
  *TP53*, *CDKN2A*, *TCF7L2* and other multi-TSS genes).
- `Genes.enhancer_to_genes(loci, prox=50_000, level='gene'|'transcript')` —
  bundles overlap / proximal / closest into one call, returns a DataFrame.
  Unlike ROSE, `closest` uses real interval-to-interval distance
  (`pyranges.nearest`), not enhancer-center → TSS.

## Performance

| Pipeline | Wallclock | Bottleneck |
|---|---|---|
| ROSE | ~5 min | `geneList.count(line[1])` O(n²) Python loop |
| Genomeblocks (this script) | <10 s | pyranges interval-tree joins |
