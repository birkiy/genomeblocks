Genomeblocks
============

Fluent building blocks for regulatory genomics, as columnar tables.

- `Loci`: interval tables with set algebra (`& | - ^`), `slop` / `sort` / `merge` / `nearest`, indexed lookups and signal extraction.
- `Genes`: GTF / GFF3 / UCSC gene models as three linked tables (genes, transcripts, features); annotation and ATAC-supported isoforms.
- `Pairs` and `Architecture`: BEDPE loops and the chromatin-contact graph over CREs — Hi-C weights from `.mcool`, O/E normalisation, hubs.
- `signal`: bigWig signal cubes, TMM normalisation, heatmaps; `motifs`: motif scanning with any engine; `Atlas`: GIGGLE-style enrichment.
- `browser`, `igv_html`, `View`: a matplotlib region view, a shareable IGV page and a one-file interactive browser.

Every table hands itself to pandas, polars, Arrow, bioframe, pyranges, pybedtools, AnnData, duckdb, seaborn and friends, and takes their frames back. The heavy work runs through swappable backends (numpy / cgranges / ncls / bioframe / pyranges / bedtools for intervals, pybigtools / pyBigWig / pure Python for bigWigs, lightmotif / MOODS / Biopython for motifs, graph-tool / scipy / igraph / networkx for graphs, ...) with the same answer whichever engine runs.

📖 **Documentation**: [birkiy.github.io/genomeblocks](https://birkiy.github.io/genomeblocks/) — start with the [quickstart](https://birkiy.github.io/genomeblocks/quickstart/).

Quick Start
-----------

```python
import genomeblocks as gb

cre = (gb.Loci.make("atac.narrowPeak")          # a table: numpy columns, row = join key
         .slop(100).sort().merge())

se = cre & gb.as_loci("H3K27ac_SE.bed")          # anything interval-like goes in

genes = gb.Genes.make("gencode.v38.annotation.gtf")   # 0-based tables, GTF or GFF3
labels = genes.annotations(se)                   # Promoter-TSS / Exonic / Intronic / ...

A = (gb.Architecture.make(cre, "loops.bedpe", r=2500)
       .add_mcool("hic.mcool", resolution=5000)
       .normalize().annotate(genes).strength())
hubs = A.prime_hubs()

S = cre.signal(["atac.bw", "h3k27ac.bw"], n_bins=200, flank=3000)   # (rows, tracks, bins)
cre.to_polars(); A.to_pandas(); genes.to_arrow()                    # and back: Loci.from_frame(df)
cre.intersect(se, backend="bioframe")             # or: with gb.use_backend(intervals="pyranges"): ...
```

Setup
-----

```bash
pip install genomeblocks                 # every default backend
pip install "genomeblocks[all]"          # + polars, pysam, logomaker and every other engine
```

The conda-only engines (graph-tool, cgranges, bedtools) are optional extras; `environment.yml` builds the full development environment:

```bash
conda env create -f environment.yml
conda activate genomeblocks
```

`gb.backends()` lists what is installed and what each family uses. See [the documentation](https://birkiy.github.io/genomeblocks/) for installation details, the concepts, module guides, interoperability, backends and the API reference.

License
-------

MIT — see [LICENSE](LICENSE).
