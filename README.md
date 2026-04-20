Genomeblocks
============

Fluent building blocks for regulatory genomics.

- `Loci`: craft and manipulate candidate regulatory element (CRE) sets.
- `Tags`: attach boolean / numeric annotations and query them with expressions.
- `Genes`: parse GTF / UCSC RefSeq and annotate CREs with promoter-anchored gene models.
- `Architecture`: chromatin-contact graphs over CREs — build from BEDPE loops, overlay mcool matrices, discover hubs and focus genes.
- `signal`: threaded bigWig extraction, TMM normalization, comparative heatmaps.
- `browser`: IGV-like, fully-vectorial region viewer.

📖 **Documentation**: [birkiy.github.io/genomeblocks](https://birkiy.github.io/genomeblocks/)

Quick Start
-----------

```python
from genomeblocks import Architecture, Genes, Loci

cre = (Loci.make("atac.narrowPeak")
           .slop(100)
           .sort()
           .merge())

se = cre.intersect(Loci.make("H3K27ac_SE.bed"))

arch = (Architecture.make(cre, "RNAP_loops.bedpe", r=2500)
                    .add_mcool(cre, "RNAP.mcool", resolution=5000)
                    .normalize(cre))

genes  = Genes.make("gencode.v38.annotation.gtf", promoter_r=1000)
counts = genes.annotations(se & cre).groupby("annotation").size()
```

Setup
-----

```bash
conda env create -f environment.yml
conda activate genomeblocks
pip install -e .
```

For a pip-only install (no `graph-tool` → no `Architecture`):

```bash
pip install genomeblocks
```

See [the documentation](https://birkiy.github.io/genomeblocks/) for installation details, a full quickstart, module guides, and an API reference.

License
-------

MIT — see [LICENSE](LICENSE).
