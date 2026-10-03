# Columnar genomeblocks (prototype)

A notebook tour of `genomeblocks.columnar`, the table-first prototype on the
`columnar-prototype` branch:

- one shared `Genome`; CREs, genes and the Architecture are numpy tables
- the row number links the tables: row `i` of the CREs is row `i` of the labels,
  the signal cube and every graph column
- the Architecture keeps all edges in one table, sorted into per-chromosome cis
  blocks plus one trans block, so `A.chrom("chr8")`, `A.cis` and `A.trans` are
  zero-copy views, and neighbours and graph-tool algorithms still see every edge

| file | what |
| --- | --- |
| `columnar_prototype.ipynb` | the tour, executed (outputs included) |
| `make_notebook.py` | rebuilds and re-executes the notebook |

## Data

The notebook uses the synthetic hg38-shaped benchmark data:

```bash
cd benchmarks
python make_data.py chromsizes peaks bw gtf pairs loops loops_trans hic_trans
```

`hic_trans` needs the `cooler` CLI. Then run `python examples/columnar_prototype/make_notebook.py`.
