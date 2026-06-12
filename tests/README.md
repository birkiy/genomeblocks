# Tests

Fast, synthetic unit tests — no external genome / bigwig / ChIP-Atlas data
required. All fixtures are tiny in-memory objects or temp files (see
[`conftest.py`](conftest.py)).

## Running

```bash
pip install -e ".[test]"     # pytest
pytest                        # from the repo root
```

graph-tool must be importable (it is the one conda-only dependency), so run in an
environment that has it.

## Coverage

| file | module under test |
| --- | --- |
| `test_locus.py` | `Locus` — uid, distance, overlap, ordering |
| `test_loci.py` | `Loci` — `make`, set algebra (`& - + \| ^`), sort/merge/slop, overlap, nearest |
| `test_bedpe.py` | `read_bedpe`, `Pair`, `Loci.pair_to_bed` |
| `test_genes.py` | GTF parsing, region-class annotation, nearest gene |
| `test_architecture.py` | build, `strength` / `elbow` / `annotate` / `prime_hubs`, `normalize`, subgraph/copy/set-ops/pickle |
| `test_signal.py` | `tmm`, `_bcast`, `_resolve_groups`, `_even_ranges`, `plan_workers` |
| `test_signal_draw.py` | `plot_heatmap` / `plot_profiles` with `groups` |
| `test_browser.py` | track detection, bigwig-list averaging, `bw_share`, scalar `bw_ymax` |
| `test_motifs.py` | `bootstrap_enrichment` LFC |
| `test_atlas.py` | `Atlas.make` / `search`, header-less `attach_meta` |
| `test_api.py` | public surface; removed names (Tags, make_spread, …) stay gone; import order |
