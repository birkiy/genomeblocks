---
title: signal
parent: API Reference
layout: default
nav_order: 6
---

# `genomeblocks.signal`
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## `signal(loci, bigwigs, *, n_bins=200, flank=3_000, agg="mean", ...)`

```python
signal(
    loci: Loci,
    bigwigs: Sequence[str],
    *,
    n_bins: int = 200,
    flank: int = 3_000,
    agg: str = "mean",               # mean / max / min / std / sum / coverage
    dtype = np.float32,
    progress: bool = True,
    workers: int = 1,                # 1 = sequential; >1 = multiprocessing
    span: bool = False,              # True: use full locus span, not center±flank
    verbose: bool = True,
    backend: str | None = None,      # 'pybigtools' | 'bigwig' | None (auto)
    exact: bool = True,              # pybigtools base-accurate binning
) -> np.ndarray                      # shape (n_loci, n_tracks, n_bins)
```

Also attached as `Loci.signal(...)`.

The default path is **sequential** — a single native pass with `pybigtools`
is fastest for typical heatmap / browser / per-locus workloads. Pass
`workers > 1` for scale (many bigWigs × many loci): extraction then runs in a
`ProcessPoolExecutor` writing into a shared-memory cube. Multiprocessing (not
threading) is used because `pybigtools` serialises concurrent Python threads;
`workers` is capped at `min(workers, n_tracks·⌈n_loci/1000⌉, cpu_count()//2)`.

---

## `tmm(cube) -> np.ndarray`

Per-track TMM normalization + library-size-per-million scaling. Input/output shape is preserved.

---

## `plot_heatmap(loci, S, *, groups=None, sets=None, samples=None, ...)`

In `genomeblocks.signal_draw` (also attached as `Loci.plot_heatmap`). Full signature:

```python
plot_heatmap(
    loci: Loci,
    S: np.ndarray,                          # (regions, tracks, bins)
    *,
    groups: dict[str, Loci] | None = None,  # row groups (None = one group)
    sets: list[str] | None = None,          # row order
    samples: list[str] | None = None,       # column labels
    colors: dict[str, tuple] | None = None,
    ymax: float | list = 10,
    ymin: float | list = 0,
    height: int = 3000,                     # flank in bp — controls x-axis labels
    cmap: str | list = "Blues",
    vmax: float | list = 10,
    profile: bool = True,                   # top average profile
    sort: str | None = "group",             # 'group' | 'global' | None
    dpi: int = 100,
) -> matplotlib.figure.Figure
```

---

## `plot_profiles(loci, S, *, groups=None, sets=None, ...)`

In `genomeblocks.signal_draw` (also attached as `Loci.plot_profiles`).

```python
plot_profiles(
    loci,
    S,
    *,
    groups: dict[str, Loci] | None = None,
    sets: list[str] | None = None,
    colors: dict | None = None,
    ylim: float | None = None,
    dpi: int = 100,
    height: int = 3000,
) -> matplotlib.figure.Figure
```

---

## `compare_heatmap(a, b, bigwigs, ...)`

In `genomeblocks.signal_draw` (re-exported as `genomeblocks.compare_heatmap`).

```python
compare_heatmap(
    a: Loci, b: Loci,
    bigwigs: Sequence[str],
    *,
    a_name: str = "A",
    b_name: str = "B",
    common_name: str = "common",
    sets: list[str] | None = None,
    samples: list[str] | dict[str, list[int]] | None = None,
    n_bins: int = 200,
    flank: int = 3_000,
    agg: str = "mean",
    normalize: bool = True,       # run tmm() before plotting
    cmap: str | list = "Blues",
    vmax: float | list = 10,
    ymax: float | list = 10,
    ymin: float | list = 0,
    profile: bool = True,
    sort: str | None = "group",
    colors: dict | None = None,
    dpi: int = 100,
    S: np.ndarray | None = None,  # pre-computed signal cube for the union
    signal_kw: dict | None = None,
) -> (fig, union_loci, S, groups)
```

Computes `a - b`, `a & b`, `b - a`; stacks into a `groups`-grouped heatmap.

---

## Backend helpers

| Symbol | Purpose |
|---|---|
| `_resolve_opener(backend)` | `(opener, backend_name)` for `'pybigtools'` / `'bigwig'` / `None` — pure, no global mutation. |
| `_open_pybigtools(path)` | Adapter returning a `_PyBigToolsHandle`. |
| `_PyBigToolsHandle` | Thin wrapper matching the internal reader interface (`chroms()`, `stats(...)`, `values(...)`, `close()`). |
| `_detect_backend()` | `(opener, backend_name)` — prefers pybigtools, falls back to pure-Python. |
