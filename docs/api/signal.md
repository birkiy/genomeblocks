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
    max_bw_parallel: int = 6,
    max_workers: int | None = None,
    span: bool = False,              # True: use full locus span, not center±flank
    verbose: bool = True,
    backend: str | None = None,      # 'pybigtools' | 'bigwig' | None (auto)
) -> np.ndarray                      # shape (n_loci, n_tracks, n_bins)
```

Also attached as `Loci.signal(...)`.

---

## `tmm(cube) -> np.ndarray`

Per-track TMM normalization + library-size-per-million scaling. Input/output shape is preserved.

---

## `plot_heatmap(loci, S, *, tags=None, groups=None, sets=None, samples=None, ...)`

Full signature:

```python
plot_heatmap(
    loci: Loci,
    S: np.ndarray,                          # (regions, tracks, bins)
    *,
    tags: Tags | None = None,
    groups: dict[str, Loci] | None = None,  # alternative to tags
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

## `plot_profiles(loci, S, *, tags=None, sets=None, ...)`

```python
plot_profiles(
    loci,
    S,
    *,
    tags: Tags | None = None,
    sets: list[str] | None = None,
    colors: dict | None = None,
    ylim: float | None = None,
    dpi: int = 100,
    height: int = 3000,
) -> matplotlib.figure.Figure
```

---

## `compare_heatmap(a, b, bigwigs, ...)`

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
) -> (fig, union_loci, S, tags)
```

Computes `a - b`, `a & b`, `b - a`; stacks into a Tags-grouped heatmap.

---

## `plan_workers(n_tracks, n_loci, *, cores=None, max_bw_parallel=6)`

Returns a list of `((t_lo, t_hi), (l_lo, l_hi))` work chunks used internally by `signal()`. Useful for custom extraction pipelines.

---

## Backend helpers

| Symbol | Purpose |
|---|---|
| `_open_pybigtools(path)` | Adapter returning a `_PyBigToolsHandle`. |
| `_PyBigToolsHandle` | Thin wrapper matching the internal reader interface (`chroms()`, `stats(...)`, `values(...)`, `close()`). |
| `_detect_backend()` | `(opener, backend_name)` — prefers pybigtools, falls back to pure-Python. |
