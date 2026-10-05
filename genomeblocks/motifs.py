"""Motifs: scanning, enrichment and archetypes.

Scanning reads each locus window (centre ± ``r``) from a FASTA (path, dict or
open pyfaidx / pysam / Biopython handle — the fasta backend) and scores a
motif library (JASPAR / TRANSFAC / uniprobe / MEME file, Biopython motifs, or
plain matrices — see :func:`~genomeblocks.load_motifs`) with the motifs
backend: MOODS by default (lightmotif when only it is installed), Biopython on request, all scoring
the same log-odds matrices so they report the same hits.

Results line up with the loci: row ``i`` of a matrix or cube is locus ``i``.
Windows that run off a chromosome end, sit on a chromosome missing from the
FASTA, or contain letters other than A C G T N count 0 (as a bigWig cube
leaves such rows at 0).
"""
from __future__ import annotations

import os
from typing import TYPE_CHECKING, Dict, List, Optional

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

from .backends.motifs import Block, Library, load_motifs, write_meme  # noqa: F401  (re-exported)


# ── shared setup ───────────────────────────────────────────────────────────

def _windows(loci, fasta, r: int, backend_fasta=None):
    """(Loci, upper-case window sequences, valid mask): windows are centre ± r."""
    from .interop import as_loci
    L = as_loci(loci)
    seqs = L.sequences(fasta, r=r, upper=True, backend=backend_fasta)
    valid = np.fromiter((len(q) == 2 * r for q in seqs), bool, len(seqs))
    if not valid.all():
        seqs = [q if ok else "" for q, ok in zip(seqs, valid)]
    return L, seqs, valid


def _thresholds(lib: Library, threshold, pvalue):
    """Per-motif score cutoffs: ``threshold`` (scalar or one per motif) or, with
    ``pvalue``, the score each motif reaches with that probability under a
    uniform background (exact score distribution on a 0.001-bit grid, the
    same whichever engine scans)."""
    if pvalue is None:
        return np.broadcast_to(np.asarray(threshold, np.float64), (len(lib),)).copy()
    from .backends.motifs import threshold_from_pvalue
    return np.array([threshold_from_pvalue(lib.logodds(i), pvalue) for i in range(len(lib))], np.float64)


def _workers(workers, n_motifs, n_seqs):
    if workers is None:
        workers = max(1, (os.cpu_count() or 2) // 2)
    if workers > 1 and (n_motifs < 16 or n_seqs < 200):     # process start-up would dominate
        workers = 1
    return workers


_STATE: dict = {}


def _init(seqs, mats, thr, both, backend):
    _STATE.update(block=Block(seqs, backend), mats=mats, thr=thr, both=both)


def _counts_task(idx):
    st = _STATE
    blk = st["block"]
    if blk.backend == "moods":
        return idx, blk.moods_counts([st["mats"][i] for i in idx], st["thr"][idx], st["both"])
    return idx, np.column_stack([blk.counts(st["mats"][i], st["thr"][i], st["both"]) for i in idx])


def _count_matrix(seqs, lib, thr, both, backend, workers, verbose, desc):
    """(windows x motifs) hit counts, motifs spread over a process pool."""
    from tqdm import tqdm
    mats = [lib.logodds(i) for i in range(len(lib))]
    n_m = len(mats)
    out = np.zeros((len(seqs), n_m), np.int64)
    workers = _workers(workers, n_m, len(seqs))
    chunk = max(1, n_m // (workers * 4)) if workers > 1 else max(1, min(64, n_m))
    batches = [np.arange(k, min(k + chunk, n_m)) for k in range(0, n_m, chunk)]
    if workers == 1:
        _init(seqs, mats, thr, both, backend)
        it = map(_counts_task, batches)
        if verbose:
            it = tqdm(it, total=len(batches), desc=desc)
        for idx, block in it:
            out[:, idx] = block
        _STATE.clear()
        return out
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers, initializer=_init,
                             initargs=(seqs, mats, thr, both, backend)) as ex:
        it = ex.map(_counts_task, batches)
        if verbose:
            it = tqdm(it, total=len(batches), desc=desc)
        for idx, block in it:
            out[:, idx] = block
    return out


def _report(verbose, valid, tag):
    if verbose and not valid.all():
        print(f"[{tag}] {int((~valid).sum())} of {len(valid)} windows run off a chromosome or "
              f"are missing from the FASTA; their rows count 0.")


# ── scanning ───────────────────────────────────────────────────────────────

def scan_motifs(loci, fasta, motifs, *, format: str = "jaspar", r: int = 250, threshold=13.0,
                pvalue: Optional[float] = None, norm: bool = True, both_strands: bool = False,
                backend: Optional[str] = None, verbose: bool = True) -> Dict[str, float]:
    """Total hits of every motif over all windows (divided by motif width when ``norm``)."""
    M = scan_motifs_matrix(loci, fasta, motifs, format=format, r=r, threshold=threshold, pvalue=pvalue,
                           norm=norm, both_strands=both_strands, backend=backend, workers=1, verbose=verbose)
    return {k: float(v) for k, v in M.sum(axis=0).items()}


def scan_motifs_matrix(loci, fasta, motifs, *, format: str = "jaspar", r: int = 250, threshold=13.0,
                       pvalue: Optional[float] = None, norm: bool = True, both_strands: bool = False,
                       workers: Optional[int] = None, backend: Optional[str] = None,
                       fasta_backend: Optional[str] = None, verbose: bool = True):
    """Hits of every motif in every window: a DataFrame (one row per locus, index
    = uid, in the loci's order; one column per motif).

    Args:
        loci: anything :func:`~genomeblocks.as_loci` takes.
        fasta: FASTA path, ``{chrom: str}`` dict, or an open pyfaidx / pysam /
            Biopython handle.
        motifs: a motif file (``format`` = 'jaspar', 'jaspar16', 'transfac',
            'uniprobe', 'meme'), a Library, Biopython motifs, or matrices.
        r: half-window around each locus centre (windows are 2r bp).
        threshold: log2-odds cutoff (scalar or one per motif); or give
            ``pvalue`` for a per-motif cutoff with that match probability.
        norm: divide counts by motif width.
        both_strands: also count reverse-strand matches.
        workers: processes over motifs (None = half the cores; small jobs stay serial).
        backend: 'moods' (default), 'lightmotif' or 'biopython' — the same hits.
    """
    import pandas as pd
    lib = load_motifs(motifs, format=format)
    L, seqs, valid = _windows(loci, fasta, r, fasta_backend)
    _report(verbose, valid, "motifs")
    thr = _thresholds(lib, threshold, pvalue)
    counts = _count_matrix(seqs, lib, thr, both_strands, backend, workers, verbose, "[motifs]")
    vals = counts / lib.widths[None, :] if norm else counts
    return pd.DataFrame(vals, index=pd.Index(L.uid, name="uid"), columns=lib.names)


def scan_motifs_matrix_masked(loci, fasta, motifs, anchors: List[str], *, format: str = "jaspar", r: int = 250,
                              window: int = 10, threshold=13.0, pvalue: Optional[float] = None,
                              anchor_threshold: Optional[float] = None,
                              norm: bool = True, both_strands: bool = False, skip_anchors: bool = True,
                              seed: Optional[int] = None, workers: Optional[int] = None,
                              backend: Optional[str] = None, verbose: bool = True):
    """:func:`scan_motifs_matrix` after hiding the matches of ``anchors`` motifs.

    Every match of a motif whose name contains one of ``anchors``
    (case-insensitive) gets ``centre ± window`` bp replaced by random bases,
    then the library is scanned on the masked windows — "which motifs set
    these regions apart *besides* the anchor?" (e.g. mask CTCF). Anchor
    motifs are left out of the result unless ``skip_anchors=False``.
    """
    import pandas as pd
    lib = load_motifs(motifs, format=format)
    anchor_idx = lib.indices(anchors, "substring")
    if not anchor_idx:
        raise ValueError(f"No motifs matched anchors {anchors!r}")
    if verbose:
        names = [lib.names[i] for i in anchor_idx]
        print(f"[mask] {len(names)} anchor PSSM(s): {', '.join(names[:5])}"
              f"{f' (+{len(names) - 5} more)' if len(names) > 5 else ''}")
    L, seqs, valid = _windows(loci, fasta, r)
    _report(verbose, valid, "mask")
    thr_all = _thresholds(lib, threshold, pvalue)           # one cutoff per library motif
    a_thr = thr_all if anchor_threshold is None else np.full(len(lib), float(anchor_threshold))
    blk = Block(seqs, backend)
    rng = np.random.default_rng(seed)
    bases = np.frombuffer(b"ACGT", np.uint8)
    masked = [bytearray(q.encode("ascii")) for q in seqs]
    for i in anchor_idx:
        m = lib.logodds(i)
        rows, pos, _ = blk.hits(m, a_thr[i], both_strands)
        c = pos + len(m) // 2
        for row, cc in zip(rows.tolist(), c.tolist()):
            seq = masked[row]
            lo, hi = max(0, cc - window), min(len(seq), cc + window + 1)
            if hi > lo:
                seq[lo:hi] = bases[rng.integers(0, 4, size=hi - lo)].tobytes()
    masked = [b.decode("ascii") for b in masked]
    keep = [i for i in range(len(lib)) if not (skip_anchors and i in set(anchor_idx))]
    sub = lib.take(keep)
    thr = np.asarray(thr_all, np.float64)[keep]
    counts = _count_matrix(masked, sub, thr, both_strands, backend, workers, verbose, "[scan]")
    vals = counts / sub.widths[None, :] if norm else counts
    return pd.DataFrame(vals, index=pd.Index(L.uid, name="uid"), columns=sub.names)


# ── positional profiles (where in the window the hits fall) ───────────────

def _profile_task(idx):
    st = _STATE
    blk, n_bins, L, norm = st["block"], st["n_bins"], st["L"], st["norm"]
    out = []
    for i in idx:
        m = st["mats"][i]
        rows, pos, _ = blk.hits(m, st["thr"][i], st["both"])
        prof = np.zeros((blk.n, n_bins), np.float32)
        b = ((pos + len(m) // 2) * n_bins) // L
        ok = (b >= 0) & (b < n_bins)
        np.add.at(prof, (rows[ok], b[ok]), 1.0)
        if norm:
            prof /= len(m)
        out.append((i, prof))
    return out


def _init_profile(seqs, mats, thr, both, backend, n_bins, L, norm):
    _init(seqs, mats, thr, both, backend)
    _STATE.update(n_bins=n_bins, L=L, norm=norm)


def scan_motifs_profile(loci, fasta, motifs, select=None, *, format: str = "jaspar", match: str = "substring",
                        r: int = 500, n_bins: int = 100, threshold=13.0, pvalue: Optional[float] = None,
                        both_strands: bool = True, smooth: float = 1.0, norm: bool = False,
                        workers: Optional[int] = None, backend: Optional[str] = None, verbose: bool = True):
    """Where the motif hits fall in each window: a ``(rows, motifs, bins)`` cube —
    the layout of a signal cube, so ``signal_draw.plot_heatmap`` (or
    :func:`~genomeblocks.plot_motif_heatmap`) draws it like one.

    Each window (centre ± ``r``) is split into ``n_bins`` bins and every hit
    is credited to the bin of its centre; both strands by default (reverse
    matches in forward coordinates, so they share the frame). ``select``
    picks motifs by name (case-insensitive substring unless
    ``match='exact'``) — usually a handful, e.g. ``['CTCF', 'GATA', 'SOX']``.
    ``smooth`` is a Gaussian sigma in bins along each row (0 keeps raw
    counts); ``norm`` divides by motif width.

    Returns:
        (M, names): the float32 cube and the motif names of its middle axis.
    """
    lib = load_motifs(motifs, format=format)
    idx = lib.indices(select, match)
    if not idx:
        raise ValueError(f"No motifs matched {select!r} (match={match!r}); check the names exist in the library.")
    lib = lib.take(idx)
    Lc, seqs, valid = _windows(loci, fasta, r)
    _report(verbose, valid, "profile")
    thr = _thresholds(lib, threshold, pvalue)
    mats = [lib.logodds(i) for i in range(len(lib))]
    if verbose:
        head = ", ".join(lib.names[:6]) + (f" (+{len(lib) - 6} more)" if len(lib) > 6 else "")
        print(f"[profile] {len(Lc)} loci x {len(lib)} motifs x {n_bins} bins | {2 * r} bp window, "
              f"{2 * r // n_bins} bp/bin: {head}")
    workers = _workers(workers, len(lib) * 4, len(seqs))
    batches = [list(range(k, min(k + 1, len(lib)))) for k in range(len(lib))]
    M = np.zeros((len(Lc), len(lib), n_bins), np.float32)
    args = (seqs, mats, thr, both_strands, backend, n_bins, 2 * r, norm)
    if workers == 1:
        _init_profile(*args)
        for b in batches:
            for i, prof in _profile_task(b):
                M[:, i, :] = prof
        _STATE.clear()
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(workers, len(batches)), initializer=_init_profile,
                                 initargs=args) as ex:
            for out in ex.map(_profile_task, batches):
                for i, prof in out:
                    M[:, i, :] = prof
    if smooth and smooth > 0:
        from scipy.ndimage import gaussian_filter1d
        M = gaussian_filter1d(M, sigma=float(smooth), axis=-1, mode="nearest")
    if verbose:
        per = M.sum(axis=(0, 2)) / max(len(Lc), 1)
        print("[profile] mean hits/locus: " + ", ".join(f"{n}={v:.2f}" for n, v in list(zip(lib.names, per))[:6]))
    return M, list(lib.names)


def bootstrap_enrichment(
    groups: Dict[str, "pd.DataFrame"],
    ref: "pd.DataFrame",
    *,
    boot: int = 100,
    sample: int = 500,
    pseudo: float = 0.1,
    seed: Optional[int] = None,
    verbose=True
):
    """Bootstrap mean motif counts from precomputed matrices and compute LFC.

    Parameters
    ----------
    groups : dict[name -> DataFrame]
        Precomputed (n_loci x n_motifs) matrices from ``scan_motifs_matrix``.
    ref : DataFrame
        Reference matrix (e.g. all CREs) sampled as the background pool.
    boot, sample : int
        Number of bootstrap iterations and per-iteration sample size.
    seed : int, optional
        RNG seed for reproducibility.

    Returns
    -------
    pandas.DataFrame with one row per motif. Columns:
        Factor, mean_ref, and for each group: LFC_<g>, mean_<g>. If exactly
        two groups are supplied, also emits LFC = LFC_<first> - LFC_<second>.
    """
    import numpy as np
    import pandas as pd
    from tqdm import tqdm


    rng = np.random.default_rng(seed)

    # Align columns across all matrices so the bootstrap math is vectorized.
    all_motifs = set(ref.columns)
    for df in groups.values():
        all_motifs.update(df.columns)
    cols = sorted(all_motifs)
    ref_aligned = ref.reindex(columns=cols, fill_value=0.0)
    groups_aligned = {k: v.reindex(columns=cols, fill_value=0.0) for k, v in groups.items()}

    def _boot_mean(mat: pd.DataFrame, n_sample: int) -> np.ndarray:
        arr = mat.to_numpy(copy=False, dtype=float)
        n_rows = arr.shape[0]
        if n_rows == 0:
            return np.zeros(arr.shape[1])
        k = min(n_sample, n_rows)
        iter_means = np.empty((boot, arr.shape[1]))
        iterator = range(boot)
        if verbose: iterator = tqdm(iterator, total=boot)

        for i in iterator:
            idx = rng.choice(n_rows, size=k, replace=False)
            iter_means[i] = arr[idx].mean(axis=0)
        return iter_means.mean(axis=0)

    mean_ref = _boot_mean(ref_aligned, sample)
    out = {'Factor': cols, 'mean_ref': mean_ref}
    lfc_by_group = {}
    for name, df in groups_aligned.items():
        m = _boot_mean(df, sample)
        lfc = np.log2((pseudo + m) / (pseudo + mean_ref))
        out[f'mean_{name}'] = m
        out[f'LFC_{name}'] = lfc
        lfc_by_group[name] = lfc

    if len(lfc_by_group) == 2:
        a, b = list(lfc_by_group.keys())
        out['LFC'] = lfc_by_group[a] - lfc_by_group[b]

    return pd.DataFrame(out)




def compare_motifs(
    mat_a,
    mat_b,
    *,
    pseudo: float = 0.1,
    alternative: str = 'two-sided',
    verbose: bool = True,
):
    """Per-motif two-sample differential test between A and B matrices.

    For each motif column, runs a Mann-Whitney U test on the per-sequence
    counts (rank-based, no normality assumption — appropriate for the
    zero-inflated count distributions ``scan_motifs_matrix`` produces) plus
    a log2 fold change of the column means.

    Parameters
    ----------
    mat_a, mat_b : DataFrame
        Outputs of :func:`scan_motifs_matrix` or
        :func:`scan_motifs_matrix_masked`. Rows are sequences, columns are
        motifs. Missing columns are filled with zero so the two matrices
        don't need identical motif sets.
    pseudo : float
        Pseudocount added to both means before the log2 ratio. Set this to
        roughly the median non-zero per-CRE mean to avoid the +1 collapse
        we discussed for ``bootstrap_enrichment``.
    alternative : {'two-sided', 'greater', 'less'}
        Passed to scipy.stats.mannwhitneyu. Default two-sided.

    Returns
    -------
    pandas.DataFrame with columns:
        Factor, mean_A, mean_B, LFC, U, p, p_adj (Benjamini-Hochberg).
        Sorted by LFC descending.
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import mannwhitneyu

    cols = sorted(set(mat_a.columns) | set(mat_b.columns))
    A = mat_a.reindex(columns=cols, fill_value=0.0).to_numpy(dtype=float)
    B = mat_b.reindex(columns=cols, fill_value=0.0).to_numpy(dtype=float)
    mean_A = A.mean(axis=0)
    mean_B = B.mean(axis=0)
    lfc = np.log2((pseudo + mean_A) / (pseudo + mean_B))

    n = len(cols)
    U = np.zeros(n)
    p = np.ones(n)
    for j in range(n):
        a_col, b_col = A[:, j], B[:, j]
        if a_col.max() == 0 and b_col.max() == 0:
            continue
        try:
            res = mannwhitneyu(a_col, b_col, alternative=alternative)
            U[j], p[j] = res.statistic, res.pvalue
        except ValueError:
            pass

    # Benjamini-Hochberg FDR
    o = np.argsort(p)
    p_adj = np.empty_like(p)
    p_adj[o] = np.minimum.accumulate(
        (p[o] * n / (np.arange(n) + 1))[::-1])[::-1].clip(0, 1)

    out = pd.DataFrame({
        'Factor': cols,
        'mean_A': mean_A, 'mean_B': mean_B,
        'LFC': lfc, 'U': U, 'p': p, 'p_adj': p_adj,
    }).sort_values('LFC', ascending=False).reset_index(drop=True)

    if verbose:
        sig = int((out['p_adj'] < 0.01).sum())
        print(f"[compare_motifs] {n} motifs | sig (p_adj<0.01): {sig}")
    return out


def compare_motifs_to_ref(
    query,
    ref,
    *,
    pseudo: float = 0.1,
    alternative: str = 'two-sided',
    verbose: bool = True,
):
    """Per-motif comparison of one or more query sets against a reference pool.

    For each query group, runs a Mann-Whitney U test of per-sequence motif
    counts (query vs ref) and reports the rank-based p-value alongside a
    log2 fold change of the means. Use when you want "is this set enriched
    above background?" with per-motif significance — unlike
    :func:`bootstrap_enrichment`, which only reports point estimates.

    Parameters
    ----------
    query : DataFrame or dict[name -> DataFrame]
        Single query matrix or a name→matrix dict for multi-group comparisons.
        Rows are sequences, columns are motifs (output of
        :func:`scan_motifs_matrix` or :func:`scan_motifs_matrix_masked`).
    ref : DataFrame
        Reference pool matrix (e.g. all CREs). **If the queries were masked,
        mask the reference with the same anchors** — otherwise queries look
        depleted across the board for reasons unrelated to biology.
    pseudo : float
        Pseudocount for the log2 fold change. Set ≈ median non-zero column
        mean to avoid the pseudocount dominating tiny rates.
    alternative : {'two-sided', 'greater', 'less'}
        Passed through to scipy.stats.mannwhitneyu.

    Returns
    -------
    pandas.DataFrame with one row per motif. Columns:
        Factor, mean_ref, and for each group g:
        mean_<g>, LFC_<g>, U_<g>, p_<g>, p_adj_<g> (BH within group).
        When exactly two groups are supplied, also emits
        LFC = LFC_<first> - LFC_<second>, matching
        :func:`bootstrap_enrichment`'s convention.
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import mannwhitneyu

    single = not isinstance(query, dict)
    groups = {'query': query} if single else dict(query)
    if not groups:
        raise ValueError("`query` must contain at least one group.")

    all_cols = set(ref.columns)
    for df in groups.values():
        all_cols.update(df.columns)
    cols = sorted(all_cols)
    R = ref.reindex(columns=cols, fill_value=0.0).to_numpy(dtype=float)
    G_aligned = {name: df.reindex(columns=cols, fill_value=0.0).to_numpy(dtype=float)
                 for name, df in groups.items()}

    n = len(cols)
    mean_ref = R.mean(axis=0)
    out = {'Factor': cols, 'mean_ref': mean_ref}
    lfc_by_group: Dict[str, "np.ndarray"] = {}

    for name, G in G_aligned.items():
        mean_g = G.mean(axis=0)
        lfc = np.log2((pseudo + mean_g) / (pseudo + mean_ref))
        U = np.zeros(n)
        p = np.ones(n)
        for j in range(n):
            g_col, r_col = G[:, j], R[:, j]
            if g_col.max() == 0 and r_col.max() == 0:
                continue
            try:
                res = mannwhitneyu(g_col, r_col, alternative=alternative)
                U[j], p[j] = res.statistic, res.pvalue
            except ValueError:
                pass
        # BH FDR within each group
        o = np.argsort(p)
        p_adj = np.empty_like(p)
        p_adj[o] = np.minimum.accumulate(
            (p[o] * n / (np.arange(n) + 1))[::-1])[::-1].clip(0, 1)

        out[f'mean_{name}']  = mean_g
        out[f'LFC_{name}']   = lfc
        out[f'U_{name}']     = U
        out[f'p_{name}']     = p
        out[f'p_adj_{name}'] = p_adj
        lfc_by_group[name] = lfc

        if verbose:
            sig = int((p_adj < 0.05).sum())
            print(f"[compare_motifs_to_ref] {name}: {n} motifs | "
                  f"sig (p_adj<0.05): {sig}")

    if len(lfc_by_group) == 2:
        a, b = list(lfc_by_group.keys())
        out['LFC'] = lfc_by_group[a] - lfc_by_group[b]

    df = pd.DataFrame(out)
    sort_col = 'LFC' if 'LFC' in df.columns else f'LFC_{next(iter(groups))}'
    df = df.sort_values(sort_col, ascending=False).reset_index(drop=True)
    return df




# ============================================================================
# Motif clustering and archetype generation
# ----------------------------------------------------------------------------
# Pairwise Sandelin-Wasserman PWM similarity → hierarchical clustering →
# consensus archetype PFM per cluster. The recipe is the same one used by
# STAMP, RSAT matrix-clustering, and Vierstra et al. 2020 (the Nature
# archetype paper) — implemented natively in numpy/scipy so it composes with
# the rest of motifs.py without external tool dependencies.
# Logos are rendered with `logomaker` (lazy import inside the plot helpers).
# ============================================================================



def _pwm_rc(pfm):
    """Reverse-complement of a (4, W) PFM (rows ACGT). Reversing both axes
    swaps A↔T / C↔G and reads the sequence 3'→5'."""
    return pfm[::-1, ::-1]


def _pfm_information_content(p):
    """Per-column information content (bits) of a (4, W) PFM."""
    import numpy as np
    safe = np.where(p > 0, p, 1.0)
    return 2.0 + (p * np.log2(safe)).sum(axis=0)


def _pwm_align(p_anchor, p_other, *, min_overlap=5):
    """Best alignment of `p_other` against `p_anchor`.

    Returns (offset, is_rc, similarity), where `offset` is the position of
    p_other's first column in p_anchor's coordinate frame (may be negative)
    and similarity is the Sandelin-Wasserman score normalized to [0, 1].
    """
    W1, W2 = p_anchor.shape[1], p_other.shape[1]
    candidates = ((False, p_other), (True, _pwm_rc(p_other)))
    best = (0, False, -1.0)
    for is_rc, p2 in candidates:
        # Offset range: enough overlap on each side.
        for o in range(-(W2 - min_overlap), W1 - min_overlap + 1):
            a = max(0, o)
            b = min(W1, o + W2)
            n = b - a
            if n < min_overlap:
                continue
            c = a - o
            diff = p_anchor[:, a:b] - p2[:, c:c + n]
            sim = (2.0 - (diff * diff).sum(axis=0)).sum() / (2.0 * n)
            if sim > best[2]:
                best = (o, is_rc, sim)
    return best


def _pwm_similarity(p1, p2, *, min_overlap=5):
    """Just the best-of-all-offsets/strands SW similarity in [0, 1]."""
    return _pwm_align(p1, p2, min_overlap=min_overlap)[2]


# ── Parallel pairwise distance computation ────────────────────────────────
_PWM_STATE: dict = {}


def _pwm_worker_init(pwms, min_overlap):
    _PWM_STATE['pwms'] = pwms
    _PWM_STATE['min_overlap'] = min_overlap


def _pwm_pairs_chunk(pairs):
    pwms = _PWM_STATE['pwms']
    min_overlap = _PWM_STATE['min_overlap']
    out = []
    for i, j in pairs:
        sim = _pwm_similarity(pwms[i], pwms[j], min_overlap=min_overlap)
        out.append((i, j, 1.0 - sim))
    return out


def _pwm_distance_matrix(pwms, *, min_overlap=5, workers=None, verbose=True):
    """Compute an (n × n) SW distance matrix from a list of PFMs."""
    import numpy as np
    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm

    n = len(pwms)
    D = np.zeros((n, n), dtype=float)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    if workers is None:
        workers = max(1, (os.cpu_count() or 2) // 2)
    if workers > 1 and len(pairs) < 1000:
        workers = 1

    if workers == 1:
        _pwm_worker_init(pwms, min_overlap)
        it = pairs
        if verbose:
            it = tqdm(it, total=len(pairs), desc='[pwm-dist]')
        for i, j in it:
            D[i, j] = D[j, i] = 1.0 - _pwm_similarity(
                pwms[i], pwms[j], min_overlap=min_overlap)
    else:
        chunk = max(1, len(pairs) // (workers * 4))
        batches = [pairs[k:k + chunk] for k in range(0, len(pairs), chunk)]
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_pwm_worker_init,
            initargs=(pwms, min_overlap),
        ) as ex:
            it = ex.map(_pwm_pairs_chunk, batches)
            if verbose:
                it = tqdm(it, total=len(batches), desc='[pwm-dist]')
            for batch_out in it:
                for i, j, d in batch_out:
                    D[i, j] = D[j, i] = d
    return D


def pwm_distance_matrix(motifs, *, format: str = "jaspar", pseudo: float = 0.01, min_overlap: int = 5,
                        workers: Optional[int] = None, verbose: bool = True):
    """Pairwise PWM distances of a motif library (1 - Sandelin-Wasserman similarity).

    ``motifs`` is anything :func:`~genomeblocks.load_motifs` takes. Returns
    ``(D, names, pwms)``: the (n x n) distances in [0, 1], the motif names, and
    the (4 x W) probability matrices in the same order."""
    lib = load_motifs(motifs, format=format)
    pwms = [lib.pfm(i, pseudo) for i in range(len(lib))]
    D = _pwm_distance_matrix(pwms, min_overlap=min_overlap, workers=workers, verbose=verbose)
    return D, list(lib.names), pwms


def cluster_motifs(D, *, cutoff: float = 0.3, linkage_method: str = 'average'):
    """Hierarchical clustering on a PWM distance matrix.

    Returns
    -------
    labels : (n,) array of 1-based cluster ids (from scipy.fcluster).
    Z : linkage matrix.
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method=linkage_method)
    labels = fcluster(Z, t=cutoff, criterion='distance')
    return labels, Z


def archetype(pwms, *, min_overlap: int = 5, weight_by_ic: bool = True):
    """Build a consensus PFM by aligning cluster members to a medoid.

    Steps:
        1. Pick the medoid (lowest mean SW distance to the rest).
        2. Align each remaining member to the medoid (best offset + strand).
        3. Average the aligned PFMs column-wise in a frame anchored on the
           medoid, expanding to fit members that extend past either edge.
           Each motif's contribution is weighted by its mean per-column
           information content if `weight_by_ic`.

    Empty columns (no member overlaps) fall back to uniform (0.25 each).
    """
    import numpy as np
    n = len(pwms)
    if n == 1:
        return pwms[0].copy()

    # Pairwise distance within the cluster (sizes are usually small)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = 1.0 - _pwm_similarity(
                pwms[i], pwms[j], min_overlap=min_overlap)
    medoid = int(D.sum(axis=1).argmin())

    p_med = pwms[medoid]
    aligned = [(0, p_med)]
    for k in range(n):
        if k == medoid:
            continue
        offset, is_rc, _ = _pwm_align(p_med, pwms[k], min_overlap=min_overlap)
        p = _pwm_rc(pwms[k]) if is_rc else pwms[k]
        aligned.append((offset, p))

    lo = min(o for o, _ in aligned)
    hi = max(o + p.shape[1] for o, p in aligned)
    Wf = hi - lo

    acc = np.zeros((4, Wf), dtype=float)
    wt = np.zeros(Wf, dtype=float)
    for o, p in aligned:
        w = float(_pfm_information_content(p).mean()) if weight_by_ic else 1.0
        a = o - lo
        acc[:, a:a + p.shape[1]] += w * p
        wt[a:a + p.shape[1]] += w

    out = np.full((4, Wf), 0.25)
    mask = wt > 0
    out[:, mask] = acc[:, mask] / wt[mask][None, :]
    out /= out.sum(axis=0, keepdims=True)
    return out


def archetype_from_names(motifs, names: List[str], *, format: str = "jaspar", match: str = "substring",
                         pseudo: float = 0.01, min_overlap: int = 5, weight_by_ic: bool = True,
                         verbose: bool = True):
    """One consensus archetype from a named subset of motifs (e.g. the
    top-enriched factors of a differential analysis), without clustering the
    whole library.

    ``names`` are matched case-insensitively as substrings of the motif names
    and descriptions (``match='substring'``) or verbatim (``'exact'``).

    Returns:
        dict with ``archetype`` (4 x W consensus), ``members`` (matched motif
        names) and ``pwms`` (their probability matrices, for a QC plot).
    """
    lib = load_motifs(motifs, format=format)
    idx = lib.indices(names, match)
    if not idx:
        raise ValueError(f"No motifs matched {names!r} (match={match!r}). Check the names exist in the "
                         f"library; try match='substring' for symbol-level lookups.")
    members = [lib.names[i] for i in idx]
    pwms = [lib.pfm(i, pseudo) for i in idx]
    if verbose:
        print(f"[archetype_from_names] {len(members)} motifs → 1 archetype: {', '.join(members[:5])}"
              f"{f' (+{len(members) - 5} more)' if len(members) > 5 else ''}")
    return {"archetype": archetype(pwms, min_overlap=min_overlap, weight_by_ic=weight_by_ic),
            "members": members, "pwms": pwms}


def build_archetypes(
    motifs,
    *,
    format: str = 'jaspar',
    cutoff: float = 0.3,
    min_overlap: int = 5,
    linkage_method: str = 'average',
    pseudo: float = 0.01,
    weight_by_ic: bool = True,
    workers: Optional[int] = None,
    name_prefix: str = 'ARCH',
    verbose: bool = True,
):
    """End-to-end: load → pairwise distance → cluster → archetype per cluster.

    Parameters
    ----------
    motifs, format
        A motif file and its format, or any motif source load_motifs takes.
    cutoff : float
        Distance cutoff for hierarchical clustering (``fcluster`` criterion
        ``'distance'``). Equivalent to "merge motifs with SW similarity ≥
        (1 − cutoff)". Vierstra used 0.7 on a related metric; for SW
        distance, 0.25–0.35 is a typical range — use ``plot_dendrogram`` to
        tune empirically for your library.
    min_overlap, linkage_method, pseudo, weight_by_ic, workers, name_prefix
        See module docstrings on individual functions.

    Returns
    -------
    dict with keys: ``D``, ``Z``, ``labels``, ``names``, ``pwms``,
    ``archetypes`` (cluster name → consensus PFM), ``members`` (cluster
    name → list of motif names).
    """
    import numpy as np
    D, names, pwms = pwm_distance_matrix(
        motifs, format=format, pseudo=pseudo,
        min_overlap=min_overlap, workers=workers, verbose=verbose)
    labels, Z = cluster_motifs(D, cutoff=cutoff, linkage_method=linkage_method)

    n_clusters = int(labels.max())
    width = max(2, len(str(n_clusters)))
    archetypes, members = {}, {}
    for cid in range(1, n_clusters + 1):
        idx = np.where(labels == cid)[0]
        cluster_pwms = [pwms[i] for i in idx]
        arch_pfm = archetype(cluster_pwms, min_overlap=min_overlap,
                             weight_by_ic=weight_by_ic)
        key = f"{name_prefix}_{cid:0{width}d}"
        archetypes[key] = arch_pfm
        members[key] = [names[i] for i in idx]

    if verbose:
        sizes = sorted([len(v) for v in members.values()], reverse=True)
        print(f"[archetypes] {len(names)} motifs → {n_clusters} archetypes "
              f"(cutoff={cutoff}, linkage={linkage_method}). "
              f"Top sizes: {sizes[:5]}")

    return {
        'D': D, 'Z': Z, 'labels': labels, 'names': names, 'pwms': pwms,
        'archetypes': archetypes, 'members': members,
    }
