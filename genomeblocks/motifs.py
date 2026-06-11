"""Motif scanning utilities."""
from typing import Dict, Iterable, List, Optional
import os

from .loci import Loci


def _parse_fasta(path):
    import gzip
    name = ""
    seq = list()
    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path, "r")
    for _, l in enumerate(f):
        l = l.rstrip("\n")
        if l[0] == ">":
            if len(seq) == 0:
                name = l[1:].split(" ")[0]
                continue
            yield (name, "".join(seq))
            seq = list()
            name = l[1:].split(" ")[0]
        else:
            seq.append(l)
    yield (name, "".join(seq))


def make_genome(path):
    return {k:v for k,v in _parse_fasta(path)}


def scan_motifs(s, genome, motif_path, motif_format='jaspar', r=250, threshold=13.0, norm=True, verbose=True):
    import lightmotif
    from tqdm import tqdm

    # ensure genome is a dict
    if not isinstance(genome, dict): genome = make_genome(genome)
    n_motif = len([_ for _ in lightmotif.load(motif_path, format=motif_format)])
    M = {}
    # extract every sequence once; stripe once; reuse across motifs
    striped = []
    for l in s:
        seq = l.sequence(genome, r=r).upper()
        if len(seq) != (2 * r):
            continue
        try:
            striped.append(lightmotif.stripe(seq))
        except ValueError:
            # skip sequences with invalid characters (e.g. N at boundaries)
            continue
    iterator = lightmotif.load(motif_path, format=motif_format)
    if verbose: iterator = tqdm(iterator, total=n_motif)
    for m in iterator:
        pssm = m.counts.normalize(0.1).log_odds()
        hits = 0
        for sseq in striped:
            for _ in lightmotif.scan(pssm, sseq, threshold=threshold):
                hits += 1
        M[m.name] = hits / len(m.counts) if norm else hits
    return M


def _extract_sequences(loci, genome, r: int):
    """Extract (uid, seq) pairs at fixed length 2*r. Skips short/invalid."""
    out = []
    for l in loci:
        seq = l.sequence(genome, r=r).upper()
        if len(seq) != (2 * r):
            continue
        out.append((l.uid, seq))
    return out


_WORKER_STATE: dict = {}


def _worker_init(seqs, motif_path, motif_format, threshold, norm):
    """Pool initializer: cache sequences (pre-striped) and motif list per worker."""
    import lightmotif
    striped = []
    for seq in seqs:
        try:
            striped.append(lightmotif.stripe(seq))
        except ValueError:
            striped.append(None)
    motifs_list = list(lightmotif.load(motif_path, format=motif_format))
    _WORKER_STATE['striped'] = striped
    _WORKER_STATE['motifs'] = motifs_list
    _WORKER_STATE['threshold'] = threshold
    _WORKER_STATE['norm'] = norm


def _scan_motif_indices(indices):
    """Worker task: scan a batch of motif indices against the cached sequences."""
    import lightmotif
    striped = _WORKER_STATE['striped']
    motifs_list = _WORKER_STATE['motifs']
    threshold = _WORKER_STATE['threshold']
    norm = _WORKER_STATE['norm']
    out = []
    for idx in indices:
        motif = motifs_list[idx]
        pssm = motif.counts.normalize(0.1).log_odds()
        width = len(motif.counts)
        counts = []
        for sseq in striped:
            if sseq is None:
                counts.append(0)
                continue
            n = 0
            for _ in lightmotif.scan(pssm, sseq, threshold=threshold):
                n += 1
            counts.append(n / width if norm else n)
        out.append((motif.name, counts))
    return out


def scan_motifs_matrix(
    s,
    genome,
    motif_path,
    motif_format: str = 'jaspar',
    r: int = 250,
    threshold: float = 13.0,
    norm: bool = True,
    workers: Optional[int] = None,
    verbose: bool = True,
):
    """Precompute a (n_loci x n_motifs) count matrix.

    Each locus's center-window sequence (length 2*r) is extracted once, then
    every motif is scanned against every sequence. Motif work is distributed
    across a process pool — motifs are independent.

    Returns
    -------
    pandas.DataFrame : rows are locus uids (order follows input ``s``,
        skipping any that failed the 2*r length check), columns are motif
        names. Cells hold hit counts (normalized by motif width if ``norm``).
    """
    import pandas as pd
    import lightmotif
    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm

    if not isinstance(genome, dict): genome = make_genome(genome)

    pairs = _extract_sequences(s, genome, r)
    if not pairs:
        return pd.DataFrame()
    uids = [u for u, _ in pairs]
    seqs = [q for _, q in pairs]

    n_motif = sum(1 for _ in lightmotif.load(motif_path, format=motif_format))

    if workers is None:
        workers = max(1, (os.cpu_count() or 2) // 2)
    # parallelism helps only when the compute dominates the ~1–2s process
    # startup + import cost; for tiny jobs, stay serial.
    if workers > 1 and (n_motif < 16 or len(seqs) < 200):
        workers = 1

    results: Dict[str, List[float]] = {}
    if workers == 1:
        _worker_init(seqs, motif_path, motif_format, threshold, norm)
        idx_iter = range(n_motif)
        if verbose: idx_iter = tqdm(idx_iter, total=n_motif, desc="[motifs]")
        for i in idx_iter:
            for name, counts in _scan_motif_indices([i]):
                results[name] = counts
    else:
        # distribute motifs in contiguous chunks so each worker scans many
        # motifs against its cached sequence list (no per-task sequence IPC).
        chunk = max(1, n_motif // (workers * 4))
        batches = [list(range(i, min(i + chunk, n_motif))) for i in range(0, n_motif, chunk)]
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_worker_init,
            initargs=(seqs, motif_path, motif_format, threshold, norm),
        ) as ex:
            it = ex.map(_scan_motif_indices, batches)
            if verbose:
                it = tqdm(it, total=len(batches), desc="[motifs]")
            for batch_out in it:
                for name, counts in batch_out:
                    results[name] = counts

    return pd.DataFrame(results, index=uids)


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


def _hit_position(hit):
    """Best-effort extraction of a 0-based start position from a lightmotif Hit.
    The Python binding's attribute name has varied across versions — try the
    common spellings before falling back to iteration."""
    for attr in ('position', 'pos', 'start'):
        v = getattr(hit, attr, None)
        if v is not None:
            return int(v)
    try:
        return int(next(iter(hit)))
    except (TypeError, ValueError, StopIteration):
        return None


def _resolve_anchors(motifs_list, anchors):
    """Resolve user-given names to motif indices by case-insensitive substring."""
    anchors_lc = [a.lower() for a in anchors]
    return [i for i, m in enumerate(motifs_list)
            if any(a in m.name.lower() for a in anchors_lc)]


def _mask_anchor_hits(seqs, motifs_list, anchor_idx, *, window, threshold, seed):
    """Replace ±window bp around each anchor-motif hit center with random
    A/C/G/T. Random fill is used (rather than 'N') because lightmotif.stripe
    rejects N — the few stray hits a random ~20-bp filler can create are
    background noise that's equal between A and B groups in the downstream
    differential test."""
    import lightmotif
    import numpy as np

    anchor_pssms = [(motifs_list[i].counts.normalize(0.1).log_odds(),
                     len(motifs_list[i].counts)) for i in anchor_idx]
    rng = np.random.default_rng(seed)
    bases = np.frombuffer(b'ACGT', dtype=np.uint8)
    out = []
    for seq in seqs:
        try:
            sseq = lightmotif.stripe(seq)
        except ValueError:
            out.append(seq)              # unstripable → scanner will skip it too
            continue
        arr = bytearray(seq.encode('ascii'))
        L = len(arr)
        for pssm, w in anchor_pssms:
            for hit in lightmotif.scan(pssm, sseq, threshold=threshold):
                pos = _hit_position(hit)
                if pos is None:
                    continue
                c = pos + w // 2
                lo, hi = max(0, c - window), min(L, c + window + 1)
                if hi <= lo:
                    continue
                arr[lo:hi] = bases[rng.integers(0, 4, size=hi - lo)].tobytes()
        out.append(arr.decode('ascii'))
    return out


def scan_motifs_matrix_masked(
    s,
    genome,
    motif_path,
    anchors: List[str],
    *,
    motif_format: str = 'jaspar',
    r: int = 250,
    window: int = 10,
    threshold: float = 13.0,
    anchor_threshold: Optional[float] = None,
    norm: bool = True,
    skip_anchors: bool = True,
    seed: Optional[int] = None,
    workers: Optional[int] = None,
    verbose: bool = True,
):
    """Scan motifs after masking out matches of anchor motifs in each sequence.

    For each locus's center window (length 2*r), find every match to any
    motif whose name matches one of ``anchors`` (case-insensitive substring
    against the library), replace a ±``window``-bp window around each
    match center with random A/C/G/T, then run the standard scanner on
    the masked sequences. Anchor motifs are excluded from the output by
    default (``skip_anchors=True``).

    Use this to ask "which motifs enrich at these regions *independent of
    the anchor*?" — e.g., mask CTCF and ask which co-factors differentiate
    your two CRE sets.

    Parameters
    ----------
    s, genome, motif_path, motif_format, r, threshold, norm, workers, verbose:
        Same semantics as :func:`scan_motifs_matrix`.
    anchors : list[str]
        Motif name fragments to mask (case-insensitive substring). E.g.
        ``['CTCF']`` matches every CTCF.* PSSM in the library.
    window : int
        Half-width of the mask on each side of the motif-match center.
        Default 10 bp (mask spans 21 bp).
    anchor_threshold : float, optional
        Threshold used when finding anchor hits to mask. Defaults to
        ``threshold`` (same as scanning).
    skip_anchors : bool
        If True (default), the anchor motifs themselves are excluded from
        the returned matrix (they'd be ~zero after masking anyway).
    seed : int, optional
        RNG seed for the random fill — set for reproducibility.

    Returns
    -------
    pandas.DataFrame
        (n_loci × n_motifs), same shape as :func:`scan_motifs_matrix`
        (minus anchor columns when ``skip_anchors``).
    """
    import pandas as pd
    import lightmotif
    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm

    if anchor_threshold is None:
        anchor_threshold = threshold
    if not isinstance(genome, dict):
        genome = make_genome(genome)

    pairs = _extract_sequences(s, genome, r)
    if not pairs:
        return pd.DataFrame()
    uids = [u for u, _ in pairs]
    seqs = [q for _, q in pairs]

    motifs_list = list(lightmotif.load(motif_path, format=motif_format))
    anchor_idx = _resolve_anchors(motifs_list, anchors)
    if not anchor_idx:
        raise ValueError(f"No motifs matched anchors {anchors!r} in {motif_path}")
    if verbose:
        names = [motifs_list[i].name for i in anchor_idx]
        head = ', '.join(names[:5])
        more = f" (+{len(names) - 5} more)" if len(names) > 5 else ''
        print(f"[mask] {len(anchor_idx)} anchor PSSM(s): {head}{more}")

    if verbose:
        seqs_iter = tqdm(seqs, total=len(seqs), desc='[mask]')
        # materialize the generator so _mask_anchor_hits sees a list
        masked_input = list(seqs_iter)
    else:
        masked_input = seqs
    masked = _mask_anchor_hits(
        masked_input, motifs_list, anchor_idx,
        window=window, threshold=anchor_threshold, seed=seed,
    )

    # Choose motifs to scan in the output
    if skip_anchors:
        anchor_set = set(anchor_idx)
        scan_idx = [i for i in range(len(motifs_list)) if i not in anchor_set]
    else:
        scan_idx = list(range(len(motifs_list)))
    n_motif_out = len(scan_idx)

    if workers is None:
        workers = max(1, (os.cpu_count() or 2) // 2)
    if workers > 1 and (n_motif_out < 16 or len(masked) < 200):
        workers = 1

    results: Dict[str, List[float]] = {}
    if workers == 1:
        _worker_init(masked, motif_path, motif_format, threshold, norm)
        idx_iter = scan_idx
        if verbose:
            idx_iter = tqdm(idx_iter, total=n_motif_out, desc='[scan]')
        for i in idx_iter:
            for name, counts in _scan_motif_indices([i]):
                results[name] = counts
    else:
        chunk = max(1, n_motif_out // (workers * 4))
        batches = [scan_idx[i:i + chunk] for i in range(0, n_motif_out, chunk)]
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_worker_init,
            initargs=(masked, motif_path, motif_format, threshold, norm),
        ) as ex:
            it = ex.map(_scan_motif_indices, batches)
            if verbose:
                it = tqdm(it, total=len(batches), desc='[scan]')
            for batch_out in it:
                for name, counts in batch_out:
                    results[name] = counts

    return pd.DataFrame(results, index=uids)


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

def _pfm_from_motif(motif, pseudo=0.01):
    """Convert a lightmotif Motif to a (4, W) PFM with rows in ACGT order.

    Adds `pseudo` to each count before column-normalizing to probabilities.
    Assumes lightmotif's CountMatrix supports base-keyed access
    (``cm[base]`` → iterable of per-position counts).
    """
    import numpy as np
    # lightmotif's CountMatrix exposes the buffer protocol — np.asarray
    # yields shape (W, 5) with column order [A, C, T, G, N]. Drop N, then
    # permute [0,1,3,2] to get ACGT, add pseudocount (avoids 0/0 → NaN on
    # any all-zero column), normalize, and transpose to (4, W).
    arr = np.asarray(motif.counts, dtype=float)
    counts = arr[:, [0, 1, 3, 2]] + pseudo
    return (counts / counts.sum(axis=1, keepdims=True)).T


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
    import numpy as np
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


def pwm_distance_matrix(
    motif_path,
    *,
    motif_format: str = 'jaspar',
    pseudo: float = 0.01,
    min_overlap: int = 5,
    workers: Optional[int] = None,
    verbose: bool = True,
):
    """Load motifs and compute pairwise PWM distance matrix.

    Returns
    -------
    D : (n, n) numpy array of distances in [0, 1] (1 − Sandelin-Wasserman).
    names : list of motif names in row/column order.
    pwms : list of (4, W) PFMs in the same order.
    """
    import lightmotif
    motifs_list = list(lightmotif.load(motif_path, format=motif_format))
    names = [m.name for m in motifs_list]
    pwms = [_pfm_from_motif(m, pseudo=pseudo) for m in motifs_list]
    D = _pwm_distance_matrix(pwms, min_overlap=min_overlap,
                             workers=workers, verbose=verbose)
    return D, names, pwms


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


def archetype_from_names(
    motif_path,
    names: List[str],
    *,
    motif_format: str = 'jaspar',
    match: str = 'substring',
    pseudo: float = 0.01,
    min_overlap: int = 5,
    weight_by_ic: bool = True,
    verbose: bool = True,
):
    """Build a single consensus archetype from a named subset of motifs.

    Use when you already have a coherent group of motifs in mind — e.g. the
    top-enriched factors from a differential analysis — and just want one
    representative PWM, without clustering the whole library.

    Parameters
    ----------
    motif_path, motif_format : str
        Motif file readable by lightmotif.
    names : list[str]
        Motif names (or fragments) to include. With ``match='substring'``
        (default), 'CTCF' matches every 'CTCF.H14CORE.*' PSSM in the
        library. With ``match='exact'``, names must equal the library's
        motif names verbatim — appropriate when you pass the full
        ``enr['Factor']`` values.
    match : {'substring', 'exact'}
        How to resolve ``names``. Substring matching is case-insensitive.
    pseudo, min_overlap, weight_by_ic
        See :func:`archetype`.

    Returns
    -------
    dict with:
        ``archetype``: consensus PFM, shape (4, W)
        ``members``  : list of full motif names that were included
        ``pwms``     : list of per-motif PFMs in member order (handy for
                       a ``plot_cluster_members`` QC plot)
    """
    import lightmotif

    motifs_list = list(lightmotif.load(motif_path, format=motif_format))

    if match == 'exact':
        wanted = set(names)
        idx = [i for i, m in enumerate(motifs_list) if m.name in wanted]
    elif match == 'substring':
        names_lc = [n.lower() for n in names]
        idx = [i for i, m in enumerate(motifs_list)
               if any(n in m.name.lower() for n in names_lc)]
    else:
        raise ValueError(f"match must be 'exact' or 'substring', got {match!r}")

    if not idx:
        raise ValueError(
            f"No motifs matched {names!r} in {motif_path} "
            f"(match={match!r}). Check the names exist in the library; "
            f"try match='substring' for symbol-level lookups.")

    matched_names = [motifs_list[i].name for i in idx]
    pwms = [_pfm_from_motif(motifs_list[i], pseudo=pseudo) for i in idx]

    if verbose:
        head = ', '.join(matched_names[:5])
        more = f' (+{len(matched_names) - 5} more)' if len(matched_names) > 5 else ''
        print(f"[archetype_from_names] {len(matched_names)} motifs → 1 archetype: "
              f"{head}{more}")

    arch_pfm = archetype(pwms, min_overlap=min_overlap, weight_by_ic=weight_by_ic)
    return {'archetype': arch_pfm, 'members': matched_names, 'pwms': pwms}


def build_archetypes(
    motif_path,
    *,
    motif_format: str = 'jaspar',
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
    motif_path, motif_format : str
        Motif file and format readable by lightmotif (e.g. 'jaspar', 'meme').
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
        motif_path, motif_format=motif_format, pseudo=pseudo,
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


def write_meme(archetypes, path, *, alphabet: str = 'ACGT',
               bg=(0.25, 0.25, 0.25, 0.25)):
    """Write an archetype dict to MEME-format so the archetypes can be fed
    back into ``scan_motifs_matrix(motif_path, motif_format='meme', ...)``.
    """
    bg_str = ' '.join(f'{a} {p:.4f}' for a, p in zip(alphabet, bg))
    with open(path, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write(f"ALPHABET= {alphabet}\n\n")
        f.write("strands: + -\n\n")
        f.write(f"Background letter frequencies\n{bg_str}\n\n")
        for name, pfm in archetypes.items():
            W = pfm.shape[1]
            f.write(f"MOTIF {name}\n")
            f.write(f"letter-probability matrix: alength= 4 w= {W}\n")
            for j in range(W):
                f.write(" ".join(f"{pfm[i, j]:.6f}" for i in range(4)) + "\n")
            f.write("\n")


# ── Plotting (logomaker — lazy imported) ──────────────────────────────────

def plot_archetype(pfm, ax=None, *, title=None, alphabet: str = 'ACGT',
                   show_xticks: bool = True, ylim=(0, 2)):
    """Plot a (4, W) PFM as an information-content sequence logo (bits)."""
    import logomaker
    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.DataFrame(pfm.T, columns=list(alphabet))
    ic_df = logomaker.transform_matrix(
        df, from_type='probability', to_type='information')

    if ax is None:
        _, ax = plt.subplots(figsize=(max(2.5, ic_df.shape[0] * 0.35), 1.4))

    logo = logomaker.Logo(ic_df, ax=ax, show_spines=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    ax.set_ylabel('bits', fontsize=8)
    ax.set_ylim(*ylim)
    if show_xticks:
        ax.set_xticks(range(ic_df.shape[0]))
        ax.set_xticklabels(range(1, ic_df.shape[0] + 1), fontsize=7)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_archetypes(archetypes, *, ncols: int = 4, figsize_per=(2.8, 1.3),
                    members=None, sort_by_size: bool = True):
    """Grid of archetype logos. Pass ``members=`` to annotate cluster size."""
    import matplotlib.pyplot as plt
    import math

    items = list(archetypes.items())
    if sort_by_size and members is not None:
        items.sort(key=lambda kv: -len(members.get(kv[0], [])))
    n = len(items)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per[0] * ncols, figsize_per[1] * nrows),
        squeeze=False,
    )
    flat = axes.flatten()
    for ax, (name, pfm) in zip(flat, items):
        title = name
        if members is not None and name in members:
            title = f"{name}  (n={len(members[name])})"
        plot_archetype(pfm, ax=ax, title=title)
    for ax in flat[n:]:
        ax.axis('off')
    plt.tight_layout()
    return fig


def plot_cluster_members(archetype_pfm, member_pwms, member_names, *,
                         ncols: int = 3, figsize_per=(2.8, 1.3),
                         archetype_title: str = 'ARCHETYPE'):
    """Archetype on top + each member underneath for QC."""
    import matplotlib.pyplot as plt
    import math

    n = len(member_pwms) + 1
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per[0] * ncols, figsize_per[1] * nrows),
        squeeze=False,
    )
    flat = axes.flatten()
    plot_archetype(archetype_pfm, ax=flat[0], title=archetype_title)
    for ax, pfm, name in zip(flat[1:], member_pwms, member_names):
        plot_archetype(pfm, ax=ax, title=name)
    for ax in flat[n:]:
        ax.axis('off')
    plt.tight_layout()
    return fig


def plot_dendrogram(Z, names=None, *, cutoff=None, ax=None,
                    color_threshold=None, leaf_font_size: int = 6,
                    no_labels: bool = False):
    """Dendrogram from a linkage matrix Z (output of ``cluster_motifs``).
    Draws a horizontal line at ``cutoff`` if provided."""
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    if ax is None:
        _, ax = plt.subplots(figsize=(14, 4))
    color_threshold = color_threshold if color_threshold is not None else cutoff
    dendrogram(
        Z, labels=names, ax=ax,
        color_threshold=color_threshold,
        leaf_font_size=leaf_font_size,
        leaf_rotation=90,
        no_labels=no_labels,
    )
    if cutoff is not None:
        ax.axhline(cutoff, color='red', ls='--', lw=0.8, alpha=0.6,
                   label=f'cutoff={cutoff}')
        ax.legend(loc='upper right', fontsize=7)
    ax.set_ylabel('SW distance', fontsize=9)
    return ax


Loci.scan_motifs = scan_motifs
Loci.scan_motifs_matrix = scan_motifs_matrix
Loci.scan_motifs_matrix_masked = scan_motifs_matrix_masked
