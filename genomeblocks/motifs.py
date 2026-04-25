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
        lfc = np.log2((1 + m) / (1 + mean_ref))
        out[f'mean_{name}'] = m
        out[f'LFC_{name}'] = lfc
        lfc_by_group[name] = lfc

    if len(lfc_by_group) == 2:
        a, b = list(lfc_by_group.keys())
        out['LFC'] = lfc_by_group[a] - lfc_by_group[b]

    return pd.DataFrame(out)


Loci.scan_motifs = scan_motifs
Loci.scan_motifs_matrix = scan_motifs_matrix
