#!/usr/bin/env python3
"""
Benchmark genomeblocks.signal.signal() under different parallelism strategies.

Default workload: 50 bigwigs × 50,000 loci, 200 bins, mean aggregation.
Bigwig paths are cycled from a pool of real on-disk bigwigs; duplicates are
fine — each process opens its own handle so this models 50 distinct tracks
from a filesystem perspective.

Usage:
    python bench_signal_parallel.py
    python bench_signal_parallel.py --n-bw 50 --n-loci 50000 \
        --workers 1 2 4 8 16 24 32 --skip-thread

What is compared, for each worker count W:
  1. Sequential (only W=1)
  2. Multiprocessing (signal's built-in ProcessPoolExecutor)
  3. Threading (ThreadPoolExecutor sharing one cube, opening one handle per
     thread) — expected to be slower than sequential because pybigtools'
     Python binding serializes across threads
"""
from __future__ import annotations
import argparse
import os
import random
import time
from pathlib import Path

import numpy as np
from tqdm import tqdm


DEFAULT_BW_POOL = [
    "/home/ualtintas/bluegill/data/bw/hg19.100way.phastCons.chr1.10kb.bw",
    "/home/ualtintas/bluegill/data/bw/A549.GR.0h-DEX.hg38.chr1.10kb.bw",
    "/home/ualtintas/bluegill/data/bw/A549.GR.4h-DEX.hg38.chr1.10kb.bw",
]


# ── helpers ──────────────────────────────────────────────────────────────────

def _fmt_time(sec: float) -> str:
    if sec < 1:
        return f"{sec*1000:6.1f} ms"
    if sec < 60:
        return f"{sec:6.2f} s "
    return f"{sec/60:5.2f} min"


def _fmt_rate(n: int, sec: float) -> str:
    return f"{n/sec:>10,.0f} rt/s"


def _build_bw_list(pool, n: int):
    """Cycle through pool to reach n paths."""
    existing = [p for p in pool if Path(p).exists()]
    if not existing:
        raise FileNotFoundError(
            f"None of the pool bigwigs exist: {pool}. "
            "Pass --bw-pool with available files."
        )
    return [existing[i % len(existing)] for i in range(n)]


def _build_loci(bw_path: str, n: int, flank: int, seed: int = 42):
    """Generate n random Locus objects fitting inside the bigwig's chroms."""
    from genomeblocks.signal import _bw_open
    from genomeblocks.locus import Locus
    from genomeblocks.loci import Loci

    h = _bw_open(bw_path)
    chroms = h.chroms()
    h.close()

    pool = [(c, sz) for c, sz in chroms.items() if sz > 2 * flank + 1000]
    if not pool:
        raise RuntimeError(f"No chrom large enough for flank={flank}")

    rng = random.Random(seed)
    loci = Loci()
    for _ in range(n):
        c, sz = rng.choice(pool)
        center = rng.randint(flank + 1, sz - flank - 1)
        loci.append(Locus(chrom=c, start=center - 100, end=center + 100))
    return loci


# ── thread-based extraction (expected to be slow) ────────────────────────────

def _thread_extract(loci, bigwigs, *, n_bins, flank, agg, exact, workers):
    """Run signal extraction with a ThreadPoolExecutor (shared cube).

    Mirrors the structure of signal._run_multiprocess but threads instead of
    processes. Each thread opens its own handles and writes its slice directly
    into a shared ndarray — no shared-memory dance needed for threads.
    """
    from concurrent.futures import ThreadPoolExecutor
    from genomeblocks.signal import (
        _bw_open, _extract_chunk, _even_ranges,
    )

    n_loci = len(loci)
    n_tracks = len(bigwigs)
    chrom_list = [l.chrom for l in loci]
    starts = np.fromiter((l.start for l in loci), dtype=np.int64, count=n_loci)
    ends = np.fromiter((l.end for l in loci), dtype=np.int64, count=n_loci)
    cube = np.zeros((n_loci, n_tracks, n_bins), dtype=np.float32)

    t_splits = min(workers, n_tracks)
    l_per_t = max(1, workers // t_splits)
    t_ranges = _even_ranges(n_tracks, t_splits)
    l_ranges = _even_ranges(n_loci, l_per_t)

    def run(t_lo, t_hi, l_lo, l_hi):
        hs = [_bw_open(p) for p in bigwigs[t_lo:t_hi]]
        try:
            ch = hs[0].chroms()
            _extract_chunk(cube, hs, ch, chrom_list, starts, ends,
                           t_lo, l_lo, l_hi, n_bins, flank, agg, False, exact)
        finally:
            for h in hs:
                h.close()

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run, t_lo, t_hi, l_lo, l_hi)
                for t_lo, t_hi in t_ranges
                for l_lo, l_hi in l_ranges]
        for f in futs:
            f.result()
    return cube


# ── backend comparison (sequential, small subset) ────────────────────────────

def _pybigwig_direct(bw_paths, loci, n_bins, flank, agg):
    """Mimic signal()'s extraction loop using pyBigWig as backend.

    Runs fully sequentially (one file at a time, one region at a time).
    Returns (cube, seconds). No signal() involvement — this is what the
    baseline Python ecosystem gives you without genomeblocks.
    """
    import pyBigWig
    n_loci = len(loci)
    n_tracks = len(bw_paths)
    cube = np.zeros((n_loci, n_tracks, n_bins), dtype=np.float32)
    t0 = time.perf_counter()
    for t, p in enumerate(bw_paths):
        f = pyBigWig.open(p)
        try:
            chroms = f.chroms()
            for i, l in enumerate(loci):
                if l.chrom not in chroms:
                    continue
                c = (l.start + l.end) // 2
                sz = chroms[l.chrom]
                L, R = max(0, c - flank), min(sz, c + flank)
                if R <= L:
                    continue
                xs = f.stats(l.chrom, L, R, type=agg, nBins=n_bins)
                for j, x in enumerate(xs):
                    if x is not None:
                        cube[i, t, j] = x
        finally:
            f.close()
    return cube, time.perf_counter() - t0


def _run_backend_compare(loci, bigwigs, *, n_bins, flank, agg,
                         n_loci_cap, n_bw_cap):
    """Compare pybigtools / pure-python / pyBigWig on a small subset (seq, 1 worker)."""
    from genomeblocks.signal import signal
    from genomeblocks.loci import Loci

    n_loci_cap = min(n_loci_cap, len(loci))
    n_bw_cap = min(n_bw_cap, len(bigwigs))
    sub_loci = Loci(loci[:n_loci_cap])
    sub_bw = bigwigs[:n_bw_cap]
    n_rt = n_loci_cap * n_bw_cap

    out = []  # (name, time, throughput)

    # pybigtools via signal()
    t0 = time.perf_counter()
    signal(sub_loci, sub_bw, n_bins=n_bins, flank=flank, agg=agg,
           workers=1, progress=False, verbose=False,
           backend='pybigtools', exact=True)
    out.append(("pybigtools (signal)", time.perf_counter() - t0, n_rt))

    # pure-python via signal()
    try:
        t0 = time.perf_counter()
        signal(sub_loci, sub_bw, n_bins=n_bins, flank=flank, agg=agg,
               workers=1, progress=False, verbose=False,
               backend='bigwig')
        out.append(("pure-python (signal)", time.perf_counter() - t0, n_rt))
    except Exception as exc:
        out.append((f"pure-python (err: {exc!r})", float('inf'), n_rt))

    # pyBigWig direct (not wired into signal(); call its loop ourselves)
    try:
        _, t = _pybigwig_direct(sub_bw, sub_loci, n_bins, flank, agg)
        out.append(("pyBigWig (direct)", t, n_rt))
    except ImportError:
        out.append(("pyBigWig (not installed)", float('inf'), n_rt))

    return out, n_loci_cap, n_bw_cap


# ── rendering smoke tests ────────────────────────────────────────────────────

def _run_render_smoke(loci, bigwigs, cube, *, n_bins, flank, out_dir):
    """Render a browser view, heatmap, and profile from the extracted cube.

    Small inputs (8 loci / 4 tracks), timed independently. Each wrapped in
    try/except so a single rendering failure does not kill the bench run.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from genomeblocks.tags import Tags

    results = []

    # 1) Browser: 1 region × 2 bigwig tracks ------------------------------
    try:
        from genomeblocks.browser import browser
        locus = loci[0]
        region = (locus.chrom, locus.start - flank, locus.end + flank)
        tracks = {f"bw_{i}": bigwigs[i] for i in range(min(2, len(bigwigs)))}
        t0 = time.perf_counter()
        fig, _ = browser(region, tracks, bw_n_bins=500)
        path = out_dir / "bench_browser.png"
        fig.savefig(path, dpi=100, bbox_inches="tight")
        plt.close(fig)
        results.append(("browser (1 loc × 2bw)", True,
                        time.perf_counter() - t0, path, None))
    except Exception as exc:
        results.append(("browser (1 loc × 2bw)", False, 0.0, None, repr(exc)))

    # 2) plot_heatmap: small slice of the already-extracted cube -----------
    try:
        from genomeblocks.signal import plot_heatmap
        n_show = min(64, cube.shape[0])
        n_tr = min(4, cube.shape[1])
        sub = cube[:n_show, :n_tr, :]
        from genomeblocks.loci import Loci
        sub_loci = Loci(loci[:n_show])
        tags = Tags.make(sub_loci, verbose=False)
        half = n_show // 2
        tags.add({
            "first_half": Loci(loci[:half]),
            "second_half": Loci(loci[half:n_show]),
        })
        t0 = time.perf_counter()
        fig = plot_heatmap(loci[:n_show], sub, tags=tags,
                           sets=["first_half", "second_half"],
                           samples=[f"bw_{i}" for i in range(n_tr)],
                           vmax=max(1.0, float(sub.max())),
                           ymax=max(1.0, float(sub.mean(axis=(0, 2)).max() * 1.5)))
        path = out_dir / "bench_heatmap.png"
        fig.savefig(path, dpi=80, bbox_inches="tight")
        plt.close(fig)
        results.append((f"heatmap ({n_show}×{n_tr})", True,
                        time.perf_counter() - t0, path, None))
    except Exception as exc:
        results.append(("heatmap", False, 0.0, None, repr(exc)))

    # 3) plot_profiles: average profile per group --------------------------
    try:
        from genomeblocks.signal import plot_profiles
        n_show = min(128, cube.shape[0])
        sub = cube[:n_show, :1, :]
        from genomeblocks.loci import Loci
        sub_loci = Loci(loci[:n_show])
        tags = Tags.make(sub_loci, verbose=False)
        half = n_show // 2
        tags.add({
            "first_half": Loci(loci[:half]),
            "second_half": Loci(loci[half:n_show]),
        })
        t0 = time.perf_counter()
        fig = plot_profiles(loci[:n_show], sub, tags=tags,
                            sets=["first_half", "second_half"])
        path = out_dir / "bench_profiles.png"
        fig.savefig(path, dpi=80, bbox_inches="tight")
        plt.close(fig)
        results.append((f"profiles ({n_show})", True,
                        time.perf_counter() - t0, path, None))
    except Exception as exc:
        results.append(("profiles", False, 0.0, None, repr(exc)))

    return results


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-bw", type=int, default=50,
                    help="Number of bigwig tracks (pool is cycled). Default 50.")
    ap.add_argument("--n-loci", type=int, default=50_000,
                    help="Number of random loci. Default 50,000.")
    ap.add_argument("--n-bins", type=int, default=200)
    ap.add_argument("--flank", type=int, default=3000)
    ap.add_argument("--agg", default="mean")
    ap.add_argument("--exact", action="store_true", default=True,
                    help="Base-accurate binning (default).")
    ap.add_argument("--no-exact", action="store_false", dest="exact",
                    help="Use pybigtools zoom interpolation (approximate, faster).")
    ap.add_argument("--workers", type=int, nargs='+',
                    default=[1, 2, 4, 8, 16, 24, 32],
                    help="Worker counts to sweep.")
    ap.add_argument("--skip-thread", action="store_true",
                    help="Skip thread-pool benchmark (it's slow by design).")
    ap.add_argument("--thread-workers", type=int, nargs='+',
                    default=[2, 4, 8],
                    help="Worker counts to try for threading. Default keeps it small.")
    ap.add_argument("--bw-pool", nargs='+', default=DEFAULT_BW_POOL,
                    help="Pool of real bigwig paths to cycle from.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-render", action="store_true",
                    help="Skip browser/heatmap/profile rendering smoke tests.")
    ap.add_argument("--render-dir", default="/tmp",
                    help="Directory to write smoke-test figures to.")
    ap.add_argument("--skip-backend-compare", action="store_true",
                    help="Skip pybigtools vs pure-python vs pyBigWig subset comparison.")
    ap.add_argument("--cmp-n-loci", type=int, default=2000,
                    help="Loci cap for backend comparison subset (pyBigWig is slow). Default 2000.")
    ap.add_argument("--cmp-n-bw", type=int, default=10,
                    help="Bigwig cap for backend comparison subset. Default 10.")
    args = ap.parse_args()

    print("=" * 78)
    print(f"SIGNAL PARALLELISM BENCHMARK")
    print(f"  bigwigs: {args.n_bw}  loci: {args.n_loci:,}  "
          f"bins: {args.n_bins}  flank: {args.flank}bp  agg: {args.agg}  "
          f"exact: {args.exact}")
    print(f"  CPU cores: {os.cpu_count()}  (signal caps workers at cpu_count()//2)")
    print("=" * 78)

    # build workload once ----------------------------------------------------
    bigwigs = _build_bw_list(args.bw_pool, args.n_bw)
    uniq = len({Path(p).name for p in bigwigs})
    print(f"\nBigwig pool: {uniq} unique file(s), cycled to {len(bigwigs)} paths")

    from genomeblocks.signal import signal, _bw_backend
    print(f"Backend: {_bw_backend}")

    print(f"\nGenerating {args.n_loci:,} random loci ...")
    t0 = time.perf_counter()
    loci = _build_loci(bigwigs[0], args.n_loci, args.flank, seed=args.seed)
    print(f"  loci built in {time.perf_counter() - t0:.2f}s")

    total_region_tracks = args.n_loci * args.n_bw
    print(f"Total region-tracks to extract: {total_region_tracks:,}\n")

    # ── warm-up (OS cache, JIT) ────────────────────────────────────────────
    print("Warm-up pass (small slice, not timed) ...")
    _ = signal(loci[:64], bigwigs[:2], n_bins=args.n_bins, flank=args.flank,
               agg=args.agg, workers=1, progress=False, verbose=False,
               exact=args.exact)

    results = []

    # ── sequential (workers=1) ─────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("SEQUENTIAL (workers=1)")
    print("-" * 78)
    t0 = time.perf_counter()
    cube_seq = signal(loci, bigwigs,
                      n_bins=args.n_bins, flank=args.flank, agg=args.agg,
                      workers=1, progress=True, verbose=False,
                      exact=args.exact)
    t_seq = time.perf_counter() - t0
    print(f"  time: {_fmt_time(t_seq)}   throughput: {_fmt_rate(total_region_tracks, t_seq)}")
    results.append(("sequential", 1, t_seq))
    baseline = t_seq

    # ── multiprocessing sweep ──────────────────────────────────────────────
    print("\n" + "=" * 78)
    print(f"MULTIPROCESSING (ProcessPoolExecutor + SharedMemory)")
    print("-" * 78)
    for w in args.workers:
        if w == 1:
            continue  # already measured as sequential
        t0 = time.perf_counter()
        cube = signal(loci, bigwigs,
                      n_bins=args.n_bins, flank=args.flank, agg=args.agg,
                      workers=w, progress=False, verbose=False,
                      exact=args.exact)
        t = time.perf_counter() - t0
        speedup = baseline / t
        results.append((f"mp w={w}", w, t))
        print(f"  workers={w:>3d}  "
              f"time: {_fmt_time(t)}   "
              f"{_fmt_rate(total_region_tracks, t)}   "
              f"speedup: {speedup:4.2f}×")

    # ── threading sweep (expected regression) ──────────────────────────────
    if not args.skip_thread:
        print("\n" + "=" * 78)
        print(f"THREADING (ThreadPoolExecutor — shown to demonstrate anti-scaling)")
        print("-" * 78)
        for w in args.thread_workers:
            t0 = time.perf_counter()
            _ = _thread_extract(loci, bigwigs,
                                n_bins=args.n_bins, flank=args.flank,
                                agg=args.agg, exact=args.exact, workers=w)
            t = time.perf_counter() - t0
            speedup = baseline / t
            results.append((f"thread w={w}", w, t))
            print(f"  threads={w:>3d}  "
                  f"time: {_fmt_time(t)}   "
                  f"{_fmt_rate(total_region_tracks, t)}   "
                  f"speedup: {speedup:4.2f}×")

    # ── backend comparison (sequential, small subset) ──────────────────────
    if not args.skip_backend_compare:
        print("\n" + "=" * 78)
        print("BACKEND COMPARISON (sequential, 1 worker, small subset)")
        print(f"  subset: {args.cmp_n_bw} BW × {args.cmp_n_loci} loci "
              f"(pyBigWig is ~50-100× slower than pybigtools; kept small).")
        print("-" * 78)
        cmp_rows, cmp_loci, cmp_bw = _run_backend_compare(
            loci, bigwigs,
            n_bins=args.n_bins, flank=args.flank, agg=args.agg,
            n_loci_cap=args.cmp_n_loci, n_bw_cap=args.cmp_n_bw,
        )
        n_rt_cmp = cmp_loci * cmp_bw
        pbt_t = next((t for n, t, _ in cmp_rows if n.startswith("pybigtools")), None)
        for name, t, _ in cmp_rows:
            if t == float('inf'):
                print(f"  {name:<24s} (skipped)")
                continue
            rel = (t / pbt_t) if pbt_t else 1.0
            print(f"  {name:<24s} {_fmt_time(t):>10s}  "
                  f"{_fmt_rate(n_rt_cmp, t):>14s}  "
                  f"{rel:6.2f}× pybigtools time")

    # ── rendering smoke tests ──────────────────────────────────────────────
    if not args.skip_render:
        print("\n" + "=" * 78)
        print("RENDERING SMOKE TESTS (small inputs, timed; output saved to --render-dir)")
        print("-" * 78)
        out_dir = Path(args.render_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        render_results = _run_render_smoke(loci, bigwigs, cube_seq,
                                           n_bins=args.n_bins,
                                           flank=args.flank,
                                           out_dir=out_dir)
        for name, ok, t, path, err in render_results:
            status = "OK  " if ok else "FAIL"
            tag = f"{status}  {_fmt_time(t)}"
            if ok:
                print(f"  {name:<20s} {tag}  -> {path}")
            else:
                print(f"  {name:<20s} {tag}  -> {err}")

    # ── summary table ──────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("-" * 78)
    print(f"{'approach':<16s}  {'workers':>7s}  {'time':>10s}  "
          f"{'throughput':>14s}  {'speedup':>7s}")
    print("-" * 78)
    for name, w, t in results:
        print(f"{name:<16s}  {w:>7d}  {_fmt_time(t):>10s}  "
              f"{_fmt_rate(total_region_tracks, t):>14s}  "
              f"{baseline/t:>6.2f}×")
    print("=" * 78)


if __name__ == "__main__":
    main()
