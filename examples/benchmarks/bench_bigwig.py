#!/usr/bin/env python3
"""
Benchmark BigWig readers: pure-Python vs pyBigWig vs pybigtools.

Usage:
    python bench_bigwig.py <bigwig_file>
    python bench_bigwig.py <bigwig_file> --n-regions 5000 --n-bins 200
    python bench_bigwig.py <bigwig_file> --threads 1 2 4 8
"""
from __future__ import annotations
import argparse
import sys
import time
import threading
from pathlib import Path

import numpy as np


# ── helpers ──────────────────────────────────────────────────────────────────

def _make_regions(chroms: dict, n: int, flank: int = 3000, seed: int = 42):
    """Generate random (chrom, start, end) regions."""
    rng = np.random.default_rng(seed)
    pool = [(c, s) for c, s in chroms.items() if s > 2 * flank]
    if not pool:
        raise RuntimeError("No chromosomes large enough for the given flank")
    regions = []
    for _ in range(n):
        c, sz = pool[rng.integers(len(pool))]
        center = rng.integers(flank, sz - flank)
        regions.append((c, int(center - flank), int(center + flank)))
    return regions


def _fmt(sec: float) -> str:
    if sec < 1:
        return f"{sec*1000:.1f} ms"
    return f"{sec:.3f} s"


def _rate(n: int, sec: float) -> str:
    return f"{n/sec:,.0f}/s"


# ── benchmarks ───────────────────────────────────────────────────────────────

def bench_open(opener, path, n=100):
    """Time: open + chroms + close."""
    # warm up
    f = opener(path); f.chroms(); f.close()
    t0 = time.perf_counter()
    for _ in range(n):
        f = opener(path)
        f.chroms()
        f.close()
    elapsed = time.perf_counter() - t0
    return elapsed / n


def bench_stats(opener, path, regions, n_bins=200, stat='mean'):
    """Time: sequential stats() calls."""
    f = opener(path)
    # warm up
    c, s, e = regions[0]
    f.stats(c, s, e, nBins=n_bins, type=stat)

    t0 = time.perf_counter()
    for c, s, e in regions:
        f.stats(c, s, e, nBins=n_bins, type=stat)
    elapsed = time.perf_counter() - t0
    f.close()
    return elapsed


def bench_stats_compat(opener, path, regions, n_bins=200, stat='mean',
                       use_new_api=False):
    """Time: sequential stats() with new-style kwargs."""
    f = opener(path)
    c, s, e = regions[0]
    if use_new_api:
        f.stats(c, s, e, n_bins=n_bins, stat=stat)
    else:
        f.stats(c, s, e, nBins=n_bins, type=stat)

    t0 = time.perf_counter()
    for c, s, e in regions:
        if use_new_api:
            f.stats(c, s, e, n_bins=n_bins, stat=stat)
        else:
            f.stats(c, s, e, nBins=n_bins, type=stat)
    elapsed = time.perf_counter() - t0
    f.close()
    return elapsed


def bench_values(opener, path, regions):
    """Time: sequential per-base values() calls."""
    f = opener(path)
    c, s, e = regions[0]
    f.values(c, s, e)

    t0 = time.perf_counter()
    for c, s, e in regions:
        f.values(c, s, e)
    elapsed = time.perf_counter() - t0
    f.close()
    return elapsed


def bench_threaded(opener, path, regions, n_bins, stat, n_threads):
    """Time: parallel stats() across threads."""
    chunk = len(regions) // n_threads
    errors = []

    def worker(my_regions):
        try:
            f = opener(path)
            for c, s, e in my_regions:
                f.stats(c, s, e, nBins=n_bins, type=stat)
            f.close()
        except Exception as exc:
            errors.append(exc)

    t0 = time.perf_counter()
    threads = []
    for i in range(n_threads):
        lo = i * chunk
        hi = lo + chunk if i < n_threads - 1 else len(regions)
        t = threading.Thread(target=worker, args=(regions[lo:hi],))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    elapsed = time.perf_counter() - t0

    if errors:
        raise errors[0]
    return elapsed


# ── correctness check ────────────────────────────────────────────────────────

def check_correctness(openers: dict, path: str, regions: list,
                      n_bins: int = 200):
    """Compare stats() output across readers."""
    ref_name = None
    ref_results = None

    for name, (opener, use_new) in openers.items():
        f = opener(path)
        results = []
        for c, s, e in regions[:50]:  # check first 50
            if use_new:
                r = f.stats(c, s, e, n_bins=n_bins, stat='mean')
            else:
                r = f.stats(c, s, e, nBins=n_bins, type='mean')
            results.append(np.array([0.0 if x is None else float(x) for x in r]))
        f.close()

        if ref_results is None:
            ref_name = name
            ref_results = results
            continue

        diffs = []
        for ref, cur in zip(ref_results, results):
            ref_a = np.array(ref, dtype=np.float64)
            cur_a = np.array(cur, dtype=np.float64)
            if ref_a.shape != cur_a.shape:
                diffs.append(np.inf)
            else:
                diffs.append(np.max(np.abs(ref_a - cur_a)))
        max_diff = max(diffs)
        mean_diff = np.mean(diffs)
        print(f"  {ref_name} vs {name}: max_diff={max_diff:.6g}  "
              f"mean_diff={mean_diff:.6g}  "
              f"{'PASS' if max_diff < 0.01 else 'WARN (>0.01)'}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Benchmark BigWig readers")
    parser.add_argument("bigwig", help="Path to a .bw file")
    parser.add_argument("--n-regions", type=int, default=1000)
    parser.add_argument("--n-bins", type=int, default=200)
    parser.add_argument("--flank", type=int, default=3000)
    parser.add_argument("--threads", type=int, nargs='+', default=[1, 2, 4])
    parser.add_argument("--stat", default="mean")
    parser.add_argument("--skip-correctness", action="store_true")
    args = parser.parse_args()

    path = str(Path(args.bigwig).resolve())
    print(f"BigWig: {path}")
    print(f"Regions: {args.n_regions}  Bins: {args.n_bins}  "
          f"Flank: {args.flank}  Stat: {args.stat}\n")

    # ── discover backends ────────────────────────────────────────────────
    backends = {}  # name → (opener, use_new_api)

    try:
        import pyBigWig
        backends['pyBigWig'] = (pyBigWig.open, False)
        print("  [+] pyBigWig available")
    except ImportError:
        print("  [-] pyBigWig not installed")

    try:
        import pybigtools
        def _pbt_open(p):
            from genomeblocks.signal import _PyBigToolsHandle
            return _PyBigToolsHandle(pybigtools.open(p, "r"))
        backends['pybigtools'] = (_pbt_open, True)
        print("  [+] pybigtools available (stats via values→bin)")
    except ImportError:
        print("  [-] pybigtools not installed")

    try:
        from genomeblocks import bigwig
        backends['pure-python'] = (bigwig.open, True)
        print("  [+] pure-python (genomeblocks.bigwig) available")
    except Exception as exc:
        print(f"  [-] pure-python import failed: {exc}")

    if not backends:
        print("\nNo backends available!")
        sys.exit(1)

    # get chroms from first available backend
    first_opener = list(backends.values())[0][0]
    f = first_opener(path)
    chroms = f.chroms()
    f.close()
    print(f"\nChromosomes: {len(chroms)}  "
          f"(largest: {max(chroms.values()):,} bp)\n")

    regions = _make_regions(chroms, args.n_regions, args.flank)
    n = args.n_regions

    # ── correctness ──────────────────────────────────────────────────────
    if not args.skip_correctness and len(backends) > 1:
        print("=" * 60)
        print("CORRECTNESS CHECK (first 50 regions, stats mean)")
        print("-" * 60)
        check_correctness(backends, path, regions, args.n_bins)
        print()

    # ── open + chroms ────────────────────────────────────────────────────
    print("=" * 60)
    print(f"OPEN + CHROMS  (100 iterations)")
    print("-" * 60)
    for name, (opener, _) in backends.items():
        try:
            t = bench_open(opener, path, n=100)
            print(f"  {name:20s}  {_fmt(t):>10s} / open")
        except Exception as e:
            print(f"  {name:20s}  ERROR: {e}")
    print()

    # ── sequential stats ─────────────────────────────────────────────────
    print("=" * 60)
    print(f"SEQUENTIAL STATS  ({n} regions, {args.n_bins} bins)")
    print("-" * 60)
    for name, (opener, use_new) in backends.items():
        try:
            t = bench_stats_compat(opener, path, regions, args.n_bins,
                                   args.stat, use_new_api=use_new)
            print(f"  {name:20s}  {_fmt(t):>10s}  "
                  f"({_rate(n, t)} regions/s)")
        except Exception as e:
            print(f"  {name:20s}  ERROR: {e}")
    print()

    # ── sequential values (base-pair) ────────────────────────────────────
    n_val = min(200, n)
    val_regions = regions[:n_val]
    print("=" * 60)
    print(f"BASE-PAIR VALUES  ({n_val} regions, {args.flank*2} bp each)")
    print("-" * 60)
    for name, (opener, _) in backends.items():
        try:
            t = bench_values(opener, path, val_regions)
            print(f"  {name:20s}  {_fmt(t):>10s}  "
                  f"({_rate(n_val, t)} regions/s)")
        except Exception as e:
            print(f"  {name:20s}  ERROR: {e}")
    print()

    # ── threaded stats ───────────────────────────────────────────────────
    print("=" * 60)
    print(f"THREADED STATS  ({n} regions, {args.n_bins} bins)")
    print("-" * 60)
    for n_t in args.threads:
        print(f"\n  --- {n_t} thread(s) ---")
        for name, (opener, _) in backends.items():
            try:
                t = bench_threaded(opener, path, regions, args.n_bins,
                                   args.stat, n_t)
                print(f"    {name:20s}  {_fmt(t):>10s}  "
                      f"({_rate(n, t)} regions/s)")
            except Exception as e:
                print(f"    {name:20s}  ERROR: {e}")
    print()


if __name__ == "__main__":
    main()
