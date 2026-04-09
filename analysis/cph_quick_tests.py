#!/usr/bin/env python3
"""
cph_quick_tests.py (terminal-friendly)

- Parses CORSIKA/GRISU-style .cph files and summarizes P-line columns
- Prints an ASCII histogram table for time column col[5] (GrOptics fPTime)

Usage:
  python cph_quick_tests.py cherenkov_hits.cph --max-photons 200000
  python cph_quick_tests.py cherenkov_hits.cph --hist-bins 80 --hist-width 60
  python cph_quick_tests.py cherenkov_hits.cph --hist-range "1549870,1553320"
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np


def parse_cph_photons(path: str | Path, max_photons: int | None = None) -> np.ndarray:
    rows = []
    n_skipped_bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == "*":
                continue
            if not line.startswith("P "):
                continue

            parts = line.split()
            try:
                vals = [float(x) for x in parts[1:]]
            except ValueError:
                n_skipped_bad += 1
                continue

            rows.append(vals)
            if max_photons is not None and len(rows) >= max_photons:
                break

    if not rows:
        raise RuntimeError("No P lines parsed. Is this the right file?")

    lengths = np.array([len(r) for r in rows], dtype=int)
    uniq, counts = np.unique(lengths, return_counts=True)
    common_len = int(uniq[np.argmax(counts)])
    filtered = [r for r in rows if len(r) == common_len]

    if len(filtered) != len(rows):
        print(f"[warn] Rows had varying lengths. Keeping {len(filtered)}/{len(rows)} with length={common_len}.")
    if n_skipped_bad:
        print(f"[warn] Skipped {n_skipped_bad} malformed P lines.")

    return np.array(filtered, dtype=np.float64)


def summarize_columns(arr: np.ndarray) -> None:
    n, m = arr.shape
    print(f"\nParsed photons: N={n}, columns per photon: M={m}\n")

    col_min = np.nanmin(arr, axis=0)
    col_max = np.nanmax(arr, axis=0)
    col_rng = col_max - col_min

    print("Per-column min/max/range (0-based column index within P numeric fields):")
    for j in range(m):
        print(f"  col[{j:2d}]  min={col_min[j]: .6g}   max={col_max[j]: .6g}   range={col_rng[j]: .6g}")


def ascii_histogram(values: np.ndarray, bins: int, width: int, value_range: tuple[float, float] | None) -> str:
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return "No finite values to histogram."

    if value_range is None:
        vmin = float(v.min())
        vmax = float(v.max())
    else:
        vmin, vmax = value_range
        if vmax <= vmin:
            raise ValueError("Invalid --hist-range: max must be > min")

    counts, edges = np.histogram(v, bins=bins, range=(vmin, vmax))
    maxc = int(counts.max()) if counts.size else 0
    total = int(counts.sum())
    lines = []
    # after qv is computed
    t01, t05, t50, t95, t99 = np.quantile(v, [0.01, 0.05, 0.5, 0.95, 0.99])
    lines.append("\nRobust spans (raw units):")
    lines.append(f"  (t99 - t01) = {t99 - t01:.6g}")
    lines.append(f"  (t95 - t05) = {t95 - t05:.6g}")
    # tail diagnostics
    tail_thr = t99
    tail_count = int(np.sum(v > tail_thr))
    lines.append(f"\nTail photons above t99 ({tail_thr:.6g}): {tail_count} / {v.size} ({tail_count/v.size:.3%})")
    lines.append(f"\nASCII histogram: time = col[5] (raw units), N={v.size}, bins={bins}")
    lines.append(f"range: [{vmin:.6g}, {vmax:.6g}]  span={vmax - vmin:.6g}")
    lines.append("")
    lines.append("bin_index | low_edge -> high_edge | count | frac   | bar")
    lines.append("-" * (len(lines[-1])))

    for i in range(bins):
        lo = edges[i]
        hi = edges[i + 1]
        c = int(counts[i])
        frac = (c / total) if total > 0 else 0.0
        bar_len = int(round((c / maxc) * width)) if maxc > 0 else 0
        bar = "#" * bar_len
        lines.append(f"{i:8d} | {lo: .6g} -> {hi: .6g} | {c:5d} | {frac:6.2%} | {bar}")

    # Basic quantiles for terminal debugging
    qs = [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0]
    qv = np.quantile(v, qs)
    lines.append("\nQuantiles (time col[5], raw units):")
    for q, val in zip(qs, qv):
        lines.append(f"  q={q:>5.2%} : {val:.6g}")

    return "\n".join(lines)


def parse_range(s: str) -> tuple[float, float]:
    # format: "min,max"
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2:
        raise ValueError("Range must be 'min,max'")
    return float(parts[0]), float(parts[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cph_file")
    ap.add_argument("--max-photons", type=int, default=None, help="limit number of P lines to read")

    # Histogram controls
    ap.add_argument("--hist-time", action="store_true", default=True,
                    help="print ASCII histogram for time col[5] (default: on)")
    ap.add_argument("--no-hist-time", action="store_false", dest="hist_time",
                    help="disable ASCII histogram output")
    ap.add_argument("--hist-bins", type=int, default=60, help="number of histogram bins")
    ap.add_argument("--hist-width", type=int, default=50, help="max bar width in characters")
    ap.add_argument("--hist-range", type=parse_range, default=None,
                    help="optional range 'min,max' in raw time units for col[5]")

    args = ap.parse_args()

    arr = parse_cph_photons(args.cph_file, max_photons=args.max_photons)
    summarize_columns(arr)

    if arr.shape[1] <= 5:
        raise RuntimeError(f"Expected at least 6 columns on P line; got M={arr.shape[1]}")

    if args.hist_time:
        t = arr[:, 5]  # GrOptics fPTime
        #print(ascii_histogram(t, bins=args.hist_bins, width=args.hist_width, value_range=args.hist_range))
        t0 = float(np.nanmin(t))
        t = t - t0
        print(f"\n[info] Relative histogram enabled: using (time - min_time). min_time={t0:.6g}")

        print(ascii_histogram(
            t, bins=args.hist_bins, width=args.hist_width, value_range=args.hist_range
        ))


if __name__ == "__main__":
    main()
