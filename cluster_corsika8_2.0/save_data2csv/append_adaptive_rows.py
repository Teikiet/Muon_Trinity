#!/usr/bin/env python3
"""
Helper script to merge adaptive-cycle CSV rows into per-(zenith, azimuth) CSVs.

One file per expected (zenith, azimuth) cell, named
``adaptive_pid{pid}_E{E}_R{r}_y{tel_y}_zen{zen:.1f}_az{az:.1f}.csv``.

Used by the adaptive controller's save-job wrap to either:
- APPEND new rows to existing cell CSV files (normal operation)
- REPLACE rows matched by key (--rerun-failed operation)

Key = (seed, energy_string, zen, az, height, tel_x, tel_z, injection_x, injection_y)
Locking: fcntl.flock on the target CSV file during read-modify-write.
"""

import argparse
import csv
import fcntl
import os
import sys
import tempfile
from pathlib import Path


def offset_key(value):
    """Match _fmt_offset() from save_CARE2csv_chunk_tree.py."""
    try:
        fval = float(value)
    except (ValueError, TypeError):
        fval = 0.0
    if abs(fval) < 1e-12:
        return "0"
    return f"{fval:.1f}"


def _row_angle(row, scan_key, fallback_key):
    """Return the expected (scanned) angle, falling back to the legacy column.

    The ``zen``/``az`` columns now hold the energy-weighted centre-of-mass
    direction of all particles at the observation level; the expected arrival is
    in ``scan_zen``/``scan_az``. Legacy CSVs have no ``scan_*`` columns (their
    ``zen``/``az`` held the scanned angles), hence the fallback.
    """
    for key in (scan_key, fallback_key):
        value = row.get(key)
        if value in (None, ""):
            continue
        try:
            return float(value)
        except (ValueError, TypeError):
            continue
    return 0.0


def row_key(row):
    """Canonical key for deduplication/replacement.

    Keyed on the *expected* (scanned) arrival angles so re-running a cell
    replaces its row regardless of the centre-of-mass direction recorded in
    ``zen``/``az``.
    """
    return (
        str(int(float(row.get("seed", 0)))),
        str(row.get("energy_string", "")).strip(),
        f"{_row_angle(row, 'scan_zen', 'zen'):.1f}",
        f"{_row_angle(row, 'scan_az', 'az'):.1f}",
        str(int(round(float(row.get("height", 0))))),
        offset_key(row.get("tel_x", 0.0)),
        offset_key(row.get("tel_z", 0.0)),
    )


def append_or_replace(target_path, rows, header, mode="append"):
    """
    Append rows to target_path, or replace matching rows by key.

    Uses flock to guard against concurrent writes.
    """
    target_path = Path(target_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.touch(exist_ok=True)

    with open(target_path, "r+", newline="", encoding="utf-8") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)

        if mode == "append":
            handle.seek(0, os.SEEK_END)
            is_empty = handle.tell() == 0
            if is_empty:
                writer = csv.DictWriter(handle, fieldnames=header)
                writer.writeheader()
            else:
                writer = csv.DictWriter(handle, fieldnames=header)
            for row in rows:
                writer.writerow({name: row.get(name, "") for name in header})
            handle.flush()
            os.fsync(handle.fileno())
        elif mode == "replace":
            handle.seek(0)
            reader = csv.DictReader(handle)
            existing_rows = list(reader) if reader.fieldnames else []

            existing = {}
            order = []
            for row in existing_rows:
                key = row_key(row)
                if key not in existing:
                    order.append(key)
                existing[key] = row

            for row in rows:
                key = row_key(row)
                if key not in existing:
                    order.append(key)
                existing[key] = row

            tmp_dir = target_path.parent
            with tempfile.NamedTemporaryFile(
                "w", delete=False, dir=str(tmp_dir), encoding="utf-8", newline=""
            ) as tmp:
                writer = csv.DictWriter(tmp, fieldnames=header)
                writer.writeheader()
                for key in order:
                    row = existing[key]
                    writer.writerow({name: row.get(name, "") for name in header})
                tmp_path_str = tmp.name

            os.replace(tmp_path_str, str(target_path))
        else:
            raise ValueError(f"Unknown mode: {mode!r}")


def split_rows_by_zen_az(rows):
    """Group rows by their expected (scanned) zenith AND azimuth.

    One output file per (zenith, azimuth) cell so the filename can carry both
    expected angles. ``_row_angle`` prefers ``scan_zen``/``scan_az`` and falls
    back to ``zen``/``az`` for CSVs written before those columns existed.
    """
    by_cell = {}
    for row in rows:
        zen_key = f"{_row_angle(row, 'scan_zen', 'zen'):.1f}"
        az_key = f"{_row_angle(row, 'scan_az', 'az'):.1f}"
        by_cell.setdefault((zen_key, az_key), []).append(row)
    return by_cell


def main():
    parser = argparse.ArgumentParser(
        description="Append/replace adaptive CSV rows to per-zenith CSVs."
    )
    parser.add_argument("--result-csv", action="append", required=True,
                        help="Path to a cycle result CSV (repeatable, e.g. one per "
                             "save batch; rows are merged into the destination).")
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--energy-str", required=True)
    parser.add_argument("--radius", type=int, required=True)
    parser.add_argument("--tel-y", type=float, required=True)
    parser.add_argument("--base-path", required=True,
                        help="Base path for simulation outputs.")
    parser.add_argument("--mode", choices=["append", "replace"],
                        default="append",
                        help="'append' for normal cycles; 'replace' for --rerun-failed.")
    args = parser.parse_args()

    result_csvs = [Path(p) for p in args.result_csv]
    header = None
    rows = []
    for rc in result_csvs:
        if not rc.exists():
            print(f"ERROR: result CSV not found: {rc}", file=sys.stderr)
            sys.exit(1)
        with rc.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if header is None:
                header = reader.fieldnames or []
            rows.extend(list(reader))

    if not rows:
        print("No rows to process.")
        return

    by_cell = split_rows_by_zen_az(rows)
    base_path = Path(args.base_path)

    for (zen_key, az_key), cell_rows in by_cell.items():
        target = (
            base_path
            / f"Muon_pid{args.pid}_E{args.energy_str}_R{args.radius}"
            / "csv_output"
            / f"adaptive_pid{args.pid}_E{args.energy_str}_R{args.radius}"
              f"_y{args.tel_y}_zen{zen_key}_az{az_key}.csv"
        )
        append_or_replace(target, cell_rows, header, mode=args.mode)
        print(f"  {args.mode} {len(cell_rows)} rows -> {target}")

    print(f"Done: {len(rows)} rows to {len(by_cell)} (zen, az) files.")


if __name__ == "__main__":
    main()