#!/usr/bin/env python3
# merge_csv_chunks.py — combine all chunk CSVs into one final CSV

import glob
import argparse
import csv
import os


def _row_key(row):
    return (
        row.get("seed", "").strip(),
        row.get("energy_string", "").strip(),
        row.get("zen", "").strip(),
        row.get("az", "").strip(),
        row.get("height", "").strip(),
        row.get("tel_x", "").strip(),
        row.get("tel_z", "").strip(),
    )


def _as_int(val):
    try:
        return int(float(val))
    except Exception:
        return 0


def _pick_row(existing, new):
    if existing is None:
        return new
    existing_deleted = _as_int(existing.get("deleted_low_pe", "0"))
    new_deleted = _as_int(new.get("deleted_low_pe", "0"))
    existing_found = _as_int(existing.get("file_found", "0"))
    new_found = _as_int(new.get("file_found", "0"))

    if existing_deleted == 1:
        return existing
    if new_deleted == 1:
        return new
    if existing_found == 1 and new_found == 0:
        return existing
    if existing_found == 0 and new_found == 1:
        return new
    return existing

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    chunk_files = sorted(glob.glob(f"{args.chunk_dir}/result_*.csv"))
    print(f"Found {len(chunk_files)} chunk files")

    if not chunk_files:
        if os.path.exists(args.output):
            print(f"No chunk files found; keeping existing {args.output}")
            return
        print("No chunk files found; nothing to merge")
        return

    chunk_header = None
    with open(chunk_files[0], "r", newline="") as f:
        reader = csv.DictReader(f)
        chunk_header = reader.fieldnames or []

    existing_rows = {}
    order = []
    existing_header = []
    if os.path.exists(args.output):
        with open(args.output, "r", newline="") as f:
            reader = csv.DictReader(f)
            existing_header = reader.fieldnames or []
            for row in reader:
                key = _row_key(row)
                if key in existing_rows:
                    existing_rows[key] = _pick_row(existing_rows[key], row)
                else:
                    existing_rows[key] = row
                    order.append(key)

    if chunk_header and existing_header and chunk_header != existing_header:
        header = chunk_header + [h for h in existing_header if h not in chunk_header]
    elif chunk_header:
        header = chunk_header
    else:
        header = existing_header

    new_rows = 0
    updated_rows = 0
    for cf in chunk_files:
        with open(cf, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = _row_key(row)
                if key in existing_rows:
                    chosen = _pick_row(existing_rows[key], row)
                    if chosen is row:
                        existing_rows[key] = row
                        updated_rows += 1
                else:
                    existing_rows[key] = row
                    order.append(key)
                    new_rows += 1

    total_rows = 0
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=header)
        writer.writeheader()
        for key in order:
            row = existing_rows.get(key, {})
            writer.writerow({h: row.get(h, "") for h in header})
            total_rows += 1

    print(
        f"Merged {total_rows} rows into {args.output} "
        f"(new={new_rows}, updated={updated_rows})"
    )

if __name__ == "__main__":
    main()
