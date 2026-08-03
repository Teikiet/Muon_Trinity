#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import sys


def _fmt(val):
    s = str(val)
    return s[:-2] if s.endswith('.0') else s


def _fmt_angle(val):
    return f"{float(val):.1f}"


def _fmt_offset(val):
    f = float(val)
    if abs(f) < 1e-12:
        return "0"
    return f"{f:.1f}"


def build_source_path(base_path, row):
    pid = _fmt(row["pid"])
    energy_string = str(row["energy_string"])
    radius = _fmt(row["radius"])
    seed = _fmt(row["seed"])
    zen = _fmt_angle(row["zen"])
    az = _fmt_angle(row["az"])
    height = _fmt(row["height"])
    tel_x = _fmt_offset(row["tel_x"])
    tel_y = _fmt_offset(row["tel_y"])
    tel_z = _fmt_offset(row["tel_z"])

    return (
        f"{base_path}/Muon_pid{pid}_E{energy_string}_R{radius}/"
        f"pdg{pid}_E{energy_string}_r{radius}_s{seed}/"
        f"zen{zen}/"
        f"az{az}/"
        f"h{height}/"
        f"x{tel_x}_y{tel_y}_z{tel_z}/"
        f"CARE/cherenkov_hits.root"
    )


def build_dest_name(row):
    pid = _fmt(row["pid"])
    energy_string = str(row["energy_string"])
    seed = _fmt(row["seed"])
    zen = _fmt_angle(row["zen"])
    az = _fmt_angle(row["az"])
    height = _fmt(row["height"])
    tel_x = _fmt_offset(row["tel_x"])
    tel_y = _fmt_offset(row["tel_y"])
    tel_z = _fmt_offset(row["tel_z"])

    return (
        f"CARE_pid{pid}_E{energy_string}_zen{zen}_az{az}_h{height}"
        f"_x{tel_x}_y{tel_y}_z{tel_z}_s{seed}.root"
    )


def main():
    parser = argparse.ArgumentParser(description="Upload triggered CARE root files to an rclone destination.")
    parser.add_argument("--csv", required=True, help="Merged CSV produced by submit_save_care2csv_tree.sh")
    parser.add_argument("--base-path", required=True, help="Simulation base path used to reconstruct CARE file paths")
    parser.add_argument("--drive-dest", required=True, help="rclone destination, for example remote:folder")
    parser.add_argument("--threshold", type=float, default=20.0, help="Upload rows with max_pe at or above this value")
    args = parser.parse_args()

    if not os.path.exists(args.csv):
        print(f"ERROR: CSV not found: {args.csv}", file=sys.stderr)
        return 1

    with open(args.csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print(f"No rows found in {args.csv}")
        return 0

    uploaded = 0
    skipped_missing = 0
    skipped_nontriggered = 0

    for row in rows:
        try:
            max_pe = float(row.get("max_pe", "nan"))
        except ValueError:
            skipped_nontriggered += 1
            continue

        if max_pe < args.threshold:
            skipped_nontriggered += 1
            continue

        source_path = build_source_path(args.base_path, row)
        if not os.path.exists(source_path):
            print(f"MISSING: {source_path}")
            skipped_missing += 1
            continue

        dest_name = build_dest_name(row)
        dest_path = f"{args.drive_dest.rstrip('/')}/{dest_name}"

        print(f"UPLOAD: {source_path} -> {dest_path}")
        subprocess.run(["rclone", "copyto", "--ignore-existing", source_path, dest_path], check=True)
        uploaded += 1

    print(f"Done: uploaded={uploaded}, missing={skipped_missing}, nontriggered_or_invalid={skipped_nontriggered}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())