#!/usr/bin/env python3
# copy_triggered_events.py — copy CARE root files for triggered events (max_pe > threshold)

import pandas as pd
import shutil
import os
import argparse

def _fmt(val):
    s = str(val)
    return s[:-2] if s.endswith('.0') else s

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy-mag", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--threshold", type=float, default=20.0)
    parser.add_argument("--base-path", default="/scratch/general/vast/u1520754/muon_sim_chain")
    parser.add_argument("--dest-path", default="/scratch/general/vast/u1520754/Triggered_Muon_Events")
    parser.add_argument("--dry-run", action="store_true", help="Print what would be copied without copying")
    args = parser.parse_args()

    pid = 13
    radius = 5
    E_mag = args.energy_mag
    seed = args.seed

    csv_path = (f"{args.base_path}/Muon_pid{pid}_E1e{E_mag}_R{radius}/"
                f"csv_output/scan_care_pid{pid}_E1e{E_mag}_R{radius}_y0_s{seed}.csv")

    if not os.path.exists(csv_path):
        print(f"CSV not found: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    triggered = df[df["max_pe"] > args.threshold]

    print(f"Found {len(triggered)} triggered events out of {len(df)} total (threshold > {args.threshold} PE)")

    os.makedirs(args.dest_path, exist_ok=True)

    copied = 0
    missing = 0
    for _, row in triggered.iterrows():
        zen = _fmt(row["zen"])
        az = _fmt(row["az"])
        h = _fmt(row["height"])
        x = _fmt(row["tel_x"])
        y = _fmt(row["tel_y"])
        z = _fmt(row["tel_z"])

        src = (f"{args.base_path}/Muon_pid{pid}_E1e{E_mag}_R{radius}/"
               f"Tilt_pdg{pid}_E1e{E_mag}_zen{zen}_az{az}_h{h}"
               f"_x{x}_y{y}_z{z}_r{radius}_s{seed}/CARE/cherenkov_hits.root")

        dst_name = (f"Event_pid{pid}_E1e{E_mag}_zen{zen}_az{az}_h{h}"
                    f"_x{x}_y{y}_z{z}_s{seed}.root")
        dst = os.path.join(args.dest_path, dst_name)

        if not os.path.exists(src):
            print(f"  MISSING: {src}")
            missing += 1
            continue

        if args.dry_run:
            print(f"  WOULD COPY: {src} -> {dst}")
        else:
            shutil.copy2(src, dst)
            copied += 1

    if args.dry_run:
        print(f"\nDry run: {copied + len(triggered) - missing} files would be copied, {missing} missing")
    else:
        print(f"\nCopied {copied} files, {missing} missing")

if __name__ == "__main__":
    main()
