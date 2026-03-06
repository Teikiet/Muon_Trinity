#!/usr/bin/env python3
# save_CARE2csv_chunk_tree.py — process one chunk of parameter combinations

import uproot
import numpy as np
import csv
import json
import argparse

DC_TO_PE       = 24.1
BASELINE_DC    = 500.0
SAMPLES_PER_NS = 10
TRIGGER_SAMPLE = 17 #17 for the for the actual CARE output, but we set to 0 for the shifted traces we saved

CSV_FIELDS = [
    "pid", "energy_magnitude", "radius", "seed",
    "zen", "az", "height",
    "tel_x", "tel_y", "tel_z", "tel_r",
    "max_pe", "time_at_max_pe_ns", "avg_pe", "total_pe", "file_found",
]

def _fmt(val):
    s = str(val)
    return s[:-2] if s.endswith('.0') else s

def _fmt_angle(val):
    f = float(val)
    return str(int(f)) if f == int(f) else f"{f:.1f}"


def make_path(base_path, pid, E_mag, zen, az, h, x, y, z, r, s):
    return (f"{base_path}/Muon_pid{pid}_E1e{E_mag}_R{r}/"
            f"pdg{pid}_E1e{E_mag}_r{r}_s{s}/"
            f"zen{_fmt_angle(zen)}/"
            f"az{_fmt_angle(az)}/"
            f"h{_fmt(h)}/"
            f"x{_fmt(x)}_y{_fmt(y)}_z{_fmt(z)}/"
            f"CARE/cherenkov_hits.root")


def read_metrics(filepath):
    try:
        with uproot.open(filepath) as f:
            tree = f['Events/T0;1']
            available = set(tree.keys())
            fadc_branches = [f"vFADCTraces{i}" for i in range(256)
                             if f"vFADCTraces{i}" in available]

            if not fadc_branches:
                return 0.0, 0.0, 0.0, 0.0, 1

            arrays = tree.arrays(fadc_branches, library="np")
            traces = np.array([arrays[b][0] for b in fadc_branches], dtype=np.float32)
            traces_pe = np.clip((traces - BASELINE_DC) / DC_TO_PE, 0.0, None)

            n_samples = traces_pe.shape[1]
            time_ns = (np.arange(n_samples) - TRIGGER_SAMPLE) * SAMPLES_PER_NS

            total_per_step = traces_pe.sum(axis=0)
            peak_step = int(np.argmax(total_per_step))
            max_pe = float(traces_pe[:, peak_step].max())
            time_at_max_pe = float(time_ns[peak_step])

            nonzero_mask = traces_pe > 0
            avg_pe = float(traces_pe[nonzero_mask].mean()) if nonzero_mask.any() else 0.0
            total_pe = float(traces_pe.sum())

            return max_pe, time_at_max_pe, avg_pe, total_pe, 1
    except Exception:
        return 0.0, 0.0, 0.0, 0.0, 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--energy-mag", type=int, required=True)
    parser.add_argument("--radius", type=int, required=True)
    parser.add_argument("--tel-y", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--base-path", required=True)
    args = parser.parse_args()

    print(f"DEBUG: base_path = {args.base_path}")
    print(f"DEBUG: pid={args.pid} E_mag={args.energy_mag} R={args.radius} seed={args.seed} tel_y={args.tel_y}")

    with open(args.chunk_file) as f:
        combos = json.load(f)

    print(f"DEBUG: Loaded {len(combos)} combos from {args.chunk_file}")
    if combos:
        # Show first combo's constructed path
        c = combos[0]
        sample_path = make_path(args.base_path, args.pid, args.energy_mag,
                                float(c[0]), float(c[1]), float(c[2]),
                                float(c[3]), args.tel_y, float(c[4]),
                                args.radius, args.seed)
        print(f"DEBUG: First combo = {c}")
        print(f"DEBUG: First path  = {sample_path}")

        import os
        # Walk up the path to find where it breaks
        parts = sample_path.split("/")
        for i in range(1, len(parts) + 1):
            partial = "/".join(parts[:i])
            exists = os.path.exists(partial)
            if not exists:
                print(f"DEBUG: PATH BREAKS AT: {partial}")
                # Show what the parent contains
                parent = "/".join(parts[:i-1])
                if os.path.isdir(parent):
                    contents = os.listdir(parent)[:20]
                    print(f"DEBUG: Parent dir '{parent}' contains: {contents}")
                break
        else:
            print(f"DEBUG: Full path exists!")

    found_count = 0
    not_found_count = 0

    with open(args.output, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDS)
        writer.writeheader()

        for combo in combos:
            zen, az, h, x_i, z_i = float(combo[0]), float(combo[1]), float(combo[2]), float(combo[3]), float(combo[4])

            path = make_path(args.base_path, args.pid, args.energy_mag,
                             zen, az, h, x_i, args.tel_y, z_i,
                             args.radius, args.seed)

            max_pe, time_at_max, avg_pe, total_pe, found = read_metrics(path)

            if found:
                found_count += 1
            else:
                not_found_count += 1
                # Print first 5 missing paths
                if not_found_count <= 5:
                    print(f"DEBUG: NOT FOUND [{not_found_count}]: {path}")

            writer.writerow({
                "pid":               args.pid,
                "energy_magnitude":  args.energy_mag,
                "radius":            args.radius,
                "seed":              args.seed,
                "zen":               zen,
                "az":                az,
                "height":            h,
                "tel_x":             x_i,
                "tel_y":             args.tel_y,
                "tel_z":             z_i,
                "tel_r":             round(np.sqrt(x_i**2 + z_i**2), 4),
                "max_pe":            round(max_pe, 4),
                "time_at_max_pe_ns": round(time_at_max, 2),
                "avg_pe":            round(avg_pe, 4),
                "total_pe":          round(total_pe, 4),
                "file_found":        found,
            })

    print(f"Wrote {len(combos)} rows to {args.output}")
    print(f"DEBUG SUMMARY: found={found_count}, not_found={not_found_count}")

if __name__ == "__main__":
    main()
