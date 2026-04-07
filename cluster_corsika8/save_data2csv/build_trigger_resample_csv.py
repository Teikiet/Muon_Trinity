#!/usr/bin/env python3
"""
build_trigger_resample_csv.py

Reads the original merged CSV, filters for max_pe >= threshold (triggered events),
then searches for all resample (tel_x, tel_z) positions for those triggered events,
reads their CARE output, and writes a combined CSV.

Usage:
  python3 build_trigger_resample_csv.py \
    --base-path /scratch/general/vast/u1520754/muon_sim_chain_tree \
    --pid 13 --energy-str 3.54813e5 --radius 5 --tel-y 0 --seed 1 \
    --pe-threshold 20 \
    --tel-x-values -5 -2 -1 -0.3 0 0.3 1 2 5 \
    --tel-z-values -5 -2 -1 -0.3 0 0.3 1 2 5
"""

import argparse
import csv
import os
import sys
import glob
import numpy as np

try:
    import uproot
except ImportError:
    print("ERROR: uproot not available. pip install uproot awkward", file=sys.stderr)
    sys.exit(1)

CSV_FIELDS = [
    "pid", "energy_string", "radius", "seed",
    "zen", "az", "height",
    "tel_x", "tel_y", "tel_z", "tel_r",
    "max_pe", "time_at_max_pe_ns", "avg_pe", "total_pe",
    "n_hit_pixels", "image_size_pe", "frac_in_brightest",
    "concentration_2", "pulse_width_ns", "rise_time_ns",
    "time_spread_ns", "time_gradient", "peak_to_charge",
    "baseline_rms_pe", "file_found",
]


def make_path(base_path, pid, energy_str, zen, az, h, tel_x, tel_y, tel_z, radius, seed):
    energy_group = f"Muon_pid{pid}_E{energy_str}_R{radius}"
    tree_root = f"pdg{pid}_E{energy_str}_r{radius}_s{seed}"

    def fmt(v):
        f = float(v)
        if f == int(f):
            return str(int(f))
        return str(v)

    zen_s = fmt(zen)
    az_s = fmt(az)
    h_s = fmt(h)
    x_s = fmt(tel_x)
    y_s = fmt(tel_y)
    z_s = fmt(tel_z)

    pos_dir = f"x{x_s}_y{y_s}_z{z_s}"
    return os.path.join(base_path, energy_group, tree_root,
                        f"zen{zen_s}", f"az{az_s}", f"h{h_s}",
                        pos_dir, "CARE")



def read_metrics(care_dir):
    """Read metrics from CARE root file. Returns tuple of metrics + found bool."""
    defaults = (0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    if not os.path.isdir(care_dir):
        return defaults

    root_files = glob.glob(os.path.join(care_dir, "*.root"))
    if not root_files:
        return defaults

    root_file = root_files[0]

    BASELINE_DC = 3800.0
    DC_TO_PE = 5.6
    SAMPLES_PER_NS = 2.0
    TRIGGER_SAMPLE = 20

    try:
        f = uproot.open(root_file)
        tree = f["Events/T0"]
        if tree.num_entries == 0:
            return defaults

        available = set(tree.keys())
        fadc_branches = [f"vFADCTraces{i}" for i in range(256)
                         if f"vFADCTraces{i}" in available]

        if not fadc_branches:
            return defaults

        arrays = tree.arrays(fadc_branches, library="np")
        traces = np.array([arrays[b][0] for b in fadc_branches], dtype=np.float32)
        traces_pe = np.clip((traces - BASELINE_DC) / DC_TO_PE, 0.0, None)

        n_pix, n_samples = traces_pe.shape
        time_ns = (np.arange(n_samples) - TRIGGER_SAMPLE) * SAMPLES_PER_NS

        pixel_max = traces_pe.max(axis=1)
        pixel_charge = traces_pe.sum(axis=1)
        pixel_peak_time = time_ns[traces_pe.argmax(axis=1)]

        pe_threshold = 1.0
        hit_mask = pixel_max >= pe_threshold
        n_hit = int(hit_mask.sum())

        summed = traces_pe.sum(axis=0)
        peak_step = int(np.argmax(summed))
        max_pe = float(traces_pe[:, peak_step].max())
        time_at_max_pe = float(time_ns[peak_step])

        total_pe = float(traces_pe.sum())
        nonzero = traces_pe > 0
        avg_pe = float(traces_pe[nonzero].mean()) if nonzero.any() else 0.0

        image_size = float(pixel_charge[hit_mask].sum()) if n_hit > 0 else 0.0

        sorted_charge = np.sort(pixel_charge)[::-1]
        frac_brightest = float(sorted_charge[0] / total_pe) if total_pe > 0 else 0.0
        conc_2 = float(sorted_charge[:2].sum() / total_pe) if total_pe > 0 else 0.0

        half_max = summed[peak_step] / 2.0
        above_half = np.where(summed >= half_max)[0]
        if len(above_half) >= 2:
            pulse_width = float((above_half[-1] - above_half[0]) * SAMPLES_PER_NS)
        else:
            pulse_width = 0.0

        peak_val = summed[peak_step]
        t10 = np.where(summed[:peak_step+1] >= 0.1 * peak_val)[0]
        t90 = np.where(summed[:peak_step+1] >= 0.9 * peak_val)[0]
        if len(t10) > 0 and len(t90) > 0:
            rise_time = float((t90[0] - t10[0]) * SAMPLES_PER_NS)
        else:
            rise_time = 0.0

        if n_hit >= 2:
            time_spread = float(np.std(pixel_peak_time[hit_mask]))
        else:
            time_spread = 0.0

        if n_hit >= 3:
            hit_idx = np.where(hit_mask)[0].astype(float)
            coeffs = np.polyfit(hit_idx, pixel_peak_time[hit_mask], 1)
            time_gradient = float(coeffs[0])
        else:
            time_gradient = 0.0

        mean_charge = image_size / n_hit if n_hit > 0 else 1.0
        peak_to_charge = float(max_pe / mean_charge) if mean_charge > 0 else 0.0

        n_baseline = min(5, n_samples)
        baseline_rms = float(np.std(traces_pe[:, :n_baseline]))

        return (max_pe, time_at_max_pe, avg_pe, total_pe, n_hit,
                image_size, frac_brightest, conc_2,
                pulse_width, rise_time, time_spread, time_gradient,
                peak_to_charge, baseline_rms, 1)

    except Exception as e:
        print(f"  WARN: Error reading {root_file}: {e}", file=sys.stderr)
        return defaults



def main():
    parser = argparse.ArgumentParser(description="Build triggered + resample CSV")
    parser.add_argument("--base-path", required=True)
    parser.add_argument("--pid", type=int, default=13)
    parser.add_argument("--energy-str", required=True)
    parser.add_argument("--radius", type=int, default=5)
    parser.add_argument("--tel-y", type=float, default=0)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--pe-threshold", type=float, default=20)
    parser.add_argument("--tel-x-values", nargs="+", type=float, required=True)
    parser.add_argument("--tel-z-values", nargs="+", type=float, required=True)
    args = parser.parse_args()

    energy_group = f"Muon_pid{args.pid}_E{args.energy_str}_R{args.radius}"
    csv_dir = os.path.join(args.base_path, energy_group, "csv_output")

    # Original merged CSV
    orig_csv = os.path.join(csv_dir,
        f"scan_care_pid{args.pid}_E{args.energy_str}_R{args.radius}_y{int(args.tel_y)}_s{args.seed}.csv")

    if not os.path.isfile(orig_csv):
        print(f"ERROR: Original CSV not found: {orig_csv}", file=sys.stderr)
        sys.exit(1)

    # Output CSV
    out_csv = os.path.join(csv_dir,
        f"scan_care_pid{args.pid}_E{args.energy_str}_R{args.radius}_y{int(args.tel_y)}_s{args.seed}_trigger.csv")

    print(f"Reading original CSV: {orig_csv}")
    print(f"PE threshold: {args.pe_threshold}")
    print(f"tel_x values: {args.tel_x_values}")
    print(f"tel_z values: {args.tel_z_values}")

    # Step 1: Read original CSV, filter triggered events
    triggered_events = []
    total_rows = 0
    with open(orig_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_rows += 1
            try:
                max_pe = float(row["max_pe"])
            except (ValueError, KeyError):
                continue
            if max_pe >= args.pe_threshold:
                triggered_events.append(row)

    print(f"Total rows in original CSV: {total_rows}")
    print(f"Triggered events (max_pe >= {args.pe_threshold}): {len(triggered_events)}")

    # Step 2: For each triggered event, collect base row + all resample positions
    all_tel_positions = []
    for tx in args.tel_x_values:
        for tz in args.tel_z_values:
            all_tel_positions.append((tx, tz))

    n_resample = len([p for p in all_tel_positions if not (p[0] == 0 and p[1] == 0)])
    print(f"Total tel positions: {len(all_tel_positions)} ({n_resample} non-base)")

    rows_written = 0
    resample_found = 0
    resample_missing = 0

    with open(out_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDS)
        writer.writeheader()

        for event in triggered_events:
            zen = event["zen"]
            az = event["az"]
            height = event["height"]

            # Strip .0 from height if present
            try:
                h_float = float(height)
                if h_float == int(h_float):
                    height = str(int(h_float))
            except ValueError:
                pass

            # Write the base event row as-is
            writer.writerow(event)
            rows_written += 1

            # Now search for all resample positions (skip base 0,0)
            for tel_x, tel_z in all_tel_positions:
                if tel_x == 0 and tel_z == 0:
                    continue

                care_dir = make_path(
                    args.base_path, args.pid, args.energy_str,
                    zen, az, height,
                    tel_x, args.tel_y, tel_z,
                    args.radius, args.seed
                )
                if tel_x == -5 and tel_z == -5 and resample_found + resample_missing <= 3:
                    print(f"  DEBUG path: {care_dir}")
                    print(f"  DEBUG exists: {os.path.isdir(care_dir)}")

                (max_pe, time_at_max, avg_pe, total_pe, n_hit,
                 image_size, frac_brightest, conc_2,
                 pulse_width, rise_time, time_spread,
                 time_gradient, peak_to_charge, baseline_rms,
                 found) = read_metrics(care_dir)

                if found:
                    resample_found += 1
                else:
                    resample_missing += 1
                    if resample_missing <= 10:
                        print(f"  MISSING: {care_dir}")

                tel_r = round(np.sqrt(tel_x**2 + tel_z**2), 4)

                writer.writerow({
                    "pid":               args.pid,
                    "energy_string":     args.energy_str,
                    "radius":            args.radius,
                    "seed":              args.seed,
                    "zen":               zen,
                    "az":                az,
                    "height":            height,
                    "tel_x":             tel_x,
                    "tel_y":             args.tel_y,
                    "tel_z":             tel_z,
                    "tel_r":             tel_r,
                    "max_pe":            round(max_pe, 4),
                    "time_at_max_pe_ns": round(time_at_max, 2),
                    "avg_pe":            round(avg_pe, 4),
                    "total_pe":          round(total_pe, 4),
                    "n_hit_pixels":      n_hit,
                    "image_size_pe":     round(image_size, 4),
                    "frac_in_brightest": round(frac_brightest, 4),
                    "concentration_2":   round(conc_2, 4),
                    "pulse_width_ns":    round(pulse_width, 2),
                    "rise_time_ns":      round(rise_time, 2),
                    "time_spread_ns":    round(time_spread, 2),
                    "time_gradient":     round(time_gradient, 4),
                    "peak_to_charge":    round(peak_to_charge, 4),
                    "baseline_rms_pe":   round(baseline_rms, 4),
                    "file_found":        found,
                })
                rows_written += 1

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Triggered base events: {len(triggered_events)}")
    print(f"Resample positions found: {resample_found}")
    print(f"Resample positions missing: {resample_missing}")
    print(f"Total rows written: {rows_written}")
    print(f"Output: {out_csv}")


if __name__ == "__main__":
    main()
