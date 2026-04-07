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
TRIGGER_SAMPLE = 0 #17 for the for the actual CARE output, but we set to 0 for the shifted traces we saved

CSV_FIELDS = [
    "pid", "energy_string", "radius", "seed",
    "zen", "az", "height",
    "tel_x", "tel_y", "tel_z", "tel_r",
    # ── Existing ──
    "max_pe", "time_at_max_pe_ns", "avg_pe", "total_pe",
    # ── Image shape (Hillas-like) ──
    "n_hit_pixels",          # number of pixels above threshold
    "image_size_pe",         # sum of PE in hit pixels (cleaned)
    "frac_in_brightest",     # fraction of total_pe in brightest pixel
    "concentration_2",       # fraction in 2 brightest pixels
    # ── Timing ──
    "pulse_width_ns",        # FWHM of summed trace
    "rise_time_ns",          # 10%→90% of peak in summed trace
    "time_spread_ns",        # RMS of per-pixel peak times (hit pixels only)
    "time_gradient",         # linear slope of peak-time vs pixel index (proxy for direction)
    # ── Trace shape ──
    "peak_to_charge",        # max_pe / (total_pe / n_hit_pixels), peakedness
    "baseline_rms_pe",       # RMS of first few samples (noise estimate)
    "file_found",
]


def _fmt(val):
    s = str(val)
    return s[:-2] if s.endswith('.0') else s

def _fmt_angle(val):
    f = float(val)
    return f"{f:.1f}"


def make_path(base_path, pid, energy_str, zen, az, h, x, y, z, r, s):
    return (f"{base_path}/Muon_pid{pid}_E{energy_str}_R{r}/"
            f"pdg{pid}_E{energy_str}_r{r}_s{s}/"
            f"zen{_fmt_angle(zen)}/"
            f"az{_fmt_angle(az)}/"
            f"h{_fmt(h)}/"
            f"x{_fmt(x)}_y{_fmt(y)}_z{_fmt(z)}/"
            f"CARE/cherenkov_hits.root")


def read_metrics(filepath, pe_threshold=1.0):
    try:
        with uproot.open(filepath) as f:
            tree = f['Events/T0;1']
            available = set(tree.keys())
            fadc_branches = [f"vFADCTraces{i}" for i in range(256)
                             if f"vFADCTraces{i}" in available]

            if not fadc_branches:
                return {k: 0.0 for k in CSV_FIELDS if k not in
                        ("pid","energy_string","radius","seed","zen","az",
                         "height","tel_x","tel_y","tel_z","tel_r","file_found")}, 1

            arrays = tree.arrays(fadc_branches, library="np")
            traces = np.array([arrays[b][0] for b in fadc_branches], dtype=np.float32)
            traces_pe = np.clip((traces - BASELINE_DC) / DC_TO_PE, 0.0, None)

            n_pix, n_samples = traces_pe.shape
            time_ns = (np.arange(n_samples) - TRIGGER_SAMPLE) * SAMPLES_PER_NS

            # ── Per-pixel quantities ─────────────────────────────────────
            pixel_max = traces_pe.max(axis=1)            # peak PE per pixel
            pixel_charge = traces_pe.sum(axis=1)         # integrated charge per pixel
            pixel_peak_time = time_ns[traces_pe.argmax(axis=1)]  # time of peak per pixel

            hit_mask = pixel_max >= pe_threshold
            n_hit = int(hit_mask.sum())

            # ── Summed trace ─────────────────────────────────────────────
            summed = traces_pe.sum(axis=0)
            peak_step = int(np.argmax(summed))
            max_pe = float(traces_pe[:, peak_step].max())
            time_at_max_pe = float(time_ns[peak_step])

            total_pe = float(traces_pe.sum())
            nonzero = traces_pe > 0
            avg_pe = float(traces_pe[nonzero].mean()) if nonzero.any() else 0.0

            # ── Image size (cleaned) ─────────────────────────────────────
            image_size = float(pixel_charge[hit_mask].sum()) if n_hit > 0 else 0.0

            # ── Concentration ────────────────────────────────────────────
            sorted_charge = np.sort(pixel_charge)[::-1]
            frac_brightest = float(sorted_charge[0] / total_pe) if total_pe > 0 else 0.0
            conc_2 = float(sorted_charge[:2].sum() / total_pe) if total_pe > 0 else 0.0

            # ── Pulse width (FWHM of summed trace) ───────────────────────
            half_max = summed[peak_step] / 2.0
            above_half = np.where(summed >= half_max)[0]
            if len(above_half) >= 2:
                pulse_width = float((above_half[-1] - above_half[0]) * SAMPLES_PER_NS)
            else:
                pulse_width = 0.0

            # ── Rise time (10%→90% of summed trace peak) ────────────────
            peak_val = summed[peak_step]
            t10 = np.where(summed[:peak_step+1] >= 0.1 * peak_val)[0]
            t90 = np.where(summed[:peak_step+1] >= 0.9 * peak_val)[0]
            if len(t10) > 0 and len(t90) > 0:
                rise_time = float((t90[0] - t10[0]) * SAMPLES_PER_NS)
            else:
                rise_time = 0.0

            # ── Time spread (RMS of hit-pixel peak times) ────────────────
            if n_hit >= 2:
                hit_times = pixel_peak_time[hit_mask]
                time_spread = float(np.std(hit_times))
            else:
                time_spread = 0.0

            # ── Time gradient (linear fit of peak time vs pixel index) ───
            if n_hit >= 3:
                hit_idx = np.where(hit_mask)[0].astype(float)
                hit_times = pixel_peak_time[hit_mask]
                coeffs = np.polyfit(hit_idx, hit_times, 1)
                time_gradient = float(coeffs[0])  # ns per pixel
            else:
                time_gradient = 0.0

            # ── Peak-to-charge ratio ─────────────────────────────────────
            mean_charge = image_size / n_hit if n_hit > 0 else 1.0
            peak_to_charge = float(max_pe / mean_charge) if mean_charge > 0 else 0.0

            # ── Baseline RMS (noise from first 5 samples) ────────────────
            n_baseline = min(5, n_samples)
            baseline_rms = float(np.std(traces_pe[:, :n_baseline]) )

            return (max_pe, time_at_max_pe, avg_pe, total_pe,
                    n_hit, image_size, frac_brightest, conc_2,
                    pulse_width, rise_time, time_spread, time_gradient,
                    peak_to_charge, baseline_rms, 1)

    except Exception:
        return (0.0, 0.0, 0.0, 0.0,
                0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--energy-str", required=True, type=str,
                    help="Energy string as it appears in directory names, e.g. '2.81838e5'")
    parser.add_argument("--radius", type=int, required=True)
    parser.add_argument("--tel-y", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--base-path", required=True)
    args = parser.parse_args()

    print(f"DEBUG: base_path = {args.base_path}")
    print(f"DEBUG: pid={args.pid} E_mag={args.energy_str} R={args.radius} seed={args.seed} tel_y={args.tel_y}")

    with open(args.chunk_file) as f:
        combos = json.load(f)

    print(f"DEBUG: Loaded {len(combos)} combos from {args.chunk_file}")
    if combos:
        # Show first combo's constructed path
        c = combos[0]
        sample_path = make_path(args.base_path, args.pid, args.energy_str,
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

            path = make_path(args.base_path, args.pid, args.energy_str,
                             zen, az, h, x_i, args.tel_y, z_i,
                             args.radius, args.seed)

            (max_pe, time_at_max, avg_pe, total_pe, n_hit, image_size, frac_brightest, conc_2, pulse_width, rise_time, time_spread, time_gradient, peak_to_charge, baseline_rms, found) = read_metrics(path)
            if found:
                found_count += 1
            else:
                not_found_count += 1
                # Print first 5 missing paths
                if not_found_count <= 5:
                    print(f"DEBUG: NOT FOUND [{not_found_count}]: {path}")
            writer.writerow({
                "pid":               args.pid,
                "energy_string":     args.energy_str,
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


            


    print(f"Wrote {len(combos)} rows to {args.output}")
    print(f"DEBUG SUMMARY: found={found_count}, not_found={not_found_count}")

if __name__ == "__main__":
    main()
