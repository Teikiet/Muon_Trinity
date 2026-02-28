#/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/analysis/save_CARE2csv.py
import uproot
import numpy as np
import csv
import os
from itertools import product

# ── Constants ─────────────────────────────────────────────────────────────────
DC_TO_PE    = 24.1
BASELINE_DC = 500.0
SAMPLES_PER_NS = 10        # 10 ns per FADC sample
TRIGGER_SAMPLE = 17        # trigger at sample 17 → t=0 at 170 ns

# ── Simulation parameters ─────────────────────────────────────────────────────
particle_id      = 13
energy_magnitude = 4
radius           = 5
y_tel            = 0
seed             = 5

ZENITHS   = [89.9, 89, 88, 87]
AZIMUTHS  = [268, 269, 270, 271, 272]
TEL_XS    = [-10, -7, -5, -2, -1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1, 2, 5, 7, 10]
TEL_ZS    = [-10, -7, -5, -2, -1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1, 2, 5, 7, 10]
HEIGHTS   = [5000, 10000, 15000, 20000, 30000]

BASE_PATH = "/scratch/general/vast/u1520754/muon_sim_chain"

# ── Output CSV name ───────────────────────────────────────────────────────────
OUT_CSV = (f"scan_care_pid{particle_id}_E1e{energy_magnitude}"
           f"_R{radius}_y{y_tel}_s{seed}.csv")

CSV_FIELDS = [
    "pid", "energy_magnitude", "radius", "seed",
    "zen", "az", "height",
    "tel_x", "tel_y", "tel_z", "tel_r",
    "max_pe",          # peak PE value (max over pixels at peak time step)
    "time_at_max_pe_ns",  # time (ns, relative to trigger) of that peak step
    "avg_pe",          # average PE over all pixels and all time steps (>0 only)
    "total_pe",        # sum of PE across all pixels and all time steps
    "file_found",      # 1 if file existed and was readable, 0 otherwise
]

def make_path(pid, E_mag, zen, az, h, x, y, z, r, s):
    return (f"{BASE_PATH}/Muon_pid{pid}_E1e{E_mag}_R{r}/"
            f"Tilt_pdg{pid}_E1e{E_mag}_zen{zen}_az{az}_h{h}"
            f"_x{x}_y{y}_z{z}_r{r}_s{s}/CARE/cherenkov_hits.root")

def read_metrics(filepath):
    """
    Open a CARE ROOT file and return:
        max_pe, time_at_max_pe_ns, avg_pe, total_pe, file_found
    All zero / 0 if the file is missing or broken.
    """
    try:
        with uproot.open(filepath) as f:
            tree = f['Events/T0;1']
            available = set(tree.keys())

            fadc_branches = [f"vFADCTraces{i}" for i in range(256)
                             if f"vFADCTraces{i}" in available]

            if not fadc_branches:
                return 0.0, 0.0, 0.0, 0.0, 1   # file found but no traces

            arrays = tree.arrays(fadc_branches, library="np")

            # Stack → shape (n_pixels, n_samples), convert to PE
            traces = np.array([arrays[b][0] for b in fadc_branches],
                              dtype=np.float32)
            traces_pe = np.clip((traces - BASELINE_DC) / DC_TO_PE, 0.0, None)

            # ── Time axis (ns, relative to trigger) ──────────────────────────
            n_samples  = traces_pe.shape[1]
            time_ns    = (np.arange(n_samples) - TRIGGER_SAMPLE) * SAMPLES_PER_NS

            # ── Peak: time step where summed-over-pixels PE is maximum ────────
            total_per_step  = traces_pe.sum(axis=0)          # (n_samples,)
            peak_step       = int(np.argmax(total_per_step))
            max_pe          = float(traces_pe[:, peak_step].max())
            time_at_max_pe  = float(time_ns[peak_step])

            # ── Average PE (over all pixels × all samples, excluding zeros) ───
            nonzero_mask = traces_pe > 0
            if nonzero_mask.any():
                avg_pe = float(traces_pe[nonzero_mask].mean())
            else:
                avg_pe = 0.0

            # ── Total PE (sum over every pixel and every sample) ──────────────
            total_pe = float(traces_pe.sum())

            return max_pe, time_at_max_pe, avg_pe, total_pe, 1

    except Exception as e:
        return 0.0, 0.0, 0.0, 0.0, 0


# ── Main scan ─────────────────────────────────────────────────────────────────
all_combos = list(product(ZENITHS, AZIMUTHS, HEIGHTS, TEL_XS, TEL_ZS))
total      = len(all_combos)

print(f"Scanning {total} parameter combinations...")
print(f"Output will be saved to: {OUT_CSV}\n")

with open(OUT_CSV, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDS)
    writer.writeheader()

    for done, (zen, az, h, x_i, z_i) in enumerate(all_combos, 1):

        if done % 50 == 0 or done == total:
            print(f"  Progress: {done}/{total}  ({100*done/total:.1f}%)", end="\r")

        path = make_path(particle_id, energy_magnitude,
                         zen, az, h,
                         x_i, y_tel, z_i,
                         radius, seed)

        max_pe, time_at_max, avg_pe, total_pe, found = read_metrics(path)
        r_off = round(np.sqrt(x_i**2 + z_i**2), 4)

        writer.writerow({
            "pid":               particle_id,
            "energy_magnitude":  energy_magnitude,
            "radius":            radius,
            "seed":              seed,
            "zen":               zen,
            "az":                az,
            "height":            h,
            "tel_x":             x_i,
            "tel_y":             y_tel,
            "tel_z":             z_i,
            "tel_r":             r_off,
            "max_pe":            round(max_pe,      4),
            "time_at_max_pe_ns": round(time_at_max, 2),
            "avg_pe":            round(avg_pe,      4),
            "total_pe":          round(total_pe,    4),
            "file_found":        found,
        })

print(f"\n\nDone. {total} rows written to '{OUT_CSV}'.")
print(f"Missing/broken files will have file_found=0 and all PE columns=0.")
