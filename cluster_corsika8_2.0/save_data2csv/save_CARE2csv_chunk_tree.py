#!/usr/bin/env python3
# save_CARE2csv_chunk_tree.py — process one chunk of parameter combinations

import uproot
import pyarrow.parquet as pq
import numpy as np
import csv
import json
import argparse
import os
import shutil

DC_TO_PE       = 24.1
BASELINE_DC    = 500.0
SAMPLES_PER_NS = 10
TRIGGER_SAMPLE = 0 #17 for the for the actual CARE output, but we set to 0 for the shifted traces we saved

CSV_FIELDS = [
    "pid", "energy_string", "radius", "seed",
    "scan_zen", "scan_az",
    "zen", "az", "height",
    "tel_x", "tel_y", "tel_z", "tel_r",
    "x_mag", "y_mag", "z_mag",
    "injection_zenith_deg", "injection_azimuth_deg",
    "injection_x", "injection_y", "injection_height_m", "injection_z_minus_R_m",
    "survival_probability",
    "file_size_MB",
    "cph_photon_count",
    "cph_max_photons_10ns",
    "correction_x_m", "correction_y_m",
    "runtime_seconds",
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
    "max_muon_energy_GeV",
    "r68_m",
    "r99_m",
    "total_particles",
    "muon_component",
    "electron_component",
    "hadronic_component",
    "gamma_component",
    "other_component",
    "deleted_low_pe",
    "file_found",
    # ── Centre-of-mass direction (energy-weighted, ALL particles) ──
    # Appended at the end so merge_csv_chunks.py keeps old merged CSVs readable.
    "com_concentration",     # |sum(w*n)| / sum(w); 1.0 = perfectly collimated
    "max_e_muon_zen",        # legacy: direction of the highest-KE muon (was "zen")
    "max_e_muon_az",         # legacy: direction of the highest-KE muon (was "az")
    "com_x_m",               # energy-weighted centroid of particle positions
    "com_y_m",
    "com_z_m",
]


def _fmt(val):
    s = str(val)
    return s[:-2] if s.endswith('.0') else s

def _fmt_angle(val):
    f = float(val)
    return f"{f:.1f}"

def _fmt_offset(val):
    f = float(val)
    if abs(f) < 1e-12:
        return "0.0"
    return f"{f:.1f}"


def make_path(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix=0.0, iy=0.0):
    # NOTE: ix/iy are retained for call-site compatibility but are no longer
    # part of the tree layout (the injection is predicted, not scanned).
    return (f"{base_path}/Muon_pid{pid}_E{energy_str}_R{r}/"
            f"pdg{pid}_E{energy_str}_r{r}_s{s}/"
            f"zen{_fmt_angle(zen)}/"
            f"az{_fmt_angle(az)}/"
            f"h{_fmt(h)}/"
            f"x{_fmt_offset(x)}_y{_fmt_offset(y)}_z{_fmt_offset(z)}/"
            f"CARE/cherenkov_hits.root")


def make_run_dir(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix=0.0, iy=0.0):
    care_path = make_path(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix, iy)
    return os.path.dirname(os.path.dirname(care_path))


def make_corsika8_path(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix=0.0, iy=0.0):
    run_dir = make_run_dir(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix, iy)
    return os.path.join(run_dir, "CORSIKA8", "cherenkov_hits.dat")


def make_cph_path(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix=0.0, iy=0.0):
    run_dir = make_run_dir(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix, iy)
    return os.path.join(run_dir, "CIO", "cherenkov_hits.cph")


def make_particles_path(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix=0.0, iy=0.0):
    run_dir = make_run_dir(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix, iy)
    return os.path.join(run_dir, "corsika8_output", "particles", "particles.parquet")


def make_correction_paths(base_path, pid, energy_str, zen, az, h, x, y, z, r, s,
                          report_name="metadata.yaml", ix=0.0, iy=0.0):
    run_dir = make_run_dir(base_path, pid, energy_str, zen, az, h, x, y, z, r, s, ix, iy)
    candidates = [
        os.path.join(run_dir, report_name),
        os.path.join(run_dir, "telescope_position_correction.yaml")
        #os.path.join(run_dir, "correction_report_firstpass.json"),
    ]
    seen = set()
    ordered = []
    for path in candidates:
        if path in seen:
            continue
        seen.add(path)
        ordered.append(path)
    return ordered


def read_directory_size_mb(dir_path):
    try:
        total_bytes = 0
        for root, _, files in os.walk(dir_path):
            for name in files:
                fp = os.path.join(root, name)
                if os.path.islink(fp):
                    continue
                try:
                    total_bytes += os.path.getsize(fp)
                except OSError:
                    continue
        return total_bytes / (1024.0 * 1024.0)
    except Exception:
        return 0.0


def remove_corsika8_tables(run_dir):
    # Remove any leftover CORSIKA 8 tables file to keep run dirs clean.
    try:
        candidate = os.path.join(run_dir, "corsika8_tables.dat")
        if os.path.isfile(candidate):
            os.remove(candidate)
            return
        for root, _, files in os.walk(run_dir):
            if "corsika8_tables.dat" in files:
                os.remove(os.path.join(root, "corsika8_tables.dat"))
                return
    except Exception:
        return


def delete_run_dir(run_dir, base_path=None, dry_run=False):
    if not run_dir or run_dir in ("/", "."):
        return False
    if base_path:
        base_root = base_path.rstrip("/") + "/"
        if not run_dir.startswith(base_root):
            return False
    if dry_run:
        print(f"DRY-RUN: would delete low-PE run dir: {run_dir}")
        return False
    try:
        shutil.rmtree(run_dir)
        print(f"Deleted low-PE run dir: {run_dir}")
        return True
    except Exception:
        return False


def _iter_cph_time_weight(filepath):
    """Yield (time_ns, photon_weight) for every 'P' line of a .cph file.

    The trailing field of a .cph 'P' line is the number of photons the bunch
    represents (VGrisu writes i_bunch.photons there).  Files without that field
    fall back to weight 1.0, so the weight sum then degrades to the bunch count.
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == "*":
                continue
            if not line.startswith("P "):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                t = float(parts[6])          # ctime [ns]
            except ValueError:
                continue
            weight = 1.0
            if len(parts) > 10:              # trailing "photons" field
                try:
                    weight = float(parts[10])
                except ValueError:
                    weight = 1.0
            if not np.isfinite(weight) or weight < 0.0:
                weight = 1.0
            yield t, weight


"""def read_cph_photon_stats_old(filepath):
    # Return total photons and the max count per 10 ns bin from .cph data.
    try:
        if not os.path.exists(filepath):
            return 0, 0

        total_photons = 0
        min_time = None
        for t in _iter_cph_times(filepath):
            total_photons += 1
            if min_time is None or t < min_time:
                min_time = t

        if total_photons == 0 or min_time is None:
            return 0, 0

        counts = {}
        max_per_10ns = 0
        for t in _iter_cph_times(filepath):
            idx = int((t - min_time) // 10.0)
            counts[idx] = counts.get(idx, 0) + 1
            if counts[idx] > max_per_10ns:
                max_per_10ns = counts[idx]

        return total_photons, max_per_10ns
    except Exception:
        return 0, 0
"""

def read_cph_photon_stats(filepath):
    # Return the total photon count (weight-summed over bunches) and the max
    # photon count in any 10 ns bin, from .cph data.
    # Single-pass: collect times/weights once, then bin with numpy.
    try:
        if not os.path.exists(filepath):
            return 0, 0

        times = []
        weights = []
        for t, w in _iter_cph_time_weight(filepath):
            times.append(t)
            weights.append(w)
        if not times:
            return 0, 0

        times = np.asarray(times, dtype=np.float64)
        weights = np.asarray(weights, dtype=np.float64)

        total_photons = float(weights.sum())
        if total_photons <= 0.0:
            return 0, 0

        min_time = times.min()
        idx = ((times - min_time) // 10.0).astype(np.int64)

        # bincount requires non-negative ints; idx is >= 0 since min_time is the min.
        counts = np.bincount(idx, weights=weights)
        max_per_10ns = float(counts.max()) if counts.size else 0.0

        return int(round(total_photons)), int(round(max_per_10ns))
    except Exception:
        return 0, 0

def _parse_metadata_recenter(lines):
    x_val = None
    y_val = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("recentered_telescope_x_m:"):
            x_val = stripped.split(":", 1)[1].strip()
        elif stripped.startswith("recentered_telescope_y_m:"):
            y_val = stripped.split(":", 1)[1].strip()
        if x_val is not None and y_val is not None:
            break
    if x_val is None or y_val is None:
        return None
    return float(x_val), float(y_val)


def _parse_telescope_yaml_offset(lines):
    x_val = None
    y_val = None
    in_offset = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("offset_local_m:"):
            in_offset = True
            continue
        if in_offset:
            if stripped.startswith("x:"):
                x_val = stripped.split(":", 1)[1].strip()
            elif stripped.startswith("y:"):
                y_val = stripped.split(":", 1)[1].strip()
            elif stripped.startswith("z:"):
                continue
            elif ":" in stripped and not stripped.startswith("x:") and not stripped.startswith("y:"):
                if x_val is not None or y_val is not None:
                    break
        if x_val is not None and y_val is not None:
            break
    if x_val is None or y_val is None:
        return None
    return float(x_val), float(y_val)


def read_correction_offsets(paths):
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        try:
            if path.endswith(".json"):
                with open(path, "r", encoding="utf-8") as f:
                    report = json.load(f)
                x = float(report["center_x_m"])
                y = float(report["center_y_m"])
                return x, y, 1
            with open(path, "r", encoding="utf-8") as f:
                lines = f.readlines()
            parsed = _parse_metadata_recenter(lines)
            if parsed is not None:
                return parsed[0], parsed[1], 1
            parsed = _parse_telescope_yaml_offset(lines)
            if parsed is not None:
                return parsed[0], parsed[1], 1
        except Exception:
            continue
    return None, None, 0


def read_duration_seconds(path):
    if not path or not os.path.exists(path):
        return None, 0
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if not stripped.startswith("duration_seconds:"):
                    continue
                raw_val = stripped.split(":", 1)[1].strip()
                if not raw_val:
                    return None, 0
                try:
                    num = float(raw_val)
                except ValueError:
                    return None, 0
                if num.is_integer():
                    num = int(num)
                return num, 1
    except Exception:
        return None, 0
    return None, 0


_INJECTION_META_KEYS = (
    "injection_zenith_deg",
    "injection_azimuth_deg",
    "injection_x_m",
    "injection_y_m",
    "injection_height_m",
    "injection_z_minus_R_m",
    "survival_probability",
)


def read_injection_metadata(path):
    """Parse injection / survival fields from a run's metadata.yaml.

    Returns a dict keyed by ``_INJECTION_META_KEYS`` with float values (or None
    when a field is missing).
    """
    out = {k: None for k in _INJECTION_META_KEYS}
    if not path or not os.path.exists(path):
        return out
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or ":" not in stripped:
                    continue
                key, raw_val = stripped.split(":", 1)
                key = key.strip()
                if key not in out:
                    continue
                raw_val = raw_val.strip()
                if not raw_val:
                    continue
                try:
                    out[key] = float(raw_val)
                except ValueError:
                    continue
    except Exception:
        return {k: None for k in _INJECTION_META_KEYS}
    return out


def delete_file_if_exists(path):
    if not path or not os.path.exists(path):
        return False
    try:
        os.remove(path)
        print(f"Deleted working file: {path}")
        return True
    except Exception:
        return False


def _weighted_radius_containment(radius, weight, frac):
    order = np.argsort(radius)
    r_sorted = np.asarray(radius)[order]
    w_sorted = np.asarray(weight)[order]
    total = w_sorted.sum()
    if total <= 0:
        return 0.0
    cum = np.cumsum(w_sorted) / total
    idx = np.searchsorted(cum, frac)
    idx = min(idx, len(r_sorted) - 1)
    return float(r_sorted[idx])


"""def read_particle_metrics_old(filepath):
    try:
        if not filepath or not os.path.exists(filepath):
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)

        table = pq.read_table(filepath)
        if table.num_rows == 0:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)

        df = table.to_pandas()
        if df.empty:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)

        energy_col = "kinetic_energy" if "kinetic_energy" in df.columns else "energy"
        if energy_col not in df.columns:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)

        if "weight" in df.columns:
            weight = df["weight"].to_numpy(dtype=float)
        else:
            weight = np.ones(len(df), dtype=float)

        if "radius" in df.columns:
            radius = df["radius"].to_numpy(dtype=float)
        elif "x" in df.columns and "y" in df.columns:
            radius = np.sqrt(df["x"].to_numpy(dtype=float) ** 2 + df["y"].to_numpy(dtype=float) ** 2)
        else:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)

        pdg = df["pdg"].to_numpy()
        energy = df[energy_col].to_numpy(dtype=float)

        total_particles = float(np.sum(weight))
        if total_particles <= 0:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)

        abs_pdg = np.abs(pdg)
        muon_mask = abs_pdg == 13
        electron_mask = abs_pdg == 11
        gamma_mask = pdg == 22
        hadronic_mask = np.isin(pdg, [111, 211, -211, 321, -321, 130, 310, 2212, -2212, 2112, -2112])
        other_mask = ~(muon_mask | electron_mask | gamma_mask | hadronic_mask)

        max_muon_energy = float(np.max(energy[muon_mask])) if np.any(muon_mask) else 0.0
        r68_m = _weighted_radius_containment(radius, weight, 0.68)
        r99_m = _weighted_radius_containment(radius, weight, 0.99)

        muon_component = float(weight[muon_mask].sum() / total_particles)
        electron_component = float(weight[electron_mask].sum() / total_particles)
        hadronic_component = float(weight[hadronic_mask].sum() / total_particles)
        gamma_component = float(weight[gamma_mask].sum() / total_particles)
        other_component = float(weight[other_mask].sum() / total_particles)

        return (
            max_muon_energy,
            r68_m,
            r99_m,
            total_particles,
            muon_component,
            electron_component,
            hadronic_component,
            gamma_component,
            other_component,
        )
    except Exception:
        return (0, 0, 0, 0, 0, 0, 0, 0, 0)"""

def read_particle_metrics(filepath):
    """Return particle metrics plus the observation-level centre-of-mass state.

    Returns a 20-tuple:
        (max_muon_energy, r68_m, r99_m, total_particles,
         muon_component, electron_component, hadronic_component,
         gamma_component, other_component,
         arrival_zen_deg, arrival_az_deg, x_mag, y_mag, z_mag,
         max_e_muon_zen_deg, max_e_muon_az_deg,
         com_x_m, com_y_m, com_z_m, com_concentration)

    ``arrival_*`` is the energy-weighted centre-of-mass (COM) direction of ALL
    particles at the observation level. With ``w_i = kinetic_energy_i *
    weight_i`` (``weight`` is the CORSIKA thinning weight) and unit direction
    vectors ``n_i``::

        n_com = sum(w_i * n_i) / |sum(w_i * n_i)|
        zen   = degrees(arccos(-n_com_z))   az = degrees(atan2(n_com_y, n_com_x))

    ``com_concentration`` is ``|sum(w_i n_i)| / sum(w_i)`` (1.0 = perfectly
    collimated). ``com_x/y/z_m`` is the energy-weighted centroid of all
    particle positions.

    ``x_mag``/``y_mag``/``z_mag`` keep their original meaning: the final
    position of the highest-kinetic-energy muon. ``max_e_muon_zen/az``
    preserves the legacy direction (that same most energetic muon), i.e. what
    ``arrival_*`` used to hold before the COM switch.

    A value is None when the parquet has no rows or lacks the needed columns.
    If the COM vector vanishes (e.g. two exactly back-to-back particles),
    ``arrival_*`` falls back to the most energetic muon direction.
    """
    empty = (0, 0, 0, 0, 0, 0, 0, 0, 0,
             None, None, None, None, None,
             None, None, None, None, None, None)
    try:
        if not filepath or not os.path.exists(filepath):
            return empty

        # First peek at schema to decide which columns actually exist,
        # so we only request columns we know are present (avoids
        # pyarrow raising on missing columns and avoids reading unused ones).
        schema_names = set(pq.ParquetFile(filepath).schema_arrow.names)

        energy_col = "kinetic_energy" if "kinetic_energy" in schema_names else "energy"
        if energy_col not in schema_names:
            return empty

        wanted = {"pdg", energy_col}
        has_weight = "weight" in schema_names
        has_radius = "radius" in schema_names
        has_xy = "x" in schema_names and "y" in schema_names
        has_xyz = all(c in schema_names for c in ("x", "y", "z"))
        has_dir = all(c in schema_names for c in ("nx", "ny", "nz"))

        if has_weight:
            wanted.add("weight")
        if has_radius:
            wanted.add("radius")
        if has_xy:
            wanted.add("x")
            wanted.add("y")
        if not (has_radius or has_xy):
            return empty
        if has_xyz:
            wanted.add("z")
        if has_dir:
            wanted.update(("nx", "ny", "nz"))

        if "pdg" not in schema_names:
            return empty

        table = pq.read_table(filepath, columns=list(wanted))
        if table.num_rows == 0:
            return empty

        df = table.to_pandas()
        if df.empty:
            return empty

        if has_weight:
            weight = df["weight"].to_numpy(dtype=float)
        else:
            weight = np.ones(len(df), dtype=float)

        if has_radius:
            radius = df["radius"].to_numpy(dtype=float)
        else:
            radius = np.sqrt(df["x"].to_numpy(dtype=float) ** 2 + df["y"].to_numpy(dtype=float) ** 2)

        pdg = df["pdg"].to_numpy()
        energy = df[energy_col].to_numpy(dtype=float)

        total_particles = float(np.sum(weight))
        if total_particles <= 0:
            return empty

        abs_pdg = np.abs(pdg)
        muon_mask = abs_pdg == 13
        electron_mask = abs_pdg == 11
        gamma_mask = pdg == 22
        hadronic_mask = np.isin(pdg, [111, 211, -211, 321, -321, 130, 310, 2212, -2212, 2112, -2112])
        other_mask = ~(muon_mask | electron_mask | gamma_mask | hadronic_mask)

        max_muon_energy = float(np.max(energy[muon_mask])) if np.any(muon_mask) else 0.0
        r68_m = _weighted_radius_containment(radius, weight, 0.68)
        r99_m = _weighted_radius_containment(radius, weight, 0.99)

        muon_component = float(weight[muon_mask].sum() / total_particles)
        electron_component = float(weight[electron_mask].sum() / total_particles)
        hadronic_component = float(weight[hadronic_mask].sum() / total_particles)
        gamma_component = float(weight[gamma_mask].sum() / total_particles)
        other_component = float(weight[other_mask].sum() / total_particles)

        # ── Legacy: arrival state of the most energetic muon ──────────────
        # Kept verbatim; exported as max_e_muon_zen/az + x_mag/y_mag/z_mag.
        arrival_zen = arrival_az = None
        x_mag = y_mag = z_mag = None
        if np.any(muon_mask):
            idx_mu = np.flatnonzero(muon_mask)
            best = int(idx_mu[int(np.argmax(energy[muon_mask]))])
            if has_dir:
                nx = float(df["nx"].to_numpy(dtype=float)[best])
                ny = float(df["ny"].to_numpy(dtype=float)[best])
                nz = float(df["nz"].to_numpy(dtype=float)[best])
                nrm = float(np.sqrt(nx * nx + ny * ny + nz * nz))
                if nrm > 0.0:
                    nx, ny, nz = nx / nrm, ny / nrm, nz / nrm
                    arrival_zen = float(np.degrees(np.arccos(np.clip(-nz, -1.0, 1.0))))
                    arrival_az = float(np.degrees(np.arctan2(ny, nx)) % 360.0)
            if has_xyz:
                x_mag = float(df["x"].to_numpy(dtype=float)[best])
                y_mag = float(df["y"].to_numpy(dtype=float)[best])
                z_mag = float(df["z"].to_numpy(dtype=float)[best])
        max_e_muon_zen, max_e_muon_az = arrival_zen, arrival_az

        # ── Energy-weighted centre-of-mass direction (ALL particles) ──────
        # w_i = kinetic_energy_i * weight_i (thinning-corrected), n_com =
        # sum(w_i n_i) / |sum(w_i n_i)|. Overrides arrival_zen/arrival_az;
        # if the COM vector vanishes it falls back to the most energetic muon.
        com_x_m = com_y_m = com_z_m = None
        com_concentration = None
        if has_dir:
            nx_all = df["nx"].to_numpy(dtype=float)
            ny_all = df["ny"].to_numpy(dtype=float)
            nz_all = df["nz"].to_numpy(dtype=float)
            nrm_all = np.sqrt(nx_all * nx_all + ny_all * ny_all + nz_all * nz_all)
            w_all = energy * weight
            good = (nrm_all > 0.0) & (w_all > 0.0)
            if np.any(good):
                wg = w_all[good]
                nxu = nx_all[good] / nrm_all[good]
                nyu = ny_all[good] / nrm_all[good]
                nzu = nz_all[good] / nrm_all[good]
                sx = float(np.sum(wg * nxu))
                sy = float(np.sum(wg * nyu))
                sz = float(np.sum(wg * nzu))
                s_mag = float(np.sqrt(sx * sx + sy * sy + sz * sz))
                w_sum = float(np.sum(wg))
                if s_mag > 0.0 and w_sum > 0.0:
                    arrival_zen = float(np.degrees(np.arccos(np.clip(-sz / s_mag, -1.0, 1.0))))
                    arrival_az = float(np.degrees(np.arctan2(sy, sx)) % 360.0)
                    com_concentration = s_mag / w_sum
                    if has_xyz:
                        com_x_m = float(np.sum(wg * df["x"].to_numpy(dtype=float)[good]) / w_sum)
                        com_y_m = float(np.sum(wg * df["y"].to_numpy(dtype=float)[good]) / w_sum)
                        com_z_m = float(np.sum(wg * df["z"].to_numpy(dtype=float)[good]) / w_sum)

        return (
            max_muon_energy,
            r68_m,
            r99_m,
            total_particles,
            muon_component,
            electron_component,
            hadronic_component,
            gamma_component,
            other_component,
            arrival_zen,
            arrival_az,
            x_mag,
            y_mag,
            z_mag,
            max_e_muon_zen,
            max_e_muon_az,
            com_x_m,
            com_y_m,
            com_z_m,
            com_concentration,
        )
    except Exception:
        return empty

def read_metrics(filepath, pe_threshold=1.0):
    try:
        with uproot.open(filepath) as f:
            tree = f['Events/T0;1']
            available = set(tree.keys())
            fadc_branches = [f"vFADCTraces{i}" for i in range(256)
                             if f"vFADCTraces{i}" in available]

            if not fadc_branches:
                return (0.0, 0.0, 0.0, 0.0,
                    0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0)

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
    parser.add_argument("--injection-x", type=float, default=0.0,
                        help="Primary injection offset in x (m). Must match the run "
                             "tree level name h*_ix*_iy*.")
    parser.add_argument("--injection-y", type=float, default=0.0,
                        help="Primary injection offset in y (m). Must match the run "
                             "tree level name h*_ix*_iy*.")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--base-path", required=True)
    parser.add_argument("--correction-report-name", default="metadata.yaml",
                        help="Correction source filename stored beside CARE/CORSIKA8/GROPT/CIO directories")
    parser.add_argument("--Max_PE_cut", type=float, default=None,
                        help="Delete run dir when max_pe < this value (PE).")
    parser.add_argument("--dry-run-delete", action="store_true",
                        help="Log low-PE deletions without removing files.")
    args = parser.parse_args()

    print(f"DEBUG: base_path = {args.base_path}")
    print(f"DEBUG: pid={args.pid} E_mag={args.energy_str} R={args.radius} seed={args.seed} tel_y={args.tel_y}")

    with open(args.chunk_file) as f:
        combos = json.load(f)

    def combo_seed(combo):
        if len(combo) == 6:
            return int(float(combo[5]))
        if args.seed is None:
            raise SystemExit("ERROR: --seed is required when chunk entries have only 5 values")
        return int(args.seed)

    print(f"DEBUG: Loaded {len(combos)} combos from {args.chunk_file}")
    if combos:
        # Show first combo's constructed path
        c = combos[0]
        first_seed = combo_seed(c)
        sample_path = make_path(args.base_path, args.pid, args.energy_str,
                                float(c[0]), float(c[1]), float(c[2]),
                                float(c[3]), args.tel_y, float(c[4]),
                                args.radius, first_seed,
                                args.injection_x, args.injection_y)
        sample_report = make_correction_paths(
            args.base_path, args.pid, args.energy_str,
            float(c[0]), float(c[1]), float(c[2]),
            float(c[3]), args.tel_y, float(c[4]),
            args.radius, first_seed,
            args.correction_report_name,
            args.injection_x, args.injection_y,
        )
        print(f"DEBUG: First combo = {c}")
        print(f"DEBUG: First path  = {sample_path}")
        print(f"DEBUG: First correction sources= {sample_report}")
        sample_run_dir = make_run_dir(
            args.base_path, args.pid, args.energy_str,
            float(c[0]), float(c[1]), float(c[2]),
            float(c[3]), args.tel_y, float(c[4]),
            args.radius, first_seed,
            args.injection_x, args.injection_y,
        )
        #sample_c8 = make_corsika8_path(
        #    args.base_path, args.pid, args.energy_str,
        #    float(c[0]), float(c[1]), float(c[2]),
        #    float(c[3]), args.tel_y, float(c[4]),
        #    args.radius, args.seed,
        #)
        sample_cph = make_cph_path(
            args.base_path, args.pid, args.energy_str,
            float(c[0]), float(c[1]), float(c[2]),
            float(c[3]), args.tel_y, float(c[4]),
            args.radius, first_seed,
            args.injection_x, args.injection_y,
        )
        print(f"DEBUG: First run dir = {sample_run_dir}")
        #print(f"DEBUG: First C8 path = {sample_c8}")
        print(f"DEBUG: First CPH path = {sample_cph}")

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
            combo_seed_value = combo_seed(combo)
            zen, az, h, x_i, z_i = float(combo[0]), float(combo[1]), float(combo[2]), float(combo[3]), float(combo[4])

            path = make_path(args.base_path, args.pid, args.energy_str,
                             zen, az, h, x_i, args.tel_y, z_i,
                             args.radius, combo_seed_value,
                             args.injection_x, args.injection_y)
            c8_path = make_corsika8_path(
                args.base_path, args.pid, args.energy_str,
                zen, az, h, x_i, args.tel_y, z_i,
                args.radius, combo_seed_value,
                args.injection_x, args.injection_y,
            )
            correction_report_path = make_correction_paths(
                args.base_path, args.pid, args.energy_str,
                zen, az, h, x_i, args.tel_y, z_i,
                args.radius, combo_seed_value,
                args.correction_report_name,
                args.injection_x, args.injection_y,
            )
            run_dir = make_run_dir(
                args.base_path, args.pid, args.energy_str,
                zen, az, h, x_i, args.tel_y, z_i,
                args.radius, combo_seed_value,
                args.injection_x, args.injection_y,
            )
            remove_corsika8_tables(run_dir)
            cph_path = make_cph_path(
                args.base_path, args.pid, args.energy_str,
                zen, az, h, x_i, args.tel_y, z_i,
                args.radius, combo_seed_value,
                args.injection_x, args.injection_y,
            )
            particles_path = make_particles_path(
                args.base_path, args.pid, args.energy_str,
                zen, az, h, x_i, args.tel_y, z_i,
                args.radius, combo_seed_value,
                args.injection_x, args.injection_y,
            )

            (max_pe, time_at_max, avg_pe, total_pe,
             n_hit, image_size, frac_brightest, conc_2,
             pulse_width, rise_time, time_spread, time_gradient,
             peak_to_charge, baseline_rms, care_found) = read_metrics(path)
            (max_muon_energy_GeV, r68_m, r99_m, total_particles,
             muon_component, electron_component, hadronic_component,
             gamma_component, other_component,
             arrival_zen, arrival_az,
             x_mag, y_mag, z_mag,
             max_e_muon_zen, max_e_muon_az,
             com_x_m, com_y_m, com_z_m, com_concentration) = read_particle_metrics(particles_path)

            correction_x, correction_y, correction_found = read_correction_offsets(correction_report_path)
            file_size_mb = read_directory_size_mb(run_dir)
            cph_photon_count, cph_max_photons_10ns = read_cph_photon_stats(cph_path)
            duration_path = os.path.join(run_dir, "metadata.yaml")
            runtime_seconds, runtime_found = read_duration_seconds(duration_path)
            injection_meta = read_injection_metadata(duration_path)
            delete_file_if_exists(c8_path)

            low_pe_candidate = (
                care_found
                and args.Max_PE_cut is not None
                and max_pe < args.Max_PE_cut
            )
            deleted_low_pe = 0
            if low_pe_candidate:
                if args.dry_run_delete:
                    print(f"DRY-RUN: max_pe {max_pe:.2f} < {args.Max_PE_cut} for {run_dir}")
                deleted_low_pe = 1 if delete_run_dir(
                    run_dir,
                    base_path=args.base_path,
                    dry_run=args.dry_run_delete,
                ) else 0

            found = 1 if (care_found and correction_found) or deleted_low_pe else 0
            if found == 1:
                found_count += 1
            else:
                not_found_count += 1
                # Print first 5 incomplete combinations for debugging
                if not_found_count <= 5:
                    if not care_found:
                        print(f"DEBUG: CARE NOT FOUND [{not_found_count}]: {path}")
                    if not correction_found:
                        print(f"DEBUG: CORRECTION SOURCE NOT FOUND [{not_found_count}]: {correction_report_path}")

            correction_x_csv = round(correction_x, 6) if correction_found else ""
            correction_y_csv = round(correction_y, 6) if correction_found else ""
            runtime_csv = runtime_seconds if runtime_found else ""
            arrival_zen_csv = round(arrival_zen, 4) if arrival_zen is not None else ""
            arrival_az_csv = round(arrival_az, 4) if arrival_az is not None else ""
            max_e_muon_zen_csv = round(max_e_muon_zen, 4) if max_e_muon_zen is not None else ""
            max_e_muon_az_csv = round(max_e_muon_az, 4) if max_e_muon_az is not None else ""
            x_mag_csv = round(x_mag, 4) if x_mag is not None else ""
            y_mag_csv = round(y_mag, 4) if y_mag is not None else ""
            z_mag_csv = round(z_mag, 4) if z_mag is not None else ""
            com_x_csv = round(com_x_m, 4) if com_x_m is not None else ""
            com_y_csv = round(com_y_m, 4) if com_y_m is not None else ""
            com_z_csv = round(com_z_m, 4) if com_z_m is not None else ""
            com_concentration_csv = round(com_concentration, 6) if com_concentration is not None else ""

            def _inj(key, nd=6):
                val = injection_meta.get(key)
                return round(val, nd) if val is not None else ""

            writer.writerow({
                "pid":               args.pid,
                "energy_string":     args.energy_str,
                "radius":            args.radius,
                "seed":              combo_seed_value,
                "scan_zen":          zen,
                "scan_az":           az,
                "zen":               arrival_zen_csv,
                "az":                arrival_az_csv,
                "com_concentration": com_concentration_csv,
                "max_e_muon_zen":    max_e_muon_zen_csv,
                "max_e_muon_az":     max_e_muon_az_csv,
                "com_x_m":           com_x_csv,
                "com_y_m":           com_y_csv,
                "com_z_m":           com_z_csv,
                "height":            h,
                "tel_x":             x_i,
                "tel_y":             args.tel_y,
                "tel_z":             z_i,
                "tel_r":             round(np.sqrt(x_i**2 + z_i**2), 4),
                "x_mag":             x_mag_csv,
                "y_mag":             y_mag_csv,
                "z_mag":             z_mag_csv,
                "injection_zenith_deg":  _inj("injection_zenith_deg"),
                "injection_azimuth_deg": _inj("injection_azimuth_deg"),
                "injection_x":           _inj("injection_x_m"),
                "injection_y":           _inj("injection_y_m"),
                "injection_height_m":    _inj("injection_height_m"),
                "injection_z_minus_R_m": _inj("injection_z_minus_R_m"),
                "survival_probability":  _inj("survival_probability"),
                "file_size_MB":      round(file_size_mb, 6),
                "cph_photon_count": cph_photon_count,
                "cph_max_photons_10ns": cph_max_photons_10ns,
                "correction_x_m":    correction_x_csv,
                "correction_y_m":    correction_y_csv,
                "runtime_seconds":   runtime_csv,
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
                "max_muon_energy_GeV": round(max_muon_energy_GeV, 6),
                "r68_m":               round(r68_m, 4),
                "r99_m":               round(r99_m, 4),
                "total_particles":     round(total_particles, 4),
                "muon_component":      round(muon_component, 6),
                "electron_component":  round(electron_component, 6),
                "hadronic_component":  round(hadronic_component, 6),
                "gamma_component":     round(gamma_component, 6),
                "other_component":     round(other_component, 6),
                "deleted_low_pe":    deleted_low_pe,
                "file_found":        found,
            })


            


    print(f"Wrote {len(combos)} rows to {args.output}")
    print(f"DEBUG SUMMARY: found={found_count}, not_found={not_found_count}")

if __name__ == "__main__":
    main()
