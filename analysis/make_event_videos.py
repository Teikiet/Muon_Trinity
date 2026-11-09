#!/usr/bin/env python3
"""
make_event_videos.py

Render time-synchronised movies of the Cherenkov-light evolution of one event:

    1. CORSIKA8/cherenkov_hits*.dat   (eventio, ground-plane photon density)
    2. CIO/cherenkov_hits.cph         (GrOptics input photon list)
    3. GROPT/cherenkov_hits.root      (tree T1, focal-plane photon density)
    4. CARE/cherenkov_hits.root       (tree Events/T0, 16x16 camera PE map)

All movies share a single time axis (see --grid / --time-bin-ns).  By default
the window is the GROPT photon span at 1 ns per frame; the CARE FADC samples
(10 ns each) are held across the 1 ns frames so every movie stays in step.

By default the script runs the CORSIKA8 -> CIO -> GrOptics -> CARE chain for the
requested event and then renders.  Alternatively, point it at an existing run
folder with --run-dir (a folder holding CORSIKA8/, CIO/, GROPT/, CARE/,
corsika8_output/) and it renders straight from those outputs.

The movies are written into an output folder named after the event id,
    pid{PDG}_E{ENERGY}_zen{ZENITH}_az{AZIMUTH}_h{HEIGHT}
(e.g. pid13_E1e4_zen90.0_az270.0_h10000).  Choose where that folder lives with
--out-root.

Examples
--------
    # run the chain (if needed) and render the movies
    python make_event_videos.py \
        --pdg 13 --energy 1e3 --zenith 90 --azimuth 267 --height 4000

    # only render from existing outputs (never submit a job)
    python make_event_videos.py \
        --pdg 13 --energy 1e3 --zenith 90 --azimuth 267 --height 4000 --no-run

    # render an existing run folder into
    #   <out-root>/pid13_E1e4_zen90.0_az270.0_h10000/
    python make_event_videos.py \
        --run-dir /path/to/.../x0.0_y0.0_z0.0 --out-root /path/to/animations
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm, LinearSegmentedColormap

# ── Reuse the project's readers (plot_GROPT / plot_CARE / read_cherenkov_hits) ─
ANALYSIS_DIR = "/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/analysis"
if ANALYSIS_DIR not in sys.path:
    sys.path.append(ANALYSIS_DIR)
from ReadTrinity import (read_cherenkov_hits, plot_GROPT, plot_CARE,  # noqa: E402
                         read_cph)

# ── Magnetic-field muon back-tracker (back_track_muon) lives here ─────────────
MAG_DIR = ("/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/"
           "cluster_corsika8_2.0")

# ── Detector / electronics constants (ReadTrinity + save_CARE2csv) ────────────
DC_TO_PE       = 24.1     # DAC -> photoelectron
BASELINE_DC    = 500.0    # pedestal
SAMPLES_PER_NS = 10.0     # 10 ns per FADC sample
TRIGGER_SAMPLE = 17       # FADC sample index that corresponds to the trigger

DEFAULT_BASE_DIR = "/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/muon_sim_chain_tree_visual"#"/scratch/general/vast/u1520754/muon_sim_chain_tree_visual"
DEFAULT_SLURM_SCRIPT = (
    "/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/cluster_corsika8_2.0/"
    "base_muon_event_generation/run_corsika8_trinity_chain_tree.slurm"
)

# CARE colour scheme (identical to the reference plotting snippet)
FUZZY_CMAP = LinearSegmentedColormap.from_list("fuzzy", [
    (0.00, "#273ca9"),
    (0.20, "#0e58c3"),
    (0.30, "#1481ca"),
    (0.45, "#07a3c0"),
    (0.60, "#3bbb98"),
    (0.75, "#9bba77"),
    (0.88, "#f3c344"),
    (1.00, "#ffff1e"),
])

if shutil.which("ffmpeg"):
    matplotlib.rcParams["animation.ffmpeg_path"] = shutil.which("ffmpeg")


# =============================================================================
# Small formatting helpers (mirror the tree-naming conventions)
# =============================================================================
def _fmt_angle(v):
    """90 -> '90.0'  (the tree uses zen{ZENITH} / az{AZIMUTH})."""
    return f"{float(v):.1f}"


def _fmt_intlike(v):
    """4000.0 -> '4000',  1.0 -> '1'  (used for h and telescope-radius levels)."""
    fv = float(v)
    return str(int(round(fv))) if abs(fv - round(fv)) < 1e-9 else f"{fv:g}"


def _fmt_offset(v):
    """0 -> '0.0',  1.5 -> '1.5'  (used for x/y/z telescope offsets)."""
    fv = float(v)
    return "0.0" if abs(fv) < 1e-12 else f"{fv:.1f}"


# =============================================================================
# Path resolution
# =============================================================================
def _numeric_match(parent, prefix, value, tol=1e-6):
    """Return parent/<prefix><num> whose <num> equals value (numeric match)."""
    if not os.path.isdir(parent):
        return None
    for name in sorted(os.listdir(parent)):
        if not name.startswith(prefix):
            continue
        rest = name[len(prefix):]
        try:
            if abs(float(rest) - float(value)) < tol:
                return os.path.join(parent, name)
        except ValueError:
            continue
    return None


def _position_match(parent, tx, ty, tz, tol=1e-6):
    """Match a 'x<a>_y<b>_z<c>' leaf directory numerically."""
    if not os.path.isdir(parent):
        return None
    pat = re.compile(r"^x([-0-9.eE+]+)_y([-0-9.eE+]+)_z([-0-9.eE+]+)$")
    for name in sorted(os.listdir(parent)):
        m = pat.match(name)
        if not m:
            continue
        try:
            a, b, c = float(m.group(1)), float(m.group(2)), float(m.group(3))
        except ValueError:
            continue
        if abs(a - tx) < tol and abs(b - ty) < tol and abs(c - tz) < tol:
            return os.path.join(parent, name)
    return None


def resolve_run_dir(base_dir, pdg, energy, tel_r, seed, zen, az, height,
                    tx, ty, tz, verbose=True):
    """Resolve the run directory from the output-tree layout."""
    root = os.path.join(base_dir, f"Muon_pid{pdg}_E{energy}_R{tel_r}")
    pdg_dir = os.path.join(root, f"pdg{pdg}_E{energy}_r{tel_r}_s{seed}")

    if not os.path.isdir(pdg_dir):
        if verbose:
            print(f"[paths] pdg directory not found: {pdg_dir}")
        return None

    zen_dir = _numeric_match(pdg_dir, "zen", zen)
    if zen_dir is None:
        return None
    az_dir = _numeric_match(zen_dir, "az", az)
    if az_dir is None:
        return None
    h_dir = _numeric_match(az_dir, "h", height)
    if h_dir is None:
        return None
    return _position_match(h_dir, tx, ty, tz)


def parse_event_from_path(run_dir):
    """Best-effort extraction of pdg/energy/zenith/azimuth/height from a run path.

    Recognises the tree components ``pdg{PDG}_E{ENERGY}_r...``, ``zen{...}``,
    ``az{...}`` and ``h{...}``.  Returns a dict; missing entries are ``None``.
    """
    out = {"pdg": None, "energy": None, "zenith": None,
           "azimuth": None, "height": None}
    for part in os.path.normpath(run_dir).split(os.sep):
        m = re.match(r"^pdg(-?\d+)_E(.+?)_r", part)
        if m:
            out["pdg"] = int(m.group(1))
            out["energy"] = m.group(2)
            continue
        m = re.match(r"^zen([-+0-9.eE]+)$", part)
        if m:
            out["zenith"] = float(m.group(1))
            continue
        m = re.match(r"^az([-+0-9.eE]+)$", part)
        if m:
            out["azimuth"] = float(m.group(1))
            continue
        m = re.match(r"^h([-+0-9.eE]+)$", part)
        if m:
            out["height"] = float(m.group(1))
            continue
    return out


def event_id_str(pdg, energy, zenith, azimuth, height):
    """Canonical event id used to name the animation output folder."""
    return (f"pid{_fmt_intlike(pdg)}_E{energy}"
            f"_zen{_fmt_angle(zenith)}_az{_fmt_angle(azimuth)}"
            f"_h{_fmt_intlike(height)}")


def find_file(*candidates):
    """Return the first existing path from the candidates."""
    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None


def find_corsika8_dat(run_dir):
    c8 = os.path.join(run_dir, "CORSIKA8")
    base = find_file(
        os.path.join(c8, "cherenkov_hits_base.dat"),
        os.path.join(c8, "cherenkov_hits.dat"),
    )
    if base:
        return base
    for root, _dirs, files in os.walk(run_dir):
        for f in sorted(files):
            if f.endswith(".dat") and f.startswith("cherenkov_hits"):
                return os.path.join(root, f)
    return None


def find_gropt_root(run_dir):
    g = os.path.join(run_dir, "GROPT")
    p = find_file(os.path.join(g, "cherenkov_hits.root"))
    if p:
        return p
    if os.path.isdir(g):
        for name in sorted(os.listdir(g)):
            sub = os.path.join(g, name)
            if os.path.isdir(sub):
                p = find_file(os.path.join(sub, "cherenkov_hits.root"))
                if p:
                    return p
        for name in sorted(os.listdir(g)):
            if name.endswith(".root"):
                return os.path.join(g, name)
    return None


def find_care_root(run_dir):
    c = os.path.join(run_dir, "CARE")
    p = find_file(os.path.join(c, "cherenkov_hits.root"))
    if p:
        return p
    if os.path.isdir(c):
        for name in sorted(os.listdir(c)):
            sub = os.path.join(c, name)
            if os.path.isdir(sub):
                p = find_file(os.path.join(sub, "cherenkov_hits.root"))
                if p:
                    return p
        for name in sorted(os.listdir(c)):
            if name.endswith(".root"):
                return os.path.join(c, name)
    return None


def find_cph_file(run_dir):
    """Locate the CIO cph photon list (GrOptics input)."""
    c = os.path.join(run_dir, "CIO")
    p = find_file(os.path.join(c, "cherenkov_hits.cph"))
    if p:
        return p
    if os.path.isdir(c):
        for name in sorted(os.listdir(c)):
            sub = os.path.join(c, name)
            if os.path.isdir(sub):
                p = find_file(os.path.join(sub, "cherenkov_hits.cph"))
                if p:
                    return p
        for name in sorted(os.listdir(c)):
            if name.endswith(".cph"):
                return os.path.join(c, name)
    return None


# =============================================================================
# Optional: run the simulation chain
# =============================================================================
def outputs_ready(run_dir):
    if not run_dir or not os.path.isdir(run_dir):
        return False
    return (find_corsika8_dat(run_dir) is not None
            and find_gropt_root(run_dir) is not None
            and find_care_root(run_dir) is not None)


def compute_injection(args):
    """Predict the muon injection (position + direction) and survival via the
    magnetic-field back-tracker ``magnetic_field_muon_correction.back_track_muon``.

    Mirrors ``adaptive_controller.compute_cell_injection``.  Returns a dict with
    the ``--injection-*`` fields, or ``None`` when the correction cannot be
    applied (non-muon primary, unparsable energy, no valid crossing, ...).
    """
    pdg = int(args.pdg)
    if pdg not in (13, -13):
        print(f"[mag] pdg={pdg} is not a muon; magnetic correction not applicable.")
        return None
    charge_sign = -1.0 if pdg == 13 else 1.0
    try:
        energy = float(args.energy)
    except (TypeError, ValueError):
        print(f"[mag] cannot parse energy '{args.energy}' as a number; "
              f"skipping magnetic correction.")
        return None

    if MAG_DIR not in sys.path:
        sys.path.append(MAG_DIR)
    try:
        from magnetic_field_muon_correction import back_track_muon
    except Exception as exc:                                    # noqa: BLE001
        print(f"[mag] could not import back_track_muon ({exc}); "
              f"skipping magnetic correction.")
        return None

    try:
        res = back_track_muon(energy, charge_sign, float(args.zenith),
                              float(args.azimuth), float(args.height),
                              H_obs_m=float(args.obs_level))
    except Exception as exc:                                    # noqa: BLE001
        print(f"[mag] back_track_muon failed ({exc}); "
              f"skipping magnetic correction.")
        return None

    if not bool(np.asarray(res["success"]).item()):
        print("[mag] back-tracking found no valid crossing; "
              "skipping magnetic correction.")
        return None
    surv = float(np.asarray(res["survival_probability"]).item())
    if not np.isfinite(surv):
        print("[mag] survival_probability is not finite; "
              "skipping magnetic correction.")
        return None

    inj = {
        "injection_zenith_deg":  float(np.asarray(res["injection_zenith_deg"]).item()),
        "injection_azimuth_deg": float(np.asarray(res["injection_azimuth_deg"]).item()),
        "injection_x_m":         float(np.asarray(res["injection_x_m"]).item()),
        "injection_y_m":         float(np.asarray(res["injection_y_m"]).item()),
        "injection_height_m":    float(np.asarray(res["injection_height_m"]).item()),
        "injection_z_minus_R_m": float(np.asarray(res["injection_z_minus_R_m"]).item()),
        "survival_probability":  surv,
    }
    print(f"[mag] injection zenith={inj['injection_zenith_deg']:.5f} deg  "
          f"azimuth={inj['injection_azimuth_deg']:.5f} deg  "
          f"x={inj['injection_x_m']:.3f} m  y={inj['injection_y_m']:.3f} m  "
          f"survival={inj['survival_probability']:.4f}")
    return inj


def submit_chain(args, output_base_dir, log_file, injection=None):
    """Submit run_corsika8_trinity_chain_tree.slurm and return the job id.

    When ``injection`` is given (from the magnetic-field back-tracker) the six
    ``--injection-*`` / ``--survival-probability`` optionals are appended, exactly
    as ``adaptive_controller`` does.  ``hadron_model`` and ``rerun_event`` are
    always passed as the two positional optionals so a rerun can be forced.
    """
    cmd = [
        "sbatch",
        f"--job-name=corsika8_trinity_pid{args.pdg}_E{args.energy}_s{args.seed}",
        args.slurm_script,
        output_base_dir,
        log_file,
        str(args.pdg),
        str(args.energy),
        _fmt_angle(args.zenith),
        _fmt_angle(args.azimuth),
        _fmt_intlike(args.height),
        _fmt_intlike(args.obs_level),
        _fmt_offset(args.tel_x),
        _fmt_offset(args.tel_y),
        _fmt_offset(args.tel_z),
        _fmt_intlike(args.cherenkov_radius),
        _fmt_intlike(args.tel_radius),
        str(int(args.seed)),
        str(args.hadron_model),
        str(args.rerun_event),
        f"--update_time_stamp={str(bool(args.update_time_stamp)).lower()}",
    ]
    if injection:
        cmd.append(f"--injection-zenith={injection['injection_zenith_deg']:.9g}")
        cmd.append(f"--injection-azimuth={injection['injection_azimuth_deg']:.9g}")
        cmd.append(f"--injection-x={injection['injection_x_m']:.9g}")
        cmd.append(f"--injection-y={injection['injection_y_m']:.9g}")
        cmd.append(f"--injection-z-minus-r={injection['injection_z_minus_R_m']:.9g}")
        cmd.append(f"--survival-probability={injection['survival_probability']:.9g}")
    print("[run] submitting:")
    print("      " + " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.stdout.strip():
        print("[run] " + res.stdout.strip())
    if res.returncode != 0:
        print("[run] " + res.stderr.strip())
        raise RuntimeError("sbatch submission failed")
    jobid = res.stdout.strip().split()[-1]
    print(f"[run] job id = {jobid}")
    return jobid


def wait_for_job(jobid, timeout):
    t_start = time.time()
    while True:
        try:
            r = subprocess.run(["squeue", "-h", "-j", str(jobid)],
                               capture_output=True, text=True)
        except FileNotFoundError:
            print("[run] squeue not available; cannot poll job state.")
            return
        if not r.stdout.strip():
            print(f"[run] job {jobid} left the queue.")
            return
        if timeout and (time.time() - t_start) > timeout:
            raise TimeoutError(f"job {jobid} still running after {timeout}s")
        time.sleep(20)


def slurm_log_paths(jobid, slurm_script):
    """Resolve the SLURM --output/--error paths from the script's #SBATCH lines."""
    out_pat = err_pat = None
    try:
        with open(slurm_script, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                s = line.strip()
                if not s.startswith("#SBATCH"):
                    continue
                if "--output=" in s:
                    out_pat = s.split("--output=", 1)[1].strip()
                elif "--error=" in s:
                    err_pat = s.split("--error=", 1)[1].strip()
    except OSError:
        pass

    def _expand(pat):
        if not pat:
            return None
        return pat.replace("%A", str(jobid)).replace("%j", str(jobid))

    return _expand(out_pat), _expand(err_pat)


def job_state(jobid):
    """Return (state, exit_code) for a SLURM job via sacct; (None, None) on failure."""
    try:
        res = subprocess.run(
            ["sacct", "-j", str(jobid), "--noheader", "--parsable2",
             "--format=State,ExitCode"],
            capture_output=True, text=True)
    except FileNotFoundError:
        return None, None
    for line in res.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("|")
        if "." in parts[0]:          # skip sub-steps (jobid.batch / jobid.extern)
            continue
        toks = parts[0].split()
        state = toks[0] if toks else parts[0]
        exit_code = parts[1] if len(parts) > 1 else ""
        return state, exit_code
    return None, None


def report_job_failure(jobid, slurm_script, tail_lines=30):
    """Print the SLURM job state, the log paths and a tail of the stderr log."""
    state, exit_code = job_state(jobid)
    out_path, err_path = slurm_log_paths(jobid, slurm_script)

    print("-" * 72)
    print(f"[run] job {jobid} did NOT produce the expected outputs.")
    if state:
        print(f"[run] final state = {state}   exit code = {exit_code}")
    if out_path:
        print(f"[run] stdout log  : {out_path}")
    if err_path:
        print(f"[run] stderr log  : {err_path}")
    if err_path and os.path.isfile(err_path):
        try:
            with open(err_path, "r", encoding="utf-8", errors="replace") as fh:
                lines = fh.readlines()
            if lines:
                n = min(tail_lines, len(lines))
                print(f"[run] ---- last {n} lines of {err_path} ----")
                for ln in lines[-n:]:
                    print("      " + ln.rstrip("\n"))
        except OSError:
            pass
    print("-" * 72)


# =============================================================================
# Data readers
# =============================================================================
def read_corsika8(dat_path):
    """Return X, Y [m], T [ns] and bunch weight from an eventio Cherenkov file."""
    X, Y, T, _cx, _cy, _wl, _zem, _th, _ph, _sl, W = \
        read_cherenkov_hits(dat_path, plot=False)
    X = np.asarray(X, dtype=float) / 100.0   # cm -> m
    Y = np.asarray(Y, dtype=float) / 100.0
    T = np.asarray(T, dtype=float)
    W = np.asarray(W, dtype=float)
    return X, Y, T, W


def read_gropt(root_path):
    """Return photonX, photonY [mm] and T [ns] from the GrOptics T1 tree."""
    _cosx, _cosy, px, py, T = plot_GROPT(root_path, plot=False)
    return (np.asarray(px, dtype=float),
            np.asarray(py, dtype=float),
            np.asarray(T, dtype=float))


def read_care(root_path):
    """Return per-sample 16x16 PE grids and the total PE per sample."""
    (_cam, _tn, _tot, _mx, arrays, pixel_config,
     min_x, _maxx, min_y, _maxy, pitch) = plot_CARE(root_path, plot=False)

    traces = {}
    for pix_id in range(256):
        name = f"vFADCTraces{pix_id}"
        if name in arrays:
            traces[pix_id] = np.clip(
                (arrays[name][0] - BASELINE_DC) / DC_TO_PE, 0.0, None)

    if not traces:
        raise RuntimeError(f"No FADC traces found in {root_path}")

    n_samples = min(len(t) for t in traces.values())
    grids = np.zeros((n_samples, 16, 16), dtype=float)
    for pix_id, tr in traces.items():
        col = int((pixel_config[pix_id][1] - min_x) / pitch)
        row = int((pixel_config[pix_id][2] - min_y) / pitch)
        if 0 <= row < 16 and 0 <= col < 16:
            grids[:n_samples, row, col] = tr[:n_samples]

    total_pe_per_sample = grids.sum(axis=(1, 2))
    return grids, total_pe_per_sample


def read_cph_photons(cph_path):
    """Return the CIO .cph photon ground positions x, y, time [ns] and weight.

    .cph columns: 0 xGrd, 1 yGrd, 2 xcos, 3 ycos, 4 hgtEmiss, 5 time,
                  6 wavelength, 7 type, 8 tel, 9 photons (bunch weight).
    Older files without the trailing weight fall back to weight 1.0, so the
    weighted histogram then degrades to a bunch count.
    """
    arr = read_cph(cph_path)
    x = np.asarray(arr[:, 0], dtype=float)
    y = np.asarray(arr[:, 1], dtype=float)
    t = np.asarray(arr[:, 5], dtype=float)
    if arr.shape[1] > 9:
        w = np.array(arr[:, 9], dtype=float)          # copy, so we can sanitise
        w[~np.isfinite(w) | (w < 0.0)] = 1.0
    else:
        w = np.ones(x.size, dtype=float)
    return x, y, t, w


# =============================================================================
# Frame builders
# =============================================================================
def _safe_range(v, pad_frac=0.05):
    v = np.asarray(v, dtype=float)
    if v.size == 0:
        return (-1.0, 1.0)
    lo, hi = float(np.min(v)), float(np.max(v))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        c = lo if np.isfinite(lo) else 0.0
        return (c - 1.0, c + 1.0)
    pad = pad_frac * (hi - lo)
    return (lo - pad, hi + pad)


def build_hist_frames(x, y, t, t0, dt, n_frames, bins, xrange, yrange,
                      weights=None):
    """Return a (n_frames, bins, bins) stack of 2D histograms (axis0 = time)."""
    H = np.zeros((n_frames, bins, bins), dtype=float)
    if x.size == 0:
        return H
    idx = np.floor((t - t0) / dt).astype(np.int64)
    valid = (idx >= 0) & (idx < n_frames)
    idx, x, y = idx[valid], x[valid], y[valid]
    w = weights[valid] if weights is not None else None
    for i in range(n_frames):
        m = idx == i
        if not m.any():
            continue
        hh, _xe, _ye = np.histogram2d(
            x[m], y[m], bins=bins, range=[xrange, yrange],
            weights=(w[m] if w is not None else None))
        H[i] = hh.T          # (rows=y, cols=x) for imshow
    return H


def care_grid_plot(grid):
    """Apply the same orientation transform as the reference CARE snippet."""
    return np.flip(np.flip(grid.T, axis=0), axis=1)


# =============================================================================
# Movie writers
# =============================================================================
def _write_movie(fig, update, n_frames, out_path, fps):
    writer = animation.FFMpegWriter(
        fps=fps, codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-loglevel", "error"])
    with writer.saving(fig, out_path, dpi=110):
        for i in range(n_frames):
            update(i)
            writer.grab_frame()
    plt.close(fig)


def write_hist_movie(H, extent, t0, dt, n_frames, out_path, fps,
                     title, xlabel, ylabel, cb_label, cmap="viridis"):
    """Render a time-stack of 2D histograms to an mp4."""
    vmax = float(H.max()) if H.size else 1.0
    norm = LogNorm(vmin=1.0, vmax=max(vmax, 1.001))
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    im = ax.imshow(np.ma.masked_where(H[0] <= 0, H[0]), origin="lower",
                   extent=extent, cmap=cmap, norm=norm,
                   interpolation="nearest", aspect="auto")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label(cb_label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    def update(i):
        im.set_data(np.ma.masked_where(H[i] <= 0, H[i]))
        ax.set_title(f"{title}\nt = {t0 + i * dt:.0f} ns")
        return (im,)

    _write_movie(fig, update, n_frames, out_path, fps)
    return out_path


def write_care_movie(grids, sample_for_frame, t0, dt, n_frames, out_path, fps,
                     title, n_pix=16):
    """Render the CARE 16x16 camera PE map over time using the fuzzy colour map."""
    vmax = max(float(grids.max()) if grids.size else 1.0, 1e-6)
    edges = np.arange(n_pix + 1)
    centers = np.arange(n_pix) + 0.5

    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    k0 = int(sample_for_frame[0])
    g0 = grids[k0] if k0 >= 0 else np.zeros((n_pix, n_pix))
    mesh = ax.pcolormesh(edges, edges, care_grid_plot(g0), cmap=FUZZY_CMAP,
                         shading="flat", vmin=0.0, vmax=vmax)
    cb = fig.colorbar(mesh, ax=ax, fraction=0.06, pad=0.01)
    cb.set_label("PE", fontsize=15)
    cb.ax.tick_params(labelsize=12)
    ax.set_xticks(centers)
    ax.set_xticklabels([str(i) for i in range(n_pix)], fontsize=12)
    ax.set_yticks(centers)
    ax.set_yticklabels([str(i) for i in range(n_pix)], fontsize=12)
    ax.set_xlabel("Pixel along camera y axis", fontsize=15)
    ax.set_ylabel("Pixel along camera x axis", fontsize=15)
    ax.set_aspect("equal")

    def update(i):
        k = int(sample_for_frame[i])
        g = grids[k] if k >= 0 else np.zeros((n_pix, n_pix))
        mesh.set_array(care_grid_plot(g).ravel())
        ax.set_title(f"{title} (t = {t0 + i * dt:.0f} ns)", fontsize=14)
        return (mesh,)

    _write_movie(fig, update, n_frames, out_path, fps)
    return out_path


def write_combined_movie(panels, grids, sample_for_frame, t_origin, dt, n_frames,
                         out_path, fps, n_pix=16):
    """Side-by-side movie of the 2D-histogram sources plus the CARE camera.

    ``panels`` is a list of dicts with keys ``H``, ``extent``, ``title``,
    ``xlabel``, ``ylabel``, ``cb`` (each ``H`` is a (n_frames, bins, bins)
    stack). The CARE camera is always appended as the last panel.
    """
    n = len(panels) + 1
    fig, axes = plt.subplots(1, n, figsize=(5.4 * n, 5.8))
    axes = np.atleast_1d(axes)

    ims = []
    for ax, p in zip(axes, panels):
        vmax = max(float(p["H"].max()) if p["H"].size else 1.0, 1.001)
        im = ax.imshow(np.ma.masked_where(p["H"][0] <= 0, p["H"][0]),
                       origin="lower", extent=p["extent"], cmap="viridis",
                       norm=LogNorm(vmin=1.0, vmax=vmax),
                       interpolation="nearest", aspect="auto")
        fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02).set_label(p["cb"])
        ax.set_xlabel(p["xlabel"])
        ax.set_ylabel(p["ylabel"])
        ax.set_title(p["title"])
        ims.append(im)

    axc = axes[-1]
    vmax_p = max(float(grids.max()) if grids.size else 1.0, 1e-6)
    edges = np.arange(n_pix + 1)
    centers = np.arange(n_pix) + 0.5
    k0 = int(sample_for_frame[0])
    g0 = grids[k0] if k0 >= 0 else np.zeros((n_pix, n_pix))
    mesh = axc.pcolormesh(edges, edges, care_grid_plot(g0), cmap=FUZZY_CMAP,
                          shading="flat", vmin=0.0, vmax=vmax_p)
    fig.colorbar(mesh, ax=axc, fraction=0.05, pad=0.02).set_label("PE")
    axc.set_xticks(centers)
    axc.set_yticks(centers)
    axc.set_xlabel("Pixel along camera y axis")
    axc.set_ylabel("Pixel along camera x axis")
    axc.set_aspect("equal")
    axc.set_title("CARE (camera)")

    def update(i):
        for im, p in zip(ims, panels):
            im.set_data(np.ma.masked_where(p["H"][i] <= 0, p["H"][i]))
        k = int(sample_for_frame[i])
        g = grids[k] if k >= 0 else np.zeros((n_pix, n_pix))
        mesh.set_array(care_grid_plot(g).ravel())
        fig.suptitle(f"t = {t_origin + i * dt:.0f} ns", fontsize=16)
        return tuple(ims) + (mesh,)

    _write_movie(fig, update, n_frames, out_path, fps)
    return out_path


# =============================================================================
# CLI / main
# =============================================================================
def build_arg_parser():
    p = argparse.ArgumentParser(
        description=("Run the CORSIKA8 -> GrOptics -> CARE chain for one event "
                     "and render time-synchronised Cherenkov-light movies."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--pdg", type=int, default=None, help="primary PDG code")
    p.add_argument("--energy", type=str, default=None,
                   help="energy string used in the tree, e.g. 1e3")
    p.add_argument("--zenith", type=float, default=None,
                   help="expected arrival zenith [deg]")
    p.add_argument("--azimuth", type=float, default=None,
                   help="expected arrival azimuth [deg]")
    p.add_argument("--height", type=float, default=None,
                   help="production height [m]")
    p.add_argument("--run-dir", dest="run_dir", default=None,
                   help="existing run folder (holding CORSIKA8/, CIO/, GROPT/, "
                        "CARE/, corsika8_output/) to render without running the "
                        "chain")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--tel-radius", dest="tel_radius", type=float, default=1.0)
    p.add_argument("--tel-x", dest="tel_x", type=float, default=0.0)
    p.add_argument("--tel-y", dest="tel_y", type=float, default=0.0)
    p.add_argument("--tel-z", dest="tel_z", type=float, default=0.0)
    p.add_argument("--obs-level", dest="obs_level", type=float, default=2944.0)
    p.add_argument("--cherenkov-radius", dest="cherenkov_radius",
                   type=float, default=15.0)
    p.add_argument("--base-dir", dest="base_dir", default=DEFAULT_BASE_DIR)
    p.add_argument("--slurm-script", dest="slurm_script",
                   default=DEFAULT_SLURM_SCRIPT)
    p.add_argument("--no-run", dest="no_run", action="store_true",
                   help="never submit a job; only render existing outputs")
    p.add_argument("--run-timeout", dest="run_timeout", type=int, default=1800,
                   help="max seconds to wait for the submitted job")
    p.add_argument("--time-bin-ns", dest="time_bin_ns", type=float, default=1.0,
                   help="time resolution (frame width) [ns]")
    p.add_argument("--fps", type=float, default=10.0, help="frames per second")
    p.add_argument("--bins", type=int, default=100,
                   help="2D histogram bins for CORSIKA8 / CIO-cph / GrOptics")
    p.add_argument("--max-frames", dest="max_frames", type=int, default=0,
                   help="cap the number of frames (0 = no cap)")
    p.add_argument("--gropt-time-offset", dest="gropt_time_offset",
                   type=float, default=None,
                   help="ns added to GROPT times before binning (default: auto)")
    p.add_argument("--care-time-offset", dest="care_time_offset",
                   type=float, default=None,
                   help="absolute time [ns] of the CARE FADC trigger sample")
    p.add_argument("--cio-time-offset", dest="cio_time_offset",
                   type=float, default=None,
                   help="ns added to CIO .cph times before binning (default: auto)")
    p.add_argument("--grid", dest="grid",
                   choices=["gropt", "cph", "care", "corsika"], default="gropt",
                   help="time window: 'gropt' (default) = first->last GROPT "
                        "photon, 'cph' = cph span, 'care' = CARE FADC readout, "
                        "'corsika' = CORSIKA span")
    p.add_argument("--out-root", dest="out_root", default=None,
                   help="parent folder for the animation folder (default: the "
                        "run directory)")
    p.add_argument("--out-dir", dest="out_dir", default=None,
                   help="exact output folder (overrides --out-root and the "
                        "event-id sub-folder)")
    p.add_argument("--combined", action="store_true",
                   help="also write a side-by-side CORSIKA8|GrOptics|CARE movie")
    p.add_argument("--magnetic-correction", dest="magnetic_correction",
                   action="store_true", default=True,
                   help="run the magnetic-field muon back-tracker and pass "
                        "--injection-* to the slurm (muons only)")
    p.add_argument("--no-magnetic-correction", dest="magnetic_correction",
                   action="store_false",
                   help="disable the magnetic-field correction")
    p.add_argument("--min-survival-probability", dest="min_survival_probability",
                   type=float, default=0.0,
                   help="abort if the back-tracked survival probability is below this")
    p.add_argument("--hadron-model", dest="hadron_model", default="SIBYLL-2.3d",
                   help="hadronic interaction model passed to the slurm")
    p.add_argument("--rerun-event", dest="rerun_event", default="false",
                   choices=["false", "detector_sim", "all"],
                   help="force a rerun of an already-existing event")
    p.add_argument("--update-time-stamp", dest="update_time_stamp",
                   action="store_true", default=False,
                   help="pass --update_time_stamp=true to the slurm")
    return p


def main():
    args = build_arg_parser().parse_args()

    tel_r = _fmt_intlike(args.tel_radius)
    seed = int(args.seed)
    dt = float(args.time_bin_ns)

    # ── 0. event identity (from --run-dir path and/or the CLI) ───────────────
    run_dir = None
    if args.run_dir:
        run_dir = os.path.abspath(args.run_dir)
        if not os.path.isdir(run_dir):
            raise SystemExit(f"[error] --run-dir not found: {run_dir}")
        derived = parse_event_from_path(run_dir)
        for key, val in derived.items():
            if getattr(args, key) is None and val is not None:
                setattr(args, key, val)

    ev_fields = (args.pdg, args.energy, args.zenith, args.azimuth, args.height)
    ev_id = (event_id_str(*ev_fields)
             if all(v is not None for v in ev_fields) else None)

    print("=" * 72)
    print("make_event_videos")
    print(f"  pdg={args.pdg}  energy={args.energy}  "
          f"zen={args.zenith}  az={args.azimuth}  height={args.height}")
    print(f"  tel_r={tel_r}  seed={seed}  dt={dt:g} ns")
    print(f"  magnetic correction={args.magnetic_correction}  "
          f"hadron_model={args.hadron_model}  rerun={args.rerun_event}")
    print("=" * 72)

    # ── 1. locate the run directory (and run the chain if needed) ─────────────
    if run_dir is not None:
        print(f"[paths] run directory : {run_dir}  (--run-dir; chain not run)")
    else:
        if not all(v is not None for v in ev_fields):
            raise SystemExit(
                "[error] --pdg/--energy/--zenith/--azimuth/--height are required "
                "unless --run-dir is given.")

        run_dir = resolve_run_dir(args.base_dir, args.pdg, args.energy, tel_r,
                                  seed, args.zenith, args.azimuth, args.height,
                                  args.tel_x, args.tel_y, args.tel_z)

        if not outputs_ready(run_dir):
            if args.no_run:
                # render whatever is available, but require at least the master file
                if run_dir is None or find_corsika8_dat(run_dir) is None:
                    raise SystemExit(
                        f"[error] no usable outputs in {run_dir!r} and "
                        f"--no-run was set.")
                print("[run] --no-run: rendering the outputs that are available.")
            else:
                output_base_dir = os.path.join(
                    args.base_dir,
                    f"Muon_pid{args.pdg}_E{args.energy}_R{tel_r}")
                log_file = os.path.join(
                    output_base_dir, "logs",
                    f"s{seed}_zen{_fmt_angle(args.zenith)}"
                    f"_az{_fmt_angle(args.azimuth)}"
                    f"_h{_fmt_intlike(args.height)}_x{_fmt_offset(args.tel_x)}"
                    f"_z{_fmt_offset(args.tel_z)}_ix0.0_iy0.0.log")
                os.makedirs(os.path.dirname(log_file), exist_ok=True)

                injection = None
                if args.magnetic_correction:
                    injection = compute_injection(args)
                    if (injection is not None
                            and args.min_survival_probability > 0.0
                            and injection["survival_probability"]
                            < args.min_survival_probability):
                        raise SystemExit(
                            f"[mag] survival probability "
                            f"{injection['survival_probability']:.4f} < "
                            f"--min-survival-probability "
                            f"{args.min_survival_probability}; aborting.")
                else:
                    print("[mag] magnetic correction disabled "
                          "(--no-magnetic-correction).")

                jobid = submit_chain(args, output_base_dir, log_file, injection)
                wait_for_job(jobid, args.run_timeout)

                run_dir = resolve_run_dir(args.base_dir, args.pdg, args.energy,
                                          tel_r, seed, args.zenith,
                                          args.azimuth, args.height, args.tel_x,
                                          args.tel_y, args.tel_z)
                if not outputs_ready(run_dir):
                    report_job_failure(jobid, args.slurm_script)
                    raise SystemExit(
                        f"[error] outputs still missing after job {jobid} "
                        f"({run_dir!r}). See the logs above.")

        print(f"[paths] run directory : {run_dir}")

    # ── 2. output folder (named after the event id) ───────────────────────────
    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_root = args.out_root or run_dir
        if ev_id is None:
            ev_id = os.path.basename(os.path.normpath(run_dir))
            print(f"[paths] event id incomplete; naming the output folder "
                  f"'{ev_id}'.")
        out_dir = os.path.join(out_root, ev_id)
    os.makedirs(out_dir, exist_ok=True)
    print(f"[paths] output folder: {out_dir}")

    c8_path = find_corsika8_dat(run_dir)
    cph_path = find_cph_file(run_dir)
    gropt_path = find_gropt_root(run_dir)
    care_path = find_care_root(run_dir)
    print(f"[paths] CORSIKA8     : {c8_path}")
    print(f"[paths] CIO cph      : {cph_path}")
    print(f"[paths] GROPT        : {gropt_path}")
    print(f"[paths] CARE         : {care_path}")

    if c8_path is None:
        raise SystemExit(
            "[error] CORSIKA8 cherenkov .dat file not found (master clock).")

    # ── 2. read the three data sources ───────────────────────────────────────
    X, Y, T, W = read_corsika8(c8_path)
    print(f"[data] CORSIKA8 bunches: {X.size}  "
          f"t=[{T.min():.1f}, {T.max():.1f}] ns")

    have_gropt = gropt_path is not None
    gx = gy = gt = None
    if have_gropt:
        gx, gy, gt = read_gropt(gropt_path)
        print(f"[data] GROPT photons   : {gx.size}  "
              f"t=[{gt.min():.1f}, {gt.max():.1f}] ns")

    have_cph = cph_path is not None
    cx = cy = ct = cw = None
    if have_cph:
        cx, cy, ct, cw = read_cph_photons(cph_path)
        print(f"[data] CIO cph photons: {cx.size} bunches "
              f"({cw.sum():.0f} photons)  t=[{ct.min():.1f}, {ct.max():.1f}] ns")

    have_care = care_path is not None
    grids = care_total = None
    if have_care:
        grids, care_total = read_care(care_path)
        print(f"[data] CARE samples    : {grids.shape[0]}  "
              f"max PE = {grids.max():.3f}")

    # ── 3. shared time grid ──────────────────────────────────────────────────
    t_c0 = float(T.min())          # CORSIKA first light (master clock)

    def _first_light_offset(t_min):
        """Shift a non-CORSIKA clock onto the CORSIKA clock (first light)."""
        if (t_min - t_c0) < 0.0 or (t_min - t_c0) > (float(T.max()) - t_c0):
            return t_c0 - t_min
        return 0.0

    g_off = 0.0
    if have_gropt:
        if args.gropt_time_offset is not None:
            g_off = float(args.gropt_time_offset)
        else:
            g_off = _first_light_offset(float(gt.min()))
            if g_off:
                print(f"[time] GROPT range does not overlap CORSIKA -> "
                      f"auto first-light shift of {g_off:.1f} ns")

    c_off = 0.0
    if have_cph:
        if args.cio_time_offset is not None:
            c_off = float(args.cio_time_offset)
        else:
            c_off = _first_light_offset(float(ct.min()))
            if c_off:
                print(f"[time] CPH range does not overlap CORSIKA -> "
                      f"auto first-light shift of {c_off:.1f} ns")

    # CARE: absolute time of the FADC trigger (auto-detected by peak alignment)
    trigger_abs = None
    if have_care:
        k_peak = int(np.argmax(care_total)) if care_total.size else 0
        if args.care_time_offset is not None:
            trigger_abs = float(args.care_time_offset)
        else:
            idx = np.floor((T - t_c0) / dt).astype(np.int64)
            idx = idx[idx >= 0]
            counts = np.bincount(idx) if idx.size else np.array([1], dtype=np.int64)
            peak_frame = int(np.argmax(counts))
            trigger_abs = t_c0 + (peak_frame - k_peak + TRIGGER_SAMPLE) * dt
        print(f"[time] CARE peak sample={k_peak}; "
              f"trigger at {trigger_abs:.1f} ns (CORSIKA clock)")

    # pick the window source (fall back when the requested one is unavailable)
    grid = args.grid
    if grid == "gropt" and not have_gropt:
        grid = "cph" if have_cph else "corsika"
        print(f"[time] no GROPT output; using the {grid} window instead.")
    if grid == "cph" and not have_cph:
        grid = "corsika"
        print("[time] no CIO cph output; using the CORSIKA window instead.")
    if grid == "care" and not have_care:
        grid = "corsika"
        print("[time] no CARE output; using the CORSIKA window instead.")

    care_dt = SAMPLES_PER_NS          # CARE FADC sample width (10 ns)

    if grid == "care":
        # frame i == CARE FADC sample i (10 ns frames); display clock = CORSIKA
        frame_dt = care_dt
        t_origin = trigger_abs - TRIGGER_SAMPLE * frame_dt
        n_frames = int(grids.shape[0])
        sample_for_frame = np.arange(n_frames, dtype=np.int64)
        sh = {"c8": 0.0, "cph": c_off, "gropt": g_off, "care": 0.0}
    else:
        # display clock == GROPT clock for gropt/cph windows, else CORSIKA
        frame_dt = dt
        if grid == "gropt":
            t_origin = float(gt.min())
            t_last = float(gt.max())
            sh = {"c8": -g_off, "cph": c_off - g_off, "gropt": 0.0, "care": -g_off}
        elif grid == "cph":
            t_origin = float(ct.min()) + (c_off - g_off)
            t_last = float(ct.max()) + (c_off - g_off)
            sh = {"c8": -g_off, "cph": c_off - g_off, "gropt": 0.0, "care": -g_off}
        else:  # corsika
            t_origin = t_c0
            t_last = float(T.max())
            sh = {"c8": 0.0, "cph": c_off, "gropt": g_off, "care": 0.0}

        n_frames = max(int(math.ceil((t_last - t_origin) / frame_dt)) + 1, 1)

        # CARE samples are held across their full 10 ns width
        sample_for_frame = np.full(n_frames, -1, dtype=np.int64)
        if have_care:
            trigger_disp = trigger_abs + sh["care"]
            for i in range(n_frames):
                k = int(math.floor((t_origin + i * frame_dt - trigger_disp)
                                   / care_dt)) + TRIGGER_SAMPLE
                if 0 <= k < grids.shape[0]:
                    sample_for_frame[i] = k

    if args.max_frames and n_frames > args.max_frames:
        print(f"[time] capping frames {n_frames} -> {args.max_frames}")
        n_frames = args.max_frames
        sample_for_frame = sample_for_frame[:n_frames]

    print(f"[time] grid = {grid.upper()} window: {n_frames} frames of "
          f"{frame_dt:g} ns ({n_frames * frame_dt:.0f} ns); "
          f"origin = {t_origin:.1f} ns ({grid} clock)")
    if n_frames <= 0:
        raise SystemExit("[error] no frames to render.")

    # ── 4. build the frame stacks ────────────────────────────────────────────
    xr = _safe_range(X)
    yr = _safe_range(Y)
    ext_c = (xr[0], xr[1], yr[0], yr[1])
    H_c = build_hist_frames(X, Y, T + sh["c8"], t_origin, frame_dt, n_frames,
                            args.bins, xr, yr, weights=W)

    if have_cph:
        cxr = _safe_range(cx)
        cyr = _safe_range(cy)
        ext_p = (cxr[0], cxr[1], cyr[0], cyr[1])
        H_p = build_hist_frames(cx, cy, ct + sh["cph"], t_origin, frame_dt,
                                n_frames, args.bins, cxr, cyr, weights=cw)
    else:
        ext_p = (-1.0, 1.0, -1.0, 1.0)
        H_p = np.zeros((n_frames, args.bins, args.bins), dtype=float)

    if have_gropt:
        gxr = _safe_range(gx)
        gyr = _safe_range(gy)
        ext_g = (gxr[0], gxr[1], gyr[0], gyr[1])
        H_g = build_hist_frames(gx, gy, gt + sh["gropt"], t_origin, frame_dt,
                                n_frames, args.bins, gxr, gyr)
    else:
        ext_g = (-1.0, 1.0, -1.0, 1.0)
        H_g = np.zeros((n_frames, args.bins, args.bins), dtype=float)

    # ── 5. render the movies ─────────────────────────────────────────────────
    written = []

    c8_out = os.path.join(out_dir, "cherenkov_corsika8.mp4")
    write_hist_movie(H_c, ext_c, t_origin, frame_dt, n_frames, c8_out, args.fps,
                     "CORSIKA 8 Cherenkov light", "x [m]", "y [m]", "photons")
    written.append(c8_out)
    print(f"[out] {c8_out}")

    if have_cph:
        cph_out = os.path.join(out_dir, "cherenkov_cph.mp4")
        write_hist_movie(H_p, ext_p, t_origin, frame_dt, n_frames, cph_out,
                         args.fps, "CIO cph photons", "cph x_grnd",
                         "cph y_grnd", "photons")
        written.append(cph_out)
        print(f"[out] {cph_out}")

    if have_gropt:
        gr_out = os.path.join(out_dir, "cherenkov_gropt.mp4")
        write_hist_movie(H_g, ext_g, t_origin, frame_dt, n_frames, gr_out,
                         args.fps, "GrOptics focal-plane photons",
                         "photon X [mm]", "photon Y [mm]", "photons")
        written.append(gr_out)
        print(f"[out] {gr_out}")

    if have_care:
        ca_out = os.path.join(out_dir, "cherenkov_care.mp4")
        write_care_movie(grids, sample_for_frame, t_origin, frame_dt, n_frames,
                         ca_out, args.fps, "CARE Camera Response")
        written.append(ca_out)
        print(f"[out] {ca_out}")

    if args.combined and have_gropt and have_care:
        panels = [dict(H=H_c, extent=ext_c, title="CORSIKA 8 (ground)",
                       xlabel="x [m]", ylabel="y [m]", cb="photons")]
        if have_cph:
            panels.append(dict(H=H_p, extent=ext_p, title="CIO cph (aperture)",
                               xlabel="cph x_grnd", ylabel="cph y_grnd",
                               cb="photons"))
        panels.append(dict(H=H_g, extent=ext_g, title="GrOptics (focal plane)",
                           xlabel="photon X [mm]", ylabel="photon Y [mm]",
                           cb="photons"))
        cb_out = os.path.join(out_dir, "cherenkov_combined.mp4")
        write_combined_movie(panels, grids, sample_for_frame, t_origin,
                             frame_dt, n_frames, cb_out, args.fps)
        written.append(cb_out)
        print(f"[out] {cb_out}")

    print("=" * 72)
    print("Done. Movies written:")
    for w in written:
        print(f"  {w}")
    print("=" * 72)


if __name__ == "__main__":
    main()
