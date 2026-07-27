#Extract data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os
import sys
import math

# ── Configuration ─────────────────────────────────────────────────────────────
PE_THRESHOLD = 20
#N_BOOTSTRAP = 1000
#CI_LEVEL = 0.99
E_MIN_GEV = 1e1
E_MAX_GEV = 1e6
SEEDS = [1,2,3,4]
PID = 13
RADIUS = 5
BASE_DIR = "/scratch/general/vast/u1520754/muon_sim_chain_tree"


# ── Helper functions ─────────────────────────────────────────────────────────
def load_scan(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    n_missing = (df["file_found"] == 0).sum()
    if n_missing > 0:
        print(f"  WARNING: {n_missing} rows had file_found=0 and will be dropped.")
        index = (df["file_found"] == 0)
        """print("zenith",set(df["zen"][index]))
        print("azimuth",set(df["az"][index]))
        print("height",set(df["height"][index]))"""
    df = df[df["file_found"] == 1].copy()
    df["r"] = np.sqrt(df["tel_x"]**2 + df["tel_z"]**2).round(4)
    df["detected"] = (df["max_pe"] >= PE_THRESHOLD).astype(float)
    return df

def slant_depth(zem, obs_height, zenith_rad):
    theta = np.pi - zenith_rad
    arg = (obs_height / zem) * np.sin(theta)
    arg = np.clip(arg, -1.0, 1.0)
    beta = np.where(np.isfinite(arg), np.arcsin(arg), 0.0)
    beta = np.clip(beta, 0.0, np.pi/2)
    #beta = np.clip(np.arcsin((obs_height / zem) * np.sin(theta)), 0, np.pi/2)
    alpha = np.pi - theta - beta
    s = zem * np.sin(alpha) / np.sin(theta)
    return np.clip(s, 0, np.inf)
def incidence_angle(Zenith_deg, Azimuth_deg):
    Zenith = Zenith_deg
    Azimuth = Azimuth_deg
    cos_Zenith = np.cos(np.radians(Zenith))
    cos_Azimuth = np.cos(np.radians(Azimuth - 180))
    cos_Incidence = np.sqrt(1 - cos_Zenith**2 - cos_Azimuth**2)
    Incidence = np.degrees(np.arccos(cos_Incidence))
    return Incidence


import re
import os

def get_available_energies(base_dir, pid, radius, e_min_GeV, e_max_GeV):
    """
    Scan BASE_DIR for folders matching the given PID and RADIUS,
    extract energy values, and filter by [e_min_GeV, e_max_GeV].
    Returns sorted list of (energy_GeV, energy_str) tuples.
    """
    pattern = re.compile(
        rf"Muon_pid{pid}_E([0-9]+\.?[0-9]*e[+-]?[0-9]+)_R{radius}$"
    )
    found = {}
    for name in os.listdir(base_dir):
        m = pattern.match(name)
        if m:
            energy_str = m.group(1)
            energy_val = float(energy_str)
            if e_min_GeV <= energy_val <= e_max_GeV:
                found[energy_val] = energy_str

    return sorted(found.items())  # sorted by energy_GeV


available = get_available_energies(BASE_DIR, PID, RADIUS, E_MIN_GEV, E_MAX_GEV)
print(f"Found {len(available)} energy files between {E_MIN_GEV:.0e} and {E_MAX_GEV:.0e} GeV:")
for e_val, e_str in available:
    print(f"  E={e_str} ({e_val:.3e} GeV)")

# ── Load only those files ────────────────────────────────────────────────────
dfs = []
missing_energies = []
for energy_val, energy_str in available:
    for seed in SEEDS:
        csv_path = (f"{BASE_DIR}/Muon_pid{PID}_E{energy_str}_R{RADIUS}/"
                    f"csv_output/scan_care_pid{PID}_E{energy_str}_R{RADIUS}_y0_s{seed}.csv")
        try:
            df_tmp = load_scan(csv_path)
            df_tmp["seed"] = seed
            df_tmp["energy_string"] = energy_str
            df_tmp["energy_GeV"] = energy_val
            dfs.append(df_tmp)
            print(f"  Loaded E={energy_str} seed={seed}: {len(df_tmp)} rows")
        except Exception as e:
            missing_energies.append(energy_str)
            print(f"  ⚠ E={energy_str} seed={seed}: {e}")


df = pd.concat(dfs, ignore_index=True)
print(f"\nMerged: {len(df)} total rows from {len(dfs)} file(s)")
print(f"Energy range: {df['energy_GeV'].min():.1f} – {df['energy_GeV'].max():.1f} GeV")
if missing_energies:
    print(f"Missing energies: {sorted(set(missing_energies))}")

# ── Derived columns ──────────────────────────────────────────────────────────
R = np.array(df["r"]) # m
X = np.array(df["tel_x"]) # m
Z = np.array(df["tel_z"]) # m
H = np.array(df["height"]) # m
Zenith = np.array(df["zen"]) #Degree
Azimuth = np.array(df["az"]) #Degree
Max_PE = np.array(df["max_pe"]) #Photo_electron
Time_at_Max_PE = np.array(df["time_at_max_pe_ns"]) #ns
Avg_PE = np.array(df["avg_pe"]) 
Total_PE = np.array(df["total_pe"])
Energy_GeV = np.array(df["energy_GeV"]) #GeV
Slant = slant_depth(H + 6371e3, 2944 + 6371e3, np.radians(Zenith))/1000 #km
cos_Zenith = np.cos(np.radians(Zenith))
cos_Azimuth = np.cos(np.radians(Azimuth - 180))
cos_Incidence = np.sqrt(1 - cos_Zenith**2 - cos_Azimuth**2)
Incidence = incidence_angle(Zenith, Azimuth) #Degree
#pulse_width = np.array(df["pulse_width_ns"])
#rise_time = np.array(df["rise_time_ns"])
x_mag = np.array(df["correction_x_m"]) #m
y_mag = np.array(df["correction_y_m"]) #m
file_size = np.array(df["file_size_MB"]) #MB
photon_count = np.array(df["cph_photon_count"])
#max_photon_10ns = np.array(df["cph_max_photons_10ns"])
E_list = np.sort(np.array(list(set(Energy_GeV))))
deleted_low_pe = np.array(df["deleted_low_pe"])
runtime_hours = np.array(df["runtime_seconds"])/3600 
File = []
Photon = []
PE = []
Time = []

Total_FIle_Size = np.nansum(file_size)/1000 - np.nansum(file_size[deleted_low_pe==1])/1000#MB to GB

# --- Coverage levels for nested bands ---
coverages = [50, 80, 95, 100]    # percent
cmap = plt.cm.Blues

# Storage for multi-coverage bands: dict[coverage] -> (lows, highs)
file_bands   = {c: ([], []) for c in coverages}
photon_bands = {c: ([], []) for c in coverages}
pe_bands     = {c: ([], []) for c in coverages}
time_bands     = {c: ([], []) for c in coverages}
for energy in E_list:
    mask = Energy_GeV == energy

    f  = np.nanmean(file_size[mask])
    p  = np.nanmean(photon_count[mask])
    pe = np.nanmean(Max_PE[mask])
    t = np.nanmean(runtime_hours[mask]) 

    File.append(f)
    Photon.append(p)
    PE.append(pe)
    Time.append(t)
    
    # Multi-coverage percentile bands
    for cov in coverages:
        lo_p = (100 - cov) / 2
        hi_p = 100 - lo_p
        file_bands[cov][0].append(np.nanpercentile(file_size[mask], lo_p))
        file_bands[cov][1].append(np.nanpercentile(file_size[mask], hi_p))
        photon_bands[cov][0].append(np.nanpercentile(photon_count[mask], lo_p))
        photon_bands[cov][1].append(np.nanpercentile(photon_count[mask], hi_p))
        pe_bands[cov][0].append(np.nanpercentile(Max_PE[mask], lo_p))
        pe_bands[cov][1].append(np.nanpercentile(Max_PE[mask], hi_p))
        time_bands[cov][0].append(np.nanpercentile(runtime_hours[mask], lo_p))
        time_bands[cov][1].append(np.nanpercentile(runtime_hours[mask], hi_p))
File   = np.array(File)
Photon = np.array(Photon)
PE     = np.array(PE)
Time   = np.array(Time)

for cov in coverages:
    file_bands[cov]   = (np.array(file_bands[cov][0]),   np.array(file_bands[cov][1]))
    photon_bands[cov] = (np.array(photon_bands[cov][0]), np.array(photon_bands[cov][1]))
    pe_bands[cov]     = (np.array(pe_bands[cov][0]),     np.array(pe_bands[cov][1]))
    time_bands[cov]   = (np.array(time_bands[cov][0]),     np.array(time_bands[cov][1]))

def plot_with_bands(E_list, mean, bands_dict, coverages, ylabel, title):
    fig, ax = plt.subplots()
    # Widest first (lightest), narrowest last (darkest)
    order = sorted(coverages, reverse=True)
    n = len(order)
    for i, cov in enumerate(order):
        lo, hi = bands_dict[cov]
        color = cmap(0.3 + 0.6 * (i / max(1, n - 1)))
        ax.fill_between(E_list, lo, hi, color=color, alpha=0.7,
                        label=f'{cov}% range')
    ax.plot(E_list, mean, 'k.-', lw=1.2, label='Mean')
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Energy (GeV)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    plt.show()


# --- Photon Count Plot ---
plot_with_bands(E_list, Photon, photon_bands, coverages,
                "Photon Count", "Average Photon Count vs Energy")

# --- PE Count Plot ---
plot_with_bands(E_list, PE, pe_bands, coverages,
                "Max PE Count per Pixel", "Average Max PE vs Energy")

plot_with_bands(E_list, File, file_bands, coverages,
                "File Size (MB)", "Average File Size vs Energy")

plot_with_bands(E_list, Time, time_bands, coverages,
                "Run time (hours))", "Average Run Time vs Energy")

print(Total_FIle_Size/1000, "TB", np.nanmean(Photon)/np.nanmean(File), "photon per MB")


# Empirical trigger efficiency: P(Max_PE >= threshold | photon_count, incidence)
# =============================================================================
from scipy.ndimage import gaussian_filter
from matplotlib.colors import SymLogNorm

# ---------- 1. Clean inputs ----------
m = (np.isfinite(photon_count) & (photon_count > 0)
     & np.isfinite(Max_PE)     & (Max_PE >= 0) & (R==0)
     & np.isfinite(Incidence))
pc_  = photon_count[m]
pe_  = Max_PE[m]
inc_ = Incidence[m]
hit_ = (pe_ >= PE_THRESHOLD).astype(float)

# ---------- 2. Define 2-D bins ----------
pc_edges  = 10**np.quantile(np.log10(pc_), np.linspace(0,1, 50))#np.logspace(np.log10(pc_.min()), np.log10(pc_.max()), 50)
inc_edges = np.quantile(inc_, np.linspace(0,1, 50)) #np.linspace(inc_.min(),inc_.max(), 50)
pc_cen    = np.sqrt(pc_edges[:-1] * pc_edges[1:])
inc_cen   = 0.5*(inc_edges[:-1] + inc_edges[1:])

# ---------- 3. Build empirical efficiency table ----------
N_tot, _, _ = np.histogram2d(pc_, inc_, bins=[pc_edges, inc_edges])
N_hit, _, _ = np.histogram2d(pc_, inc_, bins=[pc_edges, inc_edges], weights=hit_)

with np.errstate(invalid='ignore', divide='ignore'):
    eff = np.where(N_tot >= 5, N_hit / N_tot, np.nan)
    eff_err = np.where(N_tot >= 5, np.sqrt(eff*(1-eff)/N_tot), np.nan)

# Light smoothing to fill noisy bins (optional)
eff_smooth = eff.copy()
mask = np.isfinite(eff_smooth)
eff_smooth[~mask] = 0
weights = mask.astype(float)
sigma_value = 0.8
eff_smooth = gaussian_filter(eff_smooth, sigma=sigma_value) / np.clip(
    gaussian_filter(weights, sigma=sigma_value), 1e-6, None)

# ---------- 4. Lookup function with bilinear interpolation ----------
from scipy.interpolate import RegularGridInterpolator
interp = RegularGridInterpolator(
    (np.log10(pc_cen), inc_cen), eff_smooth,
    bounds_error=False, fill_value=None, method='linear')

def trigger_eff(N_gamma, inc_deg):
    N_gamma = np.asarray(N_gamma, dtype=float)
    inc_deg = np.asarray(inc_deg, dtype=float)
    N_b, inc_b = np.broadcast_arrays(N_gamma, inc_deg)
    pts = np.column_stack([
        np.log10(np.clip(N_b.ravel(), pc_.min(), pc_.max())),
        np.clip(inc_b.ravel(), inc_.min(), inc_.max()),
    ])
    vals = np.clip(interp(pts), 0.0, 1.0)
    return vals.reshape(N_b.shape)

# ---------- 5. Plot the efficiency map ----------
fig, ax = plt.subplots(figsize=(7,5))
pcm = ax.pcolormesh(pc_edges, inc_edges, eff_smooth.T, #eff_smooth.T,
                    cmap='viridis', norm=SymLogNorm(linthresh=np.min(eff_smooth[np.where(eff_smooth>0)]), linscale=0.5, base=10))
ax.set_xscale('log')
ax.set_xlabel("photon_count"); ax.set_ylabel("Incidence (deg)")
ax.set_title(f"Empirical trigger efficiency (Max_PE >= {PE_THRESHOLD})")
plt.colorbar(pcm, ax=ax, label='efficiency')
plt.tight_layout(); plt.show()

# ---------- 6. Slice validation ----------
bounds = np.linspace(inc_.min(), inc_.max(), 8)
slices = list(zip(bounds[:-1], bounds[1:]))
fig, axes = plt.subplots(1, len(slices), figsize=(4*len(slices), 4), sharey=True)
for ax, (ilo, ihi) in zip(axes, slices):
    sel = (inc_ >= ilo) & (inc_ < ihi)
    emp, err = [], []
    for i in range(len(pc_cen)):
        s2 = sel & (pc_ >= pc_edges[i]) & (pc_ < pc_edges[i+1])
        n = s2.sum()
        if n >= 5:
            p = np.mean(hit_[s2]); emp.append(p)
            err.append(np.sqrt(p*(1-p)/n) + 1.0/n)
        else:
            emp.append(np.nan); err.append(np.nan)
    emp, err = np.array(emp), np.array(err)
    pred = trigger_eff(pc_cen, 0.5*(ilo+ihi))
    ax.errorbar(pc_cen, emp, yerr=err, fmt='ko', ms=4, label='empirical')
    ax.plot(pc_cen, pred, 'r-', label='lookup')
    ax.set_xscale('log')#; ax.set_ylim(-0.05, 1.1)
    #ax.set_yscale('log')
    ax.set_xlabel("photon_count")
    ax.set_title(f"inc {ilo:.1f}-{ihi:.1f}°")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
axes[0].set_ylabel(f"Trigger eff. (Max_PE >= {PE_THRESHOLD})")
plt.tight_layout(); plt.show()


# Fitting Efficiency vs Photon Count Model (per incidence bin)
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.optimize import curve_fit
from scipy.stats import beta

# --- User settings ---
#PE_THRESHOLD = 5
n_photon_bins = 50
n_inc_bins = 20
ncols = 5

MIN_EVENTS_PER_INC_BIN = 1
MIN_EVENTS_PER_PC_BIN = 1
MIN_TRIGGERS_PER_INC = 1   # minimum total triggers in an incidence slice to attempt a fit
MAX_CHI2_DOF = 10.0
MAX_X0_ERR_LOGN = 1.0
MAX_W_ERR_LOGN = 1.0
EFF_ERR_FLOOR = 0.01
B_ERR_FLOOR = 0.01
W_ERR_FLOOR = 0.01
SATURATION_EFF = 1

B_MODEL  = "sigmoid_fixed"   # options: "const", "linear", "sigmoid"
X0_MODEL = "log_sigmoid"  # options: "linear", "quadratic", "log_sigmoid"
W_MODEL  = "sigmoid"   # options: "const", "linear", "powerlaw", "lognormal", "sigmoid"

N_BIN_VARIATIONS = 20
N_var_per_key_min = 1
NPC_RANGE  = (20, 200)
NINC_RANGE = (5, 100)

N_BOOTSTRAP = 100
CI_LEVEL = 68  # percent 68
CI_ALPHA = (100.0 - CI_LEVEL) / 2.0

_RNG = np.random.default_rng(0)

# ---- preprocess ----
mask = (
    np.isfinite(photon_count) & (photon_count > 0) &
    np.isfinite(Max_PE) & (Max_PE >= 0) & (R==0) &
    np.isfinite(Incidence)
)

pc = photon_count[mask]
inc = Incidence[mask]
pe = Max_PE[mask]
trig = (pe >= PE_THRESHOLD).astype(float)

pc_for_edges = pc[trig>0]
if pc_for_edges.size == 0:
    raise ValueError("No triggered events to define photon-count bin edges.")

#pc_edges = np.logspace(np.log10(pc_for_edges.min()), np.log10(pc_for_edges.max()), n_photon_bins + 1)
pc_edges = 10**(np.quantile(np.log10(pc_for_edges), np.linspace(0, 1, n_photon_bins + 1)))
n_photon_bins = len(pc_edges) - 1
#inc_edges = np.linspace(inc.min(), inc.max(), n_inc_bins + 1)
inc_edges = np.quantile(inc, np.linspace(0, 1, n_inc_bins + 1))

def clopper_pearson(hits, total, alpha=0.32):
    """
    Two-sided Clopper-Pearson binomial CI.
    alpha = 1 - CI_LEVEL/100  (e.g. 0.32 for 68% CI)
    Returns (lo, hi) arrays, same shape as hits/total.
    Handles hits=0 and hits=total correctly (non-zero width).
    """
    hits = np.asarray(hits, dtype=float)
    total = np.asarray(total, dtype=float)
    lo = np.where(hits > 0, beta.ppf(alpha / 2.0, hits, total - hits + 1), 0.0)
    hi = np.where(hits < total, beta.ppf(1.0 - alpha / 2.0, hits + 1, total - hits), 1.0)
    return lo, hi

from scipy.optimize import minimize

def fit_binomial_logistic(logN_bin, hits_bin, total_bin, p0_lw, bounds_lw):
    """
    Fit logistic to binomial counts via maximum likelihood.
    Returns (popt, perr, pcov, chi2_dof_like) in (b, x0, w) parameterization.
    """
    def nll(params):
        b, x0, log_w = params
        w = np.exp(log_w)
        p = logistic_logN(logN_bin, b, x0, w)
        p = np.clip(p, 1e-12, 1 - 1e-12)
        return -np.sum(hits_bin * np.log(p) + (total_bin - hits_bin) * np.log(1 - p))

    scipy_bounds = list(zip(bounds_lw[0], bounds_lw[1]))
    res_fit = minimize(nll, x0=p0_lw, method="L-BFGS-B",
                       bounds=scipy_bounds, options={"maxiter": 5000})

    if not res_fit.success:
        raise RuntimeError(f"NLL fit failed: {res_fit.message}")

    popt_lw = res_fit.x
    # Covariance from inverse Hessian (BFGS approximation)
    try:
        hess_inv = res_fit.hess_inv.todense() if hasattr(res_fit.hess_inv, "todense") else np.array(res_fit.hess_inv)
        pcov_lw = hess_inv
    except Exception:
        pcov_lw = np.eye(3) * np.nan

    # Convert to (b, x0, w) space
    popt = np.array([popt_lw[0], popt_lw[1], np.exp(popt_lw[2])])
    perr_lw = np.sqrt(np.maximum(np.diag(pcov_lw), 0.0))
    perr = np.array([perr_lw[0], perr_lw[1], np.exp(popt_lw[2]) * perr_lw[2]])

    pcov = pcov_lw.copy()
    pcov[2, :] *= np.exp(popt_lw[2])
    pcov[:, 2] *= np.exp(popt_lw[2])

    # Pseudo chi2/dof using binomial deviance
    p_fit = np.clip(logistic_logN(logN_bin, *popt), 1e-12, 1 - 1e-12)
    p_obs = np.clip(hits_bin / np.maximum(total_bin, 1), 1e-12, 1 - 1e-12)
    dev = 2.0 * np.sum(
        hits_bin * np.log(p_obs / p_fit) +
        (total_bin - hits_bin) * np.log((1 - p_obs) / (1 - p_fit))
    )
    dof = len(logN_bin) - 3
    chi2_dof = dev / dof if dof > 0 else np.nan
    return popt, perr, pcov, chi2_dof


def logistic_logN(logN, b, x0, w):
    exponent = np.clip(-(logN - x0) / w, -500, 500)
    return b / (1.0 + np.exp(exponent))

# ---- helper: compute bootstrap eff/err arrays for a slice ----
def compute_bootstrap_arrays(pc_slice, trig_slice, pc_edges_v, rng=_RNG):
    n_bins = len(pc_edges_v) - 1

    idx = np.digitize(pc_slice, pc_edges_v) - 1
    valid = (idx >= 0) & (idx < n_bins)
    idx = idx[valid]
    t = trig_slice[valid]

    total = np.bincount(idx, minlength=n_bins).astype(int)
    hits = np.bincount(idx, weights=t, minlength=n_bins)

    with np.errstate(invalid="ignore", divide="ignore"):
        p_hat = np.where(total > 0, hits / np.maximum(total, 1), np.nan)

    # Clopper-Pearson CI (replaces bootstrap)
    alpha_cp = 1.0 - CI_LEVEL / 100.0
    safe_total = np.maximum(total, 1)
    lo, hi = clopper_pearson(hits, safe_total, alpha=alpha_cp)

    # Mask empty bins
    empty = (total == 0)
    p_hat[empty] = np.nan
    lo = np.where(empty, np.nan, lo)
    hi = np.where(empty, np.nan, hi)

    err_lo = p_hat - lo
    err_hi = hi - p_hat
    err = np.maximum((hi - lo) / 2.0, EFF_ERR_FLOOR)
    err[empty] = np.nan

    # Also return hits for use in likelihood fit (Step 2)
    return p_hat, err, err_lo, err_hi, total, hits

def fit_one_inc_bin(inc_lo, inc_hi, pc_edges_v, rng):
    inc_center = 0.5 * (inc_lo + inc_hi)
    sub = (inc >= inc_lo) & (inc < inc_hi)
    n_sub = int(sub.sum())

    res = {
        "inc_lo": inc_lo,
        "inc_hi": inc_hi,
        "inc_center": inc_center,
        "n_sub": n_sub,
        "ok": False,
        "skip_reason": None,
    }

    if n_sub < MIN_EVENTS_PER_INC_BIN:
        res["skip_reason"] = "low_events"
        return res

    pc_sub = pc[sub]
    trig_sub = trig[sub]

    eff, err, err_lo, err_hi, total, hits = compute_bootstrap_arrays(
        pc_sub, trig_sub, pc_edges_v, rng=rng
    )
    pc_centers = np.sqrt(pc_edges_v[:-1] * pc_edges_v[1:])
    res.update({
        "pc_centers": pc_centers,
        "eff": eff,
        "err": err,
        "err_lo": err_lo,
        "err_hi": err_hi,
        "total": total,
        "hits": hits,
    })

    eff_valid = eff[np.isfinite(eff)]
    saturated = bool(eff_valid.size > 0 and np.all(eff_valid >= SATURATION_EFF))
    good = np.isfinite(eff) & (total >= MIN_EVENTS_PER_PC_BIN)

    # Need some triggers somewhere in this slice to constrain x0
    n_triggers = int(hits.sum())
    if n_triggers < MIN_TRIGGERS_PER_INC:
        res["skip_reason"] = "no_triggers"
        return res

    res.update({
        "good": good,
        "saturated": saturated,
        "n_triggers": n_triggers,
        "min_total": int(np.min(total[total > 0])) if np.any(total > 0) else 0,
    })

    if good.sum() < 3:
        res["skip_reason"] = "low_stats"
        return res

    if saturated:
        res["skip_reason"] = "saturated"
        return res

    # ---- Prepare data for binomial NLL fit ----
    xd = np.log10(pc_centers[good])
    hits_g = hits[good].astype(float)
    total_g = total[good].astype(float)

    # Initial guess from observed efficiencies
    yd_init = hits_g / np.maximum(total_g, 1.0)
    x_range = xd.max() - xd.min()
    b0 = min(max(np.nanmax(yd_init), 0.05), 1.0)
    x0_0 = xd[np.argmin(np.abs(yd_init - 0.5 * b0))]
    w0 = max(0.1 * x_range, 1e-3)

    p0_lw = [b0, x0_0, np.log(max(w0, 1e-3))]
    bounds_lw = (
        [0.0,  xd.min() - 2.0, np.log(1e-3)],
        [1.05, xd.max() + 3.0, np.log(5.0)],
    )

    try:
        popt, perr, pcov, chi2_dof = fit_binomial_logistic(
            xd, hits_g, total_g, p0_lw, bounds_lw
        )
    except Exception:
        res["skip_reason"] = "fit_failed"
        return res

    # Floor on w uncertainty (always positive, physically meaningful)
    perr[2] = max(perr[2], W_ERR_FLOOR)

    # Flag extrapolated x0 (true 50% point lies outside the data range)
    x0_fit = popt[1]
    x0_extrapolated = bool((x0_fit < xd.min()) or (x0_fit > xd.max()))

    # Quality cuts
    ok = True
    if not np.isfinite(chi2_dof) or chi2_dof >= MAX_CHI2_DOF:
        ok = False
    if perr[1] >= MAX_X0_ERR_LOGN or perr[2] >= MAX_W_ERR_LOGN:
        ok = False

    res.update({
        "popt": popt,
        "perr": perr,
        "chi2_dof": chi2_dof,
        "x0_extrapolated": x0_extrapolated,
        "ok": ok,
    })

    if not ok:
        res["skip_reason"] = "quality_cuts"

    return res



# ======================================================================
# Bin-variation loop (binning systematic)
# ======================================================================
agg = defaultdict(list)
agg_perr = defaultdict(list)
agg_chi2 = defaultdict(list)
agg_nsub = defaultdict(list)
agg_skip = defaultdict(lambda: defaultdict(int))

log_pc     = np.log10(pc_for_edges)
log_pc_min = np.log10(pc_for_edges.min())
log_pc_max = np.log10(pc_for_edges.max())

bin_results = [] if N_BIN_VARIATIONS == 1 else None
print("max incidence", max(inc))
for v in range(N_BIN_VARIATIONS):
    rng_v = np.random.default_rng(1000 + v)
    if N_BIN_VARIATIONS == 1:
        npc = n_photon_bins
        ninc = n_inc_bins
    else:
        npc = int(rng_v.integers(NPC_RANGE[0], NPC_RANGE[1] + 1))
        ninc = int(rng_v.integers(NINC_RANGE[0], NINC_RANGE[1] + 1))

    #pc_edges_v = np.logspace(log_pc_min, log_pc_max, npc + 1)
    pc_edges_v = 10**np.quantile(log_pc, np.linspace(0, 1, npc + 1))
    #inc_edges_v = np.linspace(inc.min(), inc.max(), ninc + 1)
    inc_edges_v = np.quantile(inc, np.linspace(0, 1, ninc + 1))
    for i in range(ninc):
        inc_lo, inc_hi = inc_edges_v[i], inc_edges_v[i + 1]
        res = fit_one_inc_bin(inc_lo, inc_hi, pc_edges_v, rng_v)

        if N_BIN_VARIATIONS == 1:
            bin_results.append(res)

        if not res.get("ok", False):
            reason = res.get("skip_reason", "unknown")
            key_r = round(res["inc_center"], 1)
            agg_skip[key_r][reason] += 1
            continue

        key = round(res["inc_center"], 1)
        agg[key].append(res["popt"])
        agg_perr[key].append(res["perr"])
        agg_chi2[key].append(res["chi2_dof"])
        agg_nsub[key].append(res["n_sub"])

# ======================================================================
# Aggregate across binning variations
# ======================================================================
keys_sorted = sorted(agg.keys())
inc_means = np.array(keys_sorted)

b_vals  = np.array([np.median([p[0] for p in agg[k]]) for k in keys_sorted])
x0_vals = np.array([np.median([p[1] for p in agg[k]]) for k in keys_sorted])
w_vals  = np.array([np.median([p[2] for p in agg[k]]) for k in keys_sorted])


def _combined_err(popts, perrs, idx):
    vals = np.array([p[idx] for p in popts])
    stat = np.array([e[idx] for e in perrs])
    syst = vals.std(ddof=1) if len(vals) > 1 else 0.0
    return np.sqrt(np.median(stat) ** 2 + syst ** 2)


b_errs  = np.array([_combined_err(agg[k], agg_perr[k], 0) for k in keys_sorted])
x0_errs = np.array([_combined_err(agg[k], agg_perr[k], 1) for k in keys_sorted])
w_errs  = np.array([_combined_err(agg[k], agg_perr[k], 2) for k in keys_sorted])

chi2_vals = np.array([np.median(agg_chi2[k]) for k in keys_sorted])
n_variations_per_key = np.array([len(agg[k]) for k in keys_sorted])
min_counts = np.array([int(np.median(agg_nsub[k])) for k in keys_sorted])

finite_mask = np.isfinite(b_vals) & np.isfinite(x0_vals) & np.isfinite(w_vals)
stats_mask = n_variations_per_key >= N_var_per_key_min #max(3, N_BIN_VARIATIONS // 10)
chi2_mask = chi2_vals < MAX_CHI2_DOF
quality_mask = finite_mask & stats_mask & chi2_mask

# ======================================================================
# Diagnostics print
# ======================================================================
print("=" * 80)
print("PER-BIN LOGISTIC FIT DIAGNOSTICS (eff vs photon count, per incidence bin)")
print("=" * 80)
print("Fit equation: eff(logN) = b / (1 + exp(-(logN - x0) / w))")
print(f"Error bars:   bootstrap {CI_LEVEL:.0f}% CI (N_boot={N_BOOTSTRAP})")
print("-" * 80)
print(f"{'inc_mean':>10} {'b(sat)':>10} {'+/-b':>10} {'x0(logN50)':>12} "
      f"{'+/-x0':>10} {'w(logN)':>10} {'+/-w':>10} {'chi2/dof':>10} {'minN':>8}")
print("-" * 80)
for i in range(len(inc_means)):
    print(f"{inc_means[i]:10.3f} "
          f"{b_vals[i]:10.3g} {b_errs[i]:10.2g} "
          f"{x0_vals[i]:12.3f} {x0_errs[i]:10.2g} "
          f"{w_vals[i]:10.3f} {w_errs[i]:10.2g} "
          f"{chi2_vals[i]:10.2f} {min_counts[i]:8d}")
print("-" * 80)
print(f"Bins fitted: {len(inc_means)} | Bins used for global fits: {int(np.sum(quality_mask))}")
print("=" * 80)
print("\nBIN STABILITY ACROSS BINNING VARIATIONS (skip_reason counts):")
print(f"{'inc_mean':>10} {'ok_count':>10} {'low_events':>12} {'low_stats':>10} "
      f"{'saturated':>10} {'fit_failed':>10} {'quality_cuts':>13} {'unknown':>8}")
print("-" * 90)
for k in keys_sorted:
    n_ok  = len(agg[k])
    skips = agg_skip.get(k, {})
    print(f"{k:10.3f} {n_ok:10d} "
          f"{skips.get('low_events', 0):12d} "
          f"{skips.get('low_stats', 0):10d} "
          f"{skips.get('saturated', 0):10d} "
          f"{skips.get('fit_failed', 0):10d} "
          f"{skips.get('quality_cuts', 0):13d} "
          f"{skips.get('unknown', 0):8d}")
print("-" * 90)
print(f"Note: ok_count + all skip counts should sum to approximately "
      f"N_BIN_VARIATIONS={N_BIN_VARIATIONS} per row (variation in ninc "
      f"means not every variation produces an inc bin at that center).")
dropped = ~quality_mask
if np.any(dropped):
    print("\nBins EXCLUDED from global fits:")
    print(f"{'inc_mean':>10} {'b/err_b':>10} {'chi2/dof':>10} {'minN':>8}  reasons")
    min_req = N_var_per_key_min #max(3, N_BIN_VARIATIONS // 10)
    for i in np.where(dropped)[0]:
        reasons = []
        if not finite_mask[i]:
            reasons.append("nonfinite")
        if not stats_mask[i]:
            reasons.append(f"minN<{min_req}")
        if not chi2_mask[i]:
            reasons.append(f"chi2>{MAX_CHI2_DOF}")
        bs = b_vals[i] / b_errs[i] if np.isfinite(b_errs[i]) and b_errs[i] > 0 else 0.0
        print(f"{inc_means[i]:10.3f} {bs:10.2f} {chi2_vals[i]:10.2f} "
              f"{min_counts[i]:8d}  " + ", ".join(reasons))

if len(inc_means) == 0:
    raise RuntimeError("No valid incidence bins available for global fits.")

# ======================================================================
# Global models: b(inc), x0(inc), w(inc)
# ======================================================================

def b_const(inc_deg, b0):
    return np.full_like(inc_deg, b0, dtype=float)


def b_linear(inc_deg, b0, s):
    return b0 + s * inc_deg

def b_sigmoid(inc_deg, b0=1.0, b_inf=0.5, inc_mid=4.75, k=4.0):
    """b at low inc -> b_inf at high inc, transition at inc_mid with steepness k."""
    return b_inf + (b0 - b_inf) / (1.0 + np.exp(k * (inc_deg - inc_mid)))

B0_FIXED = 1.0      # set to whatever value you want
B_INF_FIXED = 0.0   # set to whatever value you want

def b_sigmoid_fixed(inc_deg, inc_mid, k):
    return B_INF_FIXED + (B0_FIXED - B_INF_FIXED) / (1.0 + np.exp(k * (inc_deg - inc_mid)))


def x0_linear(inc_deg, a, s):
    return a + s * inc_deg


def x0_quadratic(inc_deg, a, s, q):
    return a + s * inc_deg + q * inc_deg**2


def x0_log_sigmoid(inc_deg, x_min, x_max, inc_mid, k):
    return x_min + (x_max - x_min) / (1.0 + np.exp(-k * (inc_deg - inc_mid)))


def w_const(inc_deg, w0):
    return np.full_like(inc_deg, w0, dtype=float)


def w_linear(inc_deg, w0, s):
    return w0 + s * inc_deg


def w_powerlaw(inc_deg, w0, k):
    return w0 * (1.0 + inc_deg) ** k


def w_lognormal(inc_deg, w0, mu, sigma):
    # Gaussian in log(inc); guard inc>0
    x = np.log(np.maximum(inc_deg, 1e-6))
    return w0 * np.exp(-((x - mu) / sigma) ** 2)

def w_sigmoid(inc_deg, w0, w_inf, inc_mid, k):
    # w0 at low inc, w_inf at high inc, transition at inc_mid with steepness k
    return w_inf + (w0 - w_inf) / (1.0 + np.exp(k * (inc_deg - inc_mid)))

inc_min = float(np.min(inc_means)) if len(inc_means) else float(np.min(inc))
inc_max = 10 #float(np.max(inc_means)) if len(inc_means) else float(np.max(inc))
x0_min = float(np.nanmin(x0_vals)) if len(x0_vals) else np.log10(pc_for_edges.min())
x0_max = float(np.nanmax(x0_vals)) if len(x0_vals) else np.log10(pc_for_edges.max())

b_model_cfg = {
    "const": {
        "fn": b_const,
        "p0": [np.nanmedian(b_vals)],
        "bounds": ([0.0], [1.05]),
        "label": "b(inc) = b0",
        "param_names": ["b0"],
    },
    "linear": {
        "fn": b_linear,
        "p0": [np.nanmedian(b_vals), 0.0],
        "bounds": ([0.0, -0.1], [1.05, 0.1]),
        "label": "b(inc) = b0 + s*inc",
        "param_names": ["b0", "s"],
    },
    "sigmoid": {
        "fn": b_sigmoid,
        "p0": [1.0, 0.0, 4.75, 4.0],
        "bounds": (
            [0.5, 0.0, inc_min, 0.0001],
            [1.05, 0.0001, inc_max, 10.0],
        ),
        "label": "b(inc) = b_inf + (b0 - b_inf) / (1 + exp(k*(inc - inc_mid)))",
        "param_names": ["b0", "b_inf", "inc_mid", "k"],
    },
    "sigmoid_fixed": {
        "fn": b_sigmoid_fixed,
        "p0": [4.75, 4.0],
        "bounds": (
            [inc_min, 0.0001],
            [inc_max, 10.0],
        ),
        "label": f"b(inc) = {B_INF_FIXED} + ({B0_FIXED}-{B_INF_FIXED}) / (1 + exp(k*(inc - inc_mid)))",
        "param_names": ["inc_mid", "k"],
    },
}

x0_model_cfg = {
    "linear": {
        "fn": x0_linear,
        "p0": [np.nanmedian(x0_vals), 0.0],
        "bounds": ([x0_min - 2.0, -1.0], [x0_max + 2.0, 1.0]),
        "label": "x0(inc) = a + s*inc",
        "param_names": ["a", "s"],
    },
    "quadratic": {
        "fn": x0_quadratic,
        "p0": [np.nanmedian(x0_vals), 0.0, 0.0],
        "bounds": ([x0_min - 2.0, -1.0, -0.1], [x0_max + 2.0, 1.0, 0.1]),
        "label": "x0(inc) = a + s*inc + q*inc^2",
        "param_names": ["a", "s", "q"],
    },
    "log_sigmoid": {
        "fn": x0_log_sigmoid,
        "p0": [x0_min, x0_max, 0.5 * (inc_min + inc_max), 1.0],
        "bounds": ([x0_min - 2.0, x0_min, inc_min, 0.01],
                   [x0_max, x0_max + 2.0, inc_max, 10.0]),
        "label": "x0(inc) = x_min + (x_max-x_min)/(1+exp(-k*(inc-inc_mid)))",
        "param_names": ["x_min", "x_max", "inc_mid", "k"],
    },
}

w_model_cfg = {
    "const": {
        "fn": w_const,
        "p0": [np.nanmedian(w_vals)],
        "bounds": ([1e-3], [5.0]),
        "label": "w(inc) = w0",
        "param_names": ["w0"],
    },
    "linear": {
        "fn": w_linear,
        "p0": [np.nanmedian(w_vals), 0.0],
        "bounds": ([1e-3, -1.0], [5.0, 1.0]),
        "label": "w(inc) = w0 + s*inc",
        "param_names": ["w0", "s"],
    },
    "powerlaw": {
        "fn": w_powerlaw,
        "p0": [np.nanmedian(w_vals), 0.0],
        "bounds": ([1e-3, -3.0], [5.0, 3.0]),
        "label": "w(inc) = w0 * (1 + inc)^k",
        "param_names": ["w0", "k"],
    },
    "lognormal": {
        "fn": w_lognormal,
        "p0": [np.nanmax(w_vals),
               np.log(max(inc_means[np.nanargmax(w_vals)], 1e-3)),
               1.0],
        "bounds": ([1e-3, np.log(max(inc_min, 1e-3)) - 5.0, 1e-2],
                   [5.0,  np.log(max(inc_max, 1e-3)) + 5.0, 10.0]),
        "label": "w(inc) = w0 * exp(-((log(inc) - mu)/sigma)^2)",
        "param_names": ["w0", "mu", "sigma"],
    },
    "sigmoid": {
        "fn": w_sigmoid,
        "p0": [
            np.nanmedian(w_vals),                 # w0: low-inc plateau
            min(0.01, np.nanmin(w_vals)),         # w_inf: high-inc asymptote
            0.5 * (inc_min + inc_max),            # inc_mid: transition point
            5.0,                                  # k: steepness
        ],
        "bounds": (
            [1e-3, 0.0, inc_min, 0.01],
            [5.0,  5.0, inc_max, 100.0],
        ),
        "label": "w(inc) = w_inf + (w0 - w_inf) / (1 + exp(k*(inc - inc_mid)))",
        "param_names": ["w0", "w_inf", "inc_mid", "k"],
    },
}

if B_MODEL not in b_model_cfg:
    raise ValueError(f"Unknown B_MODEL='{B_MODEL}'. Use 'const' or 'linear', or 'sigmoid', or 'sigmoid_fixed'.")
if X0_MODEL not in x0_model_cfg:
    raise ValueError(f"Unknown X0_MODEL='{X0_MODEL}'. Use 'linear', 'quadratic', or 'log_sigmoid'.")
if W_MODEL not in w_model_cfg:
    raise ValueError(f"Unknown W_MODEL='{W_MODEL}'. Use 'const', 'linear', 'powerlaw', 'lognormal', or 'sigmoid'.")


b_model = b_model_cfg[B_MODEL]
x0_model = x0_model_cfg[X0_MODEL]
w_model = w_model_cfg[W_MODEL]


def fit_model(fn, xdata, ydata, yerr, p0, bounds):
    sigma = np.where(yerr > 0, yerr, EFF_ERR_FLOOR)
    popt, pcov = curve_fit(
        fn, xdata, ydata, p0=p0, sigma=sigma,
        absolute_sigma=True, bounds=bounds, maxfev=20000
    )
    perr = np.sqrt(np.diag(pcov))
    residuals = ydata - fn(xdata, *popt)
    chi2 = np.sum((residuals / sigma) ** 2)
    dof = len(xdata) - len(popt)
    chi2_dof = chi2 / dof if dof > 0 else np.nan
    return popt, perr, pcov, chi2_dof


results = {}
results_meta = {}

ok_b = quality_mask & np.isfinite(b_errs)
ok_x0 = quality_mask & np.isfinite(x0_errs)
ok_w = quality_mask & np.isfinite(w_errs)

try:
    popt, perr, pcov, c2 = fit_model(
        b_model["fn"], inc_means[ok_b], b_vals[ok_b],
        np.maximum(b_errs[ok_b], B_ERR_FLOOR),
        p0=b_model["p0"], bounds=b_model["bounds"],
    )
    results["b"] = (popt, perr, pcov)
    results_meta["b"] = {"chi2_dof": c2, "model": B_MODEL}
except Exception as e:
    print(f"b model fit failed: {e}")

try:
    popt, perr, pcov, c2 = fit_model(
        x0_model["fn"], inc_means[ok_x0], x0_vals[ok_x0],
        np.maximum(x0_errs[ok_x0], W_ERR_FLOOR),
        p0=x0_model["p0"], bounds=x0_model["bounds"],
    )
    results["x0"] = (popt, perr, pcov)
    results_meta["x0"] = {"chi2_dof": c2, "model": X0_MODEL}
except Exception as e:
    print(f"x0 model fit failed: {e}")

try:
    popt, perr, pcov, c2 = fit_model(
        w_model["fn"], inc_means[ok_w], w_vals[ok_w],
        np.maximum(w_errs[ok_w], W_ERR_FLOOR),
        p0=w_model["p0"], bounds=w_model["bounds"],
    )
    results["w"] = (popt, perr, pcov)
    results_meta["w"] = {"chi2_dof": c2, "model": W_MODEL}
except Exception as e:
    print(f"w model fit failed: {e}")
# ======================================================================
# Empirical cross-group correlation of (b, x0, w) from binning variations
# ======================================================================
# For each quality-masked incidence bin, agg[k] holds a list of (b, x0, w)
# triplets — one per binning variation that passed quality cuts. The empirical
# 3x3 covariance of those triplets captures how b, x0, w co-vary in practice.
# Pooling across bins gives a representative correlation matrix that Cell 4
# uses to populate the off-diagonal blocks of the joint parameter covariance.

_emp_covs_list = []
for k, qm in zip(keys_sorted, quality_mask):
    if not qm:
        continue
    if len(agg[k]) < 4:          # need at least 4 samples for a meaningful estimate
        continue
    _triplets = np.array(agg[k])  # (n_variations, 3): columns = [b, x0, w]
    _emp_covs_list.append(np.cov(_triplets.T))   # (3, 3)

if len(_emp_covs_list) > 0:
    _EMP_JOINT_COV_3x3_raw = np.mean(_emp_covs_list, axis=0)   # (3, 3) pooled

    # Convert to correlation matrix: C_ij / (sigma_i * sigma_j)
    _std3 = np.sqrt(np.maximum(np.diag(_EMP_JOINT_COV_3x3_raw), 0.0))
    _outer_std = np.outer(_std3, _std3)
    EMP_JOINT_CORR_3x3 = np.where(
        _outer_std > 0,
        _EMP_JOINT_COV_3x3_raw / _outer_std,
        0.0,
    )
    # Clip to valid correlation range for numerical safety
    np.clip(EMP_JOINT_CORR_3x3, -1.0, 1.0, out=EMP_JOINT_CORR_3x3)
    # Force exact diagonal = 1 (rounding can push it slightly off)
    np.fill_diagonal(EMP_JOINT_CORR_3x3, 1.0)

    print("\nEmpirical cross-group correlation matrix (b, x0, w):")
    print(f"  Pooled from {len(_emp_covs_list)} quality bins "
          f"(min 4 variations each)")
    print(f"  ρ(b,  x0) = {EMP_JOINT_CORR_3x3[0, 1]:.4f}")
    print(f"  ρ(b,   w) = {EMP_JOINT_CORR_3x3[0, 2]:.4f}")
    print(f"  ρ(x0,  w) = {EMP_JOINT_CORR_3x3[1, 2]:.4f}")
else:
    EMP_JOINT_CORR_3x3 = np.eye(3, dtype=float)
    print("\nWarning: no quality bins had >= 4 binning variations. "
          "EMP_JOINT_CORR_3x3 set to identity (zero cross-group correlations).")

del _emp_covs_list
# ======================================================================
# Comparison plot: coefficients vs incidence
# ======================================================================
fig, axes = plt.subplots(2, 2, figsize=(20, 10))
axes = axes.ravel()

inc_smooth = np.linspace(inc_means.min(), inc_means.max(), 300)

labels = ["b (sat)", "x0 (logN50)", "w (w_logN)"]
yvals = [b_vals, x0_vals, w_vals]
yerrs = [b_errs, x0_errs, w_errs]
colors = ["darkorange", "green", "red"]
model_keys = ["b", "x0", "w"]
model_fns = [b_model["fn"], x0_model["fn"], w_model["fn"]]
model_labels = {
    "b": b_model["label"],
    "x0": x0_model["label"],
    "w": w_model["label"],
}
param_names = {
    "b": b_model["param_names"],
    "x0": x0_model["param_names"],
    "w": w_model["param_names"],
}

for ax, lbl, yv, ye, col, mkey, mfn in zip(
        axes[:len(labels)], labels, yvals, yerrs, colors, model_keys, model_fns):

    qm = quality_mask
    ax.errorbar(inc_means[qm], yv[qm], yerr=ye[qm], fmt="o", color=col,
                capsize=3, label="fitted value", zorder=3)
    sc = ax.scatter(inc_means[qm], yv[qm], c=chi2_vals[qm], cmap="Reds",
                    vmin=0, vmax=5, zorder=5, s=60)
    plt.colorbar(sc, ax=ax, label="chi2/dof (per-bin fit)")
    ax.scatter(inc_means[qm], yv[qm],
               facecolors="none", edgecolors="black", s=90,
               linewidths=1.2, label="used in global fit")

    if mkey in results:
        popt, perr, _ = results[mkey]
        c2 = results_meta.get(mkey, {}).get("chi2_dof", np.nan)
        ax.plot(inc_smooth, mfn(inc_smooth, *popt), "k--",
                lw=1.5, label=model_labels[mkey])
        txt = "\n".join(f"{n}={v:.3g}+/-{e:.2g}"
                        for n, v, e in zip(param_names[mkey], popt, perr))
        txt += f"\nchi2/dof={c2:.2f}"
        ax.text(0.97, 0.05, txt, transform=ax.transAxes, fontsize=8,
                ha="right", va="bottom",
                bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    ax.set_xlabel("incidence angle (deg)")
    ax.set_ylabel(lbl)
    ax.set_title(lbl, fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

for ax in axes[len(labels):]:
    ax.axis("off")

fig.suptitle("Fit coefficients vs incidence angle", fontsize=13)
plt.show()

# ======================================================================
# Physical model summary
# ======================================================================
print("\n" + "=" * 80)
print("PHYSICAL MODEL FITS (coefficient vs incidence angle)")
print("=" * 80)
model_equations = {
    "b": b_model["label"],
    "x0": x0_model["label"],
    "w": w_model["label"],
}
for key in ["b", "x0", "w"]:
    print(f"\n--- Model for parameter '{key}' ---")
    print(f"Equation:  {model_equations[key]}")
    if key in results:
        popt, perr, _ = results[key]
        c2 = results_meta.get(key, {}).get("chi2_dof", np.nan)
        for name, val, err in zip(param_names[key], popt, perr):
            print(f"  {name:>12} = {val:.4g} +/- {err:.2g}")
        print(f"  {'chi2/dof':>12} = {c2:.3f}")
    else:
        print("  FIT FAILED or not available")
print("=" * 80)

# ======================================================================
# Per-incidence-bin efficiency curves and logistic fits
# ======================================================================
if N_BIN_VARIATIONS == 1:
    nrows = int(np.ceil(n_inc_bins / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 3 * nrows), sharex=True, sharey=True
    )
    axes = np.atleast_1d(axes).ravel()

    N_smooth = np.logspace(np.log10(pc.min()), np.log10(pc.max()), 300)

    for i in range(n_inc_bins):
        ax = axes[i]
        r = bin_results[i]
        inc_lo, inc_hi = r["inc_lo"], r["inc_hi"]
        n_sub = r["n_sub"]
        inc_c = r["inc_center"]

        if n_sub == 0:
            ax.set_title(f"inc in [{inc_lo:.2f},{inc_hi:.2f}]\nN=0")
            continue

        if r.get("skip_reason") in ("low_events", "low_stats", "fit_failed", "saturated", "quality_cuts"):
            ax.set_title(f"inc in [{inc_lo:.2f},{inc_hi:.2f}]\nN={n_sub}")
            ax.text(0.5, 0.5, r.get("skip_reason", "skip"),
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="darkred")
            continue

        pc_centers = r["pc_centers"]
        eff = r["eff"]
        err_lo = r["err_lo"]
        err_hi = r["err_hi"]

        yerr_asym = np.array([
            np.where(np.isfinite(err_lo), err_lo, 0),
            np.where(np.isfinite(err_hi), err_hi, 0),
        ])

        ax.errorbar(pc_centers, eff, yerr=yerr_asym, fmt="o", capsize=2,
                    label=f"bootstrap {CI_LEVEL:.0f}% CI")

        popt = r["popt"]
        chi2_dof = r["chi2_dof"]
        ax.plot(N_smooth, logistic_logN(np.log10(N_smooth), *popt), "r-", lw=1.5)

        b, x0, w = popt
        txt = (f"b={b:.3f}\n"
               f"x0={x0:.2f}\n"
               f"w={w:.2f}\n"
               f"chi2/dof={chi2_dof:.2f}")
        ax.text(0.05, 0.95, txt, transform=ax.transAxes,
                fontsize=7, va="top", ha="left",
                bbox=dict(boxstyle="round", fc="white", alpha=0.75))

        ax.set_title(f"inc≈{inc_c:.2f}° (n={n_sub})", fontsize=9)
        ax.set_xscale("log")
        ax.set_xlabel("photon_count")
        ax.set_ylabel(f"Trigger eff. (Max_PE >= {PE_THRESHOLD})")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower right")

    for j in range(n_inc_bins, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Trigger efficiency vs photon count | bootstrap {CI_LEVEL:.0f}% CI", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

# --- Plot representative per-incidence-bin fits (median binning) ---
npc_plot = int(np.median(NPC_RANGE))
ninc_plot = int(np.median(NINC_RANGE))

#pc_edges_plot = np.logspace(log_pc_min, log_pc_max, npc_plot + 1)
pc_edges_plot = 10**np.quantile(log_pc, np.linspace(0, 1, npc_plot + 1))
#inc_edges_plot = np.linspace(inc.min(), inc.max(), ninc_plot + 1)
inc_edges_plot = np.quantile(inc, np.linspace(0, 1, ninc_plot + 1))

plot_results = []
for i in range(ninc_plot):
    res = fit_one_inc_bin(inc_edges_plot[i], inc_edges_plot[i + 1],
                          pc_edges_plot, _RNG)
    if res is not None and res.get("ok", False):
        plot_results.append(res)

n_plots = len(plot_results)
if n_plots == 0:
    print(res)
    print("No successful per-incidence-bin fits to plot.")
else:
    ncols_p = 4
    nrows_p = int(np.ceil(n_plots / ncols_p))
    fig, axes = plt.subplots(nrows_p, ncols_p,
                             figsize=(4 * ncols_p, 3 * nrows_p),
                             squeeze=False)
    axes = axes.flatten()

    N_smooth = np.logspace(np.log10(pc.min()), np.log10(pc.max()), 300)

    for ax, res in zip(axes, plot_results):
        pc_c = res["pc_centers"]
        eff = res["eff"]
        elo = res["err_lo"]
        ehi = res["err_hi"]
        popt = res["popt"]
        perr = res["perr"]
        c2 = res["chi2_dof"]
        inc_c = res["inc_center"]

        yerr = np.vstack([elo, ehi])

        ax.errorbar(pc_c, eff, yerr=yerr, fmt="o", color="C0",
                    capsize=2, ms=4, label="data")
        ax.plot(N_smooth, logistic_logN(np.log10(N_smooth), *popt),
                "r-", lw=1.5, label="fit")

        b, x0, w = popt
        be, x0e, we = perr
        txt = (f"b={b:.3f}+/-{be:.2g}\n"
               f"x0={x0:.2f}+/-{x0e:.2g}\n"
               f"w={w:.2f}+/-{we:.2g}\n"
               f"chi2/dof={c2:.2f}")
        ax.text(0.03, 0.95, txt, transform=ax.transAxes,
                fontsize=7, va="top", ha="left",
                bbox=dict(boxstyle="round", fc="white", alpha=0.75))

        ax.set_title(f"inc≈{inc_c:.2f}° (n={res['n_sub']})", fontsize=9)
        ax.set_xlabel("photon_count")
        ax.set_ylabel("efficiency")
        ax.set_xscale("log")
        ax.grid(alpha=0.3)

    for ax in axes[n_plots:]:
        ax.axis("off")

    fig.suptitle(
        f"Per-incidence-bin logistic fits "
        f"(npc={npc_plot}, ninc={ninc_plot}, {n_plots}/{ninc_plot} bins fit OK)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()


# Bootstrap with combined statistical + binning systematics
# Empirical region : bootstrap resampling (unchanged)
# Extrapolation region : Jacobian / analytic error propagation (replaces _PARAM_CUBE)
# Blending (quadrature, lo/hi separately):
#   lo = nominal - sqrt( COV^2*(nominal-lo_emp)^2 + (1-COV)^2*(z*sigma_jac)^2 )
#   hi = nominal + sqrt( COV^2*(hi_emp-nominal)^2 + (1-COV)^2*(z*sigma_jac)^2 )
# =============================================================================
import os
import warnings
import numpy as np
from joblib import Parallel, delayed
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.stats import norm as _scipy_norm
import sympy as sp

# ---------- Config ----------
N_BOOT   = N_BOOTSTRAP
N_JOBS   = os.cpu_count()
SEED     = 0
MIN_CNT  = 5

NPC_RANGE    = (5, 100)
NINC_RANGE   = (5, 100)
SMOOTH_RANGE = (0.7, 1.2)
LOG10_PC_MIN_RANGE = (np.log10(pc_.min()), 1)
MIN_LOW_N_FRACTION = 0.05

if "B_MODEL" not in globals():
    B_MODEL = "const"
if "X0_MODEL" not in globals():
    X0_MODEL = "linear"
if "W_MODEL" not in globals():
    W_MODEL = "const"

if "results" not in globals():
    raise RuntimeError("Run the per-incidence fit cell before this bootstrap cell.")

N_COMMON_PC_DATA  = 50
N_COMMON_INC_DATA = 50


# =============================================================================
# Empirical interpolator helpers  (UNCHANGED)
# =============================================================================

def build_interp_from_edges(pc_s, inc_s, hit_s,
                            pc_edges, inc_edges, smooth, min_cnt=MIN_CNT):
    pc_cen  = np.sqrt(pc_edges[:-1] * pc_edges[1:])
    inc_cen = 0.5 * (inc_edges[:-1] + inc_edges[1:])
    N_tot, _, _ = np.histogram2d(pc_s, inc_s, bins=[pc_edges, inc_edges])
    N_hit, _, _ = np.histogram2d(pc_s, inc_s, bins=[pc_edges, inc_edges],
                                 weights=hit_s)
    with np.errstate(invalid='ignore', divide='ignore'):
        e = np.where(N_tot >= min_cnt, N_hit / N_tot, np.nan)
    has_triggered = N_hit >= 1
    e = np.where(has_triggered, e, np.nan)
    mask     = np.isfinite(e)
    e_filled = np.where(mask, e, 0.0)
    w        = mask.astype(float)
    e_smooth = gaussian_filter(e_filled, sigma=smooth) / np.clip(
        gaussian_filter(w, sigma=smooth), 1e-6, None)
    return RegularGridInterpolator(
        (np.log10(pc_cen), inc_cen), e_smooth,
        bounds_error=False, fill_value=None, method='linear')


def _one_boot(seed, pc_, inc_, hit_):
    rng  = np.random.default_rng(seed)
    n    = pc_.shape[0]
    idx  = rng.integers(0, n, size=n)
    pcs, incs, hits = pc_[idx], inc_[idx], hit_[idx]
    npc    = rng.integers(NPC_RANGE[0],  NPC_RANGE[1]  + 1)
    ninc   = rng.integers(NINC_RANGE[0], NINC_RANGE[1] + 1)
    smooth = rng.uniform(*SMOOTH_RANGE)
    pcs_trig = pcs[hits > 0]
    if len(pcs_trig) < 10:
        pcs_trig = pcs
    pc_edges  = 10**np.quantile(np.log10(pcs_trig), np.linspace(0, 1, npc + 1))
    inc_edges = np.quantile(incs, np.linspace(0, 1, ninc + 1))
    return build_interp_from_edges(pcs, incs, hits, pc_edges, inc_edges, smooth)


def _eval_interp(bi, pts, shape):
    vals = bi(pts)
    return np.clip(vals, 0.0, 1.0).reshape(shape).astype(np.float32)


# =============================================================================
# Nominal interpolator  (UNCHANGED)
# =============================================================================

pc_trig_mask  = hit_ > 0
pc_edges_nom  = 10**np.quantile(np.log10(pc_[pc_trig_mask]),
                                np.linspace(0, 1, N_COMMON_PC_DATA + 1))
inc_edges_nom = np.quantile(inc_, np.linspace(0, 1, N_COMMON_INC_DATA + 1))
interp = build_interp_from_edges(pc_, inc_, hit_,
                                 pc_edges_nom, inc_edges_nom, smooth=0.8)

# =============================================================================
# Bootstrap replicas — parallel  (UNCHANGED)
# =============================================================================

ss    = np.random.SeedSequence(SEED)
seeds = [s.generate_state(1)[0] for s in ss.spawn(N_BOOT)]
print(f"Running {N_BOOT} bootstrap replicas on {N_JOBS} CPUs...")
boot_interps = Parallel(n_jobs=N_JOBS, verbose=5)(
    delayed(_one_boot)(s, pc_, inc_, hit_) for s in seeds
)
print(f"Done: {len(boot_interps)} interpolators built.")

# =============================================================================
# Extended common grid  (UNCHANGED)
# =============================================================================

_DATA_LOG_PC_LO = np.log10(pc_.min())
_DATA_LOG_PC_HI = np.log10(pc_.max())
_DATA_INC_LO    = inc_.min()
_DATA_INC_HI    = inc_.max()

N_COMMON_PC  = int(N_COMMON_PC_DATA * 2)
N_COMMON_INC = int(N_COMMON_INC_DATA * 5)
COMMON_LOG_PC = np.linspace(_DATA_LOG_PC_LO - 1.0, _DATA_LOG_PC_HI + 2.0, N_COMMON_PC)
COMMON_INC    = np.linspace(0.0, 30.0, N_COMMON_INC)

_PC2D, _INC2D = np.meshgrid(COMMON_LOG_PC, COMMON_INC, indexing='ij')
N2D      = 10.0 ** _PC2D
_INC2D_f = _INC2D.astype(float)

# =============================================================================
# Coverage mask  (UNCHANGED)
# =============================================================================

_pc_trig  = pc_[hit_ > 0]
_inc_trig = inc_[hit_ > 0]

_DATA_LOG_PC_LO = np.log10(_pc_trig.min())
_DATA_LOG_PC_HI = np.log10(_pc_trig.max())
_DATA_INC_LO    = _inc_trig.min()
_DATA_INC_HI    = _inc_trig.max()

_inside = (
    (_PC2D >= _DATA_LOG_PC_LO) & (_PC2D <= _DATA_LOG_PC_HI) &
    (_INC2D >= _DATA_INC_LO)   & (_INC2D <= _DATA_INC_HI)
).astype(np.float32)
COV = gaussian_filter(_inside, sigma=1.0)
COV = np.clip(COV, 0.0, 1.0).astype(np.float32)

# =============================================================================
# Empirical cube — parallel  (UNCHANGED)
# =============================================================================

_PC2D_flat = np.column_stack([_PC2D.ravel(), _INC2D.ravel()])
_shape_ext  = (N_COMMON_PC, N_COMMON_INC)

_emp_list = Parallel(n_jobs=N_JOBS)(
    delayed(_eval_interp)(bi, _PC2D_flat, _shape_ext)
    for bi in boot_interps
)
_EMP_CUBE = np.stack(_emp_list, axis=0)   # (N_BOOT, N_PC, N_INC), float32
del _emp_list

alpha_pct = (100.0 - CI_LEVEL) / 2.0
_EMP_LO  = np.nanpercentile(_EMP_CUBE, alpha_pct,       axis=0).astype(np.float32)
_EMP_HI  = np.nanpercentile(_EMP_CUBE, 100 - alpha_pct, axis=0).astype(np.float32)
del _EMP_CUBE

# =============================================================================
# Parametric helpers (scalar, used for nominal surface only)
# =============================================================================

_EXP_CLIP = 500.0

def _safe_exp(x):
    return np.exp(np.clip(x, -_EXP_CLIP, _EXP_CLIP))

def _b_of_inc(inc_deg, *params):
    inc_deg = np.asarray(inc_deg, dtype=float)
    if B_MODEL == "const":
        (b0,) = params
        return np.full_like(inc_deg, b0)
    if B_MODEL == "linear":
        b0, s = params
        return b0 + s * inc_deg
    if B_MODEL == "sigmoid":
        b0, b_inf, inc_mid, k = params
        return b_inf + (b0 - b_inf) / (1.0 + _safe_exp(k * (inc_deg - inc_mid)))
    if B_MODEL == "sigmoid_fixed":
        inc_mid, k = params
        return 1 / (1.0 + _safe_exp(k * (inc_deg - inc_mid)))
    raise ValueError(f"Unsupported B_MODEL: {B_MODEL}")

def _x0_of_inc(inc_deg, *params):
    inc_deg = np.asarray(inc_deg, dtype=float)
    if X0_MODEL == "linear":
        a, s = params
        return a + s * inc_deg
    if X0_MODEL == "quadratic":
        a, s, q = params
        return a + s * inc_deg + q * inc_deg**2
    if X0_MODEL == "log_sigmoid":
        x_min, x_max, inc_mid, k = params
        return x_min + (x_max - x_min) / (1.0 + _safe_exp(-k * (inc_deg - inc_mid)))
    raise ValueError(f"Unsupported X0_MODEL: {X0_MODEL}")

def _w_of_inc(inc_deg, *params):
    inc_deg = np.asarray(inc_deg, dtype=float)
    if W_MODEL == "const":
        (w0,) = params
        return np.full_like(inc_deg, w0)
    if W_MODEL == "linear":
        w0, s = params
        return w0 + s * inc_deg
    if W_MODEL == "powerlaw":
        w0, k = params
        return w0 * (1.0 + inc_deg) ** k
    if W_MODEL == "lognormal":
        w0, mu, sigma = params
        x = np.log(np.maximum(inc_deg, 1e-6))
        return w0 * _safe_exp(-((x - mu) / sigma) ** 2)
    if W_MODEL == "sigmoid":
        w0, w_inf, inc_mid, k = params
        return w_inf + (w0 - w_inf) / (1.0 + _safe_exp(k * (inc_deg - inc_mid)))
    raise ValueError(f"Unsupported W_MODEL: {W_MODEL}")

def _parametric_eff(N, inc, b_params, x0_params, w_params):
    b    = np.clip(_b_of_inc(inc, *b_params), 0.0, 1.0)
    x0   = _x0_of_inc(inc, *x0_params)
    w    = np.maximum(_w_of_inc(inc, *w_params), 1e-3)
    logN = np.log10(np.clip(N, 1e-30, None))
    return b / (1.0 + _safe_exp(-(logN - x0) / w))


# =============================================================================
# SymPy Jacobian — build once for the active model variant
# =============================================================================
# Symbolic variables
_sym_logN = sp.Symbol('logN')
_sym_inc  = sp.Symbol('inc')

def _build_b_expr(model):
    """Return (sympy_expr, [param_symbols])"""
    if model == "const":
        b0 = sp.Symbol('b0')
        return b0, [b0]
    if model == "linear":
        b0, s = sp.symbols('b0 s_b')
        return b0 + s * _sym_inc, [b0, s]
    if model == "sigmoid":
        b0, b_inf, inc_mid, k = sp.symbols('b0 b_inf inc_mid_b k_b')
        return (b_inf + (b0 - b_inf) /
                (1 + sp.exp(k * (_sym_inc - inc_mid))),
                [b0, b_inf, inc_mid, k])
    if model == "sigmoid_fixed":
        inc_mid, k = sp.symbols('inc_mid_b k_b')
        #b0_f  = sp.Symbol('B0_FIXED')
        #binf_f = sp.Symbol('B_INF_FIXED')
        return (1 /
                (1 + sp.exp(k * (_sym_inc - inc_mid))),
                [inc_mid, k])
    raise ValueError(f"Unsupported B_MODEL: {model}")

def _build_x0_expr(model):
    if model == "linear":
        a, s = sp.symbols('a_x0 s_x0')
        return a + s * _sym_inc, [a, s]
    if model == "quadratic":
        a, s, q = sp.symbols('a_x0 s_x0 q_x0')
        return a + s * _sym_inc + q * _sym_inc**2, [a, s, q]
    if model == "log_sigmoid":
        x_min, x_max, inc_mid, k = sp.symbols('x_min x_max inc_mid_x0 k_x0')
        return (x_min + (x_max - x_min) /
                (1 + sp.exp(-k * (_sym_inc - inc_mid))),
                [x_min, x_max, inc_mid, k])
    raise ValueError(f"Unsupported X0_MODEL: {model}")

def _build_w_expr(model):
    if model == "const":
        w0 = sp.Symbol('w0')
        return w0, [w0]
    if model == "linear":
        w0, s = sp.symbols('w0 s_w')
        return w0 + s * _sym_inc, [w0, s]
    if model == "powerlaw":
        w0, k = sp.symbols('w0 k_w')
        return w0 * (1 + _sym_inc)**k, [w0, k]
    if model == "lognormal":
        w0, mu, sigma = sp.symbols('w0 mu_w sigma_w')
        x = sp.log(sp.Max(_sym_inc, sp.Float(1e-6)))
        return w0 * sp.exp(-((x - mu) / sigma)**2), [w0, mu, sigma]
    if model == "sigmoid":
        w0, w_inf, inc_mid, k = sp.symbols('w0 w_inf inc_mid_w k_w')
        return (w_inf + (w0 - w_inf) /
                (1 + sp.exp(k * (_sym_inc - inc_mid))),
                [w0, w_inf, inc_mid, k])
    raise ValueError(f"Unsupported W_MODEL: {model}")

print("Building SymPy Jacobian for active model variant "
      f"({B_MODEL} / {X0_MODEL} / {W_MODEL})...")

_b_expr,  _b_syms  = _build_b_expr(B_MODEL)
_x0_expr, _x0_syms = _build_x0_expr(X0_MODEL)
_w_expr,  _w_syms  = _build_w_expr(W_MODEL)
_all_param_syms = _b_syms + _x0_syms + _w_syms

# Full parametric efficiency expression:
#   ε = b(inc) / ( 1 + exp(-(logN - x0(inc)) / w(inc)) )
# Clipping is not differentiable; we work in the unclipped domain and
# rely on the nominal model being well-behaved at the parameter optimum.
_w_safe = sp.Max(_w_expr, sp.Float(1e-3))
_eff_expr = _b_expr / (1 + sp.exp(-(_sym_logN - _x0_expr) / _w_safe))

# Jacobian row: d(eff)/d(theta_i) for each parameter in order
_jac_exprs = [sp.diff(_eff_expr, sym) for sym in _all_param_syms]

# Lambdify for fast numerical evaluation — outputs a list of arrays
_jac_func = sp.lambdify(
    [_sym_logN, _sym_inc] + _all_param_syms,
    _jac_exprs,
    modules='numpy',
)
print(f"  Parameters ({len(_all_param_syms)}): "
      f"{[str(s) for s in _all_param_syms]}")
print("  Jacobian lambdified OK.")

# =============================================================================
# Joint covariance matrix (b, x0, w correlated via EMP_JOINT_CORR_3x3)
# =============================================================================

b_mean,  b_cov  = results['b'][0],  results['b'][2]
x0_mean, x0_cov = results['x0'][0], results['x0'][2]
w_mean,  w_cov  = results['w'][0],  results['w'][2]

_nb  = len(b_mean)
_nx0 = len(x0_mean)
_nw  = len(w_mean)
_n_total = _nb + _nx0 + _nw

_joint_cov = np.zeros((_n_total, _n_total))
_joint_cov[         :_nb,          :_nb         ] = b_cov
_joint_cov[_nb      :_nb + _nx0,   _nb:_nb+_nx0 ] = x0_cov
_joint_cov[_nb+_nx0 :,             _nb+_nx0:    ] = w_cov

if "EMP_JOINT_CORR_3x3" in globals():
    _ρ_b_x0 = float(EMP_JOINT_CORR_3x3[0, 1])
    _ρ_b_w  = float(EMP_JOINT_CORR_3x3[0, 2])
    _ρ_x0_w = float(EMP_JOINT_CORR_3x3[1, 2])
    print(f"[Jacobian] Using EMP_JOINT_CORR_3x3: "
          f"ρ(b,x0)={_ρ_b_x0:.3f}  ρ(b,w)={_ρ_b_w:.3f}  ρ(x0,w)={_ρ_x0_w:.3f}")
else:
    _ρ_b_x0 = _ρ_b_w = _ρ_x0_w = 0.0
    print("[Jacobian] EMP_JOINT_CORR_3x3 not found — using block-diagonal covariance.")

_σ_b  = np.sqrt(np.maximum(np.diag(b_cov),  0.0))
_σ_x0 = np.sqrt(np.maximum(np.diag(x0_cov), 0.0))
_σ_w  = np.sqrt(np.maximum(np.diag(w_cov),  0.0))

_C_b_x0 = _ρ_b_x0 * np.outer(_σ_b,  _σ_x0)
_C_b_w  = _ρ_b_w  * np.outer(_σ_b,  _σ_w)
_C_x0_w = _ρ_x0_w * np.outer(_σ_x0, _σ_w)

_joint_cov[          :_nb,         _nb:_nb+_nx0 ] = _C_b_x0
_joint_cov[_nb       :_nb+_nx0,    :_nb         ] = _C_b_x0.T
_joint_cov[          :_nb,         _nb+_nx0:    ] = _C_b_w
_joint_cov[_nb+_nx0: ,             :_nb         ] = _C_b_w.T
_joint_cov[_nb       :_nb+_nx0,    _nb+_nx0:   ] = _C_x0_w
_joint_cov[_nb+_nx0: ,             _nb:_nb+_nx0 ] = _C_x0_w.T

_joint_cov = 0.5 * (_joint_cov + _joint_cov.T)

# PSD correction
_eigvals, _eigvecs = np.linalg.eigh(_joint_cov)
if np.any(_eigvals < 0):
    _n_neg = int((_eigvals < 0).sum())
    _eigvals = np.maximum(_eigvals, 0.0)
    _joint_cov = (_eigvecs * _eigvals) @ _eigvecs.T
    _joint_cov = 0.5 * (_joint_cov + _joint_cov.T)
    print(f"[Jacobian] PSD correction: clipped {_n_neg} negative eigenvalue(s).")

# =============================================================================
# Jacobian uncertainty surface — σ_jac(logN, inc) on the extended grid
# =============================================================================
# σ²_jac = J(θ*) · Σ_joint · J(θ*)ᵀ   where J is (1 x n_params) at each point

_joint_popt = np.concatenate([b_mean, x0_mean, w_mean])   # best-fit parameter vector

# Evaluate Jacobian at every grid point — shape (n_params, N_PC*N_INC) → (N_PC, N_INC)
_logN_flat = np.log10(np.clip(N2D.ravel(), 1e-30, None))   # (N_PC*N_INC,)
_inc_flat  = _INC2D_f.ravel()                              # (N_PC*N_INC,)

# _jac_func returns a list of n_param arrays each of shape (N_flat,)
# We need to handle SymPy Max which compiles to np.maximum correctly via lambdify.
_jac_vals = _jac_func(_logN_flat, _inc_flat, *_joint_popt)
# Stack → (n_params, N_flat)
_J = np.array([np.broadcast_to(np.asarray(jv, dtype=float), _logN_flat.shape)
               for jv in _jac_vals])   # (n_params, N_flat)

# σ²_jac = diag(J^T · Σ · J) evaluated per grid point:
#   for each column j of J:  σ²_j = j · Σ · j
# Efficient form: (J^T Σ J) diagonal = sum_i sum_k J[i,:] * Σ[i,k] * J[k,:]
#   = einsum('iN,ij,jN->N', J, Σ, J)
_ΣJ       = _joint_cov @ _J           # (n_params, N_flat)
_var_jac  = np.einsum('iN,iN->N', _J, _ΣJ)                # (N_flat,)
_var_jac  = np.maximum(_var_jac, 0.0)                     # guard float noise
_SIGMA_JAC = np.sqrt(_var_jac).reshape(N_COMMON_PC, N_COMMON_INC).astype(np.float32)

print(f"σ_jac surface: min={_SIGMA_JAC.min():.3e}  "
      f"max={_SIGMA_JAC.max():.3e}  "
      f"mean={_SIGMA_JAC.mean():.3e}")
del _J, _ΣJ, _var_jac, _logN_flat, _inc_flat, _jac_vals

# =============================================================================
# Global (nominal) surface  (UNCHANGED)
# =============================================================================

COMMON_LOG_PC_DATA = np.linspace(np.log10(pc_.min()), np.log10(pc_.max()), N_COMMON_PC_DATA)
COMMON_INC_DATA    = np.linspace(inc_.min(), inc_.max(), N_COMMON_INC_DATA)

_EMP_GLOBAL = np.clip(interp(_PC2D_flat), 0.0, 1.0).reshape(
    N_COMMON_PC, N_COMMON_INC).astype(np.float32)
_EMP_GLOBAL = np.where(np.isfinite(_EMP_GLOBAL), _EMP_GLOBAL, 0.0)

_PARAM_GLOBAL = np.clip(
    _parametric_eff(N2D, _INC2D_f, b_mean, x0_mean, w_mean),
    0.0, 1.0).astype(np.float32)

BOOT_GLOBAL = (COV * _EMP_GLOBAL + (1.0 - COV) * _PARAM_GLOBAL).astype(np.float32)

# =============================================================================
# Gaussian quantile for CI_LEVEL
# =============================================================================

_z = float(_scipy_norm.ppf(0.5 + CI_LEVEL / 200.0))
print(f"CI_LEVEL={CI_LEVEL}%  →  z = {_z:.4f}")

# =============================================================================
# Blended CI bands — quadrature, lo/hi separately
#
#   half_emp_lo  = nominal - lo_emp            (empirical lower half-width)
#   half_emp_hi  = hi_emp  - nominal           (empirical upper half-width)
#   half_jac     = z * sigma_jac               (Jacobian half-width, symmetric)
#
#   lo_total = nominal - sqrt( COV^2 * half_emp_lo^2 + (1-COV)^2 * half_jac^2 )
#   hi_total = nominal + sqrt( COV^2 * half_emp_hi^2 + (1-COV)^2 * half_jac^2 )
# =============================================================================

_nom = BOOT_GLOBAL                          # (N_PC, N_INC)
_cov = COV                                  # (N_PC, N_INC)
_half_jac = (_z * _SIGMA_JAC).astype(np.float32)

_half_emp_lo = np.maximum(_nom - _EMP_LO, 0.0)
_half_emp_hi = np.maximum(_EMP_HI - _nom, 0.0)

BOOT_LO = (_nom - np.sqrt(
    _cov**2 * _half_emp_lo**2 + (1.0 - _cov)**2 * _half_jac**2
)).astype(np.float32)

BOOT_HI = (_nom + np.sqrt(
    _cov**2 * _half_emp_hi**2 + (1.0 - _cov)**2 * _half_jac**2
)).astype(np.float32)

np.clip(BOOT_LO, 0.0, 1.0, out=BOOT_LO)
np.clip(BOOT_HI, 0.0, 1.0, out=BOOT_HI)

# Median surface: keep empirical median in data region, nominal globally elsewhere
_emp_med_list = Parallel(n_jobs=N_JOBS)(
    delayed(_eval_interp)(bi, _PC2D_flat, _shape_ext)
    for bi in boot_interps
)
_EMP_MED_CUBE = np.stack(_emp_med_list, axis=0)
del _emp_med_list
_EMP_MED = np.nanmedian(_EMP_MED_CUBE, axis=0).astype(np.float32)
del _EMP_MED_CUBE

BOOT_MED  = (COV * _EMP_MED + (1.0 - COV) * _PARAM_GLOBAL).astype(np.float32)
BOOT_MEAN = BOOT_GLOBAL.copy()   # global blended surface serves as mean

del _EMP_MED, _EMP_LO, _EMP_HI, _half_emp_lo, _half_emp_hi, _half_jac

print(f"Stored surfaces of shape {BOOT_MED.shape}, "
      f"total size: {(BOOT_MED.nbytes * 4) / 1e6:.1f} MB")

# =============================================================================
# Grid scaling constants
# =============================================================================

_LOG_PC_LO, _LOG_PC_HI = COMMON_LOG_PC[0], COMMON_LOG_PC[-1]
_INC_LO,    _INC_HI    = COMMON_INC[0],    COMMON_INC[-1]
_SCALE_PC  = (N_COMMON_PC  - 1) / (_LOG_PC_HI - _LOG_PC_LO)
_SCALE_INC = (N_COMMON_INC - 1) / (_INC_HI    - _INC_LO)

# =============================================================================
# Fast CI lookup via map_coordinates  (UNCHANGED)
# =============================================================================

def trigger_eff_ci(N_gamma, inc_deg, reducer="global"):
    N_gamma = np.asarray(N_gamma, dtype=float)
    inc_deg = np.asarray(inc_deg, dtype=float)
    N_b, inc_b = np.broadcast_arrays(N_gamma, inc_deg)

    _N_raw   = N_b.ravel()
    _inc_raw = inc_b.ravel()

    _oob_pc_lo  = _N_raw < 10**_LOG_PC_LO
    _oob_pc_hi  = _N_raw > 10**_LOG_PC_HI
    _oob_inc_lo = _inc_raw < _INC_LO
    _oob_inc_hi = _inc_raw > _INC_HI
    _oob_pc_lo_nontrivial = _oob_pc_lo & (_N_raw > 0)   # N=0 is valid, not OOB
    _n_oob_warn = int((_oob_pc_lo_nontrivial | _oob_pc_hi |
                    _oob_inc_lo | _oob_inc_hi).sum())
    if _n_oob_warn > 0:
        warnings.warn(
            f"trigger_eff_ci: {_n_oob_warn} query point(s) outside grid "
            f"[logN: {_LOG_PC_LO:.2f}–{_LOG_PC_HI:.2f}, "
            f"inc: {_INC_LO:.2f}–{_INC_HI:.2f}]. Clamped to boundary.",
            stacklevel=2,
        )

    lp = np.log10(np.clip(_N_raw, 10**_LOG_PC_LO, 10**_LOG_PC_HI))
    ii = np.clip(_inc_raw, _INC_LO, _INC_HI)
    fx = (lp - _LOG_PC_LO) * _SCALE_PC
    fy = (ii - _INC_LO)    * _SCALE_INC
    coords = np.vstack([fx, fy])

    if reducer == "median":
        cen = map_coordinates(BOOT_MED,    coords, order=1, mode='nearest')
    elif reducer == "mean":
        cen = map_coordinates(BOOT_MEAN,   coords, order=1, mode='nearest')
    elif reducer == "global":
        cen = map_coordinates(BOOT_GLOBAL, coords, order=1, mode='nearest')
    else:
        raise ValueError(f"Unknown reducer: {reducer}")

    lo = map_coordinates(BOOT_LO, coords, order=1, mode='nearest')
    hi = map_coordinates(BOOT_HI, coords, order=1, mode='nearest')

    np.clip(cen, 0.0, 1.0, out=cen)
    np.clip(lo,  0.0, 1.0, out=lo)
    np.clip(hi,  0.0, 1.0, out=hi)

    return (cen.reshape(N_b.shape),
            lo.reshape(N_b.shape),
            hi.reshape(N_b.shape))


# =============================================================================
# Sanity check
# =============================================================================

N_test   = np.logspace(2, 7, 10)
inc_test = np.full_like(N_test, 4)
nom, lo, hi = trigger_eff_ci(N_test, inc_test)
for n_, e_, l_, h_ in zip(N_test, nom, lo, hi):
    print(f"N={n_:8.2e}, inc={inc_test[0]:.0f},  "
          f"nom={e_:.3e}  CI{CI_LEVEL:g}%=[{l_:.3e}, {h_:.3e}]")

# =============================================================================
# Save surfaces to disk
# =============================================================================

_SAVE_PATH = "/scratch/general/vast/u1520754/data_Muon_Trinity/trigger_eff_surfaces.npz"
np.savez_compressed(
    _SAVE_PATH,
    BOOT_MED=BOOT_MED,
    BOOT_LO=BOOT_LO,
    BOOT_HI=BOOT_HI,
    BOOT_GLOBAL=BOOT_GLOBAL,
    SIGMA_JAC=_SIGMA_JAC,
    COMMON_LOG_PC=COMMON_LOG_PC,
    COMMON_INC=COMMON_INC,
    CI_LEVEL=np.array([CI_LEVEL]),
    N_BOOT=np.array([N_BOOT]),
    B_MODEL=np.array([B_MODEL]),
    X0_MODEL=np.array([X0_MODEL]),
    W_MODEL=np.array([W_MODEL]),
)
import os as _os
print(f"\nSurfaces saved to '{_SAVE_PATH}' "
      f"({_os.path.getsize(_SAVE_PATH) / 1e6:.1f} MB).")



# Empirical efficiency cube on (zenith, azimuth, H, E) — OPTIMIZED version
# =============================================================================
# Key optimizations vs original:
#   1. ffill_along_H        — fully vectorized, no Python loops
#   2. efficiency_grid_H_E  — bin-membership computed once per (E,H) pair via
#                             np.digitize; tasks chunked by zenith so joblib
#                             overhead is ~N_Z_BINS jobs instead of N_tasks
#   3. _process_zen_slice   — phi/trig computed once per zenith, reused across
#                             all (H, E) bins in that slice; single batched
#                             trigger_eff_ci call per (H, E) bin
#   4. extrapolate_via_slant — same chunking strategy (one job per zen slice)
#   5. build_E_slant_pc_table — inner loop replaced with np.digitize + groupby
# =============================================================================

import os
import numpy as np
from joblib import Parallel, delayed
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt

# ---------- Config ----------
BOOT_SEED   = 123
N_BOOT_EFF  = N_BOOTSTRAP
MIN_PTS_EFF = 5
N_E_BINS    = 20 #40
N_H_BINS    = 10 #20
N_Z_BINS    = 10 #30
N_AZ_BINS   = 5


N_JOBS      = os.cpu_count()

Telescope_Zenith = 90
ENERGY_MIN = 1e1
ENERGY_MAX = 1e6
ZEN_MIN    = np.nanmin(Zenith)
ZEN_MAX    = np.nanmax(Zenith)

AZIMUTH_RANGE = 10 #ZEN_MAX - ZEN_MIN
AZ_MIN = 0.0
AZ_MAX = AZIMUTH_RANGE  # degrees

H_MIN = np.nanmin(H) #5_000
H_MAX = np.nanmax(H) #50_000

OBS_HEIGHT_M   = 2_944.0
EARTH_RADIUS_M = 6_371_000.0

N_PARAMS   = 1
PARAM_ORDER = ["dummy"]
COV_FULL   = np.zeros((N_PARAMS, N_PARAMS))

# =============================================================================
# Helpers
# =============================================================================

def _phi_inc_grid(zen_deg, tele_zen_deg, phi_max_deg, n_phi):
    """Returns phi (n_phi,) in degrees AND radians, and theta_inc (n_phi,)."""
    theta_z  = np.radians(zen_deg)
    psi      = np.radians(tele_zen_deg)
    phi_rad  = np.linspace(0.0, np.radians(phi_max_deg), n_phi)
    phi_deg  = np.degrees(phi_rad)
    cos_t    = (np.cos(psi) * np.cos(theta_z)
                + np.sin(psi) * np.sin(theta_z) * np.cos(phi_rad))
    cos_t    = np.clip(cos_t, -1.0, 1.0)
    theta_inc = np.degrees(np.arccos(cos_t))   # (n_phi,)
    return phi_rad, phi_deg, theta_inc


def _phi_avg_eff(pc_sel, theta_inc, phi, reducer="mean",
                 n_pc_bins="fd", n_boot=100, rng=None):
    """
    Average trigger efficiency over phi for an array of photon counts,
    propagating BOTH model-side CI (from trigger_eff_ci) AND
    statistical sampling CI (from bootstrap over showers).

    Parameters
    ----------
    pc_sel    : (Nsh,) photon counts for showers in this bin
    theta_inc : (n_phi,) incidence angles in degrees
    phi       : (n_phi,) phi values in radians
    reducer   : "mean" | "histogram"
    n_pc_bins : int or "fd"  (only used for histogram reducer)
    n_boot    : int          number of bootstrap resamples
    rng       : np.random.Generator or None

    Returns
    -------
    eff, eff_lo, eff_hi : scalars  (nominal, 16th, 84th percentile of total CI)
    """
    npts  = len(pc_sel)
    n_phi = len(phi)
    denom = phi[-1] - phi[0]
    if rng is None:
        rng = np.random.default_rng()

    # -----------------------------------------------------------
    # Compute per-shower phi-averaged efficiency (med, lo, hi)
    # Shape: (Nsh,)
    # -----------------------------------------------------------
    pc_b  = pc_sel[:, None]                              # (Nsh, 1)
    inc_b = np.broadcast_to(theta_inc, (npts, n_phi))    # (Nsh, n_phi)
    med, lo, hi = trigger_eff_ci(pc_b, inc_b)            # each (Nsh, n_phi)
    med_phi = np.trapezoid(med, phi, axis=-1) / denom    # (Nsh,)
    lo_phi  = np.trapezoid(lo,  phi, axis=-1) / denom    # (Nsh,)
    hi_phi  = np.trapezoid(hi,  phi, axis=-1) / denom    # (Nsh,)

    if reducer == "mean":
        return _combine_uncertainty(med_phi, lo_phi, hi_phi, n_boot, rng)

    elif reducer == "histogram":
        # Build pc histogram, evaluate efficiency at bin centers, weight by PDF
        pc_pos = pc_sel[pc_sel > 0]
        if len(pc_pos) < 2 or pc_pos.min() == pc_pos.max():
            # Fall back to the per-shower mean path
            return _combine_uncertainty(med_phi, lo_phi, hi_phi, n_boot, rng)

        if n_pc_bins == "fd":
            log_pc = np.log10(pc_pos)
            edges_log = np.histogram_bin_edges(log_pc, bins='fd')
            edges = 10**edges_log
            count, _ = np.histogram(pc_pos, bins=edges)
            nb = len(count)
        else:
            edges = np.logspace(np.log10(pc_pos.min()),
                                np.log10(pc_pos.max()), n_pc_bins + 1)
            count, _ = np.histogram(pc_sel, bins=edges)
            nb = n_pc_bins

        dp     = np.diff(edges)
        total  = count.sum()
        if total == 0:
            return _combine_uncertainty(med_phi, lo_phi, hi_phi, n_boot, rng)
        pdf    = count / (total * dp)
        pc_cen = np.sqrt(edges[:-1] * edges[1:])

        pc_b  = pc_cen[:, None]
        inc_b = np.broadcast_to(theta_inc, (nb, n_phi))
        med_h, lo_h, hi_h = trigger_eff_ci(pc_b, inc_b)
        med_h_phi = np.trapezoid(med_h, phi, axis=-1) / denom    # (nb,)
        lo_h_phi  = np.trapezoid(lo_h,  phi, axis=-1) / denom
        hi_h_phi  = np.trapezoid(hi_h,  phi, axis=-1) / denom

        weights = pdf * dp                                       # (nb,)
        eff_nom = float(np.sum(weights * med_h_phi))

        # ---- model-side CI: weighted average of lo/hi over pc bins ----
        eff_model_lo = float(np.sum(weights * lo_h_phi))
        eff_model_hi = float(np.sum(weights * hi_h_phi))

        # ---- sampling CI: bootstrap the bin counts (Poisson / multinomial) ----
        if n_boot > 0 and total > 0:
            # Multinomial resample of the histogram counts
            boot_counts = rng.multinomial(int(total), count / total,
                                          size=n_boot)            # (n_boot, nb)
            boot_pdf    = boot_counts / (total * dp[None, :])     # (n_boot, nb)
            boot_w      = boot_pdf * dp[None, :]                  # (n_boot, nb)
            boot_eff    = boot_w @ med_h_phi                      # (n_boot,)
            stat_lo, stat_hi = np.percentile(boot_eff, [16, 84])
        else:
            stat_lo = stat_hi = eff_nom

        # Combine in quadrature
        sig_stat_lo = max(eff_nom - stat_lo, 0.0)
        sig_stat_hi = max(stat_hi - eff_nom, 0.0)
        sig_mod_lo  = max(eff_nom - eff_model_lo, 0.0)
        sig_mod_hi  = max(eff_model_hi - eff_nom, 0.0)
        eff_lo = eff_nom - np.hypot(sig_stat_lo, sig_mod_lo)
        eff_hi = eff_nom + np.hypot(sig_stat_hi, sig_mod_hi)
        return eff_nom, max(eff_lo, 0.0), min(eff_hi, 1.0)

    else:
        raise ValueError(f"Unknown reducer: {reducer!r}")


def _combine_uncertainty(med_phi, lo_phi, hi_phi, n_boot, rng):
    """
    Helper: combine model-side CI (lo/hi per shower) with bootstrap
    sampling CI over showers, returning (nominal, lo, hi).
    """
    npts = len(med_phi)
    eff_nom = float(np.mean(med_phi))

    # Model-side averaged CI
    eff_model_lo = float(np.mean(lo_phi))
    eff_model_hi = float(np.mean(hi_phi))

    # Bootstrap sampling CI
    if npts >= 2 and n_boot > 0:
        idx = rng.integers(0, npts, size=(n_boot, npts))
        boot_means = med_phi[idx].mean(axis=1)
        stat_lo, stat_hi = np.percentile(boot_means, [16, 84])
    else:
        stat_lo = stat_hi = eff_nom

    sig_stat_lo = max(eff_nom - stat_lo, 0.0)
    sig_stat_hi = max(stat_hi - eff_nom, 0.0)
    sig_mod_lo  = max(eff_nom - eff_model_lo, 0.0)
    sig_mod_hi  = max(eff_model_hi - eff_nom, 0.0)

    eff_lo = eff_nom - np.hypot(sig_stat_lo, sig_mod_lo)
    eff_hi = eff_nom + np.hypot(sig_stat_hi, sig_mod_hi)
    return eff_nom, max(eff_lo, 0.0), min(eff_hi, 1.0)

# =============================================================================
# 1. Forward-fill along H — fully vectorized
# =============================================================================

def ffill_along_H(EFF, *extras):
    """
    Forward-fill NaNs along the H axis (axis=2) for arrays of shape
    (n_zen, n_az, n_H, n_E[, ...]).
    """
    n_zen, n_az, n_H, n_E = EFF.shape[:4]
    trailing = EFF.shape[4:]

    def _ffill_4d(arr):
        # arr shape (n_zen, n_az, n_H, n_E)
        out = arr.copy()
        finite = np.isfinite(out)
        h_idx = np.where(finite, np.arange(n_H)[None, None, :, None], 0)
        np.maximum.accumulate(h_idx, axis=2, out=h_idx)
        z_idx = np.arange(n_zen)[:, None, None, None]
        a_idx = np.arange(n_az)[None, :, None, None]
        e_idx = np.arange(n_E)[None, None, None, :]
        out   = out[z_idx, a_idx, h_idx, e_idx]
        first_finite = np.argmax(finite, axis=2)   # (n_zen, n_az, n_E)
        all_nan      = ~finite.any(axis=2)
        h_range      = np.arange(n_H)[None, None, :, None]
        before       = h_range < first_finite[:, :, None, :]
        out[before | all_nan[:, :, None, :]] = np.nan
        return out

    def _ffill_5d(arr):
        # arr shape (n_zen, n_az, n_H, n_E, N_PARAMS)
        n_p = arr.shape[4]
        return np.stack([_ffill_4d(arr[..., p]) for p in range(n_p)], axis=-1)

    out_eff = _ffill_4d(EFF)
    out_extras = []
    for arr in extras:
        if arr.ndim == 4:
            out_extras.append(_ffill_4d(arr))
        elif arr.ndim == 5:
            out_extras.append(_ffill_5d(arr))
        else:
            raise ValueError(f"ffill_along_H: unsupported ndim={arr.ndim}")

    return (out_eff, *out_extras) if out_extras else out_eff


# =============================================================================
# 2. Per-zenith-slice worker  (used by efficiency_grid_H_E)
# =============================================================================

def _process_zen_slice_HE(k, zen_deg, bin_tasks,
                           tele_zen_deg, phi_max_deg, n_phi,
                           reducer, n_boot, rng_seed):
    if not bin_tasks:
        return []

    rng = np.random.default_rng(rng_seed)
    phi_rad, phi_deg, theta_inc = _phi_inc_grid(
        zen_deg, tele_zen_deg, phi_max_deg, n_phi
    )

    results = []
    for (j, i, pc_sel) in bin_tasks:

        # --- evaluate trigger_eff_ci once for all phi ---
        # pc_sel: (N_showers,), theta_inc: (n_phi,)
        pc_b   = pc_sel[:, None]                          # (N, n_phi)
        inc_b  = theta_inc[None, :]                       # (1, n_phi)
        med_raw, lo_raw, hi_raw = trigger_eff_ci(pc_b, inc_b, reducer=reducer)
        # med_raw: (N, n_phi)

        # --- per-phi efficiency = mean over showers ---
        med_phi = med_raw.mean(axis=0)   # (n_phi,)
        lo_phi  = lo_raw.mean(axis=0)    # (n_phi,)
        hi_phi  = hi_raw.mean(axis=0)    # (n_phi,)

        # --- bootstrap uncertainty per phi bin ---
        if n_boot > 1 and len(pc_sel) >= 2:
            idx = rng.integers(0, len(pc_sel), size=(n_boot, len(pc_sel)))
            boot = med_raw[idx].mean(axis=1)   # (n_boot, n_phi)
            stat_lo = np.percentile(boot, 16, axis=0)  # (n_phi,)
            stat_hi = np.percentile(boot, 84, axis=0)
        else:
            stat_lo = med_phi.copy()
            stat_hi = med_phi.copy()

        # Combine model + stat uncertainty in quadrature per phi
        sig_stat_lo = np.maximum(med_phi - stat_lo, 0.0)
        sig_stat_hi = np.maximum(stat_hi - med_phi, 0.0)
        sig_mod_lo  = np.maximum(med_phi - lo_phi,  0.0)
        sig_mod_hi  = np.maximum(hi_phi  - med_phi, 0.0)
        eff_lo = med_phi - np.hypot(sig_stat_lo, sig_mod_lo)
        eff_hi = med_phi + np.hypot(sig_stat_hi, sig_mod_hi)

        results.append((
            k, j, i,
            med_phi,                        # (n_phi,)
            np.clip(eff_lo, 0.0, 1.0),
            np.clip(eff_hi, 0.0, 1.0),
        ))

    return results

# =============================================================================
# 3. Empirical efficiency cube  (parallel over zenith slices)
# =============================================================================

def efficiency_grid_H_E(
    E_edges, H_edges, Zen_grid, Az_grid,
    Energy_GeV, H, photon_count,
    tele_zen_deg=90.0, phi_max_deg=3.0, n_phi=50,
    min_pts=5, reducer='mean',
    n_boot=200, rng=None, n_jobs=1,
):
    n_zen = len(Zen_grid)
    n_az  = len(Az_grid)
    n_H   = len(H_edges) - 1
    n_E   = len(E_edges) - 1
    E_cen = np.sqrt(E_edges[:-1] * E_edges[1:])
    H_cen = 0.5 * (H_edges[:-1] + H_edges[1:])

    EFF    = np.full((n_zen, n_az, n_H, n_E), np.nan)
    EFF_LO = np.full((n_zen, n_az, n_H, n_E), np.nan)
    EFF_HI = np.full((n_zen, n_az, n_H, n_E), np.nan)
    COUNTS = np.zeros((n_zen, n_az, n_H, n_E), dtype=int)
    GRAD   = np.zeros((n_zen, n_az, n_H, n_E, N_PARAMS))
    pc_in_bin = [[[None] * n_E for _ in range(n_H)] for _ in range(n_zen)]

    # ------------------------------------------------------------------
    # Bin membership computed once with np.digitize (vectorized)
    # ------------------------------------------------------------------
    E_idx = np.digitize(Energy_GeV, E_edges) - 1   # 0-based; -1 or n_E = outside
    H_idx = np.digitize(H,          H_edges) - 1

    valid = (E_idx >= 0) & (E_idx < n_E) & (H_idx >= 0) & (H_idx < n_H)
    E_idx_v = E_idx[valid]
    H_idx_v = H_idx[valid]
    pc_v    = photon_count[valid]

    # Count hits per (H, E) bin  — same for every zenith
    flat_idx = H_idx_v * n_E + E_idx_v                # (N_valid,)
    counts_HE = np.bincount(flat_idx, minlength=n_H * n_E).reshape(n_H, n_E)

    # Build per-(j, i) pc arrays once
    _order  = np.argsort(flat_idx)
    _fi_s   = flat_idx[_order]
    _pc_s   = pc_v[_order]
    _splits = np.searchsorted(_fi_s, np.arange(n_H * n_E + 1))
    pc_HE   = {}                                       # {(j,i): pc_array}
    for flat in range(n_H * n_E):
        j, i = divmod(flat, n_E)
        sl   = slice(_splits[flat], _splits[flat + 1])
        cnt  = _splits[flat + 1] - _splits[flat]
        COUNTS[:, :, j, i] = cnt          # broadcast over both zen and az axes
        if cnt >= min_pts:
            pc_HE[(j, i)] = _pc_s[sl]
            for k in range(n_zen):
                pc_in_bin[k][j][i] = _pc_s[sl]

    if not pc_HE:
        return (EFF, EFF_LO, EFF_HI, COUNTS, E_cen, H_cen, pc_in_bin, GRAD)

    # ------------------------------------------------------------------
    # Build per-zenith task lists and dispatch one job per zenith slice
    # ------------------------------------------------------------------
    zen_tasks = []
    for k, zen_deg in enumerate(Zen_grid):
        bin_tasks = [(j, i, arr) for (j, i), arr in pc_HE.items()]
        zen_tasks.append((k, zen_deg, bin_tasks))

    n_total_bins = len(Zen_grid) * len(pc_HE)
    print(f"efficiency_grid_H_E: {n_total_bins} bins "
          f"({len(Zen_grid)} zen × {len(pc_HE)} HE), {n_jobs} CPUs")

    #slice_results = Parallel(n_jobs=n_jobs, verbose=5)(
    #    delayed(_process_zen_slice_HE)(k, zen_deg, bin_tasks,
    #                                   tele_zen_deg, phi_max_deg, n_phi, reducer)
    #    for k, zen_deg, bin_tasks in zen_tasks
    #)
    # Extract a deterministic seed from rng
    rng_seed = int(rng.integers(0, 2**31 - 1)) if rng is not None else None

    slice_results = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_process_zen_slice_HE)(k, zen_deg, bin_tasks,
                                   tele_zen_deg, phi_max_deg, n_phi, reducer,
                                   n_boot, rng_seed)
    for k, zen_deg, bin_tasks in zen_tasks
    )

    for slice_out in slice_results:
        for k, j, i, med_phi, lo_phi, hi_phi in slice_out:
            # med_phi: (n_phi,) → maps to azimuth axis
            for a in range(n_az):
                EFF[k, a, j, i]    = med_phi[a]
                EFF_LO[k, a, j, i] = lo_phi[a]
                EFF_HI[k, a, j, i] = hi_phi[a]
            COUNTS[k, :, j, i] = COUNTS[k, 0, j, i]  # same for all az

    return (EFF, EFF_LO, EFF_HI, COUNTS, E_cen, H_cen, pc_in_bin, GRAD)


# =============================================================================
# 4. Build (E, Slant) photon-count lookup — vectorized inner loop
# =============================================================================

def build_E_slant_pc_table(E_edges, H_emit, Energy_GeV, photon_count,
                           Zenith_deg, n_slant_bins):
    E_edges      = np.asarray(E_edges,      dtype=float)
    H_emit       = np.asarray(H_emit,       dtype=float)
    Energy_GeV   = np.asarray(Energy_GeV,   dtype=float)
    photon_count = np.asarray(photon_count, dtype=float)
    Zenith_rad   = np.radians(np.asarray(Zenith_deg, dtype=float))

    Slant_km = slant_depth(
        H_emit + EARTH_RADIUS_M,
        OBS_HEIGHT_M + EARTH_RADIUS_M,
        Zenith_rad,
    ) / 1000.0

    S_edges = np.linspace(np.nanmin(Slant_km), np.nanmax(Slant_km), n_slant_bins + 1)
    S_cen   = 0.5 * (S_edges[:-1] + S_edges[1:])
    n_E     = len(E_edges) - 1

    # Bin membership — vectorized
    finite_mask = np.isfinite(photon_count) & (photon_count > 0)
    E_idx = np.digitize(Energy_GeV, E_edges) - 1       # 0-based
    S_idx = np.digitize(Slant_km,   S_edges) - 1

    valid = (finite_mask
             & (E_idx >= 0) & (E_idx < n_E)
             & (S_idx >= 0) & (S_idx < n_slant_bins))
    E_idx_v  = E_idx[valid]
    S_idx_v  = S_idx[valid]
    pc_v     = photon_count[valid]

    flat_idx = E_idx_v * n_slant_bins + S_idx_v
    order    = np.argsort(flat_idx)
    fi_s     = flat_idx[order]
    pc_s     = pc_v[order]
    splits   = np.searchsorted(fi_s, np.arange(n_E * n_slant_bins + 1))

    pc_ES = [[None] * n_slant_bins for _ in range(n_E)]
    for flat in range(n_E * n_slant_bins):
        i, s = divmod(flat, n_slant_bins)
        sl = slice(splits[flat], splits[flat + 1])
        arr = pc_s[sl]
        pc_ES[i][s] = arr if len(arr) > 0 else np.empty(0)

    return S_edges, S_cen, pc_ES


# =============================================================================
# 5. Slant-based extrapolation — parallel over zenith slices
# =============================================================================

def _process_zen_slice_slant(k, zen, j_tasks, tele_zen_deg, phi_max_deg, n_phi,
                              reducer, n_boot, rng_seed):
    if not j_tasks:
        return []
    rng = np.random.default_rng(rng_seed + 1000 + k) if rng_seed is not None else np.random.default_rng()
    phi_rad, phi_deg, theta_inc = _phi_inc_grid(zen, tele_zen_deg, phi_max_deg, n_phi)

    out = []
    for j, i, pc in j_tasks:
        if len(pc) == 0:
            out.append((k, j, i, np.full(n_phi, np.nan),
                                  np.full(n_phi, np.nan),
                                  np.full(n_phi, np.nan)))
            continue

        pc_b  = pc[:, None]
        inc_b = theta_inc[None, :]
        med_raw, lo_raw, hi_raw = trigger_eff_ci(pc_b, inc_b, reducer=reducer)

        med_phi = med_raw.mean(axis=0)
        lo_phi  = lo_raw.mean(axis=0)
        hi_phi  = hi_raw.mean(axis=0)

        if n_boot > 1 and len(pc) >= 2:
            idx  = rng.integers(0, len(pc), size=(n_boot, len(pc)))
            boot = med_raw[idx].mean(axis=1)
            stat_lo = np.percentile(boot, 16, axis=0)
            stat_hi = np.percentile(boot, 84, axis=0)
        else:
            stat_lo = med_phi.copy()
            stat_hi = med_phi.copy()

        sig_stat_lo = np.maximum(med_phi - stat_lo, 0.0)
        sig_stat_hi = np.maximum(stat_hi - med_phi, 0.0)
        sig_mod_lo  = np.maximum(med_phi - lo_phi,  0.0)
        sig_mod_hi  = np.maximum(hi_phi  - med_phi, 0.0)
        eff_lo = np.clip(med_phi - np.hypot(sig_stat_lo, sig_mod_lo), 0.0, 1.0)
        eff_hi = np.clip(med_phi + np.hypot(sig_stat_hi, sig_mod_hi), 0.0, 1.0)

        out.append((k, j, i, med_phi, eff_lo, eff_hi))
    return out


def extrapolate_via_slant(EFF, EFF_LO, EFF_HI,
                          E_cen, H_cen, zenith_deg, S_edges, pc_ES,
                          tele_zen_deg=90.0, phi_max_deg=5.0, n_phi=20,
                          min_pts=5, reducer='mean', n_boot=200,
                          rng=None, n_jobs=N_JOBS):
    n_zen, n_az, n_H, n_E = EFF.shape   # <-- was (n_zen, n_H, n_E)

    zen_tasks = []
    for k in range(n_zen):
        zen = zenith_deg[k]
        j_tasks = []
        for j in range(n_H):
            slant_km = slant_depth(
                H_cen[j] + EARTH_RADIUS_M,
                OBS_HEIGHT_M + EARTH_RADIUS_M,
                np.deg2rad(zen),
            ) / 1000.0
            si = int(np.clip(np.searchsorted(S_edges, slant_km) - 1,
                             0, len(S_edges) - 2))
            for i in range(n_E):
                # Check if ANY azimuth slice is still NaN
                if np.all(np.isfinite(EFF[k, :, j, i])):
                    continue
                pc = pc_ES[i][si]
                if pc is None or len(pc) < min_pts:
                    continue
                j_tasks.append((j, i, pc))
        if j_tasks:
            zen_tasks.append((k, zen, j_tasks))

    if not zen_tasks:
        return EFF, EFF_LO, EFF_HI

    n_total = sum(len(t[2]) for t in zen_tasks)
    print(f"extrapolate_via_slant: {n_total} bins "
          f"across {len(zen_tasks)} zenith slices, {n_jobs} CPUs")

    rng_seed = int(rng.integers(0, 2**31 - 1)) if rng is not None else None
    slice_results = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_process_zen_slice_slant)(k, zen, j_tasks,
                                          tele_zen_deg, phi_max_deg, n_phi, reducer,
                                          n_boot, rng_seed)
        for k, zen, j_tasks in zen_tasks
    )

    for slice_out in slice_results:
        for k, j, i, med_phi, lo_phi, hi_phi in slice_out:
            if np.any(np.isfinite(med_phi)):
                EFF[k, :, j, i]    = med_phi    # (n_phi,) → az axis
                EFF_LO[k, :, j, i] = lo_phi
                EFF_HI[k, :, j, i] = hi_phi

    return EFF, EFF_LO, EFF_HI


# =============================================================================
# 6. Build efficiency interpolator (logic unchanged, calls optimized helpers)
# =============================================================================

def build_efficiency_interpolator(EFF, EFF_LO, EFF_HI,
                                  E_cen, H_cen, zenith_deg, S_edges, pc_ES,
                                  tele_zen_deg=90.0, phi_max_deg=5.0, n_phi=20,
                                  min_pts=3, reducer='mean', n_boot=200,
                                  rng=None,
                                  H_extrap_range=(0.0, 120_000.0),
                                  n_H_extra=4,
                                  zen_extrap_range=None,
                                  n_zen_extra=2,
                                  n_jobs=N_JOBS):
    H_cen      = np.asarray(H_cen,      dtype=float)
    zenith_deg = np.asarray(zenith_deg, dtype=float)
    n_az       = EFF.shape[1]            # <-- read from data

    # --- H / zenith extension (unchanged logic) ---
    if n_H_extra > 0:
        H_lo      = np.linspace(H_extrap_range[0], H_cen[0],   n_H_extra + 1)[:-1]
        H_hi      = np.linspace(H_cen[-1], H_extrap_range[1],  n_H_extra + 1)[1:]
        H_cen_ext = np.concatenate([H_lo, H_cen, H_hi])
    else:
        H_cen_ext = H_cen.copy()

    if zen_extrap_range is not None and n_zen_extra > 0:
        Z_lo    = np.linspace(zen_extrap_range[0], zenith_deg[0],  n_zen_extra + 1)[:-1]
        Z_hi    = np.linspace(zenith_deg[-1], zen_extrap_range[1], n_zen_extra + 1)[1:]
        zen_ext = np.concatenate([Z_lo, zenith_deg, Z_hi])
    else:
        zen_ext = zenith_deg.copy()

    z_off     = int(np.searchsorted(zen_ext, zenith_deg[0]))
    h_off     = int(np.searchsorted(H_cen_ext, H_cen[0]))
    n_zen_ext = len(zen_ext)
    n_H_ext   = len(H_cen_ext)
    n_E       = EFF.shape[3]             # <-- was shape[2]

    def _embed(arr_old):
        # arr_old shape: (n_zen, n_az, n_H, n_E)
        out = np.full((n_zen_ext, n_az, n_H_ext, n_E), np.nan)
        out[z_off:z_off + len(zenith_deg),
            :,
            h_off:h_off + len(H_cen),
            :] = arr_old
        return out

    EFF    = _embed(EFF)
    EFF_LO = _embed(EFF_LO)
    EFF_HI = _embed(EFF_HI)

    # fill_mask must be computed AFTER embed (on the extended arrays)
    fill_mask = ~np.isfinite(EFF)        # shape (n_zen_ext, n_az, n_H_ext, n_E)

    H_cen      = H_cen_ext
    zenith_deg = zen_ext

    # --- EDT fill (unchanged logic, now operates on 4-D arrays) ---
    EFF_filled    = EFF.copy()
    EFF_LO_filled = EFF_LO.copy()
    EFF_HI_filled = EFF_HI.copy()

    nan_mask = ~np.isfinite(EFF_filled)
    if np.any(nan_mask):
        _, idx = distance_transform_edt(nan_mask, return_indices=True)
        EFF_filled[nan_mask]    = EFF_filled[tuple(idx)][nan_mask]
        EFF_LO_filled[nan_mask] = EFF_LO_filled[tuple(idx)][nan_mask]
        EFF_HI_filled[nan_mask] = EFF_HI_filled[tuple(idx)][nan_mask]

    original_nan = ~np.isfinite(EFF)
    edge_mask = np.zeros_like(EFF_filled, dtype=bool)
    edge_mask[:, :, :2, :]  = True      # <-- H axis is now axis 2
    edge_mask[:2, :, :, :]  = True
    edge_mask[-2:, :, :, :] = True
    reopen = edge_mask & original_nan
    EFF_filled[reopen]    = np.nan
    EFF_LO_filled[reopen] = np.nan
    EFF_HI_filled[reopen] = np.nan

    EFF_filled, EFF_LO_filled, EFF_HI_filled = extrapolate_via_slant(
        EFF_filled, EFF_LO_filled, EFF_HI_filled,
        E_cen, H_cen, zenith_deg, S_edges, pc_ES,
        tele_zen_deg=tele_zen_deg, phi_max_deg=phi_max_deg, n_phi=n_phi,
        min_pts=min_pts, reducer=reducer, n_boot=n_boot, rng=rng,
        n_jobs=n_jobs,
    )

    remaining = ~np.isfinite(EFF_filled)
    if np.any(remaining):
        _, idx2 = distance_transform_edt(remaining, return_indices=True)
        EFF_filled[remaining]    = EFF_filled[tuple(idx2)][remaining]
        EFF_LO_filled[remaining] = EFF_LO_filled[tuple(idx2)][remaining]
        EFF_HI_filled[remaining] = EFF_HI_filled[tuple(idx2)][remaining]

    EFF_filled[~np.isfinite(EFF_filled)]       = 0.0
    EFF_LO_filled[~np.isfinite(EFF_LO_filled)] = 0.0
    EFF_HI_filled[~np.isfinite(EFF_HI_filled)] = 0.0

    log10_E = np.log10(E_cen)
    # Axis order: (zen, az, H, log10E) — consistent throughout
    rgi_eff  = RegularGridInterpolator(
        (zenith_deg, Az_grid, H_cen, log10_E), EFF_filled,
        method='linear', bounds_error=False, fill_value=None)
    rgi_lo   = RegularGridInterpolator(
        (zenith_deg, Az_grid, H_cen, log10_E), EFF_LO_filled,
        method='linear', bounds_error=False, fill_value=None)
    rgi_hi   = RegularGridInterpolator(
        (zenith_deg, Az_grid, H_cen, log10_E), EFF_HI_filled,
        method='linear', bounds_error=False, fill_value=None)
    rgi_mask = RegularGridInterpolator(
        (zenith_deg, Az_grid, H_cen, log10_E),   # <-- az axis added
        fill_mask.astype(float),
        method='nearest', bounds_error=False, fill_value=1.0)

    def _make_pts(E_GeV, H_m, zenith_deg_q, az_deg_q):
        E_b, H_b, Z_b, A_b = np.broadcast_arrays(
            np.asarray(E_GeV,        float),
            np.asarray(H_m,          float),
            np.asarray(zenith_deg_q, float),
            np.asarray(az_deg_q,     float),
        )
        return np.column_stack([
            Z_b.ravel(), A_b.ravel(), H_b.ravel(), np.log10(E_b.ravel())
        ]), E_b.shape

    def eff_func(E_GeV, H_m, zenith_deg_q, az_deg_q):
        pts, shp = _make_pts(E_GeV, H_m, zenith_deg_q, az_deg_q)
        return np.clip(rgi_eff(pts), 0.0, 1.0).reshape(shp)

    def eff_sigma_func(E_GeV, H_m, zenith_deg_q, az_deg_q):
        pts, shp = _make_pts(E_GeV, H_m, zenith_deg_q, az_deg_q)
        med = np.clip(rgi_eff(pts), 0.0, 1.0).reshape(shp)
        lo  = np.clip(rgi_lo(pts),  0.0, 1.0).reshape(shp)
        hi  = np.clip(rgi_hi(pts),  0.0, 1.0).reshape(shp)
        return med, lo, hi

    def eff_fill_mask_func(E_GeV, H_m, zenith_deg_q, az_deg_q):  # <-- az added
        pts, shp = _make_pts(E_GeV, H_m, zenith_deg_q, az_deg_q)
        return (rgi_mask(pts).reshape(shp) >= 0.5)

    extra = {
        "EFF_filled":    EFF_filled,
        "EFF_LO_filled": EFF_LO_filled,
        "EFF_HI_filled": EFF_HI_filled,
        "fill_mask":     fill_mask,
        "H_cen_ext":     H_cen,
        "zen_ext":       zenith_deg,
    }
    return eff_func, eff_sigma_func, eff_fill_mask_func, extra

# =============================================================================
# Driver
# =============================================================================

rng     = np.random.default_rng(BOOT_SEED)
E_edges = np.logspace(np.log10(ENERGY_MIN), np.log10(ENERGY_MAX), N_E_BINS + 1) #equal_count_bins(Energy_GeV, N_E_BINS) #np.logspace(np.log10(ENERGY_MIN), np.log10(ENERGY_MAX), N_E_BINS + 1)
H_edges = np.logspace(np.log10(H_MIN), np.log10(H_MAX), N_H_BINS + 1) #equal_count_bins(H, N_H_BINS)
Zen_grid = np.linspace(ZEN_MIN, ZEN_MAX, N_Z_BINS)
Az_grid = np.linspace(AZ_MIN, AZ_MAX, N_AZ_BINS)
REDUCER = "median" # mean, global
photon_mask = (photon_count >= 0)&(R==0) #np.where(photon_count >= 0) 
(EFF, EFF_LO, EFF_HI, COUNTS, E_cen, H_cen, pc_in_bin, GRAD) = efficiency_grid_H_E(
    E_edges, H_edges, Zen_grid, Az_grid,
    Energy_GeV[photon_mask], H[photon_mask], photon_count[photon_mask],
    tele_zen_deg=Telescope_Zenith, phi_max_deg=AZIMUTH_RANGE, n_phi=N_AZ_BINS,
    min_pts=MIN_PTS_EFF, reducer=REDUCER,
    n_boot=N_BOOT_EFF, rng=rng, n_jobs=N_JOBS,
)
filled_mask = ~np.isnan(EFF)  # before ffill
print(f"Filled bins: {filled_mask.sum()} / {EFF.size}  ({100*filled_mask.mean():.1f}%)")
sparse_mask = COUNTS < MIN_PTS_EFF
EFF[sparse_mask] = np.nan  # exclude from interpolation

EFF, EFF_LO, EFF_HI, GRAD = ffill_along_H(EFF, EFF_LO, EFF_HI, GRAD)

S_edges, S_cen, pc_ES = build_E_slant_pc_table(
    E_edges, H, Energy_GeV, photon_count, Zenith,
    n_slant_bins=N_H_BINS,
)

eff_func, eff_sigma_func, eff_fill_mask_func, eff_extra = build_efficiency_interpolator(
    EFF, EFF_LO, EFF_HI, E_cen, H_cen, Zen_grid, S_edges, pc_ES,
    tele_zen_deg=Telescope_Zenith, phi_max_deg=AZIMUTH_RANGE, n_phi=N_AZ_BINS,
    min_pts=MIN_PTS_EFF, reducer=REDUCER,
    n_boot=N_BOOT_EFF, rng=rng, n_jobs=N_JOBS,
)


# Effective area from digitized trigger-rate data
# =========================================================================

# --- Digitized data ---
R_data = np.array([
    0.08165819247995529, 0.3077888370767397, 0.48994972995969954, 0.684673244738047, 
    0.9108038893348314, 1.0929647822177913, 1.2939696079438325,1.4949750087498417, 
    2.0917084625202556, 2.2738690678632314, 2.908291537479744,5.100502125323036, 
    5.301507526129045
])
R_data = R_data - R_data[1]

trigger_data = np.array([
    0.013083803640247065, 0.014237288135593221,0.014336158192090397, 0.013709980161850063, 
    0.012721279596878311, 0.012424669427386787, 0.011864406779661017, 0.011172316384180791, 
    0.010150658127951758, 0.009392655367231639, 0.008634651097874185, 0.0033615819209039553, 
    0.001252353043205994
])

trigger_rate_max = np.array([
    0.01858757062146893, 0.016016949152542374, 0.015918079096045198, 0.015160074826687743, 
    0.014138418079096047, 0.013841807909604521, 0.013775894035727291, 0.013380413809738591, 
    0.011040488636426335, 0.011139358692923509, 0.010776836158192091, 0.003855932203389831, 
    0.0018126156909317622
])
trigger_rate_err = trigger_rate_max - trigger_data

trigger_rate_err = abs(max(trigger_data)*trigger_rate_err - trigger_data*max(trigger_rate_err))/(max(trigger_data))**2 #+ 1e-10
trigger_data = trigger_data/max(trigger_data)

trigger_rate_err = trigger_rate_err[1:]
trigger_data = trigger_data[1:]
R_data = R_data[1:]
print(len(trigger_rate_err), len(trigger_data))
# --- Weighted linear fit ---
fit_weights = 1.0 / trigger_rate_err**2

(p, cov_fit) = np.polyfit(
    R_data,
    trigger_data,
    deg=1,
    w=np.sqrt(fit_weights),
    cov=True
)

m_fit, b_fit = p
m_err = np.sqrt(cov_fit[0, 0])
b_err = np.sqrt(cov_fit[1, 1])

print(f"Slope     m = {m_fit:.6e} ± {m_err:.6e}")
print(f"Intercept b = {b_fit:.6e} ± {b_err:.6e}")

# --- Weighted R² ---
y_fit        = m_fit * R_data + b_fit
weighted_mean = np.average(trigger_data, weights=fit_weights)
ss_res = np.sum(fit_weights * (trigger_data - y_fit)**2)
ss_tot = np.sum(fit_weights * (trigger_data - weighted_mean)**2)
r_squared = 1 - ss_res / ss_tot
print(f"Weighted R² = {r_squared:.5f}")

# --- Central effective area ---
r_grid = np.linspace(0, 6, 10000)
trigger_rate_r      = np.clip(m_fit * r_grid + b_fit, 0, None)
trigger_rate_r_norm = trigger_rate_r / trigger_rate_r.max()
effective_area_central = 2 * np.pi * np.trapezoid(
    trigger_rate_r_norm * r_grid, x=r_grid
)
print(f"Effective Area (central) = {effective_area_central * 1e4:.4f} cm²")

# --- Monte Carlo uncertainty propagation ---
N_mc = 1000
samples = np.random.multivariate_normal(
    mean=[m_fit, b_fit],
    cov=cov_fit,
    size=N_mc
)

area_samples = []
for m_i, b_i in samples:
    trigger_i = np.clip(m_i * r_grid + b_i, 0, None)
    if trigger_i.max() <= 0:
        continue
    trigger_i_norm = trigger_i 
    area_i = 2 * np.pi * np.trapezoid(trigger_i_norm * r_grid, x=r_grid)
    area_samples.append(area_i)

area_samples = np.array(area_samples) * 1e4 #convert to cm^2

effective_area_mean = np.mean(area_samples) 
effective_area_std  = np.std(area_samples) 
raw_area = np.pi*500**2

print(f"Effective Area (MC mean) = {effective_area_mean * 1e4:.2f} ± "
      f"{effective_area_std * 1e4:.2f} cm²")
print("Raw Area", raw_area, " cm²")
# --- Fit plot ---
plt.figure(figsize=(8, 6))
plt.errorbar(R_data, trigger_data, yerr=trigger_rate_err,
             fmt='o', color="black", capsize=4, label='Data')
plt.plot(r_grid, (m_fit * r_grid + b_fit), '--',
         label=(
                f'm = {m_fit:.3e} ± {m_err:.1e}\n'
                f'b = {b_fit:.3e} ± {b_err:.1e}\n'
                f'R² = {r_squared:.4f}'))


plt.xlabel("Distance from Detector Center (m)")
plt.ylim(0,1.01)
plt.ylabel("Detector Trigger Efficiency")
plt.grid(True)
plt.legend()
plt.savefig("/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/plot/eff_r.png",transparent=True)

plt.show()

# Muon Flux Calculation
import pickle
import numpy as np
import matplotlib.pyplot as plt

def convert_altitude_to_depth(altitude_m, zenith_deg):
    altitude_m  = np.asarray(altitude_m, dtype=float)
    altitude_km = altitude_m * 1e-3
    altitude_cm = altitude_m * 100.0

    zenith_deg  = np.asarray(zenith_deg, dtype=float)
    zenith_deg  = np.where(zenith_deg >= 90.0, 89.9, zenith_deg)

    depth = np.empty_like(altitude_m)

    m1 = (altitude_km >   0.0) & (altitude_km <   4.0)
    m2 = (altitude_km >=  4.0) & (altitude_km <  10.0)
    m3 = (altitude_km >= 10.0) & (altitude_km <  40.0)
    m4 = (altitude_km >= 40.0) & (altitude_km <= 100.0)
    m5 =  altitude_km > 100.0

    depth[m1] = -186.555305 + 1222.6562 * np.exp(-altitude_cm[m1] / 994186.38)
    depth[m2] =  -94.919    + 1144.9069 * np.exp(-altitude_cm[m2] / 878153.55)
    depth[m3] =    0.61289  + 1305.5948 * np.exp(-altitude_cm[m3] / 636143.04)
    depth[m4] =    0.0      +  540.1778 * np.exp(-altitude_cm[m4] / 772170.16)
    depth[m5] =    0.01128292 - altitude_cm[m5] / 1e9

    depth /= np.cos(np.deg2rad(zenith_deg))

    return depth


def eff_sigma_func_X(E_mesh, H_mesh, Z_mesh, A_mesh):
    X_mesh = convert_altitude_to_depth(H_mesh, Z_mesh)
    return eff_sigma_func(E_mesh, H_mesh, Z_mesh, A_mesh)

# =========================================================================
# 2. Load the precomputed MCEq grid
# =========================================================================
"""with open('/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/data/'
          'mceq_grid_80_90deg_Nzen100_no_decay_no_E_loss.pkl', 'rb') as f:"""
with open('/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/data/'
          'mceq_grid_80_90deg_Nzen100_no_decay_no_E_loss.pkl', 'rb') as f:          
    grid_data = pickle.load(f)

zen_centers_all = grid_data['zen_centers']
cos_edges       = grid_data['cos_edges']
sin_edges       = np.sqrt(1-cos_edges**2)
theta_edges     = np.arccos(cos_edges)
int_sin2_edges  = (theta_edges - sin_edges*cos_edges)/2
E_GeV           = grid_data['E_GeV']

min_E           = ENERGY_MIN #grid_data['min_E'] #
max_E           = 1e6#grid_data['max_E'] #
E_mask          = (E_GeV>=min_E)&(E_GeV<=max_E)
E_GeV           = E_GeV[E_mask]

zen_min_use, zen_max_use = 85, 90.0
tele_azimuth = 270
n_azimuth = 5

zen_mask = (zen_centers_all >= zen_min_use) & (zen_centers_all <= zen_max_use)
if not np.any(zen_mask):
    raise ValueError(
        f"No zenith bins in [{zen_min_use}, {zen_max_use}] deg. "
        f"Available centers: {zen_centers_all}"
    )

zen_centers = zen_centers_all[zen_mask]

min_az, max_az = tele_azimuth - AZIMUTH_RANGE, tele_azimuth + AZIMUTH_RANGE
phi = np.linspace(0, AZIMUTH_RANGE, n_azimuth)

dphi = np.radians(max_az - min_az)

selected_entries = [e for e, keep in zip(grid_data['per_zenith'], zen_mask) if keep]

# Azimuth offsets (same for every zenith). 
# phi already defined above as: phi = np.linspace(0, AZIMUTH_RANGE, n_azimuth)
n_A = len(phi)

flux_slices = []   # will hold one (n_H, n_E) -> (n_H, n_E, n_A) array per zenith
zen_list    = []   # keep track of the zenith value for each slice (for bookkeeping)
X_list      = []
for entry in selected_entries:
    theta  = entry['theta']      # zenith (deg) for this slice
    H_grid = entry['H_grid']     # emission heights (cm)
    d_CDF  = entry['d_CDF']      # muon flux, shape (n_H, n_E)  <-- adjust if different
    X_grid = entry['X_grid']
    # --- the flux for this zenith: shape (n_H, n_E) ---
    d_X = np.gradient(X_grid)[:,None]
    #d_CDF = d_CDF/d_X
    flux_HE = np.asarray(d_CDF)[:, E_mask]

    # --- broadcast over azimuth: (n_H, n_E) -> (n_H, n_E, n_A) ---
    flux_HEA = np.broadcast_to(flux_HE[:, :, None], flux_HE.shape + (n_A,))

    flux_slices.append(flux_HEA)
    zen_list.append(theta)
    X_list.append(X_grid)
X_arr = np.array(X_list)
# Stack along a new leading zenith axis -> (Z, H, E, A)
flux_tensor = np.stack(flux_slices, axis=0)

zen_arr = np.array(zen_list)   # (Z,)

print("flux_tensor shape (Z, H, E, A):", flux_tensor.shape)
print("zenith values:", zen_arr)

# Common axes (built once)
H_m      = np.asarray(selected_entries[0]['H_grid']) / 100.0   # (n_H,)
E_axis   = np.asarray(E_GeV)                                   # (n_E,)
phi_axis = np.asarray(phi)                                     # (n_A,)

# 3D meshes for ONE zenith slice — built once, reused every iteration
H_mesh, E_mesh, A_mesh = np.meshgrid(H_m, E_axis, phi_axis, indexing='ij')
# each shape (n_H, n_E, n_A) ≈ 18 MB

Z_mesh = np.empty_like(H_mesh)   # reused buffer for the scalar zenith

eff_slices = []
eff_slices_lo = []
eff_slices_hi = []
for theta in zen_arr:
    print(theta)
    Z_mesh.fill(theta)                                   # constant zenith
    #eff_HEA = eff_func(E_mesh, H_mesh, Z_mesh, A_mesh)   # (n_H, n_E, n_A)
    eff_HEA, eff_lo, eff_hi = eff_sigma_func(E_mesh, H_mesh, Z_mesh, A_mesh)
    eff_slices.append(eff_HEA.astype(np.float32))        # downcast to save RAM if neccesary
    eff_slices_lo.append(eff_lo.astype(np.float32))        # downcast to save RAM if neccesary
    eff_slices_hi.append(eff_hi.astype(np.float32))        # downcast to save RAM if neccesary

eff_tensor = np.stack(eff_slices, axis=0)   # (Z, H, E, A)
eff_tensor_lo = np.stack(eff_slices_lo, axis=0)
eff_tensor_hi = np.stack(eff_slices_hi, axis=0)
print("eff_tensor shape:", eff_tensor.shape, "dtype:", eff_tensor.dtype)
assert eff_tensor.shape == flux_tensor.shape

flux_eff_tensor = eff_tensor*flux_tensor
flux_eff_tensor_lo = eff_tensor_lo*flux_tensor
flux_eff_tensor_hi = eff_tensor_hi*flux_tensor
# =========================================================================
# Extrapolate flux grid to 91.8 degrees by copying the 90-degree slice
# =========================================================================

# Find the index of the 90-degree zenith slice
idx_90 = np.argmin(np.abs(zen_arr - 90.0))
print(f"Using zenith slice at {zen_arr[idx_90]:.2f}° as the basis for extrapolation")

# Define the new zenith angles to extrapolate to (above 90 degrees)
zen_extrap = np.degrees(np.arccos(np.linspace(np.cos(np.radians(91.8)), np.cos(np.radians(90)), 10))) #np.array([91.8, 91.0, 90.9])  # can add more values here if needed

# --- Flux extrapolation ---
# Copy the 90° flux slice for each new zenith angle
flux_extrap_slices = []
X_extrap = []
for z_new in zen_extrap:
    # Copy the 90° slice: shape (H, E, A)
    flux_extrap_slices.append(flux_tensor[idx_90].copy())
    X_extrap.append(X_arr[idx_90].copy())
flux_extrap_tensor = np.stack(flux_extrap_slices, axis=0)  # (Z_new, H, E, A)
print(f"flux_extrap_tensor shape: {flux_extrap_tensor.shape}")
# --- Efficiency extrapolation ---
# Recompute efficiency at the new zenith angles
eff_extrap_slices    = []
eff_extrap_slices_lo = []
eff_extrap_slices_hi = []

for z_new in zen_extrap:
    print(f"Computing efficiency at zenith = {z_new}°")
    Z_mesh_extrap = np.full_like(H_mesh, z_new)
    eff_HEA, eff_lo, eff_hi = eff_sigma_func(E_mesh, H_mesh, Z_mesh_extrap, A_mesh)
    eff_extrap_slices.append(eff_HEA.astype(np.float32))
    eff_extrap_slices_lo.append(eff_lo.astype(np.float32))
    eff_extrap_slices_hi.append(eff_hi.astype(np.float32))

eff_extrap_tensor    = np.stack(eff_extrap_slices,    axis=0)  # (Z_new, H, E, A)
eff_extrap_tensor_lo = np.stack(eff_extrap_slices_lo, axis=0)
eff_extrap_tensor_hi = np.stack(eff_extrap_slices_hi, axis=0)
print(f"eff_extrap_tensor shape: {eff_extrap_tensor.shape}")

# --- Concatenate original + extrapolated tensors ---
zen_arr_full          = np.concatenate([zen_extrap, zen_arr         ])
X_arr_full            = np.concatenate([np.array(X_extrap),   X_arr])
flux_tensor_full      = np.concatenate([flux_extrap_tensor, flux_tensor],      axis=0)
eff_tensor_full       = np.concatenate([eff_extrap_tensor, eff_tensor       ],       axis=0)
eff_tensor_lo_full    = np.concatenate([eff_extrap_tensor_lo, eff_tensor_lo    ],    axis=0)
eff_tensor_hi_full    = np.concatenate([eff_extrap_tensor_hi, eff_tensor_hi    ],    axis=0)

flux_eff_tensor_full    = eff_tensor_full    * flux_tensor_full
flux_eff_tensor_lo_full = eff_tensor_lo_full * flux_tensor_full
flux_eff_tensor_hi_full = eff_tensor_hi_full * flux_tensor_full

print(f"\nFinal combined tensor shapes:")
print(f"  zen_arr_full:          {zen_arr_full}")
print(f"  flux_tensor_full:      {flux_tensor_full.shape}")
print(f"  eff_tensor_full:       {eff_tensor_full.shape}")
print(f"  flux_eff_tensor_full:  {flux_eff_tensor_full.shape}")


# Replace the original tensors with the full ones for downstream use
zen_arr          = zen_arr_full
X_arr            = X_arr_full
flux_tensor      = flux_tensor_full
eff_tensor       = eff_tensor_full
eff_tensor_lo    = eff_tensor_lo_full
eff_tensor_hi    = eff_tensor_hi_full
flux_eff_tensor    = flux_eff_tensor_full
flux_eff_tensor_lo = flux_eff_tensor_lo_full
flux_eff_tensor_hi = flux_eff_tensor_hi_full

print("\nDone. Downstream variables updated to include 91.8° extrapolation.")

def Muon_Rate_Calculation(flux_eff_tensor, effective_area_mean, E_axis, phi_axis,zen_arr):    
    R_zen_az = np.sum(np.trapezoid(flux_eff_tensor, x = E_axis, axis = 2), axis = 1)*effective_area_mean
    R_zen = 2*np.trapezoid(R_zen_az*np.cos(np.radians(phi_axis)), x = np.radians(phi_axis), axis = 1)

    zen_arr = zen_arr-Telescope_Zenith
    zen_arr_pos = zen_arr[zen_arr >=0]
    R_zen_pos = R_zen[zen_arr >=0]
    R_pos = np.trapezoid(R_zen_pos*np.cos(np.radians(zen_arr_pos))**2, x = -(np.radians(zen_arr_pos)), axis = 0) * 3600 *24

    zen_arr_neg = zen_arr[zen_arr < 0]
    R_zen_neg = R_zen[zen_arr < 0]
    R_neg = np.trapezoid(R_zen_neg*np.cos(np.radians(zen_arr_neg))**2, x = -(np.radians(zen_arr_neg)), axis = 0) * 3600 *24

    R = R_neg + R_pos
    return R, R_neg, R_pos


Rate = Muon_Rate_Calculation(flux_eff_tensor, effective_area_mean , E_axis, phi_axis, zen_arr)[0]  #effective_area_mean
Rate_raw = Muon_Rate_Calculation(flux_tensor_full, raw_area, E_axis, phi_axis, zen_arr)[0]
Rate_lo = Muon_Rate_Calculation(flux_eff_tensor_lo, effective_area_mean - effective_area_std, E_axis, phi_axis, zen_arr)[0]
Rate_hi = Muon_Rate_Calculation(flux_eff_tensor_hi, effective_area_mean + effective_area_std, E_axis, phi_axis, zen_arr)[0]
print("Raw rate", Rate_raw)
print(Rate/3, Rate_lo/3, Rate_hi/3, "muons per 8 hours")

def Muon_Rate_per_E_Calculation(flux_eff_tensor, effective_area_mean, phi_axis,zen_arr):    
    R_zen_az_E = np.sum(flux_eff_tensor, axis = 1)
    R_zen_E = 2*np.trapezoid(R_zen_az_E*np.cos(np.radians(phi_axis)), x = np.radians(phi_axis), axis = 2)
    zen_arr = zen_arr-Telescope_Zenith
    zen_arr_pos = zen_arr[zen_arr >=0][:, None]
    R_zen_pos = R_zen_E[zen_arr >=0]
    R_pos = np.trapezoid(R_zen_pos*np.cos(np.radians(zen_arr_pos))**2, x = -(np.radians(zen_arr_pos)), axis = 0) * 3600 *24

    zen_arr_neg = zen_arr[zen_arr < 0][:, None]
    R_zen_neg = R_zen_E[zen_arr < 0]
    R_neg = np.trapezoid(R_zen_neg*np.cos(np.radians(zen_arr_neg))**2, x = -(np.radians(zen_arr_neg)), axis = 0) * 3600 *24

    R_E = R_neg + R_pos
    return R_E*effective_area_mean
"""def Muon_Rate_per_E_Calculation(flux_eff_tensor, effective_area_mean, phi_axis, zen_arr, X_arr):
    # flux_eff_tensor: (Z, H, E, A)
    # X_arr:           (Z, H)
    n_Z, n_H, n_E, n_A = flux_eff_tensor.shape

    # Integrate over H using X as variable, per zenith slice -> (Z, E, A)
    R_ZEA = np.empty((n_Z, n_E, n_A))
    for iz in range(n_Z):
        X_slice     = X_arr[iz]           # (n_H,)
        flux_slice  = flux_eff_tensor[iz] # (n_H, n_E, n_A)

        sort_idx    = np.argsort(X_slice)
        X_sorted    = X_slice[sort_idx]
        flux_sorted = flux_slice[sort_idx]

        R_ZEA[iz]   = np.trapezoid(flux_sorted, x=X_sorted, axis=0)  # (n_E, n_A)

    # Integrate over azimuth with cos(phi) weighting: (Z, E, A) -> (Z, E)
    R_zen_E = 2 * np.trapezoid(
        R_ZEA * np.cos(np.radians(phi_axis)),
        x=np.radians(phi_axis),
        axis=2
    )

    # Zenith offset from telescope pointing
    zen_offset = zen_arr - Telescope_Zenith   # (Z,)

    zen_pos    = zen_offset[zen_offset >= 0][:, None]   # (Z+, 1)
    R_zen_pos  = R_zen_E[zen_offset >= 0]               # (Z+, E)
    R_pos      = np.trapezoid(
        R_zen_pos * np.cos(np.radians(zen_pos))**2,
        x=-np.radians(zen_pos),
        axis=0
    ) * 3600 * 24                                       # (E,)

    zen_neg    = zen_offset[zen_offset < 0][:, None]    # (Z-, 1)
    R_zen_neg  = R_zen_E[zen_offset < 0]                # (Z-, E)
    R_neg      = np.trapezoid(
        R_zen_neg * np.cos(np.radians(zen_neg))**2,
        x=-np.radians(zen_neg),
        axis=0
    ) * 3600 * 24                                       # (E,)

    R_E = (R_pos + R_neg) * effective_area_mean
    return R_E


R_E     = Muon_Rate_per_E_Calculation(flux_eff_tensor,    effective_area_mean,                  phi_axis, zen_arr, X_arr)
R_E_lo  = Muon_Rate_per_E_Calculation(flux_eff_tensor_lo, effective_area_mean - effective_area_std, phi_axis, zen_arr, X_arr)
R_E_hi  = Muon_Rate_per_E_Calculation(flux_eff_tensor_hi, effective_area_mean + effective_area_std, phi_axis, zen_arr, X_arr)
R_E_raw = Muon_Rate_per_E_Calculation(flux_tensor,        raw_area,                             phi_axis, zen_arr, X_arr)
R_eff   = R_E / R_E_raw"""


R_E = Muon_Rate_per_E_Calculation(flux_eff_tensor, effective_area_mean, phi_axis, zen_arr)
R_E_lo = Muon_Rate_per_E_Calculation(flux_eff_tensor_lo, effective_area_mean - effective_area_std, phi_axis, zen_arr)
R_E_hi = Muon_Rate_per_E_Calculation(flux_eff_tensor_hi, effective_area_mean + effective_area_std, phi_axis, zen_arr)
R_E_raw = Muon_Rate_per_E_Calculation(flux_tensor, raw_area, phi_axis, zen_arr)
R_eff = R_E/R_E_raw

fig, ax1 = plt.subplots(figsize=(7, 5))

# ── left axis : muon rate ──────────────────────────────────────────────────
color_rate = 'blue'

err_lo = R_E - R_E_lo  # lower error bar lengths (always positive)
err_lo = np.clip(err_lo, np.min(err_lo[err_lo>0]), None)
err_hi = R_E_hi - R_E   # upper error bar lengths (always positive)
E_data_mask = E_axis<= max(Energy_GeV)
ax1.errorbar(
    np.log10(E_axis[E_data_mask]), R_E[E_data_mask],
    yerr=[err_lo[E_data_mask], err_hi[E_data_mask]],
    fmt='--o',
    color=color_rate,
    ecolor=color_rate,
    elinewidth=1.5,
    capsize=3,
    markersize=4,
    linewidth=2,
    label='Detected Muon Rate',
)
#ax1.plot(np.log10(E_axis), R_E, color=color_rate, linewidth=1.5, linestyle='--', label='Extrapolate Muon Rate')
ax1.plot(np.log10(E_axis), R_E_raw,
         color='navy', linewidth=1.5, linestyle='--', label='Raw Muon Rate (no eff.)')

ax1.set_yscale('log')
ax1.set_xlabel(r'$\log_{10}(E\;/\;\mathrm{GeV})$')
ax1.set_ylabel(r'Muon Rate $(\frac{d^{2}N}{dtdE})$ [day$^{-1}$ GeV$^{-1}$]', color=color_rate, fontsize=17)
ax1.tick_params(axis='y', labelcolor=color_rate)
ax1.legend(loc='upper left', fontsize=10)

# ── right axis : efficiency ────────────────────────────────────────────────
color_eff = 'black'
ax2 = ax1.twinx()

# Efficiency error bars propagated from lo/hi rate ratios
R_eff_lo = R_E_lo / R_E_raw
R_eff_hi = R_E_hi / R_E_raw

eff_err_lo = R_eff - R_eff_lo
eff_err_lo = np.clip(eff_err_lo, np.min(eff_err_lo[eff_err_lo>0]), None)

eff_err_hi = R_eff_hi - R_eff

#ax2.plot(np.log10(E_axis), R_eff, color=color_eff, linewidth=1.5, linestyle='--', label='Raw Muon Rate (no eff.)')

ax2.errorbar(
    np.log10(E_axis[E_data_mask]), R_eff[E_data_mask],
    yerr=[eff_err_lo[E_data_mask], eff_err_hi[E_data_mask]],
    fmt='--s',
    color=color_eff,
    ecolor=color_eff,
    elinewidth=1.5,
    capsize=3,
    markersize=4,
    linewidth=2,
    label='Efficiency',
)

ax2.set_ylabel(r'Detector Trigger Efficiency $(\epsilon)$', color=color_eff, fontsize=17)
ax2.tick_params(axis='y', labelcolor=color_eff)
ax2.set_yscale("log")
ax2.legend(loc='upper right', fontsize=10)

# ── cosmetics ─────────────────────────────────────────────────────────────
ax1.grid(True, which='both', alpha=0.3)
ax1.text(1.8, 1e-12, 
        r"$\frac{dN}{dt}\ = \int_{E_{\min}}^{E_{\max}} \frac{d^{2}N(E)}{dtdE} \epsilon(E) \, dE = 0.3$"+" "+r"$ day^{-1}$",
        fontsize=15,
        color = color_rate
)
fig.tight_layout()
ax1.set_ylim(1e-13, 1e8)
ax2.set_ylim(1e-9, 1e-1)
plt.xlim(1.4, 6)
plt.savefig("/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/plot/Muon_rate.png",transparent=True)
plt.show()
