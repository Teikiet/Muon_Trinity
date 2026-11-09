"""
Landau distribution fitting worker module.
Must be a real .py file (not defined in the notebook) so that
ProcessPoolExecutor workers can pickle/import the functions.
"""
import numpy as np
from numba import njit, vectorize, float64
from scipy.optimize import minimize, approx_fprime, curve_fit

# ---------------------------------------------------------------------
# CERNLIB Landau PDF rational-approximation coefficients
# ---------------------------------------------------------------------
_p1 = (0.4259894875, -0.1249762550, 0.03984243700, -0.006298287635, 0.001511162253)
_q1 = (1.0, -0.3388260629, 0.09594393323, -0.01608042283, 0.003778942063)
_p2 = (0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411, 0.0001283617211)
_q2 = (1.0, 0.7428795082, 0.3153932961, 0.06694219548, 0.008790609714)
_p3 = (0.1788357263, 0.09202856865, 0.008682256834, 0.0006367807209, 0.0)
_q3 = (1.0, 0.7803488730, 0.3345163797, 0.07135755393, 0.008670274108)
_p4 = (0.9874054407, 118.6723273, 849.2794360, -743.7792444, 427.0262186)
_q4 = (1.0, 106.8615961, 337.6496214, 2016.712389, 1597.063511)
_p5 = (1.003675074, 167.5702434, 4789.711289, 21217.86767, -22324.94910)
_q5 = (1.0, 156.9424537, 3745.310488, 22759.60148, 12279.90862)
_p6 = (1.000827619, 664.9143136, 62972.92665, 475554.6998, -5743609.109)
_q6 = (1.0, 651.4101098, 56974.73333, 165917.4725, -2815759.939)
_a1 = (0.04166666667, -0.01996527778, 0.02709538966)
_a2 = (-1.845568670, -4.284640743)


@njit(cache=True, fastmath=True)
def _landau_pdf_scalar(x):
    if x < -5.5:
        u = np.exp(x + 1.0)
        if u < 1e-10:
            return 0.0
        ue = np.exp(-1.0 / u)
        us = np.sqrt(u)
        denom = _a1[0] + (_a1[1] + _a1[2] * u) * u
        return 0.3989422803 * (ue / us) * (1.0 + denom * u)
    elif x < -1.0:
        u = np.exp(-x - 1.0)
        num = _p1[0] + (_p1[1] + (_p1[2] + (_p1[3] + _p1[4] * x) * x) * x) * x
        den = _q1[0] + (_q1[1] + (_q1[2] + (_q1[3] + _q1[4] * x) * x) * x) * x
        return num / den * np.exp(-u) * np.sqrt(u)
    elif x < 1.0:
        num = _p2[0] + (_p2[1] + (_p2[2] + (_p2[3] + _p2[4] * x) * x) * x) * x
        den = _q2[0] + (_q2[1] + (_q2[2] + (_q2[3] + _q2[4] * x) * x) * x) * x
        return num / den
    elif x < 5.0:
        num = _p3[0] + (_p3[1] + (_p3[2] + _p3[3] * x) * x) * x
        den = _q3[0] + (_q3[1] + (_q3[2] + _q3[3] * x) * x) * x
        return num / den
    elif x < 12.0:
        u = 1.0 / x
        num = _p4[0] + (_p4[1] + (_p4[2] + (_p4[3] + _p4[4] * u) * u) * u) * u
        den = _q4[0] + (_q4[1] + (_q4[2] + (_q4[3] + _q4[4] * u) * u) * u) * u
        return num / den * (u * u)
    elif x < 50.0:
        u = 1.0 / x
        num = _p5[0] + (_p5[1] + (_p5[2] + (_p5[3] + _p5[4] * u) * u) * u) * u
        den = _q5[0] + (_q5[1] + (_q5[2] + (_q5[3] + _q5[4] * u) * u) * u) * u
        return num / den * (u * u)
    elif x < 300.0:
        u = 1.0 / x
        num = _p6[0] + (_p6[1] + (_p6[2] + (_p6[3] + _p6[4] * u) * u) * u) * u
        den = _q6[0] + (_q6[1] + (_q6[2] + (_q6[3] + _q6[4] * u) * u) * u) * u
        return num / den * (u * u)
    else:
        u = 1.0 / (x - x * np.log(x) / (x + 1.0))
        return 1.0 - _a2[0] * u - _a2[1] * u * u


@vectorize([float64(float64)], cache=True, fastmath=True)
def _landau_pdf_vec(x):
    return _landau_pdf_scalar(x)


def landau_pdf(x, mpv, eta):
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    if eta <= 0:
        return np.zeros_like(x)
    xi = (x - mpv) / eta
    return _landau_pdf_vec(xi) / eta


def landau_cdf_numeric(x, mpv, eta, x_lo=None, n=600):
    """Numeric CDF via cumulative trapezoid on a fixed grid (coarser
    grid than original for speed; still plenty accurate for Landau)."""
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    if x_lo is None:
        x_lo = mpv - 30 * eta
    x_hi = max(np.max(x) + 1e-6, x_lo + 1e-6)
    grid = np.linspace(x_lo, x_hi, n)
    pdf_grid = landau_pdf(grid, mpv, eta)
    cdf_grid = np.concatenate(([0.0], np.cumsum(0.5 * (pdf_grid[1:] + pdf_grid[:-1]) * np.diff(grid))))
    return np.interp(x, grid, cdf_grid)


def chi_square_reduced(observed_counts, expected_counts, n_params, min_count=5):
    """
    Proper reduced chi-square using raw (unnormalized) counts.
    observed_counts, expected_counts: raw counts per bin (NOT densities).
    Bins with expected_counts < min_count are excluded (low-count bins
    violate the Gaussian/chi2 approximation and blow up the statistic).
    """
    if expected_counts.size == 0 or np.max(expected_counts) <= 0:
        return np.nan
    m = expected_counts >= min_count
    ndof = int(np.sum(m)) - n_params
    if ndof < 1:
        return np.nan
    return np.sum((observed_counts[m] - expected_counts[m]) ** 2 / expected_counts[m]) / ndof


def censored_landau_nll(params, pe_mid, n_left, n_right, x_lo, x_hi):
    mpv, log_eta = params
    eta = np.exp(log_eta)
    if eta <= 0 or not np.isfinite(eta):
        return np.inf

    nll = 0.0
    if len(pe_mid) > 0:
        pdf_vals = landau_pdf(pe_mid, mpv, eta)
        pdf_vals = np.clip(pdf_vals, 1e-300, None)
        nll -= np.sum(np.log(pdf_vals))

    if n_left > 0 or n_right > 0:
        grid_lo = mpv - 30 * eta
        grid_hi = max(x_hi, mpv + 30 * eta)
        grid = np.linspace(grid_lo, grid_hi, 600)
        pdf_grid = landau_pdf(grid, mpv, eta)
        cdf_grid = np.concatenate(([0.0], np.cumsum(0.5 * (pdf_grid[1:] + pdf_grid[:-1]) * np.diff(grid))))
        denom = cdf_grid[-1] if cdf_grid[-1] > 0 else 1.0
        cdf_grid = cdf_grid / denom

        if n_left > 0:
            F_left = np.interp(x_lo, grid, cdf_grid)
            F_left = np.clip(F_left, 1e-300, 1.0)
            nll -= n_left * np.log(F_left)
        if n_right > 0:
            F_right = np.interp(x_hi, grid, cdf_grid)
            surv = np.clip(1.0 - F_right, 1e-300, 1.0)
            nll -= n_right * np.log(surv)
    return nll


def fit_censored_landau(pe_raw, pe_min, pe_sat):
    """Censored MLE for Landau(mpv, eta)."""
    left = pe_raw <= pe_min
    right = pe_raw >= pe_sat
    mid = ~left & ~right

    n_left, n_right = int(left.sum()), int(right.sum())
    pe_mid = pe_raw[mid]
    n_mid = len(pe_mid)

    out = dict(ok=False, mpv=np.nan, eta=np.nan,
               mpv_err=np.nan, eta_err=np.nan,
               n_left=n_left, n_right=n_right, n_mid=n_mid)
    if n_mid < 2:
        return out

    mpv0 = np.median(pe_mid)
    eta0 = max(0.5 * np.std(pe_mid), 1e-2)

    res = minimize(censored_landau_nll, x0=[mpv0, np.log(eta0)],
                   args=(pe_mid, n_left, n_right, pe_min, pe_sat),
                   method='Nelder-Mead',
                   options={'xatol': 1e-6, 'fatol': 1e-6, 'maxiter': 5000})
    if not res.success:
        return out

    mpv_fit, log_eta_fit = res.x
    eta_fit = np.exp(log_eta_fit)

    mpv_err = eta_err = np.nan
    try:
        eps = 1e-4
        f = lambda p: censored_landau_nll(p, pe_mid, n_left, n_right, pe_min, pe_sat)
        hess = np.zeros((2, 2))
        g0 = approx_fprime(res.x, f, eps)
        for k in range(2):
            xp = res.x.copy()
            xp[k] += eps
            hess[:, k] = (approx_fprime(xp, f, eps) - g0) / eps
        cov = np.linalg.inv(hess)
        if cov[0, 0] > 0:
            mpv_err = np.sqrt(cov[0, 0])
        if cov[1, 1] > 0:
            eta_err = eta_fit * np.sqrt(cov[1, 1])
    except Exception:
        pass

    out.update(ok=True, mpv=mpv_fit, eta=eta_fit,
               mpv_err=mpv_err, eta_err=eta_err)
    return out


def fit_curvefit_landau(pe_mid, lin_counts, lin_centers):
    """Least-squares Landau fit to the linear-domain histogram, built
    ONLY from uncensored data."""
    out = dict(ok=False, mpv=np.nan, eta=np.nan, mpv_err=np.nan, eta_err=np.nan)
    if len(pe_mid) < 10 or lin_counts.size == 0:
        return out

    mpv0 = lin_centers[np.argmax(lin_counts)] if lin_counts.max() > 0 else np.median(pe_mid)
    eta0 = max(0.5 * np.std(pe_mid), 1e-2)
    A0 = max(lin_counts.max(), 1e-9)

    def _model(x, mpv, eta, A):
        eta = abs(eta)
        if eta <= 0:
            return np.zeros_like(x)
        return A * landau_pdf(x, mpv, eta)

    try:
        popt, pcov = curve_fit(_model, lin_centers, lin_counts,
                                p0=[mpv0, eta0, A0], maxfev=20000)
        mpv_fit, eta_fit, _A = popt
        eta_fit = abs(eta_fit)
        if not np.isfinite(mpv_fit) or eta_fit <= 0:
            return out
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else [np.nan] * 3
        mpv_err, eta_err = perr[0], perr[1]

        out.update(ok=True, mpv=mpv_fit, eta=eta_fit,
                   mpv_err=mpv_err, eta_err=eta_err)
    except Exception:
        pass
    return out


def fit_one_cell(i, j, pe_raw, pe_min_valid, pe_saturation, n_lin_bins, n_disp_bins=40):
    """
    Do all the work for a single (pc_bin=i, inc_bin=j) cell:
    build histograms, run both fits, compute chi2, build display curves.
    Returns a plain-dict payload (picklable) with everything the parent
    process needs to populate grids/panel_data/masks.
    """
    n_raw = len(pe_raw)
    left_mask = pe_raw <= pe_min_valid
    right_mask = pe_raw >= pe_saturation
    mid_mask = ~left_mask & ~right_mask
    pe_mid = pe_raw[mid_mask]
    n_left, n_right, n_mid = int(left_mask.sum()), int(right_mask.sum()), len(pe_mid)

    if n_mid < 2:
        return dict(i=i, j=j, status='nodata', n_raw=n_raw,
                   n_left=n_left, n_right=n_right, n_mid=n_mid)

        # Shared linear-domain histogram (chi2 reference for both methods)
    lin_edges = np.linspace(0.0, pe_saturation, n_lin_bins + 1)
    lin_counts_chi2, _ = np.histogram(pe_mid, bins=lin_edges)          # raw counts
    lin_centers_chi2 = 0.5 * (lin_edges[1:] + lin_edges[:-1])
    lin_widths_chi2 = np.diff(lin_edges)

    # Display histogram (finer binning, density for plotting)
    disp_edges = np.linspace(0.0, pe_saturation, n_disp_bins + 1)
    disp_counts, _ = np.histogram(pe_mid, bins=disp_edges, density=True)
    disp_centers = 0.5 * (disp_edges[1:] + disp_edges[:-1])
    disp_widths = np.diff(disp_edges)

    # Method A: censored MLE
    mle = fit_censored_landau(pe_raw, pe_min_valid, pe_saturation)
    if mle['ok']:
        pdf_mle = landau_pdf(lin_centers_chi2, mle['mpv'], mle['eta'])
        exp_mle = pdf_mle * n_mid * lin_widths_chi2
        mle['chi2'] = chi_square_reduced(lin_counts_chi2, exp_mle, n_params=2)
    else:
        mle['chi2'] = np.nan

    # Method B: curve_fit on censored-out histogram
    cf = fit_curvefit_landau(pe_mid, lin_counts_chi2, lin_centers_chi2)  # see below
    if cf['ok']:
        pdf_cf = landau_pdf(lin_centers_chi2, cf['mpv'], cf['eta'])
        exp_cf = pdf_cf * n_mid * lin_widths_chi2
        cf['chi2'] = chi_square_reduced(lin_counts_chi2, exp_cf, n_params=3)
    else:
        cf['chi2'] = np.nan

    # Winner
    c_mle = mle['chi2'] if (mle['ok'] and np.isfinite(mle['chi2'])) else np.inf
    c_cf = cf['chi2'] if (cf['ok'] and np.isfinite(cf['chi2'])) else np.inf
    if np.isinf(c_mle) and np.isinf(c_cf):
        win = -1
    elif c_mle <= c_cf:
        win = 0
    else:
        win = 1

    # Display curves
    xs = np.linspace(max(pe_min_valid, 1e-6), pe_saturation, 400)
    curve_mle = landau_pdf(xs, mle['mpv'], mle['eta']) if mle['ok'] else np.array([])
    curve_cf = landau_pdf(xs, cf['mpv'], cf['eta']) if cf['ok'] else np.array([])

    return dict(
        i=i, j=j, status='ok', n_raw=n_raw,
        n_left=n_left, n_right=n_right, n_mid=n_mid,
        lin_x=disp_centers, lin_y=disp_counts, lin_w=disp_widths,
        curve_x=xs, curve_mle=curve_mle, curve_cf=curve_cf,
        mle=mle, cf=cf, winner=win,
    )