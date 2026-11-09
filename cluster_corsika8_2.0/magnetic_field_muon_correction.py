#Define Backtrack muon
import numpy as np
import shlex

# ===============================================================
# Constants
# ===============================================================
C          = 299792458.0
E_CHARGE   = 1.602176634e-19
MEV_TO_J   = 1.602176634e-13
M_MU_MEV   = 105.6583755
R_EARTH    = 6371e3
TAU_MU_S   = 2.1969811e-6


# ===============================================================
# Exact vectorized constant-field helix
# Unique name: unaffected by an older function named _helix.
# ===============================================================
def _muon_helix_vectorized(x0, u0, pc, charge, B, s):
    """
    x0, u0: (..., 3)
    pc: momentum*c [MeV]
    charge: charge in units of elementary charge
    B: shared magnetic field, shape (3,), [T]
    s: path length [m]
    """
    x0 = np.asarray(x0, dtype=float)
    u0 = np.asarray(u0, dtype=float)
    B = np.asarray(B, dtype=float)

    pc, charge, s = np.broadcast_arrays(
        np.asarray(pc, dtype=float),
        np.asarray(charge, dtype=float),
        np.asarray(s, dtype=float),
    )

    shape = np.broadcast_shapes(
        x0.shape[:-1], u0.shape[:-1], pc.shape
    )

    x0 = np.broadcast_to(x0, shape + (3,))
    u0 = np.broadcast_to(u0, shape + (3,))
    pc = np.broadcast_to(pc, shape)
    charge = np.broadcast_to(charge, shape)
    s = np.broadcast_to(s, shape)

    Bmag = np.linalg.norm(B)

    if Bmag == 0.0:
        return x0 + s[..., None] * u0, u0.copy()

    bhat = B / Bmag

    u_parallel = np.sum(
        u0 * bhat, axis=-1, keepdims=True
    ) * bhat

    u_perp = u0 - u_parallel
    w = np.cross(bhat, u_perp)

    angle = (
        -charge * E_CHARGE * Bmag * C / MEV_TO_J
    ) * s / pc

    # Stable at angle = 0.
    sinc = np.sinc(angle / np.pi)
    cosc = (
        0.5 * angle
        * np.sinc(angle / (2.0 * np.pi))**2
    )

    x = x0 + s[..., None] * (
        u_parallel
        + sinc[..., None] * u_perp
        + cosc[..., None] * w
    )

    u = (
        u_parallel
        + np.cos(angle)[..., None] * u_perp
        + np.sin(angle)[..., None] * w
    )

    u /= np.linalg.norm(u, axis=-1, keepdims=True)

    return x, u


# ===============================================================
# Shared numerical implementation for BOTH scalar and grid calls
# ===============================================================
def _muon_backtrack_core(
    E_tot_GeV,
    charge_sign,
    zenith_deg,
    azimuth_deg,
    h_prod_m,
    H_obs_m=3000.0,
    B_vec=(28.6e-6, 0.0, -46.4e-6),
    rest_mass_MeV=M_MU_MEV,
    lifetime_s=TAU_MU_S,
    pdg_abs=13,
    R_earth=R_EARTH,
    scan_step_m=500.0,
    max_path_m=2_000_000.0,
    path_tol_m=1e-6,
    max_scan_rotation_rad=0.02,
):
    """
    Constant-energy, uniform-field geometric backtracking.

    Angles describe FORWARD travel:
        u = (sin(zen)*cos(az), sin(zen)*sin(az), -cos(zen))

    Production altitude is spherical:
        norm(x_prod) = R_earth + h_prod_m

    Backtracking reverses velocity AND charge, keeping B unchanged.

    Scalar inputs return scalars and shape-(3,) vectors.
    Array inputs follow NumPy broadcasting.

    Non-muon usage:
      Pass ``B_vec=(0.0, 0.0, 0.0)`` to disable magnetic deflection.  With B = 0
      the trajectory is a straight line, so the injection position/direction are
      independent of ``rest_mass_MeV``; only beta/gamma/times/survival/PDG depend
      on it.  Set ``rest_mass_MeV`` to the primary's rest mass (0.0 for a photon),
      ``lifetime_s=None`` for a stable particle (survival = 1, decay = 0) and
      ``pdg_abs`` to the primary's |PDG|.  A neutral primary must still pass
      ``charge_sign = +/-1`` (the value is irrelevant while B = 0).

    Search limitations:
      - First sign-changing crossing detected by the scan.
      - Tangencies / multiple crossings within one interval can be missed.
      - No ground rejection or material interactions.
      - Equal detector/production heights return s = 0.
    """
    arrays = np.broadcast_arrays(
        *[
            np.asarray(v, dtype=float)
            for v in (
                E_tot_GeV,
                charge_sign,
                zenith_deg,
                azimuth_deg,
                h_prod_m,
                H_obs_m,
            )
        ]
    )

    shape = arrays[0].shape
    E, q, zen, azi, hp, ho = [
        a.ravel() for a in arrays
    ]

    B = np.asarray(B_vec, dtype=float)

    if B.shape != (3,) or not np.all(np.isfinite(B)):
        raise ValueError("B_vec must contain three finite numbers.")

    if E.size == 0:
        raise ValueError("Muon inputs cannot be empty.")

    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("All muon inputs must be finite.")

    if rest_mass_MeV < 0.0:
        raise ValueError("rest_mass_MeV must be >= 0.")

    if np.any(E <= rest_mass_MeV / 1000.0):
        raise ValueError("Total energy must exceed the particle rest energy.")

    if np.any((q != -1) & (q != 1)):
        raise ValueError("charge_sign must be -1 or +1.")

    if np.any((zen < 0) | (zen > 180)):
        raise ValueError("zenith_deg must be in [0, 180].")

    if np.any(hp < 0) or np.any(ho < 0):
        raise ValueError("Production/detector heights must be >= 0.")

    controls = (
        R_earth,
        scan_step_m,
        max_path_m,
        path_tol_m,
        max_scan_rotation_rad,
    )
    if any(not np.isfinite(v) or v <= 0 for v in controls):
        raise ValueError("Earth radius and search controls must be positive.")

    E_MeV = E * 1000.0
    pc = np.sqrt(
        np.maximum(E_MeV - rest_mass_MeV, 0.0) * (E_MeV + rest_mass_MeV)
    )

    th = np.radians(zen)
    ph = np.radians(azi)

    u_tel = np.stack(
        [
            np.sin(th) * np.cos(ph),
            np.sin(th) * np.sin(ph),
            -np.cos(th),
        ],
        axis=-1,
    )

    x_tel = np.zeros((E.size, 3))
    x_tel[:, 2] = R_earth + ho

    target_r = R_earth + hp

    def backward(s, indices):
        return _muon_helix_vectorized(
            x_tel[indices],
            -u_tel[indices],
            pc[indices],
            -q[indices],
            B,
            s,
        )

    def residual(s, indices):
        x, _ = backward(s, indices)
        return np.linalg.norm(x, axis=-1) - target_r[indices]

    kb = E_CHARGE * np.linalg.norm(B) * C / MEV_TO_J

    if kb == 0.0:
        step = np.full(E.size, scan_step_m)
    else:
        step = np.minimum(
            scan_step_m,
            max_scan_rotation_rad * pc / kb,
        )

    lo = np.zeros(E.size)
    hi = np.zeros(E.size)
    flo = ho - hp

    found = (ho == hp)
    active = ~found

    # -----------------------------------------------------------
    # Bracket each trajectory independently.
    # No Python loop over muons.
    # -----------------------------------------------------------
    while np.any(active):
        idx = np.flatnonzero(active)

        trial = np.minimum(
            lo[idx] + step[idx], max_path_m
        )

        if np.any(trial <= lo[idx]):
            raise RuntimeError(
                "Scan cannot advance at floating-point precision."
            )

        ft = residual(trial, idx)

        crossed = (
            (ft == 0.0)
            | (np.signbit(ft) != np.signbit(flo[idx]))
        )

        hit = idx[crossed]
        hi[hit] = trial[crossed]
        found[hit] = True
        active[hit] = False

        miss = idx[~crossed]
        lo[miss] = trial[~crossed]
        flo[miss] = ft[~crossed]
        active[miss] = trial[~crossed] < max_path_m

    # -----------------------------------------------------------
    # Bisection, including exact-root and roundoff handling.
    # -----------------------------------------------------------
    active = found & ((hi - lo) > path_tol_m)

    while np.any(active):
        idx = np.flatnonzero(active)
        mid = lo[idx] + 0.5 * (hi[idx] - lo[idx])

        stalled = (mid == lo[idx]) | (mid == hi[idx])
        if np.any(stalled):
            active[idx[stalled]] = False

        idx = idx[~stalled]
        mid = mid[~stalled]

        if idx.size == 0:
            continue

        fm = residual(mid, idx)
        exact = (fm == 0.0)

        root = idx[exact]
        lo[root] = mid[exact]
        hi[root] = mid[exact]
        flo[root] = 0.0

        remaining = idx[~exact]
        middle = mid[~exact]
        fmiddle = fm[~exact]

        left = (
            np.signbit(fmiddle)
            != np.signbit(flo[remaining])
        )

        hi[remaining[left]] = middle[left]
        lo[remaining[~left]] = middle[~left]
        flo[remaining[~left]] = fmiddle[~left]

        active[idx] = (
            (hi[idx] - lo[idx]) > path_tol_m
        )

    s = np.where(
        found, lo + 0.5 * (hi - lo), np.nan
    )

    x_prod = np.full((E.size, 3), np.nan)
    u_inj = np.full((E.size, 3), np.nan)

    position_error = np.full(E.size, np.nan)
    direction_error = np.full(E.size, np.nan)
    altitude_error = np.full(E.size, np.nan)

    idx = np.flatnonzero(found)

    if idx.size:
        xp, ub = backward(s[idx], idx)
        x_prod[idx] = xp
        u_inj[idx] = -ub

        xr, ur = _muon_helix_vectorized(
            x_prod[idx],
            u_inj[idx],
            pc[idx],
            q[idx],
            B,
            s[idx],
        )

        position_error[idx] = np.linalg.norm(
            xr - x_tel[idx], axis=-1
        )
        direction_error[idx] = np.linalg.norm(
            ur - u_tel[idx], axis=-1
        )
        altitude_error[idx] = (
            np.linalg.norm(xp, axis=-1) - target_r[idx]
        )

    theta = np.degrees(
        np.arccos(np.clip(-u_inj[:, 2], -1.0, 1.0))
    )

    phi = np.degrees(
        np.arctan2(u_inj[:, 1], u_inj[:, 0])
    ) % 360.0

    beta = pc / E_MeV
    if rest_mass_MeV > 0.0:
        gamma = E_MeV / rest_mass_MeV
    else:
        gamma = np.full_like(E_MeV, np.inf)

    flight_time = s / (beta * C)
    proper_time = s * rest_mass_MeV / (pc * C)
    if lifetime_s is not None and np.isfinite(lifetime_s) and lifetime_s > 0.0:
        exponent = proper_time / lifetime_s
    else:
        # Stable particle: no decay -> survival = 1, decay = 0.
        exponent = np.zeros_like(proper_time)

    result = dict(
        success=found,
        pdg=(-pdg_abs * q).astype(int),
        charge_sign=q,

        detector_zenith_deg=zen,
        detector_azimuth_deg=azi % 360.0,
        production_altitude_m=hp,

        injection_zenith_deg=theta,
        injection_azimuth_deg=phi,
        injection_x_m=x_prod[:, 0],
        injection_y_m=x_prod[:, 1],

        # Preserve YOUR existing command-height convention.
        injection_height_m=hp,

        # Distinct quantity: Cartesian z minus Earth radius.
        injection_z_minus_R_m=x_prod[:, 2] - R_earth,

        injection_total_energy_GeV=E,
        observation_level_m=ho,

        position_earth_m=x_prod,
        direction_unit=u_inj,
        path_length_m=s,

        beta=beta,
        gamma=gamma,
        flight_time_s=flight_time,
        proper_time_s=proper_time,
        survival_probability=np.exp(-exponent),
        decay_probability=-np.expm1(-exponent),

        forward_position_error_m=position_error,
        forward_direction_error=direction_error,
        production_altitude_error_m=altitude_error,
    )

    vector_keys = {"position_earth_m", "direction_unit"}

    for key, value in result.items():
        trailing = (3,) if key in vector_keys else ()
        value = value.reshape(shape + trailing)
        result[key] = (
            value.item() if value.ndim == 0 else value
        )

    return result


# ===============================================================
# Public scalar/broadcast interface
# ===============================================================
def back_track_muon(
    E_tot_GeV,
    charge_sign,
    zenith_deg,
    azimuth_deg,
    h_prod_m,
    **kwargs,
):
    return _muon_backtrack_core(
        E_tot_GeV,
        charge_sign,
        zenith_deg,
        azimuth_deg,
        h_prod_m,
        **kwargs,
    )


# Explicit alternative name.
back_track_muon_vectorized = back_track_muon


# ===============================================================
# Cartesian-product grid interface
# Calls the private core, NOT the public scalar function.
# ===============================================================
def back_track_muon_grid(
    E_tot_GeV,
    charge_sign,
    zenith_deg,
    azimuth_deg,
    h_prod_m,
    **kwargs,
):
    """
    Scalar-field output shape:
        (energy, charge, zenith, azimuth, production_height)

    Vector fields have an additional final dimension of 3.
    """
    inputs = [
        np.asarray(v, dtype=float)
        for v in (
            E_tot_GeV,
            charge_sign,
            zenith_deg,
            azimuth_deg,
            h_prod_m,
        )
    ]

    if any(a.ndim > 1 for a in inputs):
        raise ValueError("Grid inputs must be scalars or 1-D arrays.")

    if any(a.size == 0 for a in inputs):
        raise ValueError("Grid inputs cannot be empty.")

    E, Q, Z, A, H = np.meshgrid(
        *(np.atleast_1d(a) for a in inputs),
        indexing="ij",
        sparse=True,
    )

    return _muon_backtrack_core(
        E, Q, Z, A, H, **kwargs
    )


# ===============================================================
# Extract one result, preserving vector fields
# ===============================================================
def select_muon(res, idx):
    """
    Example:
        one = select_muon(res, (0, 0, 3, 5, 0))
    """
    shape = np.asarray(res["success"]).shape

    if not isinstance(idx, tuple):
        idx = (idx,)

    if len(idx) != len(shape):
        raise ValueError(
            f"Expected {len(shape)} indices for shape {shape}."
        )

    if not all(isinstance(i, (int, np.integer)) for i in idx):
        raise ValueError("Use integer indices, not slices.")

    selected = {}

    for key, value in res.items():
        element = np.asarray(value)[idx]
        selected[key] = (
            element.item()
            if element.ndim == 0
            else element.copy()
        )

    if not selected["success"]:
        raise ValueError(f"No successful crossing at index {idx}.")

    return selected


# ===============================================================
# Command generation
# ===============================================================
def c8_commands(
    res,
    out="backtrack_muon",
    extra="",
    executable="corsika-install/bin/muon_cherenkov_3d",
    print_commands=True,
    height_mode="spherical",
):
    """
    Return a flat list of commands.

    height_mode:
        "spherical": use h_prod_m (your existing convention).
        "z_minus_R": use computed Cartesian z - R.

    Choose the mode matching your executable's interpretation of
    --injection-height; these quantities are not interchangeable.

    --observation-level is added unless already present in extra.
    """
    ok = np.asarray(res["success"]).ravel()

    if not np.all(ok):
        raise ValueError(
            "Some trajectories have no crossing. "
            "Select successful elements before generating commands."
        )

    if height_mode not in ("spherical", "z_minus_R"):
        raise ValueError(
            "height_mode must be 'spherical' or 'z_minus_R'."
        )

    height_key = (
        "injection_height_m"
        if height_mode == "spherical"
        else "injection_z_minus_R_m"
    )

    def flat(key):
        return np.asarray(res[key]).ravel()

    pdg = flat("pdg")
    energy = flat("injection_total_energy_GeV")
    height = flat(height_key)
    x = flat("injection_x_m")
    y = flat("injection_y_m")
    zen = flat("injection_zenith_deg")
    azi = flat("injection_azimuth_deg")
    obs = flat("observation_level_m")

    extra = extra.strip()
    tokens = shlex.split(extra)

    has_observation = any(
        token == "--observation-level"
        or token.startswith("--observation-level=")
        for token in tokens
    )

    commands = []

    for i in range(ok.size):
        output = out if ok.size == 1 else f"{out}_{i}"

        cmd = (
            f"{shlex.quote(executable)} "
            f"--pdg {int(pdg[i])} "
            f"-E {energy[i]:.15g} "
            f"--injection-height {height[i]:.15g} "
            f"--injection-x {x[i]:.15g} "
            f"--injection-y {y[i]:.15g} "
            f"--injection-zenith {zen[i]:.15g} "
            f"--injection-azimuth {azi[i]:.15g} "
        )

        if not has_observation:
            cmd += f"--observation-level {obs[i]:.15g} "

        cmd += f"-f {shlex.quote(output)}"

        if extra:
            cmd += " " + extra

        commands.append(cmd)

    if print_commands:
        print("\n".join(commands))

    return commands


def c8_command(
    res,
    pdg=None,
    out="backtrack_muon",
    extra="",
    **kwargs,
):
    """
    Generate one command. PDG is inferred from the result.
    Optional pdg is checked for consistency.
    """
    if np.asarray(res["success"]).size != 1:
        raise ValueError(
            "Select one muon first using select_muon()."
        )

    actual_pdg = int(np.asarray(res["pdg"]).item())

    if pdg is not None and pdg != actual_pdg:
        raise ValueError(
            f"Requested PDG {pdg}, but result has PDG {actual_pdg}."
        )

    return c8_commands(
        res,
        out=out,
        extra=extra,
        **kwargs,
    )[0]