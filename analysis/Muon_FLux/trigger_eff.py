"""
trigger_eff.py  -  fast trigger-efficiency lookup  eff(incidence, photon_count)

Build once:
    from trigger_eff import build_efficiency_table
    build_efficiency_table("landau_gev_gauss_2d_fit_results_weighted_tested.npz",
                           "trigger_eff_table.npz")
Use:
    from trigger_eff import TriggerEfficiency
    teff = TriggerEfficiency("trigger_eff_table.npz")
    eff, err = teff(inc, pc, mode="hybrid")          # any shapes, broadcast
"""
import numpy as np

_CL68 = 0.6827
SOURCE_NAMES = {0: "none (0)", 1: "data", 2: "data, 0 triggers (0 ± UL)", 3: "GEV⊗Gauss fit"}


# =============================================================================
# 1. BUILD (run once)
# =============================================================================
def build_efficiency_table(fit_npz, out_npz, n_min_trig=1, pe_threshold=20.0):
    d = np.load(fit_npz, allow_pickle=True)

    inc_edges = np.asarray(d["inc_edges"], float)
    pc_edges  = np.asarray(d["pc_edges"],  float)
    if np.any(pc_edges <= 0):
        raise ValueError("pc_edges must be > 0 for log interpolation")

    n_events   = np.asarray(d["n_events"], float)
    n_total    = np.asarray(d["n_total"],  float)
    frac_valid = np.asarray(d["frac_valid"]).astype(bool)

    eff_d = np.asarray(d["efficiency_data"],     float)
    err_d = np.asarray(d["efficiency_data_err"], float)
    eff_g = np.asarray(d["efficiency_gev"],      float)
    err_g = np.asarray(d["efficiency_gev_err"],  float)
    ok_g  = np.asarray(d["success_gev"]).astype(bool)

    with np.errstate(invalid="ignore"):
        n_trig = eff_d * n_total

    has_ev    = frac_valid & (n_events > 0)
    has_trig  = has_ev & np.isfinite(eff_d) & (n_trig >= n_min_trig)
    zero_trig = has_ev & ~has_trig

    # one-sided 68% Clopper-Pearson upper limit for k = 0
    ul = np.where(has_ev, 1.0 - (1.0 - _CL68) ** (1.0 / np.maximum(n_events, 1.0)), 0.0)

    # ---------- data-only table ----------
    E_dat = np.zeros_like(eff_d)
    S_dat = np.zeros_like(eff_d)
    src_dat = np.zeros(eff_d.shape, np.int8)

    E_dat[has_trig] = eff_d[has_trig]
    S_dat[has_trig] = np.nan_to_num(err_d[has_trig], nan=0.0)
    src_dat[has_trig] = 1
    S_dat[zero_trig] = ul[zero_trig]
    src_dat[zero_trig] = 2

    # ---------- hybrid table ----------
    E_hyb, S_hyb, src_hyb = E_dat.copy(), S_dat.copy(), src_dat.copy()
    use_g = zero_trig & ok_g & np.isfinite(eff_g)
    g_err = np.where(np.isfinite(err_g), err_g, eff_g)   # unknown err -> 100 %
    E_hyb[use_g] = eff_g[use_g]
    S_hyb[use_g] = g_err[use_g]
    src_hyb[use_g] = 3

    np.clip(E_dat, 0, 1, out=E_dat)
    np.clip(E_hyb, 0, 1, out=E_hyb)

    np.savez(out_npz,
             inc_edges=inc_edges, pc_edges=pc_edges,
             eff_data=E_dat, err_data=S_dat, src_data=src_dat,
             eff_hybrid=E_hyb, err_hybrid=S_hyb, src_hybrid=src_hyb,
             n_events=n_events, n_min_trig=n_min_trig, pe_threshold=pe_threshold)

    print(f"saved {out_npz}  ({E_dat.shape[0]} inc x {E_dat.shape[1]} pc bins)")
    for k, name in SOURCE_NAMES.items():
        print(f"  hybrid  {name:28s}: {(src_hyb == k).sum()}")
    return out_npz


# =============================================================================
# 2. FAST LOOKUP
# =============================================================================
class TriggerEfficiency:
    def __init__(self, table_npz):
        t = np.load(table_npz)                       # no pickle -> fast & safe
        self.inc_edges = t["inc_edges"]
        self.pc_edges  = t["pc_edges"]
        self.lpc_edges = np.log10(self.pc_edges)
        self.inc_c = 0.5 * (self.inc_edges[:-1] + self.inc_edges[1:])
        self.lpc_c = 0.5 * (self.lpc_edges[:-1] + self.lpc_edges[1:])
        self.pe_threshold = float(t["pe_threshold"])
        self._tab = {
            "hybrid": (t["eff_hybrid"], t["err_hybrid"], t["src_hybrid"]),
            "data":   (t["eff_data"],   t["err_data"],   t["src_data"]),
        }

    @staticmethod
    def _bracket(c, v):
        n = c.size
        if n == 1:
            k = np.zeros(v.size, np.intp)
            return k, k, np.zeros(v.size)
        k = np.clip(np.searchsorted(c, v) - 1, 0, n - 2)
        t = np.clip((v - c[k]) / (c[k + 1] - c[k]), 0.0, 1.0)  # clamp = hold edge value
        return k, k + 1, t

    def __call__(self, inc, pc, mode="hybrid", method="linear", return_source=False):
        """
        inc : incidence angle [deg], any shape
        pc  : photon count, any shape (broadcast with inc)
        mode   : "hybrid" | "data"
        method : "linear" (bilinear in inc, log10 pc) | "nearest" (exact bin)
        returns eff, err  (+ source code array if return_source)
        """
        E, S, SRC = self._tab[mode]
        inc, pc = np.broadcast_arrays(np.asarray(inc, float), np.asarray(pc, float))
        shape = inc.shape
        x = inc.ravel()
        with np.errstate(divide="ignore", invalid="ignore"):
            y = np.log10(pc.ravel())

        ie, le = self.inc_edges, self.lpc_edges
        inside = (np.isfinite(x) & np.isfinite(y)
                  & (x >= ie[0]) & (x <= ie[-1]) & (y >= le[0]) & (y <= le[-1]))

        eff = np.zeros(x.size)
        err = np.zeros(x.size)
        src = np.zeros(x.size, np.int8)
        xi, yi = x[inside], y[inside]

        # nearest bin (always needed for source; is the answer for "nearest")
        ni = np.clip(np.searchsorted(ie, xi, "right") - 1, 0, E.shape[0] - 1)
        nj = np.clip(np.searchsorted(le, yi, "right") - 1, 0, E.shape[1] - 1)

        if method == "nearest":
            eff[inside] = E[ni, nj]
            err[inside] = S[ni, nj]
        elif method == "linear":
            i0, i1, t = self._bracket(self.inc_c, xi)
            j0, j1, u = self._bracket(self.lpc_c, yi)
            w00 = (1 - t) * (1 - u); w10 = t * (1 - u)
            w01 = (1 - t) * u;       w11 = t * u
            eff[inside] = (w00 * E[i0, j0] + w10 * E[i1, j0]
                           + w01 * E[i0, j1] + w11 * E[i1, j1])
            err[inside] = (w00 * S[i0, j0] + w10 * S[i1, j0]
                           + w01 * S[i0, j1] + w11 * S[i1, j1])
        else:
            raise ValueError("method must be 'linear' or 'nearest'")

        src[inside] = SRC[ni, nj]
        out = (eff.reshape(shape), err.reshape(shape))
        if return_source:
            out += (src.reshape(shape),)
        return out if shape else tuple(float(o) if o.dtype.kind == "f" else int(o)
                                       for o in out)