// map_b.cpp -- coordinate definition, mapping test against the REAL geometry table,
// and validation against an independent reference for the ACTUAL configured ray
// (dispersion included). Diagnostic only; no production code is touched.
//
// Coordinate (unambiguous):
//     p       = n_phase(r_e, lambda) * r_e * sin(theta_e)     [Bouguer invariant]
//     p_crit  = n_phase(r_obs, lambda) * r_obs
//     xi      = p / p_crit          <- THE TABLE AXIS (impact-parameter ratio)
//     mu_obs  = |cos theta_obs| = sqrt(1 - xi^2)   (direction cosine; NOT the axis)
// xi and mu_obs are different quantities and are never conflated below.
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0;   // earth radius [m]
static constexpr double H_TOP = 60000.0; // table ceiling above ground [m]

// ---- the optical configuration (what "active immutable configuration" means) ----
// Refractivity shape N(h) = (n-1) at f = 1. Two terms: free troposphere + boundary layer.
static double Nshape(double h) {
  return 2.20e-4 * std::exp(-h / 7500.0) + 0.55e-4 * std::exp(-h / 1600.0);
}
// Phase dispersion F(lambda), normalised to 1 at 400 nm (representative Cauchy form).
static double Fdisp(double lam_nm) {
  double const l = lam_nm / 1000.0;
  return (1.0 + 0.0032 / (l * l)) / (1.0 + 0.0032 / 0.16);
}
// The RI table the geometry route uses; production interpolates ln N on a uniform grid.
static RefractiveIndexTable makeRI() {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= H_TOP + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(std::log(Nshape(h)));
  }
  return RefractiveIndexTable(alt, lnN, R);
}
// Phase index actually used by the ray: n = 1 + N(z) * f.  Production's f convention
// f = (n-1)/N(r_emit) gives exactly f = F(lambda) for this separable model, so passing
// fPhase = Fdisp(lambda) into trace() reproduces the production call site exactly.
static double nPhase(RefractiveIndexTable const& ri, double zAbs, double f) {
  return 1.0 + ri.N(zAbs) * f;
}
static double bCrit(RefractiveIndexTable const& ri, double f) { return nPhase(ri, R, f) * R; }

// ---- extinction profile, reproduced exactly as the production callback does -------
static AtmosphericAbsorption g_ab;
static double alphaOf(double lam, double h) {
  double const dh = 50.0;
  auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h - dh * 0.5), hi = std::min(hr.second, h + dh * 0.5);
  double const a = hi - lo;
  if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam, hi) - g_ab.getOpticalDepth(lam, lo)) / a);
}

// ---- INDEPENDENT reference along the actual ray ----------------------------------
// tau = int alpha(z) / sqrt(1 - (p/(n_phase(z,lam) z))^2) dz, p fixed by the actual
// photon, n_phase carrying the actual dispersion. No table is consulted.
static double tauRef(RefractiveIndexTable const& ri, double p, double f, double lam,
                     double hTop, double dr) {
  double acc = 0.0, z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zm = z + 0.5 * st;
    double const s = p / (nPhase(ri, zm, f) * zm);
    if (s >= 1.0) return acc;
    acc += alphaOf(lam, zm - R) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
    z += st;
  }
  return acc;
}
// Same integrator, integrand 1 (geometric path length) -- cross-checks the traced ray.
static double lenRef(RefractiveIndexTable const& ri, double p, double f, double hTop, double dr) {
  double acc = 0.0, z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zm = z + 0.5 * st;
    double const s = p / (nPhase(ri, zm, f) * zm);
    if (s >= 1.0) return acc;
    acc += 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
    z += st;
  }
  return acc;
}

// ---- candidate A: cumulative optical depth on (h, xi, lambda) --------------------
struct TableA {
  std::vector<double> h, xi, lam;
  int nH = 0, nX = 0, nL = 0;
  std::vector<float> a;
  bool configured = false;                 // true: built with the configured dispersion
  double fOf(double lamNm) const { return configured ? Fdisp(lamNm) : 1.0; }
  double bytes() const { return double(a.size()) * 4.0; }

  double lookup(double hq, double xq, double lq) const {
    auto cell = [](std::vector<double> const& v, double q, int& i) {
      int const n = int(v.size());
      i = 0;
      while (i + 2 < n && v[i + 1] < q) ++i;
      double const d = v[i + 1] - v[i];
      double const t = (d > 0.0) ? (q - v[i]) / d : 0.0;
      return std::min(1.0, std::max(0.0, t));   // clamp, never extrapolate
    };
    int i0, x0, l0;
    double const th = cell(h, hq, i0);
    double const tx = cell(xi, xq, x0);
    double const tl = cell(lam, lq, l0);
    double acc = 0.0;
    for (int p = 0; p < 2; ++p)
      for (int q = 0; q < 2; ++q)
        for (int c = 0; c < 2; ++c) {
          double const w = (p ? th : 1.0 - th) * (q ? tx : 1.0 - tx) * (c ? tl : 1.0 - tl);
          if (w == 0.0) continue;
          acc += w * double(a[std::size_t(i0 + p) * nX * nL + (x0 + q) * nL + (l0 + c)]);
        }
    return acc;
  }
};

static std::vector<double> makeAxesXi() {
  std::vector<double> v;
  for (int i = 0; i < 48; ++i) v.push_back(2.0e-4 * std::pow(0.5 / 2.0e-4, double(i) / 47.0));
  for (int i = 1; i <= 208; ++i) {
    double const u = std::sqrt(0.5) * double(i) / 208.0;
    v.push_back(1.0 - u * u);
  }
  std::sort(v.begin(), v.end());
  return v;
}

static TableA buildTable(RefractiveIndexTable const& ri, bool configured, double* secs) {
  TableA t;
  t.configured = configured;
  t.h = {0.0, 25.0, 50.0, 100.0, 250.0, 500.0};
  for (int i = 1; i <= 50; ++i) t.h.push_back(1000.0 * i);
  t.xi = makeAxesXi();
  for (double w = 300.0; w <= 900.0 + 1e-9; w += 5.0) t.lam.push_back(w);
  t.nH = int(t.h.size());
  t.nX = int(t.xi.size());
  t.nL = int(t.lam.size());
  t.a.assign(std::size_t(t.nH) * t.nX * t.nL, 0.0f);
  auto const t0 = std::chrono::steady_clock::now();
  for (int ix = 0; ix < t.nX; ++ix) {
    for (int il = 0; il < t.nL; ++il) {
      double const lam = t.lam[il];
      double const f = t.fOf(lam);
      double const b = t.xi[ix] * bCrit(ri, f);   // xi normalised to THIS f's b_crit
      double acc = 0.0, z = R;
      int ih = 0;
      while (ih < t.nH) {
        double const zTop = R + t.h[ih];
        while (z < zTop - 1e-12) {
          double const st = std::min(5.0, zTop - z);
          double const zm = z + 0.5 * st;
          double const s = b / (nPhase(ri, zm, f) * zm);
          if (s >= 1.0) break;
          acc += alphaOf(lam, zm - R) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
          z += st;
        }
        t.a[std::size_t(ih) * t.nX * t.nL + ix * t.nL + il] = float(acc);
        ++ih;
      }
    }
  }
  *secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return t;
}

int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  RefractiveIndexTable ri = makeRI();
  std::printf("=== coordinate + mapping + ACTUAL-RAY validation (dispersion included) ===\n");
  std::printf("configuration: R=%.1f m  H_TOP=%.0f m  RI grid [%.1f, %.1f] km, %zu nodes\n", R,
              H_TOP, ri.zMin() / 1000.0, ri.zMax() / 1000.0, ri.size());
  std::printf("  N(0)=%.6e  N(10km)=%.6e\n", Nshape(0.0), Nshape(10000.0));
  {
    bool mono = true;
    double prev = -1.0;
    for (double h = 0.0; h <= H_TOP; h += 25.0) {
      double const nr = nPhase(ri, R + h, 1.0) * (R + h);
      if (nr < prev) mono = false;
      prev = nr;
    }
    std::printf("  n*r strictly increasing on [0,%.0f m]: %s -> min(n*r) at the detector\n", H_TOP,
                mono ? "YES" : "NO");
  }
  std::printf("  dispersion F(lambda): 300nm %.6f  400nm %.6f  500nm %.6f  700nm %.6f\n",
              Fdisp(300.0), Fdisp(400.0), Fdisp(500.0), Fdisp(700.0));
  std::printf("  b_crit(400nm)=%.3f m   b_crit(700nm)=%.3f m   spread %.3f m\n",
              bCrit(ri, Fdisp(400.0)), bCrit(ri, Fdisp(700.0)),
              bCrit(ri, Fdisp(400.0)) - bCrit(ri, Fdisp(700.0)));

  double gsec = 0.0;
  CurvedRayTable tab;
  {
    auto const t0 = std::chrono::steady_clock::now();
    tab = CurvedRayTable(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
    gsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  std::printf("\ngeometry table (REAL production class): valid=%d hNodes=%zu muNodes=%zu "
              "memory=%.1f MB build=%.1f s\n",
              int(tab.valid()), tab.hNodes(), tab.muNodes(), tab.memoryMB(), gsec);
  std::printf("  muMin=%.2e (more grazing than this is REJECTED, never clamped)\n", tab.muMin());

  double s1 = 0.0, s2 = 0.0;
  TableA A1 = buildTable(ri, false, &s1);
  std::printf("\nA1 (built at f=1)           %d x %d x %d = %.2f MB  build %.1f s\n", A1.nH, A1.nX,
              A1.nL, A1.bytes() / 1048576.0, s1);
  TableA A2 = buildTable(ri, true, &s2);
  std::printf("A2 (configured dispersion)  %d x %d x %d = %.2f MB  build %.1f s\n", A2.nH, A2.nX,
              A2.nL, A2.bytes() / 1048576.0, s2);

  // ---- the mapping test and the actual-ray validation -----------------------------
  double worst1 = 0.0, worst2 = 0.0, maxDxi = 0.0, maxDL = 0.0, maxFsec = 0.0;
  char i1[200] = {0}, i2[200] = {0}, imap[200] = {0};
  int n = 0, nrej = 0;
  for (double h_e : {1500.0, 6000.0, 18000.0, 32000.0, 46000.0}) {
    for (double lam : {310.0, 400.0, 500.0, 620.0, 690.0, 720.0, 850.0}) {
      double const f = Fdisp(lam);
      double const bc = bCrit(ri, f);
      double const re = R + h_e;
      for (double xt : {0.15, 0.35, 0.6, 0.8, 0.92, 0.97, 0.995, 0.9995}) {
        double const p = xt * bc;                       // the ACTUAL photon's invariant
        double const nE = nPhase(ri, re, f);
        double const sinT = p / (nE * re);
        if (sinT >= 1.0) { ++nrej; continue; }
        double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
        CurvedRayTable::Result const r =
            tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
        if (r.status != CurvedRayTable::Status::Hit) { ++nrej; continue; }
        // (1) mapping: the geometry route's own mu must reproduce xi = p/p_crit
        double const muTrace = -std::cos(r.theta_obs);
        double const xiTrace = std::sqrt(std::max(0.0, 1.0 - muTrace * muTrace));
        double const dXi = std::abs(xiTrace - xt);
        // (2) same ray: the traced geometric path must match the independent integral
        double const Lref = lenRef(ri, p, f, h_e, 0.5);
        double const dL = std::abs(Lref - r.geometric_path_m);
        // (3) optical depth of the ACTUAL ray vs both tables
        double const tRef = tauRef(ri, p, f, lam, h_e, 0.5);
        double const tRefF1 = tauRef(ri, p, 1.0, lam, h_e, 0.5);
        double const tA1 = A1.lookup(h_e, xt, lam);
        double const tA2 = A2.lookup(h_e, xt, lam);
        double const dT1 = std::abs(std::exp(-tA1) - std::exp(-tRef));
        double const dT2 = std::abs(std::exp(-tA2) - std::exp(-tRef));
        double const dF = std::abs(std::exp(-tRef) - std::exp(-tRefF1));
        ++n;
        if (dXi > maxDxi) { maxDxi = dXi; std::snprintf(imap, sizeof(imap),
            "h=%.0f lam=%.0f xi=%.4f  xi_trace=%.6f  dxi=%.2e  dL=%.2e m", h_e, lam, xt, xiTrace, dXi, dL); }
        if (dL > maxDL) maxDL = dL;
        if (dF > maxFsec) maxFsec = dF;
        if (dT1 > worst1) { worst1 = dT1; std::snprintf(i1, sizeof(i1),
            "h=%.0f lam=%.0f xi=%.4f tau=%.4f", h_e, lam, xt, tRef); }
        if (dT2 > worst2) { worst2 = dT2; std::snprintf(i2, sizeof(i2),
            "h=%.0f lam=%.0f xi=%.4f tau=%.4f", h_e, lam, xt, tRef); }
      }
    }
  }
  std::printf("\nrays traced: %d accepted, %d inadmissible\n", n, nrej);
  std::printf("MAPPING TEST (geometry route vs independent xi = p/p_crit)\n");
  std::printf("  worst |xi_trace - xi| = %.3e   [%s]\n", maxDxi, imap);
  std::printf("  worst |L_trace  - L_ref| = %.3e m  (same physical ray)\n", maxDL);
  std::printf("ACTUAL-RAY VALIDATION (reference = true ray, dispersion included)\n");
  std::printf("  A1 f=1 table       worst |dT| = %.4e  [%s]\n", worst1, i1);
  std::printf("  A2 configured      worst |dT| = %.4e  [%s]\n", worst2, i2);
  std::printf("  (dispersion-only sensitivity of the reference itself: %.4e)\n", maxFsec);
  std::printf("  gate: A1 %s   A2 %s\n", (worst1 <= 1e-3) ? "PASS" : "FAIL",
              (worst2 <= 1e-3) ? "PASS" : "FAIL");
  {
    std::ifstream st("/proc/self/status");
    std::string l;
    while (std::getline(st, l))
      if (l.rfind("VmHWM", 0) == 0 || l.rfind("VmRSS", 0) == 0) std::printf("  %s\n", l.c_str());
  }
  return 0;
}
