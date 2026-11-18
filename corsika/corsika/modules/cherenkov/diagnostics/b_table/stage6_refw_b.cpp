// stage6_refw_b.cpp -- Stage 6, Task 2: a reference integrator in w = sqrt(z - R).
//
// Stage 5's reference is a midpoint march in z, and it is reference-limited below
// mu = 1e-3 (at mu = 2e-4 it moves 15.4 ns between dr = 0.1 and 0.02).  This program
// integrates in w instead, where the near-ground resolution is naturally fine, and
// covers BOTH L and Sigma:
//
//   z = R + w^2,  dz = 2 w dw
//   L     = int 2 w dw / sqrt(1 - s^2)
//   Sigma = int 2 w N(z) dw / sqrt(1 - s^2),   s = p / (n(z) z),  n = 1 + N(z) f
//
// Validations required by the brief:
//   (a) agrees with the Stage 5 z-march reference to 1 cm on 5 rays with mu >= 1e-2
//   (b) converges on its own (dw halved twice) at mu = 2e-4, 5e-4, 1e-3 and at
//       h = 59900 / 60000
//   (c) DL increases monotonically with h at 59900 -> 60000
//
// Diagnostics only; the stall guard is carried as required.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage6_refw_b.cpp -o stage6_refw_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0;

static double Nshape(double h) {
  return 2.20e-4 * std::exp(-h / 7500.0) + 0.55e-4 * std::exp(-h / 1600.0);
}
static double Fdisp(double lam_nm) {
  double const l = lam_nm / 1000.0;
  return (1.0 + 0.0032 / (l * l)) / (1.0 + 0.0032 / 0.16);
}
static RefractiveIndexTable makeRI() {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= H_TOP + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(std::log(Nshape(h)));
  }
  return RefractiveIndexTable(alt, lnN, R);
}
static double nPhase(RefractiveIndexTable const& ri, double zAbs, double f) {
  return 1.0 + ri.N(zAbs) * f;
}
static double bCrit(RefractiveIndexTable const& ri, double f) { return nPhase(ri, R, f) * R; }

// ---- the z-march reference used in Stages 3 and 5 (for validation (a)) -------
static double refZ(RefractiveIndexTable const& ri, double p, double f, double hTop, double dr) {
  double acc = 0.0, z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zn = z + st;
    if (!(zn > z)) break; // stall guard
    double const zm = z + 0.5 * st;
    double const s = p / (nPhase(ri, zm, f) * zm);
    if (s >= 1.0) return acc;
    acc += 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
    z = zn;
  }
  return acc;
}

// ---- the w-substitution reference, both L and Sigma -------------------------
static void refW(RefractiveIndexTable const& ri, double p, double f, double hTop, double dw,
                 double& L, double& S, bool* turned) {
  L = 0.0;
  S = 0.0;
  if (turned) *turned = false;
  double const wTop = std::sqrt(hTop);
  for (double w = 0.0; w < wTop - 1e-12;) {
    double const st = std::min(dw, wTop - w);
    double const wn = w + st;
    if (!(wn > w)) break; // stall guard
    double const wm = w + 0.5 * st;
    double const z = R + wm * wm;
    double const s = p / (nPhase(ri, z, f) * z);
    if (s >= 1.0) {
      if (turned) *turned = true;
      return;
    }
    double const jac = 2.0 * wm; // dz/dw at the midpoint
    double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s));
    L += inv * jac * st;
    S += ri.N(z) * inv * jac * st;
    w = wn;
  }
}

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 6, Task 2: w = sqrt(z-R) reference for L and Sigma ===\n\n");
  if (!tab.valid()) return 1;
  double const f = Fdisp(400.0);

  // ---- (a) agreement with the z-march reference, 5 rays with mu >= 1e-2 -----
  std::printf("(a) agreement with the Stage 5 z-march reference (dr=0.01) on 5 rays, mu >= 1e-2:\n");
  std::printf("    %-9s %-8s %-16s %-16s %s\n", "h[m]", "mu", "L_w[m]", "L_z[m]", "dL[mm]");
  double const ha[5] = {1000.0, 5000.0, 20000.0, 46000.0, 59900.0};
  double const mua[5] = {0.1, 0.1, 0.5, 0.1, 0.3};
  double worstA = 0.0;
  for (int q = 0; q < 5; ++q) {
    double const xi = std::sqrt(std::max(0.0, 1.0 - mua[q] * mua[q]));
    double const p = xi * bCrit(ri, f);
    double Lw = 0.0, Sw = 0.0;
    bool t1 = false;
    refW(ri, p, f, ha[q], 1.0e-3, Lw, Sw, &t1);
    double const Lz = refZ(ri, p, f, ha[q], 0.01);
    std::printf("    %-9.0f %-8.2f %-16.6f %-16.6f %.3f\n", ha[q], mua[q], Lw, Lz, (Lw - Lz) * 1e3);
    worstA = std::max(worstA, std::fabs(Lw - Lz) * 1e3);
  }
  std::printf("    worst |L_w - L_z| = %.3f mm  -> %s (limit 10 mm)\n\n", worstA,
              worstA < 10.0 ? "PASS" : "FAIL");

  // ---- (b) self-convergence at the hard rays -------------------------------
  std::printf("(b) self-convergence of the w-reference (dw halved twice):\n");
  std::printf("    %-9s %-10s %-16s %-16s %-16s %s\n", "h[m]", "mu", "L(dw=4e-3)", "L(dw=2e-3)",
              "L(dw=1e-3)", "drift[mm]");
  double const hb[5] = {46000.0, 46000.0, 46000.0, 59900.0, 60000.0};
  double const mub[5] = {2.0e-4, 5.0e-4, 1.0e-3, 1.412e-2, 1.412e-2};
  for (int q = 0; q < 5; ++q) {
    double const xi = std::sqrt(std::max(0.0, 1.0 - mub[q] * mub[q]));
    double const p = xi * bCrit(ri, f);
    double L1 = 0.0, S1 = 0.0, L2 = 0.0, S2 = 0.0, L3 = 0.0, S3 = 0.0;
    bool t1 = false, t2 = false, t3 = false;
    refW(ri, p, f, hb[q], 4.0e-3, L1, S1, &t1);
    refW(ri, p, f, hb[q], 2.0e-3, L2, S2, &t2);
    refW(ri, p, f, hb[q], 1.0e-3, L3, S3, &t3);
    std::printf("    %-9.0f %-10.2e %-16.6f %-16.6f %-16.6f %.3f%s\n", hb[q], mub[q], L1, L2, L3,
                (L3 - L1) * 1e3, t3 ? "  (turned)" : "");
  }

  // ---- (c) monotonicity of DL across 59900 -> 60000 ------------------------
  std::printf("\n(c) DL(h) must increase monotonically at the top of the column (mu=1.412e-2):\n");
  double const mu = 1.412e-2;
  double const xi = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  double const p = xi * bCrit(ri, f);
  double prev = -1.0;
  bool mono = true;
  double const hc[6] = {59600.0, 59700.0, 59800.0, 59900.0, 59950.0, 60000.0};
  for (int q = 0; q < 6; ++q) {
    double Lw = 0.0, Sw = 0.0;
    bool tt = false;
    refW(ri, p, f, hc[q], 1.0e-3, Lw, Sw, &tt);
    double const DL = Lw - tab.vacuumL(hc[q], mu);
    bool const bad = (prev > 0.0 && DL <= prev);
    std::printf("    h=%-8.0f DL_ref = %.6f m  %s\n", hc[q], DL, bad ? "NON-MONOTONIC" : "ok");
    if (bad) mono = false;
    prev = DL;
  }
  std::printf("    -> %s\n", mono ? "MONOTONIC PASS" : "NON-MONOTONIC FAIL");
  std::printf("\n=== end of Task 2 ===\n");
  return 0;
}
