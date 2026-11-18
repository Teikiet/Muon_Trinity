// path_b.cpp -- independent trajectory cross-check of the geometry table:
// the traced geometric path vs an independent Bouguer integral of the same ray.
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
// independent geometric path along the ray: L = int dz / sqrt(1 - (p/(n(z) z))^2)
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

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("geometry table: valid=%d hNodes=%zu muNodes=%zu memory=%.1f MB\n", int(tab.valid()),
              tab.hNodes(), tab.muNodes(), tab.memoryMB());
  std::printf("independent trajectory cross-check (traced L vs Bouguer integral)\n");
  double worst = 0.0, worstRel = 0.0;
  char info[220] = {0};
  for (double h_e : {1000.0, 1500.0, 5000.0, 6000.0, 18000.0, 32000.0, 46000.0})
    for (double lam : {400.0, 690.0}) {
      double const f = Fdisp(lam);
      double const re = R + h_e;
      for (double xt : {0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999}) {
        double const p = xt * bCrit(ri, f);
        double const nE = nPhase(ri, re, f);
        double const sinT = p / (nE * re);
        if (sinT >= 1.0) continue;
        double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
        CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
        if (r.status != CurvedRayTable::Status::Hit) continue;
        double const Lref = lenRef(ri, p, f, h_e, 0.1);
        double const dL = std::abs(Lref - r.geometric_path_m);
        if (dL > worst) {
          worst = dL;
          worstRel = dL / Lref;
          std::snprintf(info, sizeof(info), "h=%.0f lam=%.0f xi=%.4f  L=%.3f m", h_e, lam, xt, Lref);
        }
      }
    }
  std::printf("  worst |L_trace - L_ref| = %.4e m (relative %.2e)  [%s]\n", worst, worstRel, info);
  std::printf("\n  breakdown at h=1500 m, lambda=690 nm:\n");
  for (double xt : {0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999}) {
    double const f = Fdisp(690.0), re = R + 1500.0;
    double const p = xt * bCrit(ri, f);
    double const sinT = p / (nPhase(ri, re, f) * re);
    if (sinT >= 1.0) continue;
    double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
    CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
    if (r.status != CurvedRayTable::Status::Hit) continue;
    double const Lref = lenRef(ri, p, f, 1500.0, 0.1);
    std::printf("    xi=%.4f  mu=%.6f  L_trace=%.6f  L_ref=%.6f  dL=%.3e m\n", xt,
                std::sqrt(std::max(0.0, 1.0 - xt * xt)), r.geometric_path_m, Lref,
                std::abs(Lref - r.geometric_path_m));
  }
  return 0;
}
