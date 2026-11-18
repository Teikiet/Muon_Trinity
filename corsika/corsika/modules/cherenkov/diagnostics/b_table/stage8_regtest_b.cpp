// stage8_regtest_b.cpp -- Stage 8, Task 1.3 companion: run the SUBSTANCE of the new
// Catch2 regression test standalone.
//
// The Catch2 suite itself cannot be built in this environment (tests/modules/
// testCherenkovCurvedRayTable.cpp includes Cherenkov.hpp, which needs spdlog, and no
// spdlog is present in the repo or the conan cache).  This harness exercises the same two
// assertions against the production table, which is what the regression is for:
//
//   (1) both bracketing NODE values on the Stage 3 ray, against the w-reference: <= 0.03 m
//   (2) the total arrival time on that ray, against the w-reference: <= 0.1 ns
//
// Diagnostics only.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage8_regtest_b.cpp -o stage8_regtest_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0, C_LIGHT = 2.99792458e8;

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
static void refW(RefractiveIndexTable const& ri, double p, double f, double hTop, double dw,
                 double& L, double& S) {
  L = 0.0;
  S = 0.0;
  double const wTop = std::sqrt(hTop);
  for (double w = 0.0; w < wTop - 1e-12;) {
    double const st = std::min(dw, wTop - w);
    double const wn = w + st;
    if (!(wn > w)) break; // stall guard
    double const wm = w + 0.5 * st;
    double const z = R + wm * wm;
    double const s = p / ((1.0 + ri.N(z) * f) * z);
    if (s >= 1.0) return;
    double const jac = 2.0 * wm;
    double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s));
    L += inv * jac * st;
    S += ri.N(z) * inv * jac * st;
    w = wn;
  }
}

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 8: standalone run of the new regression assertions ===\n\n");
  if (!tab.valid()) return 1;
  double const f = Fdisp(400.0);
  double const hE = 46000.0, xiRay = 0.9999;
  double const muRay = std::sqrt(std::max(0.0, 1.0 - xiRay * xiRay));
  std::size_t j = 0;
  double tq = 0.0;
  tab.muCell(muRay, j, tq);
  std::size_t const i = tab.hCell(hE);
  int nFail = 0;

  std::printf("(1) bracketing node values vs the w-reference (limit 0.03 m)\n");
  for (int which = 0; which < 2; ++which) {
    double const muN = tab.muNode(j + std::size_t(which));
    double const xiN = std::sqrt(std::max(0.0, 1.0 - muN * muN));
    double const pN = xiN * (1.0 + ri.N(R) * f) * R;
    double Lr = 0.0, Sr = 0.0;
    refW(ri, pN, f, hE, 1.0e-3, Lr, Sr);
    double const DLref = Lr - tab.vacuumL(hE, muN);
    double const err = std::abs(tab.storedDL(i, j + std::size_t(which)) - DLref);
    bool const ok = err <= 0.03;
    if (!ok) ++nFail;
    std::printf("    j=%zu mu=%.9e  storedDL=%.6f  DL_ref=%.6f  |err|=%.6f m  %s\n",
                j + std::size_t(which), muN, tab.storedDL(i, j + std::size_t(which)), DLref, err,
                ok ? "PASS" : "FAIL");
  }

  std::printf("\n(2) total arrival time on the Stage 3 ray (limit 0.1 ns)\n");
  double const re = R + hE;
  double const p = xiRay * (1.0 + ri.N(R) * f) * R;
  double const sinT = p / ((1.0 + ri.N(re) * f) * re);
  double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
  CurvedRayTable::Result const res = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
  if (res.status != CurvedRayTable::Status::Hit) {
    std::printf("    FAIL: status is not Hit\n");
    ++nFail;
  } else {
    double Lr = 0.0, Sr = 0.0;
    refW(ri, p, f, hE, 1.0e-3, Lr, Sr);
    double const tRef = (Lr + f * Sr) / C_LIGHT;
    double const dt = std::abs(res.time_s - tRef) * 1e9;
    bool const ok = dt <= 0.1;
    if (!ok) ++nFail;
    std::printf("    t_prod = %.12e s   t_ref = %.12e s   |dt| = %.6f ns  %s\n", res.time_s, tRef,
                dt, ok ? "PASS" : "FAIL");
  }
  std::printf("\n(3) mu = 1.0 (radial emission) must be a Hit with zero displacement\n");
  {
    double const hs2[6] = {100.0, 1000.0, 5000.0, 20000.0, 46000.0, 50000.0};
    for (int q = 0; q < 6; ++q) {
      double const re = R + hs2[q];
      CurvedRayTable::Result const rr = tab.trace(0.0, 0.0, re, 0.0, 0.0, -1.0, 1.0, 1.0);
      // d = r_obs*Phi is not exposed in Result, but Phi = 0 puts the impact point on the
      // emission radial, so the transverse coordinates must vanish.  This is well
      // conditioned, unlike R*acos(z/R), which is singular at z/R = 1.
      double const dtr = std::sqrt(rr.x * rr.x + rr.y * rr.y);
      bool const ok = rr.hit && rr.status == CurvedRayTable::Status::Hit &&
                      dtr <= 1e-9 && std::abs(rr.dx) <= 1e-12 &&
                      std::abs(rr.dy) <= 1e-12 && rr.dz <= -1.0 + 1e-12;
      if (!ok) ++nFail;
      std::printf("    h=%-7.0f status=%d hit=%d  transverse=%.3e m  dir=(%.3e,%.3e,%.3f)  %s\n",
                  hs2[q], int(rr.status), int(rr.hit), dtr, rr.dx, rr.dy, rr.dz, ok ? "PASS" : "FAIL");
    }
  }

  std::printf("\n%s (%d failing assertion(s))\n", nFail == 0 ? "ALL ASSERTIONS PASS" : "FAILURES",
              nFail);
  std::printf("\n=== end of standalone regression run ===\n");
  return nFail == 0 ? 0 : 1;
}
