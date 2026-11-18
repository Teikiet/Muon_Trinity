// stage3Cp_b.cpp -- Stage 3, Step C': rebuild ONE column with a scaled substep count.
//
// Step C established (measurement, log_stage3C_b.txt):
//   * the stored NODE values are wrong by ~2.15 m  -> the cause is (b), not (b')
//   * the error accumulates up the column, fastest in the lowest few km
//   * the stored table is a STAIRCASE in mu at the nSubBase transitions, with a
//     factor 2.05 jump at mu = 1e-3 where nSubBase goes 8 -> 4
//
// This program reproduces buildColumn() locally (the production header is NOT
// touched; buildColumn is private and has no substep override) and re-runs the
// failing column at nSubBase x1 (self-check), x4 and x16.
//
// Predictions under test:
//   4  RK4 error ~ nSub^-4  =>  x4 should cut 2.15 m to ~8.4 mm, x16 to ~0.03 mm
//   3  the staircase / the mu = 1e-3 jump should disappear once both sides of the
//      threshold are built with the SAME substep count
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage3Cp_b.cpp -o stage3Cp_b
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

// ---- the h grid, exactly as CurvedRayTable builds it (lines 306-315) --------
static std::vector<double> makeZNodes(double hTop, double zBreak, double dhLo, double dhHi) {
  std::vector<double> z;
  z.push_back(0.0);
  for (double h = 0.0;;) {
    double const step = (h < zBreak) ? dhLo : dhHi;
    double nh = h + step;
    if (nh > hTop) nh = hTop;
    z.push_back(nh);
    if (nh >= hTop - 1e-9) break;
    h = nh;
  }
  return z;
}

// ---- the build's own substep rule (lines 589-601) --------------------------
static int nSubOf(RefractiveIndexTable const& ri, double mu, double f) {
  int nSubBase = 4;
  if (mu < 1.0e-3) {
    double const r0 = R;
    double const Nv = ri.N(r0);
    double const dN = ri.logSlope(ri.cellIndex(r0)) * Nv;
    double const a1 = 1.0 / r0 + dN / (1.0 + Nv);
    double const af = 1.0 / r0 + (f * dN) / (1.0 + f * Nv);
    double const stiff = (af > 0.0 && a1 > 0.0) ? std::max(1.0, a1 / af) : 1.0;
    double const k = std::min(256.0, std::ceil(1.0e-3 / mu) * stiff);
    nSubBase = static_cast<int>(4.0 * k);
  }
  return nSubBase;
}

// ---- local copy of buildColumn (CurvedRayTable.hpp:549-633) ----------------
// Carries only theta and L: dtheta/dh does not depend on phi or Sigma, so those
// states are irrelevant to the path length.  Everything else is verbatim.
static void buildColumnCopy(RefractiveIndexTable const& ri, double rObs,
                            std::vector<double> const& z, double mu, double f, int nSubBase,
                            int firstCellBoost, std::vector<double>& lOut) {
  std::size_t const n_h = z.size();
  lOut.assign(n_h, 0.0);
  double th = std::acos(std::min(1.0, std::max(0.0, mu)));
  double L = 0.0;
  auto rhs = [&](double hh, double t, double& dth, double& dl) {
    double const r = rObs + hh;
    double const Nz = ri.N(r);
    double const n = 1.0 + Nz * f;
    double const dN = ri.logSlope(ri.cellIndex(r)) * Nz;
    double const npn = (f * dN) / n;
    double const k = 1.0 / r + npn;
    double const c = std::cos(t), s = std::sin(t);
    double const cc = (c > 1e-12) ? c : 1e-12;
    dth = -s * k / cc;
    dl = 1.0 / cc;
  };
  for (std::size_t i = 0; i + 1 < n_h; ++i) {
    double const hA = z[i], hB = z[i + 1];
    int const nSub = (i == 0) ? firstCellBoost * nSubBase : nSubBase;
    double const dh = (hB - hA) / double(nSub);
    double hh = hA;
    for (int m = 0; m < nSub; ++m) {
      double k1t, k1l, k2t, k2l, k3t, k3l, k4t, k4l;
      rhs(hh, th, k1t, k1l);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k1t, k2t, k2l);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k2t, k3t, k3l);
      rhs(hh + dh, th + dh * k3t, k4t, k4l);
      th += dh / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
      L += dh / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
      hh += dh;
    }
    lOut[i + 1] = L;
  }
}

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 3, Step C': one-column rebuild with scaled substeps ===\n\n");
  std::printf("table: valid=%d hNodes=%zu muNodes=%zu memory=%.1f MB  muBreak=%.9f\n",
              int(tab.valid()), tab.hNodes(), tab.muNodes(), tab.memoryMB(), tab.muBreak());
  if (!tab.valid()) return 1;

  std::vector<double> const z = makeZNodes(H_TOP, 30000.0, 10.0, 100.0);
  std::printf("h grid: %zu nodes (expect %zu)\n", z.size(), tab.hNodes());

  double const hE = 46000.0, lam = 400.0, f = Fdisp(lam);
  double const xiRay = 0.9999;
  double const muRay = std::sqrt(std::max(0.0, 1.0 - xiRay * xiRay));
  std::size_t j = 0;
  double t = 0.0;
  tab.muCell(muRay, j, t);
  double const muJ = tab.muNode(j);
  std::size_t const iE = tab.hCell(hE);
  std::printf("ray mu=%.9f -> cell j=%zu (t=%.6f), column node mu=%.9e, hNode(hCell)=%.1f\n\n", muRay,
              j, t, muJ, tab.hNode(iE));

  // ---- reference at the node mu, at several heights ------------------------
  double const xiN = std::sqrt(std::max(0.0, 1.0 - muJ * muJ));
  double const pN = xiN * bCrit(ri, f);
  double const hs[5] = {1000.0, 5000.0, 20000.0, 46000.0, 60000.0};
  double dref[5];
  std::printf("reference DL_ref(h) at the column's node mu (dr=0.1 m):\n");
  for (int q = 0; q < 5; ++q) {
    dref[q] = lenRef(ri, pN, f, hs[q], 0.1) - tab.vacuumL(hs[q], muJ);
    std::printf("   h=%-8.0f DL_ref = %+.6f m\n", hs[q], dref[q]);
  }

  // ---- 1. SELF-CHECK: reproduce the stored column -------------------------
  std::printf("\n1. self-check: local copy with the build's own nSubBase (%d), boost x4\n",
              nSubOf(ri, muJ, f));
  std::vector<double> lCopy;
  buildColumnCopy(ri, R, z, muJ, f, nSubOf(ri, muJ, f), 4, lCopy);
  std::printf("   %-8s %-18s %-18s %s\n", "h[m]", "copy DL[m]", "storedDL[m]", "copy-stored");
  double maxSelf = 0.0;
  for (int q = 0; q < 5; ++q) {
    std::size_t const i = tab.hCell(hs[q]);
    double const dc = lCopy[i] - tab.vacuumL(tab.hNode(i), muJ);
    double const ds = tab.storedDL(i, j);
    maxSelf = std::max(maxSelf, std::fabs(dc - ds));
    std::printf("   %-8.0f %+-18.6f %+-18.6f %+.3e\n", tab.hNode(i), dc, ds, dc - ds);
  }
  std::printf("   -> max |copy - stored| = %.3e m  %s\n", maxSelf,
              (maxSelf < 0.01) ? "(COPY VERIFIED)" : "(COPY MISMATCH -- distrust the scaling below)");

  // ---- 2. SCALE TEST at h = 46000 -----------------------------------------
  std::printf("\n2. substep scaling at h = 46000 m (prediction 4)\n");
  std::printf("   %-12s %-8s %-18s %-18s %-12s %s\n", "scale", "nSubBase", "DL_copy[m]", "DL_ref[m]",
              "error[m]", "observed order");
  int const base = nSubOf(ri, muJ, f);
  double err1 = 0.0;
  int const scales[3] = {1, 4, 16};
  for (int si = 0; si < 3; ++si) {
    int const sc = scales[si];
    std::vector<double> lc;
    buildColumnCopy(ri, R, z, muJ, f, base * sc, 4, lc);
    std::size_t const i = tab.hCell(hE);
    double const dc = lc[i] - tab.vacuumL(tab.hNode(i), muJ);
    double const err = dc - dref[3];
    if (si == 0) err1 = err;
    double order = 0.0;
    if (si > 0 && err != 0.0 && err1 != 0.0) order = std::log(std::fabs(err1 / err)) / std::log(double(sc));
    std::printf("   x%-11d %-8d %+-18.6f %+-18.6f %+-12.6f %s\n", sc, base * sc, dc, dref[3], err,
                (si == 0) ? "-" : "");
    if (si > 0)
      std::printf("   %-12s %-8s %-18s %-18s %-12s p = %.2f\n", "", "", "", "", "", order);
  }

  // ---- 3. where the error lives, per scale --------------------------------
  std::printf("\n3. error vs height, per scale (error = DL_copy - DL_ref)\n");
  std::printf("   %-10s %-14s %-14s %-14s %-14s %s\n", "h[m]", "scale x1", "scale x4", "scale x16",
              "", "");
  std::vector<double> l1, l4, l16;
  buildColumnCopy(ri, R, z, muJ, f, base, 4, l1);
  buildColumnCopy(ri, R, z, muJ, f, base * 4, 4, l4);
  buildColumnCopy(ri, R, z, muJ, f, base * 16, 4, l16);
  for (int q = 0; q < 5; ++q) {
    std::size_t const i = tab.hCell(hs[q]);
    double const e1 = l1[i] - tab.vacuumL(tab.hNode(i), muJ) - dref[q];
    double const e4 = l4[i] - tab.vacuumL(tab.hNode(i), muJ) - dref[q];
    double const e16 = l16[i] - tab.vacuumL(tab.hNode(i), muJ) - dref[q];
    std::printf("   %-10.0f %+-14.6f %+-14.6f %+-14.6f\n", tab.hNode(i), e1, e4, e16);
  }

  // ---- 4. the mu = 1e-3 threshold, causally ---------------------------------
  std::printf("\n4. threshold test at h = 1000 m (prediction 3): does the jump survive\n");
  std::printf("   when BOTH sides are built with the same substep count?\n");
  std::printf("   %-16s %-8s %-18s %-18s %s\n", "muNode", "nSubBase", "DL_copy[m]", "DL_ref[m]",
              "error[m]");
  for (int q = 0; q < 2; ++q) {
    double const muT = (q == 0) ? 9.90e-4 : 1.01e-3;
    double const xiT = std::sqrt(std::max(0.0, 1.0 - muT * muT));
    double const pT = xiT * bCrit(ri, f);
    double const drefT = lenRef(ri, pT, f, 1000.0, 0.05) - tab.vacuumL(1000.0, muT);
    std::size_t const iT = tab.hCell(1000.0);
    std::vector<double> lc;
    // (a) the build's own rule
    buildColumnCopy(ri, R, z, muT, f, nSubOf(ri, muT, f), 4, lc);
    double const eA = lc[iT] - tab.vacuumL(1000.0, muT) - drefT;
    std::printf("   %-16.6e %-8d %+-18.6f %+-18.6f %+.6f   <- build's own rule\n", muT,
                nSubOf(ri, muT, f), lc[iT] - tab.vacuumL(1000.0, muT), drefT, eA);
    // (b) common nSubBase = 16 on both sides
    buildColumnCopy(ri, R, z, muT, f, 16, 4, lc);
    double const eB = lc[iT] - tab.vacuumL(1000.0, muT) - drefT;
    std::printf("   %-16.6e %-8d %+-18.6f %+-18.6f %+.6f   <- forced nSubBase=16\n", muT, 16,
                lc[iT] - tab.vacuumL(1000.0, muT), drefT, eB);
  }

  std::printf("\n=== end of Step C' ===\n");
  return 0;
}
