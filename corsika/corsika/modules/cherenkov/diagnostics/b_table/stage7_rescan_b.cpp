// stage7_rescan_b.cpp -- Stage 7, Task C: re-scan the Stage 5 grid with the V4 rule.
//
// For every (mu, lambda, h) of the Stage 5 grid (plus the mu in [2e-4, 1e-3) bands):
//   - the column is integrated with the V4 rule (the refractive-index slope frozen to
//     the substep's start cell), nSubBase and the first-cell boost unchanged;
//   - the ground truth is the w = sqrt(z-R) reference from Stage 6 Task 2 (stall guard);
//   - mu = 2e-4 is treated as ELIGIBLE via an explicit tolerance, because the table's
//     own trace() rejects it through round-off in reconstructing mu from p/pCrit.
//
// The column is evaluated directly at the queried (mu, h) rather than through the
// table's (h, mu) Hermite interpolation.  That interpolation was measured in Stage 3 to
// contribute 2.3 mm = 0.008 ns on the failing ray, so bypassing it cannot hide a
// 0.1 ns requirement.
//
// Diagnostics only; no production source touched.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage7_rescan_b.cpp -o stage7_rescan_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>

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
static double nPhase(RefractiveIndexTable const& ri, double zAbs, double f) {
  return 1.0 + ri.N(zAbs) * f;
}
static double bCrit(RefractiveIndexTable const& ri, double f) { return nPhase(ri, R, f) * R; }

// w-substitution reference, both L and Sigma (Stage 6 Task 2), stall guard carried
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
    double const s = p / (nPhase(ri, z, f) * z);
    if (s >= 1.0) return;
    double const jac = 2.0 * wm;
    double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s));
    L += inv * jac * st;
    S += ri.N(z) * inv * jac * st;
    w = wn;
  }
}

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

// V4 column: slope frozen to the substep's start cell, for the L and Sigma integrands
static void buildColumnV4(RefractiveIndexTable const& ri, double rObs, std::vector<double> const& z,
                          double mu, double f, int nSubBase, int firstCellBoost, bool v4,
                          std::vector<double>& lOut, std::vector<double>& sOut) {
  std::size_t const n_h = z.size();
  lOut.assign(n_h, 0.0);
  sOut.assign(n_h, 0.0);
  double th = std::acos(std::min(1.0, std::max(0.0, mu)));
  double L = 0.0, Sig = 0.0, stepSlope = 0.0;
  auto rhs = [&](double hh, double t, double& dth, double& dl, double& ds) {
    double const r = rObs + hh;
    double const Nz = ri.N(r);
    double const dN = v4 ? (stepSlope * Nz) : (ri.logSlope(ri.cellIndex(r)) * Nz);
    double const k = 1.0 / r + (f * dN) / (1.0 + Nz * f);
    double const c = std::cos(t), s = std::sin(t);
    double const cc = (c > 1e-12) ? c : 1e-12;
    dth = -s * k / cc;
    dl = 1.0 / cc;
    ds = Nz / cc;
  };
  for (std::size_t i = 0; i + 1 < n_h; ++i) {
    double const hA = z[i], hB = z[i + 1];
    int const nSub = (i == 0) ? firstCellBoost * nSubBase : nSubBase;
    double const dh = (hB - hA) / double(nSub);
    double hh = hA;
    for (int m = 0; m < nSub; ++m) {
      stepSlope = ri.logSlope(ri.cellIndex(R + hh));
      double k1t, k1l, k1s, k2t, k2l, k2s, k3t, k3l, k3s, k4t, k4l, k4s;
      rhs(hh, th, k1t, k1l, k1s);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k1t, k2t, k2l, k2s);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k2t, k3t, k3l, k3s);
      rhs(hh + dh, th + dh * k3t, k4t, k4l, k4s);
      th += dh / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
      L += dh / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
      Sig += dh / 6.0 * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
      hh += dh;
    }
    lOut[i + 1] = L;
    sOut[i + 1] = Sig;
  }
}

static int bandOf(double mu) {
  if (mu < 5.0e-4) return 0;
  if (mu < 1.0e-3) return 1;
  if (mu < 2.0e-3) return 2;
  if (mu < 5.0e-3) return 3;
  if (mu < 1.0e-2) return 4;
  if (mu < 2.0e-2) return 5;
  if (mu < 5.0e-2) return 6;
  if (mu < 1.0e-1) return 7;
  return 8;
}

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 7, Task C: re-scan the Stage 5 grid with the V4 rule ===\n\n");
  if (!tab.valid()) return 1;
  std::vector<double> const z = makeZNodes(H_TOP, 30000.0, 10.0, 100.0);
  double const hs[10] = {100.0, 300.0, 1000.0, 3000.0, 5000.0,
                         10000.0, 20000.0, 32000.0, 46000.0, 59900.0};
  double const mus[13] = {2.0e-4, 3.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2,
                          2.0e-2, 5.0e-2, 1.0e-1, 0.3, 0.6, 1.0};
  double const lams[5] = {350.0, 400.0, 500.0, 690.0, 700.0};

  // ---- self-check: V4 OFF must reproduce the stored column -----------------
  double const f0 = Fdisp(400.0);
  double const muRay = std::sqrt(std::max(0.0, 1.0 - 0.9999 * 0.9999));
  std::size_t j = 0;
  double tq = 0.0;
  tab.muCell(muRay, j, tq);
  double const muJ = tab.muNode(j);
  std::size_t const iE = tab.hCell(46000.0);
  std::vector<double> lo, so;
  buildColumnV4(ri, R, z, muJ, f0, nSubOf(ri, muJ, f0), 4, false, lo, so);
  double const dOff = lo[iE] - tab.vacuumL(46000.0, muJ);
  std::printf("self-check V4 OFF: DL_copy = %.6f, storedDL = %.6f, diff = %.3e m  -> %s\n\n", dOff,
              tab.storedDL(iE, j), dOff - tab.storedDL(iE, j),
              std::fabs(dOff - tab.storedDL(iE, j)) < 2e-3 ? "MATCHES" : "MISMATCH");

  // ---- the scan ------------------------------------------------------------
  std::printf("mu = 2.0e-4 is treated as ELIGIBLE via an explicit tolerance: the table's own\n");
  std::printf("trace() rejects it through round-off when reconstructing mu from p/pCrit.\n\n");
  int tot[9][10] = {{0}}, fail[9][10] = {{0}};
  double worst[9][10];
  double cellUnc[9][10];
  double unc[9] = {0.0};
  for (int b = 0; b < 9; ++b)
    for (int q = 0; q < 10; ++q) {
      worst[b][q] = 0.0;
      cellUnc[b][q] = 0.0;
    }
  double topDt[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double topMu[5] = {0.0}, topH[5] = {0.0}, topL[5] = {0.0};
  int bandTot[9] = {0}, bandFail[9] = {0};
  double bandWorst[9] = {0.0};
  for (int il = 0; il < 5; ++il) {
    double const f = Fdisp(lams[il]);
    for (int im = 0; im < 13; ++im) {
      double const mu = mus[im];
      std::vector<double> lv, sv;
      buildColumnV4(ri, R, z, mu, f, nSubOf(ri, mu, f), 4, true, lv, sv);
      double const xi = std::sqrt(std::max(0.0, 1.0 - mu * mu));
      double const p = xi * bCrit(ri, f);
      int const b = bandOf(mu);
      for (int ih = 0; ih < 10; ++ih) {
        double const h = hs[ih];
        std::size_t const i = tab.hCell(h);
        double const vac = tab.vacuumL(h, mu);
        double Lr = 0.0, Sr = 0.0, Lr2 = 0.0, Sr2 = 0.0;
        refW(ri, p, f, h, 1.0e-3, Lr, Sr);
        refW(ri, p, f, h, 2.0e-3, Lr2, Sr2);
        // NOTE: buildColumnV4 returns the TOTAL geometric path, not a residual, so the
        // vacuum chord must NOT be added again here (an earlier revision did, which is
        // why a first run reported ~1e6 ns -- the size of vac/c itself).
        double const tV4 = (lv[i] + f * sv[i]) / C_LIGHT;
        double const tR = (Lr + f * Sr) / C_LIGHT;
        double const tR2 = (Lr2 + f * Sr2) / C_LIGHT;
        double const dt = std::fabs(tV4 - tR) * 1e9;
        double const u = std::fabs(tR2 - tR) * 1e9;
        ++tot[b][ih];
        if (dt > 0.1) ++fail[b][ih];
        if (dt > worst[b][ih]) worst[b][ih] = dt;
        if (u > unc[b]) unc[b] = u;
        ++bandTot[b];
        if (dt > 0.1) ++bandFail[b];
        if (dt > bandWorst[b]) bandWorst[b] = dt;
        if (u > cellUnc[b][ih]) cellUnc[b][ih] = u;
        if (dt > topDt[4]) {
          int q = 4;
          while (q > 0 && dt > topDt[q - 1]) {
            topDt[q] = topDt[q - 1];
            topMu[q] = topMu[q - 1];
            topH[q] = topH[q - 1];
            topL[q] = topL[q - 1];
            --q;
          }
          topDt[q] = dt;
          topMu[q] = mu;
          topH[q] = h;
          topL[q] = lams[il];
        }
      }
    }
  }

  // ---- report --------------------------------------------------------------
  char const* bl[9] = {"2e-4-5e-4", "5e-4-1e-3", "1e-3-2e-3", "2e-3-5e-3", "5e-3-1e-2",
                       "1e-2-2e-2", "2e-2-5e-2", "5e-2-1e-1", "1e-1-1.0 "};
  std::printf("fails / total (worst |dt| in ns) per mu band x h bin:\n");
  std::printf("%-11s", "h[m]");
  for (int b = 0; b < 9; ++b) std::printf(" %-13s", bl[b]);
  std::printf("\n");
  for (int ih = 0; ih < 10; ++ih) {
    std::printf("%-11.0f", hs[ih]);
    for (int b = 0; b < 9; ++b) {
      char buf[32];
      if (tot[b][ih] == 0) std::snprintf(buf, sizeof(buf), "-");
      else std::snprintf(buf, sizeof(buf), "%d/%-2d w=%.3f", fail[b][ih], tot[b][ih], worst[b][ih]);
      std::printf(" %-13s", buf);
    }
    std::printf("\n");
  }
  std::printf("\nper band: fails/total, worst |dt|, reference uncertainty, verdict\n");
  int gf = 0, gt = 0;
  for (int b = 0; b < 9; ++b) {
    char const* v = (unc[b] > 0.03) ? "NOT CHARACTERISED" : (bandFail[b] == 0 ? "PASS" : "FAIL");
    std::printf("   %-11s %4d/%-4d  worst %.4f ns   ref unc %.5f ns   %s\n", bl[b], bandFail[b],
                bandTot[b], bandWorst[b], unc[b], v);
    gf += bandFail[b];
    gt += bandTot[b];
  }
  std::printf("\noverall: %d of %d rays over 0.1 ns  (pass fraction %.4f)\n", gf, gt,
              gt ? 1.0 - double(gf) / double(gt) : 0.0);

  // per-cell reference uncertainty
  std::printf("\nreference uncertainty per cell (ns):\n");
  std::printf("%-11s", "h[m]");
  for (int b = 0; b < 9; ++b) std::printf(" %-9s", bl[b]);
  std::printf("\n");
  for (int ih = 0; ih < 10; ++ih) {
    std::printf("%-11.0f", hs[ih]);
    for (int b = 0; b < 9; ++b) {
      char buf[24];
      if (tot[b][ih] == 0) std::snprintf(buf, sizeof(buf), "-");
      else std::snprintf(buf, sizeof(buf), "%.5f", cellUnc[b][ih]);
      std::printf(" %-9s", buf);
    }
    std::printf("\n");
  }
  double worstCellUnc = 0.0;
  for (int b = 0; b < 9; ++b)
    for (int ih = 0; ih < 10; ++ih)
      if (cellUnc[b][ih] > worstCellUnc) worstCellUnc = cellUnc[b][ih];
  std::printf("worst single-cell reference uncertainty = %.5f ns (limit 0.03)\n", worstCellUnc);

  std::printf("\n5 worst rays overall (|dt|, mu, h, lambda):\n");
  for (int q = 0; q < 5; ++q)
    std::printf("   %2d. %.4f ns   mu=%.2e  h=%.0f m  lambda=%.0f nm\n", q + 1, topDt[q], topMu[q],
                topH[q], topL[q]);

  std::printf("\n=== end of Task C ===\n");
  return 0;
}
