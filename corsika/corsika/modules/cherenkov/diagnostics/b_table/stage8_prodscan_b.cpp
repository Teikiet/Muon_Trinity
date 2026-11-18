// stage8_prodscan_b.cpp -- Stage 8, Task 1.4: the Stage 7 grid through the PRODUCTION
// lookup path, after the V4 fix was applied to CurvedRayTable::buildColumn.
//
// Unlike Stage 7 Task C, this goes through trace(): the (h, mu) Hermite interpolation and
// the f-axis interpolation are exercised, so it tests the production table, not a copy of
// the column builder.  The observable is the production Result.time_s, and the ground
// truth is the Stage 6 Task 2 w-reference.
//
// Diagnostics only.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage8_prodscan_b.cpp -o stage8_prodscan_b
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

// Displacement reference: the project defines the transverse displacement as
//   d = r_obs * Phi      (CurvedRayTable.hpp:775 in trace(); the build stores
//                         dOut[i+1] = r_obs_ * phi, same file, buildColumn)
// with dPhi/dh = tan(theta)/r (the build's own RHS, dph = s/(r*cc)).  Along the
// Bouguer ray s = p/(n(z) z) and tan(theta) = s/sqrt(1-s^2), so
//   Phi = int tan(theta)/r dz ,  dz = 2w dw  in the w variable.
static void refPhi(RefractiveIndexTable const& ri, double p, double f, double hTop, double dw,
                   double& Phi) {
  Phi = 0.0;
  double const wTop = std::sqrt(hTop);
  for (double w = 0.0; w < wTop - 1e-12;) {
    double const st = std::min(dw, wTop - w);
    double const wn = w + st;
    if (!(wn > w)) break; // stall guard
    double const wm = w + 0.5 * st;
    double const z = R + wm * wm;
    double const s = p / (nPhase(ri, z, f) * z);
    if (s >= 1.0) return;
    double const tanT = s / std::sqrt(std::max(1e-300, 1.0 - s * s));
    Phi += tanT / z * (2.0 * wm) * st;
    w = wn;
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
  std::printf("=== Stage 8, Task 1.4: Stage 7 grid through the PRODUCTION lookup path ===\n\n");
  std::printf("table: valid=%d hNodes=%zu muNodes=%zu memory=%.1f MB  muBreak=%.9f\n",
              int(tab.valid()), tab.hNodes(), tab.muNodes(), tab.memoryMB(), tab.muBreak());
  if (!tab.valid()) return 1;
  double const hs[10] = {100.0, 300.0, 1000.0, 3000.0, 5000.0,
                         10000.0, 20000.0, 32000.0, 46000.0, 59900.0};
  double const mus[13] = {2.0e-4, 3.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2,
                          2.0e-2, 5.0e-2, 1.0e-1, 0.3, 0.6, 1.0};
  double const lams[5] = {350.0, 400.0, 500.0, 690.0, 700.0};

  int tot[9][10] = {{0}}, fail[9][10] = {{0}};
  double worst[9][10], cellUnc[9][10], unc[9] = {0.0};
  for (int b = 0; b < 9; ++b)
    for (int q = 0; q < 10; ++q) {
      worst[b][q] = 0.0;
      cellUnc[b][q] = 0.0;
    }
  int bandTot[9] = {0}, bandFail[9] = {0};
  double bandWorst[9] = {0.0};
  // displacement (1 cm target) accumulators
  int dFail[9][10] = {{0}};
  double dWorst[9][10], dCellUnc[9][10], dUnc[9] = {0.0};
  for (int b = 0; b < 9; ++b)
    for (int q = 0; q < 10; ++q) {
      dWorst[b][q] = 0.0;
      dCellUnc[b][q] = 0.0;
    }
  int dBandTot[9] = {0}, dBandFail[9] = {0};
  double dBandWorst[9] = {0.0};
  double dTop[5] = {0.0, 0.0, 0.0, 0.0, 0.0}, dTopMu[5] = {0.0}, dTopH[5] = {0.0};
  // rejection listing (Task C)
  int nRejList = 0;
  double rejMu[100] = {0.0}, rejH[100] = {0.0}, rejLam[100] = {0.0};
  int rejSt[100] = {0};
  double topDt[5] = {0.0, 0.0, 0.0, 0.0, 0.0}, topMu[5] = {0.0}, topH[5] = {0.0}, topL[5] = {0.0};
  int nHit = 0, nRej = 0, rejStat[9] = {0};

  for (int il = 0; il < 5; ++il) {
    double const f = Fdisp(lams[il]);
    for (int im = 0; im < 13; ++im) {
      double const mu = mus[im];
      double const xi = std::sqrt(std::max(0.0, 1.0 - mu * mu));
      double const p = xi * bCrit(ri, f);
      int const b = bandOf(mu);
      for (int ih = 0; ih < 10; ++ih) {
        double const h = hs[ih];
        double const re = R + h;
        double const nE = nPhase(ri, re, f);
        double const sinT = p / (nE * re);
        if (sinT >= 1.0) continue;
        double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
        CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
        int const st = int(r.status);
        if (r.status != CurvedRayTable::Status::Hit) {
          ++nRej;
          if (st >= 0 && st <= 8) ++rejStat[st];
          if (nRejList < 100) {
            rejMu[nRejList] = mu;
            rejH[nRejList] = h;
            rejLam[nRejList] = lams[il];
            rejSt[nRejList] = st;
            ++nRejList;
          }
          continue;
        }
        ++nHit;
        double Lr = 0.0, Sr = 0.0, Lr2 = 0.0, Sr2 = 0.0;
        refW(ri, p, f, h, 1.0e-3, Lr, Sr);
        refW(ri, p, f, h, 2.0e-3, Lr2, Sr2);
        double const tP = r.time_s; // the production observable
        double const tR = (Lr + f * Sr) / C_LIGHT;
        double const tR2 = (Lr2 + f * Sr2) / C_LIGHT;
        double const dt = std::fabs(tP - tR) * 1e9;
        double const u = std::fabs(tR2 - tR) * 1e9;
        ++tot[b][ih];
        if (dt > 0.1) ++fail[b][ih];
        if (dt > worst[b][ih]) worst[b][ih] = dt;
        if (u > cellUnc[b][ih]) cellUnc[b][ih] = u;
        if (u > unc[b]) unc[b] = u;
        ++bandTot[b];
        if (dt > 0.1) ++bandFail[b];
        if (dt > bandWorst[b]) bandWorst[b] = dt;
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

        // displacement: d = r_obs * Phi, Phi being the angle between the arrival
        // radial (res.x,res.y,res.z)/r_obs and the emission radial (0,0,1)
        double const dTab = R * std::acos(std::min(1.0, std::max(-1.0, r.z / R)));
        double Phi1 = 0.0, Phi2 = 0.0;
        refPhi(ri, p, f, h, 1.0e-3, Phi1);
        refPhi(ri, p, f, h, 2.0e-3, Phi2);
        double const dRef = R * Phi1, dRef2 = R * Phi2;
        double const de = std::fabs(dTab - dRef);
        double const du = std::fabs(dRef2 - dRef);
        if (de > 0.01) ++dFail[b][ih];
        if (de > dWorst[b][ih]) dWorst[b][ih] = de;
        if (du > dCellUnc[b][ih]) dCellUnc[b][ih] = du;
        if (du > dUnc[b]) dUnc[b] = du;
        ++dBandTot[b];
        if (de > 0.01) ++dBandFail[b];
        if (de > dBandWorst[b]) dBandWorst[b] = de;
        if (de > dTop[4]) {
          int q = 4;
          while (q > 0 && de > dTop[q - 1]) {
            dTop[q] = dTop[q - 1];
            dTopMu[q] = dTopMu[q - 1];
            dTopH[q] = dTopH[q - 1];
            --q;
          }
          dTop[q] = de;
          dTopMu[q] = mu;
          dTopH[q] = h;
        }
      }
    }
  }

  char const* bl[9] = {"2e-4-5e-4", "5e-4-1e-3", "1e-3-2e-3", "2e-3-5e-3", "5e-3-1e-2",
                       "1e-2-2e-2", "2e-2-5e-2", "5e-2-1e-1", "1e-1-1.0 "};
  std::printf("\nfails / total (worst |dt| ns) per mu band x h bin, PRODUCTION path:\n");
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
  std::printf("\nper band: fails/total, worst |dt|, ref uncertainty, verdict\n");
  int gf = 0, gt = 0;
  for (int b = 0; b < 9; ++b) {
    char const* v = (unc[b] > 0.03) ? "NOT CHARACTERISED" : (bandFail[b] == 0 ? "PASS" : "FAIL");
    std::printf("   %-11s %4d/%-4d  worst %.4f ns   ref unc %.5f ns   %s\n", bl[b], bandFail[b],
                bandTot[b], bandWorst[b], unc[b], v);
    gf += bandFail[b];
    gt += bandTot[b];
  }
  std::printf("\nproduction hits %d, rejections %d (", nHit, nRej);
  for (int s = 0; s < 9; ++s)
    if (rejStat[s]) std::printf("status %d: %d  ", s, rejStat[s]);
  std::printf(")\n");
  std::printf("overall: %d of %d production rays over 0.1 ns  (pass fraction %.4f)\n", gf, gt,
              gt ? 1.0 - double(gf) / double(gt) : 0.0);
  double worstCellUnc = 0.0;
  for (int b = 0; b < 9; ++b)
    for (int ih = 0; ih < 10; ++ih)
      if (cellUnc[b][ih] > worstCellUnc) worstCellUnc = cellUnc[b][ih];
  std::printf("worst single-cell reference uncertainty = %.5f ns (limit 0.03)\n", worstCellUnc);
  std::printf("\n5 worst production rays (|dt|, mu, h, lambda):\n");
  for (int q = 0; q < 5; ++q)
    std::printf("   %2d. %.4f ns   mu=%.2e  h=%.0f m  lambda=%.0f nm\n", q + 1, topDt[q], topMu[q],
                topH[q], topL[q]);

  // ---- Task B: displacement against the 1 cm target ------------------------
  std::printf("\nDISPLACEMENT (target 1 cm): fails / total (worst |dd| mm) per mu band x h bin\n");
  std::printf("%-11s", "h[m]");
  for (int b = 0; b < 9; ++b) std::printf(" %-13s", bl[b]);
  std::printf("\n");
  for (int ih = 0; ih < 10; ++ih) {
    std::printf("%-11.0f", hs[ih]);
    for (int b = 0; b < 9; ++b) {
      char buf[32];
      if (tot[b][ih] == 0) std::snprintf(buf, sizeof(buf), "-");
      else std::snprintf(buf, sizeof(buf), "%d/%-2d w=%.2f", dFail[b][ih], tot[b][ih],
                         dWorst[b][ih] * 1e3);
      std::printf(" %-13s", buf);
    }
    std::printf("\n");
  }
  std::printf("\nper band: fails/total, worst |dd|, ref uncertainty, verdict (unc limit 3 mm)\n");
  int df = 0, dt2 = 0;
  for (int b = 0; b < 9; ++b) {
    char const* v = (dUnc[b] > 3.0e-3) ? "NOT CHARACTERISED"
                                       : (dBandFail[b] == 0 ? "PASS" : "FAIL");
    std::printf("   %-11s %4d/%-4d  worst %.3f mm   ref unc %.4f mm   %s\n", bl[b], dBandFail[b],
                dBandTot[b], dBandWorst[b] * 1e3, dUnc[b] * 1e3, v);
    df += dBandFail[b];
    dt2 += dBandTot[b];
  }
  std::printf("overall displacement: %d of %d rays over 1 cm  (pass fraction %.4f)\n", df, dt2,
              dt2 ? 1.0 - double(df) / double(dt2) : 0.0);
  std::printf("5 worst displacement rays (|dd| mm, mu, h):\n");
  for (int q = 0; q < 5; ++q)
    std::printf("   %2d. %.4f mm   mu=%.2e  h=%.0f m\n", q + 1, dTop[q] * 1e3, dTopMu[q], dTopH[q]);

  // ---- Task C: the rejection listing ---------------------------------------
  std::printf("\nrejections (%d listed of %d):\n", nRejList, nRej);
  for (int q = 0; q < nRejList; ++q)
    std::printf("   mu=%.2e  h=%-7.0f lambda=%-5.0f status=%d\n", rejMu[q], rejH[q], rejLam[q],
                rejSt[q]);

  std::printf("\n=== end of Task 1.4 ===\n");
  return 0;
}
