// stage5_b.cpp -- Stage 5: log-spaced mu x h x lambda scan of the TOTAL L and t.
//
// Stage 3 diagnosed ONE ray (h=46000, lambda=400, xi=0.9999) and showed the error is
// in the stored node values (b), that the build integration is first order (p = 1.00),
// and that a second defect appears near the ray's turning point.  Stage 5 maps how
// widely that extends.
//
// For every (h, lambda, mu) the program reports
//     L_trace, L_ref, dL        (geometric path)
//     Sigma_trace, Sigma_ref, dSigma
//     t_trace, t_ref, |dt|      (the quantity the 0.1 ns requirement applies to)
// plus the table's own status and the nSubBase the build used for that column.
//
// The reference is an independent Bouguer midpoint march of the SAME ray:
//     L_ref = int dz / sqrt(1 - s^2),  Sigma_ref = int N dz / sqrt(1 - s^2),
//     s = p / (n(z) z),   n = 1 + N(z) f,
// which matches the build's own dL/dh = 1/cos(theta), dSigma/dh = N/cos(theta).
//
// Scan-only: no production source is touched, no table is rebuilt, no population run.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage5_b.cpp -o stage5_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0;
static constexpr double C_LIGHT = 2.99792458e8; // m/s

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

// independent reference for BOTH L and Sigma along the same ray
//
// NOTE: the archived integrators (path_b.cpp:37-44, combined_b.cpp) use the bare
// form `while (z < zTop - 1e-12) { st = min(dr, zTop - z); ...; z += st; }`.
// That form can STALL: if the remaining gap lands in (1e-12, ULP(z)/2), adding it
// to z leaves z unchanged and the loop never exits.  It was hit during Stage 5 at
// (h=100 m, mu=3e-4).  The guard below breaks out instead; the residual gap is
// below one ULP (~1e-9 m) and cannot affect any reported value.  Stage 3's runs
// did not encounter the window, so its numbers are unaffected.
static void refLS(RefractiveIndexTable const& ri, double p, double f, double hTop, double dr,
                  double& L, double& S, bool* turned) {
  L = 0.0;
  S = 0.0;
  if (turned) *turned = false;
  double z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zn = z + st;
    if (!(zn > z)) break; // ULP stall guard (see note above)
    double const zm = z + 0.5 * st;
    double const nz = nPhase(ri, zm, f);
    double const s = p / (nz * zm);
    if (s >= 1.0) {
      if (turned) *turned = true;
      return;
    }
    double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s));
    L += inv * st;
    S += ri.N(zm) * inv * st;
    z = zn;
  }
}

// the build's own substep rule, so the scan can attribute each ray to a plateau
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

// ray construction identical to path_b.cpp: p = xi * bCrit, mu = sqrt(1 - xi^2)
struct Ray {
  double p, sinT, cosT, mu;
};
static Ray makeRay(RefractiveIndexTable const& ri, double hE, double f, double mu) {
  double const re = R + hE;
  double const pCrit = nPhase(ri, R, f) * R;
  double const xi = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  double const p = xi * pCrit;
  double const nE = nPhase(ri, re, f);
  double const sinT = p / (nE * re);
  double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
  return {p, sinT, cosT, mu};
}

static char const* statusName(int s) {
  switch (s) {
    case 0: return "Hit";
    case 1: return "InvalidTable";
    case 2: return "InvalidInput";
    case 3: return "EmissionBelowObs";
    case 4: return "AboveTable";
    case 5: return "Upgoing";
    case 6: return "BeyondCritical";
    case 7: return "DispersionOutOfRange";
    case 8: return "Grazing";
    default: return "?";
  }
}

struct Rec {
  double h, lam, mu, dt, dL, dS, Ltr, Lref, Str, Sref;
  int nsub, status;
  bool turned;
};

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 5: log-spaced mu x h x lambda scan of total L and t ===\n");
  std::printf("=== scan-only; no production change, no table rebuild ===\n\n");
  std::printf("table: valid=%d hNodes=%zu muNodes=%zu memory=%.1f MB  muMin=%.3e muBreak=%.9f\n",
              int(tab.valid()), tab.hNodes(), tab.muNodes(), tab.memoryMB(), tab.muMin(),
              tab.muBreak());
  std::printf("requirement under test: |dt| <= 0.1 ns over the eligible domain (mu >= 2e-4)\n");
  if (!tab.valid()) return 1;

  double const hs[10] = {100.0, 300.0, 1000.0, 3000.0, 5000.0, 10000.0, 20000.0, 32000.0, 46000.0, 59900.0};
  double const mus[15] = {1.0e-4, 1.5e-4, 2.0e-4, 3.0e-4, 5.0e-4, 1.0e-3, 2.0e-3,
                          5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1, 0.3, 0.6, 1.0};
  double const lams[5] = {350.0, 400.0, 500.0, 690.0, 700.0};
  double const dr = 0.1;

  std::vector<Rec> recs;
  int nRej[9] = {0};
  int nRay = 0;
  std::printf("\n--- per-ray scan (dr=%.2f m for the reference) ---\n", dr);
  std::printf("%-9s %-6s %-10s %-5s %-6s %-15s %-15s %-11s %-11s %s\n", "h[m]", "lam", "mu", "nSub",
              "stat", "L_trace[m]", "L_ref[m]", "dL[m]", "dt[ns]", "flag");
  for (int ih = 0; ih < 10; ++ih) {
    for (int il = 0; il < 5; ++il) {
      for (int im = 0; im < 15; ++im) {
        double const hE = hs[ih], lam = lams[il], mu = mus[im];
        double const f = Fdisp(lam);
        Ray const ry = makeRay(ri, hE, f, mu);
        double const re = R + hE;
        CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, ry.sinT, 0.0, -ry.cosT, f, f);
        ++nRay;
        int const st = int(r.status);
        if (st < 0 || st > 8) continue;
        ++nRej[st];
        if (r.status != CurvedRayTable::Status::Hit) {
          std::printf("%-9.0f %-6.0f %-10.2e %-5d %-6s (rejected)\n", hE, lam, mu, nSubOf(ri, mu, f),
                      statusName(st));
          continue;
        }
        double Lref = 0.0, Sref = 0.0;
        bool turned = false;
        refLS(ri, ry.p, f, hE, dr, Lref, Sref, &turned);
        double const Ltr = r.geometric_path_m, Str = r.sigma_m;
        double const ttr = r.time_s, tref = (Lref + f * Sref) / C_LIGHT;
        double const dt = (ttr - tref) * 1e9;
        double const dL = Ltr - Lref, dS = Str - Sref;
        std::printf("%-9.0f %-6.0f %-10.2e %-5d %-6s %-15.4f %-15.4f %+-11.4f %+-11.4f %s\n", hE, lam,
                    mu, nSubOf(ri, mu, f), statusName(st), Ltr, Lref, dL, dt,
                    turned ? "(ref turned early)" : "");
        Rec rc;
        rc.h = hE;
        rc.lam = lam;
        rc.mu = mu;
        rc.dt = dt;
        rc.dL = dL;
        rc.dS = dS;
        rc.Ltr = Ltr;
        rc.Lref = Lref;
        rc.Str = Str;
        rc.Sref = Sref;
        rc.nsub = nSubOf(ri, mu, f);
        rc.status = st;
        rc.turned = turned;
        recs.push_back(rc);
      }
    }
  }
  std::printf("rays attempted: %d, hits: %zu, rejected: %d\n", nRay, recs.size(), nRay - int(recs.size()));
  for (int s = 0; s < 9; ++s)
    if (nRej[s]) std::printf("   rejected status %-20s %d\n", statusName(s), nRej[s]);
  std::printf("   (mu = 1.0e-4 and 1.5e-4 are BELOW muMin=2.0e-4 by construction: the\n");
  std::printf("    contract is to fall back to the exact tracer, not a table failure)\n");

  // ---- summaries -----------------------------------------------------------
  int nHit = 0, nOver = 0, nElig = 0, nEligOver = 0;
  double worst = 0.0;
  Rec wr = recs.empty() ? Rec{} : recs[0];
  for (Rec const& rc : recs) {
    ++nHit;
    if (std::fabs(rc.dt) > std::fabs(worst)) {
      worst = rc.dt;
      wr = rc;
    }
    if (std::fabs(rc.dt) > 0.1) ++nOver;
    if (rc.mu >= 2.0e-4) {
      ++nElig;
      if (std::fabs(rc.dt) > 0.1) ++nEligOver;
    }
  }
  std::printf("\n=== SUMMARY ===\n");
  std::printf("hits compared: %d\n", nHit);
  std::printf("eligible (mu >= 2e-4): %d,  of which |dt| > 0.1 ns: %d  (%.1f %%)\n", nElig, nEligOver,
              nElig ? 100.0 * nEligOver / nElig : 0.0);
  std::printf("worst |dt| overall = %+.4f ns at h=%.0f m  lambda=%.0f nm  mu=%.4e  (nSubBase=%d)\n",
              wr.dt, wr.h, wr.lam, wr.mu, wr.nsub);
  std::printf("   at that ray: dL=%+.4f m  dSigma=%+.6f m  L_ref=%.4f m  Sigma_ref=%.4f m\n", wr.dL,
              wr.dS, wr.Lref, wr.Sref);

  double const mEd[11] = {0.0, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1, 1.01};
  std::printf("\nworst |dt| by mu band (all lambda, all h):\n");
  std::printf("   %-24s %-13s %-11s %s\n", "mu band", "worst|dt|[ns]", "at mu", "over 0.1ns / hits");
  for (int b = 0; b + 1 < 11; ++b) {
    double w = 0.0, wmu = 0.0;
    int over = 0, n = 0;
    for (Rec const& rc : recs)
      if (rc.mu >= mEd[b] && rc.mu < mEd[b + 1]) {
        ++n;
        if (std::fabs(rc.dt) > 0.1) ++over;
        if (std::fabs(rc.dt) > std::fabs(w)) {
          w = rc.dt;
          wmu = rc.mu;
        }
      }
    if (n)
      std::printf("   [%-9.1e,%-9.1e) %-13.4f %-11.2e %d / %d\n", mEd[b], mEd[b + 1], w, wmu, over, n);
  }

  std::printf("\nworst |dt| by h (all lambda, mu >= 2e-4):\n");
  std::printf("   %-10s %-13s %-11s %-7s %s\n", "h[m]", "worst|dt|[ns]", "at mu", "lam", "over 0.1ns / hits");
  for (int ih = 0; ih < 10; ++ih) {
    double w = 0.0, wmu = 0.0, wlam = 0.0;
    int over = 0, n = 0;
    for (Rec const& rc : recs)
      if (rc.h == hs[ih] && rc.mu >= 2.0e-4) {
        ++n;
        if (std::fabs(rc.dt) > 0.1) ++over;
        if (std::fabs(rc.dt) > std::fabs(w)) {
          w = rc.dt;
          wmu = rc.mu;
          wlam = rc.lam;
        }
      }
    if (n) std::printf("   %-10.0f %-13.4f %-11.2e %-7.0f %d / %d\n", hs[ih], w, wmu, wlam, over, n);
  }

  std::printf("\nworst |dt| by lambda (all h, mu >= 2e-4):\n");
  for (int il = 0; il < 5; ++il) {
    double w = 0.0, wmu = 0.0, wh = 0.0;
    int over = 0, n = 0;
    for (Rec const& rc : recs)
      if (rc.lam == lams[il] && rc.mu >= 2.0e-4) {
        ++n;
        if (std::fabs(rc.dt) > 0.1) ++over;
        if (std::fabs(rc.dt) > std::fabs(w)) {
          w = rc.dt;
          wmu = rc.mu;
          wh = rc.h;
        }
      }
    if (n) std::printf("   lam=%-6.0f %-13.4f %-11.2e (h=%-8.0f) %d / %d\n", lams[il], w, wmu, wh, over, n);
  }

  int nTurn = 0;
  double wTurn = 0.0;
  for (Rec const& rc : recs)
    if (rc.turned) {
      ++nTurn;
      if (std::fabs(rc.dt) > std::fabs(wTurn)) wTurn = rc.dt;
    }
  std::printf("\nrays where the reference turned before hTop (turning-point region): %d\n", nTurn);
  std::printf("   worst |dt| among them = %+.4f ns\n", wTurn);

  std::printf("\nreference self-convergence on the WORST ray of each mu band (dr=%.2f -> dr=%.3f):\n",
              dr, dr / 5.0);
  std::printf("   %-23s %-9s %-6s %-11s %-13s %-13s %-13s %s\n", "mu band", "h[m]", "lam", "mu",
              "worst|dt|[ns]", "dL_ref[mm]", "d(dt)[ns]", "verdict");
  for (int b = 0; b + 1 < 11; ++b) {
    Rec const* best = nullptr;
    for (Rec const& rc : recs)
      if (rc.mu >= mEd[b] && rc.mu < mEd[b + 1])
        if (!best || std::fabs(rc.dt) > std::fabs(best->dt)) best = &rc;
    if (!best) continue;
    double const f = Fdisp(best->lam);
    Ray const ry = makeRay(ri, best->h, f, best->mu);
    double L2 = 0.0, S2 = 0.0;
    bool tn = false;
    refLS(ri, ry.p, f, best->h, dr / 5.0, L2, S2, &tn);
    double const dLmm = (L2 - best->Lref) * 1e3;
    double const ddt = ((L2 + f * S2) - (best->Lref + f * best->Sref)) / C_LIGHT * 1e9;
    double const refShift = std::fabs(ddt);
    char const* verdict = (refShift > 0.1)     ? "REFERENCE-LIMITED"
                          : (refShift > 0.01) ? "ref ~10x margin"
                                              : "ref converged";
    std::printf("   [%-9.1e,%-9.1e) %-9.0f %-6.0f %-11.2e %-13.4f %-13.3f %-13.6f %s\n", mEd[b],
                mEd[b + 1], best->h, best->lam, best->mu, best->dt, dLmm, ddt, verdict);
  }
  std::printf("\ninterpretation: a band is only readable where the reference has converged,\n");
  std::printf("   i.e. |d(dt)| of the reference refinement is well below the reported |dt|.\n");

  std::printf("\n=== end of Stage 5 ===\n");
  return 0;
}
