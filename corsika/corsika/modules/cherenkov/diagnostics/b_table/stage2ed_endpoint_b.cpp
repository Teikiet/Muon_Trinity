// ============================================================================
// stage2ed_endpoint_b.cpp -- Checkpoint 2e-d: how much of `db` comes from the
// diagnostic's final step extending below the ground.
//
//   ./stage2ed_endpoint_b <profile .dat>
//
// Diagnostic only. Same three reference rays, same physical inputs, same
// compiler settings. Nothing here touches production, the lookup table, the
// atmosphere, or stage9b_abs_b.cpp.
//
//   db = tauIntS(400 m, mode 0) - tauIntS(1 m, mode 0)          (stage9b:422)
//
// ENDPOINT-TREATMENT EXPERIMENT (NOT a physically validated reference):
// all steps before the final one are kept bit-identical; only the final
// ground-crossing step is re-run at a shortened length, located by bisection
// with the SAME step function, and its tau uses the SAME quadrature rule with
// the sample position recomputed for the shortened step.
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include "tau_march_b.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
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
static double bCrit(RefractiveIndexTable const& ri, double f) {
  return (1.0 + ri.N(R) * f) * R;
}

struct Profile {
  std::vector<double> lam, hgt;
  std::vector<std::vector<double> > tau;
};
static bool readProfile(std::string const& fn, Profile& p) {
  std::ifstream f(fn.c_str());
  if (!f) return false;
  int nl = 0, nh = 0;
  f >> nl >> nh;
  if (nl <= 0 || nh <= 0) return false;
  p.lam.resize(nl);
  p.hgt.resize(nh);
  for (int i = 0; i < nl; ++i) f >> p.lam[i];
  for (int i = 0; i < nh; ++i) f >> p.hgt[i];
  p.tau.assign(nl, std::vector<double>(nh, 0.0));
  for (int i = 0; i < nl; ++i)
    for (int j = 0; j < nh; ++j) f >> p.tau[i][j];
  return true;
}
static int nearestLam(Profile const& p, double lam) {
  int best = 0;
  double bd = 1e30;
  for (std::size_t i = 0; i < p.lam.size(); ++i) {
    double const d = std::fabs(p.lam[i] - lam);
    if (d < bd) {
      bd = d;
      best = int(i);
    }
  }
  return best;
}
static double alphaRef(Profile const& p, int il, double h) {
  std::vector<double> const& H = p.hgt;
  std::vector<double> const& T = p.tau[il];
  if (h <= H.front()) return (T[1] - T[0]) / (H[1] - H[0]);
  if (h >= H.back())
    return (T[H.size() - 1] - T[H.size() - 2]) / (H.back() - H[H.size() - 2]);
  std::size_t j = 0;
  while (j + 2 < H.size() && H[j + 1] <= h) ++j;
  return (T[j + 1] - T[j]) / (H[j + 1] - H[j]);
}
// production rule (Cherenkov.inl:96-113), via the production class
static double alphaProd(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                        double hHi) {
  double const hl = std::max(hLo, h - 25.0), hh = std::min(hHi, h + 25.0);
  double const adh = hh - hl;
  if (adh <= 0.0) return 0.0;
  double const t1 = abs.getOpticalDepth(lam, hh);
  double const t0 = abs.getOpticalDepth(lam, hl);
  return std::max(0.0, (t1 - t0) / adh);
}
// --- stage9b_abs_b.cpp:166-202, VERBATIM (the unchanged baseline) ------------
static double tauIntS(AtmosphericAbsorption const& abs, Profile const& prof,
                      RefractiveIndexTable const& ri, double lam, int ilp, double p, double f,
                      double hTop, double ds, int mode, double hLo, double hHi, long* nSteps,
                      std::vector<double>* rec = nullptr) {
  double z = R + hTop;
  double const nTop = 1.0 + ri.N(z) * f;
  double const sinT = std::min(1.0, p / (nTop * z));
  double th = std::acos(-std::sqrt(std::max(0.0, 1.0 - sinT * sinT))); // downgoing
  double tau = 0.0;
  long n = 0;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const nn = 1.0 + ri.N(rr) * f;
    double const np = f * ri.dNdz(rr);
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c;
    dt = -s * (1.0 / rr + (nn > 1e-12 ? np / nn : 0.0));
  };
  while (z > R + 1e-9 && n < 200000000L) {
    double dr1, dt1, dr2, dt2, dr3, dt3, dr4, dt4;
    der(z, th, dr1, dt1);
    der(z + 0.5 * ds * dr1, th + 0.5 * ds * dt1, dr2, dt2);
    der(z + 0.5 * ds * dr2, th + 0.5 * ds * dt2, dr3, dt3);
    der(z + ds * dr3, th + ds * dt3, dr4, dt4);
    double const zPrev = z;
    z += ds / 6.0 * (dr1 + 2 * dr2 + 2 * dr3 + dr4);
    th += ds / 6.0 * (dt1 + 2 * dt2 + 2 * dt3 + dt4);
    double const h = std::max(0.0, 0.5 * (zPrev + z) - R);
    tau += (mode == 0 ? alphaProd(abs, lam, h, hLo, hHi) : alphaRef(prof, ilp, h)) * ds;
    if (rec) {
      rec->push_back(h);
      rec->push_back(ds);
    }
    ++n;
  }
  if (nSteps) *nSteps = n;
  return tau;
}


// ============================================================================
// gate plumbing
// ============================================================================
static int g_gates = 0, g_gates_fail = 0;

static void gate(char const* label, double residual, double tol) {
  ++g_gates;
  bool const ok = (residual <= tol);
  if (!ok) ++g_gates_fail;
  std::printf("   %-46s residual=%.17g  tol=%.17g  %s\n", label, residual, tol,
              ok ? "PASS" : "FAIL");
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <profile .dat>\n", argv[0]);
    return 2;
  }
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[1]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  Profile prof;
  if (!readProfile(argv[1], prof)) {
    std::printf("STOP: cannot re-read the native knots\n");
    return 1;
  }
  auto const hr = abs.getHeightRange();
  RefractiveIndexTable ri = makeRI();

  struct C {
    double lam, mu, h;
  };
  C const cs[3] = {{300.0, 0.6, 46000.0}, {480.0, 0.1, 1000.0}, {480.0, 0.1, 975.0}};
  double const dsBig = 8000.0 / 20.0;  // 400 m, as stage9b_abs_b.cpp:404
  double const dsOne = 1.0;

  std::printf("=== Checkpoint 2e-d: endpoint contribution to db ===\n");
  std::printf("profile            : %s\n", argv[1]);
  std::printf("R                  : %.17g m ;  table height range [%.17g, %.17g] m\n", R,
              hr.first, hr.second);
  std::printf("step sizes         : 400 m = %.17g m , 1 m = %.17g m\n", dsBig, dsOne);
  std::printf("extinction         : alphaProd (Cherenkov.inl:96-113, dh = 50 m, clamped)\n");
  std::printf("experiment         : ENDPOINT-TREATMENT EXPERIMENT, not a validated reference\n\n");

  // ---------------------------------------------------------------- section 3
  // Reproduction of the unchanged baseline via the verbatim copy.
  std::printf("--- section 3: baseline reproduction (2e-c full-precision values) ---\n");
  double tRuleRef[3], tProdRef[3], dbRef[3];
  long nRuleRef[3], nProdRef[3];
  {
    // 2e-c values, read from cp2ea_audit.out
    double const expRule[3] = {4.5627185999106041, 1.0360241577028855, 1.0132745869998552};
    double const expProd[3] = {4.6390883962468301, 1.0423360298413382, 1.0442899999999999};
    double const expDb[3] = {0.076369796336225981, 0.0063118721384527454,
                             0.031015413000144765};
    long const expN[3] = {76217, 9951, 9703};
    long const expNP[3] = {191, 25, 25};
    for (int q = 0; q < 3; ++q) {
      double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
      int const ilq = nearestLam(prof, lm);
      double const ff = Fdisp(lm);
      double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
      long n2 = 0, n3 = 0;
      tRuleRef[q] = tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsOne, 0, hr.first, hr.second, &n2);
      tProdRef[q] = tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsBig, 0, hr.first, hr.second, &n3);
      dbRef[q] = tProdRef[q] - tRuleRef[q];
      nRuleRef[q] = n2;
      nProdRef[q] = n3;
      std::printf("ray %d: tau_rule1m=%.17g  tau_prod_curved=%.17g  db=%.17g  n=%ld/%ld\n",
                  q, tRuleRef[q], tProdRef[q], dbRef[q], n2, n3);
      gate("3 tau_rule1m reproduces 2e-c", std::fabs(tRuleRef[q] - expRule[q]),
           1e-12 * std::max(1.0, std::fabs(expRule[q])));
      gate("3 tau_prod_curved reproduces 2e-c", std::fabs(tProdRef[q] - expProd[q]),
           1e-12 * std::max(1.0, std::fabs(expProd[q])));
      gate("3 db reproduces 2e-c", std::fabs(dbRef[q] - expDb[q]),
           1e-12 * std::max(1.0, std::fabs(expDb[q])));
      gate("3 step count (1 m)", std::fabs(double(n2 - expN[q])),
           1e-12 * std::max(1.0, std::fabs(double(expN[q]))));
      gate("3 step count (400 m)", std::fabs(double(n3 - expNP[q])),
           1e-12 * std::max(1.0, std::fabs(double(expNP[q]))));
    }
  }
  if (g_gates_fail != 0) {
    std::printf("\nSTOP: baseline reproduction FAILED (%d gate(s)) -- not continuing.\n",
                g_gates_fail);
    return 1;
  }
  std::printf("baseline OK (all gates pass)\n\n");


  // ------------------------------------------------------------ sections 4-6
  double tauFull[3][2] = {{0}}, tauClip[3][2] = {{0}};
  double kConst = 1e-4;  // section 6A test coefficient [1/m]

  for (int q = 0; q < 3; ++q) {
    double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
    double const ff = Fdisp(lm);
    double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);

    // the unchanged diagnostic's propagation + extinction callables
    refmarch::NFn Nof = [&](double rr) { return ri.N(rr) * ff; };
    refmarch::DerFn der = [&](double rr, double tt, double& dr, double& dt) {
      double const nn = 1.0 + ri.N(rr) * ff;
      double const np = ff * ri.dNdz(rr);
      double const s = std::sin(tt), c = std::cos(tt);
      dr = c;
      dt = -s * (1.0 / rr + (nn > 1e-12 ? np / nn : 0.0));
    };
    refmarch::AlphaFn alpha = [&](double h) {
      return alphaProd(abs, lm, h, hr.first, hr.second);
    };
    refmarch::AlphaFn alphaK = [&](double) { return kConst; };

    for (int wi = 0; wi < 2; ++wi) {
      double const ds = (wi == 0) ? dsBig : dsOne;
      std::printf("--- ray %d  lambda=%.17g nm  mu=%.17g  h=%.17g m  ds=%.17g m ---\n", q, lm,
                  mu, hh, ds);

      refmarch::MarchResult const full = refmarch::march(R, hh, ds, p, Nof, der, alpha, false);
      refmarch::MarchResult const clip = refmarch::march(R, hh, ds, p, Nof, der, alpha, true);
      refmarch::FinalStepRecord const& rf = full.fin;
      refmarch::FinalStepRecord const& rc = clip.fin;

      // helper must reproduce the verbatim baseline bit for bit
      double const base = (wi == 0) ? tProdRef[q] : tRuleRef[q];
      double const baseN = double((wi == 0) ? nProdRef[q] : nRuleRef[q]);
      gate("helper(FULL) == tauIntS baseline, tau", std::fabs(full.tau - base), 0.0);
      gate("helper(FULL) == tauIntS baseline, steps", std::fabs(double(full.n) - baseN), 0.0);
      if (!rf.valid) {
        std::printf("STOP: no ground-crossing final step was reached (step limit).\n");
        return 1;
      }

      std::printf("   steps=%ld  nominal path=%.17g m  tau_before_final=%.17g\n", rf.nSteps,
                  double(rf.nSteps) * ds, rf.tauBefore);
      std::printf("   tau_total_full=%.17g  tau_total_clipped=%.17g\n", full.tau, clip.tau);
      std::printf("   final step index=%ld\n", rf.stepIndex);
      std::printf("   pre  : z=%.17g  h=%.17g  th=%.17g\n", rf.zPre, rf.hPre, rf.thPre);
      std::printf("   full : z=%.17g  h=%.17g  th=%.17g\n", rf.zPostFull, rf.hPostFull,
                  rf.thPostFull);
      std::printf("   nominal final-step length=%.17g m\n", rf.dsNominal);
      std::printf("   extinction samples in that step: 1   h=%.17g  alpha=%.17g\n",
                  rf.hSampleFull, rf.alphaSampleFull);
      std::printf("   final-step tau (full)=%.17g\n", rf.dTauFull);
      std::printf("   brackets r=R: h_pre>0:%s  h_post<0:%s\n",
                  (rf.hPre > 0.0) ? "yes" : "NO", (rf.hPostFull < 0.0) ? "yes" : "NO");
      if (!(rf.hPre > 0.0 && rf.hPostFull < 0.0)) {
        std::printf("   STOP: endpoints do not bracket r=R; termination is the loop condition\n");
        std::printf("         while (z > R + 1e-9) of stage9b_abs_b.cpp:183, not a crossing.\n");
        return 1;
      }
      std::printf("   clipped: bracketed=%s monotonic_down=%s worst_backward=%.17g m\n",
                  rc.bracketed ? "yes" : "NO", rc.monotonicDown ? "yes" : "NO",
                  rc.worstBackward);
      if (!rc.bracketed || !rc.monotonicDown) {
        std::printf("   STOP: single downward crossing within the final step not established.\n");
        return 1;
      }
      std::printf("   clipped: bisect_iters=%ld  bracket_width=%.17g m  ds_hit=%.17g m\n",
                  rc.bisectIters, rc.bracketWidth, rc.dsHit);
      std::printf("   clipped: z=%.17g  h=%.17g  th=%.17g  |z-R|=%.17g m\n", rc.zPostClip,
                  rc.hPostClip, rc.thPostClip, std::fabs(rc.zPostClip - R));
      std::printf("   clipped: sample h=%.17g  alpha=%.17g  tau=%.17g\n", rc.hSampleClip,
                  rc.alphaSampleClip, rc.dTauClip);
      std::printf("   delta = tau_full_final - tau_clipped_final = %.17g\n",
                  rf.dTauFull - rc.dTauClip);
      gate("5 |r_end - R| <= 1e-5 m", std::fabs(rc.zPostClip - R), 1e-5);
      gate("5 pre-state identical in both runs",
           std::fabs(rc.zPre - rf.zPre) + std::fabs(rc.thPre - rf.thPre), 0.0);


      // ---------- A. constant-extinction check (test calculation only) ------
      // Extinction never enters the propagation, so the trajectory and ds_hit
      // must be unchanged; only the tau increments use k instead of alphaProd.
      refmarch::MarchResult const fullK = refmarch::march(R, hh, ds, p, Nof, der, alphaK, false);
      refmarch::MarchResult const clipK = refmarch::march(R, hh, ds, p, Nof, der, alphaK, true);
      double const dtFullK = fullK.fin.dTauFull;
      double const dtClipK = clipK.fin.dTauClip;
      double const remK = dtFullK - dtClipK;
      double const expFullK = kConst * ds;
      double const expClipK = kConst * clipK.fin.dsHit;
      double const expRemK = kConst * (ds - clipK.fin.dsHit);
      std::printf("   A k=%.17g 1/m: full=%.17g expected=%.17g\n", kConst, dtFullK, expFullK);
      std::printf("                   clipped=%.17g expected=%.17g\n", dtClipK, expClipK);
      std::printf("                   removed(calc)=%.17g expected=%.17g\n", remK, expRemK);
      gate("6A full increment == k*ds", std::fabs(dtFullK - expFullK),
           1e-12 * std::max(1.0, std::fabs(dtFullK)));
      gate("6A clipped increment == k*ds_hit", std::fabs(dtClipK - expClipK),
           1e-12 * std::max(1.0, std::fabs(dtClipK)));
      gate("6A removed == k*(ds - ds_hit)", std::fabs(remK - expRemK),
           1e-12 * std::max(1.0, std::fabs(remK)));
      gate("6A ds_hit unchanged by the coefficient", std::fabs(clipK.fin.dsHit - rc.dsHit), 0.0);

      // ---------- B. back-of-envelope atmospheric check (informational) -----
      double const dsRemoved = ds - rc.dsHit;
      double const aGround = alphaProd(abs, lm, 0.0, hr.first, hr.second);
      double const dtEstimate = aGround * dsRemoved;
      double const dtActual = rf.dTauFull - rc.dTauClip;
      std::printf("   B delta_s=%.17g m  alphaProd(0,lambda)=%.17g 1/m  estimate=%.17g\n",
                  dsRemoved, aGround, dtEstimate);
      std::printf("     actual removed tau=%.17g  estimate-actual=%.17g  (not a gate)\n",
                  dtActual, dtEstimate - dtActual);

      tauFull[q][wi] = full.tau;
      tauClip[q][wi] = clip.tau;
      std::printf("\n");
      std::fflush(stdout);
    }
  }

  // ---------- C. decomposition identity ------------------------------------
  std::printf("--- section 6C: decomposition identity ---\n");
  double d400[3], d1[3], dbClipped[3];
  for (int q = 0; q < 3; ++q) {
    d400[q] = tauFull[q][0] - tauClip[q][0];
    d1[q] = tauFull[q][1] - tauClip[q][1];
    dbClipped[q] = tauClip[q][0] - tauClip[q][1];
    double const scale =
        std::max(std::max(1.0, std::fabs(dbRef[q])),
                 std::max(std::fabs(tauFull[q][0]), std::fabs(tauFull[q][1])));
    std::printf("ray %d: db=%.17g  db_clipped=%.17g  delta_400=%.17g  delta_1=%.17g\n", q,
                dbRef[q], dbClipped[q], d400[q], d1[q]);
    gate("6C db == db_clipped + delta_400 - delta_1",
         std::fabs(dbRef[q] - (dbClipped[q] + d400[q] - d1[q])), 1e-12 * scale);
  }

  std::printf("\n=== summary ===\n");
  std::printf("%-4s %-22s %-22s %-22s %-22s %-22s %-22s\n", "ray", "db_original", "delta_400",
              "delta_1", "endpoint_difference", "db_clipped", "endpoint/db_original");
  for (int q = 0; q < 3; ++q) {
    double const ep = d400[q] - d1[q];
    double const ratio = (dbRef[q] != 0.0) ? ep / dbRef[q] : 0.0;
    std::printf("%-4d %-22.17g %-22.17g %-22.17g %-22.17g %-22.17g %-22.17g\n", q, dbRef[q],
                d400[q], d1[q], ep, dbClipped[q], ratio);
  }
  std::printf("\ngates: %d evaluated, %d failed\n", g_gates, g_gates_fail);
  std::printf("%s\n", (g_gates_fail == 0)
                          ? "ALL GATES PASS"
                          : "STOP: a mandatory numerical gate FAILED -- reported above");
  std::printf("=== end of 2e-d (endpoint-treatment experiment, not a validated reference) ===\n");
  return (g_gates_fail == 0) ? 0 : 1;
}

