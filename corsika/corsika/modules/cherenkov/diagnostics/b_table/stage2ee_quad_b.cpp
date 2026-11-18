// ============================================================================
// stage2ee_quad_b.cpp -- Checkpoint 2e-e: isolate extinction quadrature on the
// existing (endpoint-clipped) diagnostic paths.
//
//   ./stage2ee_quad_b <profile .dat>
//
// Diagnostic experiment, NOT production validation. It does not change the
// propagated path: every accepted step of the 2e-d clipped march keeps its
// pre-step state, its accepted length and its post-step state, and interior
// points are RK4Step(y_i, u) evaluated from the ORIGINAL pre-step state.
// Only the extinction integral is recomputed, with composite midpoints on
// panels of at most q = 1.0 m and q = 0.5 m.
//
// Decomposes (T400 - T1) = (C400 - C1) + DQ, with C = T - Q.
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include "tau_march_b.hpp"
#include "tau_path_b.hpp"

#include <algorithm>
#include <chrono>
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

static int g_gates = 0, g_gates_fail = 0;
static void gate(char const* label, double residual, double tol) {
  ++g_gates;
  bool const ok = (residual <= tol);
  if (!ok) ++g_gates_fail;
  std::printf("   %-52s residual=%.17g  tol=%.17g  %s\n", label, residual, tol,
              ok ? "PASS" : "FAIL");
}

static double sumLengths(refpath::Path const& p) {
  double s = 0.0;
  for (std::size_t i = 0; i < p.steps.size(); ++i) s += p.steps[i].L;
  return s;
}

int main(int argc, char** argv) {
  auto const tStart = std::chrono::steady_clock::now();
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
  double const dsBig = 8000.0 / 20.0;  // 400 m
  double const dsOne = 1.0;
  double const kConst = 1e-4;  // sanity check A [1/m]

  std::printf("=== Checkpoint 2e-e: extinction quadrature on the fixed clipped paths ===\n");
  std::printf("profile   : %s\n", argv[1]);
  std::printf("R         : %.17g m ; table height range [%.17g, %.17g] m\n", R, hr.first,
              hr.second);
  std::printf("paths     : endpoint-clipped 2e-d trails, ds = %.17g m and %.17g m\n", dsBig,
              dsOne);
  std::printf("integrand : alphaProd (Cherenkov.inl:96-113), unchanged\n");
  std::printf("panels    : composite midpoint, q = 1 m and q = 0.5 m, per accepted step\n");
  std::printf("experiment: DIAGNOSTIC ONLY, not production validation\n\n");

  double T400[3], T1[3], Q400_1[3], Q400_h[3], Q1_1[3], Q1_h[3];

  // ---------------------------------------------------------------- section 2
  std::printf("--- section 2: reproduce the endpoint-clipped baseline ---\n");
  double const expT400[3] = {4.5920502978309203, 1.0371636655991936, 1.0132736724165459};
  double const expT1[3] = {4.5625115139592811, 1.0359709986856456, 1.0132743713910122};
  for (int q = 0; q < 3; ++q) {
    double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
    double const ff = Fdisp(lm);
    double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
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
    refmarch::MarchResult const m400 = refmarch::march(R, hh, dsBig, p, Nof, der, alpha, true);
    refmarch::MarchResult const m1 = refmarch::march(R, hh, dsOne, p, Nof, der, alpha, true);
    T400[q] = m400.tau;
    T1[q] = m1.tau;
    std::printf("ray %d: T400=%.17g  T1=%.17g\n", q, T400[q], T1[q]);
    gate("2 T400 reproduces 2e-d", std::fabs(T400[q] - expT400[q]),
         1e-12 * std::max(1.0, std::fabs(expT400[q])));
    gate("2 T1 reproduces 2e-d", std::fabs(T1[q] - expT1[q]),
         1e-12 * std::max(1.0, std::fabs(expT1[q])));
  }
  if (g_gates_fail != 0) {
    std::printf("\nSTOP: clipped-baseline reproduction FAILED -- not continuing.\n");
    return 1;
  }
  std::printf("clipped baseline OK\n\n");


  // ------------------------------------------------------------ sections 3-5
  for (int q = 0; q < 3; ++q) {
    double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
    int const ilq = nearestLam(prof, lm);
    double const ff = Fdisp(lm);
    double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);
    refmarch::NFn Nof = [&](double rr) { return ri.N(rr) * ff; };
    refmarch::DerFn der = [&](double zz, double tt, double& dz, double& dth) {
      double const nTop = 1.0 + Nof(zz);
      double const slope = ri.dNdz(zz) * ff;
      double const s = std::sin(tt), c = std::cos(tt);
      dz = c;
      dth = -s * (1.0 / zz + (nTop > 1e-12 ? slope / nTop : 0.0));
    };
    refmarch::AlphaFn alpha = [&](double h) {
      return alphaProd(abs, lm, h, hr.first, hr.second);
    };
    refmarch::AlphaFn alphaK = [&](double) { return kConst; };

    // ---- section 3: the fixed clipped paths -------------------------------
    refpath::Path const pathB = refpath::buildPath(R, hh, dsBig, p, Nof, der, alpha, true);
    refpath::Path const pathO = refpath::buildPath(R, hh, dsOne, p, Nof, der, alpha, true);
    if (!pathB.ok || !pathO.ok) {
      std::printf("STOP: path recording failed for ray %d (no bracketed single crossing).\n", q);
      return 1;
    }
    double const sB = sumLengths(pathB);
    double const sO = sumLengths(pathO);
    std::printf("--- ray %d  lambda=%.17g nm  mu=%.17g  h=%.17g m ---\n", q, lm, mu, hh);
    std::printf("   400 m path: %zu accepted steps, sum L = %.17g m\n", pathB.steps.size(), sB);
    std::printf("     1 m path: %zu accepted steps, sum L = %.17g m\n", pathO.steps.size(), sO);
    gate("3 path(400) tau == T400 (bit-identical)", std::fabs(pathB.tau - T400[q]), 0.0);
    gate("3 path(1) tau == T1 (bit-identical)", std::fabs(pathO.tau - T1[q]), 0.0);
    refpath::EndpointResiduals const eB = refpath::endpointResiduals(pathB, der);
    refpath::EndpointResiduals const eO = refpath::endpointResiduals(pathO, der);
    std::printf("   endpoint reproduction: 400 m max|dr|=%.17g m  max|dtheta|=%.17g rad\n",
                eB.maxRadius, eB.maxAngle);
    std::printf("   endpoint reproduction:   1 m max|dr|=%.17g m  max|dtheta|=%.17g rad\n",
                eO.maxRadius, eO.maxAngle);
    gate("3 endpoints 400 m: |dr| <= 1e-8 m", eB.maxRadius, 1e-8);
    gate("3 endpoints 400 m: |dtheta| <= 1e-13 rad", eB.maxAngle, 1e-13);
    gate("3 endpoints 1 m: |dr| <= 1e-8 m", eO.maxRadius, 1e-8);
    gate("3 endpoints 1 m: |dtheta| <= 1e-13 rad", eO.maxAngle, 1e-13);


    // ---- section 4: composite-midpoint re-integration of extinction only ---
    refpath::Quad const Qb1 = refpath::compositeMidpoint(pathB, 1.0, der, alpha, R);
    refpath::Quad const Qbh = refpath::compositeMidpoint(pathB, 0.5, der, alpha, R);
    refpath::Quad const Qo1 = refpath::compositeMidpoint(pathO, 1.0, der, alpha, R);
    refpath::Quad const Qoh = refpath::compositeMidpoint(pathO, 0.5, der, alpha, R);
    Q400_1[q] = Qb1.value();
    Q400_h[q] = Qbh.value();
    Q1_1[q] = Qo1.value();
    Q1_h[q] = Qoh.value();
    std::printf("   Q400_1    = %.17g   panels=%ld  sum L=%.17g m\n", Q400_1[q], Qb1.panels,
                Qb1.totalLength);
    std::printf("   Q400_half = %.17g   panels=%ld  sum L=%.17g m\n", Q400_h[q], Qbh.panels,
                Qbh.totalLength);
    std::printf("   Q1_1      = %.17g   panels=%ld  sum L=%.17g m\n", Q1_1[q], Qo1.panels,
                Qo1.totalLength);
    std::printf("   Q1_half   = %.17g   panels=%ld  sum L=%.17g m\n", Q1_h[q], Qoh.panels,
                Qoh.totalLength);

    // ---- section 5A: constant-extinction bookkeeping ----------------------
    refpath::Quad const Kb1 = refpath::compositeMidpoint(pathB, 1.0, der, alphaK, R);
    refpath::Quad const Kbh = refpath::compositeMidpoint(pathB, 0.5, der, alphaK, R);
    refpath::Quad const Ko1 = refpath::compositeMidpoint(pathO, 1.0, der, alphaK, R);
    refpath::Quad const Koh = refpath::compositeMidpoint(pathO, 0.5, der, alphaK, R);
    double const expB = kConst * sB;
    double const expO = kConst * sO;
    std::printf("   A 400 m: expected k*sumL=%.17g  q=1: %.17g  q=0.5: %.17g\n", expB,
                Kb1.value(), Kbh.value());
    std::printf("   A   1 m: expected k*sumL=%.17g  q=1: %.17g  q=0.5: %.17g\n", expO,
                Ko1.value(), Koh.value());
    gate("5A 400 m q=1 constitutive", std::fabs(Kb1.value() - expB),
         1e-12 * std::max(1.0, std::fabs(expB)));
    gate("5A 400 m q=0.5 constitutive", std::fabs(Kbh.value() - expB),
         1e-12 * std::max(1.0, std::fabs(expB)));
    gate("5A 1 m q=1 constitutive", std::fabs(Ko1.value() - expO),
         1e-12 * std::max(1.0, std::fabs(expO)));
    gate("5A 1 m q=0.5 constitutive", std::fabs(Koh.value() - expO),
         1e-12 * std::max(1.0, std::fabs(expO)));

    // ---- section 5B: independent coefficient bounds -----------------------
    std::vector<double> slopes;
    for (std::size_t k = 0; k + 1 < prof.hgt.size(); ++k) {
      slopes.push_back((prof.tau[std::size_t(ilq)][k + 1] - prof.tau[std::size_t(ilq)][k]) /
                       (prof.hgt[k + 1] - prof.hgt[k]));
    }
    double rawMin = slopes[0], rawMax = slopes[0], minSpacing = 1e30;
    for (std::size_t k = 0; k < slopes.size(); ++k) {
      rawMin = std::min(rawMin, slopes[k]);
      rawMax = std::max(rawMax, slopes[k]);
    }
    for (std::size_t k = 0; k + 1 < prof.hgt.size(); ++k) {
      minSpacing = std::min(minSpacing, prof.hgt[k + 1] - prof.hgt[k]);
    }
    double const aMin = std::max(0.0, rawMin);
    double const aMax = std::max(0.0, rawMax);
    double worstRow = 0.0;
    for (std::size_t k = 0; k < prof.hgt.size(); k += 25) {
      double const direct = abs.getOpticalDepth(lm, prof.hgt[k]);
      worstRow = std::max(worstRow, std::fabs(direct - prof.tau[std::size_t(ilq)][k]));
    }
    std::printf("   B slopes: raw min=%.17g  raw max=%.17g -> alpha in [%.17g, %.17g] 1/m\n",
                rawMin, rawMax, aMin, aMax);
    std::printf("     support: %zu native segments, min spacing=%.17g m, window dh=50 m\n",
                slopes.size(), minSpacing);
    std::printf("     rule (tau(h+25)-tau(h-25))/50 clamped to [%.17g, %.17g] m\n", hr.first,
                hr.second);
    std::printf("     row check: max |tau_row(knot) - getOpticalDepth| = %.17g\n", worstRow);


    // per-step bound check on the four production-coefficient quads
    auto worstExcess = [&](refpath::Path const& pa, refpath::Quad const& qu) {
      double worst = 0.0;
      for (std::size_t i = 0; i < qu.perStep.size(); ++i) {
        double const L = pa.steps[i].L;
        double const c = qu.perStep[i];
        double const tol = 1e-12 * std::max(1.0, std::fabs(c));
        double const lo = aMin * L, hi = aMax * L;
        double const e = std::max(lo - c, c - hi);  // > 0 means outside the bounds
        worst = std::max(worst, std::max(0.0, e - tol));
      }
      return worst;
    };
    auto meanRange = [&](refpath::Path const& pa, refpath::Quad const& qu, double& lo,
                         double& hi) {
      lo = 1e30;
      hi = -1e30;
      for (std::size_t i = 0; i < qu.perStep.size(); ++i) {
        double const m = qu.perStep[i] / pa.steps[i].L;
        lo = std::min(lo, m);
        hi = std::max(hi, m);
      }
    };
    double mMinB = 0.0, mMaxB = 0.0, mMinO = 0.0, mMaxO = 0.0;
    meanRange(pathB, Qbh, mMinB, mMaxB);
    meanRange(pathO, Qoh, mMinO, mMaxO);
    std::printf("   B observed mean alpha per step: 400 m [%.17g, %.17g]  1 m [%.17g, %.17g]\n",
                mMinB, mMaxB, mMinO, mMaxO);
    double const wx = std::max(std::max(worstExcess(pathB, Qb1), worstExcess(pathB, Qbh)),
                               std::max(worstExcess(pathO, Qo1), worstExcess(pathO, Qoh)));
    gate("5B all step contributions within [alpha_min,alpha_max]*L", wx, 0.0);

    // ---- section 5C: refinement agreement ---------------------------------
    double const dQB = Qbh.value() - Qb1.value();
    double const dQO = Qoh.value() - Qo1.value();
    std::printf("   C 400 m: |Q_half-Q_1| = %.17g  relative = %.17g\n", std::fabs(dQB),
                std::fabs(dQB) / std::max(1.0, std::fabs(Qbh.value())));
    std::printf("   C   1 m: |Q_half-Q_1| = %.17g  relative = %.17g\n", std::fabs(dQO),
                std::fabs(dQO) / std::max(1.0, std::fabs(Qoh.value())));
    gate("5C refinement agreement 400 m", std::fabs(dQB),
         1e-6 * std::max(1.0, std::fabs(Qbh.value())));
    gate("5C refinement agreement 1 m", std::fabs(dQO),
         1e-6 * std::max(1.0, std::fabs(Qoh.value())));
    std::printf("\n");
  }


  // ---------------------------------------------------------------- section 6
  std::printf("--- section 6: decomposition (T400 - T1) = (C400 - C1) + DQ ---\n");
  double C400[3], C1[3], DQ[3], decRes[3];
  for (int q = 0; q < 3; ++q) {
    C400[q] = T400[q] - Q400_h[q];
    C1[q] = T1[q] - Q1_h[q];
    DQ[q] = Q400_h[q] - Q1_h[q];
    double const lhs = T400[q] - T1[q];
    double const rhs = (C400[q] - C1[q]) + DQ[q];
    decRes[q] = std::fabs(lhs - rhs);
    double const tol =
        1e-12 * std::max(std::max(1.0, std::fabs(T400[q])),
                         std::max(std::fabs(T1[q]),
                                  std::max(std::fabs(Q400_h[q]), std::fabs(Q1_h[q]))));
    std::printf("ray %d: T400-T1=%.17g  C400-C1=%.17g  DQ=%.17g  residual=%.17g\n", q, lhs,
                C400[q] - C1[q], DQ[q], decRes[q]);
    gate("6 decomposition identity", decRes[q], tol);
  }

  std::printf("\n=== table (17 significant digits) ===\n");
  std::printf("ray T400 T1 Q400_1 Q400_half Q1_1 Q1_half C400 C1 C400_minus_C1 DQ "
              "decomposition_residual\n");
  for (int q = 0; q < 3; ++q) {
    std::printf("%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", q,
                T400[q], T1[q], Q400_1[q], Q400_h[q], Q1_1[q], Q1_h[q], C400[q], C1[q],
                C400[q] - C1[q], DQ[q], decRes[q]);
  }

  auto const tEnd = std::chrono::steady_clock::now();
  double const secs = std::chrono::duration<double>(tEnd - tStart).count();
  std::printf("\ngates: %d evaluated, %d failed\n", g_gates, g_gates_fail);
  std::printf("%s\n", (g_gates_fail == 0) ? "ALL GATES PASS"
                                          : "STOP: a mandatory gate FAILED -- reported above");
  std::printf("runtime: %.3f s\n", secs);
  std::printf("=== end of 2e-e (diagnostic experiment, not production validation) ===\n");
  return (g_gates_fail == 0) ? 0 : 1;
}

