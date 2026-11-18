// ============================================================================
// stage2e_audit_b.cpp -- Checkpoint 2e-c: audit of the `db` bookkeeping.
//
//   ./stage2e_audit_b <profile .dat>
//
// One bounded diagnostic checkpoint. It does NOT diagnose physics and does NOT
// change any production code. It
//   1. recovers the fixed inputs of the existing 2e-b harness and names the
//      source of every reported quantity,
//   2. runs the three existing reference rays through (a) the diagnostic
//      integrations exactly as stage9b_abs_b.cpp computes them and (b) the
//      instrumented extracted tracer,
//   3. evaluates the mandatory sanity identities (tau accounting, db identity,
//      straight-line near-plane check, main/curved contribution fraction).
//
// The diagnostic operands are copied VERBATIM from stage2eb_cmp_b.cpp (which in
// turn takes them from stage9b_abs_b.cpp): same RefractiveIndexTable, same
// Fdisp, same alphaProd rule, same RK4 path-length march. Nothing is retuned.
//
// QUANTITY -> SOURCE (see the printed manifest):
//   tau_ref          tauIntS(ds = 1 m,  mode 1)  : stage9b_abs_b.cpp:413
//   tau_rule1m       tauIntS(ds = 1 m,  mode 0)  : stage9b_abs_b.cpp:414-415
//   tau_prod_curved  tauIntS(ds = 400 m, mode 0) : stage9b_abs_b.cpp:416-417
//                    (ds proxy: Cherenkov.inl:608-615 with H_scale ~ 8 km)
//   db_direct        tau_prod_curved - tau_rule1m: stage9b_abs_b.cpp:422
//   tau_tracer       refdiag::RefTracer::trace(), instrumented extraction of
//                    CherenkovProcess::rayTraceSphericalToPlane
//                    (Cherenkov.inl:358-719)
//   path figures     n_steps * ds  (NOMINAL step x count), NOT an accumulated
//                    segment length, except the tracer's geometric_path.
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include "ref_tracer_b.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0;

// --- diagnostic refractive model, identical to stage2eb_cmp_b.cpp:33-48 ------
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
// exact slope of the native piecewise-linear tau(h) (stage2eb_cmp_b.cpp:81-90)
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
// production rule, Cherenkov.inl:96-113 (dh = 50 m, centred), via the class
static double alphaProd(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                        double hHi) {
  double const hl = std::max(hLo, h - 25.0), hh = std::min(hHi, h + 25.0);
  double const adh = hh - hl;
  if (adh <= 0.0) return 0.0;
  double const t1 = abs.getOpticalDepth(lam, hh);
  double const t0 = abs.getOpticalDepth(lam, hl);
  return std::max(0.0, (t1 - t0) / adh);
}
// diagnostic march: RK4 in path length s from the emission down to r = R
// (stage2eb_cmp_b.cpp:104-135, verbatim)
static double tauIntS(AtmosphericAbsorption const& abs, Profile const& prof,
                      RefractiveIndexTable const& ri, double lam, int ilp, double p, double f,
                      double hTop, double ds, int mode, double hLo, double hHi, long* nSteps) {
  double z = R + hTop;
  double const nTop = 1.0 + ri.N(z) * f;
  double const sinT = std::min(1.0, p / (nTop * z));
  double th = std::acos(-std::sqrt(std::max(0.0, 1.0 - sinT * sinT)));
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
    ++n;
  }
  if (nSteps) *nSteps = n;
  return tau;
}

// ============================================================================
// audit plumbing
// ============================================================================
static int g_gates = 0, g_fails = 0;

static void gate(char const* name, double residual, double tolerance) {
  ++g_gates;
  bool const ok = (residual <= tolerance);
  if (!ok) ++g_fails;
  std::printf("   %-34s residual=%.17g  tol=%.17g  %s\n", name, residual, tolerance,
              ok ? "PASS" : "FAIL");
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <profile .dat>\n", argv[0]);
    return 2;
  }
  std::string const profArg = argv[1];
  char resolved[PATH_MAX];
  std::string const profResolved =
      (realpath(argv[1], resolved) != nullptr) ? std::string(resolved) : profArg;

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

  // The three existing reference rays, unchanged (stage2eb_cmp_b.cpp:160).
  struct C {
    double lam, mu, h;
  };
  C const cs[3] = {{300.0, 0.6, 46000.0}, {480.0, 0.1, 1000.0}, {480.0, 0.1, 975.0}};

  double const dsRef = 1.0;        // tau_ref / tau_rule1m step        [m]
  double const dsTr = 8000.0 / 20.0;  // tau_prod_curved step (400 m)  [m]
  double const tauStart = 0.0;     // tracer entry accumulated_tau
  double const survivalThresholdIn = 0.0;  // disables early termination

  std::printf("=== Checkpoint 2e-c: audit of the db bookkeeping ===\n");
  std::printf("profile argument : %s\n", profArg.c_str());
  std::printf("profile resolved : %s\n", profResolved.c_str());
  std::printf("profile knots    : %zu wavelengths, %zu heights, h=[%.17g, %.17g] m\n",
              prof.lam.size(), prof.hgt.size(), prof.hgt.front(), prof.hgt.back());
  std::printf("table height range (AtmosphericAbsorption::getHeightRange) = [%.17g, %.17g] m\n",
              hr.first, hr.second);
  std::printf("earth radius R   : %.17g m  (both implementations)\n", R);
  std::printf("diagnostic ds    : ref=%.17g m, rule1m=%.17g m, prod_curved=%.17g m\n",
              dsRef, dsRef, dsTr);
  std::printf("tau_start        : %.17g  (refdiag::RefTracer::trace tau_start argument)\n",
              tauStart);
  std::printf("survivalThreshold: %.17g  (trace survivalThreshold_in argument)\n",
              survivalThresholdIn);
  std::printf("tracer plane     : point (0, 0, %.17g), normal (0, 0, 1)\n", R);
  std::printf("tracer callables : f_phase == f_group == Fdisp(lambda) of the same\n");
  std::printf("                   RefractiveIndexTable used by the diagnostic;\n");
  std::printf("                   extinction = alphaProd(Cherenkov.inl:96-113).\n");
  std::printf("mu convention    : mu = cos(angle between the emission direction and the\n");
  std::printf("                   OUTWARD radial unit vector); dz = -mu, so mu > 0 is\n");
  std::printf("                   downward. Horizontal part = sqrt(1-mu^2) along +x.\n");
  std::printf("path figures     : n_steps * ds unless stated otherwise (NOMINAL, not an\n");
  std::printf("                   accumulated length); tracer path is geometric_path.\n");
  std::printf("units            : length m, tau dimensionless, extinction 1/m, lambda nm.\n\n");


  for (int q = 0; q < 3; ++q) {
    double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
    int const ilq = nearestLam(prof, lm);
    double const ff = Fdisp(lm);
    double const p = std::sqrt(std::max(0.0, 1.0 - mu * mu)) * bCrit(ri, ff);

    // What the diagnostic's initial angle actually is: p = sqrt(1-mu^2)*bCrit is a
    // GROUND-level impact parameter, so the march starts at
    //   sin(theta_top) = p / (n(R+h) * (R+h))
    // and the implied local |cos| equals mu only if n(R)R == n(R+h)(R+h).
    double const zTop = R + hh;
    double const nTopDiag = 1.0 + ri.N(zTop) * ff;
    double const sinTopDiag = std::min(1.0, p / (nTopDiag * zTop));
    double const muDiagImplied = std::sqrt(std::max(0.0, 1.0 - sinTopDiag * sinTopDiag));

    // --- the diagnostic columns, exactly as stage9b_abs_b.cpp computes them ---
    long n1 = 0, n2 = 0, n3 = 0;
    double const tauRef =
        tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsRef, 1, hr.first, hr.second, &n1);
    double const tauRule1m =
        tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsRef, 0, hr.first, hr.second, &n2);
    double const tauProdCurved =
        tauIntS(abs, prof, ri, lm, ilq, p, ff, hh, dsTr, 0, hr.first, hr.second, &n3);
    // stage9b_abs_b.cpp:422, reproduced as its own statement
    double const dbDirect = tauProdCurved - tauRule1m;
    // the identity this audit has to verify
    double const dbRecomputed = tauProdCurved - tauRule1m;

    // --- the extracted tracer, same callables and same fixed inputs ----------
    refdiag::RefTracer tr;
    tr.earthRadiusM_ = R;
    tr.enableAbsorption_ = true;
    tr.refractiveIndexFunc_ = [&](double h_m, double lam_nm) {
      double const f = Fdisp(lam_nm);
      double const n = 1.0 + ri.N(R + h_m) * f;
      return std::make_pair(n, n);  // equal-f convention, as the diagnostic uses
    };
    tr.extinctionCoeffFunc_ = [&](double h_m, double lam_nm) {
      return alphaProd(abs, lam_nm, h_m, hr.first, hr.second);
    };

    double const sinT = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    double const dirNorm = std::sqrt(sinT * sinT + mu * mu);
    double const x0z = R + hh;  // emission point z (x = y = 0)
    // n . u for the flat plane normal (0,0,1)
    double const denNU = -mu;
    refdiag::RefTracer::Result const r =
        tr.trace(0.0, 0.0, x0z, sinT, 0.0, -mu, 0.0, 0.0, R, 0.0, 0.0, 1.0, lm,
                 survivalThresholdIn, tauStart);

    double const tauTr = r.accumulated_tau;
    double const tMinusDiag = tauTr - tauProdCurved;
    double const eTracerDiag =
        std::fabs(tMinusDiag) /
        std::max(std::max(std::fabs(tauTr), std::fabs(tauProdCurved)), 1e-30);

    // callable-identity probes (must be bit-identical, same functions/table)
    double const probe = 1000.0, probeF = Fdisp(lm);
    auto const riProbe = tr.refractiveIndexFunc_(probe, lm);
    double const diagProbe = nPhase(ri, R + probe, probeF);
    double const kextProbeTr = tr.extinctionCoeffFunc_(probe, lm);
    double const kextProbeDiag = alphaProd(abs, lm, probe, hr.first, hr.second);
    bool const fPhaseEqFGroup = (riProbe.first == riProbe.second);
    bool const riSame = (riProbe.first == diagProbe);
    bool const kextSame = (kextProbeTr == kextProbeDiag);

    std::printf("--- ray %d: lambda=%.17g nm  mu=%.17g  h=%.17g m ---\n", q, lm, mu, hh);
    std::printf("   nearest native knot: lambda=%.17g nm  |dlambda|=%.17g nm\n",
                prof.lam[std::size_t(ilq)], std::fabs(prof.lam[std::size_t(ilq)] - lm));
    std::printf("   Fdisp(lambda)=%.17g\n", ff);
    std::printf("   diagnostic angle: sin(theta_top)=%.17g  implied |cos|=%.17g  "
                "(ray mu=%.17g)\n",
                sinTopDiag, muDiagImplied, mu);
    std::printf("   emission=(0, 0, %.17g)  direction=(%.17g, 0, %.17g)  |u|=%.17g  n.u=%.17g\n",
                x0z, sinT, -mu, dirNorm, denNU);
    std::printf("   callables: tracerRI(h=%.17g)=(%.17g, %.17g)  f_phase==f_group:%s  "
                "diag nPhase=%.17g  RI identical:%s\n",
                probe, riProbe.first, riProbe.second, fPhaseEqFGroup ? "yes" : "NO",
                diagProbe, riSame ? "yes" : "NO");
    std::printf("              kext tracer=%.17g  alphaProd=%.17g  identical:%s\n",
                kextProbeTr, kextProbeDiag, kextSame ? "yes" : "NO");
    std::printf("   tau_start=%.17g\n", tauStart);
    std::printf("   tau_ref=%.17g\n", tauRef);
    std::printf("   tau_rule1m=%.17g\n", tauRule1m);
    std::printf("   tau_prod_curved=%.17g\n", tauProdCurved);
    std::printf("   tau_tracer=%.17g\n", tauTr);
    std::printf("   db_direct=%.17g  db_recomputed=%.17g\n", dbDirect, dbRecomputed);
    std::printf("   tracer_minus_diag=%.17g  E_tracer_diag=%.17g\n", tMinusDiag, eTracerDiag);
    std::printf("   termination=%s  legacy_branch=%d  hit=%d  steps=%d\n",
                refdiag::terminationName(r.termination), r.branch, int(r.hit), r.steps);
    std::printf("   main_advances=%d  turning_events=%d\n", r.main_advances, r.turning);
    std::printf("   tau_main=%.17g  tau_near=%.17g\n", r.tau_main, r.tau_near);
    std::printf("   geometric_path_tracer=%.17g m  optical_path_tracer=%.17g m\n",
                r.geometric_path, r.optical_path);
    std::printf("   impact=(%.17g, %.17g, %.17g)\n", r.x_impact, r.y_impact, r.z_impact);
    if (r.near_applicable) {
      std::printf("   near_segment_length=%.17g m  near_tau_substeps=%d\n", r.near_seg,
                  r.near_sub);
    } else {
      std::printf("   near_segment_length=not_applicable  near_tau_substeps=not_applicable\n");
    }
    std::printf("   path_ref_nominal=%.17g m  (%ld steps x %.17g m)\n", double(n1) * dsRef,
                n1, dsRef);
    std::printf("   path_rule1m_nominal=%.17g m  (%ld steps x %.17g m)\n", double(n2) * dsRef,
                n2, dsRef);
    std::printf("   path_prod_curved_nominal=%.17g m  (%ld steps x %.17g m)\n",
                double(n3) * dsTr, n3, dsTr);
    std::printf("   path_straight_line h/mu=%.17g m\n", hh / mu);


    // ---------------- mandatory sanity checks -------------------------------
    std::printf("   checks:\n");
    // A. tau accounting: tau_tracer - tau_start - tau_main - tau_near
    {
      double const resid = std::fabs(tauTr - tauStart - r.tau_main - r.tau_near);
      double const tol = 1e-12 * std::max(1.0, std::fabs(tauTr));
      gate("A tau accounting", resid, tol);
    }
    // B. db accounting: db_direct vs tau_prod_curved - tau_rule1m
    {
      double const resid = std::fabs(dbDirect - (tauProdCurved - tauRule1m));
      double const scale =
          std::max(std::max(1.0, std::fabs(dbDirect)),
                   std::max(std::fabs(tauProdCurved), std::fabs(tauRule1m)));
      gate("B db accounting", resid, 1e-12 * scale);
    }
    // C. straight-line back-of-envelope, only for a ray that reaches the
    //    near-plane branch before ANY completed main advance
    bool const nearFirst =
        (r.main_advances == 0) &&
        (r.termination == refdiag::Termination::NearPlane ||
         r.termination == refdiag::Termination::AbsorptionNearPlane);
    if (nearFirst) {
      double const nDotDelta = (x0z - R) * 1.0;  // n . (x0 - x_plane), n = (0,0,1)
      double const s = -nDotDelta / denNU;
      double const normDev = std::fabs(dirNorm - 1.0);
      std::printf("      C preconditions: |u|-1=%.17g (need <=1e-12)  n.u=%.17g (need !=0)\n",
                  normDev, denNU);
      std::printf("      C s=%.17g m  h/mu=%.17g m  |s-h/mu|=%.17g\n", s, hh / mu,
                  std::fabs(s - hh / mu));
      std::printf("      C s_near(recorded)=%.17g m  |s-s_near|=%.17g\n", r.near_seg,
                  std::fabs(s - r.near_seg));
      gate("C straight-line s vs s_near", std::fabs(s - r.near_seg),
           1e-9 * std::max(1.0, std::fabs(s)));
    } else {
      std::printf("      C straight-line check: NOT APPLICABLE "
                  "(main_advances=%d, termination=%s)\n",
                  r.main_advances, refdiag::terminationName(r.termination));
    }
    // D. contribution fraction (descriptive, not a gate)
    if (tauTr != 0.0) {
      std::printf("      D f_main = tau_main/tau_tracer = %.17g\n", r.tau_main / tauTr);
    } else {
      std::printf("      D f_main = not_applicable (tau_tracer == 0)\n");
    }
    std::printf("\n");
    std::fflush(stdout);
  }

  std::printf("=== summary: %d gate(s) evaluated, %d failed ===\n", g_gates, g_fails);
  std::printf("%s\n", (g_fails == 0)
                          ? "ALL APPLICABLE GATES PASS"
                          : "STOP: a gate FAILED -- report the failed identity only");
  std::printf("=== end of 2e-c audit ===\n");
  return (g_fails == 0) ? 0 : 1;
}

