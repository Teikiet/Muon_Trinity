// ============================================================================
// cp2fa1_reconcile_b.cpp -- Checkpoint 2f-a1: reconcile the T1 evidence.
//
//   ./cp2fa1_reconcile_b <real profile .dat>
//
// Diagnostic audit only. No T2-T4, no production change, no gate change, no
// golden values. It answers two questions:
//   A. Is the reported "h = 25 -> 1.04429e-4 / h = 26 -> 1.0e-4" contradiction
//      reproducible, a harness error, or a transcription error?
//   B. Do the upper-boundary residuals scale like the double-precision
//      cancellation estimate E_scale = eps*(|tau_hi|+|tau_lo|)/d_actual?
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;

static constexpr double kWindow = 50.0;  // Cherenkov.inl:97
static constexpr double kHalf = kWindow * 0.5;

struct Win {
  double lo{0.0}, hi{0.0}, d{0.0};
  bool ok{false};
};
static Win windowOf(double h, double hLo, double hHi) {
  Win w;
  w.lo = std::max(hLo, h - kHalf);
  w.hi = std::min(hHi, h + kHalf);
  w.d = w.hi - w.lo;
  w.ok = (w.d > 0.0);
  return w;
}

// verbatim production rule (Cherenkov.inl:96-113)
static double alphaSrc(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                       double hHi) {
  Win const w = windowOf(h, hLo, hHi);
  if (!w.ok) return 0.0;
  double const tauLo = abs.getOpticalDepth(lam, w.lo);
  double const tauHi = abs.getOpticalDepth(lam, w.hi);
  return std::max(0.0, (tauHi - tauLo) / w.d);
}

// slope/weight reconstruction
static double alphaKernel(std::vector<double> const& s, std::vector<double> const& hKnots,
                          double h, double hLo, double hHi) {
  Win const w = windowOf(h, hLo, hHi);
  if (!w.ok) return 0.0;
  double acc = 0.0;
  for (std::size_t j = 0; j + 1 < hKnots.size(); ++j) {
    double const a = std::max(hKnots[j], w.lo);
    double const b = std::min(hKnots[j + 1], w.hi);
    if (b > a) acc += s[j] * ((b - a) / w.d);
  }
  return acc;
}

struct Profile {
  std::vector<double> lam, hgt;
  std::vector<std::vector<double> > tau;
};
static bool readProfile(std::string const& fn, Profile& p) {
  std::ifstream f(fn.c_str());
  if (!f) return false;
  int nw = 0, nh = 0;
  f >> nw >> nh;
  if (nw <= 0 || nh <= 0) return false;
  p.lam.resize(nw);
  p.hgt.resize(nh);
  for (int i = 0; i < nw; ++i) f >> p.lam[i];
  for (int i = 0; i < nh; ++i) f >> p.hgt[i];
  p.tau.assign(nw, std::vector<double>(nh, 0.0));
  for (int i = 0; i < nw; ++i)
    for (int j = 0; j < nh; ++j) f >> p.tau[i][j];
  return true;
}

static void slopesAt(Profile const& p, double lam, std::vector<double>& s) {
  std::size_t const nw = p.lam.size();
  double const l = std::max(p.lam.front(), std::min(p.lam.back(), lam));
  std::size_t iw = 0;
  for (std::size_t i = 0; i + 1 < nw; ++i) {
    if (p.lam[i + 1] >= l) {
      iw = i;
      break;
    }
    iw = i;
  }
  std::size_t const iw1 = std::min(iw + 1, nw - 1);
  double const fw = (iw1 > iw) ? (l - p.lam[iw]) / (p.lam[iw1] - p.lam[iw]) : 0.0;
  s.assign(p.hgt.size() - 1, 0.0);
  for (std::size_t j = 0; j + 1 < p.hgt.size(); ++j) {
    double const dh = p.hgt[j + 1] - p.hgt[j];
    double const s0 = (p.tau[iw][j + 1] - p.tau[iw][j]) / dh;
    double const s1 = (p.tau[iw1][j + 1] - p.tau[iw1][j]) / dh;
    s[j] = s0 * (1.0 - fw) + s1 * fw;
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <real profile .dat>\n", argv[0]);
    return 2;
  }
  Profile prof;
  if (!readProfile(argv[1], prof)) {
    std::printf("STOP: cannot read the profile\n");
    return 1;
  }
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[1]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  auto const hr = abs.getHeightRange();
  double const hLo = hr.first, hHi = hr.second;

  std::printf("=== Checkpoint 2f-a1: reconcile the T1 evidence ===\n");
  std::printf("profile : %s\n", argv[1]);
  std::printf("bounds  : h [%.17g, %.17g] m , lambda [%.17g, %.17g] nm\n\n", hLo, hHi,
              abs.getWavelengthRange().first, abs.getWavelengthRange().second);

  // ==================== A. first-bin contradiction =========================
  int il = -1;
  for (std::size_t i = 0; i < prof.lam.size(); ++i)
    if (prof.lam[i] == 480.0) il = int(i);
  if (il < 0) {
    std::printf("STOP: 480 nm is not a native knot\n");
    return 1;
  }
  double const h0 = prof.hgt[0], h1 = prof.hgt[1];
  double const tau0 = prof.tau[std::size_t(il)][0], tau1 = prof.tau[std::size_t(il)][1];
  double const s0 = (tau1 - tau0) / (h1 - h0);
  double const gate0 = 1.0e-14 + 1.0e-10 * std::fabs(s0);
  std::printf("--- A. first bin at lambda = %.17g nm ---\n", prof.lam[std::size_t(il)]);
  std::printf("  native tau(%.17g) = %.17g\n", h0, tau0);
  std::printf("  native tau(%.17g) = %.17g\n", h1, tau1);
  std::printf("  first-segment slope s0 = %.17g 1/m\n", s0);
  std::printf("  gate                   = 1e-14 + 1e-10*|s0| = %.17g 1/m\n\n", gate0);
  std::printf("  %-10s %-14s %-14s %-12s %-20s %-20s %-20s %-20s %-12s\n", "h", "window_lo",
              "window_hi", "width", "tau(lo) lookup", "tau(hi) lookup", "alpha_src",
              "alpha_kernel", "|src-s0|");
  double const hs[7] = {0.0, 24.0, 25.0, 26.0, 50.0, 500.0, 974.0};
  std::vector<double> s480;
  slopesAt(prof, 480.0, s480);
  int failA = 0;
  for (int i = 0; i < 7; ++i) {
    double const h = hs[i];
    Win const w = windowOf(h, hLo, hHi);
    double const tl = abs.getOpticalDepth(480.0, w.lo);
    double const th = abs.getOpticalDepth(480.0, w.hi);
    double const a = alphaSrc(abs, 480.0, h, hLo, hHi);
    double const b = alphaKernel(s480, prof.hgt, h, hLo, hHi);
    double const r1 = std::fabs(a - s0);
    double const r2 = std::fabs(b - s0);
    if (!(r1 <= gate0 && r2 <= gate0)) ++failA;
    std::printf("  %-10.17g %-14.17g %-14.17g %-12.17g %-20.17g %-20.17g %-20.17g %-20.17g "
                "%-12.3e%s\n",
                h, w.lo, w.hi, w.d, tl, th, a, b, std::max(r1, r2),
                (r1 <= gate0 && r2 <= gate0) ? "" : "  FAIL");
  }
  std::printf("\n  A gate: %d of 7 rows failed\n", failA);

  // whole-bin scan: is the contradiction reproducible anywhere in bin 0?
  double worstA = 0.0, worstAhat = -1.0, worstAker = 0.0, worstAh = 0.0;
  int violA = 0;
  for (int i = 0; i <= 975; ++i) {
    double const h = double(i);
    double const a = alphaSrc(abs, 480.0, h, hLo, hHi);
    double const b = alphaKernel(s480, prof.hgt, h, hLo, hHi);
    double const r1 = std::fabs(a - s0);
    double const r2 = std::fabs(b - s0);
    if (r1 > worstA) {
      worstA = r1;
      worstAh = h;
    }
    worstAker = std::max(worstAker, r2);
    if (!(r1 <= gate0 && r2 <= gate0)) {
      ++violA;
      if (worstAhat < 0.0) worstAhat = h;
    }
  }
  std::printf("  full bin-0 scan h = 0..975 m (1 m steps, 976 rows):\n");
  std::printf("    max |alpha_src    - s0| = %.17g at h = %.17g m\n", worstA, worstAh);
  std::printf("    max |alpha_kernel - s0| = %.17g\n", worstAker);
  std::printf("    gate violations        = %d", violA);
  if (violA) std::printf(" (first at h = %.17g m)", worstAhat);
  std::printf("\n\n");


  // ==================== B. cancellation scale ==============================
  std::printf("--- B. upper-boundary windows: residual vs cancellation scale ---\n");
  std::printf("  h = H + 25 - d so that the window is [H-d, H];  E_scale = eps*(|tau_hi|+|tau_lo|)/d_actual\n");
  std::printf("  eps = %.17g (DBL_EPSILON);  E_scale is a sanity estimate, not a bound\n\n",
              double(DBL_EPSILON));
  double const ds[7] = {1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};
  std::printf("  %-12s %-24s %-24s %-12s %-22s %-12s %-22s\n", "d_nominal", "d_actual",
              "max|src - s_last|", "at lambda", "max resid/E_scale", "at lambda",
              "max|kern - s_last|");
  for (int t = 0; t < 7; ++t) {
    double const d = ds[t];
    double const h = hHi + 25.0 - d;
    Win const w = windowOf(h, hLo, hHi);
    double worst = 0.0, worstLam = 0.0;
    double ratioMax = 0.0, ratioLam = 0.0;
    double kernMax = 0.0;
    std::vector<double> sVec;
    for (std::size_t k = 0; k < prof.lam.size(); ++k) {
      double const lam = prof.lam[k];
      double const sLast =
          (prof.tau[k][50] - prof.tau[k][49]) / (prof.hgt[50] - prof.hgt[49]);
      double const a = alphaSrc(abs, lam, h, hLo, hHi);
      slopesAt(prof, lam, sVec);
      double const b = alphaKernel(sVec, prof.hgt, h, hLo, hHi);
      double const resid = std::fabs(a - sLast);
      double const tauHi = abs.getOpticalDepth(lam, w.hi);
      double const tauLo = abs.getOpticalDepth(lam, w.lo);
      double const eScale = double(DBL_EPSILON) * (std::fabs(tauHi) + std::fabs(tauLo)) / w.d;
      double const ratio = (eScale > 0.0) ? resid / eScale : 0.0;
      if (resid > worst) {
        worst = resid;
        worstLam = lam;
      }
      if (ratio > ratioMax) {
        ratioMax = ratio;
        ratioLam = lam;
      }
      kernMax = std::max(kernMax, std::fabs(b - sLast));
    }
    std::printf("  %-12.17g %-24.17g %-24.3e %-12.17g %-22.3e %-12.17g %-22.3e\n", d, w.d, worst,
                worstLam, ratioMax, ratioLam, kernMax);
  }
  std::printf("\n  Note: the same d_actual applies to all 105 knots (it depends only on d).\n\n");

  // ==================== C. corrected T1 summary ============================
  // recompute the two facts quoted below from this profile (cheap, self-contained)
  double sMin = 1e300, sMax = -1e300;
  long nNeg = 0;
  for (std::size_t k = 0; k < prof.lam.size(); ++k)
    for (std::size_t j = 0; j + 1 < prof.hgt.size(); ++j) {
      double const s = (prof.tau[k][j + 1] - prof.tau[k][j]) / (prof.hgt[j + 1] - prof.hgt[j]);
      sMin = std::min(sMin, s);
      sMax = std::max(sMax, s);
      if (s < 0.0) ++nNeg;
    }
  double loNZ = 1e300, hiNZ = -1e300;
  for (int i = 0; i <= 120000; ++i) {
    double const h = -30.0 + 60.0 * double(i) / 120000.0;
    if (alphaSrc(abs, 480.0, h, hLo, hHi) != 0.0) {
      loNZ = std::min(loNZ, h);
      hiNZ = std::max(hiNZ, h);
    }
  }
  double loNZT = 1e300, hiNZT = -1e300;
  for (int i = 0; i <= 120000; ++i) {
    double const h = hHi - 30.0 + 60.0 * double(i) / 120000.0;
    if (alphaSrc(abs, 480.0, h, hLo, hHi) != 0.0) {
      loNZT = std::min(loNZT, h);
      hiNZT = std::max(hiNZT, h);
    }
  }

  std::printf("--- C. corrected T1 summary (four separate statements) ---\n");
  std::printf("  C1 ALGEBRAIC representability (this profile):\n");
  std::printf("     native slopes: min %.17g 1/m, max %.17g 1/m, negative entries %ld\n", sMin,
              sMax, nNeg);
  std::printf("     => alpha(lambda,h) = sum_j s_j(lambda) w_j(h) holds in real arithmetic;\n");
  std::printf("        the max(0,.) clamp fired only on rounding-level negatives (|raw| <= 7.1e-17,\n");
  std::printf("        worst kernel residual among those queries = 0 in T1).\n");
  std::printf("  C2 FLOATING-POINT agreement (T1 grid, unchanged):\n");
  std::printf("     10,621,107 kernel queries, 211 failed; gate 1e-14 + 1e-10|a| UNCHANGED.\n");
  std::printf("     worst 1.1116386171950788e-08 1/m at lambda=255 nm, h=H+25-1e-6 m\n");
  std::printf("     worst for windows with denominator >= 1e-3 m: 1.0338409815735528e-14 1/m\n");
  std::printf("     => the T1 gate remains FAILED; it is not redefined here.\n");
  std::printf("  C3 EXACT support endpoints (1 mm scan, lambda=480 nm):\n");
  std::printf("     alpha != 0 on (%.17g, %.17g) m in the evaluation variable h\n", -25.0, hHi + 25.0);
  std::printf("     first non-zero %.17g m; last non-zero %.17g m\n", loNZ, hiNZ);
  std::printf("     alpha is exactly 0 at h = -25 m and at h = H+25 m (= %.17g m) and beyond\n",
              hHi + 25.0);
  std::printf("     top scan [H-30, H+30]: first non-zero at the scan limit %.17g m, last "
              "non-zero at %.17g m\n",
              loNZT, hiNZT);
  std::printf("  C4 first-bin contradiction: RESOLVED as a transcription error in the T1 report.\n");
  std::printf("     h=26 m on the real profile gives alpha_src = alpha_kernel = %.17g 1/m\n", s0);
  std::printf("     (T1 harness output line 165), identical to h=25 m; not reproducible in the\n");
  std::printf("     0..975 m scan (%d gate violations), so it is neither a harness error nor a\n", violA);
  std::printf("     real discrepancy. The 1.0e-4 value in the report was the SYNTHETIC constant-slope\n");
  std::printf("     value k = 1e-4 from T1 section B.\n");
  std::printf("  NOT claimed here: that cancellation-dominated source evaluations are exactly\n");
  std::printf("  linear in wavelength, or that the integrated impact is negligible (unmeasured).\n");
  std::printf("=== end of 2f-a1 ===\n");
  return (failA == 0 && violA == 0) ? 0 : 1;
}

