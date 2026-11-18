// ============================================================================
// cp2fa_kernel_b.cpp -- Checkpoint 2f-a (T1): source-faithful extinction kernel.
//
//   ./cp2fa_kernel_b <real profile .dat>
//
// Diagnostic only. Establishes whether production extinction
// (Cherenkov.inl:96-113) is representable as
//     alpha(lambda,h) = sum_j s_j(lambda) * w_j(h),
// with s_j the native profile slopes and w_j wavelength-independent altitude
// weights of the +-25 m window. No ray integration, no reachability code, no
// compression, no production edits.
//
// Two INDEPENDENT implementations are compared:
//   alphaSrc()    : verbatim production rule, calls AtmosphericAbsorption::
//                   getOpticalDepth (the production lookup), including the
//                   max(0, ...) clamp and the actual_dh <= 0 guard.
//   alphaKernel() : the slope/weight reconstruction, no getOpticalDepth call.
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;

// --- production constants, Cherenkov.inl:97 ---------------------------------
static constexpr double kWindow = 50.0;   // double const dh = 50.0;
static constexpr double kHalf = kWindow * 0.5;

// --- the production window, Cherenkov.inl:99-103 ----------------------------
struct Win {
  double lo{0.0}, hi{0.0}, d{0.0};
  bool ok{false};
};
static Win windowOf(double h, double hLo, double hHi) {
  Win w;
  w.lo = std::max(hLo, h - kHalf);
  w.hi = std::min(hHi, h + kHalf);
  w.d = w.hi - w.lo;
  w.ok = (w.d > 0.0);  // "if (actual_dh <= 0.0) return 0.0;"
  return w;
}

// --- source-faithful callable: Cherenkov.inl:96-113, statement for statement --
static double alphaSrc(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                       double hHi, double* raw = nullptr) {
  Win const w = windowOf(h, hLo, hHi);
  if (!w.ok) {
    if (raw) *raw = 0.0;
    return 0.0;
  }
  double const tauLo = abs.getOpticalDepth(lam, w.lo);
  double const tauHi = abs.getOpticalDepth(lam, w.hi);
  double const a = (tauHi - tauLo) / w.d;
  if (raw) *raw = a;
  return std::max(0.0, a);
}

// --- the reconstruction under test: alpha = sum_j s_j * (overlap_j / Delta) ---
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

// --- native profile, as stored -------------------------------------------------
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

// --- slopes of the interpolated tau at a query wavelength ----------------------
// tau(lambda, .) is bilinear (AtmosphericAbsorption.hpp:140-158), hence linear in
// lambda between knots and linear in h between knots, so the segment slope is
//     s_j(lambda) = (1-fw) s_j(lambda_i) + fw s_j(lambda_i+1).
// The lambda clamp copies getOpticalDepth's own clamp (line 112).
static void slopesAt(Profile const& p, double lam, std::vector<double>& s, double& lamUsed) {
  std::size_t const nw = p.lam.size();
  double const l = std::max(p.lam.front(), std::min(p.lam.back(), lam));
  lamUsed = l;
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

// --- synthetic cumulative profile writer (diagnostic input only) ---------------
static bool writeSynthetic(std::string const& fn, bool zeroSlope, double k, double constVal) {
  std::vector<double> lam{300.0, 480.0, 700.0};
  std::vector<double> hgt;
  for (int i = 0; i <= 50; ++i) hgt.push_back(1000.0 * double(i));
  std::ofstream f(fn.c_str());
  if (!f) return false;
  f.precision(17);
  f << lam.size() << " " << hgt.size() << "\n";
  for (std::size_t i = 0; i < lam.size(); ++i) f << lam[i] << " ";
  f << "\n";
  for (std::size_t j = 0; j < hgt.size(); ++j) f << hgt[j] << " ";
  f << "\n";
  for (std::size_t i = 0; i < lam.size(); ++i) {
    for (std::size_t j = 0; j < hgt.size(); ++j)
      f << (zeroSlope ? constVal : k * hgt[j]) << " ";
    f << "\n";
  }
  return true;
}

// --- the SAME interpolation and window rule evaluated in long double ----------
// Purpose: separate a kernel-representation error from a rounding error inside the
// production expression. The bracket search is copied from getOpticalDepth; only
// the arithmetic is widened.
static long double alphaLD(Profile const& p, double lam, double h, double hLo, double hHi) {
  Win const w = windowOf(h, hLo, hHi);
  if (!w.ok) return 0.0L;
  auto tauLD = [&](long double hl) -> long double {
    long double const l = std::max((long double)p.lam.front(),
                                  std::min((long double)p.lam.back(), (long double)lam));
    long double const hh = std::max((long double)p.hgt.front(),
                                   std::min((long double)p.hgt.back(), hl));
    std::size_t iw = 0;
    for (std::size_t i = 0; i + 1 < p.lam.size(); ++i) {
      if ((long double)p.lam[i + 1] >= l) {
        iw = i;
        break;
      }
      iw = i;
    }
    std::size_t ih = 0;
    for (std::size_t i = 0; i + 1 < p.hgt.size(); ++i) {
      if ((long double)p.hgt[i + 1] >= hh) {
        ih = i;
        break;
      }
      ih = i;
    }
    std::size_t const iw1 = std::min(iw + 1, p.lam.size() - 1);
    std::size_t const ih1 = std::min(ih + 1, p.hgt.size() - 1);
    long double const fw = (iw1 > iw) ? (l - (long double)p.lam[iw]) /
                                            ((long double)p.lam[iw1] - (long double)p.lam[iw])
                                      : 0.0L;
    long double const fh = (ih1 > ih) ? (hh - (long double)p.hgt[ih]) /
                                            ((long double)p.hgt[ih1] - (long double)p.hgt[ih])
                                      : 0.0L;
    long double const t0 = (long double)p.tau[iw][ih] * (1.0L - fh) +
                           (long double)p.tau[iw][ih1] * fh;
    long double const t1 = (long double)p.tau[iw1][ih] * (1.0L - fh) +
                           (long double)p.tau[iw1][ih1] * fh;
    return t0 * (1.0L - fw) + t1 * fw;
  };
  long double const d = (long double)w.hi - (long double)w.lo;
  return (tauLD((long double)w.hi) - tauLD((long double)w.lo)) / d;
}

static int g_points = 0, g_fail = 0;

int main(int argc, char** argv) {
  auto const tStart = std::chrono::steady_clock::now();
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
  auto const wr = abs.getWavelengthRange();
  double const hLo = hr.first, hHi = hr.second;

  std::printf("=== Checkpoint 2f-a: source-faithful extinction kernel ===\n");
  std::printf("profile : %s\n\n", argv[1]);

  // ======================= A. profile facts ================================
  std::printf("--- A. profile facts (source: AtmosphericAbsorption) ---\n");
  std::printf("getHeightRange      = [%.17g, %.17g] m\n", hr.first, hr.second);
  std::printf("getWavelengthRange  = [%.17g, %.17g] nm\n", wr.first, wr.second);
  std::printf("knots               = %zu wavelengths, %zu heights\n", prof.lam.size(),
              prof.hgt.size());
  double dmin = 1e30, dmax = -1e30;
  for (std::size_t i = 0; i + 1 < prof.hgt.size(); ++i) {
    double const d = prof.hgt[i + 1] - prof.hgt[i];
    dmin = std::min(dmin, d);
    dmax = std::max(dmax, d);
  }
  std::printf("height spacing      = min %.17g m, max %.17g m\n", dmin, dmax);
  dmin = 1e30;
  dmax = -1e30;
  for (std::size_t i = 0; i + 1 < prof.lam.size(); ++i) {
    double const d = prof.lam[i + 1] - prof.lam[i];
    dmin = std::min(dmin, d);
    dmax = std::max(dmax, d);
  }
  std::printf("wavelength spacing  = min %.17g nm, max %.17g nm\n", dmin, dmax);

  // the reconstruction reads the native rows; require that the production lookup
  // returns exactly those rows at the knots (it must: fw = fh = 1 there)
  double worstRow = 0.0;
  for (std::size_t k = 0; k < prof.lam.size(); ++k)
    for (std::size_t j = 0; j < prof.hgt.size(); ++j)
      worstRow = std::max(worstRow, std::fabs(abs.getOpticalDepth(prof.lam[k], prof.hgt[j]) -
                                              prof.tau[k][j]));
  std::printf("row identity        = max |getOpticalDepth(knot) - native row| = %.17g\n",
              worstRow);
  std::printf("\n");


  // ======================= B. constant-slope sanity ========================
  double const kSlope = 1.0e-4;
  std::string const fnK = "/tmp/cp2fa_synth_const.dat";
  if (!writeSynthetic(fnK, false, kSlope, 0.0)) {
    std::printf("STOP: cannot write the synthetic profile\n");
    return 1;
  }
  AtmosphericAbsorption absK;
  absK.load(fnK);
  auto const hrK = absK.getHeightRange();
  double const HK = hrK.second;

  std::printf("--- B. constant-slope profile k = %.17g 1/m (h in [%.17g, %.17g]) ---\n", kSlope,
              hrK.first, hrK.second);
  std::printf("%-12s %-12s %-12s %-10s %-18s %-18s %-12s\n", "h", "window_lo", "window_hi",
              "denomin.", "expected", "computed", "residual");
  double probes[18];
  double const offs[9] = {-26.0, -25.0, -24.0, -1.0, 0.0, 1.0, 24.0, 25.0, 26.0};
  for (int i = 0; i < 9; ++i) probes[i] = offs[i];
  for (int i = 0; i < 9; ++i) probes[9 + i] = HK + offs[i];
  double worstBsyn = 0.0;
  for (int i = 0; i < 18; ++i) {
    double const h = probes[i];
    Win const w = windowOf(h, hrK.first, hrK.second);
    double const expected = w.ok ? kSlope : 0.0;  // branch-derived expectation
    double const got = alphaSrc(absK, 480.0, h, hrK.first, hrK.second);
    double const resid = std::fabs(got - expected);
    worstBsyn = std::max(worstBsyn, resid);
    ++g_points;
    if (!(resid <= 1.0e-14)) ++g_fail;
    std::printf("%-12.17g %-12.17g %-12.17g %-10.17g %-18.17g %-18.17g %-12.3e%s\n", h, w.lo,
                w.hi, w.d, expected, got, resid, (resid <= 1.0e-14) ? "" : "  FAIL");
  }
  std::printf("worst residual = %.17g  (gate 1e-14 absolute)\n\n", worstBsyn);

  // ======================= C. zero-slope sanity ============================
  std::string const fnZ = "/tmp/cp2fa_synth_zero.dat";
  if (!writeSynthetic(fnZ, true, 0.0, 0.25)) {
    std::printf("STOP: cannot write the synthetic profile\n");
    return 1;
  }
  AtmosphericAbsorption absZ;
  absZ.load(fnZ);
  auto const hrZ = absZ.getHeightRange();
  double const HZ = hrZ.second;
  std::printf("--- C. zero-slope profile (tau = 0.25 constant) ---\n");
  std::printf("%-12s %-12s %-12s %-10s %-18s %-18s %-12s\n", "h", "window_lo", "window_hi",
              "denomin.", "expected", "computed", "residual");
  double worstC = 0.0;
  for (int i = 0; i < 18; ++i) {
    double const h = (i < 9) ? probes[i] : HZ + offs[i - 9];
    Win const w = windowOf(h, hrZ.first, hrZ.second);
    double const got = alphaSrc(absZ, 480.0, h, hrZ.first, hrZ.second);
    double const resid = std::fabs(got);
    worstC = std::max(worstC, resid);
    ++g_points;
    if (!(resid <= 1.0e-14)) ++g_fail;
    std::printf("%-12.17g %-12.17g %-12.17g %-10.17g %-18.17g %-18.17g %-12.3e%s\n", h, w.lo,
                w.hi, w.d, 0.0, got, resid, (resid <= 1.0e-14) ? "" : "  FAIL");
  }
  std::printf("worst residual = %.17g  (gate 1e-14 absolute)\n\n", worstC);

  // ======================= D. native slopes of the real profile ============
  double sMin = 1e300, sMax = -1e300;
  std::size_t sMinK = 0, sMinJ = 0;
  long nNegative = 0;
  for (std::size_t k = 0; k < prof.lam.size(); ++k) {
    for (std::size_t j = 0; j + 1 < prof.hgt.size(); ++j) {
      double const s = (prof.tau[k][j + 1] - prof.tau[k][j]) / (prof.hgt[j + 1] - prof.hgt[j]);
      if (s < sMin) {
        sMin = s;
        sMinK = k;
        sMinJ = j;
      }
      sMax = std::max(sMax, s);
      if (s < 0.0) ++nNegative;
    }
  }
  std::printf("--- D. native slopes of the real profile (%zu bins x %zu knots) ---\n",
              prof.hgt.size() - 1, prof.lam.size());
  std::printf("min slope = %.17g 1/m at lambda=%.17g nm, bin [%.17g, %.17g] m\n", sMin,
              prof.lam[sMinK], prof.hgt[sMinJ], prof.hgt[sMinJ + 1]);
  std::printf("max slope = %.17g 1/m\n", sMax);
  std::printf("negative-slope entries = %ld of %zu\n",
              nNegative, (prof.hgt.size() - 1) * prof.lam.size());
  std::printf("\n");


  // ======================= E. fixed real-profile test set =================
  std::vector<double> heights;
  heights.reserve(51000);
  for (int i = 0; i <= int(hHi + 26.0); ++i) heights.push_back(double(i));  // 1 m grid
  for (std::size_t j = 0; j < prof.hgt.size(); ++j) {
    double const base[3] = {prof.hgt[j] - 25.0, prof.hgt[j], prof.hgt[j] + 25.0};
    for (int t = 0; t < 3; ++t) {
      heights.push_back(base[t]);
      heights.push_back(base[t] - 1.0e-6);
      heights.push_back(base[t] + 1.0e-6);
    }
  }
  double const below[4] = {-26.0, -25.0, -24.0, -1.0};
  for (int t = 0; t < 4; ++t) heights.push_back(below[t]);
  std::sort(heights.begin(), heights.end());
  heights.erase(std::unique(heights.begin(), heights.end()), heights.end());

  struct Mid {
    int m, a, b;
  };
  std::vector<double> lams;
  std::vector<Mid> mids;
  for (std::size_t k = 0; k < prof.lam.size(); ++k) lams.push_back(prof.lam[k]);
  for (std::size_t k = 0; k + 1 < prof.lam.size(); ++k) {
    Mid mid;
    mid.m = int(lams.size());
    mid.a = int(k);
    mid.b = int(k + 1);
    mids.push_back(mid);
    lams.push_back(0.5 * (prof.lam[k] + prof.lam[k + 1]));
  }
  int const idx150 = int(lams.size());
  lams.push_back(150.0);
  int const idx900 = int(lams.size());
  lams.push_back(900.0);

  std::vector<std::vector<double> > slopes(lams.size());
  std::vector<double> lamUsedV(lams.size(), 0.0);
  for (std::size_t i = 0; i < lams.size(); ++i) {
    slopesAt(prof, lams[i], slopes[i], lamUsedV[i]);
  }
  std::printf("    lambda clamp (getOpticalDepth semantics): 150 nm -> %.17g nm ; "
              "900 nm -> %.17g nm\n",
              lamUsedV[std::size_t(idx150)], lamUsedV[std::size_t(idx900)]);

  std::printf("--- E. fixed test set: %zu heights x %zu wavelengths = %zu queries ---\n",
              heights.size(), lams.size(), heights.size() * lams.size());
  std::printf("    heights: 1 m grid 0..%.17g m; every knot shift {-25,0,+25} m with +-1e-6 m;\n",
              hHi + 26.0);
  std::printf("    plus below-profile points -26,-25,-24,-1 m\n");
  std::printf("    wavelengths: 105 native knots, 104 adjacent midpoints, 150 nm, 900 nm\n");
  std::printf("    gate: |a-b| <= 1e-14 + 1e-10*|a|\n\n");

  double worst = 0.0, worstLam = 0.0, worstH = 0.0, worstA = 0.0, worstB = 0.0;
  double worstLin = 0.0, worstLinLam = 0.0, worstLinH = 0.0;
  long nFailE = 0, nClamp = 0, nFailLin = 0;
  double worstRaw = 0.0, worstRawLam = 0.0, worstRawH = 0.0;
  double failHmin = 1e300, failHmax = -1e300, worstClampResid = 0.0;
  double worstCond = 0.0, worstCondLam = 0.0, worstCondH = 0.0, worstCondD = 0.0;
  double worstIll = 0.0, worstIllLam = 0.0, worstIllH = 0.0, worstIllD = 0.0;
  std::vector<double> aRow(lams.size()), bRow(lams.size());

  for (std::size_t hi = 0; hi < heights.size(); ++hi) {
    double const h = heights[hi];
    for (std::size_t i = 0; i < lams.size(); ++i) {
      double raw = 0.0;
      double const a = alphaSrc(abs, lams[i], h, hLo, hHi, &raw);
      double const b = alphaKernel(slopes[i], prof.hgt, h, hLo, hHi);
      aRow[i] = a;
      bRow[i] = b;
      double const resid = std::fabs(a - b);
      double const gate = 1.0e-14 + 1.0e-10 * std::fabs(a);
      ++g_points;
      if (!(resid <= gate)) {
        ++g_fail;
        ++nFailE;
        failHmin = std::min(failHmin, h);
        failHmax = std::max(failHmax, h);
        if (nFailE <= 12)
          std::printf("    FAIL lam=%.17g h=%.17g a=%.17g b=%.17g |a-b|=%.3e gate=%.3e\n",
                      lams[i], h, a, b, resid, gate);
      }
      if (resid > worst) {
        worst = resid;
        worstLam = lams[i];
        worstH = h;
        worstA = a;
        worstB = b;
      }
      {
        Win const wq = windowOf(h, hLo, hHi);
        if (wq.ok) {
          if (wq.d >= 1.0e-3) {
            if (resid > worstCond) {
              worstCond = resid;
              worstCondLam = lams[i];
              worstCondH = h;
              worstCondD = wq.d;
            }
          } else if (resid > worstIll) {
            worstIll = resid;
            worstIllLam = lams[i];
            worstIllH = h;
            worstIllD = wq.d;
          }
        }
      }
      if (raw < 0.0) {
        ++nClamp;
        worstClampResid = std::max(worstClampResid, resid);
        if (raw < worstRaw) {
          worstRaw = raw;
          worstRawLam = lams[i];
          worstRawH = h;
        }
      }
    }
    // wavelength linearity at fixed height: alpha(mid) vs 0.5*(alpha(k)+alpha(k+1))
    for (std::size_t mi = 0; mi < mids.size(); ++mi) {
      Mid const& m = mids[mi];
      double const lhs = aRow[std::size_t(m.m)];
      double const rhs = 0.5 * (aRow[std::size_t(m.a)] + aRow[std::size_t(m.b)]);
      double const resid = std::fabs(lhs - rhs);
      double const gate = 1.0e-14 + 1.0e-10 * std::fabs(lhs);
      ++g_points;
      if (!(resid <= gate)) {
        ++g_fail;
        ++nFailLin;
      }
      if (resid > worstLin) {
        worstLin = resid;
        worstLinLam = lams[std::size_t(m.m)];
        worstLinH = h;
      }
      // and the reconstruction must be linear too
      double const lhsB = bRow[std::size_t(m.m)];
      double const rhsB = 0.5 * (bRow[std::size_t(m.a)] + bRow[std::size_t(m.b)]);
      double const residB = std::fabs(lhsB - rhsB);
      ++g_points;
      if (!(residB <= 1.0e-14 + 1.0e-10 * std::fabs(lhsB))) {
        ++g_fail;
        ++nFailLin;
      }
    }
  }
  std::printf("    kernel residual: worst %.17g at lam=%.17g nm, h=%.17g m  (a=%.17g b=%.17g)\n",
              worst, worstLam, worstH, worstA, worstB);
  std::printf("    lambda linearity: worst %.17g at lam=%.17g nm, h=%.17g m\n", worstLin,
              worstLinLam, worstLinH);
  std::printf("    worst residual, windows with denominator >= 1e-3 m: %.17g "
              "(lam=%.17g, h=%.17g, Delta=%.17g)\n",
              worstCond, worstCondLam, worstCondH, worstCondD);
  std::printf("    worst residual, windows with denominator <  1e-3 m: %.17g "
              "(lam=%.17g, h=%.17g, Delta=%.17g)\n",
              worstIll, worstIllLam, worstIllH, worstIllD);
  std::printf("    clamp-active queries (raw window average < 0): %ld\n", nClamp);
  if (nClamp) {
    std::printf("      most negative raw average %.17g at lam=%.17g nm, h=%.17g m\n", worstRaw,
                worstRawLam, worstRawH);
    std::printf("      worst kernel residual among clamp-active queries: %.17g\n",
                worstClampResid);
  }
  std::printf("    FAILED queries: %ld", nFailE);
  if (nFailE) std::printf("  (failing heights span [%.17g, %.17g] m)", failHmin, failHmax);
  std::printf("\n");

  // high-precision cross-check at the worst query: separates a representation
  // error from a rounding error inside the production expression
  {
    long double const aLd = alphaLD(prof, worstLam, worstH, hLo, hHi);
    Win const wq = windowOf(worstH, hLo, hHi);
    std::printf("    long-double cross-check at lam=%.17g nm, h=%.17g m, window=[%.17g, "
                "%.17g], denominator=%.17g:\n",
                worstLam, worstH, wq.lo, wq.hi, wq.d);
    std::printf("      alpha_long_double = %.21Lg\n", aLd);
    std::printf("      alpha_src  (double) = %.17g   |src  - ld| = %.3Le\n", worstA,
                std::fabs((long double)worstA - aLd));
    std::printf("      alpha_kern (double) = %.17g   |kern - ld| = %.3Le\n", worstB,
                std::fabs((long double)worstB - aLd));
  }

  // Delta-scaling probe: if the source discrepancy is rounding inside the
  // production expression, |alpha_src - exact| must grow like 1/Delta while the
  // reconstruction stays exactly at the last-bin slope.
  {
    double const lamProbe = 255.0;
    std::vector<double> sProbe;
    double lu = 0.0;
    slopesAt(prof, lamProbe, sProbe, lu);
    double const sLast = sProbe.back();
    std::printf("\n    Delta-scaling probe at lam=%.17g nm (exact value = last-bin slope "
                "%.17g 1/m):\n",
                lamProbe, sLast);
    std::printf("      %-12s %-18s %-20s %-20s %-14s %-14s\n", "delta", "Delta",
                "alpha_src", "alpha_kernel", "|src-exact|", "|src-exact|*Delta");
    double const deltas[4] = {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};
    for (int t = 0; t < 4; ++t) {
      double const h = hHi + 25.0 - deltas[t];
      Win const w = windowOf(h, hLo, hHi);
      double const a = alphaSrc(abs, lamProbe, h, hLo, hHi);
      double const b = alphaKernel(sProbe, prof.hgt, h, hLo, hHi);
      double const err = std::fabs(a - sLast);
      std::printf("      %-12.17g %-18.17g %-20.17g %-20.17g %-14.3e %-14.3e\n", deltas[t],
                  w.d, a, b, err, err * w.d);
    }
  }
  std::printf("\n");


  // ======================= F. boundary / support tables ====================
  int lam480 = -1;
  for (std::size_t i = 0; i < lams.size(); ++i)
    if (lams[i] == 480.0) lam480 = int(i);
  std::printf("--- F. boundary tables at lambda = 480 nm (source vs reconstruction) ---\n");

  std::printf("  ground band: h, window, denominator, alpha_src, alpha_kernel\n");
  std::printf("  %-16s %-14s %-14s %-12s %-18s %-18s\n", "h", "lo", "hi", "denomin.",
              "alpha_src", "alpha_kernel");
  std::vector<double> hGround;
  for (int i = -30; i <= 30; ++i) hGround.push_back(double(i));
  double const eps = 1.0e-6;
  double const edgeG[6] = {-25.0 - eps, -25.0 + eps, -eps, eps, 25.0 - eps, 25.0 + eps};
  for (int i = 0; i < 6; ++i) hGround.push_back(edgeG[i]);
  std::sort(hGround.begin(), hGround.end());
  for (std::size_t i = 0; i < hGround.size(); ++i) {
    double const h = hGround[i];
    Win const w = windowOf(h, hLo, hHi);
    double const a = alphaSrc(abs, 480.0, h, hLo, hHi);
    double const b = alphaKernel(slopes[std::size_t(lam480)], prof.hgt, h, hLo, hHi);
    std::printf("  %-16.17g %-14.17g %-14.17g %-12.17g %-18.17g %-18.17g\n", h, w.lo, w.hi, w.d,
                a, b);
  }

  std::printf("\n  top band: h, window, denominator, alpha_src, alpha_kernel\n");
  std::printf("  %-16s %-14s %-14s %-12s %-18s %-18s\n", "h", "lo", "hi", "denomin.",
              "alpha_src", "alpha_kernel");
  std::vector<double> hTop;
  for (int i = -30; i <= 30; ++i) hTop.push_back(hHi + double(i));
  double const edgeT[6] = {hHi - 25.0 - eps, hHi - 25.0 + eps, hHi - eps, hHi + eps,
                           hHi + 25.0 - eps, hHi + 25.0 + eps};
  for (int i = 0; i < 6; ++i) hTop.push_back(edgeT[i]);
  std::sort(hTop.begin(), hTop.end());
  for (std::size_t i = 0; i < hTop.size(); ++i) {
    double const h = hTop[i];
    Win const w = windowOf(h, hLo, hHi);
    double const a = alphaSrc(abs, 480.0, h, hLo, hHi);
    double const b = alphaKernel(slopes[std::size_t(lam480)], prof.hgt, h, hLo, hHi);
    std::printf("  %-16.17g %-14.17g %-14.17g %-12.17g %-18.17g %-18.17g\n", h, w.lo, w.hi, w.d,
                a, b);
  }

  // empirical support: where the source value is non-zero (1 mm scan)
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
  std::printf("\n  empirical support (1 mm scan, non-zero source value):\n");
  std::printf("    scan [-30, 30] m    : first non-zero at h = %.17g m; still non-zero at the "
              "scan limit %.17g m\n",
              loNZ, hiNZ);
  std::printf("    scan [H-30, H+30] m : non-zero at the scan limit %.17g m; last non-zero at "
              "h = %.17g m\n",
              loNZT, hiNZT);
  std::printf("\n");

  // ======================= G. summary ======================================
  auto const tEnd = std::chrono::steady_clock::now();
  double const secs = std::chrono::duration<double>(tEnd - tStart).count();
  std::printf("--- G. summary ---\n");
  std::printf("queries evaluated      : %d\n", g_points);
  std::printf("kernel queries failed  : %ld (of %zu kernel queries)\n", nFailE,
              heights.size() * lams.size());
  std::printf("lambda-linearity fails : %ld (of %zu checks: 2 per midpoint per height)\n",
              nFailLin, 2 * mids.size() * heights.size());
  std::printf("failed queries total   : %d\n", g_fail);
  std::printf("worst kernel residual  : %.17g at lam=%.17g nm, h=%.17g m\n", worst, worstLam,
              worstH);
  std::printf("worst lambda-linearity : %.17g\n", worstLin);
  std::printf("min native slope       : %.17g 1/m ; negative entries: %ld\n", sMin, nNegative);
  std::printf("clamp-active queries   : %ld\n", nClamp);
  std::printf("runtime                : %.3f s\n", secs);
  bool const pass = (g_fail == 0);
  std::printf("VERDICT: %s\n", pass ? "PASS -- production extinction equals "
                                    "sum_j s_j(lambda) w_j(h) on the tested set"
                                  : "FAIL -- see the failing queries above (report and stop)");
  std::printf("=== end of 2f-a ===\n");
  return pass ? 0 : 1;
}

