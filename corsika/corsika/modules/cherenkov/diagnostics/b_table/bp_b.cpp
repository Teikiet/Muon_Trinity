// bp_b.cpp -- callback breakpoint verification (h >= 0 only) and the lambda edge
// behaviour of AtmosphericAbsorption. No table build; fast.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
static AtmosphericAbsorption g_ab;
static std::vector<double> g_knots;

static void parseKnots(const char* fn) {
  std::ifstream f(fn);
  size_t nw = 0, nh = 0;
  f >> nw >> nh;
  std::vector<double> wl(nw);
  for (size_t i = 0; i < nw; ++i) f >> wl[i];
  g_knots.assign(nh, 0.0);
  for (size_t i = 0; i < nh; ++i) f >> g_knots[i];
}
// production callback, verbatim
static double alphaOf(double lam, double h) {
  double const dh = 50.0;
  auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h - dh * 0.5), hi = std::min(hr.second, h + dh * 0.5);
  double const a = hi - lo;
  if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam, hi) - g_ab.getOpticalDepth(lam, lo)) / a);
}

int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  parseKnots(argv[1]);
  auto const wr = g_ab.getWavelengthRange();
  auto const hr = g_ab.getHeightRange();
  std::printf("profile: %zu lambda knots [%.1f, %.1f] nm, %zu height knots [%.1f, %.1f] m\n",
              size_t(1), wr.first, wr.second, g_knots.size(), hr.first, hr.second);
  {
    std::ifstream f(argv[1]);
    size_t nw = 0, nh = 0;
    f >> nw >> nh;
    std::printf("  native lambda knots: %zu, spacing %.2f nm\n", nw, (wr.second - wr.first) / (nw - 1));
  }

  // ---- breakpoint verification, restricted to the callback's domain [0, hmax] --------
  std::vector<double> bp;
  for (double k : g_knots) {
    bp.push_back(k - 25.0);
    bp.push_back(k + 25.0);
  }
  bp.push_back(0.0);
  bp.push_back(25.0);
  bp.push_back(hr.second - 25.0);
  bp.push_back(hr.second);
  std::sort(bp.begin(), bp.end());
  bp.erase(std::unique(bp.begin(), bp.end()), bp.end());
  std::vector<double> bpd;
  for (double v : bp)
    if (v >= 0.0 && v <= hr.second) bpd.push_back(v);
  double worst = 0.0, arg = 0.0;
  int nseg = 0;
  for (size_t i = 0; i + 1 < bpd.size(); ++i) {
    double const a = bpd[i], b = bpd[i + 1];
    if (b - a < 1.0) continue;
    ++nseg;
    double const ya = alphaOf(690.0, a), yb = alphaOf(690.0, b);
    for (int q = 1; q <= 3; ++q) {
      double const h = a + (b - a) * double(q) / 4.0;
      double const lin = ya + (yb - ya) * (h - a) / (b - a);
      double const scale = std::max(1e-30, std::abs(ya) + std::abs(yb));
      double const res = std::abs(alphaOf(690.0, h) - lin) / scale;
      if (res > worst) {
        worst = res;
        arg = h;
      }
    }
  }
  std::printf("\nbreakpoint set {knot +/- 25 m} within [0, hmax]: %zu points, %d segments tested\n",
              bpd.size(), nseg);
  std::printf("  alpha collinearity between consecutive breakpoints: worst relative residual %.3e at h=%.1f m\n",
              worst, arg);
  std::printf("  -> %s\n", (worst < 1e-6) ? "set is COMPLETE (alpha linear between predicted breakpoints)"
                                          : "NOT verified: residual above noise");

  // ---- lambda edge behaviour ---------------------------------------------------------
  std::printf("\nlambda edge behaviour of getOpticalDepth (h=18000 m):\n");
  double const h0 = 18000.0;
  for (double lam : {295.0, 299.9, 300.0, 305.0, 310.0, 695.0, 700.0, 700.1, 720.0, 850.0, 900.0, 1000.0}) {
    double const t = g_ab.getOpticalDepth(lam, h0);
    double const tc = g_ab.getOpticalDepth(std::min(std::max(lam, wr.first), wr.second), h0);
    std::printf("  lambda=%7.1f  tau=%12.8f   clamped(%6.1f)=%12.8f   |diff|=%9.3e\n", lam, t,
                std::min(std::max(lam, wr.first), wr.second), tc, std::abs(t - tc));
  }
  std::printf("\n  alpha(295) vs alpha(300) at h=18000: %.8e vs %.8e\n", alphaOf(295.0, h0),
              alphaOf(300.0, h0));
  std::printf("  alpha(720) vs alpha(700) at h=18000: %.8e vs %.8e\n", alphaOf(720.0, h0),
              alphaOf(700.0, h0));
  return 0;
}
