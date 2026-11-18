// PHASE 2: reference validation. Analytic checks on the integrator, breakpoint
// verification against the callback, and grazing admissibility along the whole path.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static double nP(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
// generic b-integral with an injected alpha
static double bInt(double b, double f, double rTop, std::function<double(double)> const& al, double dr) {
  double acc = 0.0, r = R;
  while (r < rTop - 1e-12) {
    double const st = std::min(dr, rTop - r);
    double const rm = r + 0.5*st;
    double const s = b/(nP(rm,f)*rm);
    if (s >= 1.0) return acc;
    acc += al(rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st;
    r += st;
  }
  return acc;
}
// production callback alpha, with the quadrature split at native knots +/-25 m
static double alphaCB(double lam, double h) {
  double const dh = 50.0; auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h-dh*0.5), hi = std::min(hr.second, h+dh*0.5);
  double const a = hi-lo; if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam,hi) - g_ab.getOpticalDepth(lam,lo))/a);
}
static std::vector<double> g_hk;
static double splitInt(double b, double f, double lam, double rTop, double dr) {
  std::vector<double> cuts; for (double h : g_hk) { double c1 = R+h-25.0, c2 = R+h+25.0; if (c1 > R && c1 < rTop) cuts.push_back(c1); if (c2 > R && c2 < rTop) cuts.push_back(c2); }
  cuts.push_back(rTop); std::sort(cuts.begin(), cuts.end());
  double acc = 0.0, r = R;
  for (double c : cuts) {
    if (c <= R) continue;
    while (r < c - 1e-12) {
      double const st = std::min(dr, c - r);
      double const rm = r + 0.5*st;
      double const s = b/(nP(rm,f)*rm);
      if (s >= 1.0) return acc;
      acc += alphaCB(lam, rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st;
      r += st;
    }
  }
  return acc;
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  auto const hr = g_ab.getHeightRange();
  for (double h = hr.first; h <= hr.second+1e-9; h += 1000.0) g_hk.push_back(h);
  double const f = 1.0;
  // ---- analytic checks (metric stated: T unless noted) -----------------------
  std::printf("ANALYTIC CHECKS (integrator, synthetic alpha; metric = transmission T)"); puts("");
  {
    auto zero = [](double) { return 0.0; };
    double const t = std::exp(-bInt(0.3*R, f, R+20000.0, zero, 5.0));
    std::printf("  zero extinction        T = %.15f   (expect 1)            |dT| = %.2e", t, std::abs(t-1.0)); puts("");
  }
  {
    auto c = [](double) { return 1e-4; };
    double const t = std::exp(-bInt(0.3*R, f, R, c, 5.0));
    std::printf("  zero path (h=0)        T = %.15f   (expect 1)            |dT| = %.2e", t, std::abs(t-1.0)); puts("");
  }
  {
    // constant alpha, RADIAL ray: path length = h, so tau = alpha*h exactly
    auto c = [](double) { return 1e-4; };
    for (double h : {1000.0, 20000.0, 50000.0}) {
      double const tau = bInt(0.0, f, R+h, c, 5.0);
      double const want = 1e-4*h;
      std::printf("  const alpha radial h=%6.0f  tau=%.9f  alpha*L=%.9f  |dT|=%.2e", h, tau, want,
                  std::abs(std::exp(-tau)-std::exp(-want))); puts("");
    }
  }
  {
    // piecewise-constant layers with independently known thicknesses (radial ray)
    std::vector<double> edges = {0.0, 1000.0, 5000.0, 20000.0};
    std::vector<double> a = {3e-4, 1e-4, 2e-5};
    auto piece = [&](double h) { for (std::size_t k = 0; k+1 < edges.size(); ++k) if (h >= edges[k] && h < edges[k+1]) return a[k]; return 0.0; };
    double want = 0.0; for (std::size_t k = 0; k+1 < edges.size(); ++k) want += a[k]*(edges[k+1]-edges[k]);
    double const tau = bInt(0.0, f, R+edges.back(), piece, 1.0);
    std::printf("  piecewise-constant radial  tau=%.9f  sum(a_j*dh_j)=%.9f  |dT|=%.2e", tau, want,
                std::abs(std::exp(-tau)-std::exp(-want))); puts("");
  }
  {
    // no-refraction limit: n == 1 -> straight-line geometry. Compare the b-integral
    // with n=1 against the chord integral along the same straight path.
    auto c = [](double) { return 1e-4; };
    double const hE = 20000.0, ct = 0.5;                 // cos(theta) of the straight ray
    double const sE = std::sqrt(1.0-ct*ct);              // sin(theta_e)
    double const bq = (R+hE)*sE;                         // impact parameter of the straight line
    double const L = std::sqrt((R+hE)*(R+hE)-bq*bq) - std::sqrt(R*R-bq*bq);  // SPHERICAL chord
    double const want = 1e-4*L;
    // with n=1 the conserved b-integral reduces to int alpha dr / sqrt(1-(b/r)^2), b=r_e sin(theta_e)
    double const b = bq;
    double acc = 0.0; double const dr = 0.5;
    for (double r = R; r < R+hE; r += dr) { double const s = b/r; acc += c(0.0)/std::sqrt(std::max(1e-300,1.0-s*s))*dr; }
    std::printf("  no-refraction straight  tau=%.6f  alpha*L=%.6f  |dT|=%.2e  (L=%.1f m)", acc, want,
                std::abs(std::exp(-acc)-std::exp(-want)), L); puts("");
  }
  // ---- breakpoint verification: does the +/-25 m split capture the callback? ---
  std::printf("BREAKPOINT VERIFICATION (metric = |d tau|, then |dT|)"); puts("");
  for (double const ratio : {0.9, 0.999, 0.9999}) {
    double const b = ratio*(nP(R,f)*R);
    double const a1 = splitInt(b, f, 400.0, R+20000.0, 25.0);   // split at knots +/-25 m
    double const a2 = splitInt(b, f, 400.0, R+20000.0, 0.5);    // same split, finer
    double const a3 = splitInt(b, f, 400.0, R+20000.0, 0.01);   // same split, finest
    double const u1 = bInt(b, f, R+20000.0, [](double h){ return alphaCB(400.0,h); }, 0.01); // no split, uniform
    std::printf("  ratio=%.4f  split(25)=%.9f split(0.5)=%.9f split(0.01)=%.9f | unsplit(0.01)=%.9f  |dT split-vs-unsplit|=%.2e",
                ratio, a1, a2, a3, u1, std::abs(std::exp(-a3)-std::exp(-u1))); puts("");
  }
  // ---- grazing admissibility along the WHOLE path (not just the detector) -----
  std::printf("GRAZING ADMISSIBILITY: min of n(r) r over [R, R+h]"); puts("");
  for (double h : {1000.0, 20000.0, 50000.0}) {
    double mn = 1e300; double rAt = R;
    for (double r = R; r <= R+h+1e-9; r += 1.0) { double const v = nP(r,f)*r; if (v < mn) { mn = v; rAt = r; } }
    std::printf("  h=%6.0f m  min(n r) = %.3f at r-R = %8.1f m   (= n(R) R ? %s)",
                h, mn, rAt-R, (std::abs(mn - nP(R,f)*R) < 1e-6) ? "yes" : "NO"); puts("");
  }
  return 0;
}
