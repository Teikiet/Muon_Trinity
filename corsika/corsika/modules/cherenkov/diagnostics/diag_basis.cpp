// Diagnostic: native-altitude-basis weights for the CORRECTED extinction callback.
//   tau(lambda) = sum_k s_k(lambda) * W_k,  W_k = int_ray phi_k(h(s)) ds
// phi_k is the exact trapezoidal basis implied by the callback's 25 m smoothing of a
// piecewise-linear column. Verified BEFORE any coarsening. Production routing untouched.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static std::vector<double> g_hk;                  // native knots (0..50 km, 1 km)
static double alphaCB(double lam, double h) {
  double const dh = 50.0;
  auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h - dh*0.5), hi = std::min(hr.second, h + dh*0.5);
  double const a = hi - lo; if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam, hi) - g_ab.getOpticalDepth(lam, lo)) / a);
}
// native slope per interval k, at lambda (column is piecewise linear on the knots)
static double slopeK(int k, double lam) {
  double const t0 = g_ab.getOpticalDepth(lam, g_hk[k]);
  double const t1 = g_ab.getOpticalDepth(lam, g_hk[k+1]);
  return (t1 - t0) / (g_hk[k+1] - g_hk[k]);
}
// exact trapezoidal basis for interval k
// Exact basis for the ACTUAL callback. Interior: linear ramps over the 25 m skirts.
// At the profile boundaries the clamped denominator (h_hi - h_lo) exactly cancels the
// skirt truncation for a linear column, so alpha_cb is CONSTANT there and the boundary
// basis saturates at 1.0 (one-sided skirt). Verified pointwise before use.
static double phiK(int k, double h) {
  double const a = g_hk[k], b = g_hk[k+1], w = 25.0;
  bool const first = (k == 0), last = (k + 2 == int(g_hk.size()));
  if (h <= a - w || h >= b + w) return 0.0;
  if (h < a + w) return first ? 1.0 : (h - (a - w)) / (2.0*w);
  if (h > b - w) return last ? 1.0 : ((b + w) - h) / (2.0*w);
  return 1.0;
}
struct Out { bool hit{false}; double tau{0.0}; double L{0.0}; std::vector<double> W; long steps{0}; };
static Out traceWeights(double h0, double theta0, double fPhase, double lam, double ds) {
  Out o; o.W.assign(int(g_hk.size())-1, 0.0);
  double r = R + h0, th = theta0;
  if (r <= R) { o.hit = true; return o; }
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const h = rr - R;
    double const Nv = N0 * std::exp(-h / HSC) * fPhase;
    double const np = -Nv / HSC, nn = 1.0 + Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c; dt = -s * (1.0/rr + (nn > 1e-12 ? np/nn : 0.0));
  };
  for (long i = 0; i < 400000000L; ++i) {
    double dr1,dt1,dr2,dt2,dr3,dt3,dr4,dt4;
    der(r, th, dr1, dt1);
    der(r+0.5*ds*dr1, th+0.5*ds*dt1, dr2, dt2);
    der(r+0.5*ds*dr2, th+0.5*ds*dt2, dr3, dt3);
    der(r+ds*dr3, th+ds*dt3, dr4, dt4);
    double const rPrev = r;
    r += ds/6.0*(dr1+2*dr2+2*dr3+dr4);
    th += ds/6.0*(dt1+2*dt2+2*dt3+dt4);
    double const hMid = std::max(0.0, 0.5*(rPrev + r) - R);
    for (std::size_t k = 0; k + 1 < g_hk.size(); ++k) o.W[k] += phiK(int(k), hMid) * ds;
    o.tau += alphaCB(lam, hMid) * ds;
    o.L += ds; ++o.steps;
    if (r <= R) { o.hit = true; return o; }
  }
  return o;
}

// basis tau with weights accumulated at step ds
static double tauBasis(double h0, double theta0, double fPhase, double lam, double ds) {
  Out const o = traceWeights(h0, theta0, fPhase, lam, ds);
  double t = 0.0;
  for (std::size_t k = 0; k + 1 < g_hk.size(); ++k) t += slopeK(int(k), lam) * o.W[k];
  return t;
}
// direct callback integration at step ds (INDEPENDENT of the basis)
static double tauDirect(double h0, double theta0, double fPhase, double lam, double ds) {
  return traceWeights(h0, theta0, fPhase, lam, ds).tau;
}
static double thetaFor(double h0, double mu, double f) {
  double const nObs = 1.0 + N0*f, nEmit = 1.0 + N0*std::exp(-h0/HSC)*f;
  double const p = nObs*R*std::sqrt(std::max(0.0, 1.0-mu*mu));
  double const s0 = p/(nEmit*(R+h0));
  if (s0 >= 1.0) return -1.0;
  return std::acos(-std::sqrt(1.0 - s0*s0));
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  auto const hr = g_ab.getHeightRange();
  for (double h = hr.first; h <= hr.second + 1e-9; h += 1000.0) g_hk.push_back(h);
  std::printf("knots=%zu", g_hk.size()); puts("");
  struct Probe { double h, mu, f, lam; const char* tag; };
  Probe const probes[] = {
    {7500.0,  0.700, 1.0347, 300.0, "previous worst"},
    {25000.0, 0.010, 1.0000, 400.0, "grazing"},
    {50000.0, 0.030, 1.0000, 400.0, "top boundary"},
    {500.0,   0.300, 1.0000, 400.0, "ground boundary"},
    {20000.0, 0.500, 0.9713, 400.0, "f endpoint low"},
    {20000.0, 0.500, 1.0347, 400.0, "f endpoint high"},
    {45000.0, 0.010, 1.0347, 690.0, "grazing + high f"},
  };
  // ---- (5) step refinement: reference uncertainty ----------------------------
  std::printf("step refinement (tau_basis and tau_direct vs tau_converged):"); puts("");
  double worstRef = 0.0;
  for (Probe const& pr : probes) {
    double const th0 = thetaFor(pr.h, pr.mu, pr.f);
    if (th0 < 0.0) { std::printf("   %-18s cannot launch", pr.tag); puts(""); continue; }
    double const L = traceWeights(pr.h, th0, pr.f, pr.lam, 1.0).L;
    double const dsc = std::max(0.05, L / 2.0e6);      // bounded converged reference
    double const conv = tauDirect(pr.h, th0, pr.f, pr.lam, dsc);
    double const b2 = tauBasis(pr.h, th0, pr.f, pr.lam, 2.0), b1 = tauBasis(pr.h, th0, pr.f, pr.lam, 1.0);
    double const b5 = tauBasis(pr.h, th0, pr.f, pr.lam, 0.5), b25 = tauBasis(pr.h, th0, pr.f, pr.lam, 0.25);
    double const d2 = tauDirect(pr.h, th0, pr.f, pr.lam, 2.0), d1 = tauDirect(pr.h, th0, pr.f, pr.lam, 1.0);
    double const d5 = tauDirect(pr.h, th0, pr.f, pr.lam, 0.5), d25 = tauDirect(pr.h, th0, pr.f, pr.lam, 0.25);
    std::printf("   %-18s conv=%.6f (ds=%.3f m)  |b(2)-c|=%.2e |b(1)-c|=%.2e |b(.5)-c|=%.2e |b(.25)-c|=%.2e",
                pr.tag, conv, dsc, std::abs(b2-conv), std::abs(b1-conv), std::abs(b5-conv), std::abs(b25-conv)); puts("");
    std::printf("   %-18s                      |d(2)-c|=%.2e |d(1)-c|=%.2e |d(.5)-c|=%.2e |d(.25)-c|=%.2e",
                "", std::abs(d2-conv), std::abs(d1-conv), std::abs(d5-conv), std::abs(d25-conv)); puts("");
    double const dT = std::abs(std::exp(-b1) - std::exp(-conv));
    std::printf("   %-18s combined |dT| (basis ds=1.0 vs converged) = %.4e", pr.tag, dT); puts("");
    worstRef = std::max(worstRef, dT);
  }
  std::printf("WORST combined |dT| over probes = %.4e   (budget 1e-3)", worstRef); puts("");
  { std::ifstream f("/proc/self/status"); std::string l;
    while (std::getline(f, l)) if (l.rfind("VmHWM",0)==0) { std::printf("VmHWM: %s", l.c_str()); puts(""); } }
  return 0;
}
