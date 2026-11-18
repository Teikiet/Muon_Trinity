// Diagnostic: per-photon DIRECT native-profile integration of tau along the curved ray.
// Production routing untouched. Purpose: (1) establish a converged reference with explicit
// convergence checks (incl. grazing rays and profile boundaries), (2) validate the combined
// |dT| <= 1e-3 criterion on the frozen dataset, (3) measure FULL-route throughput
// (trajectory + extinction lookup), (4) report allocated memory and initialization cost.
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

static double alphaCB(double lam, double h) {
  double const dh = 50.0;
  auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h - dh * 0.5);
  double const hi = std::min(hr.second, h + dh * 0.5);
  double const a = hi - lo;
  if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam, hi) - g_ab.getOpticalDepth(lam, lo)) / a);
}
struct Out { bool hit{false}; double tau{0.0}; double L{0.0}; long steps{0}; };
// full route: integrate the Bouguer trajectory AND the extinction together
static Out directTau(double h0, double theta0, double fPhase, double lam, double ds) {
  Out o;
  double r = R + h0, th = theta0;
  if (r <= R) { o.hit = true; return o; }
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const h = rr - R;
    double const Nv = N0 * std::exp(-h / HSC) * fPhase;
    double const np = -Nv / HSC, nn = 1.0 + Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c; dt = -s * (1.0 / rr + (nn > 1e-12 ? np / nn : 0.0));
  };
  for (long i = 0; i < 400000000L; ++i) {
    double dr1,dt1,dr2,dt2,dr3,dt3,dr4,dt4;
    der(r, th, dr1, dt1);
    der(r + 0.5*ds*dr1, th + 0.5*ds*dt1, dr2, dt2);
    der(r + 0.5*ds*dr2, th + 0.5*ds*dt2, dr3, dt3);
    der(r + ds*dr3, th + ds*dt3, dr4, dt4);
    double const rPrev = r;
    r += ds/6.0*(dr1+2*dr2+2*dr3+dr4);
    th += ds/6.0*(dt1+2*dt2+2*dt3+dt4);
    double const hMid = 0.5*(rPrev + r) - R;
    o.tau += alphaCB(lam, std::max(0.0, hMid)) * ds;
    o.L += ds; ++o.steps;
    if (r <= R) { o.hit = true; return o; }
  }
  return o;
}
static double thetaFor(double h0, double mu, double f) {
  double const nObs = 1.0 + N0 * f, nEmit = 1.0 + N0 * std::exp(-h0 / HSC) * f;
  double const p = nObs * R * std::sqrt(std::max(0.0, 1.0 - mu * mu));
  double const s0 = p / (nEmit * (R + h0));
  if (s0 >= 1.0) return -1.0;
  return std::acos(-std::sqrt(1.0 - s0 * s0));
}
int main(int argc, char** argv) {
  auto const t0 = std::chrono::steady_clock::now();
  g_ab.load(argv[1]);
  auto const t1 = std::chrono::steady_clock::now();
  if (!g_ab.isLoaded()) { std::printf("load failed"); return 1; }
  std::printf("init: profile load %.1f ms", std::chrono::duration<double,std::milli>(t1-t0).count());
  puts("");
  double const heldH[] = {2500.0, 7500.0, 15000.0, 25000.0, 45000.0};
  double const heldMu[] = {0.999, 0.7, 0.3, 0.12, 0.05, 0.01};
  double const heldF[] = {0.9713, 1.0, 1.0347};
  double const heldLam[] = {300.0, 302.5, 400.0, 402.5, 500.0, 502.5, 690.0, 700.0};
  // ---- (1) convergence of the direct integration, incl. grazing ---------------
  std::printf("convergence (tau at ds, and |tau - tau_converged| with ds/16):"); puts("");
  for (double h0 : {25000.0, 45000.0}) for (double mu : {0.7, 0.12, 0.01}) {
    double const th0 = thetaFor(h0, mu, 1.0);
    if (th0 < 0.0) { std::printf("   h=%.0f mu=%.3f  (cannot launch)", h0, mu); puts(""); continue; }
    double t1m = directTau(h0, th0, 1.0, 400.0, 1.0).tau;
    double t2m = directTau(h0, th0, 1.0, 400.0, 0.5).tau;
    double t3m = directTau(h0, th0, 1.0, 400.0, 0.25).tau;
    double t4m = directTau(h0, th0, 1.0, 400.0, 0.125).tau;
    double const conv = directTau(h0, th0, 1.0, 400.0, 0.125/16.0).tau;
    std::printf("   h=%6.0f mu=%.3f  tau(1)=%.6f tau(.5)=%.6f tau(.25)=%.6f tau(.125)=%.6f  conv=%.6f  |d(0.125-conv)|=%.2e",
                h0, mu, t1m, t2m, t3m, t4m, conv, std::abs(t4m - conv)); puts("");
  }
  // ---- (2) frozen dataset against the converged reference ---------------------
  double worst = 0.0; char info[160] = {0}; int nray = 0;
  for (double h0 : heldH) for (double mu : heldMu) for (double f : heldF) {
    double const th0 = thetaFor(h0, mu, f);
    if (th0 < 0.0) continue;
    for (double lam : heldLam) {
      double const tA = directTau(h0, th0, f, lam, 1.0).tau;   // candidate production step
      double const tR = directTau(h0, th0, f, lam, 0.25).tau;  // reference step
      double const dT = std::abs(std::exp(-tA) - std::exp(-tR));
      ++nray;
      if (dT > worst) { worst = dT; std::snprintf(info, sizeof(info), "h=%.0f mu=%.3f f=%.4f lam=%.1f tau=%.4f", h0, mu, f, lam, tR); }
    }
  }
  std::printf("frozen dataset: %d rays  WORST |dT| (ds=1.0 vs ds=0.25) = %.4e  [%s]  (budget 1e-3)",
              nray, worst, info); puts("");
  // ---- (3) FULL-route throughput: trajectory + extinction lookup --------------
  for (double mu : {0.9, 0.3, 0.05}) {
    double const th0 = thetaFor(20000.0, mu, 1.0);
    if (th0 < 0.0) continue;
    long const n = 20;
    auto const a0 = std::chrono::steady_clock::now();
    double acc = 0.0;
    for (long k = 0; k < n; ++k) acc += directTau(20000.0, th0, 1.0, 400.0, 0.5).tau;
    auto const a1 = std::chrono::steady_clock::now();
    double const us = std::chrono::duration<double,std::micro>(a1-a0).count()/double(n);
    Out const one = directTau(20000.0, th0, 1.0, 400.0, 0.5);
    std::printf("   full route mu=%.2f : %.1f us/photon  (steps=%ld, tau=%.4f, acc=%.3e)",
                mu, us, one.steps, one.tau, acc); puts("");
  }
  // ---- (4) allocated memory + init -------------------------------------------
  {
    std::ifstream f("/proc/self/status"); std::string l;
    while (std::getline(f, l)) if (l.rfind("VmHWM",0)==0 || l.rfind("VmRSS",0)==0) { std::printf("   %s", l.c_str()); puts(""); }
  }
  return 0;
}
