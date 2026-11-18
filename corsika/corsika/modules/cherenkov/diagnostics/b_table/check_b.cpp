// Verify the conserved-impact-parameter formulation against an independent, converged
// trajectory reference. Diagnostic only; production untouched.
//   tau(r_e,b,lam) = int_{r_obs}^{r_e} alpha(r,lam) / sqrt(1 - (b/(n(r,lam) r))^2) dr
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static double alphaCB(double lam, double h) {
  double const dh = 50.0; auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h-dh*0.5), hi = std::min(hr.second, h+dh*0.5);
  double const a = hi-lo; if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam,hi) - g_ab.getOpticalDepth(lam,lo))/a);
}
// phase index n(r) = 1 + N(r)*f
static double nPhase(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
// tau via the conserved-b 1-D integral in r (Simpson, fine)
static double tauByB(double rEmit, double b, double f, double lam, int nstep) {
  double acc = 0.0; double const dr = (rEmit - R)/nstep;
  for (int i = 0; i < nstep; ++i) {
    double const r = R + (i+0.5)*dr;
    double const s = b/(nPhase(r,f)*r);
    double const den = std::sqrt(std::max(1e-300, 1.0 - s*s));
    acc += alphaCB(lam, r-R)/den * dr;
  }
  return acc;
}
// independent trajectory reference: RK4 along the Bouguer ray, converged step
static double tauByTraj(double hEmit, double theta0, double f, double lam, double ds) {
  double r = R + hEmit, th = theta0, acc = 0.0;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const Nv = N0*std::exp(-(rr-R)/HSC)*f, np = -Nv/HSC, nn = 1.0+Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c; dt = -s*(1.0/rr + (np/nn));
  };
  for (long i = 0; i < 400000000L; ++i) {
    double a1,b1,a2,b2,a3,b3,a4,b4;
    der(r,th,a1,b1); der(r+0.5*ds*a1, th+0.5*ds*b1, a2,b2);
    der(r+0.5*ds*a2, th+0.5*ds*b2, a3,b3); der(r+ds*a3, th+ds*b3, a4,b4);
    double const rp = r;
    r += ds/6.0*(a1+2*a2+2*a3+a4); th += ds/6.0*(b1+2*b2+2*b3+b4);
    acc += alphaCB(lam, std::max(0.0, 0.5*(rp+r)-R))*ds;
    if (r <= R) return acc;
  }
  return acc;
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  struct C { double h, mu, f, lam; };
  C const cases[] = {
    {20000.0, 0.900, 1.0000, 400.0}, {20000.0, 0.500, 1.0000, 400.0},
    {20000.0, 0.100, 1.0000, 400.0}, {20000.0, 0.010, 1.0000, 400.0},
    {45000.0, 0.030, 1.0347, 690.0}, {2500.0,  0.700, 0.9713, 300.0},
  };
  std::printf("%8s %7s %8s %6s %14s %14s %10s %12s", "h", "mu", "f", "lam", "tau_by_b", "tau_traj", "|diff|", "b/bcrit"); puts("");
  double worst = 0.0;
  for (C const& c : cases) {
    double const rObs = R, rEmit = R + c.h;
    double const nObs = nPhase(rObs, c.f), nEmit = nPhase(rEmit, c.f);
    double const b = nObs*rObs*std::sqrt(std::max(0.0, 1.0 - c.mu*c.mu));
    double const s0 = b/(nEmit*rEmit);
    if (s0 >= 1.0) { std::printf("   %6.0f %7.3f  cannot launch", c.h, c.mu); puts(""); continue; }
    double const th0 = std::acos(-std::sqrt(1.0 - s0*s0));
    double const tb = tauByB(rEmit, b, c.f, c.lam, 400000);
    double const tt = tauByTraj(c.h, th0, c.f, c.lam, 0.05);
    double const d = std::abs(tb - tt);
    worst = std::max(worst, d);
    std::printf("%8.0f %7.3f %8.4f %6.0f %14.6f %14.6f %10.2e %12.6f",
                c.h, c.mu, c.f, c.lam, tb, tt, d, b/(nObs*rObs)); puts("");
  }
  std::printf("worst |tau_by_b - tau_traj| = %.3e   (a 1-D r-integral vs the RK4 trajectory)", worst); puts("");
  return 0;
}
