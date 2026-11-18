// Section 2-5 diagnostics: guard derivation on DEVELOPMENT probes, then held-out
// validation with interior f. Production untouched.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
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
static double nP(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
// b_crit per configuration: min of n(r) r over the ground-bound route (protocol section 2)
static double bCritOf(double f) { double m = 1e300; for (double r = R; r <= R+50000.0; r += 25.0) m = std::min(m, nP(r,f)*r); return m; }
static double Aint(double b, double f, double lam, double rTop, double dr) {
  double acc = 0.0, r = R;
  while (r < rTop - 1e-12) {
    double const st = std::min(dr, rTop - r);
    double const rm = r + 0.5*st;
    double const s = b/(nP(rm,f)*rm);
    if (s >= 1.0) return acc;
    acc += alphaCB(lam, rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st;
    r += st;
  }
  return acc;
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  double const f = 1.0, bCrit = bCritOf(f);
  // ---- (1) reference convergence on difficult DEVELOPMENT cases ----------------
  std::printf("(1) reference convergence (development probes, incl. near-critical):"); puts("");
  double worstRef = 0.0;
  struct Dev { double h, ratio, lam; };
  Dev const devs[] = {{20000.0,0.9,400.0},{20000.0,0.999,400.0},{20000.0,0.9999,400.0},
                      {20000.0,0.99999,400.0},{20000.0,0.999999,400.0},{5000.0,0.99999,300.0},
                      {40000.0,0.9999,690.0},{1000.0,0.999999,500.0}};
  for (Dev const& d : devs) {
    double const b = d.ratio*bCrit;
    double const a1 = Aint(b,f,d.lam,R+d.h,2.0), a2 = Aint(b,f,d.lam,R+d.h,0.5), a3 = Aint(b,f,d.lam,R+d.h,0.05);
    double const dT = std::abs(std::exp(-a2)-std::exp(-a3));
    worstRef = std::max(worstRef, dT);
    std::printf("    h=%6.0f ratio=%.6f lam=%3.0f  tau=%.6f  |dT(0.5 vs 0.05)|=%.2e", d.h, d.ratio, d.lam, a3, dT); puts("");
  }
  std::printf("    worst reference uncertainty = %.2e   (budget 1e-5)", worstRef); puts("");
  // ---- (2) floating-point sensitivity on development rays ---------------------
  std::printf("(2) float sensitivity (rel 1e-13 on b, h, lambda):"); puts("");
  double worstSens = 0.0;
  for (Dev const& d : devs) {
    double const b = d.ratio*bCrit, h = d.h, lam = d.lam;
    double const t0 = Aint(b,f,lam,R+h,0.05);
    double const t1 = Aint(b*(1+1e-13), f, lam, R+h, 0.05);
    double const t2 = Aint(b, f, lam*(1+1e-13), R+h*(1+1e-13), 0.05);
    double const s = std::max(std::abs(std::exp(-t0)-std::exp(-t1)), std::abs(std::exp(-t0)-std::exp(-t2)));
    worstSens = std::max(worstSens, s);
  }
  std::printf("    worst |dT| under perturbation = %.2e   (budget 1e-5)", worstSens); puts("");
  // ---- (3) reachability / root conditioning -----------------------------------
  std::printf("(3) root conditioning at the endpoint (d(n r)/dr at r_obs):"); puts("");
  {
    double const dn = (nP(R+0.5,f)*(R+0.5) - nP(R,f)*R)/0.5;
    std::printf("    d(nr)/dr at r_obs = %.6e  (positive -> the endpoint root is a simple,", dn); puts("");
    std::printf("      integrable sqrt singularity; tau at b=b_crit = %.6f is FINITE)", Aint(bCrit,f,400.0,R+20000.0,0.01)); puts("");
  }
  // ---- (4) guard: acceptance ceiling = the table's top b node ------------------
  double const muTop = 1.0 - std::pow(std::sqrt(0.5)/216.0, 2.0);
  double const bTop = muTop*bCrit;
  std::printf("(4) guard: table b range [0, %.10f] = [0, %.8f b_crit]; b > bTop -> REJECT", bTop, muTop); puts("");
  std::printf("    the equality probe b = b_crit is %.2e above the top node -> rejected (no clamping)", bCrit-bTop); puts("");
  // ---- (5) held-out evaluation, INCLUDING INTERIOR f --------------------------
  std::printf("(5) held-out (frozen set, interior f included):"); puts("");
  double const heldH[] = {1500.0, 6000.0, 18000.0, 32000.0, 46000.0};
  double const heldR[] = {0.15, 0.35, 0.6, 0.8, 0.92, 0.97, 0.995, 0.9995};
  double const heldL[] = {310.0,355.0,400.0,445.0,490.0,535.0,580.0,625.0,670.0,695.0};
  double const heldM[] = {312.5,357.5,447.5,537.5,627.5,697.5};
  double const heldF[] = {0.9713, 0.99, 1.0, 1.01, 1.0347};
  int nAcc = 0, nRej = 0; double worstHeld = 0.0; char info[180] = {0};
  for (double h : heldH) for (double rr : heldR) for (double fv : heldF) {
    double const bc = bCritOf(fv);
    if (rr*bc > muTop*bc) { ++nRej; continue; }              // above the table's b range
    if (h > 50000.0 || h < 0.0) { ++nRej; continue; }
    for (double lam : heldL) for (double lm : {0.0, 1.0}) {
      double const L = (lm == 0.0) ? lam : lam + 2.5;         // knot or adjacent-knot midpoint
      if (L < 300.0 || L > 700.0) { ++nRej; continue; }
      double const b = rr*bc;
      if (b > nP(R+h,fv)*(R+h)) { ++nRej; continue; }         // launchability (independent)
      double const tR = Aint(b, fv, L, R+h, 0.5);             // converged reference (dr=0.5 justified by check 1)
      double const tA = Aint(b, f, L, R+h, 5.0);              // table-equivalent: f=1 construction
      double const dT = std::abs(std::exp(-tA)-std::exp(-tR));
      ++nAcc;
      if (dT > worstHeld) { worstHeld = dT; std::snprintf(info, sizeof(info), "h=%.0f ratio=%.4f f=%.4f lam=%.1f tau=%.4f", h, rr, fv, L, tR); }
    }
  }
  std::printf("    accepted %d, rejected %d", nAcc, nRej); puts("");
  std::printf("    WORST combined |dT| = %.4e  [%s]", worstHeld, info); puts("");
  std::printf("    acceptance: every eligible held-out ray <= 1e-3 ? %s", (worstHeld <= 1e-3) ? "PASS" : "FAIL"); puts("");
  return 0;
}
