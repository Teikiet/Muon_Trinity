// Standalone diagnostic: curved-ray layer/basis-weight accumulation.
// Production routing is NOT touched; this mimics the build pass in a separate program.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0;
static AtmosphericAbsorption g_ab;
static int g_ncomp = 25;
static double g_hmax = 50000.0;

// corrected production callback: alpha = (tau(h+25)-tau(h-25))/50
static double alphaCB(double lam, double h) {
  double const dh = 50.0;
  auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h - dh * 0.5);
  double const hi = std::min(hr.second, h + dh * 0.5);
  double const a = hi - lo;
  if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam, hi) - g_ab.getOpticalDepth(lam, lo)) / a);
}
// alpha_j(lambda): thickness-weighted mean over component j, at the native lambda knots
static std::vector<double> g_lams;
static std::vector<std::vector<double>> g_alphaJ;  // [j][lambda index]
static void buildAlphaJ() {
  auto const wr = g_ab.getWavelengthRange();
  int const nl = 81;
  g_lams.resize(nl);
  for (int i = 0; i < nl; ++i) g_lams[i] = wr.first + (wr.second - wr.first) * double(i) / double(nl - 1);
  double const dcomp = g_hmax / g_ncomp;
  g_alphaJ.assign(g_ncomp, std::vector<double>(nl, 0.0));
  for (int j = 0; j < g_ncomp; ++j) {
    double const h0 = j * dcomp, h1 = (j + 1) * dcomp;
    int const ns = 10;                       // sub-samples inside the component
    for (int i = 0; i < nl; ++i) {
      double acc = 0.0;
      for (int k = 0; k < ns; ++k) acc += alphaCB(g_lams[i], h0 + (k + 0.5) * (h1 - h0) / ns);
      g_alphaJ[j][i] = acc / ns;
    }
  }
}
static double alphaJAt(int j, double lam) {
  int const nl = int(g_lams.size());
  if (lam <= g_lams[0]) return g_alphaJ[j][0];
  if (lam >= g_lams[nl-1]) return g_alphaJ[j][nl-1];
  int i = 0; while (i < nl-2 && g_lams[i+1] < lam) ++i;
  double const f = (lam - g_lams[i]) / (g_lams[i+1] - g_lams[i]);
  return g_alphaJ[j][i] * (1.0 - f) + g_alphaJ[j][i+1] * f;
}
// Bouguer ray trace with fine RK4, accumulating per-component path lengths
struct RayRes { bool hit{false}; double L{0.0}; std::vector<double> lj; };
static RayRes traceRay(double h0, double theta0, double fPhase, double ds) {
  RayRes o; o.lj.assign(g_ncomp, 0.0);
  double const rObs = R;
  double r = R + h0, th = theta0;
  if (r <= rObs) { o.hit = true; return o; }
  double const dcomp = g_hmax / g_ncomp;
  // exponential atmosphere: N(h) = N0 exp(-h/H), phase index n = 1 + N*fPhase
  double const N0 = 3.27e-4, H = 8000.0;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const h = rr - R;
    double const Nv = N0 * std::exp(-h / H) * fPhase;
    double const np = -Nv / H;                 // dN/dh
    double const nn = 1.0 + Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c;
    dt = -s * (1.0 / rr + (nn > 1e-12 ? np / nn : 0.0));
  };
  for (long i = 0; i < 200000000L; ++i) {
    double dr1,dt1,dr2,dt2,dr3,dt3,dr4,dt4;
    der(r, th, dr1, dt1);
    der(r + 0.5*ds*dr1, th + 0.5*ds*dt1, dr2, dt2);
    der(r + 0.5*ds*dr2, th + 0.5*ds*dt2, dr3, dt3);
    der(r + ds*dr3, th + ds*dt3, dr4, dt4);
    double const rPrev = r, thPrev = th;
    r += ds/6.0*(dr1+2*dr2+2*dr3+dr4);
    th += ds/6.0*(dt1+2*dt2+2*dt3+dt4);
    double const hMid = 0.5*(rPrev + r) - R;
    int const j = std::min(g_ncomp-1, std::max(0, int(hMid / dcomp)));
    o.lj[j] += ds; o.L += ds;
    if (r <= rObs) { o.hit = true; return o; }
  }
  return o;
}

// independent reference: integrate the corrected callback along the curved ray
static double tauReference(double h0, double theta0, double fPhase, double lam, double ds) {
  double r = R + h0, th = theta0;
  double const N0 = 3.27e-4, H = 8000.0;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const h = rr - R;
    double const Nv = N0 * std::exp(-h / H) * fPhase;
    double const np = -Nv / H, nn = 1.0 + Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c; dt = -s * (1.0 / rr + (nn > 1e-12 ? np / nn : 0.0));
  };
  double acc = 0.0;
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
    if (hMid >= 0.0) acc += alphaCB(lam, hMid) * ds;
    if (r <= R) return acc;
  }
  return acc;
}
int main(int argc, char** argv) {
  std::string prof = argv[1];
  g_ab.load(prof);
  if (!g_ab.isLoaded()) { std::printf("profile load failed"); return 1; }
  buildAlphaJ();
  double const dcomp = g_hmax / g_ncomp;
  std::printf("components=%d  alphaJ knots=%zu  alphaJ bytes=%.1f kB", g_ncomp,
              g_lams.size(), double(g_ncomp * g_lams.size() * 8) / 1024.0);

  // ---- constant-extinction analytic check: tau = alpha * L -----------------
  {
    double const alpha = 1.0e-4;
    for (double const ct : {1.0, 0.5, 0.2}) {
      double const th0 = std::acos(ct);
      RayRes const rr = traceRay(20000.0, th0, 1.0, 0.5);
      double tauSum = 0.0;
      for (int j = 0; j < g_ncomp; ++j) tauSum += alpha * rr.lj[j];
      std::printf("  constant-extinction: cos=%.2f  sum(alpha*l_j)=%.6f  alpha*L=%.6f  |d|=%.2e",
                  ct, tauSum, alpha * rr.L, std::abs(tauSum - alpha * rr.L));
    }
  }

  // ---- held-out validation: curved rays, off-grid heights, grazing, f range --
  // FROZEN dataset - never used for tuning.
  double const heldH[] = {2500.0, 7500.0, 15000.0, 25000.0, 45000.0};
  double const heldMu[] = {0.999, 0.7, 0.3, 0.12, 0.05, 0.01};
  double const heldF[] = {0.9713, 1.0, 1.0347};
  double const heldLam[] = {300.0, 302.5, 400.0, 402.5, 500.0, 502.5, 690.0, 700.0};
  double worst = 0.0, worstTau = 0.0;
  char info[160] = {0};
  int nray = 0, nmiss = 0, nCoarseRef = 0;
  for (double h0 : heldH)
    for (double mu : heldMu)
      for (double f : heldF) {
        // emission zenith from the Bouguer invariant at the observation level
        double const N0 = 3.27e-4, H = 8000.0;
        double const nObs = 1.0 + N0 * f, nEmit = 1.0 + N0 * std::exp(-h0 / H) * f;
        double const p = nObs * R * std::sqrt(std::max(0.0, 1.0 - mu * mu));
        double const s0 = p / (nEmit * (R + h0));
        if (s0 >= 1.0) { ++nmiss; continue; }
        double const th0 = std::acos(-std::sqrt(1.0 - s0 * s0));
        RayRes const rr = traceRay(h0, th0, f, 0.5);
        if (!rr.hit) { ++nmiss; continue; }
        for (double lam : heldLam) {
          double tauL = 0.0;
          for (int j = 0; j < g_ncomp; ++j) tauL += alphaJAt(j, lam) * rr.lj[j];
          // Reference quadrature: 0.05 m for short paths, bounded to ~1e5 steps for
          // long grazing paths. Rays where this is NOT finer than the 0.5 m layer
          // accumulation are counted and reported, not hidden.
          double const dsRef = std::max(0.05, rr.L / 1.0e5);
          if (dsRef > 0.5) ++nCoarseRef;
          double const tauR = tauReference(h0, th0, f, lam, dsRef);
          double const dT = std::abs(std::exp(-tauL) - std::exp(-tauR));
          ++nray;
          if (dT > worst) {
            worst = dT; worstTau = tauR;
            std::snprintf(info, sizeof(info), "h=%.0f mu=%.3f f=%.4f lam=%.1f tau=%.4f", h0, mu, f, lam, tauR);
          }
        }
      }
  std::printf("held-out rays=%d (missed=%d, coarse-ref=%d)  WORST |dT| = %.4e  [%s] (budget 1e-3)",
              nray, nmiss, nCoarseRef, worst, info);
  std::printf("worst-case tau there = %.4f (transmission %.4f)", worstTau, std::exp(-worstTau));

  // ---- cost: l_j table memory, build time, per-photon throughput ------------
  {
    int const nH = 976, nMu = 200;   // 4x coarsened (h, mu) grid
    double const bytes = double(nH) * nMu * g_ncomp * 4.0;
    std::printf("l_j table: %d x %d x %d floats = %.1f MB", nH, nMu, g_ncomp, bytes / 1048576.0);
    // per-photon evaluation cost (25 lambda-interpolations + 25 multiply-adds)
    std::vector<double> lj(g_ncomp, 100.0);
    long const n = 200000;
    auto const t0 = std::chrono::steady_clock::now();
    double acc = 0.0;
    for (long k = 0; k < n; ++k) {
      double const lam = 300.0 + 400.0 * double(k % 997) / 997.0;
      double t = 0.0;
      for (int j = 0; j < g_ncomp; ++j) t += alphaJAt(j, lam) * lj[j];
      acc += t;
    }
    auto const t1 = std::chrono::steady_clock::now();
    double const us = std::chrono::duration<double, std::micro>(t1 - t0).count() / double(n);
    std::printf("per-photon layer evaluation: %.4f us  (acc=%.6e)", us, acc);
  }
  // peak RSS
  {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
      if (line.rfind("VmHWM", 0) == 0) std::printf("peak RSS: %s", line.c_str() + 6);
  }
  return 0;
}
