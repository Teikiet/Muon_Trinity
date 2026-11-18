// Diagnostic-only coarsening/compression study. Production routing untouched.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
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
static std::vector<double> g_hk;
static double alphaCB(double lam, double h) {
  double const dh = 50.0; auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h-dh*0.5), hi = std::min(hr.second, h+dh*0.5);
  double const a = hi-lo; if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam,hi) - g_ab.getOpticalDepth(lam,lo))/a);
}
static double slopeK(int k, double lam) {
  return (g_ab.getOpticalDepth(lam,g_hk[k+1]) - g_ab.getOpticalDepth(lam,g_hk[k]))/(g_hk[k+1]-g_hk[k]);
}
static double phiK(int k, double h) {
  double const a = g_hk[k], b = g_hk[k+1], w = 25.0;
  bool const first = (k==0), last = (k+2==int(g_hk.size()));
  if (h <= a-w || h >= b+w) return 0.0;
  if (h < a+w) return first ? 1.0 : (h-(a-w))/(2.0*w);
  if (h > b-w) return last ? 1.0 : ((b+w)-h)/(2.0*w);
  return 1.0;
}
static double thetaFor(double h0, double mu, double f) {
  double const nObs = 1.0 + N0*f, nEmit = 1.0 + N0*std::exp(-h0/HSC)*f;
  double const p = nObs*R*std::sqrt(std::max(0.0,1.0-mu*mu));
  double const s0 = p/(nEmit*(R+h0)); if (s0 >= 1.0) return -1.0;
  return std::acos(-std::sqrt(1.0-s0*s0));
}
// accumulate W_k (native 50) and the direct tau, at step ds
static void traceW(double h0, double th0, double f, double lam, double ds, std::vector<double>& W, double& tau, double& L) {
  W.assign(50, 0.0); tau = 0.0; L = 0.0;
  double r = R + h0, th = th0;
  if (r <= R) return;
  auto der = [&](double rr, double tt, double& dr, double& dt) {
    double const h = rr-R, Nv = N0*std::exp(-h/HSC)*f, np = -Nv/HSC, nn = 1.0+Nv;
    double const s = std::sin(tt), c = std::cos(tt);
    dr = c; dt = -s*(1.0/rr + (nn>1e-12?np/nn:0.0));
  };
  for (long i = 0; i < 400000000L; ++i) {
    double a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,a4,b4,c4,d4;
    der(r,th,a1,b1); der(r+0.5*ds*a1, th+0.5*ds*b1, a2,b2);
    der(r+0.5*ds*a2, th+0.5*ds*b2, a3,b3); der(r+ds*a3, th+ds*b3, a4,b4);
    double const rp = r;
    r += ds/6.0*(a1+2*a2+2*a3+a4); th += ds/6.0*(b1+2*b2+2*b3+b4);
    double const hm = std::max(0.0, 0.5*(rp+r)-R);
    for (int k = 0; k < 50; ++k) W[k] += phiK(k, hm)*ds;
    tau += alphaCB(lam, hm)*ds; L += ds;
    if (r <= R) return;
  }
}
static double tauDirectConv(double h0, double th0, double f, double lam) {
  double W[1]; double tau = 0.0, L = 0.0;
  std::vector<double> v; traceW(h0, th0, f, lam, 1.0, v, tau, L);
  double const dsc = std::max(0.05, L/2.0e6);
  std::vector<double> v2; traceW(h0, th0, f, lam, dsc, v2, tau, L);
  return tau;
}

// trilinear interpolation of the stored W grid; returns merged-component weights for a
// candidate merge factor mm (50/mm components, s_G = unweighted mean of s_k over the group)
static void interpW(std::vector<float> const& Wg, std::vector<double> const& gH,
                    std::vector<double> const& gMu, double h, double mu, double f,
                    int mm, std::vector<double>& Wout) {
  int const nF = 3; double const gF[3] = {0.97, 1.0, 1.03};
  int i0 = 0; while (i0+2 < int(gH.size()) && gH[i0+1] < h) ++i0;
  int m0 = 0; while (m0+2 < int(gMu.size()) && gMu[m0+1] > mu) ++m0;   // mu DESCENDING
  int f0 = 0; if (f > gF[1]) f0 = 1;
  double th = (gH[i0+1] > gH[i0]) ? (h-gH[i0])/(gH[i0+1]-gH[i0]) : 0.0;
  double tm = (gMu[m0] > gMu[m0+1]) ? (gMu[m0]-mu)/(gMu[m0]-gMu[m0+1]) : 0.0;
  double tf = (gF[f0+1] > gF[f0]) ? (f-gF[f0])/(gF[f0+1]-gF[f0]) : 0.0;
  th = std::min(1.0, std::max(0.0, th)); tm = std::min(1.0, std::max(0.0, tm));
  int const nMu = int(gMu.size());
  int const nM = 50/mm;
  Wout.assign(nM, 0.0);
  for (int dh = 0; dh < 2; ++dh) for (int dm = 0; dm < 2; ++dm) for (int df = 0; df < 2; ++df) {
    double const w = (dh? th : 1.0-th) * (dm? tm : 1.0-tm) * (df? tf : 1.0-tf);
    if (w == 0.0) continue;
    std::size_t const base = ((( std::size_t(i0+dh))*nMu + (m0+dm))*nF + (f0+df))*50;
    for (int g = 0; g < nM; ++g) { double acc = 0.0;
      for (int k = g*mm; k < (g+1)*mm; ++k) acc += double(Wg[base+k]);
      Wout[g] += w*acc; }
  }
}
// merged s_G(lambda): unweighted mean of the native slopes over the group (lambda-only)
static double sGroup(int g, int mm, double lam) {
  double acc = 0.0;
  for (int k = g*mm; k < (g+1)*mm; ++k) acc += slopeK(k, lam);
  return acc/mm;
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  auto const hr = g_ab.getHeightRange();
  for (double h = hr.first; h <= hr.second+1e-9; h += 1000.0) g_hk.push_back(h);
  // ---- build the W storage grid (frozen): 21 h x 11 mu x 3 f ----------------
  std::vector<double> gH, gMu; double const gF[3] = {0.97, 1.0, 1.03};
  for (int i = 0; i < 21; ++i) gH.push_back(2500.0*i);
  for (double m : {0.999,0.9,0.7,0.5,0.35,0.2,0.12,0.06,0.03,0.015,0.005}) gMu.push_back(m);
  std::size_t const nG = gH.size()*gMu.size()*3;
  std::printf("W grid: %zu h x %zu mu x 3 f = %zu rays", gH.size(), gMu.size(), nG); puts("");
  auto const t0 = std::chrono::steady_clock::now();
  std::vector<float> Wgrid(nG*50);
  for (std::size_t ih = 0; ih < gH.size(); ++ih) for (std::size_t im = 0; im < gMu.size(); ++im) for (int iF = 0; iF < 3; ++iF) {
    double const th0 = thetaFor(gH[ih], gMu[im], gF[iF]);
    std::size_t const base = ((ih*gMu.size()+im)*3+iF)*50;
    if (th0 < 0.0) { for (int k = 0; k < 50; ++k) Wgrid[base+k] = 0.f; continue; }
    std::vector<double> W; double tau=0.0, L=0.0;
    traceW(gH[ih], th0, gF[iF], 400.0, 2.0, W, tau, L);
    for (int k = 0; k < 50; ++k) Wgrid[base+k] = float(W[k]);
  }
  auto const t1 = std::chrono::steady_clock::now();
  std::printf("build time: %.1f s   W grid bytes: %.1f MB", std::chrono::duration<double>(t1-t0).count(), double(Wgrid.size())*4/1048576.0); puts("");
  { std::ifstream f("/proc/self/status"); std::string l;
    while (std::getline(f,l)) if (l.rfind("VmHWM",0)==0) { std::printf("peak after W grid: %s", l.c_str()); puts(""); } }
  // ---- evaluation: isolated and combined errors on dev then held-out ----------
  struct S { double h, mu, f, lam; };
  std::vector<S> dev, held;
  for (double h : {4000.0,10000.0,30000.0}) for (double m : {0.8,0.4,0.2,0.08}) for (double f : {0.99,1.01}) for (double l : {310.0,430.0,600.0}) dev.push_back({h,m,f,l});
  for (double h : {3500.0,12000.0,35000.0}) for (double m : {0.85,0.45,0.15,0.03}) for (double f : {0.99,1.02}) for (double l : {305.0,420.0,610.0,512.5}) held.push_back({h,m,f,l});
  for (int pass = 0; pass < 2; ++pass) {
    std::vector<S> const& set = pass ? held : dev;
    std::printf("%s set: %zu rays", pass ? "HELD-OUT" : "dev", set.size()); puts("");
    for (int mm : {1, 2, 5, 10}) {   // merge factor; M = 50/mm components
      double worstMerge = 0.0, worstComb = 0.0;
      for (S const& q : set) {
        double const th0 = thetaFor(q.h, q.mu, q.f);
        if (th0 < 0.0) continue;
        std::vector<double> W; double tau=0.0, L=0.0; traceW(q.h, th0, q.f, q.lam, 2.0, W, tau, L);
        // isolated: merge only (native weights, unmerged tau vs merged tau)
        double tNative = 0.0; for (int k = 0; k < 50; ++k) tNative += slopeK(k, q.lam)*W[k];
        double tMerged = 0.0; for (int g = 0; g < 50/mm; ++g) { double wg = 0.0; for (int k=g*mm;k<(g+1)*mm;++k) wg += W[k]; tMerged += sGroup(g, mm, q.lam)*wg; }
        worstMerge = std::max(worstMerge, std::abs(std::exp(-tNative)-std::exp(-tMerged)));
        // combined: merged + trilinear interpolation + float32 storage
        std::vector<double> Wi; interpW(Wgrid, gH, gMu, q.h, q.mu, q.f, mm, Wi);
        double tComb = 0.0; for (int g = 0; g < 50/mm; ++g) tComb += sGroup(g, mm, q.lam)*Wi[g];
        double const tRef = tauDirectConv(q.h, th0, q.f, q.lam);
        worstComb = std::max(worstComb, std::abs(std::exp(-tComb)-std::exp(-tRef)));
      }
      std::printf("   M=%-3d merge-only |dT|=%.3e   COMBINED |dT|=%.3e", 50/mm, worstMerge, worstComb); puts("");
    }
  }
  // ---- lookup timing (stored weights + interpolation + dot product) ----------
  {
    std::vector<double> Wi; long const n = 200000;
    auto const a0 = std::chrono::steady_clock::now(); double acc = 0.0;
    for (long k = 0; k < n; ++k) {
      double const mu = 0.005 + 0.99*double(k%997)/997.0;
      interpW(Wgrid, gH, gMu, 20000.0, mu, 1.0, 10, Wi);
      for (int g = 0; g < 10; ++g) acc += Wi[g];
    }
    auto const a1 = std::chrono::steady_clock::now();
    std::printf("lookup (M=10, incl. interpolation): %.3f us/photon (acc=%.3e)",
                std::chrono::duration<double,std::micro>(a1-a0).count()/double(n), acc); puts("");
  }
  return 0;
}
