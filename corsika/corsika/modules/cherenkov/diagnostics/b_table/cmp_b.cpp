// PHASE 3: bounded architecture comparison with ACTUAL lookup outputs.
// Candidates: (A) cumulative-b tau table, (B) direct tau(h,mu,lambda) table.
// Lambda axis spans 300-900 nm reproducing the CURRENT clamping semantics.
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
using corsika::cherenkov::AtmosphericAbsorption;
static constexpr double R = 6371000.0, N0 = 3.27e-4, HSC = 8000.0;
static AtmosphericAbsorption g_ab;
static double alphaCB(double lam, double h) {   // reproduces production, incl. clamping
  double const dh = 50.0; auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h-dh*0.5), hi = std::min(hr.second, h+dh*0.5);
  double const a = hi-lo; if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam,hi) - g_ab.getOpticalDepth(lam,lo))/a);
}
static double nP(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
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
  double const f = 1.0, bCrit = nP(R,f)*R;
  // ---- axes (frozen): h 1 km knots + near-ground; b graded; lambda 300..900 clamped
  std::vector<double> gH; for (double h : {0.0,25.0,50.0,100.0,250.0,500.0}) gH.push_back(h);
  for (int i = 1; i <= 50; ++i) gH.push_back(1000.0*i);
  std::vector<double> gMu;
  for (int i = 0; i < 40; ++i) gMu.push_back(2.0e-4*std::pow(0.5/2.0e-4, double(i)/39.0));
  for (int i = 1; i <= 216; ++i) { double const u = std::sqrt(0.5)*double(i)/216.0; gMu.push_back(1.0-u*u); }
  std::sort(gMu.begin(), gMu.end());
  std::vector<double> gL; for (double w = 300.0; w <= 900.0+1e-9; w += 5.0) gL.push_back(w);   // incl. clamped region
  int const nH=int(gH.size()), nMu=int(gMu.size()), nL=int(gL.size());
  std::printf("CANDIDATE A (cumulative-b): axes %d h x %d mu x %d lambda = %d nodes = %.2f MB",
              nH, nMu, nL, nH*nMu*nL, double(nH)*nMu*nL*4/1048576.0); puts("");
  auto const t0 = std::chrono::steady_clock::now();
  std::vector<float> A(std::size_t(nH)*nMu*nL);
  for (int im=0; im<nMu; ++im) { double const b=gMu[im]*bCrit;
    for (int il=0; il<nL; ++il) { double acc=0.0, r=R; int ih=0;
      while (ih<nH) { double const rTop=R+gH[ih];
        while (r < rTop-1e-12) { double const st=std::min(5.0, rTop-r); double const rm=r+0.5*st;
          double const s=b/(nP(rm,f)*rm); if (s>=1.0) break;
          acc += alphaCB(gL[il], rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*st; r+=st; }
        A[std::size_t(ih)*nMu*nL + im*nL + il]=float(acc); ++ih; } } }
  auto const t1 = std::chrono::steady_clock::now();
  std::printf("  build %.1f s   table bytes %.2f MB", std::chrono::duration<double>(t1-t0).count(), double(A.size())*4/1048576.0); puts("");
  auto lookupA = [&](double h, double mu, double lam) -> double {
    int i0=0; while (i0+2<nH && gH[i0+1]<h) ++i0;
    int m0=0; while (m0+2<nMu && gMu[m0+1]<mu) ++m0;
    int l0=0; while (l0+2<nL && gL[l0+1]<lam) ++l0;
    double const th=(gH[i0+1]>gH[i0])?std::min(1.0,std::max(0.0,(h-gH[i0])/(gH[i0+1]-gH[i0]))):0.0;
    double const tm=(gMu[m0+1]>gMu[m0])?std::min(1.0,std::max(0.0,(mu-gMu[m0])/(gMu[m0+1]-gMu[m0]))):0.0;
    double const tl=(gL[l0+1]>gL[l0])?std::min(1.0,std::max(0.0,(lam-gL[l0])/(gL[l0+1]-gL[l0]))):0.0;
    double acc=0.0;
    for (int a=0;a<2;++a) for (int q=0;q<2;++q) for (int c=0;c<2;++c) {
      double const w=(a?th:1-th)*(q?tm:1-tm)*(c?tl:1-tl); if (w==0.0) continue;
      acc += w*double(A[std::size_t(i0+a)*nMu*nL + (m0+q)*nL + (l0+c)]); }
    return acc; };
  // ---- REGRESSION set (the old held-out set; no longer fresh) ----------------
  struct Q { double h, mu, lam; };
  std::vector<Q> reg;
  for (double h : {1500.0, 6000.0, 18000.0, 32000.0, 46000.0})
    for (double r : {0.15,0.35,0.6,0.8,0.92,0.97,0.995,0.9995})
      for (double lam : {310.0,312.5,355.0,357.5,400.0,402.5,445.0,447.5,490.0,492.5,535.0,537.5,580.0,582.5,625.0,627.5,670.0,672.5,695.0,697.5,720.0,850.0}) reg.push_back({h,r,lam});
  double worst=0.0, worstTau=0.0; char info[170]={0}; int nAcc=0, nRej=0;
  for (Q const& q : reg) {
    double const b = q.mu*bCrit;
    if (b >= nP(R+q.h,f)*(R+q.h)) { ++nRej; continue; }
    double const tR = Aint(b, f, q.lam, R+q.h, 0.5);      // converged reference
    double const tA = lookupA(q.h, q.mu, q.lam);          // ACTUAL interpolated lookup
    double const dT = std::abs(std::exp(-tA)-std::exp(-tR));
    ++nAcc;
    if (dT>worst) { worst=dT; worstTau=tR; std::snprintf(info,sizeof(info),"h=%.0f ratio=%.4f lam=%.1f tau=%.4f",q.h,q.mu,q.lam,tR); }
  }
  std::printf("  ACTUAL lookup accuracy: %d rays (%d rejected)  WORST |dT| = %.4e  [%s]", nAcc, nRej, worst, info); puts("");
  std::printf("  gate 2 (accuracy <= 1e-3): %s", (worst<=1e-3)?"PASS":"FAIL"); puts("");
  // ---- CONFIGURATION SENSITIVITY: table built at f=1, photon has f != 1 ------
  std::printf("  CONFIGURATION SENSITIVITY (table built at f=1; reference at true f):"); puts("");
  for (double fTrue : {0.9713, 1.01, 1.0347}) {
    double worstC=0.0; char ic[170]={0};
    for (Q const& q : reg) {
      double const bc = nP(R,fTrue)*R;
      double const b = q.mu*bc;
      if (b >= nP(R+q.h,fTrue)*(R+q.h)) continue;
      double const tR = Aint(b, fTrue, q.lam, R+q.h, 0.5);
      double const tA = lookupA(q.h, q.mu, q.lam);
      double const dT = std::abs(std::exp(-tA)-std::exp(-tR));
      if (dT>worstC) { worstC=dT; std::snprintf(ic,sizeof(ic),"h=%.0f ratio=%.4f lam=%.1f",q.h,q.mu,q.lam); }
    }
    std::printf("    f=%.4f  WORST |dT| = %.4e  [%s]  %s", fTrue, worstC, ic, (worstC<=1e-3)?"PASS":"FAIL"); puts("");
  }
  // ---- lookup cost ---------------------------------------------------------
  long const n=500000; auto const u0=std::chrono::steady_clock::now(); double acc2=0.0;
  for (long k=0;k<n;++k) acc2 += lookupA(1000.0+30000.0*double(k%991)/991.0, 0.02+0.9*double(k%997)/997.0, 300.0+600.0*double(k%983)/983.0);
  auto const u1=std::chrono::steady_clock::now();
  std::printf("  lookup cost %.4f us/call (acc=%.3e)", std::chrono::duration<double,std::micro>(u1-u0).count()/double(n), acc2); puts("");
  { std::ifstream st("/proc/self/status"); std::string l;
    while (std::getline(st,l)) if (l.rfind("VmHWM",0)==0) { std::printf("  peak RSS (this program): %s", l.c_str()); puts(""); } }
  return 0;
}
