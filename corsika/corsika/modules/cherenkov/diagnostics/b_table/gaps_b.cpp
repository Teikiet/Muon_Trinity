// Gap-fixing diagnostics for the cumulative-b design. Production untouched.
//  (1) f_phase dependence: table built at f=1 vs a converged reference at the TRUE f.
//  (2) explicit fallback exercise with reasons and counts.
//  (3) combined diagnostic route timing: geometry lookup + absorption lookup + dispatch.
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
static double nPhase(double r, double f) { return 1.0 + N0*std::exp(-(r-R)/HSC)*f; }
static double Aintegral(double b, double f, double lam, double rTop, double dr) {
  double acc = 0.0, r = R;
  while (r < rTop - 1e-12) {
    double const step = std::min(dr, rTop - r);
    double const rm = r + 0.5*step;
    double const s = b/(nPhase(rm,f)*rm);
    if (s >= 1.0) return acc;
    acc += alphaCB(lam, rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*step;
    r += step;
  }
  return acc;
}
int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  auto const hr = g_ab.getHeightRange();
  for (double h = hr.first; h <= hr.second+1e-9; h += 1000.0) g_hk.push_back(h);
  // frozen axes (unchanged)
  std::vector<double> gH;
  for (double h : {0.0,25.0,50.0,100.0,250.0,500.0}) gH.push_back(h);
  for (int i = 1; i <= 50; ++i) gH.push_back(1000.0*i);
  std::vector<double> gMu;
  for (int i = 0; i < 40; ++i) gMu.push_back(2.0e-4*std::pow(0.5/2.0e-4, double(i)/39.0));
  for (int i = 1; i <= 216; ++i) { double const u = std::sqrt(0.5)*double(i)/216.0; gMu.push_back(1.0-u*u); }
  std::sort(gMu.begin(), gMu.end());
  std::vector<double> gL; for (double w = 300.0; w <= 700.0+1e-9; w += 5.0) gL.push_back(w);
  int const nH = int(gH.size()), nMu = int(gMu.size()), nL = int(gL.size());
  double const fBuild = 1.0, bCrit = nPhase(R,fBuild)*R;
  std::printf("axes: %d h x %d mu x %d lam = %.2f MB   (built at f=1)", nH, nMu, nL,
              double(nH)*nMu*nL*4/1048576.0); puts("");
  std::vector<float> A(std::size_t(nH)*nMu*nL);
  for (int im = 0; im < nMu; ++im) for (int il = 0; il < nL; ++il) {
    double const b = gMu[im]*bCrit; double acc = 0.0, r = R; int ih = 0;
    while (ih < nH) {
      double const rTop = R + gH[ih];
      while (r < rTop - 1e-12) {
        double const step = std::min(5.0, rTop - r);
        double const rm = r + 0.5*step;
        double const s = b/(nPhase(rm,fBuild)*rm);
        if (s >= 1.0) break;
        acc += alphaCB(gL[il], rm-R)/std::sqrt(std::max(1e-300,1.0-s*s))*step;
        r += step;
      }
      A[std::size_t(ih)*nMu*nL + im*nL + il] = float(acc); ++ih;
    }
  }
  auto lookup = [&](double h, double mu, double lam) -> double {
    int i0=0; while (i0+2<nH && gH[i0+1]<h) ++i0;
    int m0=0; while (m0+2<nMu && gMu[m0+1]<mu) ++m0;
    int l0=0; while (l0+2<nL && gL[l0+1]<lam) ++l0;
    double const th=(gH[i0+1]>gH[i0])?std::min(1.0,std::max(0.0,(h-gH[i0])/(gH[i0+1]-gH[i0]))):0.0;
    double const tm=(gMu[m0+1]>gMu[m0])?std::min(1.0,std::max(0.0,(mu-gMu[m0])/(gMu[m0+1]-gMu[m0]))):0.0;
    double const tl=(gL[l0+1]>gL[l0])?std::min(1.0,std::max(0.0,(lam-gL[l0])/(gL[l0+1]-gL[l0]))):0.0;
    double acc=0.0;
    for (int a=0;a<2;++a) for (int q=0;q<2;++q) for (int c=0;c<2;++c) {
      double const w=(a?th:1-th)*(q?tm:1-tm)*(c?tl:1-tl);
      if (w==0.0) continue;
      acc += w*double(A[std::size_t(i0+a)*nMu*nL + (q? m0+1:m0)*nL + (c? l0+1:l0)]);
    }
    return acc;
  };
  // ---- (1) f_phase dependence ---------------------------------------------
  {
    std::printf("(1) f_phase dependence (table built at f=1; reference uses the true f):"); puts("");
    for (double f : {0.9713, 1.0, 1.0347}) {
      double worst = 0.0; char info[140] = {0};
      for (double h : {5000.0, 20000.0, 40000.0})
        for (double mu : {0.9, 0.5, 0.15, 0.04})
          for (double lam : {305.0, 400.0, 500.0, 690.0}) {
            double const bTrue = mu*nPhase(R,f)*R;
            if (bTrue >= nPhase(R+h,f)*(R+h)) continue;
            double const tR = Aintegral(bTrue, f, lam, R+h, 0.5);
            double const tA = lookup(h, mu, lam);
            double const dT = std::abs(std::exp(-tA) - std::exp(-tR));
            if (dT > worst) { worst = dT; std::snprintf(info, sizeof(info), "h=%.0f mu=%.2f lam=%.0f tau=%.4f", h, mu, lam, tR); }
          }
      std::printf("    f=%.4f : WORST |dT| = %.4e   [%s]", f, worst, info); puts("");
    }
  }
  // ---- (2) explicit fallback exercise -------------------------------------
  {
    std::printf("(2) fallback exercise:"); puts("");
    int nBeyond = 0, nUpgoing = 0, nAboveTop = 0, nBelowObs = 0, nLamOut = 0, nNoTable = 0, nOk = 0;
    for (double h : {0.0, 1000.0, 20000.0, 50000.0, 60000.0})
      for (double mu : {0.9, 0.1, 0.005, 0.0, -0.5})
        for (double lam : {290.0, 300.0, 700.0, 710.0}) {
          if (h > 50000.0) { ++nAboveTop; continue; }
          if (h < 0.0) { ++nBelowObs; continue; }
          if (lam < 300.0 || lam > 700.0) { ++nLamOut; continue; }
          double const f = 1.0;
          if (mu < 0.0) { ++nUpgoing; continue; }                       // emission outward
          double const b = mu*nPhase(R,f)*R;
          if (b > bCrit) { ++nBeyond; continue; }
          // reaches the observation sphere iff b <= n(r)r for all r in [R, R+h] (n*r monotone)
          if (b > nPhase(R+h,f)*(R+h)) { ++nBeyond; continue; }
          ++nOk;
        }
    std::printf("    beyond-critical %d | upgoing %d | above-top %d | below-obs %d | lambda-out %d | no-table %d | accepted %d",
                nBeyond, nUpgoing, nAboveTop, nBelowObs, nLamOut, nNoTable, nOk); puts("");
  }
  // ---- (3) combined diagnostic route timing --------------------------------
  {
    std::vector<double> alt, lnN;
    for (double h = 0.0; h <= 120.0+1e-9; h += 0.01) { alt.push_back(h); lnN.push_back(std::log(N0*std::exp(-h*1000.0/HSC))); }
    corsika::cherenkov::RefractiveIndexTable ri(alt, lnN, R);
    corsika::cherenkov::CurvedRayTable geom(ri, R, R, 120000.0, 800, 10.0, 100.0, 30000.0);
    long const n = 200000;
    auto const t0 = std::chrono::steady_clock::now();
    double acc = 0.0; int used = 0, fell = 0;
    for (long k = 0; k < n; ++k) {
      double const h = 1000.0 + 30000.0*double(k%991)/991.0;
      double const mu = 0.02 + 0.9*double(k%997)/997.0;
      double const lam = 300.0 + 400.0*double(k%983)/983.0;
      // dispatch: eligibility, then geometry lookup, then absorption lookup
      bool eligible = (h <= 50000.0) && (lam >= 300.0 && lam <= 700.0) && (mu > 0.0);
      if (!eligible) { ++fell; continue; }
      double const r0 = R + h;
      double const s0 = mu*nPhase(R,1.0)*R/(nPhase(r0,1.0)*r0);
      if (s0 >= 1.0) { ++fell; continue; }
      double const th0 = std::acos(-std::sqrt(1.0 - s0*s0));
      auto const rs = geom.trace(0.0, 0.0, r0, std::sin(th0), 0.0, std::cos(th0), 1.0, 1.0);
      if (!rs.hit) { ++fell; continue; }
      acc += lookup(h, mu, lam);
      ++used;
    }
    auto const t1 = std::chrono::steady_clock::now();
    std::printf("(3) COMBINED route (eligibility + geometry lookup + absorption lookup): %.4f us/photon",
                std::chrono::duration<double,std::micro>(t1-t0).count()/double(n)); puts("");
    std::printf("    accepted %d, fell back %d  (acc=%.6e)", used, fell, acc); puts("");
    std::ifstream st("/proc/self/status"); std::string l;
    while (std::getline(st,l)) if (l.rfind("VmHWM",0)==0) { std::printf("    TOTAL peak RSS: %s", l.c_str()); puts(""); }
  }
  return 0;
}
