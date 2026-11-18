// reaxis_b.cpp -- diagnostics-only re-axis experiment.
//
// Compares axis variants for the cumulative-xi absorption table using MEASURED error
// and MEASURED cost, on a development ray set (never the frozen final-validation set):
//   - grazing axis: targeted refinement vs a different coordinate vs the current one
//   - height axis:  short-path grading (tau ~ sqrt(h) near the ground for grazing rays)
//   - lambda axis:  the profile's native knots over its DECLARED domain
// Reference is computed ONCE per ray and cached, so many variants share it.
// No production file is touched; no geometry code is modified.
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0;
static double Nshape(double h) {
  return 2.20e-4 * std::exp(-h / 7500.0) + 0.55e-4 * std::exp(-h / 1600.0);
}
static double Fdisp(double lam_nm) {
  double const l = lam_nm / 1000.0;
  return (1.0 + 0.0032 / (l * l)) / (1.0 + 0.0032 / 0.16);
}
static RefractiveIndexTable makeRI() {
  std::vector<double> alt, lnN;
  for (double h = 0.0; h <= H_TOP + 1e-9; h += 100.0) {
    alt.push_back(h / 1000.0);
    lnN.push_back(std::log(Nshape(h)));
  }
  return RefractiveIndexTable(alt, lnN, R);
}
static double nPhase(RefractiveIndexTable const& ri, double zAbs, double f) {
  return 1.0 + ri.N(zAbs) * f;
}
static double bCrit(RefractiveIndexTable const& ri, double f) { return nPhase(ri, R, f) * R; }

static AtmosphericAbsorption g_ab;
static double alphaOf(double lam, double h) {
  double const dh = 50.0;
  auto const hr = g_ab.getHeightRange();
  double const lo = std::max(hr.first, h - dh * 0.5), hi = std::min(hr.second, h + dh * 0.5);
  double const a = hi - lo;
  if (a <= 0.0) return 0.0;
  return std::max(0.0, (g_ab.getOpticalDepth(lam, hi) - g_ab.getOpticalDepth(lam, lo)) / a);
}
// Reference: w = sqrt(z - R) substitution (dz = 2w dw) removes the endpoint stiffness.
static double tauRefW(RefractiveIndexTable const& ri, double p, double f, double lam, double hTop,
                      long nW) {
  if (hTop <= 0.0) return 0.0;
  double const wTop = std::sqrt(hTop), dw = wTop / double(nW);
  double acc = 0.0;
  for (long k = 0; k < nW; ++k) {
    double const w = (double(k) + 0.5) * dw;
    double const z = R + w * w;
    double const s = p / (nPhase(ri, z, f) * z);
    if (s >= 1.0) return acc;
    acc += alphaOf(lam, w * w) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * (2.0 * w * dw);
  }
  return acc;
}
// ---- axis builders -----------------------------------------------------------------
struct Axes {
  std::vector<double> h, xi, lam;
  std::string name;
};
static std::vector<double> hCurrent() {
  std::vector<double> v{0, 25, 50, 100, 250, 500};
  for (int i = 1; i <= 50; ++i) v.push_back(1000.0 * i);
  return v;
}
static std::vector<double> hGraded() {  // geometric to 500 m, then 1 km  (61 nodes)
  std::vector<double> v{0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 500};
  for (int i = 1; i <= 50; ++i) v.push_back(1000.0 * i);
  return v;
}
static std::vector<double> hGraded2() {  // geometric to 2 km, then 2 km  (37 nodes)
  std::vector<double> v{0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2000};
  for (int i = 2; i <= 50; i += 2) v.push_back(1000.0 * i);
  return v;
}
static void xiLog(std::vector<double>& v, int n) {
  for (int i = 0; i < n; ++i) v.push_back(2.0e-4 * std::pow(0.5 / 2.0e-4, double(i) / double(n - 1)));
}
static std::vector<double> xiCurrent() {  // 48 log + 208 uniform in u = sqrt(1-xi)
  std::vector<double> v;
  xiLog(v, 48);
  for (int i = 1; i <= 208; ++i) {
    double const u = std::sqrt(0.5) * double(i) / 208.0;
    v.push_back(1.0 - u * u);
  }
  std::sort(v.begin(), v.end());
  return v;
}
static std::vector<double> xiRefined(int n) {  // targeted refinement of the same coordinate
  std::vector<double> v;
  xiLog(v, 48);
  for (int i = 1; i <= n; ++i) {
    double const u = std::sqrt(0.5) * double(i) / double(n);
    v.push_back(1.0 - u * u);
  }
  std::sort(v.begin(), v.end());
  return v;
}
static std::vector<double> xiLogOneMinus(int n) {  // uniform in y = -log10(1 - xi)
  std::vector<double> v;
  xiLog(v, 48);
  double const y0 = -std::log10(0.5), y1 = -std::log10(1.0e-9);
  for (int i = 0; i < n; ++i) {
    double const y = y0 + (y1 - y0) * double(i) / double(n - 1);
    v.push_back(1.0 - std::pow(10.0, -y));
  }
  std::sort(v.begin(), v.end());
  return v;
}
static std::vector<double> xiSinTheta(int n) {  // uniform in sin(theta_obs) = sqrt(1-xi^2)
  std::vector<double> v;
  xiLog(v, 48);
  double const s0 = std::sqrt(0.75);
  for (int i = 0; i < n; ++i) {
    double const s = s0 + (1.0 - s0) * double(i) / double(n - 1);
    v.push_back(std::sqrt(std::max(0.0, 1.0 - s * s)));
  }
  std::sort(v.begin(), v.end());
  return v;
}
static std::vector<double> lamNative() {  // the profile's own knots over its declared domain
  std::vector<double> v;
  for (double w = 180.0; w <= 700.0 + 1e-9; w += 5.0) v.push_back(w);
  return v;
}
static std::vector<double> lamWide() {  // the phase-3 axis (300-900, clamped above 700)
  std::vector<double> v;
  for (double w = 300.0; w <= 900.0 + 1e-9; w += 5.0) v.push_back(w);
  return v;
}

// ---- table -------------------------------------------------------------------------
struct Table {
  Axes ax;
  int nH = 0, nX = 0, nL = 0;
  std::vector<float> a;
  double secs = 0.0;
  double bytes() const { return double(a.size()) * 4.0; }
  double lookup(double hq, double xq, double lq) const {
    auto cell = [](std::vector<double> const& v, double q, int& i) {
      int const n = int(v.size());
      i = 0;
      while (i + 2 < n && v[i + 1] < q) ++i;
      double const d = v[i + 1] - v[i];
      double const t = (d > 0.0) ? (q - v[i]) / d : 0.0;
      return std::min(1.0, std::max(0.0, t));  // clamp: query-time semantics preserved
    };
    int i0, x0, l0;
    double const th = cell(ax.h, hq, i0);
    double const tx = cell(ax.xi, xq, x0);
    double const tl = cell(ax.lam, lq, l0);
    double acc = 0.0;
    for (int p = 0; p < 2; ++p)
      for (int q = 0; q < 2; ++q)
        for (int c = 0; c < 2; ++c) {
          double const w = (p ? th : 1.0 - th) * (q ? tx : 1.0 - tx) * (c ? tl : 1.0 - tl);
          if (w == 0.0) continue;
          acc += w * double(a[std::size_t(i0 + p) * nX * nL + (x0 + q) * nL + (l0 + c)]);
        }
    return acc;
  }
};
static Table build(Axes ax, RefractiveIndexTable const& ri, double* secs) {
  Table t;
  t.ax = std::move(ax);
  t.nH = int(t.ax.h.size());
  t.nX = int(t.ax.xi.size());
  t.nL = int(t.ax.lam.size());
  t.a.assign(std::size_t(t.nH) * t.nX * t.nL, 0.0f);
  auto const t0 = std::chrono::steady_clock::now();
  for (int ix = 0; ix < t.nX; ++ix) {
    for (int il = 0; il < t.nL; ++il) {
      double const lam = t.ax.lam[il];
      double const f = Fdisp(lam);
      double const b = t.ax.xi[ix] * bCrit(ri, f);
      double acc = 0.0;
      double zPrev = R;
      for (int ih = 0; ih < t.nH; ++ih) {
        double const zTop = R + t.ax.h[ih];
        constexpr int M = 60;
        double const w0 = std::sqrt(std::max(0.0, zPrev - R));
        double const w1 = std::sqrt(std::max(0.0, zTop - R));
        double const dw = (w1 - w0) / double(M);
        for (int k = 0; k < M; ++k) {
          double const w = w0 + (double(k) + 0.5) * dw;
          double const z = R + w * w;
          double const s = b / (nPhase(ri, z, f) * z);
          if (s >= 1.0) break;
          acc += alphaOf(lam, w * w) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * (2.0 * w * dw);
        }
        t.a[std::size_t(ih) * t.nX * t.nL + ix * t.nL + il] = float(acc);
        zPrev = zTop;
      }
    }
  }
  t.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  *secs = t.secs;
  return t;
}
static double hwm_kB() {
  std::ifstream st("/proc/self/status");
  std::string l;
  double v = 0.0;
  while (std::getline(st, l))
    if (l.rfind("VmHWM", 0) == 0) std::sscanf(l.c_str() + 6, "%lf", &v);
  return v;
}

// ---- development ray set ------------------------------------------------------------
// Axis-INDEPENDENT, so every variant is measured on the same rays. This is a development
// set: it is NOT the frozen final-validation set and is not used to tune anything that
// the final validation will then be judged on without fresh evidence.
struct Ray {
  double h, lam, xi;
  int fam;
};
static std::vector<Ray> devRays() {
  std::vector<Ray> v;
  std::vector<double> x1;
  for (double x = 0.5; x <= 0.95001; x += 0.05) x1.push_back(x);
  for (double x : {0.96, 0.97, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9999, 0.99995})
    x1.push_back(x);
  for (double h : {375.0, 1500.0, 7000.0})
    for (double lam : {350.0, 500.0, 690.0})
      for (double x : x1) v.push_back({h, lam, x, 1});
  double const offs[] = {-1e-3, 1e-3, 1.0, 10.0, 24.999, 25.0, 25.001, 26.0, 100.0};
  for (double k : {0.0, 25.0, 50.0, 100.0, 250.0, 500.0, 1000.0, 5000.0, 10000.0, 20000.0, 30000.0,
                   40000.0, 50000.0})
    for (double o : offs) {
      double const h = k + o;
      if (h < 0.0) continue;
      for (double lam : {400.0, 690.0})
        for (double x : {0.5, 0.9, 0.99}) v.push_back({h, lam, x, 2});
    }
  for (double h : {0.0, 1e-6, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0})
    for (double lam : {350.0, 500.0, 690.0})
      for (double x : {0.3, 0.8, 0.99, 0.9995}) v.push_back({h, lam, x, 3});
  for (double lam : {180.0, 180.001, 182.5, 200.0, 295.0, 300.0, 400.0, 695.0, 697.5, 699.999, 700.0,
                     700.001, 750.0, 850.0, 900.0})
    for (double h : {1500.0, 18000.0})
      for (double x : {0.5, 0.9, 0.995}) v.push_back({h, lam, x, 4});
  for (double x : {0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999})
    for (double h : {1500.0, 18000.0, 46000.0})
      for (double lam : {400.0, 690.0}) v.push_back({h, lam, x, 5});
  for (double d : {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 2e-8})
    for (double h : {1500.0, 18000.0, 46000.0})
      for (double lam : {400.0, 690.0}) v.push_back({h, lam, 1.0 - d, 6});
  return v;
}

int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  RefractiveIndexTable ri = makeRI();
  auto const wr = g_ab.getWavelengthRange();
  auto const hr = g_ab.getHeightRange();
  std::printf("=== re-axis experiment (diagnostics only; no production change) ===\n");
  std::printf("profile DECLARED domain: lambda [%.1f, %.1f] nm, h [%.1f, %.1f] m\n", wr.first,
              wr.second, hr.first, hr.second);
  std::printf("query-time clamping preserved: the lookup clamps, it never extrapolates\n");
  auto const rays = devRays();
  std::printf("development rays: %zu (F1..F6), axis-independent, NOT the frozen set\n", rays.size());

  std::vector<double> ref(rays.size(), -1.0);
  auto const r0 = std::chrono::steady_clock::now();
  for (size_t i = 0; i < rays.size(); ++i) {
    Ray const& y = rays[i];
    if (y.h <= 0.0) { ref[i] = 0.0; continue; }
    double const f = Fdisp(y.lam);
    double const p = y.xi * bCrit(ri, f);
    if (p >= nPhase(ri, R + y.h, f) * (R + y.h)) continue;  // inadmissible
    ref[i] = tauRefW(ri, p, f, y.lam, y.h, 60000);
  }
  std::printf("references computed in %.1f s\n",
              std::chrono::duration<double>(std::chrono::steady_clock::now() - r0).count());
  {
    double worst = 0.0;
    char info[200] = {0};
    int n = 0;
    for (size_t i = 0; i < rays.size(); i += 17) {
      Ray const& y = rays[i];
      if (ref[i] < 0.0 || y.h <= 0.0) continue;
      double const f = Fdisp(y.lam);
      double const p = y.xi * bCrit(ri, f);
      double const hi = tauRefW(ri, p, f, y.lam, y.h, 240000);
      ++n;
      if (std::abs(hi - ref[i]) > worst) {
        worst = std::abs(hi - ref[i]);
        std::snprintf(info, sizeof(info), "h=%.6g lam=%.6g xi=%.8f", y.h, y.lam, y.xi);
      }
    }
    std::printf("REFERENCE UNCERTAINTY (nW 60000 vs 240000, %d rays): worst |d tau| = %.3e  [%s]\n",
                n, worst, info);
  }

  std::vector<Axes> variants;
  variants.push_back({hCurrent(), xiCurrent(), lamWide(), "V0 phase-3 axes (reproduce)"});
  variants.push_back({hCurrent(), xiCurrent(), lamNative(), "V1 native lambda"});
  variants.push_back({hGraded(), xiCurrent(), lamNative(), "V2 graded h"});
  variants.push_back({hGraded(), xiRefined(416), lamNative(), "V3 refined xi (464 nodes)"});
  variants.push_back({hGraded(), xiLogOneMinus(208), lamNative(), "V4 log10(1-xi) coordinate"});
  variants.push_back({hGraded(), xiSinTheta(208), lamNative(), "V5 sin(theta) coordinate"});
  variants.push_back({hGraded2(), xiCurrent(), lamNative(), "V6 graded h, fewer nodes"});

  std::printf("\n%-30s %4s %4s %4s %8s %8s  %s\n", "variant", "h", "xi", "lam", "MB", "build_s",
              "worst |dT| by family (F1..F6)");
  double bestAll = 1e9;
  size_t bestIdx = 0;
  for (size_t vi = 0; vi < variants.size(); ++vi) {
    Axes const& ax = variants[vi];
    double s = 0.0;
    Table T = build(ax, ri, &s);
    double fw[7] = {0, 0, 0, 0, 0, 0, 0};
    char fi[7][220];
    for (int k = 0; k < 7; ++k) fi[k][0] = 0;
    for (size_t i = 0; i < rays.size(); ++i) {
      Ray const& y = rays[i];
      if (ref[i] < 0.0) continue;
      double const t = T.lookup(y.h, y.xi, y.lam);
      double const d = std::abs(std::exp(-t) - std::exp(-ref[i]));
      if (d > fw[y.fam]) {
        fw[y.fam] = d;
        std::snprintf(fi[y.fam], sizeof(fi[y.fam]), "h=%.6g lam=%.6g xi=%.8f", y.h, y.lam, y.xi);
      }
    }
    double all = 0.0;
    int wf = 1;
    for (int k = 1; k <= 6; ++k)
      if (fw[k] > all) {
        all = fw[k];
        wf = k;
      }
    std::printf("%-30s %4d %4d %4d %8.2f %8.1f  ", ax.name.c_str(), T.nH, T.nX, T.nL,
                T.bytes() / 1048576.0, s);
    for (int k = 1; k <= 6; ++k) std::printf("F%d=%.2e ", k, fw[k]);
    std::printf("\n    -> worst F%d %.3e [%s]  %s\n", wf, all, fi[wf], (all <= 1e-3) ? "PASS" : "FAIL");
    long const n = 300000;
    double acc = 0.0;
    auto const u0 = std::chrono::steady_clock::now();
    for (long k = 0; k < n; ++k)
      acc += T.lookup(1000.0 + 30000.0 * double(k % 991) / 991.0, 0.05 + 0.94 * double(k % 997) / 997.0,
                      180.0 + 520.0 * double(k % 983) / 983.0);
    auto const u1 = std::chrono::steady_clock::now();
    std::printf("    lookup %.4f us/call (acc=%.3e)\n",
                std::chrono::duration<double, std::micro>(u1 - u0).count() / double(n), acc);
    if (all < bestAll) {
      bestAll = all;
      bestIdx = vi;
    }
  }
  std::printf("\nbest variant by worst-case error: %s (%.3e)\n", variants[bestIdx].name.c_str(),
              bestAll);

  // ---- mapping check + actual-configured-ray comparison, for the best variant ---------
  {
    double s = 0.0;
    Table T = build(variants[bestIdx], ri, &s);
    CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
    std::printf("\n--- best variant %s: mapping + actual-configured-ray validation ---\n",
                variants[bestIdx].name.c_str());
    double worstMap = 0.0, worstRay = 0.0, maxFsec = 0.0;
    char im[220] = {0}, ir[220] = {0};
    int n = 0, nrej = 0;
    for (double h_e : {1500.0, 6000.0, 18000.0, 32000.0, 46000.0})
      for (double lam : {310.0, 400.0, 500.0, 620.0, 690.0})
        for (double xt : {0.15, 0.35, 0.6, 0.8, 0.92, 0.97, 0.995, 0.9995}) {
          double const f = Fdisp(lam), re = R + h_e, p = xt * bCrit(ri, f);
          double const sinT = p / (nPhase(ri, re, f) * re);
          if (sinT >= 1.0) { ++nrej; continue; }
          double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
          CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
          if (r.status != CurvedRayTable::Status::Hit) { ++nrej; continue; }
          double const mu = -std::cos(r.theta_obs);
          double const xiT = std::sqrt(std::max(0.0, 1.0 - mu * mu));
          double const tRef = tauRefW(ri, p, f, lam, h_e, 60000);
          double const tTab = T.lookup(h_e, xt, lam);
          double const dT = std::abs(std::exp(-tTab) - std::exp(-tRef));
          double const tF1 = tauRefW(ri, p, 1.0, lam, h_e, 60000);
          double const dF = std::abs(std::exp(-tRef) - std::exp(-tF1));
          ++n;
          if (std::abs(xiT - xt) > worstMap) {
            worstMap = std::abs(xiT - xt);
            std::snprintf(im, sizeof(im), "h=%.0f lam=%.0f xi=%.4f", h_e, lam, xt);
          }
          if (dT > worstRay) {
            worstRay = dT;
            std::snprintf(ir, sizeof(ir), "h=%.0f lam=%.0f xi=%.4f", h_e, lam, xt);
          }
          maxFsec = std::max(maxFsec, dF);
        }
    std::printf("  rays %d (inadmissible %d)\n", n, nrej);
    std::printf("  mapping |xi_trace - xi_independent| worst = %.3e  [%s]\n", worstMap, im);
    std::printf("  actual configured ray: worst |dT| = %.4e  [%s]  %s\n", worstRay, ir,
                (worstRay <= 1e-3) ? "PASS" : "FAIL");
    std::printf("  dispersion-only sensitivity of the reference: %.4e\n", maxFsec);
  }

  // ---- geometry table: effect on transmission / endpoints / timing (SEPARATE) --------
  std::printf("\n--- geometry table assessment (tracked separately; no geometry change) ---\n");
  {
    CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
    double wTime = 0.0, wEnd = 0.0, wTrans = 0.0;
    char it[220] = {0}, ie[220] = {0}, iT[220] = {0};
    for (double h_e : {1000.0, 1500.0, 5000.0, 6000.0, 18000.0, 32000.0, 46000.0})
      for (double lam : {400.0, 690.0})
        for (double xt : {0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999}) {
          double const f = Fdisp(lam), re = R + h_e, p = xt * bCrit(ri, f);
          double const sinT = p / (nPhase(ri, re, f) * re);
          if (sinT >= 1.0) continue;
          double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
          CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
          if (r.status != CurvedRayTable::Status::Hit) continue;
          double Lr = 0.0, Sr = 0.0, z = R;
          while (z < re - 1e-12) {
            double const st = std::min(0.5, re - z);
            double const zm = z + 0.5 * st;
            double const s2 = p / (nPhase(ri, zm, f) * zm);
            if (s2 >= 1.0) break;
            double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s2 * s2));
            Lr += inv * st;
            Sr += ri.N(zm) * inv * st;
            z += st;
          }
          double const tIndep = (Lr + f * Sr) / 299792458.0;
          double const dTns = std::abs(tIndep - r.time_s) * 1e9;
          double const muA = -std::cos(r.theta_obs);
          double const xiA = std::sqrt(std::max(0.0, 1.0 - muA * muA));
          double const dEnd = std::abs(xiA - xt) * R;
          double const alpha = alphaOf(lam, 0.5 * h_e);
          double const dTr = 1.0 - std::exp(-alpha * std::abs(Lr - r.geometric_path_m));
          if (dTns > wTime) {
            wTime = dTns;
            std::snprintf(it, sizeof(it), "h=%.0f lam=%.0f xi=%.4f", h_e, lam, xt);
          }
          if (dEnd > wEnd) {
            wEnd = dEnd;
            std::snprintf(ie, sizeof(ie), "h=%.0f lam=%.0f xi=%.4f", h_e, lam, xt);
          }
          if (dTr > wTrans) {
            wTrans = dTr;
            std::snprintf(iT, sizeof(iT), "h=%.0f lam=%.0f xi=%.4f", h_e, lam, xt);
          }
        }
    std::printf("  timing   : worst |d t| = %.4f ns   [%s]   (stated level 0.1 ns)\n", wTime, it);
    std::printf("  endpoint : worst |d sin(theta)|*R = %.4e m  [%s]   (stated level 0.01 m)\n", wEnd, ie);
    std::printf("  transmission: worst 1-exp(-alpha*dL) = %.4e  [%s]   (absorption tolerance 1e-3)\n",
                wTrans, iT);
    std::printf("  NOTE: transmission is a FIRST-ORDER alpha*dL estimate, not an integration along\n");
    std::printf("        the traced trajectory: it bounds the effect, it does not measure it.\n");
  }
  std::printf("\npeak RSS at exit: %.1f MB\n", hwm_kB() / 1024.0);
  return 0;
}
