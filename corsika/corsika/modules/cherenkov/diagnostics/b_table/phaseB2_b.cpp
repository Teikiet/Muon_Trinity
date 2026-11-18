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

// ---------- phase 0d / C / B helpers ------------------------------------------------
static std::vector<double> g_hknots;  // the profile's native height knots
static void parseKnots(char const* fn) {
  std::ifstream f(fn);
  size_t nw = 0, nh = 0;
  f >> nw >> nh;
  std::vector<double> wl(nw);
  for (size_t i = 0; i < nw; ++i) f >> wl[i];
  g_hknots.assign(nh, 0.0);
  for (size_t i = 0; i < nh; ++i) f >> g_hknots[i];
}

// exact local extinction: slope of the piecewise-linear cumulative column inside its cell
static double alphaExact(double lam, double h) {
  size_t k = 0;
  while (k + 2 < g_hknots.size() && g_hknots[k + 1] <= h) ++k;
  double const h0 = g_hknots[k], h1 = g_hknots[k + 1];
  double const t0 = g_ab.getOpticalDepth(lam, h0), t1 = g_ab.getOpticalDepth(lam, h1);
  return (h1 > h0) ? std::max(0.0, (t1 - t0) / (h1 - h0)) : 0.0;
}

// table with selectable STORAGE; lookup returns the TRANSMISSION either way
struct Tab {
  Axes ax;
  int nH = 0, nX = 0, nL = 0;
  std::vector<float> a;
  bool storeT = false;
  double secs = 0.0;
  double mib() const { return double(a.size()) * 4.0 / 1048576.0; }
  double lookupT(double hq, double xq, double lq) const {
    auto cell = [](std::vector<double> const& v, double q, int& i) {
      int const n = int(v.size());
      i = 0;
      while (i + 2 < n && v[i + 1] < q) ++i;
      double const d = v[i + 1] - v[i];
      double const t = (d > 0.0) ? (q - v[i]) / d : 0.0;
      return std::min(1.0, std::max(0.0, t));
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
    return storeT ? std::min(1.0, std::max(0.0, acc)) : std::exp(-acc);
  }
};

static Tab buildTab(Axes ax, bool storeT, RefractiveIndexTable const& ri, double* secs) {
  Tab t;
  t.ax = std::move(ax);
  t.storeT = storeT;
  t.nH = int(t.ax.h.size());
  t.nX = int(t.ax.xi.size());
  t.nL = int(t.ax.lam.size());
  t.a.assign(std::size_t(t.nH) * t.nX * t.nL, 0.0f);
  auto const t0 = std::chrono::steady_clock::now();
  for (int ix = 0; ix < t.nX; ++ix)
    for (int il = 0; il < t.nL; ++il) {
      double const lam = t.ax.lam[il], f = Fdisp(lam), b = t.ax.xi[ix] * bCrit(ri, f);
      double acc = 0.0, zPrev = R;
      for (int ih = 0; ih < t.nH; ++ih) {
        double const zTop = R + t.ax.h[ih];
        constexpr int M = 60;
        double const w0 = std::sqrt(std::max(0.0, zPrev - R));
        double const w1 = std::sqrt(std::max(0.0, zTop - R));
        double const dw = (w1 - w0) / double(M);
        for (int k = 0; k < M; ++k) {
          double const w = w0 + (double(k) + 0.5) * dw, z = R + w * w;
          double const s = b / (nPhase(ri, z, f) * z);
          if (s >= 1.0) break;
          acc += alphaOf(lam, w * w) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * (2.0 * w * dw);
        }
        t.a[std::size_t(ih) * t.nX * t.nL + ix * t.nL + il] = float(storeT ? std::exp(-acc) : acc);
        zPrev = zTop;
      }
    }
  t.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  *secs = t.secs;
  return t;
}

// the SAME quadrature as the build, evaluated directly at one point: the correct
// null-check reference (per-h-cell w-substitution, M substeps per cell)
static double tauBuild(RefractiveIndexTable const& ri, double p, double f, double lam,
                       Axes const& ax, double hTop) {
  double acc = 0.0, zPrev = R;
  for (size_t ih = 0; ih < ax.h.size() && ax.h[ih] <= hTop + 1e-12; ++ih) {
    double const zTop = R + ax.h[ih];
    constexpr int M = 60;
    double const w0 = std::sqrt(std::max(0.0, zPrev - R));
    double const w1 = std::sqrt(std::max(0.0, zTop - R));
    double const dw = (w1 - w0) / double(M);
    for (int k = 0; k < M; ++k) {
      double const w = w0 + (double(k) + 0.5) * dw, z = R + w * w;
      double const s = p / (nPhase(ri, z, f) * z);
      if (s >= 1.0) break;
      acc += alphaOf(lam, w * w) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * (2.0 * w * dw);
    }
    zPrev = zTop;
  }
  return acc;
}

// ---- frozen sets (phase 0d) --------------------------------------------------------
struct Ray2 {
  double h, lam, xi;
};
static std::vector<double> hRefined() {  // V3's h axis plus local nodes in [200, 600] m
  std::vector<double> v = hGraded();
  for (double x : {275.0, 300.0, 350.0, 400.0, 450.0, 550.0}) v.push_back(x);
  std::sort(v.begin(), v.end());
  return v;
}
static std::vector<double> xiRefinedLocal(int extra) {  // V3's xi axis plus local nodes near 1
  std::vector<double> v = xiRefined(416);
  for (int i = 1; i <= extra; ++i) {
    double const u = 0.03 * double(i) / double(extra);
    v.push_back(1.0 - u * u);
  }
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}
static std::vector<double> lamRefined() {  // V3's lambda axis plus the 687.5 nm neighbourhood
  std::vector<double> v = lamNative();
  for (double x : {682.5, 687.5, 692.5}) v.push_back(x);
  std::sort(v.begin(), v.end());
  return v;
}
static std::vector<Ray2> denseGrid() {  // TUNING grid, disposable
  std::vector<Ray2> v;
  for (double h : {256.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0})
    for (double xi : {0.99985, 0.99988, 0.99990, 0.99991, 0.99992, 0.99994, 0.9999555381,
                      0.99996, 0.99998})
      for (double lam : {682.5, 685.0, 687.5, 690.0, 692.5}) v.push_back({h, lam, xi});
  return v;
}
static std::vector<Ray2> setI() {  // HELD-OUT, in-band, fresh, off the tuning grid
  std::vector<Ray2> v;
  for (double h : {262.0, 288.0, 312.0, 341.0, 366.0, 402.0, 431.0, 463.0, 487.0})
    for (double xi : {0.9998620, 0.9998850, 0.9999025, 0.9999310, 0.9999475, 0.9999640, 0.9999880})
      for (double lam : {681.5, 684.0, 686.5, 688.5, 691.0, 693.5}) v.push_back({h, lam, xi});
  return v;
}
static std::vector<Ray2> setII() {  // HELD-OUT, out-of-band
  std::vector<Ray2> v;
  for (double h : {100.0, 1000.0, 5000.0, 20000.0, 40000.0})
    for (double xi : {0.30, 0.80, 0.99, 0.999})
      for (double lam : {255.0, 405.0, 555.0, 688.0}) v.push_back({h, lam, xi});
  return v;
}
static void writeSet(char const* path, std::vector<Ray2> const& v, char const* name) {
  std::ofstream o(path);
  o << "# " << name << " -- frozen held-out set, " << v.size() << " rays\n# h_m lambda_nm xi\n";
  o.precision(12);
  for (Ray2 const& r : v) o << r.h << " " << r.lam << " " << r.xi << "\n";
}

static std::vector<double> hRefinedN(int extra) {  // n extra h nodes evenly in [200, 600] m
  std::vector<double> v = hGraded();
  for (int i = 1; i <= extra; ++i) v.push_back(200.0 + 400.0 * double(i) / double(extra + 1));
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}
static std::vector<Ray2> setIII() {  // FRESH in-band acceptance set (frozen, never used to select)
  std::vector<Ray2> v;
  for (double h : {270.0, 295.0, 320.0, 355.0, 390.0, 420.0, 455.0, 480.0})
    for (double xi : {0.9998700, 0.9998950, 0.9999200, 0.9999425, 0.9999700, 0.9999900})
      for (double lam : {683.5, 686.5, 688.5, 691.5}) v.push_back({h, lam, xi});
  return v;
}
static std::vector<Ray2> bandScan() {  // low-h x high-xi x ALL 105 native lambda knots
  std::vector<Ray2> v;
  for (double h : {100.0, 150.0, 200.0, 256.0, 300.0, 400.0, 500.0, 700.0, 1000.0})
    for (double xi : {0.9995, 0.9997, 0.9998, 0.9999, 0.99995, 0.99998, 0.999995})
      for (int iw = 0; iw < 105; ++iw) v.push_back({h, 180.0 + 5.0 * iw, xi});
  return v;
}

int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  parseKnots(argv[1]);
  RefractiveIndexTable ri = makeRI();
  std::printf("=== phase B follow-up: V3 direct, fresh acceptance set, band scan ===\n");

  auto const tun = denseGrid();
  auto const s1 = setI();
  auto const s2 = setII();
  auto const s3 = setIII();
  writeSet("/tmp/heldout_setIII_acceptance.txt", s3, "Set III fresh in-band acceptance");
  std::printf("\n--- V3 size reconciliation ---\n");
  std::printf("  phase-3 axes (candidate A / V0): 56 x 256 x 121 = %zu nodes = %.2f MiB\n",
              size_t(56) * 256 * 121, double(56 * 256 * 121) * 4.0 / 1048576.0);
  std::printf("  V3 axes                        : 61 x 464 x 105 = %zu nodes = %.2f MiB\n",
              size_t(61) * 464 * 105, double(61 * 464 * 105) * 4.0 / 1048576.0);
  std::printf("  +lambda only (V3 + 3 lambda)   : 61 x 464 x 108 = %zu nodes = %.2f MiB\n",
              size_t(61) * 464 * 108, double(61 * 464 * 108) * 4.0 / 1048576.0);
  std::printf("  => the 6.62 MiB figure is the PHASE-3 axes, not V3; no contradiction.\n");

  std::printf("\n--- Set III (fresh in-band acceptance, frozen before selection) ---\n");
  std::printf("  %zu rays -> /tmp/heldout_setIII_acceptance.txt\n", s3.size());
  {
    int clash = 0;
    for (Ray2 const& a : s3)
      for (Ray2 const& b : tun)
        if (a.h == b.h && a.lam == b.lam && a.xi == b.xi) ++clash;
    int clash1 = 0;
    for (Ray2 const& a : s3)
      for (Ray2 const& b : s1)
        if (a.h == b.h && a.lam == b.lam && a.xi == b.xi) ++clash1;
    std::printf("  clashes vs tuning grid %d, vs Set I %d (both must be 0)\n", clash, clash1);
  }

  // references
  struct Grid {
    const char* name;
    std::vector<Ray2> rays;
    std::vector<double> tref;
  };
  std::vector<Grid> grids;
  grids.push_back({"tuning", tun, {}});
  grids.push_back({"SetI", s1, {}});
  grids.push_back({"SetII", s2, {}});
  grids.push_back({"SetIII", s3, {}});
  for (Grid& g : grids) {
    g.tref.assign(g.rays.size(), -1.0);
    for (size_t i = 0; i < g.rays.size(); ++i) {
      Ray2 const& r = g.rays[i];
      if (r.h <= 0.0) { g.tref[i] = 0.0; continue; }
      double const f = Fdisp(r.lam), p = r.xi * bCrit(ri, f);
      if (p >= nPhase(ri, R + r.h, f) * (R + r.h)) continue;
      g.tref[i] = tauRefW(ri, p, f, r.lam, r.h, 60000);
    }
    std::printf("  reference cached for %-7s (%zu rays)\n", g.name, g.rays.size());
  }

  std::printf("\n--- V3 direct vs the passing candidate, tau storage, on all four grids ---\n");
  struct Cfg {
    const char* name;
    Axes ax;
  };
  std::vector<Cfg> cfgs;
  cfgs.push_back({"V3 (61x464x105)", {hGraded(), xiRefined(416), lamNative(), "V3"}});
  cfgs.push_back({"+h+xi (67x664x105)", {hRefinedN(6), xiRefinedLocal(200), lamNative(), "cand"}});
  for (Cfg const& c : cfgs) {
    double s = 0.0;
    Tab T = buildTab(c.ax, false, ri, &s);
    std::printf("\n  %s   %d x %d x %d = %.2f MiB   build %.1f s\n", c.name, T.nH, T.nX, T.nL,
                T.mib(), s);
    for (Grid const& g : grids) {
      double w = 0.0;
      char loc[200] = {0};
      for (size_t i = 0; i < g.rays.size(); ++i) {
        if (g.tref[i] < 0.0) continue;
        Ray2 const& r = g.rays[i];
        double const d = std::abs(T.lookupT(r.h, r.xi, r.lam) - std::exp(-g.tref[i]));
        if (d > w) {
          w = d;
          std::snprintf(loc, sizeof(loc), "h=%.6g lam=%.6g xi=%.10f", r.h, r.lam, r.xi);
        }
      }
      std::printf("    %-7s worst |dT| = %.4e  %s   [%s]\n", g.name, w, (w <= 1e-3) ? "PASS" : "FAIL",
                  loc);
    }
  }

  // ---------------- band scan: other low-h, high-xi failure regions ------------------
  std::printf("\n--- band scan: low h x high xi x ALL 105 native lambda knots ---\n");
  {
    auto const bs = bandScan();
    std::printf("  %zu rays\n", bs.size());
    std::vector<double> tref(bs.size(), -1.0);
    for (size_t i = 0; i < bs.size(); ++i) {
      Ray2 const& r = bs[i];
      double const f = Fdisp(r.lam), p = r.xi * bCrit(ri, f);
      if (p >= nPhase(ri, R + r.h, f) * (R + r.h)) continue;
      tref[i] = tauRefW(ri, p, f, r.lam, r.h, 60000);
    }
    for (Cfg const& c : cfgs) {
      double s = 0.0;
      Tab T = buildTab(c.ax, false, ri, &s);
      std::vector<double> perLam(105, 0.0);
      double worst = 0.0;
      char loc[200] = {0};
      for (size_t i = 0; i < bs.size(); ++i) {
        if (tref[i] < 0.0) continue;
        Ray2 const& r = bs[i];
        double const d = std::abs(T.lookupT(r.h, r.xi, r.lam) - std::exp(-tref[i]));
        int const iw = int((r.lam - 180.0) / 5.0 + 0.5);
        if (iw >= 0 && iw < 105) perLam[iw] = std::max(perLam[iw], d);
        if (d > worst) {
          worst = d;
          std::snprintf(loc, sizeof(loc), "h=%.6g lam=%.6g xi=%.10f", r.h, r.lam, r.xi);
        }
      }
      std::printf("  %-20s overall worst %.4e %s  [%s]\n", c.name, worst,
                  (worst <= 1e-3) ? "PASS" : "FAIL", loc);
      // top-5 worst wavelength knots: are they absorption bands?
      std::vector<int> idx(105);
      for (int i = 0; i < 105; ++i) idx[i] = i;
      std::sort(idx.begin(), idx.end(),
                [&](int a, int b) { return perLam[a] > perLam[b]; });
      std::printf("    worst lambda knots:");
      for (int k = 0; k < 5; ++k)
        std::printf("  %.0f nm:%.3e", 180.0 + 5.0 * idx[k], perLam[idx[k]]);
      std::printf("\n");
      int nBad = 0;
      for (int i = 0; i < 105; ++i)
        if (perLam[i] > 1e-3) ++nBad;
      std::printf("    lambda knots with any ray over the gate: %d of 105\n", nBad);
    }
  }

  // ---------------- refinement-level sweep ------------------------------------------
  std::printf("\n--- refinement-level sweep (two levels per axis, so the coupling rests on more) ---\n");
  struct Sweep {
    const char* name;
    int nh, nxi;
  };
  Sweep const sw[] = {{"h+0   xi+0    (=V3)", 0, 0},   {"h+0   xi+400", 0, 400},
                      {"h+12  xi+0", 12, 0},          {"h+12  xi+400", 12, 400},
                      {"h+6   xi+200  (cand)", 6, 200}};
  std::printf("%-22s %5s %5s %7s  %-11s %-11s\n", "config", "h", "xi", "MiB", "SetI", "SetIII");
  for (Sweep const& q : sw) {
    Axes ax{q.nh ? hRefinedN(q.nh) : hGraded(),
            q.nxi ? xiRefinedLocal(q.nxi) : xiRefined(416), lamNative(), "sw"};
    double s = 0.0;
    Tab T = buildTab(ax, false, ri, &s);
    double w1 = 0.0, w3 = 0.0;
    for (Grid const& g : grids) {
      bool const isI = (std::string(g.name) == "SetI");
      bool const isIII = (std::string(g.name) == "SetIII");
      if (!isI && !isIII) continue;
      for (size_t i = 0; i < g.rays.size(); ++i) {
        if (g.tref[i] < 0.0) continue;
        Ray2 const& r = g.rays[i];
        double const d = std::abs(T.lookupT(r.h, r.xi, r.lam) - std::exp(-g.tref[i]));
        if (isI) w1 = std::max(w1, d);
        else w3 = std::max(w3, d);
      }
    }
    std::printf("%-22s %5d %5d %7.2f  %-11.3e %-11.3e\n", q.name, T.nH, T.nX, T.mib(), w1, w3);
  }
  std::printf("\npeak RSS at exit: %.1f MiB\n", hwm_kB() / 1024.0);
  return 0;
}

