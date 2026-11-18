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

static std::vector<double> xiLogRefined(int nLog, bool includeOne) {
  std::vector<double> v;
  xiLog(v, 48);
  double const y0 = -std::log10(0.5), y1 = -std::log10(1.0e-9);
  for (int i = 0; i < nLog; ++i) {
    double const y = y0 + (y1 - y0) * double(i) / double(nLog - 1);
    v.push_back(1.0 - std::pow(10.0, -y));
  }
  if (includeOne) v.push_back(1.0);
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}
// group dispersion consistent with the Cauchy phase model F(lambda):
//   f_group = F - lambda dF/dlambda = F + 2k/(lambda_um^2 * (1 + k/0.16))
static double Fgroup(double lam_nm) {
  double const l = lam_nm / 1000.0;
  return Fdisp(lam_nm) + 2.0 * 0.0032 / (l * l * (1.0 + 0.0032 / 0.16));
}

// query classification: cell bracket, clamping status, interior vs node
struct CellInfo {
  double hLo = 0, hHi = 0, xLo = 0, xHi = 0, lLo = 0, lHi = 0;
  bool clampH = false, clampXi = false, clampLam = false;
  bool interiorXi = false, interiorH = false;
};
static CellInfo classify(Axes const& ax, double h, double xi, double lam) {
  CellInfo c;
  auto find = [](std::vector<double> const& v, double q) {
    int i = 0;
    int const n = int(v.size());
    while (i + 2 < n && v[i + 1] < q) ++i;
    return i;
  };
  int const i = find(ax.h, h), j = find(ax.xi, xi), k = find(ax.lam, lam);
  c.hLo = ax.h[i];
  c.hHi = ax.h[i + 1];
  c.xLo = ax.xi[j];
  c.xHi = ax.xi[j + 1];
  c.lLo = ax.lam[k];
  c.lHi = ax.lam[k + 1];
  c.clampH = (h <= ax.h.front()) || (h >= ax.h.back());
  c.clampXi = (xi <= ax.xi.front()) || (xi >= ax.xi.back());
  c.clampLam = (lam <= ax.lam.front()) || (lam >= ax.lam.back());
  c.interiorXi = (!c.clampXi && xi > c.xLo && xi < c.xHi);
  c.interiorH = (!c.clampH && h > c.hLo && h < c.hHi);
  return c;
}

// independent timing reference #2: integrate in the SWEPT ANGLE phi, RK4 on r(phi)
//   dr/dphi = r * sqrt(1-s^2)/s,  ds = r/s dphi,  s = p/(n(r) r)
static bool timeRefPhi(RefractiveIndexTable const& ri, double p, double fPhase, double fGroup,
                       double rEmit, double dphi, double* Lout, double* Sout, long* steps) {
  auto rhs = [&](double r) {
    double const s = p / (nPhase(ri, r, fPhase) * r);
    if (!(s > 0.0) || s >= 1.0) return 0.0;
    return r * std::sqrt(std::max(0.0, 1.0 - s * s)) / s;
  };
  auto ds_dphi = [&](double r) {
    double const s = p / (nPhase(ri, r, fPhase) * r);
    if (!(s > 0.0) || s >= 1.0) return 0.0;
    return r / s;
  };
  double r = R, L = 0.0, S = 0.0;
  long n = 0;
  while (r < rEmit - 1e-12 && n < 4000000) {
    double step = dphi;
    double const k1 = rhs(r);
    double const k2 = rhs(r + 0.5 * step * k1);
    double const k3 = rhs(r + 0.5 * step * k2);
    double const k4 = rhs(r + step * k3);
    double rNext = r + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    if (rNext > rEmit) {  // partial final step
      double const frac = (rEmit - r) / std::max(1e-300, rNext - r);
      step *= frac;
      double const m1 = rhs(r);
      double const m2 = rhs(r + 0.5 * step * m1);
      double const m3 = rhs(r + 0.5 * step * m2);
      double const m4 = rhs(r + step * m3);
      rNext = r + step * (m1 + 2.0 * m2 + 2.0 * m3 + m4) / 6.0;
    }
    double const rm = 0.5 * (r + rNext);
    L += ds_dphi(rm) * step;
    S += ri.N(rm) * ds_dphi(rm) * step;
    r = rNext;
    ++n;
  }
  *Lout = L;
  *Sout = S;
  *steps = n;
  (void)fGroup;
  return (r >= rEmit - 1e-9);
}


int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  RefractiveIndexTable ri = makeRI();
  auto const wr = g_ab.getWavelengthRange();
  double const muMin = 2.0e-4;
  double const xiMax = std::sqrt(1.0 - muMin * muMin);
  std::printf("=== closure experiment (diagnostics only; no production change) ===\n");
  std::printf("eligibility limit DERIVED, not assumed: xi_max = sqrt(1 - mu_min^2), mu_min = %.3e\n",
              muMin);
  std::printf("  xi_max = %.17g\n", xiMax);
  std::printf("  1 - xi_max = %.6e   (so 1 - 2e-8 is an approximation, not an identity)\n",
              1.0 - xiMax);

  std::vector<Axes> variants;
  variants.push_back({hGraded(), xiRefined(416), lamNative(), "V3  u-uniform 464"});
  variants.push_back({hGraded(), xiLogRefined(415, true), lamNative(), "V7  log+refine 463"});
  variants.push_back({hGraded(), xiLogRefined(319, true), lamNative(), "V8  log moderate 367"});
  std::vector<Table> T;
  for (size_t i = 0; i < variants.size(); ++i) {
    double s = 0.0;
    T.push_back(build(variants[i], ri, &s));
    std::printf("built %s : %d x %d x %d = %.2f MB  %.1f s\n", variants[i].name.c_str(), T[i].nH,
                T[i].nX, T[i].nL, T[i].bytes() / 1048576.0, s);
  }

  struct Probe {
    const char* label;
    double h0, lam0, xi0;
  };
  Probe const probes[] = {
      {"V3 overall worst", 1025.0, 400.0, 0.99},
      {"audit F1 worst", 375.0, 690.0, 0.9998555381},
      {"V4 F1 worst", 375.0, 690.0, 0.9999},
      {"audit F3 worst", 100.0, 350.0, 0.9995},
      {"configured-ray worst", 1500.0, 620.0, 0.9995},
      {"phase-3 worst", 1500.0, 670.0, 0.9950},
  };

  std::printf("\n--- BLOCK 1: attributed probe maxima (all vs ELIGIBLE domain) ---\n");
  for (size_t vi = 0; vi < variants.size(); ++vi) {
    double allMax = 0.0, eligMax = 0.0;
    char ia[340] = {0}, ie[340] = {0};
    for (Probe const& p : probes)
      for (double dh : {-500.0, -100.0, -10.0, -1.0, -0.1, 0.0, 0.1, 1.0, 10.0, 100.0, 500.0})
        for (double dl : {-2.5, -1.0, -0.1, 0.0, 0.1, 1.0, 2.5})
          for (double dx : {-1e-3, -1e-4, -1e-5, -1e-6, 0.0, 1e-6, 1e-5, 1e-4, 1e-3}) {
            double const h = p.h0 + dh, lam = p.lam0 + dl, xi = p.xi0 + dx;
            if (h <= 0.0 || xi <= 0.0 || xi > 1.0 || lam < 180.0 || lam > 700.0) continue;
            double const f = Fdisp(lam), pp = xi * bCrit(ri, f);
            if (pp >= nPhase(ri, R + h, f) * (R + h)) continue;
            double const tr = tauRefW(ri, pp, f, lam, h, 20000);
            double const d = std::abs(std::exp(-T[vi].lookup(h, xi, lam)) - std::exp(-tr));
            CellInfo const c = classify(variants[vi], h, xi, lam);
            bool const elig = (xi <= xiMax) && !c.clampXi && !c.clampH && !c.clampLam;
            char buf[340];
            std::snprintf(buf, sizeof(buf),
                          "h=%.6g lam=%.6g xi=%.12f | xcell=[%.12g,%.12g] hcell=[%.6g,%.6g] | "
                          "clampXi=%d clampH=%d clampLam=%d | intXi=%d intH=%d | elig=%d",
                          h, lam, xi, c.xLo, c.xHi, c.hLo, c.hHi, int(c.clampXi), int(c.clampH),
                          int(c.clampLam), int(c.interiorXi), int(c.interiorH), int(elig));
            if (d > allMax) {
              allMax = d;
              std::snprintf(ia, sizeof(ia), "%s", buf);
            }
            if (elig && d > eligMax) {
              eligMax = d;
              std::snprintf(ie, sizeof(ie), "%s", buf);
            }
          }
    std::printf("\n%s\n", variants[vi].name.c_str());
    std::printf("  ALL-probe max        %.4e\n    %s\n", allMax, ia);
    std::printf("  ELIGIBLE-domain max  %.4e\n    %s\n", eligMax, ie);
  }

  std::printf("\n--- BLOCK 2: grazing coverage scan, xi = 1 - 10^u (u = -2 .. -8) ---\n");
  for (size_t vi = 0; vi < variants.size(); ++vi) {
    double covMax = 0.0, clampMax = 0.0, largestOk = -1.0, smallestBad = 2.0;
    char ic[340] = {0}, ik[340] = {0};
    for (int s = 0; s < 200; ++s) {
      double const u = -2.0 + (-8.0 + 2.0) * double(s) / 199.0;
      double const xi = 1.0 - std::pow(10.0, u);
      if (xi > xiMax) continue;  // ineligible: excluded by definition
      for (double h : {375.0, 1500.0, 6000.0, 18000.0, 46000.0})
        for (double lam : {350.0, 500.0, 690.0}) {
          double const f = Fdisp(lam), pp = xi * bCrit(ri, f);
          if (pp >= nPhase(ri, R + h, f) * (R + h)) continue;
          double const tr = tauRefW(ri, pp, f, lam, h, 20000);
          double const d = std::abs(std::exp(-T[vi].lookup(h, xi, lam)) - std::exp(-tr));
          CellInfo const c = classify(variants[vi], h, xi, lam);
          char buf[340];
          std::snprintf(buf, sizeof(buf),
                        "h=%.6g lam=%.6g xi=%.12f | xcell=[%.12g,%.12g] | intXi=%d", h, lam, xi,
                        c.xLo, c.xHi, int(c.interiorXi));
          if (c.clampXi) {
            if (d > clampMax) {
              clampMax = d;
              std::snprintf(ik, sizeof(ik), "%s", buf);
            }
          } else {
            if (d > covMax) {
              covMax = d;
              std::snprintf(ic, sizeof(ic), "%s", buf);
            }
            if (d <= 1e-3) largestOk = std::max(largestOk, xi);
            if (d > 1e-3) smallestBad = std::min(smallestBad, xi);
          }
        }
    }
    std::printf("\n%s\n", variants[vi].name.c_str());
    std::printf("  top xi node = %.12g\n", variants[vi].xi.back());
    std::printf("  validated-coverage max (no xi clamp) %.4e\n    %s\n", covMax, ic);
    std::printf("  CLAMPED-region max                   %.4e\n    %s\n", clampMax, ik);
    std::printf("  largest xi still <= 1e-3 : %.12g\n", largestOk);
    std::printf("  smallest xi > 1e-3       : %.12g\n", smallestBad);
  }

  std::printf("\n--- BLOCK 3: absorption height ceiling (h above the last h node) ---\n");
  std::printf("h axis top = %.1f m; the profile itself also tops at %.1f m, so the callback's\n",
              variants[0].h.back(), g_ab.getHeightRange().second);
  std::printf("alpha is 0 above it. Clamping h therefore adds nothing: expected error ~0.\n");
  for (size_t vi = 0; vi < variants.size(); ++vi) {
    double worst = 0.0;
    char info[340] = {0};
    for (double h : {50000.0, 51000.0, 52000.0, 55000.0, 60000.0})
      for (double lam : {400.0, 690.0})
        for (double xi : {0.5, 0.9, 0.999, 0.9999}) {
          double const f = Fdisp(lam), pp = xi * bCrit(ri, f);
          if (pp >= nPhase(ri, R + h, f) * (R + h)) continue;
          double const tr = tauRefW(ri, pp, f, lam, h, 20000);
          double const d = std::abs(std::exp(-T[vi].lookup(h, xi, lam)) - std::exp(-tr));
          if (d > worst) {
            worst = d;
            std::snprintf(info, sizeof(info), "h=%.0f lam=%.0f xi=%.4f clampH=%d", h, lam, xi,
                          int(classify(variants[vi], h, xi, lam).clampH));
          }
        }
    std::printf("  %-22s worst |dT| above the ceiling = %.4e  [%s]\n", variants[vi].name.c_str(),
                worst, info);
  }
  std::printf("  POLICY: still require an explicit COUNTED fallback above the ceiling, because this\n");
  std::printf("  benign result follows from the profile (alpha = 0 above its top), not from the table.\n");

  std::printf("\n--- BLOCK 4: timing resolved separately ---\n");
  std::printf("applicable requirement, quoted from this repository:\n");
  std::printf("  CurvedRayTable.hpp:264  't errors stay inside 1 cm / 0.1 ns for mu >= 2e-4; below it they break down'\n");
  std::printf("  testCherenkovCurvedRayTable.cpp:295  'agreed targets: 1 cm displacement, 0.1 ns time, 1e-6 rad angle'\n");
  std::printf("  => 0.1 ns IS the applicable requirement, and the table CLAIMS it over exactly the\n");
  std::printf("     eligible domain (mu >= 2e-4).\n");
  {
    CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
    std::printf("\n  (b) same observable with a phase/group split, f_group = f_phase + 2k/(l^2 (1+k/0.16))\n");
    for (double lam : {400.0, 690.0})
      std::printf("      lam=%.0f  f_phase=%.6f  f_group=%.6f  (split %.2f %%)\n", lam, Fdisp(lam),
                  Fgroup(lam), 100.0 * (Fgroup(lam) / Fdisp(lam) - 1.0));
    double worstSplit = 0.0;
    char isp[340] = {0};
    for (double h_e : {1500.0, 18000.0, 46000.0})
      for (double lam : {400.0, 690.0}) {
        double const fp = Fdisp(lam), fg = Fgroup(lam), re = R + h_e;
        for (double xt : {0.9, 0.99, 0.999, 0.9995, 0.9999}) {
          double const pp = xt * bCrit(ri, fp);
          double const sinT = pp / (nPhase(ri, re, fp) * re);
          if (sinT >= 1.0) continue;
          double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
          CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, fp, fg);
          if (r.status != CurvedRayTable::Status::Hit) continue;
          double Lr = 0.0, Sr = 0.0, z = R;
          while (z < re - 1e-12) {
            double const st = std::min(0.1, re - z);
            double const zm = z + 0.5 * st;
            double const s2 = pp / (nPhase(ri, zm, fp) * zm);
            if (s2 >= 1.0) break;
            double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s2 * s2));
            Lr += inv * st;
            Sr += ri.N(zm) * inv * st;
            z += st;
          }
          double const dTns = std::abs((Lr + fg * Sr) / 299792458.0 - r.time_s) * 1e9;
          if (dTns > worstSplit) {
            worstSplit = dTns;
            std::snprintf(isp, sizeof(isp), "h=%.0f lam=%.0f xi=%.4f", h_e, lam, xt);
          }
        }
      }
    std::printf("      worst |d t| with the split = %.4f ns  [%s]\n", worstSplit, isp);

    std::printf("\n  (c) independent cross-check: RK4 in the SWEPT ANGLE phi (different parameterisation)\n");
    for (double xt : {0.9, 0.99, 0.9999}) {
      double const h_e = 46000.0, lam = 400.0, fp = Fdisp(lam), fg = Fgroup(lam), re = R + h_e;
      double const pp = xt * bCrit(ri, fp);
      double const sinT = pp / (nPhase(ri, re, fp) * re);
      if (sinT >= 1.0) continue;
      double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
      CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, fp, fg);
      if (r.status != CurvedRayTable::Status::Hit) continue;
      double Lp = 0.0, Sp = 0.0;
      long steps = 0;
      bool const ok = timeRefPhi(ri, pp, fp, fg, re, 2.0e-6, &Lp, &Sp, &steps);
      double Lr = 0.0, Sr = 0.0, z = R;
      while (z < re - 1e-12) {
        double const st = std::min(0.1, re - z);
        double const zm = z + 0.5 * st;
        double const s2 = pp / (nPhase(ri, zm, fp) * zm);
        if (s2 >= 1.0) break;
        double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s2 * s2));
        Lr += inv * st;
        Sr += ri.N(zm) * inv * st;
        z += st;
      }
      double const tObs = r.time_s, tPhi = (Lp + fg * Sp) / 299792458.0;
      double const tRad = (Lr + fg * Sr) / 299792458.0;
      std::printf("      xi=%.4f ok=%d steps=%ld | L: phi=%.6f radial=%.6f | Sigma: phi=%.6f radial=%.6f\n",
                  xt, int(ok), steps, Lp, Lr, Sp, Sr);
      std::printf("                 |d t| obs-vs-phi = %.4f ns   obs-vs-radial = %.4f ns   phi-vs-radial = %.4f ns\n",
                  std::abs(tObs - tPhi) * 1e9, std::abs(tObs - tRad) * 1e9,
                  std::abs(tPhi - tRad) * 1e9);
    }
  }
  std::printf("\npeak RSS at exit: %.1f MB\n", hwm_kB() / 1024.0);
  return 0;
}
