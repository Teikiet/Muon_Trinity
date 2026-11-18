// audit_b.cpp -- boundary accuracy expansion, callback-level breakpoint analysis,
// admissibility argument, candidate-B elimination estimate, performance completion.
// Diagnostic only; no production code is touched.
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

static constexpr double R = 6371000.0;
static constexpr double H_TOP = 60000.0;

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

// ---- reference integrators -------------------------------------------------------
// (1) plain uniform-step radial integral (what the table build used)
static double tauRefU(RefractiveIndexTable const& ri, double p, double f, double lam, double hTop,
                      double dr) {
  double acc = 0.0, z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zm = z + 0.5 * st;
    double const s = p / (nPhase(ri, zm, f) * zm);
    if (s >= 1.0) return acc;
    acc += alphaOf(lam, zm - R) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
    z += st;
  }
  return acc;
}
// (2) w = sqrt(z - R) substitution: dz = 2w dw removes the sqrt(1-xi^2) stiffness at
//     the detector, so this is the trustworthy reference even at xi -> 1.
static double tauRefW(RefractiveIndexTable const& ri, double p, double f, double lam, double hTop,
                      long nW) {
  double const wTop = std::sqrt(hTop);
  double const dw = wTop / double(nW);
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

// ---- candidate A table with selectable BUILD quadrature ---------------------------
struct TableA {
  std::vector<double> h, xi, lam;
  int nH = 0, nX = 0, nL = 0;
  std::vector<float> a;
  double bytes() const { return double(a.size()) * 4.0; }
  double lookup(double hq, double xq, double lq) const {
    auto cell = [](std::vector<double> const& v, double q, int& i) {
      int const n = int(v.size());
      i = 0;
      while (i + 2 < n && v[i + 1] < q) ++i;
      double const d = v[i + 1] - v[i];
      double const t = (d > 0.0) ? (q - v[i]) / d : 0.0;
      return std::min(1.0, std::max(0.0, t));
    };
    int i0, x0, l0;
    double const th = cell(h, hq, i0);
    double const tx = cell(xi, xq, x0);
    double const tl = cell(lam, lq, l0);
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

static std::vector<double> makeAxesXi() {
  std::vector<double> v;
  for (int i = 0; i < 48; ++i) v.push_back(2.0e-4 * std::pow(0.5 / 2.0e-4, double(i) / 47.0));
  for (int i = 1; i <= 208; ++i) {
    double const u = std::sqrt(0.5) * double(i) / 208.0;
    v.push_back(1.0 - u * u);
  }
  std::sort(v.begin(), v.end());
  return v;
}

// quad = 0: uniform dr = 5 m (the phase-3 build)
// quad = 1: w = sqrt(z-R) substitution inside each h cell (resolves xi -> 1)
static TableA buildTable(RefractiveIndexTable const& ri, int quad, double* secs) {
  TableA t;
  t.h = {0.0, 25.0, 50.0, 100.0, 250.0, 500.0};
  for (int i = 1; i <= 50; ++i) t.h.push_back(1000.0 * i);
  t.xi = makeAxesXi();
  for (double w = 300.0; w <= 900.0 + 1e-9; w += 5.0) t.lam.push_back(w);
  t.nH = int(t.h.size());
  t.nX = int(t.xi.size());
  t.nL = int(t.lam.size());
  t.a.assign(std::size_t(t.nH) * t.nX * t.nL, 0.0f);
  auto const t0 = std::chrono::steady_clock::now();
  for (int ix = 0; ix < t.nX; ++ix) {
    for (int il = 0; il < t.nL; ++il) {
      double const lam = t.lam[il];
      double const f = Fdisp(lam);
      double const b = t.xi[ix] * bCrit(ri, f);
      double acc = 0.0;
      double zPrev = R;
      for (int ih = 0; ih < t.nH; ++ih) {
        double const zTop = R + t.h[ih];
        if (quad == 0) {
          double z = zPrev;
          while (z < zTop - 1e-12) {
            double const st = std::min(5.0, zTop - z);
            double const zm = z + 0.5 * st;
            double const s = b / (nPhase(ri, zm, f) * zm);
            if (s >= 1.0) break;
            acc += alphaOf(lam, zm - R) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
            z += st;
          }
        } else {
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
        }
        t.a[std::size_t(ih) * t.nX * t.nL + ix * t.nL + il] = float(acc);
        zPrev = zTop;
      }
    }
  }
  *secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return t;
}

// ---- profile knots, read straight from the .dat header -----------------------------
static std::vector<double> g_knots;
static void parseKnots(const char* fn) {
  std::ifstream f(fn);
  size_t nw = 0, nh = 0;
  f >> nw >> nh;
  std::vector<double> wl(nw);
  for (size_t i = 0; i < nw; ++i) f >> wl[i];
  g_knots.assign(nh, 0.0);
  for (size_t i = 0; i < nh; ++i) f >> g_knots[i];
}
static double hwm_kB() {
  std::ifstream st("/proc/self/status");
  std::string l;
  double v = 0.0;
  while (std::getline(st, l))
    if (l.rfind("VmHWM", 0) == 0) std::sscanf(l.c_str() + 6, "%lf", &v);
  return v;
}

// ---- one ray: w-substituted reference vs the two candidate builds ------------------
struct Ev {
  double tRef = 0.0, tU = 0.0, tW = 0.0;
  bool ok = false;
};
static Ev evalRay(RefractiveIndexTable const& ri, TableA const& AU, TableA const& AW, double h_e,
                  double lam, double xt, long nW) {
  Ev e;
  if (h_e <= 0.0) { e.ok = true; return e; }  // zero path -> tau = 0 by construction
  double const f = Fdisp(lam);
  double const p = xt * bCrit(ri, f);
  if (p >= nPhase(ri, R + h_e, f) * (R + h_e)) return e;  // not launchable
  e.tRef = tauRefW(ri, p, f, lam, h_e, nW);
  e.tU = AU.lookup(h_e, xt, lam);
  e.tW = AW.lookup(h_e, xt, lam);
  e.ok = true;
  return e;
}
struct Fam {
  const char* name = "";
  double wU = 0.0, wW = 0.0;
  char iU[200] = {0}, iW[200] = {0};
  int n = 0;
};
static void addRay(Fam& F, RefractiveIndexTable const& ri, TableA const& AU, TableA const& AW,
                   double h_e, double lam, double xt, long nW) {
  Ev const e = evalRay(ri, AU, AW, h_e, lam, xt, nW);
  if (!e.ok) return;
  ++F.n;
  double const dU = std::abs(std::exp(-e.tU) - std::exp(-e.tRef));
  double const dW = std::abs(std::exp(-e.tW) - std::exp(-e.tRef));
  if (dU > F.wU) {
    F.wU = dU;
    std::snprintf(F.iU, sizeof(F.iU), "h=%.6g lam=%.6g xi=%.10f", h_e, lam, xt);
  }
  if (dW > F.wW) {
    F.wW = dW;
    std::snprintf(F.iW, sizeof(F.iW), "h=%.6g lam=%.6g xi=%.10f", h_e, lam, xt);
  }
}

int main(int argc, char** argv) {
  g_ab.load(argv[1]);
  parseKnots(argv[1]);
  RefractiveIndexTable ri = makeRI();
  long const NW = 60000;  // w-steps for the reference integral
  std::printf("=== boundary / breakpoint / admissibility / performance audit ===\n");
  std::printf("profile: %zu height knots, [%.1f, %.1f] m\n", g_knots.size(), g_knots.front(),
              g_knots.back());

  double gsec = 0.0;
  CurvedRayTable tab;
  {
    auto const t0 = std::chrono::steady_clock::now();
    tab = CurvedRayTable(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
    gsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  std::printf("geometry table: %.1f MB  build %.1f s   peak RSS now %.1f MB\n", tab.memoryMB(),
              gsec, hwm_kB() / 1024.0);

  double sU = 0.0, sW = 0.0;
  TableA AU = buildTable(ri, 0, &sU);
  TableA AW = buildTable(ri, 1, &sW);
  std::printf("A (uniform dr=5 m build): %.2f MB  %.1f s\n", AU.bytes() / 1048576.0, sU);
  std::printf("A (w-substituted build) : %.2f MB  %.1f s\n", AW.bytes() / 1048576.0, sW);
  std::printf("peak RSS with geometry + both absorption tables: %.1f MB\n", hwm_kB() / 1024.0);

  std::vector<Fam> fams;
  {
    Fam F;
    F.name = "F1 cell interiors (mid h, mid xi)";
    for (size_t i = 0; i + 1 < AU.h.size(); i += 4)
      for (size_t j = 0; j + 1 < AU.xi.size(); j += 4)
        for (double lam : {350.0, 500.0, 690.0})
          addRay(F, ri, AU, AW, 0.5 * (AU.h[i] + AU.h[i + 1]), lam, 0.5 * (AU.xi[j] + AU.xi[j + 1]), NW);
    fams.push_back(F);
  }
  {
    Fam F;
    F.name = "F2 h knots + neighbourhoods";
    double const offs[] = {-1e-3, 1e-3, 1.0, 10.0, 24.999, 25.0, 25.001, 26.0, 100.0};
    for (double k : {0.0, 25.0, 50.0, 100.0, 250.0, 500.0, 1000.0, 5000.0, 10000.0, 20000.0,
                     30000.0, 40000.0, 50000.0})
      for (double o : offs) {
        double const h = k + o;
        if (h < 0.0) continue;
        for (double lam : {400.0, 690.0})
          for (double xt : {0.5, 0.9, 0.99}) addRay(F, ri, AU, AW, h, lam, xt, NW);
      }
    for (double k : g_knots)
      for (double o : {-25.0, -24.999, 25.0, 25.001}) {
        double const h = k + o;
        if (h < 0.0 || h > 50000.0) continue;
        for (double lam : {400.0, 690.0})
          for (double xt : {0.5, 0.9}) addRay(F, ri, AU, AW, h, lam, xt, NW);
      }
    fams.push_back(F);
  }
  {
    Fam F;
    F.name = "F3 short paths / ground boundary";
    for (double h : {0.0, 1e-6, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0})
      for (double lam : {350.0, 500.0, 690.0})
        for (double xt : {0.3, 0.8, 0.99, 0.9995}) addRay(F, ri, AU, AW, h, lam, xt, NW);
    fams.push_back(F);
  }
  {
    Fam F;
    F.name = "F4 lambda endpoints / clamp edge";
    for (double lam : {295.0, 299.999, 300.0, 300.001, 400.0, 695.0, 697.5, 699.999, 700.0, 700.001,
                       705.0, 750.0, 850.0, 900.0, 900.001, 1000.0})
      for (double h : {1500.0, 18000.0})
        for (double xt : {0.5, 0.9, 0.995}) addRay(F, ri, AU, AW, h, lam, xt, NW);
    fams.push_back(F);
  }
  {
    Fam F;
    F.name = "F5 near-critical xi";
    for (double xt : {0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999})
      for (double h : {1500.0, 18000.0, 46000.0})
        for (double lam : {400.0, 690.0}) addRay(F, ri, AU, AW, h, lam, xt, NW);
    fams.push_back(F);
  }
  {
    Fam F;
    F.name = "F6 declared grazing guard (xi = 1 - d)";
    for (double d : {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 2e-8})
      for (double h : {1500.0, 18000.0, 46000.0})
        for (double lam : {400.0, 690.0}) addRay(F, ri, AU, AW, h, lam, 1.0 - d, NW);
    fams.push_back(F);
  }

  std::printf("\n--- boundary accuracy vs the w-substituted reference (|dT|) ---\n");
  std::printf("%-38s %6s  %-11s %-11s  %s\n", "family", "rays", "uniform", "w-build", "worst ray (w-build)");
  double worstU = 0.0, worstW = 0.0;
  for (Fam const& F : fams) {
    std::printf("%-38s %6d  %-11.3e %-11.3e  %s\n", F.name, F.n, F.wU, F.wW, F.iW);
    worstU = std::max(worstU, F.wU);
    worstW = std::max(worstW, F.wW);
  }
  std::printf("worst over all families: uniform %.3e   w-build %.3e   (gate 1e-3: %s)\n", worstU,
              worstW, (worstW <= 1e-3) ? "PASS" : "FAIL");

  // ---- callback-level breakpoint analysis ------------------------------------------
  // From the callback source: alpha = (tau(hi) - tau(lo)) / (hi - lo) with
  // hi = min(hmax, h+25), lo = max(0, h-25), and tau piecewise LINEAR in h between the
  // profile's native knots. Hence alpha is piecewise linear with breakpoints exactly at
  // {knot +/- 25 m} plus the boundary set {0, 25, hmax-25, hmax}. Verified numerically:
  // alpha must be collinear between consecutive predicted breakpoints.
  std::printf("\n--- callback breakpoints ---\n");
  {
    double const hmax = g_ab.getHeightRange().second;
    std::vector<double> bp;
    for (double k : g_knots) {
      bp.push_back(k - 25.0);
      bp.push_back(k + 25.0);
    }
    bp.push_back(0.0);
    bp.push_back(25.0);
    bp.push_back(hmax - 25.0);
    bp.push_back(hmax);
    std::sort(bp.begin(), bp.end());
    bp.erase(std::unique(bp.begin(), bp.end()), bp.end());
    double worstLin = 0.0, argLin = 0.0;
    for (size_t i = 0; i + 1 < bp.size(); ++i) {
      double const a = bp[i], b = bp[i + 1];
      if (b - a < 1.0) continue;
      double const ya = alphaOf(690.0, a), yb = alphaOf(690.0, b);
      for (int q = 1; q <= 3; ++q) {
        double const h = a + (b - a) * double(q) / 4.0;
        double const lin = ya + (yb - ya) * (h - a) / (b - a);
        double const res = std::abs(alphaOf(690.0, h) - lin);
        double const scale = std::max(1e-30, std::abs(ya) + std::abs(yb));
        if (res / scale > worstLin) {
          worstLin = res / scale;
          argLin = h;
        }
      }
    }
    std::printf("predicted breakpoint set: %zu points from %zu profile knots (+/-25 m)\n", bp.size(),
                g_knots.size());
    std::printf("alpha piecewise LINEAR between them: worst relative residual %.3e at h=%.1f m\n",
                worstLin, argLin);
  }

  // ---- admissibility ---------------------------------------------------------------
  std::printf("\n--- admissibility of the ground-bound route ---\n");
  {
    double maxG = 0.0, argG = 0.0;
    for (double h = 0.0; h <= H_TOP; h += 10.0) {
      double const r = R + h;
      double const N = nPhase(ri, r, 1.0) - 1.0;
      double const dNdh = (nPhase(ri, r + 1.0, 1.0) - nPhase(ri, r - 1.0, 1.0)) / 2.0;
      double const G = r * std::abs(dNdh) - N;
      if (G > maxG) {
        maxG = G;
        argG = h;
      }
    }
    std::printf("n*r increasing  <=>  f * G < 1 with G = r|dN/dh| - N\n");
    std::printf("  max G = %.6f at h=%.0f m  -> holds for ALL f < %.3f\n", maxG, argG, 1.0 / maxG);
    std::printf("  dispersion spans f in [%.4f, %.4f]  ->  %s\n", Fdisp(700.0), Fdisp(300.0),
                (Fdisp(300.0) < 1.0 / maxG) ? "ADMISSIBLE" : "VIOLATED");
    double margin = 1e300;
    for (double h = 0.0; h + 25.0 <= H_TOP; h += 25.0) {
      double lo = 1e300;
      for (double z = h; z <= h + 25.0 + 1e-9; z += 2.5)
        lo = std::min(lo, nPhase(ri, R + z, 1.0) * (R + z));
      margin = std::min(margin, lo - nPhase(ri, R, 1.0) * R);
    }
    std::printf("  conservative interval-minimum margin over 25 m sub-intervals: %.4e m (>= 0 required)\n",
                margin);
  }

  // ---- performance -----------------------------------------------------------------
  std::printf("\n--- performance ---\n");
  {
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> uh(0.0, 40000.0), ul(300.0, 700.0), uu(0.0, 1.0);
    long const NPH = 20000;
    tab.resetCounters();
    long hit = 0, other = 0;
    double acc = 0.0;
    auto const t0 = std::chrono::steady_clock::now();
    for (long k = 0; k < NPH; ++k) {
      double const h_e = uh(rng);
      double const lam = ul(rng);
      double const f = Fdisp(lam);
      double const c = -uu(rng);
      double const s = std::sqrt(std::max(0.0, 1.0 - c * c));
      double const ph = 6.283185307179586 * uu(rng);
      CurvedRayTable::Result const r =
          tab.trace(0.0, 0.0, R + h_e, s * std::cos(ph), s * std::sin(ph), c, f, f);
      if (r.status == CurvedRayTable::Status::Hit) {
        ++hit;
        double const mu = -std::cos(r.theta_obs);
        double const xt = std::sqrt(std::max(0.0, 1.0 - mu * mu));
        acc += AW.lookup(h_e, xt, lam);  // the absorption lookup the caller performs
      } else {
        ++other;
      }
    }
    auto const t1 = std::chrono::steady_clock::now();
    double const usPer = std::chrono::duration<double, std::micro>(t1 - t0).count() / double(NPH);
    std::printf("combined geometry trace + absorption lookup: %.4f us/photon (%ld photons, acc=%.3e)\n",
                usPer, NPH, acc);
    std::printf("  hit %ld, rejected %ld (%.3f)  [upgoing %zu beyondCrit %zu dispersion %zu grazing %zu]\n",
                hit, other, double(other) / double(NPH), tab.rejectedUpgoing(),
                tab.rejectedBeyondCrit(), tab.rejectedDispersion(), tab.rejectedGrazing());
    std::printf("  NOTE: directions sampled isotropically downgoing; a real Cherenkov cone is\n"
                "        narrower, so this reject fraction is an upper bound.\n");

    long const NL = 2000;
    double acc2 = 0.0;
    auto const t2 = std::chrono::steady_clock::now();
    for (long k = 0; k < NL; ++k) {
      double const h_e = uh(rng);
      double const lam = ul(rng);
      double const f = Fdisp(lam);
      double const p = 0.5 * bCrit(ri, f);
      double z = R;
      while (z < R + h_e - 1e-12) {
        double const st = std::min(50.0, R + h_e - z);
        double const zm = z + 0.5 * st;
        double const s = p / (nPhase(ri, zm, f) * zm);
        if (s >= 1.0) break;
        acc2 += alphaOf(lam, zm - R) / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
        z += st;
      }
    }
    auto const t3 = std::chrono::steady_clock::now();
    double const usLegacy = std::chrono::duration<double, std::micro>(t3 - t2).count() / double(NL);
    std::printf("legacy per-photon tau integration (50 m sub-steps): %.2f us/photon (acc=%.3e)\n",
                usLegacy, acc2);
    std::printf("construction amortisation: %.1f s / (%.2f - %.4f) us = %.3g photons to break even\n",
                sU, usLegacy, usPer, sU * 1e6 / std::max(1e-9, usLegacy - usPer));
    std::printf("peak RSS at exit (geometry + both absorption tables + reference): %.1f MB\n",
                hwm_kB() / 1024.0);
  }
  return 0;
}
