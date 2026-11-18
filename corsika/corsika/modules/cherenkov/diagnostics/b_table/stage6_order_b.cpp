// stage6_order_b.cpp -- Stage 6, Task 1: WHY is the column integration first order?
//
// Step C' measured p = 1.00 (nSub x4 -> error /4.0) on the Stage 3 failing ray.
// This program runs the SAME local column copy in four variants at nSub x1/x4/x16/x64
// and reports the error against a matching reference and the observed order.
//
//   V0  baseline            table N and dN/dz, exactly as the build uses them
//   V1  smooth atmosphere   analytic N(h) = a1 exp(-h/H1) + a2 exp(-h/H2) and its exact
//                           derivative -- the same two-term profile the table samples at
//                           100 m and stores as a piecewise-exponential interpolant.
//                           Reference integrated on the SAME analytic profile.
//   V2  knot-aligned steps  table N, but each h cell is split at the refractive-index
//                           knots so no substep spans one; V0 also reports how many
//                           substeps crossed a knot.
//   V3  w-substitution      table N, variable changed to w = sqrt(h), substeps uniform
//                           in w, so the near-ground stiffness is resolved.
//
// Diagnostics only.  buildColumn is private and has no override, so the copy lives here.
// No production source is touched.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage6_order_b.cpp -o stage6_order_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0, C_LIGHT = 2.99792458e8;

// the profile the table samples (100 m knots, log-linear interpolation)
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
static double bCrit(RefractiveIndexTable const& ri, double f) { return (1.0 + ri.N(R) * f) * R; }

// ---- atmosphere accessors, switchable ---------------------------------------
// atmo = 0 : the table (piecewise exponential, slope jumps at every 100 m knot)
// atmo = 1 : the analytic two-term profile and its exact derivative (smooth)
static constexpr double A1 = 2.20e-4, H1 = 7500.0, A2 = 0.55e-4, H2 = 1600.0;
static double Nof(RefractiveIndexTable const& ri, int atmo, double hh) {
  if (atmo == 1) return A1 * std::exp(-hh / H1) + A2 * std::exp(-hh / H2);
  return ri.N(R + hh);
}
static double dNof(RefractiveIndexTable const& ri, int atmo, double hh) {
  if (atmo == 1) return -(A1 / H1) * std::exp(-hh / H1) - (A2 / H2) * std::exp(-hh / H2);
  double const z = R + hh;
  return ri.logSlope(ri.cellIndex(z)) * ri.N(z);
}

// ---- reference for L and Sigma on the SAME atmosphere -----------------------
// (with the stall guard required by the brief; the archived bare form can hang)
static void refLS(RefractiveIndexTable const& ri, int atmo, double p, double f, double hTop,
                  double dr, double& L, double& S, bool* turned) {
  L = 0.0;
  S = 0.0;
  if (turned) *turned = false;
  double z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zn = z + st;
    if (!(zn > z)) break; // stall guard
    double const zm = z + 0.5 * st;
    double const nz = 1.0 + Nof(ri, atmo, zm - R) * f;
    double const s = p / (nz * zm);
    if (s >= 1.0) {
      if (turned) *turned = true;
      return;
    }
    double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s));
    L += inv * st;
    S += Nof(ri, atmo, zm - R) * inv * st;
    z = zn;
  }
}

static int nSubOf(RefractiveIndexTable const& ri, double mu, double f) {
  int nSubBase = 4;
  if (mu < 1.0e-3) {
    double const r0 = R;
    double const Nv = ri.N(r0);
    double const dN = ri.logSlope(ri.cellIndex(r0)) * Nv;
    double const a1 = 1.0 / r0 + dN / (1.0 + Nv);
    double const af = 1.0 / r0 + (f * dN) / (1.0 + f * Nv);
    double const stiff = (af > 0.0 && a1 > 0.0) ? std::max(1.0, a1 / af) : 1.0;
    double const k = std::min(256.0, std::ceil(1.0e-3 / mu) * stiff);
    nSubBase = static_cast<int>(4.0 * k);
  }
  return nSubBase;
}

static std::vector<double> makeZNodes(double hTop, double zBreak, double dhLo, double dhHi) {
  std::vector<double> z;
  z.push_back(0.0);
  for (double h = 0.0;;) {
    double const step = (h < zBreak) ? dhLo : dhHi;
    double nh = h + step;
    if (nh > hTop) nh = hTop;
    z.push_back(nh);
    if (nh >= hTop - 1e-9) break;
    h = nh;
  }
  return z;
}

struct Cfg {
  int atmo = 0;                 // 0 table, 1 analytic
  int splitKnots = 0;           // 1 = no substep may span a 100 m refractive-index knot
  int wsub = 0;                 // 1 = substeps uniform in w = sqrt(h)
  int fixedSlope = 0;           // 1 = use the substep's START cell slope throughout the step
  long* knotsCrossed = nullptr; // counted for atmo == 0
};

// ---- local copy of buildColumn, with the variants --------------------------
static void buildColumnVar(RefractiveIndexTable const& ri, double rObs, std::vector<double> const& z,
                           double mu, double f, int nSubBase, int firstCellBoost, Cfg const& cfg,
                           std::vector<double>& lOut) {
  std::size_t const n_h = z.size();
  lOut.assign(n_h, 0.0);
  double th = std::acos(std::min(1.0, std::max(0.0, mu)));
  double L = 0.0;
  double stepSlope = 0.0; // used when cfg.fixedSlope
  auto rhs = [&](double hh, double t, double& dth, double& dl) {
    double const r = rObs + hh;
    double const Nz = Nof(ri, cfg.atmo, hh);
    double const n = 1.0 + Nz * f;
    double const dN = cfg.fixedSlope ? (stepSlope * Nz) : dNof(ri, cfg.atmo, hh);
    double const npn = (f * dN) / n;
    double const k = 1.0 / r + npn;
    double const c = std::cos(t), s = std::sin(t);
    double const cc = (c > 1e-12) ? c : 1e-12;
    dth = -s * k / cc;
    dl = 1.0 / cc;
  };
  auto sweep = [&](double ha, double hb, int n) { // RK4 over [ha,hb] with n substeps
    double const dh = (hb - ha) / double(n);
    double hh = ha;
    for (int m = 0; m < n; ++m) {
      if (cfg.fixedSlope) stepSlope = ri.logSlope(ri.cellIndex(R + hh));
      double k1t, k1l, k2t, k2l, k3t, k3l, k4t, k4l;
      rhs(hh, th, k1t, k1l);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k1t, k2t, k2l);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k2t, k3t, k3l);
      rhs(hh + dh, th + dh * k3t, k4t, k4l);
      th += dh / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
      L += dh / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
      hh += dh;
    }
  };
  for (std::size_t i = 0; i + 1 < n_h; ++i) {
    double const hA = z[i], hB = z[i + 1];
    int const nSub = (i == 0) ? firstCellBoost * nSubBase : nSubBase;
    if (cfg.knotsCrossed && cfg.atmo == 0) {
      double const dh = (hB - hA) / double(nSub);
      for (int m = 0; m < nSub; ++m) {
        double const lo = hA + m * dh, hi = lo + dh;
        if (std::floor((lo + 1e-9) / 100.0) != std::floor((hi - 1e-9) / 100.0)) ++(*cfg.knotsCrossed);
      }
    }
    if (cfg.wsub) {
      // uniform in w = sqrt(h); dh/dw = 2w so d/dw = 2w * d/dh
      double const wa = std::sqrt(std::max(0.0, hA)), wb = std::sqrt(hB);
      double const dw = (wb - wa) / double(nSub);
      double w = wa;
      for (int m = 0; m < nSub; ++m) {
        double const jac = 2.0 * w;
        double k1t, k1l, k2t, k2l, k3t, k3l, k4t, k4l;
        rhs(w * w, th, k1t, k1l);
        double const w2 = w + 0.5 * dw;
        rhs(w2 * w2, th + 0.5 * dw * jac * k1t, k2t, k2l);
        rhs(w2 * w2, th + 0.5 * dw * jac * k2t, k3t, k3l);
        double const w3 = w + dw;
        rhs(w3 * w3, th + dw * jac * k3t, k4t, k4l);
        th += dw * jac / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
        L += dw * jac / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
        w = w3;
      }
    } else if (cfg.splitKnots) {
      double h = hA;
      while (h < hB - 1e-12) {
        double const nextKnot = (std::floor(h / 100.0) + 1.0) * 100.0;
        double const segEnd = std::min(hB, nextKnot);
        sweep(h, segEnd, nSub);
        h = segEnd;
      }
    } else {
      sweep(hA, hB, nSub);
    }
    lOut[i + 1] = L;
  }
}

struct Var {
  char const* name;
  Cfg cfg;
};

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 6, Task 1: why is the column integration first order? ===\n\n");
  std::printf("table: valid=%d hNodes=%zu memory=%.1f MB\n", int(tab.valid()), tab.hNodes(),
              tab.memoryMB());
  if (!tab.valid()) return 1;

  std::vector<double> const z = makeZNodes(H_TOP, 30000.0, 10.0, 100.0);
  double const hE = 46000.0, lam = 400.0, xi = 0.9999, f = Fdisp(lam);
  double const muRay = std::sqrt(std::max(0.0, 1.0 - xi * xi));
  std::size_t j = 0;
  double t = 0.0;
  tab.muCell(muRay, j, t);
  double const muJ = tab.muNode(j);
  std::size_t const iE = tab.hCell(hE);
  int const base = nSubOf(ri, muJ, f);
  std::printf("ray: h=%.0f lam=%.0f xi=%.6f  -> cell j=%zu (t=%.6f), node mu=%.9e\n", hE, lam, xi, j,
              t, muJ);
  std::printf("     build's own nSubBase = %d, first-cell boost x4\n\n", base);

  // ---- STOP CONDITION: V0 must reproduce the stored column -----------------
  std::vector<double> lv;
  Cfg c0;
  buildColumnVar(ri, R, z, muJ, f, base, 4, c0, lv);
  double const dl0 = lv[iE] - tab.vacuumL(hE, muJ);
  double const stored = tab.storedDL(iE, j);
  std::printf("STOP CHECK  V0 DL_copy = %.6f m   storedDL = %.6f m   diff = %.3e m\n", dl0, stored,
              dl0 - stored);
  if (std::fabs(dl0 - stored) > 0.01) {
    std::printf("V0 does NOT reproduce the stored column -- STOPPING.\n");
    return 1;
  }
  std::printf("            (reproduces; %s)\n\n", std::fabs(dl0 - stored) < 2e-3 ? "within 2e-3 m" : "within 1e-2 m");

  // ---- knot crossings in V0 -------------------------------------------------
  long cross = 0;
  Cfg ck;
  ck.knotsCrossed = &cross;
  std::vector<double> lk;
  buildColumnVar(ri, R, z, muJ, f, base, 4, ck, lk);
  std::printf("V0 substeps that straddle a 100 m refractive-index knot: %ld  (of ~%d)\n\n", cross,
              int(z.size() - 1) * base);

  // ---- references ----------------------------------------------------------
  double Lr0 = 0.0, Sr0 = 0.0, Lr1 = 0.0, Sr1 = 0.0;
  bool tn = false;
  double const pN = std::sqrt(std::max(0.0, 1.0 - muJ * muJ)) * bCrit(ri, f); // NODE's xi
  refLS(ri, 0, pN, f, hE, 0.01, Lr0, Sr0, &tn);
  refLS(ri, 1, pN, f, hE, 0.01, Lr1, Sr1, &tn);
  double const vac = tab.vacuumL(hE, muJ);
  std::printf("references at dr=0.01 m:\n");
  std::printf("   table atmosphere   DL_ref = %.6f m   Sigma_ref = %.6f m\n", Lr0 - vac, Sr0);
  std::printf("   analytic profile   DL_ref = %.6f m   Sigma_ref = %.6f m\n", Lr1 - vac, Sr1);
  std::printf("   (difference between the two references: %.6f m)\n\n", (Lr0 - Lr1));

  Var vars[5];
  vars[0].name = "V0 baseline (table)";
  vars[1].name = "V1 smooth analytic profile";
  vars[2].name = "V2 knot-aligned substeps";
  vars[3].name = "V3 w=sqrt(h) substitution";
  vars[4].name = "V4 slope frozen per substep";
  vars[1].cfg.atmo = 1;
  vars[2].cfg.splitKnots = 1;
  vars[3].cfg.wsub = 1;
  vars[4].cfg.fixedSlope = 1;

  int const scales[4] = {1, 4, 16, 64};
  std::printf("%-28s %-6s %-14s %-13s %-12s %s\n", "variant", "scale", "nSubBase", "DL_copy[m]",
              "error[m]", "p (order)");
  double bestErr = 1e30;
  int bestVar = -1, bestScale = -1;
  for (int v = 0; v < 5; ++v) {
    double e[4];
    for (int s = 0; s < 4; ++s) {
      std::vector<double> lo;
      buildColumnVar(ri, R, z, muJ, f, base * scales[s], 4, vars[v].cfg, lo);
      double const dl = lo[iE] - vac;
      double const ref = (vars[v].cfg.atmo == 1) ? (Lr1 - vac) : (Lr0 - vac);
      e[s] = dl - ref;
      std::printf("%-28s x%-5d %-14d %-13.6f %+-12.6f %s\n", vars[v].name, scales[s], base * scales[s],
                  dl, e[s], "-");
      if (s > 0 && e[s] != 0.0 && e[s - 1] != 0.0) {
        double const pp = std::log(std::fabs(e[s - 1] / e[s])) / std::log(double(scales[s]) / double(scales[s - 1]));
        std::printf("%-28s %-6s %-14s %-13s %-12s p = %.3f\n", "", "", "", "", "", pp);
      }
      if (std::fabs(e[s]) < bestErr) {
        bestErr = std::fabs(e[s]);
        bestVar = v;
        bestScale = s;
      }
    }
    std::printf("\n");
  }

  // ---- smallest nSub that reaches 0.03 m (about 0.1 ns) --------------------
  std::printf("smallest nSub reaching |error| < 0.03 m (about 0.1 ns) per variant:\n");
  for (int v = 0; v < 5; ++v) {
    double const e1 = [&] {
      std::vector<double> lo;
      buildColumnVar(ri, R, z, muJ, f, base, 4, vars[v].cfg, lo);
      double const ref = (vars[v].cfg.atmo == 1) ? (Lr1 - vac) : (Lr0 - vac);
      return lo[iE] - vac - ref;
    }();
    double const need = 0.03 / std::fabs(e1); // ~ 1/|e1| scaling if first order
    if (std::fabs(e1) < 0.03)
      std::printf("   %-28s error(x1) = %+-9.6f m  -> ALREADY below 0.03 m at x1\n", vars[v].name, e1);
    else
      std::printf("   %-28s error(x1) = %+-9.6f m  -> needs about x%.0f  (nSubBase ~ %.0f)\n",
                  vars[v].name, e1, need, base * need);
  }
  std::printf("\nbest single result: %s at x%d, |error| = %.6f m\n",
              bestVar >= 0 ? vars[bestVar].name : "-", bestScale >= 0 ? scales[bestScale] : 0,
              bestErr);
  std::printf("\n=== end of Task 1 ===\n");
  return 0;
}
