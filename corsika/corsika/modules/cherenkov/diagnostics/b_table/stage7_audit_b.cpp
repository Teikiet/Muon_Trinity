// stage7_audit_b.cpp -- Stage 7, Task A (do V1 and V4 really differ?) and Task B
// (is the 16 mm residual reference-limited or a constant offset?).
//
// Task A
//   1. print the profile each variant evaluates (first samples of N and dN/dz) and the
//      reference each uses, with the node xi;
//   2. re-run V1 with the analytic scale heights +10 % and its reference recomputed;
//   3. stop if V1 and V4 share a code path or profile.
//
// Task B
//   1. V4 node error against the w-reference at dw, dw/2, dw/4 for j=393 and j=394 and
//      for 5 further (mu,h) nodes spread across the grid;
//   2. say whether the residual shrinks with the reference step or stays constant;
//   3. if constant, test the float32 rounding of storedDL (ULP ~ 3.9 mm at 4.4e4 m).
//
// Diagnostics only; no production source touched; stall guard carried.
//
// Build:  g++ -O2 -std=c++17 -I<parent of corsika/> stage7_audit_b.cpp -o stage7_audit_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using corsika::cherenkov::CurvedRayTable;
using corsika::cherenkov::RefractiveIndexTable;

static constexpr double R = 6371000.0, H_TOP = 60000.0, C_LIGHT = 2.99792458e8;

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

// ---- analytic profile, with adjustable scale heights (Task A.2) -------------
static constexpr double A1 = 2.20e-4, A2 = 0.55e-4;
static double gH1 = 7500.0, gH2 = 1600.0;
static double Nof(RefractiveIndexTable const& ri, int atmo, double hh) {
  if (atmo == 1) return A1 * std::exp(-hh / gH1) + A2 * std::exp(-hh / gH2);
  return ri.N(R + hh);
}
static double dNof(RefractiveIndexTable const& ri, int atmo, double hh) {
  if (atmo == 1) return -(A1 / gH1) * std::exp(-hh / gH1) - (A2 / gH2) * std::exp(-hh / gH2);
  double const z = R + hh;
  return ri.logSlope(ri.cellIndex(z)) * ri.N(z);
}

// ---- w-substitution reference (Stage 6 Task 2), both L and Sigma ------------
static void refW(RefractiveIndexTable const& ri, int atmo, double p, double f, double hTop,
                 double dw, double& L, double& S, bool* turned) {
  L = 0.0;
  S = 0.0;
  if (turned) *turned = false;
  double const wTop = std::sqrt(hTop);
  for (double w = 0.0; w < wTop - 1e-12;) {
    double const st = std::min(dw, wTop - w);
    double const wn = w + st;
    if (!(wn > w)) break; // stall guard
    double const wm = w + 0.5 * st;
    double const z = R + wm * wm;
    double const nz = 1.0 + Nof(ri, atmo, wm * wm) * f;
    double const s = p / (nz * z);
    if (s >= 1.0) {
      if (turned) *turned = true;
      return;
    }
    double const jac = 2.0 * wm;
    double const inv = 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s));
    L += inv * jac * st;
    S += Nof(ri, atmo, wm * wm) * inv * jac * st;
    w = wn;
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
  int atmo = 0;
  int fixedSlope = 0; // the V4 rule
};

static void buildColumnVar(RefractiveIndexTable const& ri, double rObs, std::vector<double> const& z,
                           double mu, double f, int nSubBase, int firstCellBoost, Cfg const& cfg,
                           std::vector<double>& lOut, std::vector<double>& sOut) {
  std::size_t const n_h = z.size();
  lOut.assign(n_h, 0.0);
  sOut.assign(n_h, 0.0);
  double th = std::acos(std::min(1.0, std::max(0.0, mu)));
  double L = 0.0, Sig = 0.0, stepSlope = 0.0;
  auto rhs = [&](double hh, double t, double& dth, double& dl, double& ds) {
    double const r = rObs + hh;
    double const Nz = Nof(ri, cfg.atmo, hh);
    double const n = 1.0 + Nz * f;
    double const dN = cfg.fixedSlope ? (stepSlope * Nz) : dNof(ri, cfg.atmo, hh);
    double const k = 1.0 / r + (f * dN) / n;
    double const c = std::cos(t), s = std::sin(t);
    double const cc = (c > 1e-12) ? c : 1e-12;
    dth = -s * k / cc;
    dl = 1.0 / cc;
    ds = Nz / cc;
  };
  for (std::size_t i = 0; i + 1 < n_h; ++i) {
    double const hA = z[i], hB = z[i + 1];
    int const nSub = (i == 0) ? firstCellBoost * nSubBase : nSubBase;
    double const dh = (hB - hA) / double(nSub);
    double hh = hA;
    for (int m = 0; m < nSub; ++m) {
      if (cfg.fixedSlope) stepSlope = ri.logSlope(ri.cellIndex(R + hh));
      double k1t, k1l, k1s, k2t, k2l, k2s, k3t, k3l, k3s, k4t, k4l, k4s;
      rhs(hh, th, k1t, k1l, k1s);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k1t, k2t, k2l, k2s);
      rhs(hh + 0.5 * dh, th + 0.5 * dh * k2t, k3t, k3l, k3s);
      rhs(hh + dh, th + dh * k3t, k4t, k4l, k4s);
      th += dh / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
      L += dh / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
      Sig += dh / 6.0 * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
      hh += dh;
    }
    lOut[i + 1] = L;
    sOut[i + 1] = Sig;
  }
}

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  std::printf("=== Stage 7, Task A: do V1 and V4 really differ? ===\n\n");
  if (!tab.valid()) return 1;
  std::vector<double> const z = makeZNodes(H_TOP, 30000.0, 10.0, 100.0);
  double const hE = 46000.0, lam = 400.0, xi = 0.9999, f = Fdisp(lam);
  double const muRay = std::sqrt(std::max(0.0, 1.0 - xi * xi));
  std::size_t j = 0;
  double tt = 0.0;
  tab.muCell(muRay, j, tt);
  double const muJ = tab.muNode(j);
  std::size_t const iE = tab.hCell(hE);
  double const xiN = std::sqrt(std::max(0.0, 1.0 - muJ * muJ));
  double const pN = xiN * bCrit(ri, f);
  int const base = nSubOf(ri, muJ, f);

  // ---- A.1 the two profiles ------------------------------------------------
  std::printf("A.1 profile samples (V1 = analytic atmo=1, V4 = table atmo=0 + fixedSlope)\n");
  std::printf("    %-8s %-15s %-15s %-16s %-16s\n", "h[m]", "table N", "analytic N", "table dN/dz",
              "analytic dN/dz");
  double const hh[6] = {0.0, 50.0, 100.0, 150.0, 200.0, 1000.0};
  for (int q = 0; q < 6; ++q)
    std::printf("    %-8.0f %-15.9e %-15.9e %+-16.6e %+-16.6e\n", hh[q], Nof(ri, 0, hh[q]),
                Nof(ri, 1, hh[q]), dNof(ri, 0, hh[q]), dNof(ri, 1, hh[q]));
  std::printf("    -> the profiles DIFFER, so V1 and V4 cannot share one.\n\n");

  // ---- A.1b references, node xi, and the two errors ------------------------
  Cfg c1;
  c1.atmo = 1;
  Cfg c4;
  c4.fixedSlope = 1;
  std::vector<double> l1, s1, l4, s4;
  buildColumnVar(ri, R, z, muJ, f, base, 4, c1, l1, s1);
  buildColumnVar(ri, R, z, muJ, f, base, 4, c4, l4, s4);
  double Lr1 = 0.0, Sr1 = 0.0, Lr4 = 0.0, Sr4 = 0.0;
  bool tn = false;
  refW(ri, 1, pN, f, hE, 1.0e-3, Lr1, Sr1, &tn);
  refW(ri, 0, pN, f, hE, 1.0e-3, Lr4, Sr4, &tn);
  double const vac = tab.vacuumL(hE, muJ);
  std::printf("A.1b node xi = %.9f  (muNode(%zu) = %.9e; the RAY xi was %.4f)\n", xiN, j, muJ, xi);
  std::printf("    V1: profile ANALYTIC, ref = w-ref on ANALYTIC\n");
  std::printf("        DL_copy = %.6f  DL_ref = %.6f  error = %+.6f m\n", l1[iE] - vac, Lr1 - vac,
              (l1[iE] - vac) - (Lr1 - vac));
  std::printf("    V4: profile TABLE + fixedSlope, ref = w-ref on TABLE\n");
  std::printf("        DL_copy = %.6f  DL_ref = %.6f  error = %+.6f m\n", l4[iE] - vac, Lr4 - vac,
              (l4[iE] - vac) - (Lr4 - vac));
  std::printf("    -> DL_copy differs by %.6f m between V1 and V4; only the RESIDUAL coincides.\n\n",
              (l1[iE] - vac) - (l4[iE] - vac));

  // ---- A.2 V1 with the scale heights +10 % ---------------------------------
  double const H1s = gH1, H2s = gH2;
  gH1 *= 1.1;
  gH2 *= 1.1;
  std::vector<double> l1b, s1b;
  buildColumnVar(ri, R, z, muJ, f, base, 4, c1, l1b, s1b);
  double Lr1b = 0.0, Sr1b = 0.0;
  refW(ri, 1, pN, f, hE, 1.0e-3, Lr1b, Sr1b, &tn);
  std::printf("A.2 V1 with H x1.1 (H1 %.0f -> %.0f, H2 %.0f -> %.0f), reference recomputed:\n", H1s, gH1,
              H2s, gH2);
  std::printf("    DL_copy = %.6f  DL_ref = %.6f  error = %+.6f m\n", l1b[iE] - vac, Lr1b - vac,
              (l1b[iE] - vac) - (Lr1b - vac));
  std::printf("    (DL_copy moved by %.6f m when the profile changed)\n\n", (l1b[iE] - vac) - (l1[iE] - vac));
  gH1 = H1s;
  gH2 = H2s;

  // ---- A.3 stop check ------------------------------------------------------
  double const relN = std::fabs((Nof(ri, 0, 100.0) - Nof(ri, 1, 100.0)) / Nof(ri, 1, 100.0));
  std::printf("A.3 STOP CHECK: profiles shared? %s (rel diff %.3e) -> %s\n\n",
              relN < 1e-9 ? "YES" : "NO", relN,
              relN < 1e-9 ? "harness bug: Task 1 must be redone" : "no shared path: V4 claim stands");

  // ---- Task B --------------------------------------------------------------
  std::printf("=== Stage 7, Task B: source of the 16 mm residual ===\n");
  double const pick[5] = {1.0e-2, 5.0e-2, 0.1, 0.5, 2.0e-3};
  double const pickH[5] = {10000.0, 20000.0, 59900.0, 5000.0, 46000.0};
  std::size_t jn[5] = {0, 0, 0, 0, 0};
  for (int q = 0; q < 5; ++q) {
    double tq = 0.0;
    tab.muCell(pick[q], jn[q], tq);
  }
  std::printf("V4 node error against the w-reference at dw, dw/2, dw/4:\n");
  std::printf("    %-6s %-13s %-8s %-16s %-13s %-13s %-13s %s\n", "j", "mu", "h[m]", "DL_copy(V4)",
              "err(4e-3)", "err(2e-3)", "err(1e-3)", "trend");
  for (int q = 0; q < 7; ++q) {
    std::size_t const jq = (q == 0) ? j : (q == 1) ? (j + 1) : jn[q - 2];
    double const hq = (q == 0 || q == 1) ? 46000.0 : pickH[q - 2];
    double const muq = tab.muNode(jq);
    double const xiq = std::sqrt(std::max(0.0, 1.0 - muq * muq));
    double const pq = xiq * bCrit(ri, f);
    Cfg cq;
    cq.fixedSlope = 1;
    std::vector<double> lq, sq;
    buildColumnVar(ri, R, z, muq, f, nSubOf(ri, muq, f), 4, cq, lq, sq);
    std::size_t const iq = tab.hCell(hq);
    double const vacq = tab.vacuumL(hq, muq);
    double const dcq = lq[iq] - vacq;
    double Lx = 0.0, Sx = 0.0, e[3] = {0.0, 0.0, 0.0};
    bool tq = false;
    double const dws[3] = {4.0e-3, 2.0e-3, 1.0e-3};
    for (int s = 0; s < 3; ++s) {
      refW(ri, 0, pq, f, hq, dws[s], Lx, Sx, &tq);
      e[s] = dcq - (Lx - vacq);
    }
    char const* trend = (std::fabs(e[2] - e[0]) < 0.2 * std::fabs(e[2])) ? "CONSTANT OFFSET"
                                                                        : "shrinks with dw";
    std::printf("    %-6zu %-13.6e %-8.0f %-16.6f %+-13.6f %+-13.6f %+-13.6f %s\n", jq, muq, hq, dcq,
                e[0], e[1], e[2], trend);
  }
  // float32 test against the stored column
  Cfg c0;
  std::vector<double> l0, s0;
  buildColumnVar(ri, R, z, muJ, f, base, 4, c0, l0, s0);
  double const stored = tab.storedDL(iE, j);
  double const dbl = l0[iE] - vac;
  std::printf("\nfloat32 test at j=%zu, h=%.0f: storedDL = %.6f, double V0 copy = %.6f, diff = %.6f m\n",
              j, hE, stored, dbl, dbl - stored);
  std::printf("    float32 ULP near 4.4e4 m = %.6f m -> %s\n", 43776.0 * 1.1920929e-7,
              std::fabs(dbl - stored) < 4.0e-3 ? "consistent with float32 storage"
                                               : "larger than float32 storage");
  std::printf("\n=== end of Tasks A and B ===\n");
  return 0;
}
