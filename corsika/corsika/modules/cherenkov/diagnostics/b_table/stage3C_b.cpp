// stage3C_b.cpp -- Stage 3, Step C: localise the geometry table's path-length error.
//
// Single ray: h = 46000 m, lambda = 400 nm (f_phase = f_group = Fdisp(400) = 1 exactly),
// xi = 0.9999  =>  mu = sqrt(1 - xi^2) = 0.0141421 (the table's own mu convention,
// see CurvedRayTable.hpp:682).  One table, f = 1, so the f-axis contributes nothing
// (the Lagrange weights all carry a factor t = (fPhase-1)/eps and vanish at t = 0).
//
// Prints, all from the table ACTUALLY BUILT:
//   1. build parameters: muMin, muBreak, sTot, deltaS, fLo, fHi, hNodes, memory
//   2. the h-node drift at h = 46000 (the one open condition on "(b'') closed")
//   3. the s-cell (j, t) for the ray's mu and the two bracketing node mu values
//   4. node value vs converged reference at those two nodes      -> (b)  node/integration
//   5. interpolated value vs converged reference                 -> (b') s-interpolation
//   6. the per-node error profile up the whole column            -> prediction 2
//   7. a mu-threshold scan over the stored table                 -> prediction 3
//
// The reference integrator (lenRef) and the ray construction are IDENTICAL to
// path_b.cpp, so these numbers are directly comparable with
// log_path_b_trajectory_crosscheck.txt.
//
// Build (see README.md):  g++ -O2 -std=c++17 -I<parent of corsika/> stage3C_b.cpp -o stage3C_b
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

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

// independent geometric path along the ray (byte-identical to path_b.cpp), plus a
// flag telling whether the ray turned before reaching hTop (s >= 1).
static double lenRef(RefractiveIndexTable const& ri, double p, double f, double hTop, double dr,
                     bool* turned = nullptr) {
  if (turned) *turned = false;
  double acc = 0.0, z = R;
  double const zTop = R + hTop;
  while (z < zTop - 1e-12) {
    double const st = std::min(dr, zTop - z);
    double const zm = z + 0.5 * st;
    double const s = p / (nPhase(ri, zm, f) * zm);
    if (s >= 1.0) {
      if (turned) *turned = true;
      return acc;
    }
    acc += 1.0 / std::sqrt(std::max(1e-300, 1.0 - s * s)) * st;
    z += st;
  }
  return acc;
}

// ray construction identical to path_b.cpp.  p = xi * bCrit, and the table's own
// mu = sqrt(1 - (p/pCrit)^2) = sqrt(1 - xi^2) (CurvedRayTable.hpp:682).
struct Ray {
  double p, sinT, cosT, mu;
};
static Ray makeRay(RefractiveIndexTable const& ri, double hE, double f, double xi) {
  double const re = R + hE;
  double const pCrit = nPhase(ri, R, f) * R;
  double const p = xi * pCrit;
  double const nE = nPhase(ri, re, f);
  double const sinT = p / (nE * re);
  double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
  double const mu = std::sqrt(std::max(0.0, 1.0 - (p / pCrit) * (p / pCrit)));
  return {p, sinT, cosT, mu};
}

// The build's own substep rule (CurvedRayTable.hpp:589-601), replicated so the
// diagnostic can print the nSubBase the build actually used for each column.
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

int main() {
  RefractiveIndexTable ri = makeRI();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);

  std::printf("=== Stage 3, Step C: path-length error localisation ===\n");
  std::printf("=== one ray, one table, f = 1 exactly ===\n\n");

  // ---- 1. build parameters from the table actually built ------------------
  std::printf("1. build parameters (from the table actually built)\n");
  std::printf("   valid=%d  hNodes=%zu  muNodes=%zu  memory=%.1f MB\n", int(tab.valid()),
              tab.hNodes(), tab.muNodes(), tab.memoryMB());
  std::printf("   muMin    = %.6e\n", tab.muMin());
  std::printf("   muBreak  = %.9f\n", tab.muBreak());
  std::printf("   sTot     = %.9f\n", tab.sTot());
  std::printf("   deltaS   = %.12e\n", tab.deltaS());
  std::printf("   fLo, fHi = %.6f, %.6f\n", tab.fLo(), tab.fHi());
  std::printf("   buildFailCol=%d  buildFailReason=%d\n", tab.buildFailCol(),
              tab.buildFailReason());
  if (!tab.valid()) {
    std::printf("TABLE INVALID -- aborting\n");
    return 1;
  }

  // ---- the ray -------------------------------------------------------------
  double const hE = 46000.0, lam = 400.0, xi = 0.9999;
  double const f = Fdisp(lam);
  Ray const ry = makeRay(ri, hE, f, xi);
  double const mu = ry.mu;
  std::printf("\n   ray: h=%.0f m  lambda=%.0f nm  f=Fdisp(%.0f)=%.9f  xi=%.6f\n", hE, lam, lam, f,
              xi);
  std::printf("        mu = sqrt(1-xi^2) = %.9f\n", mu);
  std::printf("        nSubBase the build used for this column = %d\n", nSubOf(ri, mu, f));

  // ---- 2. h axis and the node drift ---------------------------------------
  std::size_t const iH = tab.hCell(hE);
  std::printf("\n2. h axis at h = %.0f m\n", hE);
  std::printf("   hCell(h)   = %zu\n", iH);
  std::printf("   hNode(i)   = %.12f    drift = %+.3e m\n", tab.hNode(iH), tab.hNode(iH) - hE);
  std::printf("   hNode(i+1) = %.12f    offset= %+.3e m\n", tab.hNode(iH + 1),
              tab.hNode(iH + 1) - hE);

  // ---- 3. s cell -----------------------------------------------------------
  std::size_t j = 0;
  double t = 0.0;
  tab.muCell(mu, j, t);
  std::printf("\n3. s axis at mu = %.9f\n", mu);
  std::printf("   cell j = %zu   local position t = %.9f\n", j, t);
  std::printf("   nodeS(j)=%.9f  nodeS(j+1)=%.9f\n", tab.nodeSOf(j), tab.nodeSOf(j + 1));
  std::printf("   muNode(j)=%.9e  muNode(j+1)=%.9e\n", tab.muNode(j), tab.muNode(j + 1));
  std::printf("   mu is %s a node\n",
              (std::fabs(t) < 1e-12 || std::fabs(t - 1.0) < 1e-12) ? "ON" : "BETWEEN");

  // ---- 4. traced (interpolated) value vs the converged reference ----------
  double const re = R + hE;
  CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, ry.sinT, 0.0, -ry.cosT, f, f);
  std::printf("\n4. trace() on this ray:  status=%d  hit=%d\n", int(r.status), int(r.hit));
  if (r.status != CurvedRayTable::Status::Hit) {
    std::printf("   NOT A HIT -- aborting Step C\n");
    return 1;
  }
  double const vacL = tab.vacuumL(hE, mu);
  double const DL_trace = r.geometric_path_m - vacL;
  std::printf("   geometric_path_m  = %.6f m\n", r.geometric_path_m);
  std::printf("   vacuumL(h, mu)    = %.6f m\n", vacL);
  std::printf("   DL (interpolated) = %.6f m\n", DL_trace);
  std::printf("   sigma_m           = %.6f m\n", r.sigma_m);
  std::printf("   optical_path_m    = %.6f m\n", r.optical_path_m);
  std::printf("   time_s            = %.12e s  = %.4f ns\n", r.time_s, r.time_s * 1e9);

  bool turned = false;
  double const Lref01 = lenRef(ri, ry.p, f, hE, 0.1);
  double const Lref001 = lenRef(ri, ry.p, f, hE, 0.01, &turned);
  std::printf("\n   independent Bouguer reference (same ray, same f, turned=%d)\n", int(turned));
  std::printf("   L_ref(dr=0.1)  = %.6f m   DL_ref = %.6f m\n", Lref01, Lref01 - vacL);
  std::printf("   L_ref(dr=0.01) = %.6f m   DL_ref = %.6f m\n", Lref001, Lref001 - vacL);
  std::printf("   reference self-convergence |dr0.01-dr0.1| = %.3e m\n", std::fabs(Lref001 - Lref01));
  std::printf("   err_interp = DL_interp - DL_ref(dr=0.01) = %+.6f m\n", DL_trace - (Lref001 - vacL));

  // ---- 5. node values vs the reference at the two bracketing node mu ------
  std::printf("\n5. node comparison (separates node-value error from interpolation error)\n");
  for (int which = 0; which < 2; ++which) {
    std::size_t const jn = j + std::size_t(which);
    if (jn >= tab.muNodes()) continue;
    double const muN = tab.muNode(jn);
    double const xiN = std::sqrt(std::max(0.0, 1.0 - muN * muN));
    Ray const rn = makeRay(ri, hE, f, xiN);
    double const vacN = tab.vacuumL(hE, muN);
    double const LrN = lenRef(ri, rn.p, f, hE, 0.1);
    double const DLrefN = LrN - vacN;
    double const storedN = tab.storedDL(iH, jn);
    std::printf("   --- node j=%zu : mu=%.9e  nSubBase=%d ---\n", jn, muN, nSubOf(ri, muN, f));
    std::printf("       storedDL(i,j) = %+.6f m\n", storedN);
    std::printf("       DL_ref        = %+.6f m\n", DLrefN);
    std::printf("       node error    = %+.6f m   (stored - ref)\n", storedN - DLrefN);
    std::printf("       storedDLh(i,j)= %+.6e\n", tab.storedDLh(iH, jn));
  }

  // ---- 6. per-node error profile up the column (prediction 2) -------------
  double const muP = tab.muNode(j);
  std::printf("\n6. per-node error profile at mu = muNode(j) = %.9e (prediction 2)\n", muP);
  std::printf("   reference evaluated every 10th node at dr = 1.0 m\n");
  std::printf("   %-6s %-15s %-18s %-18s %s\n", "i", "hNode[m]", "storedDL[m]", "DL_ref[m]",
              "diff[m]");
  double const xiP = std::sqrt(std::max(0.0, 1.0 - muP * muP));
  Ray const rp = makeRay(ri, hE, f, xiP);
  double maxAbs = 0.0, hAtMax = 0.0;
  std::size_t iAtMax = 0, nRows = 0;
  for (std::size_t i = 0; i < tab.hNodes(); ++i) {
    if (!((i % 10 == 0) || (i + 1 == tab.hNodes()))) continue;
    double const hh = tab.hNode(i);
    double const sd = tab.storedDL(i, j);
    bool tk = false;
    double const vacI = tab.vacuumL(hh, muP);
    double const LrI = lenRef(ri, rp.p, f, hh, 1.0, &tk);
    double const diff = sd - (LrI - vacI);
    if (std::fabs(diff) > maxAbs) {
      maxAbs = std::fabs(diff);
      iAtMax = i;
      hAtMax = hh;
    }
    std::printf("   %-6zu %-15.3f %+-18.6f %+-18.6f %+.6f%s\n", i, hh, sd, LrI - vacI, diff,
                tk ? "  (ref turned early)" : "");
    ++nRows;
  }
  std::printf("   -> max |node error| = %.6f m at i=%zu (h=%.1f m); %zu rows printed\n", maxAbs,
              iAtMax, hAtMax, nRows);

  // ---- 7. mu-threshold scan on the stored table (prediction 3) ------------
  std::printf("\n7. mu-threshold scan on the STORED table at h = 1000 m (prediction 3)\n");
  std::printf("   %-6s %-15s %-10s %-18s %-18s %s\n", "j", "muNode", "nSubBase", "storedDL[m]",
              "DL_ref[m]", "diff[m]");
  double const hScan = 1000.0;
  std::size_t const iScan = tab.hCell(hScan);
  for (std::size_t jj = 0; jj + 1 < tab.muNodes(); ++jj) {
    double const muJ = tab.muNode(jj);
    if (muJ > 2.5e-3) continue;
    if (muJ < tab.muMin()) break;
    double const xiJ = std::sqrt(std::max(0.0, 1.0 - muJ * muJ));
    Ray const rj = makeRay(ri, hScan, f, xiJ);
    bool tk = false;
    double const vacJ = tab.vacuumL(hScan, muJ);
    double const LrJ = lenRef(ri, rj.p, f, hScan, 0.1, &tk);
    double const sdJ = tab.storedDL(iScan, jj);
    std::printf("   %-6zu %-15.9e %-10d %+-18.6f %+-18.6f %+.6f%s\n", jj, muJ, nSubOf(ri, muJ, f),
                sdJ, LrJ - vacJ, sdJ - (LrJ - vacJ), tk ? "  (ref turned early)" : "");
  }
  std::printf("\n=== end of Step C ===\n");
  return 0;
}
