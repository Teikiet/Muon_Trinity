// stage9_abs_b.cpp -- Task E: absorption accuracy along the corrected curved ray.
//
// Step 0 identifications (all cited in README Section 22):
//   * stored quantity   : cumulative optical depth tau(lambda,h), AtmosphericAbsorption.hpp:41
//                         (ecoeff_, [wavelength_idx][height_idx]); looked up by
//                         getOpticalDepth (line 107, bilinear 140-154).
//   * production path   : extinctionCoeffFunc_ in Cherenkov.inl:96-113 -- the local
//                         extinction is a +/-(dh/2) finite difference of getOpticalDepth
//                         with dh = 50 m, returning max(0, (tau_hi - tau_lo)/actual_dh).
//                         This program calls that same rule through the production class.
//   * Sigma in CurvedRayTable.hpp is NOT absorption: it is int N ds, the refractivity
//                         column used for timing (lines 253-257, RHS ds = Nz/cc).
//
// The scan compares, along the SAME corrected curved path:
//   tau_prod = int alpha_prod(h) ds      (production rule above)
//   tau_ref  = int alpha_ref(h)  ds      (exact slope of the profile's native tau knots)
// and reports |dtau| against the 1e-3 gate of spec section 7.1.
//
// Usage: ./stage9_abs_b <absorption profile .dat>
// Build: g++ -O2 -std=c++17 -I<parent of corsika/> stage9_abs_b.cpp -o stage9_abs_b
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

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

// ---- the profile's own native knots (read independently, for alpha_ref) -----
struct Profile {
  std::vector<double> lam, hgt;                 // native knots
  std::vector<std::vector<double> > tau;        // [ilam][ihgt], cumulative from ground up
};
static bool readProfile(std::string const& fn, Profile& p) {
  std::ifstream f(fn.c_str());
  if (!f) return false;
  int nl = 0, nh = 0;
  f >> nl >> nh;
  if (nl <= 0 || nh <= 0) return false;
  p.lam.resize(nl);
  p.hgt.resize(nh);
  for (int i = 0; i < nl; ++i) f >> p.lam[i];
  for (int i = 0; i < nh; ++i) f >> p.hgt[i];
  p.tau.assign(nl, std::vector<double>(nh, 0.0));
  for (int i = 0; i < nl; ++i)
    for (int j = 0; j < nh; ++j) f >> p.tau[i][j];
  return true;
}
// exact slope of the native piecewise-linear tau(h), no interpolation smoothing
static double alphaRef(Profile const& p, int il, double h) {
  std::vector<double> const& H = p.hgt;
  std::vector<double> const& T = p.tau[il];
  if (h <= H.front()) return (T[1] - T[0]) / (H[1] - H[0]);
  if (h >= H.back()) return (T[H.size() - 1] - T[H.size() - 2]) / (H.back() - H[H.size() - 2]);
  std::size_t j = 0;
  while (j + 2 < H.size() && H[j + 1] <= h) ++j;
  return (T[j + 1] - T[j]) / (H[j + 1] - H[j]);
}
static int nearestLam(Profile const& p, double lam) {
  int best = 0;
  double bd = 1e30;
  for (std::size_t i = 0; i < p.lam.size(); ++i) {
    double const d = std::fabs(p.lam[i] - lam);
    if (d < bd) {
      bd = d;
      best = int(i);
    }
  }
  return best;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <absorption profile .dat>\n", argv[0]);
    return 2;
  }
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[1]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  Profile prof;
  if (!readProfile(argv[1], prof)) {
    std::printf("STOP: cannot re-read the profile's native knots\n");
    return 1;
  }
  auto wr = abs.getWavelengthRange();
  auto hr = abs.getHeightRange();
  std::printf("=== Task E: absorption accuracy on the corrected curved ray ===\n\n");
  std::printf("profile: %s\n", argv[1]);
  std::printf("  table  lambda range [%.1f, %.1f] nm, h range [%.1f, %.1f] m\n", wr.first, wr.second,
              hr.first, hr.second);
  std::printf("  native knots: %zu lambda, %zu h (h spacing %.1f m)\n", prof.lam.size(),
              prof.hgt.size(), prof.hgt.size() > 1 ? prof.hgt[1] - prof.hgt[0] : 0.0);

  RefractiveIndexTable ri = makeRI();
  auto const t0 = std::chrono::steady_clock::now();
  CurvedRayTable tab(ri, R, R, H_TOP, 800, 10.0, 100.0, 30000.0);
  auto const t1 = std::chrono::steady_clock::now();
  double const buildS = std::chrono::duration<double>(t1 - t0).count();
  std::printf("geometry table build: %.1f s (valid=%d)\n\n", buildS, int(tab.valid()));
  if (!tab.valid()) return 1;

  std::vector<double> lams; // the 121 nodes of section 12.2: 300..900 nm at 5 nm
  for (double l = 300.0; l <= 900.0 + 1e-9; l += 5.0) lams.push_back(l);
  double const mus[8] = {2.0e-4, 1.0e-3, 5.0e-3, 2.0e-2, 1.0e-1, 0.3, 0.6, 1.0};
  double const hs[6] = {100.0, 1000.0, 5000.0, 20000.0, 46000.0, 59900.0};
  double const dr = 1.0;

  int tot[8][6] = {{0}}, fail[8][6] = {{0}}, excl[8][6] = {{0}};
  double worst[8][6];
  for (int a = 0; a < 8; ++a)
    for (int b = 0; b < 6; ++b) worst[a][b] = 0.0;
  double gWorst = 0.0, gLam = 0, gMu = 0, gH = 0;
  double refDiffWorst = 0.0;
  int nRay = 0, nHit = 0, nExcl = 0;
  double tLookup = 0.0, tExact = 0.0;
  double worstRec[10] = {0.0}, worstRecLam[10] = {0.0}, worstRecMu[10] = {0.0}, worstRecH[10] = {0.0};

  for (std::size_t il = 0; il < lams.size(); ++il) {
    double const lam = lams[il];
    int const ilp = nearestLam(prof, lam);
    double const f = Fdisp(lam);
    for (int a = 0; a < 8; ++a) {
      double const mu = mus[a];
      double const xi = std::sqrt(std::max(0.0, 1.0 - mu * mu));
      double const p = xi * bCrit(ri, f);
      for (int b = 0; b < 6; ++b) {
        double const hE = hs[b];
        ++nRay;
        double const re = R + hE;
        double const nE = nPhase(ri, re, f);
        double const sinT = p / (nE * re);
        if (sinT >= 1.0) {
          ++excl[a][b];
          ++nExcl;
          continue;
        }
        double const cosT = std::sqrt(std::max(0.0, 1.0 - sinT * sinT));
        CurvedRayTable::Result const r = tab.trace(0.0, 0.0, re, sinT, 0.0, -cosT, f, f);
        if (r.status != CurvedRayTable::Status::Hit) {
          ++excl[a][b]; // mu=1.0 degenerate (line 780) or the grazing endpoint
          ++nExcl;
          continue;
        }
        ++nHit;
        double tauProdW = 0.0, tauRefW = 0.0, tauProdZ = 0.0, tauRefZ = 0.0;
        { // (a) w = sqrt(z-R), the substituted integral of sections 6.1 / 9.2
          double const wTop = std::sqrt(hE);
          for (double w = 0.0; w < wTop - 1e-12;) {
            double const st = std::min(dr, wTop - w);
            double const wn = w + st;
            if (!(wn > w)) break;
            double const wm = w + 0.5 * st;
            double const z = R + wm * wm;
            double const h = wm * wm;
            double const s = p / (nPhase(ri, z, f) * z);
            if (s >= 1.0) break;
            double const ds = (2.0 * wm) * st / std::sqrt(std::max(1e-300, 1.0 - s * s));
            double const hl = std::max(hr.first, h - 25.0), hh = std::min(hr.second, h + 25.0);
            double const adh = hh - hl;
            auto const q0 = std::chrono::steady_clock::now();
            double ap = 0.0;
            if (adh > 0.0)
              ap = std::max(0.0,
                            (abs.getOpticalDepth(lam, hh) - abs.getOpticalDepth(lam, hl)) / adh);
            auto const q1 = std::chrono::steady_clock::now();
            tLookup += std::chrono::duration<double>(q1 - q0).count();
            auto const q2 = std::chrono::steady_clock::now();
            double const ar = alphaRef(prof, ilp, h);
            auto const q3 = std::chrono::steady_clock::now();
            tExact += std::chrono::duration<double>(q3 - q2).count();
            tauProdW += ap * ds;
            tauRefW += ar * ds;
            w = wn;
          }
        }
        { // (b) z-march, the independent Bouguer parameterisation of section 6.5
          double z = R;
          double const zTop = R + hE;
          while (z < zTop - 1e-12) {
            double const st = std::min(dr, zTop - z);
            double const zn = z + st;
            if (!(zn > z)) break;
            double const zm = z + 0.5 * st;
            double const h = zm - R;
            double const s = p / (nPhase(ri, zm, f) * zm);
            if (s >= 1.0) break;
            double const ds = st / std::sqrt(std::max(1e-300, 1.0 - s * s));
            double const hl = std::max(hr.first, h - 25.0), hh = std::min(hr.second, h + 25.0);
            double const adh = hh - hl;
            double ap = 0.0;
            if (adh > 0.0)
              ap = std::max(0.0,
                            (abs.getOpticalDepth(lam, hh) - abs.getOpticalDepth(lam, hl)) / adh);
            tauProdZ += ap * ds;
            tauRefZ += alphaRef(prof, ilp, h) * ds;
            z = zn;
          }
        }
        double const dW = std::fabs(tauProdW - tauRefW);
        double const dZ = std::fabs(tauProdZ - tauRefZ);
        double const d = std::max(dW, dZ);
        double const refDiff = std::fabs(tauRefW - tauRefZ);
        if (refDiff > refDiffWorst) refDiffWorst = refDiff;
        ++tot[a][b];
        if (d > 1e-3) ++fail[a][b];
        if (d > worst[a][b]) worst[a][b] = d;
        if (d > gWorst) {
          gWorst = d;
          gLam = lam;
          gMu = mu;
          gH = hE;
        }
        if (d > worstRec[9]) {
          int q = 9;
          while (q > 0 && d > worstRec[q - 1]) {
            worstRec[q] = worstRec[q - 1];
            worstRecLam[q] = worstRecLam[q - 1];
            worstRecMu[q] = worstRecMu[q - 1];
            worstRecH[q] = worstRecH[q - 1];
            --q;
          }
          worstRec[q] = d;
          worstRecLam[q] = lam;
          worstRecMu[q] = mu;
          worstRecH[q] = hE;
        }
      }
    }
  }

  // ---- report --------------------------------------------------------------
  char const* ml[8] = {"2e-4", "1e-3", "5e-3", "2e-2", "1e-1", "0.3", "0.6", "1.0"};
  std::printf("fails / total (worst |dtau|) per mu x h; E = excluded (rejection/fallback)\n");
  std::printf("%-9s", "mu \\ h");
  for (int b = 0; b < 6; ++b) std::printf(" %-16s", "");
  std::printf("\n");
  std::printf("%-9s", "");
  for (int b = 0; b < 6; ++b) std::printf(" h=%-14.0f", hs[b]);
  std::printf("\n");
  for (int a = 0; a < 8; ++a) {
    std::printf("%-9s", ml[a]);
    for (int b = 0; b < 6; ++b) {
      char buf[32];
      if (tot[a][b] == 0) std::snprintf(buf, sizeof(buf), "E(%d)", excl[a][b]);
      else std::snprintf(buf, sizeof(buf), "%d/%-2d w=%.1e", fail[a][b], tot[a][b], worst[a][b]);
      std::printf(" %-16s", buf);
    }
    std::printf("\n");
  }
  std::printf("\nrays %d, hits %d, excluded %d (fraction %.3f)\n", nRay, nHit, nExcl,
              nRay ? double(nExcl) / double(nRay) : 0.0);
  std::printf("worst |dtau| = %.4e at lambda=%.0f nm, mu=%.2e, h=%.0f m   (gate 1e-3)\n", gWorst,
              gLam, gMu, gH);
  std::printf("archived value on file: 5.775e-07 (a REFERENCE resolution study, not a table error)\n");
  std::printf("worst disagreement between the two references (tau_ref in w vs in z) = %.4e\n",
              refDiffWorst);
  std::printf("  -> if this exceeds 1e-3 the two references disagree by more than the gate\n");
  std::printf("\n10 worst rays:\n");
  for (int q = 0; q < 10; ++q)
    std::printf("   %2d. %.4e  lambda=%.0f mu=%.2e h=%.0f\n", q + 1, worstRec[q], worstRecLam[q],
                worstRecMu[q], worstRecH[q]);
  std::printf("\ncost: production lookup %.3f us/step, exact slope %.3f us/step (ratio %.2f)\n",
              nHit ? tLookup / double(nHit) * 1e6 : 0.0, nHit ? tExact / double(nHit) * 1e6 : 0.0,
              tExact > 0.0 ? tLookup / tExact : 0.0);
  std::printf("geometry table build %.1f s\n", buildS);
  std::printf("\n=== end of Task E ===\n");
  return 0;
}
