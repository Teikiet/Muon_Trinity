// cp2fh2_floor_probe_b.cpp -- DIAGNOSTIC ONLY (no production edit, no gate)
// (1) grid/monotonicity of the loaded cumulative-tau profile across ALL lambdas
// (2) bounded probe of the RAW +-25 m window representation (covered-length
//     normalised, NO floor) to quantify whether max(0,.) can matter.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

int main(int argc, char** argv) {
  char const* fn = (argc > 1) ? argv[1] : "abstable_MODTRAN_bottomup.dat";
  std::ifstream f(fn);
  if (!f) { std::printf("cannot open %s\n", fn); return 2; }
  int nw = 0, nh = 0;
  f >> nw >> nh;
  std::vector<double> lam(nw), hgt(nh);
  for (int i = 0; i < nw; ++i) f >> lam[i];
  for (int j = 0; j < nh; ++j) f >> hgt[j];
  std::vector<std::vector<double> > tau(nw, std::vector<double>(nh));
  for (int i = 0; i < nw; ++i)
    for (int j = 0; j < nh; ++j) f >> tau[i][j];

  double const zmin = hgt.front(), zmax = hgt.back();
  std::printf("GRID nw=%d nh=%d | lambda %.0f..%.0f nm | height %.0f..%.0f m | dz=%.0f m\n", nw,
              nh, lam.front(), lam.back(), zmin, zmax, hgt[1] - hgt[0]);

  long nD = 0, nNeg = 0, nZero = 0;
  double mostNeg = 0.0;
  int mLam = -1, mK = -1;
  for (int i = 0; i < nw; ++i)
    for (int j = 0; j + 1 < nh; ++j) {
      double const d = tau[i][j + 1] - tau[i][j];
      ++nD;
      if (d < 0.0) {
        ++nNeg;
        if (d < mostNeg) { mostNeg = d; mLam = i; mK = j; }
      }
      if (d == 0.0) ++nZero;
    }
  std::printf("MONOTONICITY: %ld adjacent-height differences ; %ld negative ; %ld exactly zero ; "
              "most negative = %.17g at (lam=%.0f nm, k=%d, h=%.0f..%.0f m)\n",
              nD, nNeg, nZero, mostNeg, (mLam >= 0 ? lam[mLam] : 0.0), mK,
              (mK >= 0 ? hgt[mK] : 0.0), (mK >= 0 ? hgt[mK + 1] : 0.0));
  std::printf("MONOTONICITY ESTABLISHED for THIS profile: %s (statement is per-profile only)\n",
              nNeg == 0 ? "YES" : "NO");

  // production height interpolation: linear in tau, clamped at both ends
  auto tauAt = [&](int i, double h) -> double {
    if (h <= hgt[0]) return tau[i][0];
    if (h >= hgt[nh - 1]) return tau[i][nh - 1];
    int j = 0;
    while (j + 2 < nh && hgt[j + 1] < h) ++j;
    double const dh = hgt[j + 1] - hgt[j];
    double const t = (dh > 0.0) ? (h - hgt[j]) / dh : 0.0;
    return tau[i][j] * (1.0 - t) + tau[i][j + 1] * t;
  };

  long nTot = 0, nZeroLen = 0, nNegWin = 0;
  double mostNegWin = 0.0, wLam = -1.0, wH = 0.0, maxFloor = 0.0;
  auto probe = [&](int i, double h) {
    double const lo = std::max(zmin, h - 25.0);
    double const hi = std::min(zmax, h + 25.0);
    double const d = hi - lo;
    ++nTot;
    if (!(d > 0.0)) { ++nZeroLen; return; }
    double const raw = (tauAt(i, hi) - tauAt(i, lo)) / d;
    if (raw < 0.0) {
      ++nNegWin;
      if (raw < mostNegWin) { mostNegWin = raw; wLam = lam[i]; wH = h; }
    }
    maxFloor = std::max(maxFloor, std::fabs(std::max(0.0, raw) - raw));
  };
  for (int i = 0; i < nw; ++i)
    for (double h = -100.0; h <= 50100.0 + 1e-9; h += 1.0) probe(i, h);
  double const bnd[6] = {zmin - 25.0, zmin, zmin + 25.0, zmax - 25.0, zmax, zmax + 25.0};
  for (int i = 0; i < nw; ++i)
    for (int b = 0; b < 6; ++b)
      for (int s = -1; s <= 1; s += 2) probe(i, bnd[b] + double(s) * 1.0e-6);

  std::printf("FLOOR PROBE: %ld evaluations (105 lambdas x [-100,50100] m at 1 m + 12 boundary "
              "offsets)\n", nTot);
  std::printf("  nonpositive covered-length cases : %ld\n", nZeroLen);
  std::printf("  strictly negative raw windows    : %ld\n", nNegWin);
  std::printf("  most negative raw window         : %.17g m^-1 at (lam=%.0f nm, h=%.6f m)\n",
              mostNegWin, wLam, wH);
  std::printf("  max change due to max(0,.)       : %.17g m^-1\n", maxFloor);
  std::printf("RESULT: floor is %s for this profile\n",
              (maxFloor == 0.0) ? "INACTIVE (mathematically inert)" : "ACTIVE (changes values)");
  return 0;
}
