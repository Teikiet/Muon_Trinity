// ============================================================================
// stage2e_invariance_b.cpp -- Checkpoint 2e-c instrumentation invariance gate.
//
//   g++ ... -DREF_TRACER_HEADER='"ref_tracer_b.hpp.pre2eaudit"' \
//       stage2e_invariance_b.cpp -o inv_base
//   g++ ... -DREF_TRACER_HEADER='"ref_tracer_b.hpp"' \
//       stage2e_invariance_b.cpp -o inv_new
//   ./inv_base > a.txt ; ./inv_new > b.txt ; diff a.txt b.txt
//
// Prints, with 17 significant digits, only the tracer outputs that the audit
// must not perturb: success status, impact coordinates, optical path, geometric
// path and accumulated tau. Identical text for both headers (same compiler,
// same flags, same inputs) is the invariance gate.
// The same three reference rays and the same callables as stage2e_audit_b.cpp
// are used; no physics is recomputed here.
// ============================================================================
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

// Header selected on the command line: baseline snapshot or the instrumented copy.
#include REF_TRACER_HEADER

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using corsika::cherenkov::AtmosphericAbsorption;
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
static double alphaProd(AtmosphericAbsorption const& abs, double lam, double h, double hLo,
                        double hHi) {
  double const hl = std::max(hLo, h - 25.0), hh = std::min(hHi, h + 25.0);
  double const adh = hh - hl;
  if (adh <= 0.0) return 0.0;
  double const t1 = abs.getOpticalDepth(lam, hh);
  double const t0 = abs.getOpticalDepth(lam, hl);
  return std::max(0.0, (t1 - t0) / adh);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::printf("usage: %s <profile .dat>\n", argv[0]);
    return 2;
  }
  AtmosphericAbsorption abs;
  try {
    abs.load(argv[1]);
  } catch (std::exception const& e) {
    std::printf("STOP: profile load failed: %s\n", e.what());
    return 1;
  }
  auto const hr = abs.getHeightRange();
  RefractiveIndexTable ri = makeRI();

  struct C {
    double lam, mu, h;
  };
  C const cs[3] = {{300.0, 0.6, 46000.0}, {480.0, 0.1, 1000.0}, {480.0, 0.1, 975.0}};

  std::printf("=== 2e-c instrumentation invariance input dump ===\n");
  for (int q = 0; q < 3; ++q) {
    double const lm = cs[q].lam, mu = cs[q].mu, hh = cs[q].h;
    refdiag::RefTracer tr;
    tr.earthRadiusM_ = R;
    tr.enableAbsorption_ = true;
    tr.refractiveIndexFunc_ = [&](double h_m, double lam_nm) {
      double const f = Fdisp(lam_nm);
      double const n = 1.0 + ri.N(R + h_m) * f;
      return std::make_pair(n, n);
    };
    tr.extinctionCoeffFunc_ = [&](double h_m, double lam_nm) {
      return alphaProd(abs, lam_nm, h_m, hr.first, hr.second);
    };
    double const sinT = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    refdiag::RefTracer::Result const r =
        tr.trace(0.0, 0.0, R + hh, sinT, 0.0, -mu, 0.0, 0.0, R, 0.0, 0.0, 1.0, lm, 0.0, 0.0);
    std::printf("ray %d hit=%d\n", q, int(r.hit));
    std::printf("  x=%.17g y=%.17g z=%.17g\n", r.x_impact, r.y_impact, r.z_impact);
    std::printf("  optical=%.17g geometric=%.17g tau=%.17g steps=%d\n", r.optical_path,
                r.geometric_path, r.accumulated_tau, r.steps);
  }
  std::printf("=== end ===\n");
  return 0;
}
