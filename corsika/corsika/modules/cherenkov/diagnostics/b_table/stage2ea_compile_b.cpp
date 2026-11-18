// ============================================================================
// stage2ea_compile_b.cpp -- COMPILE-ONLY harness for ref_tracer_b.hpp.
//
// Checkpoint 2e-a. This program is built but NEVER executed. Its only purpose
// is to force the compiler to check the extracted tracer's body and its
// callable interface. No physics is exercised here.
// ============================================================================

#include "ref_tracer_b.hpp"

#include <utility>

namespace {

refdiag::RefTracer make_tracer() {
  refdiag::RefTracer t;
  t.earthRadiusM_ = 6371000.0;
  // Toy, monotone, physically-shaped profiles only: this harness never runs.
  t.refractiveIndexFunc_ = [](double z_m, double lam_nm) {
    (void)lam_nm;
    double const N = 2.9e-4 * std::exp(-z_m / 8000.0);
    return std::make_pair(1.0 + N, 1.0 + 1.035 * N);
  };
  t.extinctionCoeffFunc_ = [](double h_m, double lam_nm) {
    return 1.0e-4 * std::exp(-h_m / 1200.0) + 1.0e-7 * (lam_nm / 500.0);
  };
  return t;
}

}  // namespace

// External linkage (not static) so -Wall -Wextra does not flag it as unused
// while still guaranteeing the body is compiled at -O0.
double stage2ea_probe() {
  refdiag::RefTracer t = make_tracer();

  double const R = 6371000.0;
  refdiag::RefTracer::Result const r = t.trace(
      /*x0,y0,z0_altitude*/ 0.0, 0.0, R + 46000.0,
      /*dx,dy,dz        */ 0.0, 0.0, -1.0,
      /*plane_x,y,z     */ 0.0, 0.0, R,
      /*norm_x,y,z      */ 0.0, 0.0, 1.0,
      /*lambdaNm        */ 400.0,
      /*survivalThreshold*/ 0.0,
      /*tau_start       */ 0.0);

  double out = r.optical_path + r.geometric_path + r.accumulated_tau;
  out += r.x_impact + r.y_impact + r.z_impact;
  out += static_cast<double>(r.steps);
  if (r.hit) out += 1.0;
  return out;
}

int main() {
  // Deliberately does not call stage2ea_probe(): compile-only harness.
  double (*p)() = &stage2ea_probe;
  (void)p;
  return 0;
}
