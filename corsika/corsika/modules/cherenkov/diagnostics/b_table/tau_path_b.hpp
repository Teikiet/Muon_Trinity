#pragma once
// ============================================================================
// tau_path_b.hpp -- DIAGNOSTIC-ONLY helper for Checkpoint 2e-e.
//
// NOT production code, NOT in any build. It does not modify tau_march_b.hpp,
// stage9b_abs_b.cpp or any other existing file; it INCLUDES tau_march_b.hpp to
// reuse the identical step function (refmarch::rk4Step) and the identical
// clip-final endpoint treatment of checkpoint 2e-d.
//
// It records, for every ACCEPTED step of a march, the pre-step state, the
// accepted length (shortened for the clipped final step) and the post-step
// state. Interior points are defined by
//     y_i(u) = RK4Step(y_i, u),   0 <= u <= L_i,
// always starting from the ORIGINAL pre-step state -- no chaining of substeps,
// no replacement of post-step states, no change of the trajectory.
//
// This is a numerical interpolating path for the experiment; it is NOT an exact
// physical ray and NOT an arc-length-validated interpolant.
// ============================================================================
#include "tau_march_b.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace refpath {

struct Step {
  long index = 0;      // 1-based accepted-step index
  double zPre = 0.0;   // radius from the Earth centre [m]
  double thPre = 0.0;  // angle to the outward radial [rad]
  double L = 0.0;      // accepted length [m] (shortened for the final step)
  double zPost = 0.0;
  double thPost = 0.0;
};

struct Path {
  std::vector<Step> steps;
  double tau = 0.0;
  bool ok = false;
  bool clipFinal = false;
};

// Same loop, same arithmetic as refmarch::march(); adds one record per accepted
// step. clipFinal=true is the exact 2e-d endpoint treatment.
inline Path buildPath(double const R, double const hTop, double const ds, double const p,
                      refmarch::NFn const& Nof, refmarch::DerFn const& der,
                      refmarch::AlphaFn const& alpha, bool const clipFinal) {
  Path out;
  out.clipFinal = clipFinal;

  double z = R + hTop;
  double const nTop = 1.0 + Nof(z);
  double const sinT = std::min(1.0, p / (nTop * z));
  double th = std::acos(-std::sqrt(std::max(0.0, 1.0 - sinT * sinT)));
  double tau = 0.0;
  long n = 0;

  while (z > R + 1e-9 && n < 200000000L) {
    double zNext, thNext;
    refmarch::rk4Step(z, th, ds, der, zNext, thNext);
    double const zPrev = z;
    double const thPrev = th;
    bool const isFinal = !(zNext > R + 1e-9);

    if (!isFinal) {
      z = zNext;
      th = thNext;
      double const h = std::max(0.0, 0.5 * (zPrev + z) - R);
      tau += alpha(h) * ds;
      ++n;
      out.steps.push_back(Step{n, zPrev, thPrev, ds, zNext, thNext});
      continue;
    }

    // ---- final step: identical treatment to refmarch::march() --------------
    if (!clipFinal) {
      double const hFull = std::max(0.0, 0.5 * (zPrev + zNext) - R);
      tau += alpha(hFull) * ds;
      ++n;
      out.steps.push_back(Step{n, zPrev, thPrev, ds, zNext, thNext});
      out.tau = tau;
      out.ok = true;
      return out;
    }

    bool const bracketed = (zPrev > R) && (zNext <= R);
    double rPrev = zPrev;
    bool mono = true;
    int const K = 64;
    for (int k = 1; k <= K; ++k) {
      double zz, tt;
      refmarch::rk4Step(zPrev, thPrev, ds * (static_cast<double>(k) / static_cast<double>(K)),
                        der, zz, tt);
      if (zz - rPrev > 0.0) mono = false;
      rPrev = zz;
    }
    if (!bracketed || !mono) {  // caller reports; nothing is invented
      out.tau = tau;
      out.ok = false;
      return out;
    }

    double lo = 0.0, hi = ds;
    long it = 0;
    while ((hi - lo) > 1e-7 && it < 80) {
      double const mid = 0.5 * (lo + hi);
      double zz, tt;
      refmarch::rk4Step(zPrev, thPrev, mid, der, zz, tt);
      if (zz > R) {
        lo = mid;
      } else {
        hi = mid;
      }
      ++it;
    }
    double zHi, thHi;
    refmarch::rk4Step(zPrev, thPrev, hi, der, zHi, thHi);
    double const aClip = alpha(std::max(0.0, 0.5 * (zPrev + zHi) - R));
    tau += aClip * hi;
    ++n;
    out.steps.push_back(Step{n, zPrev, thPrev, hi, zHi, thHi});
    z = zHi;
    th = thHi;
    out.tau = tau;
    out.ok = true;
    return out;
  }

  out.tau = tau;
  return out;
}

// ============================================================================
// Composite-midpoint re-integration of extinction ON THE FIXED PATH.
//
//   per accepted step:  N = ceil(L/q),  w = L/N,
//                       u_j = (j + 0.5) * w,
//                       h_j = max(0, radius(y_i(u_j)) - R),
//                       contribution_i = w * sum_j alphaProd(h_j, lambda)
//
// Interior evaluations start from the original pre-step state. Totals use
// Neumaier compensated summation.
// ============================================================================
struct Quad {
  double sum = 0.0;   // running compensated sum
  double comp = 0.0;  // compensation term
  long panels = 0;
  double totalLength = 0.0;
  std::vector<double> perStep;  // compensated contribution per accepted step
  std::vector<long> panelsPerStep;
  double value() const { return sum + comp; }
};

inline void neumaierAdd(double& sum, double& c, double const x) {
  double const t = sum + x;
  if (std::fabs(sum) >= std::fabs(x)) {
    c += (sum - t) + x;
  } else {
    c += (x - t) + sum;
  }
  sum = t;
}

inline Quad compositeMidpoint(Path const& path, double const q, refmarch::DerFn const& der,
                              refmarch::AlphaFn const& alpha, double const R) {
  Quad res;
  res.perStep.reserve(path.steps.size());
  res.panelsPerStep.reserve(path.steps.size());
  for (std::size_t i = 0; i < path.steps.size(); ++i) {
    Step const& s = path.steps[i];
    long const N = std::max(1L, static_cast<long>(std::ceil(s.L / q)));
    double const w = s.L / static_cast<double>(N);
    double ssum = 0.0, scomp = 0.0;
    for (long j = 0; j < N; ++j) {
      double const u = (static_cast<double>(j) + 0.5) * w;
      double zz, tt;
      refmarch::rk4Step(s.zPre, s.thPre, u, der, zz, tt);
      double const h = std::max(0.0, zz - R);
      neumaierAdd(ssum, scomp, w * alpha(h));
    }
    double const contribution = ssum + scomp;
    res.perStep.push_back(contribution);
    res.panelsPerStep.push_back(N);
    neumaierAdd(res.sum, res.comp, contribution);
    res.panels += N;
    res.totalLength += s.L;
  }
  return res;
}

// Endpoint reproduction: y_i(0) vs the recorded pre-state and y_i(L_i) vs the
// recorded post-step state.
struct EndpointResiduals {
  double maxRadius = 0.0;  // [m]
  double maxAngle = 0.0;   // [rad]
};

inline EndpointResiduals endpointResiduals(Path const& path, refmarch::DerFn const& der) {
  EndpointResiduals out;
  for (std::size_t i = 0; i < path.steps.size(); ++i) {
    Step const& s = path.steps[i];
    double z0, t0, zL, tL;
    refmarch::rk4Step(s.zPre, s.thPre, 0.0, der, z0, t0);
    refmarch::rk4Step(s.zPre, s.thPre, s.L, der, zL, tL);
    out.maxRadius = std::max(out.maxRadius, std::fabs(z0 - s.zPre));
    out.maxRadius = std::max(out.maxRadius, std::fabs(zL - s.zPost));
    out.maxAngle = std::max(out.maxAngle, std::fabs(t0 - s.thPre));
    out.maxAngle = std::max(out.maxAngle, std::fabs(tL - s.thPost));
  }
  return out;
}

}  // namespace refpath

