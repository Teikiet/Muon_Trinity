#pragma once
// ============================================================================
// tau_march_b.hpp -- DIAGNOSTIC-ONLY helper for Checkpoint 2e-d.
//
// NOT production code, NOT in any build, does not touch stage9b_abs_b.cpp.
//
// It re-expresses exactly ONE algorithm: the path-length RK4 march of
// stage9b_abs_b.cpp:166-202 (tauIntS). Differences are confined to
//   (a) recording the final (ground-crossing) step, and
//   (b) an ENDPOINT-TREATMENT EXPERIMENT that re-runs only that final step at a
//       shorter length, found by bisection with the SAME step function.
//
// The propagation step is the same four-stage combination with the same
// derivative expression; the step length is simply a parameter. No new
// quadrature rule is introduced -- tau still receives ONE midpoint sample of the
// same callable, multiplied by the length of the step actually taken.
//
// THIS IS AN ENDPOINT-TREATMENT EXPERIMENT, NOT A PHYSICALLY VALIDATED
// REFERENCE.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <functional>

namespace refmarch {

// exactly tauIntS:176-182, supplied by the caller (captures ri and f)
using DerFn = std::function<void(double, double, double&, double&)>;
// h -> extinction coefficient [1/m]; the caller passes alphaProd or a constant
using AlphaFn = std::function<double(double)>;
// r -> N(r)*f, i.e. the refractivity used by the step (ri.N(r)*f in tauIntS)
using NFn = std::function<double(double)>;

// Same four-stage combination as tauIntS:184-191, with the step length as a
// parameter. zNext/thNext are written exactly where tauIntS would write z/th.
inline void rk4Step(double const z, double const th, double const s, DerFn const& der,
                    double& zNext, double& thNext) {
  double dr1, dt1, dr2, dt2, dr3, dt3, dr4, dt4;
  der(z, th, dr1, dt1);
  der(z + 0.5 * s * dr1, th + 0.5 * s * dt1, dr2, dt2);
  der(z + 0.5 * s * dr2, th + 0.5 * s * dt2, dr3, dt3);
  der(z + s * dr3, th + s * dt3, dr4, dt4);
  zNext = z + s / 6.0 * (dr1 + 2 * dr2 + 2 * dr3 + dr4);
  thNext = th + s / 6.0 * (dt1 + 2 * dt2 + 2 * dt3 + dt4);
}

struct FinalStepRecord {
  bool valid = false;          // a final ground step was reached
  bool terminatedByStepLimit = false;  // no final step: n hit the step limit
  long stepIndex = 0;          // 1-based index of the final step
  // state immediately before the final step
  double zPre = 0.0, thPre = 0.0, hPre = 0.0, tauBefore = 0.0;
  // state after the FULL final step (as tauIntS performs it)
  double zPostFull = 0.0, thPostFull = 0.0, hPostFull = 0.0;
  double dsNominal = 0.0;
  double hSampleFull = 0.0;    // the single extinction sample altitude
  double alphaSampleFull = 0.0;// the coefficient used
  double dTauFull = 0.0;       // final-step tau contribution, full length
  double tauTotalFull = 0.0;   // tau returned by the march
  long nSteps = 0;             // step count (nominal path = nSteps*ds)
  // shortened (clipped) final step
  bool bracketed = false;      // zPre > R and z(ds) <= R
  bool monotonicDown = false;  // r(s) non-increasing on a 64-point grid
  double worstBackward = 0.0;  // worst r(s_{i+1}) - r(s_i) of that grid [m]
  long bisectIters = 0;
  double bracketWidth = 0.0;   // final [lo,hi] width [m]
  double dsHit = 0.0;          // chosen shortened length (below-ground side)
  double zPostClip = 0.0, thPostClip = 0.0, hPostClip = 0.0;
  double hSampleClip = 0.0, alphaSampleClip = 0.0;
  double dTauClip = 0.0;
  double tauTotalClip = 0.0;
};

struct MarchResult {
  double tau = 0.0;
  long n = 0;
  FinalStepRecord fin;
};

// ============================================================================
// The march. State and update order are those of stage9b_abs_b.cpp:170-199.
//
//   z     radius from the Earth centre [m]      (tauIntS calls it "z")
//   th    angle between the ray and the outward radial [rad], downgoing
//   tau   accumulated optical depth [dimensionless]
//   n     completed steps
//
// Per step: the FULL length ds is applied through rk4Step, then
//   tau += alpha(max(0, 0.5*(zPrev+z) - R)) * ds
// exactly as tauIntS:189-193 (one midpoint sample, clamped at ground level).
//
// clipFinal == false reproduces tauIntS bit for the same inputs.
// ============================================================================
inline MarchResult march(double const R, double const hTop, double const ds,
                         double const p, NFn const& Nof, DerFn const& der,
                         AlphaFn const& alpha, bool const clipFinal) {
  MarchResult out;

  // --- tauIntS:170-175, verbatim (ri.N(z)*f -> Nof(z)) ----------------------
  double z = R + hTop;
  double const nTop = 1.0 + Nof(z);
  double const sinT = std::min(1.0, p / (nTop * z));
  double th = std::acos(-std::sqrt(std::max(0.0, 1.0 - sinT * sinT)));  // downgoing
  double tau = 0.0;
  long n = 0;

  while (z > R + 1e-9 && n < 200000000L) {
    double zNext, thNext;
    rk4Step(z, th, ds, der, zNext, thNext);
    double const zPrev = z;
    bool const isFinal = !(zNext > R + 1e-9);

    if (!isFinal) {
      // identical arithmetic to tauIntS:189-193
      z = zNext;
      th = thNext;
      double const h = std::max(0.0, 0.5 * (zPrev + z) - R);
      tau += alpha(h) * ds;
      ++n;
      continue;
    }

    // ---------------- the final (ground-crossing) step ---------------------
    FinalStepRecord& r = out.fin;
    r.valid = true;
    r.stepIndex = n + 1;
    r.zPre = zPrev;
    r.thPre = th;
    r.hPre = zPrev - R;
    r.tauBefore = tau;
    r.dsNominal = ds;
    r.zPostFull = zNext;
    r.thPostFull = thNext;
    r.hPostFull = zNext - R;
    r.hSampleFull = std::max(0.0, 0.5 * (zPrev + zNext) - R);
    r.alphaSampleFull = alpha(r.hSampleFull);
    r.dTauFull = r.alphaSampleFull * ds;
    // tau after the FULL final step is a report value in both modes; only the
    // FULL run actually returns it.
    r.tauTotalFull = r.tauBefore + r.dTauFull;
    ++n;
    out.fin.nSteps = n;

    if (!clipFinal) {
      tau += r.dTauFull;
      z = zNext;
      th = thNext;
      break;
    }

    // ---- ENDPOINT-TREATMENT EXPERIMENT: rerun the final step, shortened ----
    // (a) bracketing + monotonicity on the same step function
    r.bracketed = (zPrev > R) && (zNext <= R);
    double rPrev = zPrev;
    double worstBack = 0.0;
    bool mono = true;
    int const K = 64;
    for (int k = 1; k <= K; ++k) {
      double zz, tt;
      rk4Step(zPrev, th, ds * (static_cast<double>(k) / static_cast<double>(K)), der, zz, tt);
      double const d = zz - rPrev;
      if (d > 0.0) mono = false;
      worstBack = std::max(worstBack, d);
      rPrev = zz;
    }
    r.monotonicDown = mono;
    r.worstBackward = worstBack;
    if (!r.bracketed || !mono) {
      z = zNext;
      th = thNext;
      break;  // caller reports the failure and stops
    }

    // (b) bisection for the first ground intersection in [0, ds]
    double lo = 0.0, hi = ds;
    long it = 0;
    while ((hi - lo) > 1e-7 && it < 80) {
      double const mid = 0.5 * (lo + hi);
      double zz, tt;
      rk4Step(zPrev, th, mid, der, zz, tt);
      if (zz > R) {
        lo = mid;
      } else {
        hi = mid;
      }
      ++it;
    }
    r.bisectIters = it;
    r.bracketWidth = hi - lo;
    r.dsHit = hi;
    double zHi, thHi;
    rk4Step(zPrev, th, hi, der, zHi, thHi);
    r.zPostClip = zHi;
    r.thPostClip = thHi;
    r.hPostClip = zHi - R;
    // the same quadrature rule, sample position recomputed for the shorter step
    r.hSampleClip = std::max(0.0, 0.5 * (zPrev + zHi) - R);
    r.alphaSampleClip = alpha(r.hSampleClip);
    r.dTauClip = r.alphaSampleClip * hi;
    tau += r.dTauClip;
    out.fin.tauTotalClip = tau;

    z = zHi;
    th = thHi;
    break;
  }

  if (!out.fin.valid) {
    out.fin.terminatedByStepLimit = true;
    out.fin.nSteps = n;
    out.fin.tauTotalFull = tau;
  }
  out.tau = tau;
  out.n = n;
  return out;
}

}  // namespace refmarch
