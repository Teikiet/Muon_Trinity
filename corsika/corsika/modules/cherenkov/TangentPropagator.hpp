#pragma once

/*
 * TangentPropagator.hpp
 *
 * Geometry for the regularised (turning-point) propagation mode of the Cherenkov plane
 * tracer: q(r) = n(r) r, its cell-exact derivatives, the turning radius r_t of a ray of
 * invariant p, and the regularised coefficient
 *
 *     B(v) = ( q(r_t + v^2)^2 - p^2 ) / v^2 ,   r = r_t + v^2 ,   B(0) = 2 p q'(r_t)
 *
 * which appears in
 *
 *     dPhi/dv = 2 p / [ (r_t + v^2) sqrt(B) ],   ds/dv = 2 q / sqrt(B).
 *
 * Why bother: at tangency the ordinary (r, Phi) equations carry a 1/sqrt(r - r_t)
 * singularity, so a step near the perigee is either tiny or inaccurate.  In v the two
 * equations are finite and smooth, which is the whole content of the mode.
 *
 * Numerical care (this is where a naive implementation fails):
 *  * The refractive-index table stores a PIECEWISE-CONSTANT log-slope b_i, so q is only
 *    C1 with a kink in q'' at every 10 m node.  Everything here is evaluated inside the
 *    bracketing cell (n = 1 + N_i exp(b_i (r - r_i))), so q' and q'' are exact for the
 *    model the integrator actually sees, and no finite difference is involved.
 *  * B(v) is a difference of two large numbers (q^2 ~ 4e13) over v^2, so it is evaluated
 *    from the local series 2 p q' + (q'^2 + p q'') v^2 close to the turning point
 *    (v^2 < 1 m) and from the direct form beyond.  The two agree to O(v^4) there.
 *  * The turning-radius solve is a bracketed bisection on the monotone q, and it REFUSES
 *    the ray when q'(r_t) <= 0 (a degenerate turning point or a duct), because the
 *    regularisation assumes a simple turning point.
 *
 * Header-only; depends on RefractiveIndexTable only.
 */

#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace corsika::cherenkov {

  class TangentGeometry {
  public:
    TangentGeometry(RefractiveIndexTable const& ri, double const fPhase)
        : ri_(&ri), fPhase_(fPhase) {}

    bool valid() const { return ri_ != nullptr && ri_->valid(); }
    double dispersionFactor() const { return fPhase_; }

    /// Refractivity N(r) of the model (absolute z = radius).
    double N(double const r) const { return ri_->N(r); }
    /// Phase index n(r) = 1 + N(r) f, the same decomposition the tracers use.
    double n(double const r) const { return 1.0 + fPhase_ * ri_->N(r); }
    /// q(r) = n(r) r.
    double q(double const r) const { return n(r) * r; }

    /// q'(r) = d/dr[(1 + N f) r], exact inside the stored cell.
    double qPrime(double const r) const {
      double const Nv = ri_->N(r);
      double const b = ri_->logSlope(ri_->cellIndex(r)); // dN/dr of this cell
      return 1.0 + fPhase_ * Nv + r * fPhase_ * b * Nv;
    }

    /// q''(r), exact inside the stored cell (N'' = b^2 N there).
    double qSecond(double const r) const {
      double const Nv = ri_->N(r);
      double const b = ri_->logSlope(ri_->cellIndex(r));
      return fPhase_ * b * Nv * (2.0 + r * b);
    }

    /// Invariant p = n r sin(psi) = q(r) sqrt(1 - mu_r^2).
    static double invariant(double const qr, double const muRadial) {
      double const s2 = std::max(0.0, 1.0 - muRadial * muRadial);
      return qr * std::sqrt(s2);
    }

    /**
     * Turning radius: the root of q(r) = p.  Monotone q makes this a bracketed
     * bisection.  Returns false (and leaves rt untouched) when the bracket does not
     * contain a sign change, or when q'(rt) <= 0 there: a degenerate turning point is
     * outside this mode's contract.
     */
    bool turningRadius(double const p, double const rLo, double const rHi,
                       double& rt) const {
      if (!valid() || !(p > 0.0) || !(rHi > rLo)) return false;
      double fLo = q(rLo) - p;
      double fHi = q(rHi) - p;
      if (!(fLo <= 0.0 && fHi >= 0.0)) return false;
      double lo = rLo, hi = rHi;
      for (int it = 0; it < 90; ++it) {
        double const mid = 0.5 * (lo + hi);
        double const fm = q(mid) - p;
        if (fm < 0.0) {
          lo = mid;
        } else {
          hi = mid;
        }
        if (!(hi - lo > 1.0e-9)) break;
      }
      double const rts = 0.5 * (lo + hi);
      if (!(qPrime(rts) > 0.0)) return false; // ducting / degenerate
      rt = rts;
      return true;
    }

    /**
     * B(v) of a ray with invariant p turning at r_t.  Series below V2_SERIES_M, direct
     * form above; both are exact for the stored piecewise-exponential model.
     */
    double B(double const rt, double const p, double const v) const {
      double constexpr kV2Series_m = 1.0;
      double const v2 = v * v;
      if (v2 < kV2Series_m) {
        double const qp = qPrime(rt);
        double const qpp = qSecond(rt);
        return 2.0 * p * qp + (qp * qp + p * qpp) * v2;
      }
      double const qr = q(rt + v2);
      return (qr * qr - p * p) / v2;
    }

    /// sin(psi) = p / q(r); 1 when p exceeds q (should not happen on a valid ray).
    double sinPsi(double const r, double const p) const {
      double const qr = q(r);
      return (qr > 0.0) ? std::min(1.0, p / qr) : 1.0;
    }

  private:
    RefractiveIndexTable const* ri_{nullptr};
    double fPhase_{1.0};
  };

} // namespace corsika::cherenkov
