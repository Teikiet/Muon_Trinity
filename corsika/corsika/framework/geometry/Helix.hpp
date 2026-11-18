/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <cmath>

namespace corsika {

  /*!
   * Defines a helical path
   *
   * A Helix is defined by the cyclotron frequency \f$ \omega_c \f$, the initial
   * Point r0 and
   * the velocity vectors \f$ \vec{v}_{\parallel} \f$ and \f$ \vec{v}_{\perp} \f$
   * denoting the projections of the initial velocity \f$ \vec{v}_0 \f$ parallel
   * and perpendicular to the axis \f$ \vec{B} \f$, respectively, i.e.
   * \f{align*}{
        \vec{v}_{\parallel} &= \frac{\vec{v}_0 \cdot \vec{B}}{\vec{B}^2} \vec{B} \\
        \vec{v}_{\perp} &= \vec{v}_0 - \vec{v}_{\parallel}
     \f}
   */

  class Helix {

  public:
    Helix(Point const& pR0, FrequencyType pOmegaC, VelocityVector const& pvPar,
          VelocityVector const& pvPerp)
        : r0_(pR0)
        , omegaC_(pOmegaC)
        , vPar_(pvPar)
        , vPerp_(pvPerp)
        , uPerp_(vPerp_.cross(vPar_.normalized()))
        , radius_(pvPar.getNorm() / abs(pOmegaC)) {}

    LengthType getRadius() const;

    Point getPosition(TimeType const t) const;

    VelocityVector getVelocity(TimeType const t) const;

    Point getPositionFromArclength(LengthType const l) const;

    LengthType getArcLength(TimeType const t1, TimeType const t2) const;

    TimeType getTimeFromArclength(LengthType const l) const;

  private:
    Point r0_;             ///! origin of helix, but this is in the center of the
                           ///! "cylinder" on which the helix rotates
    FrequencyType omegaC_; ///! speed of angular rotation
    VelocityVector vPar_;  ///! speed along direction of "cylinder"
    VelocityVector vPerp_, uPerp_;
    LengthType radius_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/Helix.inl>
