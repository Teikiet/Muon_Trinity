/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/utility/QuarticSolver.hpp>

namespace corsika {

  inline Line LeapFrogTrajectory::getLine() const {
    auto D = getPosition(1) - getPosition(0);
    auto d = D.getNorm();
    auto v = initialVelocity_;
    if (d > 1_um) { // if trajectory is ultra-short, we do not
                    // re-calculate velocity, just use initial
                    // value. Otherwise, this is numerically unstable
      v = D / d * getVelocity(0).getNorm();
    }
    return Line(getPosition(0), v);
  }

  inline Point LeapFrogTrajectory::getPosition(double const u) const {
    if (u == 0) return initialPosition_;
    Point const position = initialPosition_ + initialVelocity_ * timeStep_ * u / 2;
    VelocityVector const velocity =
        initialVelocity_ + initialVelocity_.cross(magneticfield_) * timeStep_ * u * k_;
    return position + velocity * timeStep_ * u / 2;
  }

  inline VelocityVector LeapFrogTrajectory::getVelocity(double const u) const {
    return getDirection(u) * initialVelocity_.getNorm();
  }

  inline DirectionVector LeapFrogTrajectory::getDirection(double const u) const {
    if (u == 0) return initialDirection_;

    return (initialDirection_ +
            initialDirection_.cross(magneticfield_) * timeStep_ * u * k_)
        .normalized();
  }

  inline TimeType LeapFrogTrajectory::getDuration(double const u) const {
    TimeType const step = timeStep_ * u;
    double const correction = 1;
    // the eventual (delta-L to L) correction factor is:
    //    (initialDirection_ + initialDirection_.cross(magneticfield_) * step *
    //    k_).getNorm();
    return step / 2 * (correction + 1);
  }

  template <typename Particle>
  inline TimeType LeapFrogTrajectory::getTime(Particle const& particle,
                                              double const u) const {
    return particle.getTime() + getDuration(u);
  }

  inline LengthType LeapFrogTrajectory::getLength(double const u) const {
    return getDuration(u) * initialVelocity_.getNorm();
  }

  inline void LeapFrogTrajectory::setLength(LengthType const limit) {
    if (initialVelocity_.getNorm() == SpeedType::zero()) setDuration(0_s);
    setDuration(limit / initialVelocity_.getNorm());
  }

  inline void LeapFrogTrajectory::setDuration(TimeType const limit) {
    /*
    initial attempt to calculate delta-L from assumed full-leap-frog-length L:

    Note: often return 0. Not good enough yet.

    LengthType const L = initialVelocity_.getNorm() * limit; // distance
    double const a = (initialVelocity_.cross(magneticfield_) * k_).getSquaredNorm() / 4 /
                     square(1_m) * static_pow<4>(1_s);
    double const d = L * initialVelocity_.getNorm() / square(1_m) * 1_s;
    double const e = -square(L) / square(1_m);
    std::vector<double> solutions = solve_quartic_real(a, 0, 0, d, e);
    CORSIKA_LOG_DEBUG("setDuration limit={} L={} solution={}", limit, L,
                      fmt::join(solutions, ", "));
    */
    timeStep_ = limit;
  }

} // namespace corsika
