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

namespace corsika {

  inline VelocityVector StraightTrajectory::getVelocity(double const u) const {
    return initialVelocity_ * (1 - u) + finalVelocity_ * u;
  }

  inline TimeType StraightTrajectory::getDuration(double const u) const {
    return u * timeStep_;
  }

  template <typename Particle>
  inline TimeType StraightTrajectory::getTime(Particle const& particle,
                                              double const u) const {
    return particle.getTime() + getDuration(u); // timeStep_ * u;
  }

  inline LengthType StraightTrajectory::getLength(double const u) const {
    if (timeLength_ == 0_s) return 0_m;
    if (timeStep_ == std::numeric_limits<TimeType::value_type>::infinity() * 1_s)
      return std::numeric_limits<LengthType::value_type>::infinity() * 1_m;
    return getDistance(u) * timeStep_ / timeLength_;
  }

  inline void StraightTrajectory::setLength(LengthType const limit) {
    setDuration(line_.getTimeFromArclength(limit));
  }

  inline void StraightTrajectory::setDuration(TimeType const limit) {
    if (timeStep_ == 0_s) {
      timeLength_ = 0_s;
      setFinalVelocity(getVelocity(0));
      timeStep_ = limit;
    } else {
      // for infinite steps there can't be a difference between
      // curved and straight trajectory, this is fundamentally
      // undefined: assume they are the same (which, i.e. is always correct for a
      // straight line trajectory).
      //
      // Final note: only straight-line trajectories should have
      // infinite steps! Everything else is ill-defined.
      if (timeStep_ == std::numeric_limits<TimeType::value_type>::infinity() * 1_s ||
          timeLength_ == std::numeric_limits<TimeType::value_type>::infinity() * 1_s) {
        timeLength_ = limit;
        timeStep_ = limit;
        // ...and don't touch velocity
      } else {
        const double scale = limit / timeStep_;
        timeLength_ *= scale;
        setFinalVelocity(getVelocity(scale));
        timeStep_ = limit;
      }
    }
  }

  inline LengthType StraightTrajectory::getDistance(double const u) const {
    if (!(0 <= u && u <= 1)) {
      CORSIKA_LOG_ERROR("argument must be 0 <= {} <= 1", u);
      throw std::runtime_error("value out of range");
    }

    return line_.getArcLength(0 * second, u * timeLength_);
  }

} // namespace corsika
