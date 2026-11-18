/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Vector.hpp>

namespace corsika {

  inline Point Line::getPosition(TimeType const t) const {
    return start_point_ + velocity_ * t;
  }

  inline VelocityVector const& Line::getVelocity(TimeType const) const {
    return velocity_;
  }

  inline Point Line::getPositionFromArclength(LengthType const l) const {
    return start_point_ + velocity_.normalized() * l;
  }

  inline LengthType Line::getArcLength(TimeType const t1, TimeType const t2) const {
    return velocity_.getNorm() * (t2 - t1);
  }

  inline TimeType Line::getTimeFromArclength(LengthType const t) const {
    return t / velocity_.getNorm();
  }

  inline Point const& Line::getStartPoint() const { return start_point_; }

  inline DirectionVector Line::getDirection() const { return velocity_.normalized(); }

  inline VelocityVector const& Line::getVelocity() const { return velocity_; }

} // namespace corsika
