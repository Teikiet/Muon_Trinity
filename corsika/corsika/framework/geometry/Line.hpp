/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>

namespace corsika {

  /**
   *
   * A Line describes a movement in three dimensional space. It
   * consists of a Point `$\vec{p_0}$` and and a speed-Vector
   * `$\vec{v}$`, so that it can return GetPosition as
   * `$\vec{p_0}*\vec{v}*t$` for any value of time `$t$`.
   *
   **/

  class Line {

  public:
    Line(Point const& pR0, VelocityVector const& pV0)
        : start_point_(pR0)
        , velocity_(pV0) {}

    Point getPosition(TimeType const t) const;

    VelocityVector const& getVelocity(TimeType const) const;

    Point getPositionFromArclength(LengthType const l) const;

    LengthType getArcLength(TimeType const t1, TimeType const t2) const;

    TimeType getTimeFromArclength(LengthType const t) const;

    Point const& getStartPoint() const;

    DirectionVector getDirection() const;

    VelocityVector const& getVelocity() const;

  private:
    Point start_point_;
    VelocityVector velocity_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/Line.inl>
