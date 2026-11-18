/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>

namespace corsika {

  /**
   *
   * A Trajectory is a description of a momvement of an object in
   * three-dimensional space that describes the trajectory (connection
   * between two Points in space), as well as the direction of motion
   * at any given point.
   *
   * A Trajectory has a start `0` and an end `1`, where
   * e.g. getPosition(0) returns the start point and getDirection(1)
   * the direction of motion at the end. Values outside 0...1 are not
   * defined.
   *
   * A Trajectory has a length in [m], getLength, a duration in [s], getDuration.
   *
   * Note: so far it is assumed that the speed (d|vec{r}|/dt) between
   * start and end does not change and is constant for the entire
   * Trajectory.
   *
   **/

  class BaseTrajectory {

  public:
    virtual Point getPosition(double const u) const = 0;

    virtual VelocityVector getVelocity(double const u) const = 0;

    virtual DirectionVector getDirection(double const u) const = 0;

    virtual TimeType getDuration(double const u = 1) const = 0;

    virtual LengthType getLength(double const u = 1) const = 0;

    virtual void setLength(LengthType const limit) = 0;

    virtual void setDuration(TimeType const limit) = 0;
  };

} // namespace corsika
