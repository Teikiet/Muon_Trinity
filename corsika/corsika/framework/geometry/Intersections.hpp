/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>

#include <map> // for pair

namespace corsika {

  /**
   *
   * Container to store and return a list of intersections of a
   * trajectory with a geometric volume objects in space.
   *
   **/

  class Intersections {

    Intersections(const Intersections&) = delete;
    Intersections(Intersections&&) = delete;
    Intersections& operator=(const Intersections&) = delete;

  public:
    Intersections()
        : has_intersections_(false) {}

    Intersections(TimeType&& t1, TimeType&& t2)
        : has_intersections_(true)
        , intersections_(std::make_pair(t1, t2)) {}

    Intersections(TimeType&& t)
        : has_intersections_(true)
        , intersections_(std::make_pair(
              t, std::numeric_limits<TimeType::value_type>::infinity() * second)) {}

    bool hasIntersections() const { return has_intersections_; }

    ///! where did the trajectory currently enter the volume
    TimeType getEntry() const {
      if (has_intersections_)
        return intersections_.first;
      else
        return std::numeric_limits<TimeType::value_type>::infinity() * second;
    }

    ///! where did the trajectory currently exit the volume
    TimeType getExit() const {
      if (has_intersections_)
        return intersections_.second;
      else
        return std::numeric_limits<TimeType::value_type>::infinity() * second;
    }

  private:
    bool has_intersections_;
    std::pair<TimeType, TimeType> intersections_;
  };

} // namespace corsika
