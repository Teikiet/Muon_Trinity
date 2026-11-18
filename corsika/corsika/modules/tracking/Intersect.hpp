/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/geometry/Intersections.hpp>

#include <limits>

namespace corsika {

  /**
   *
   * This is a CRTP class to provide a generic volume-tree
   * intersection for the purpose of tracking.
   *
   * It return the closest distance in time to the next geometric
   * intersection, as well as a pointer to the corresponding new
   * volume.
   *
   * User may provide an optional global step-length limit as
   * parameter. This may be needd for (simpler) algorithms in magnetic
   * fields, where tracking errors grow linearly with step-length.
   * Obviously, in the case of the step-length limit, the returend
   * "next" volume is just the current one.
   *
   **/

  template <typename TDerived>
  class Intersect {

  protected:
    /**
     * Determines next intersection with any of the geometry volumes.
     *
     *
     */
    template <typename TParticle>
    auto nextIntersect(TParticle const& particle,
                       TimeType const step_limit =
                           std::numeric_limits<TimeType::value_type>::infinity() *
                           second) const;
  };

} // namespace corsika

#include <corsika/detail/modules/tracking/Intersect.inl>
