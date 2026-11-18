/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Path.hpp>
#include <corsika/framework/geometry/Vector.hpp>

namespace corsika {

  /**
   * Store the photon signal path between two points.
   *
   * This is basically a container class
   */
  struct SignalPath final : public Path {

    // TODO: discuss if we need average refractivity or average refractive index
    TimeType const propagation_time_;       ///< The total propagation time.
    double const average_refractive_index_; ///< The average refractive index.
    double const refractive_index_source_;  ///< The refractive index at the source.
    double const
        refractive_index_destination_; ///< The refractive index at the destination point.
    Vector<dimensionless_d> const emit_;    ///< The (unit-length) emission vector.
    Vector<dimensionless_d> const receive_; ///< The (unit-length) receive vector.
    LengthType const
        R_distance_; ///< The distance from the point of emission to an observer. TODO:
                     ///< optical path, not geometrical! (probably)

    /**
     * Create a new SignalPath instance.
     */
    SignalPath(TimeType const propagation_time, double const average_refractive_index,
               double const refractive_index_source,
               double const refractive_index_destination,
               Vector<dimensionless_d> const& emit,
               Vector<dimensionless_d> const& receive, LengthType const R_distance,
               std::deque<Point> const& points);

  }; // END: class SignalPath final

} // namespace corsika

#include <corsika/detail/modules/radio/propagators/SignalPath.inl>