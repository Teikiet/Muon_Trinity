/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {

  inline SignalPath::SignalPath(
      TimeType const propagation_time, double const average_refractive_index,
      double const refractive_index_source, double const refractive_index_destination,
      Vector<dimensionless_d> const& emit, Vector<dimensionless_d> const& receive,
      LengthType const R_distance, std::deque<Point> const& points)
      : Path(points)
      , propagation_time_(propagation_time)
      , average_refractive_index_(average_refractive_index)
      , refractive_index_source_(refractive_index_source)
      , refractive_index_destination_(refractive_index_destination)
      , emit_(emit)
      , receive_(receive)
      , R_distance_(R_distance) {}

} // namespace corsika