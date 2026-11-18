/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Helix.hpp>
#include <corsika/framework/geometry/StraightTrajectory.hpp>

#include <corsika/framework/geometry/LeapFrogTrajectory.hpp>
// #include <corsika/modules/TrackingLine.hpp>
// #include <corsika/modules/TrackingCurved.hpp> // simple leap-frog implementation
// #include <corsika/modules/TrackingLeapFrog.hpp> // more complete leap-frog
//  implementation

namespace corsika::setup::testing {

  template <typename TTrack>
  TTrack make_track(
      Line const line,
      TimeType const tEnd = std::numeric_limits<TimeType::value_type>::infinity() * 1_s);

  template <>
  inline StraightTrajectory make_track<StraightTrajectory>(Line const line,
                                                           TimeType const tEnd) {
    return StraightTrajectory(line, tEnd);
  }

  template <>
  inline LeapFrogTrajectory make_track<LeapFrogTrajectory>(Line const line,
                                                           TimeType const tEnd) {

    auto const k = square(0_m) / (square(1_s) * 1_V);
    return LeapFrogTrajectory(line.getStartPoint(), line.getVelocity(),
                              MagneticFieldVector{0_T, 0_T, 0_T}, k, tEnd);
  }

} // namespace corsika::setup::testing
