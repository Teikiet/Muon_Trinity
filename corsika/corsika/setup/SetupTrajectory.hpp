/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/StraightTrajectory.hpp>
#include <corsika/framework/geometry/LeapFrogTrajectory.hpp>
#include <corsika/modules/Tracking.hpp>

namespace corsika::setup {

  /**
    \file SetupTrajectory.hpp

    Note/Warning:     Tracking and Trajectory must fit together !

    tracking_leapfrog_curved::Tracking is the result of the Bachelor
    thesis of Andre Schmidt, KIT. This is a leap-frog algorithm with
    an analytical, precise calculation of volume intersections. This
    algorithm needs a LeapFrogTrajectory.

    tracking_leapfrog_straight::Tracking is a more simple and direct
    leap-frog implementation. The two halve steps are coded explicitly
    as two straight segments. Intersections with other volumes are
    calculate only on the straight segments. This algorithm is based
    on LineTrajectory.

    tracking_line::Tracking is a pure straight tracker. It is based on
    LineTrajectory.
   */

  /**
     The default tracking algorithm.
   */

  typedef corsika::tracking_leapfrog_curved::Tracking Tracking;
  // typedef corsika::tracking_leapfrog_straight::Tracking Tracking;
  // typedef corsika::tracking_line::Tracking Tracking;

  /**
   The default trajectory.
  */
  /// definition of Trajectory base class, to be used in tracking and cascades
  // typedef StraightTrajectory Trajectory;
  typedef corsika::LeapFrogTrajectory Trajectory;

} // namespace corsika::setup
