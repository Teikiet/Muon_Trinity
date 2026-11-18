/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/StraightTrajectory.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/modules/tracking/TrackingStraight.hpp>

#include <type_traits>
#include <utility>

namespace corsika {

  namespace tracking_leapfrog_straight {

    /**
     *
     * The class tracking_leapfrog_straight::Tracking inherits from
     * tracking_line::Tracking and adds a (two-step) Leap-Frog
     * algorithms with two halve-steps and magnetic deflection.
     *
     * The two halve steps are implemented as two straight explicit
     * `tracking_line::Tracking`s and all geometry intersections are,
     * thus, based on those two straight line elements.
     *
     * As a precaution for numerical instability, the steplength is
     * limited to correspond to a straight line distance to the next
     * volume intersection. In typical situations this leads to about
     * (at least) one full leap-frog step to the next volume boundary.
     *
     **/

    class Tracking : public tracking_line::Tracking {

    public:
      /**
       * \param firstFraction fraction of first leap-frog halve step
       * relative to full linear step to next volume boundary. This
       * should not be less than 0.5, otherwise you risk that
       * particles will never travel from one volume to the next
       * one. A crossing must be possible (even likely). However, if
       * firstFraction is too big (~1) the resulting calculated
       * numerical error will be largest.
       *
       */
      Tracking(double const firstFraction = 0.55)
          : firstFraction_(firstFraction) {}

      template <typename Particle>
      auto getTrack(Particle& particle);

      static std::string getName() { return "LeapFrogStraight"; }
      static std::string getVersion() { return "1.0.0"; }

    protected:
      double firstFraction_;
    };

  } // namespace tracking_leapfrog_straight

} // namespace corsika

#include <corsika/detail/modules/tracking/TrackingLeapFrogStraight.inl>
