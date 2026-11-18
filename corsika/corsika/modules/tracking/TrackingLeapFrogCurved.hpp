/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/tracking/TrackingStraight.hpp> // for neutral particles
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/SeparationPlane.hpp>
#include <corsika/framework/geometry/LeapFrogTrajectory.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/Intersections.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/QuarticSolver.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/modules/tracking/Intersect.hpp>

#include <type_traits>
#include <utility>

namespace corsika {

  namespace tracking_leapfrog_curved {

    /**
     * \file TrackingLeapFrogCurved.hpp The leap-frog tracking.
     */

    /**
     * The class tracking_leapfrog_curved::Tracking is based on the
     * Bachelor thesis of Andre Schmidt (KIT). It implements a
     * two-step leap-frog algorithm, but with analytically exact geometric
     * intersections between leap-frog steps and geometric volumes
     * (spheres, planes).
     *
     * Note that leap-frog times and length always reflect the actual properties of
     * the final step. The internal steplength is slightly shorter, because the second
     * halve steps of the algorithm is slightly longer than the first one (in principle
     * violating |v|=const).
     */

    class Tracking : public Intersect<Tracking> {

      using Intersect<Tracking>::nextIntersect;

    public:
      Tracking(double maxDeflection = 0.2)
          : maxMagneticDeflectionAngle_(maxDeflection)
          , straightTracking_{tracking_line::Tracking()} {}

      template <typename TParticle>
      auto getTrack(TParticle const& particle);

      /**
       *  find intersection of Sphere with Track
       *
       * Returns intersection of particle assuming a curved leap-frog step, with a sphere.
       * Entry and exit times are calculated, where the velocity is constant and the
       * steplength is the geometric steplength of the leap-frog.
       *
       * @param particle Current particle state
       * @param sphere Sphere object
       */
      template <typename TParticle>
      static Intersections intersect(TParticle const& particle, Sphere const& sphere);

      /**
       *  find intersection of any Volume node with particle
       *
       * The intersection time(s) of a particle, assuming a curved leap-frog
       * step, are calculated for any volume type.
       */
      template <typename TParticle, typename TBaseNodeType>
      static Intersections intersect(TParticle const& particle,
                                     TBaseNodeType const& node);

      /**
       *  find intersection of Plane with Track
       *
       * Intersection times of particle are caculated with a plane, assuming a curved leap
       * frog trajectory. The intersection time is assuming constant velocity (no change)
       * along the geometric leap-frog step.
       *
       * @tparam TParticle Type of particle object on stack.
       * @param particle Particle initial state.
       * @param plane Plane.
       * @return Intersections in time units.
       */
      template <typename TParticle>
      static Intersections intersect(TParticle const& particle, Plane const& plane);

      /**
       *  find intersection of Separation with Track
       *
       * Intersection times of particle are caculated by
       * and using the intersect-with-plane routine
       *
       * @tparam TParticle Type of particle object on stack.
       * @param particle Particle initial state.
       * @param sepPlane Separation.
       * @return Intersections in time units.
       */
      template <typename TParticle>
      static Intersections intersect(TParticle const& particle,
                                     SeparationPlane const& sepPlane);

      static std::string getName() { return "LeapFrog-curved"; }
      static std::string getVersion() { return "1.0.0"; }

    protected:
      /**
       * Use internally stored class tracking_line::Tracking to
       * perform a straight line tracking, if no magnetic bendig was
       * detected.
       */
      template <typename TParticle>
      auto getLinearTrajectory(TParticle& particle);

    protected:
      double maxMagneticDeflectionAngle_;
      tracking_line::Tracking
          straightTracking_; ///! we want this for neutral and B=0T tracks

    }; // namespace tracking_leapfrog_curved

  } // namespace tracking_leapfrog_curved

} // namespace corsika

#include <corsika/detail/modules/tracking/TrackingLeapFrogCurved.inl>
