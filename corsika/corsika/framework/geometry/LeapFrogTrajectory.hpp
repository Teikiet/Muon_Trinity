/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {

  /**
   * The LeapFrogTrajectory stores information on one leap-frog step.
   *
   * The leap-frog algorithm uses two half-steps and is used in magnetic
   * field tracking. The LeapFrogTrajectory will solve the leap-frog
   * algorithm equation for a given constant @f$k@f$ that has to be
   * specified during construction (essentially fixing the magnetic
   * field). Thus, different steps (length) can be dynamically
   * generated here. The velocity vector will correctly point into the
   * direction as calculated by the algorithm for any steplength, or
   * intermediate position.
   *
   * One complete leap-frog step is
   * @f{eqnarray*}{
   * \vec{x}(t_{i+0.5}) &=& \vec{x}(t_{i}) + \vec{v(t_{i})} * \Delta t / 2 \\
   * \vec{v}(t_{i+1}) &=& \vec{v}(t_{i}) + \vec{v}(t_{i})\times\vec{B}(x_{i}, t_{i}) *
   * \Delta t \\
   * \vec{x}(t_{i+1}) &=& \vec{x}(t_{i+0.5}) + \vec{v}(t_{i+1}) * \Delta t /2 \\
   * @f}
   *
   * The volocity update has the characteristics @f$|\vec{v}(t_{i+1})|>1@f$, thus final
   * velocities are renormalised. The full leap-frog steplength is thus
   * @f[ L = |\vec{v}(t_{i+1})| \cdot \Delta t / 2  + |\vec{v}(t_{i+1})| \cdot \Delta t /
   * 2 @f]
   */

  class LeapFrogTrajectory : public BaseTrajectory {

  public:
    LeapFrogTrajectory() = delete;
    LeapFrogTrajectory(LeapFrogTrajectory const&) = default;
    LeapFrogTrajectory(LeapFrogTrajectory&&) = default;
    LeapFrogTrajectory& operator=(LeapFrogTrajectory const&) = delete;

    LeapFrogTrajectory(Point const& pos, VelocityVector const& initialVelocity,
                       MagneticFieldVector const& Bfield,
                       decltype(1 / (tesla * second)) const k,
                       TimeType const timeStep) // leap-from total length
        : initialPosition_(pos)
        , initialVelocity_(initialVelocity)
        , initialDirection_(initialVelocity.normalized())
        , magneticfield_(Bfield)
        , k_(k)
        , timeStep_(timeStep) {}

    Line getLine() const;

    Point getPosition(double const u) const;

    VelocityVector getVelocity(double const u) const;

    DirectionVector getDirection(double const u) const;

    ///! duration along potentially bend trajectory
    TimeType getDuration(double const u = 1) const;

    ///! time at the start (u=0) or at the end (u=1) of the track of a particle
    template <typename Particle>
    TimeType getTime(Particle const& particle, double const u) const;

    ///! total length along potentially bend trajectory
    LengthType getLength(double const u = 1) const;

    ///! set new duration along potentially bend trajectory.
    void setLength(LengthType const limit);

    ///! set new duration along potentially bend trajectory.
    //   Scale other properties by "limit/timeLength_"
    void setDuration(TimeType const limit);

  private:
    Point initialPosition_;
    VelocityVector initialVelocity_;
    DirectionVector initialDirection_;
    MagneticFieldVector magneticfield_;
    decltype(1 / (tesla * second)) k_;
    TimeType timeStep_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/LeapFrogTrajectory.inl>
