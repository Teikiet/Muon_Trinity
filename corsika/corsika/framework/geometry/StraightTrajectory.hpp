/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {

  /**
   *
   * This implements a straight trajectory between two points.
   *
   * Note: so far it is assumed that the speed (d|vec{r}|/dt) between
   * start and end does not change and is constant for the entire
   * Trajectory.
   *
   **/

  class StraightTrajectory : public BaseTrajectory {

  public:
    StraightTrajectory() = delete;
    StraightTrajectory(StraightTrajectory const&) = default;
    StraightTrajectory(StraightTrajectory&&) = default;
    StraightTrajectory& operator=(StraightTrajectory const&) = delete;

    /**
     * \param theLine The geometric Line object that represents a straight-line
     * connection
     *
     * \param timeLength The time duration to traverse the straight trajectory
     * in units of TimeType
     */
    StraightTrajectory(Line const& theLine, TimeType timeLength)
        : line_(theLine)
        , timeLength_(timeLength)
        , timeStep_(timeLength)
        , initialVelocity_(theLine.getVelocity(TimeType::zero()))
        , finalVelocity_(theLine.getVelocity(timeLength)) {}

    /**
     * \param theLine The geometric Line object that represents a straight-line
     * connection
     *
     * \param timeLength The time duration to traverse the straight trajectory
     * in units of TimeType
     *
     * \param timeStep Time duration to folow eventually curved
     * trajectory in units of TimesType
     *
     * \param initialV Initial velocity vector at
     * start of trajectory
     *
     * \param finalV Final velocity vector at start of trajectory
     */
    StraightTrajectory(Line const& theLine,
                       TimeType const timeLength, // length of theLine (straight)
                       TimeType const timeStep,   // length of bend step (curved)
                       VelocityVector const& initialV, VelocityVector const& finalV)
        : line_(theLine)
        , timeLength_(timeLength)
        , timeStep_(timeStep)
        , initialVelocity_(initialV)
        , finalVelocity_(finalV) {}

    Line const& getLine() const { return line_; }

    Point getPosition(double const u) const { return line_.getPosition(timeLength_ * u); }

    VelocityVector getVelocity(double const u) const;

    DirectionVector getDirection(double const u) const {
      return getVelocity(u).normalized();
    }

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

  protected:
    ///! total length along straight trajectory
    LengthType getDistance(double const u) const;

    void setFinalVelocity(VelocityVector const& v) { finalVelocity_ = v; }

  private:
    Line line_;
    TimeType timeLength_; ///! length of straight step (shortest connecting line)
    TimeType timeStep_;   ///! length of bend step (curved)
    VelocityVector initialVelocity_;
    VelocityVector finalVelocity_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/StraightTrajectory.inl>
