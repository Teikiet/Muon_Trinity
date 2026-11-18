/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/StraightTrajectory.hpp>

namespace corsika {

  struct DeltaParticleState {
  public:
    DeltaParticleState(HEPEnergyType dEkin, TimeType dT, LengthVector const& ds,
                       DirectionVector const& du)
        : delta_Ekin_{dEkin}
        , delta_time_{dT}
        , displacement_{ds}
        , delta_direction_{du} {}

    HEPEnergyType delta_Ekin_ = HEPEnergyType::zero();
    TimeType delta_time_ = TimeType::zero();
    LengthVector displacement_;
    DirectionVector delta_direction_;
  };

  template <typename TParticle>
  class Step {
  public:
    template <typename TTrajectory>
    Step(TParticle const& particle, TTrajectory const& track)
        : particlePreStep_{particle} //~ , track_{track}
        , diff_{HEPEnergyType::zero(), track.getDuration(1),
                track.getPosition(1) - particle.getPosition(),
                track.getDirection(1) - particle.getDirection()}

    {}

    HEPEnergyType const& add_dEkin(HEPEnergyType dEkin) {
      diff_.delta_Ekin_ += dEkin;
      return diff_.delta_Ekin_;
    }

    TimeType const& add_dt(TimeType dt) {
      diff_.delta_time_ += dt;
      return diff_.delta_time_;
    }

    LengthVector const& add_displacement(LengthVector const& dis) {
      diff_.displacement_ += dis;
      return diff_.displacement_;
    }

    DirectionVector const& add_dU(DirectionVector const& du) {
      diff_.delta_direction_ += du;
      return diff_.delta_direction_;
    }

    // getters for difference

    HEPEnergyType getDiffEkin() const { return diff_.delta_Ekin_; }

    TimeType getDiffT() const { return diff_.delta_time_; };

    DirectionVector const& getDiffDirection() const { return diff_.delta_direction_; }

    LengthVector const& getDisplacement() const { return diff_.displacement_; }

    //! alias for getDisplacement()
    LengthVector const& getDiffPosition() const { return getDisplacement(); }

    // getters for absolute
    TParticle const& getParticlePre() const { return particlePreStep_; }

    HEPEnergyType getEkinPre() const { return getParticlePre().getKineticEnergy(); }

    HEPEnergyType getEkinPost() const { return getEkinPre() + getDiffEkin(); }

    TimeType getTimePre() const { return getParticlePre().getTime(); }

    TimeType getTimePost() const { return getTimePre() + getDiffT(); }

    DirectionVector const getDirectionPre() const {
      return getParticlePre().getDirection();
    }

    VelocityVector getVelocityVector() const { return getDisplacement() / getDiffT(); }

    StraightTrajectory getStraightTrack() const {
      Line const line(getPositionPre(), getVelocityVector());
      StraightTrajectory track(line, getDiffT());
      return track;
    }

    DirectionVector getDirectionPost() const {
      return (getDirectionPre() + getDiffDirection()).normalized();
    }

    Point const& getPositionPre() const { return getParticlePre().getPosition(); }

    Point getPositionPost() const { return getPositionPre() + getDisplacement(); }

  private:
    TParticle const& particlePreStep_;
    //~ TTrajectory const& track_; // TODO: perhaps remove
    DeltaParticleState diff_;
  };

} // namespace corsika
