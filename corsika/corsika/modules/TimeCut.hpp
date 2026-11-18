/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/core/Step.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

namespace corsika {

  /**
   * Simple TimeCut process, removes particles older than specified cut time.
   */
  class TimeCut : public ContinuousProcess<TimeCut> {

  public:
    TimeCut(TimeType const time);

    template <typename Particle>
    ProcessReturn doContinuous(
        Step<Particle> const& step,
        const bool limitFlag = false); // this is not used for TimeCut
    template <typename Particle, typename Track>
    LengthType getMaxStepLength(Particle const&, Track const&) {
      return meter * std::numeric_limits<double>::infinity();
    }

  private:
    TimeType time_;
  };

} // namespace corsika

#include <corsika/detail/modules/TimeCut.inl>
