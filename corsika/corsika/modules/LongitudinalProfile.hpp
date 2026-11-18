/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/ShowerAxis.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Step.hpp>

#include <corsika/modules/writers/LongitudinalProfileWriterParquet.hpp>

#include <array>
#include <fstream>
#include <limits>
#include <string>

namespace corsika {

  /**
   * \class LongitudinalProfile
   *
   * \todo test missing
   *
   * is a ContinuousProcess, which is constructed from an environment::ShowerAxis
   * object, and a dX in units of g/cm2  (GrammageType).
   *
   * LongitudinalProfile does then convert each single Track of the
   * simulation into a projected grammage range and counts for
   * different particle species when they cross dX (default: 10g/cm2)
   * boundaries.
   */

  template <typename TOutput>
  class LongitudinalProfile : public ContinuousProcess<LongitudinalProfile<TOutput>>,
                              public TOutput {

  public:
    template <typename... TArgs>
    LongitudinalProfile(TArgs&&... args);

    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle> const&, bool);

    template <typename TParticle, typename TTrack>
    LengthType getMaxStepLength(TParticle const&, TTrack const&) {
      return meter * std::numeric_limits<double>::infinity();
    }

    YAML::Node getConfig() const;
  };

} // namespace corsika

#include <corsika/detail/modules/LongitudinalProfile.inl>
