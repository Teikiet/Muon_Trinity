/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/sibyll/ParticleConversion.hpp>
#include <corsika/modules/sibyll/HadronInteractionModel.hpp>
#include <corsika/modules/sibyll/Decay.hpp>
#include <corsika/modules/sibyll/NuclearInteractionModel.hpp>

#include <corsika/modules/sibyll/InteractionModel.hpp>

#include <corsika/framework/process/InteractionProcess.hpp>

/**
 * @file Sibyll.hpp
 *
 * Includes all the parts of the Sibyll model. Defines the InteractionProcess<TModel>
 * classes needed for the ProcessSequence.
 */

namespace corsika::sibyll {
  /**
   * @brief sibyll::Interaction is the process for ProcessSequence.
   *
   * The sibyll::InteractionModel is wrapped as an InteractionProcess here in order
   * to provide all the functions for ProcessSequence.
   */
  struct Interaction : public InteractionModel, public InteractionProcess<Interaction> {
    Interaction(std::set<Code> const& nuccomp, std::set<Code> const& stablehad)
        : InteractionModel{nuccomp, stablehad} {}
  };

  /**
   * @brief sibyll::NuclearInteraction is the process for ProcessSequence.
   *
   * The sibyll::NuclearInteractionModel is wrapped as an InteractionProcess here in order
   * to provide all the functions for ProcessSequence.
   */
  template <class TNucleonModel>
  class NuclearInteraction
      : public NuclearInteractionModel<TNucleonModel>,
        public InteractionProcess<NuclearInteraction<TNucleonModel>> {
  public:
    NuclearInteraction(TNucleonModel& model, std::set<Code> const& nuccomp)
        : NuclearInteractionModel<TNucleonModel>{model, nuccomp} {}
  };

} // namespace corsika::sibyll
