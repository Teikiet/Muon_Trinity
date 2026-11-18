/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/pythia8/Decay.hpp>
#include <corsika/modules/pythia8/NeutrinoInteraction.hpp>
#include <corsika/modules/pythia8/InteractionModel.hpp>

namespace corsika::pythia8 {
  /**
   * pythia8::Interaction is the process for ProcessSequence.
   *
   * The pythia8::InteractionModel is wrapped as an InteractionProcess here in order
   * to provide all the functions for ProcessSequence.
   */
  class Interaction : public InteractionModel, public InteractionProcess<Interaction> {
  public:
    Interaction(std::set<Code> const& stableList = {})
        : InteractionModel{stableList} {};
  };
} // namespace corsika::pythia8
