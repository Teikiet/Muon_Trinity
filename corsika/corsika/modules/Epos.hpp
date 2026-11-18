/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/epos/ParticleConversion.hpp>
#include <corsika/modules/epos/InteractionModel.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>

/**
 * @file Epos.hpp
 *
 * Includes all the parts of the EPOS model. Defines the InteractionProcess<TModel>
 * classes needed for the ProcessSequence.
 */

namespace corsika::epos {
  /**
   * epos::Interaction is the process for ProcessSequence.
   *
   * The epos::InteractionModel is wrapped as an InteractionProcess here in order
   * to provide all the functions for ProcessSequence.
   */
  class Interaction : public InteractionModel, public InteractionProcess<Interaction> {
  public:
    Interaction(std::set<Code> const& stableList)
        : InteractionModel{stableList} {};
  };
} // namespace corsika::epos
