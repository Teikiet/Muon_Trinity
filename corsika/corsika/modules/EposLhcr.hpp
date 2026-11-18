/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/epos-lhcr/ParticleConversion.hpp>
#include <corsika/modules/epos-lhcr/InteractionModel.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>

/**
 * @file EposLhcr.hpp
 *
 * Includes all the parts of the EPOS-LHCR model. Defines the InteractionProcess<TModel>
 * classes needed for the ProcessSequence.
 */

namespace corsika::EPOS_LHCR {
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
} // namespace corsika::EPOS_LHCR
