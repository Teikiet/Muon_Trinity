/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/fluka/InteractionModel.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>

/**
 * @file FLUKA.hpp
 *
 * Includes all the parts of the FLUKA model. Defines the InteractionProcess<TModel>
 * classes needed for the ProcessSequence.
 */

namespace corsika::fluka {
  /**
   * fluka::Interaction is the process for ProcessSequence.
   *
   * The fluka::InteractionModel is wrapped as an InteractionProcess here in order
   * to provide all the functions for ProcessSequence.
   */
  class Interaction : public fluka::InteractionModel,
                      public corsika::InteractionProcess<Interaction> {
  public:
    template <typename TEnvironment>
    Interaction(TEnvironment const& env)
        : fluka::InteractionModel{env} {}
  };
} // namespace corsika::fluka
