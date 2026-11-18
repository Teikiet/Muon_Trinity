/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/qgsjetII/InteractionModel.hpp>

#include <corsika/framework/process/InteractionProcess.hpp>

/**
 * @file QGSJetII.hpp
 *
 * Includes all the parts of the QGSJetII model. Defines the InteractionProcess<TModel>
 * classes needed for the ProcessSequence.
 */

namespace corsika::qgsjetII {
  /**
   * @brief qgsjetII::Interaction is the process for ProcessSequence.
   *
   * The qgsjetII::InteractionModel is wrapped as an InteractionProcess here in order
   * to provide all the functions for ProcessSequence.
   */
  class Interaction : public InteractionModel, public InteractionProcess<Interaction> {};
} // namespace corsika::qgsjetII
