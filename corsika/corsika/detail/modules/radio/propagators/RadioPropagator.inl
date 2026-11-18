/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/propagators/RadioPropagator.hpp>

namespace corsika {

  template <typename TImpl, typename TEnvironment>
  inline RadioPropagator<TImpl, TEnvironment>::RadioPropagator(TEnvironment const& env)
      : env_(env) {}

} // namespace corsika
