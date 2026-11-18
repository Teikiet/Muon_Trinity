/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <vector>

#include <corsika/framework/geometry/Path.hpp>
#include <corsika/modules/radio/propagators/SignalPath.hpp>

namespace corsika {

  /**
   * Radio propagators are used to calculate the propagation
   * paths from particles to observers. Any class that wants
   * to be used as a RadioPropagator must implement the
   * following methods:
   *
   *  SignalPathCollection Propagate(Particle const& particle,
   *                                 Point const& source,
                                     Point const& destination) const
   */
  template <typename TImpl, typename TEnvironment>
  class RadioPropagator {

  protected:
    // Since we typically know roughly how many paths will
    // be computed for a given propagator, we use a std::vector here.
    using SignalPathCollection = std::vector<SignalPath> const;

    TEnvironment const& env_; ///< The environment.

  public:
    /**
     * Construct a new RadioPropagator instance.
     */
    RadioPropagator(TEnvironment const& env);

  }; // class RadioPropagator

} // namespace corsika

#include <corsika/detail/modules/radio/propagators/RadioPropagator.inl>
