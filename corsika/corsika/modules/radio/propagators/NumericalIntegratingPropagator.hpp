/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/media/Environment.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/modules/radio/propagators/RadioPropagator.hpp>

namespace corsika {

  /**
   * This class implements a basic propagator that uses
   * the straight-line (vector) between the particle
   * location and the observer as the trajectory.
   * To calculate the time delay of the signal, a basic
   * numerical integration scheme based on Simpson's rule
   * takes place. This propagator is slow and not
   * recommended for big showers simulations.
   *
   * This is what is used in ZHAireS and CoREAS in C7.
   */
  template <typename TEnvironment>
  class NumericalIntegratingPropagator final
      : public RadioPropagator<NumericalIntegratingPropagator<TEnvironment>,
                               TEnvironment> {

    using Base =
        RadioPropagator<NumericalIntegratingPropagator<TEnvironment>, TEnvironment>;
    using SignalPathCollection = typename Base::SignalPathCollection;

  public:
    /**
     * Construct a new StraightPropagator with a given environment.
     *
     */
    NumericalIntegratingPropagator(TEnvironment const& env, LengthType const stepsize);

    /**
     * Return the collection of paths from `start` to `end`.
     * or from 'source' which is the emission point to 'destination'
     * which is the location of the observer
     */
    template <typename Particle>
    SignalPathCollection propagate(Particle const& particle, Point const& source,
                                   Point const& destination) const;

  private:
    LengthType const stepsize_;

  }; // End: StraightPropagator

  template <typename TEnvironment>
  NumericalIntegratingPropagator<TEnvironment>
  make_numerical_integrating_radio_propagator(TEnvironment const& env,
                                              LengthType const stepsize) {
    return NumericalIntegratingPropagator<TEnvironment>(env, stepsize);
  }

} // namespace corsika

#include <corsika/detail/modules/radio/propagators/NumericalIntegratingPropagator.inl>
