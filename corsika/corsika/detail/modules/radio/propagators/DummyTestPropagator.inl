/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/propagators/DummyTestPropagator.hpp>

namespace corsika {

  template <typename TEnvironment>
  inline DummyTestPropagator<TEnvironment>::DummyTestPropagator(TEnvironment const& env)
      : RadioPropagator<DummyTestPropagator, TEnvironment>(env) {
    rindex.reserve(2);
  }

  template <typename TEnvironment>
  template <typename Particle>
  inline typename DummyTestPropagator<TEnvironment>::SignalPathCollection
  DummyTestPropagator<TEnvironment>::propagate(Particle const& particle,
                                               Point const& source,
                                               Point const& destination) {

    /**
     * This is the simplest case of straight propagator
     * where no integration takes place.
     * This can be used for fast tests and checks of the radio module.
     *
     */

    // these are used for the direction of emission and reception of signal at the
    // observer
    auto const emit_{(destination - source).normalized()};
    auto const receive_{-emit_};

    // the geometrical distance from the point of emission to an observer
    auto const distance_{(destination - source).getNorm()};

    // get the universe for this environment
    auto const* const universe{Base::env_.getUniverse().get()};

    // clear the refractive index vector and points deque for this signal propagation.
    rindex.clear();
    points.clear();

    // get and store the refractive index of the first point 'source'.
    auto const* const nodeSource{particle.getNode()};
    auto const ri_source{nodeSource->getModelProperties().getRefractiveIndex(source)};
    rindex.push_back(ri_source);
    points.push_back(source);

    // add the refractive index of last point 'destination' and store it.
    auto const* const node{universe->getContainingNode(destination)};
    auto const ri_destination{node->getModelProperties().getRefractiveIndex(destination)};
    rindex.push_back(ri_destination);
    points.push_back(destination);

    // compute the average refractive index.
    auto const averageRefractiveIndex_ = (ri_source + ri_destination) * 0.5;

    // compute the total time delay.
    TimeType const time = averageRefractiveIndex_ * (distance_ / constants::c);

    return std::vector<SignalPath>(
        1, SignalPath(time, averageRefractiveIndex_, ri_source, ri_destination, emit_,
                      receive_, distance_, points));

  } // END: propagate()

} // namespace corsika
