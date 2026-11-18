/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/RadioProcess.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  /**
   * A concrete implementation of the ZHS algorithm.
   */
  template <typename TRadioDetector, typename TPropagator>
  class ZHS final : public RadioProcess<TRadioDetector, ZHS<TRadioDetector, TPropagator>,
                                        TPropagator> {

  public:
    // an identifier for which algorithm was used
    static constexpr auto algorithm = "ZHS";

    /**
     * Construct a new ZHS instance.
     *
     * This forwards the detector and other arguments to
     * the RadioProcess parent.
     *
     */
    template <typename... TArgs>
    ZHS(TRadioDetector& detector, TArgs&&... args)
        : RadioProcess<TRadioDetector, ZHS, TPropagator>(detector, args...) {}

    /**
     * Simulate the radio emission from a particle across a track.
     *
     * This must be provided by the TRadioImpl.
     *
     * @param particle    The current particle.
     * @param track       The current track.
     *
     */
    template <typename Particle>
    ProcessReturn simulate(Step<Particle> const& step);
    static constexpr auto emConstant_{1.0 / (4.0 * M_PI) / (constants::epsilonZero) /
                                      constants::c};

  private:
    using Base =
        RadioProcess<TRadioDetector, ZHS<TRadioDetector, TPropagator>, TPropagator>;
    using Base::observers_;

  }; // END: class ZHS

} // namespace corsika

#include <corsika/detail/modules/radio/ZHS.inl>
