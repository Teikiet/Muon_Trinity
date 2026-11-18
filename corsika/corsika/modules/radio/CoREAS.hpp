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

  template <typename TRadioDetector, typename TPropagator>
  class CoREAS final
      : public RadioProcess<TRadioDetector, CoREAS<TRadioDetector, TPropagator>,
                            TPropagator> {

  public:
    // an identifier for which algorithm was used
    static constexpr auto algorithm = "CoREAS";

    /**
     * Construct a new CoREAS instance.
     *
     * This forwards the detector and other arguments to
     * the RadioProcess parent.
     *
     */
    template <typename... TArgs>
    CoREAS(TRadioDetector& detector, TArgs&&... args)
        : RadioProcess<TRadioDetector, CoREAS, TPropagator>(detector, args...) {}

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

    using Base =
        RadioProcess<TRadioDetector, CoREAS<TRadioDetector, TPropagator>, TPropagator>;
    using Base::observers_;

  }; // end of class CoREAS

} // namespace corsika

#include <corsika/detail/modules/radio/CoREAS.inl>
