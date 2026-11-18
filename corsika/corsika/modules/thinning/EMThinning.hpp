/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <random>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/process/SecondariesProcess.hpp>

namespace corsika {

  /**
   * This process implements thinning for EM splitting processes (1 -> 2).
   */

  class EMThinning : public SecondariesProcess<EMThinning> {
  public:
    /**
     * Construct a new EMThinning process.
     *
     * @param threshold: thinning applied below this energy
     * @param maxWeight: maximum allowed weight
     */
    EMThinning(HEPEnergyType threshold, double maxWeight,
               bool const eraseParticles = true);

    /**
     * Apply thinning to secondaries. Only EM primaries with two EM secondaries are
     * considered.
     *
     * If the maximum weight is still out of reach, Hillas thinning is applied (i.e.
     * one of the two secondaries is kept, the other one discarded). If the acceptance
     * probabilities would lead to a weight factor exceeding the maximum weight, we resort
     * to statistical thinning (i.e. the secondaries are kept/discared randomly each on is
     * own). In that case, acceptance probabilities can be assigned without constraints
     * (sum does not need to be 1) and we increase the acceptance probability such that
     * the maximum weight is not exceeded.
     */
    template <typename TStackView>
    void doSecondaries(TStackView&);

  private:
    default_prng_type& rng_ = RNGManager<>::getInstance().getRandomStream("thinning");
    std::uniform_real_distribution<double> uniform_{};
    HEPEnergyType const threshold_;
    double const maxWeight_;
    bool const eraseParticles_;
  };
} // namespace corsika

#include <corsika/detail/modules/thinning/EMThinning.inl>
