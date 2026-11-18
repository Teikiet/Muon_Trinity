/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <Pythia8/Pythia.h>

#include <corsika/framework/random/RNGManager.hpp>

namespace corsika::pythia8 {

  class Random : public Pythia8::RndmEngine {
    double flat();

  private:
    std::uniform_real_distribution<double> Dist_;
    default_prng_type& RNG_ = RNGManager<>::getInstance().getRandomStream("pythia");
  };

} // namespace corsika::pythia8

#include <corsika/detail/modules/pythia8/Random.inl>
