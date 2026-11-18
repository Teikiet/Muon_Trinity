/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/modules/LongitudinalProfile.hpp>

#include <cmath>
#include <iomanip>
#include <limits>
#include <utility>

namespace corsika {

  template <typename TOutput>
  template <typename... TArgs>
  inline LongitudinalProfile<TOutput>::LongitudinalProfile(TArgs&&... args)
      : TOutput(std::forward<TArgs>(args)...) {}

  template <typename TOutput>
  template <typename TParticle>
  inline ProcessReturn LongitudinalProfile<TOutput>::doContinuous(
      Step<TParticle> const& step, bool const) {

    auto const pid = step.getParticlePre().getPID();
    this->write(step.getPositionPre(), step.getPositionPost(), pid,
                step.getParticlePre().getWeight()); // weight hardcoded so far
    return ProcessReturn::Ok;
  }

  template <typename TOutput>
  inline YAML::Node LongitudinalProfile<TOutput>::getConfig() const {
    YAML::Node node;
    node["type"] = "LongitudinalProfile";

    return node;
  }

} // namespace corsika
