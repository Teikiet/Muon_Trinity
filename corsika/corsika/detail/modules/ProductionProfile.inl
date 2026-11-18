/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/modules/ProductionProfile.hpp>

#include <cmath>
#include <iomanip>
#include <limits>
#include <utility>

namespace corsika {

  template <typename TOutput>
  template <typename... TArgs>
  inline ProductionProfile<TOutput>::ProductionProfile(TArgs&&... args)
      : TOutput(std::forward<TArgs>(args)...) {}

  template <typename TOutput>
  template <typename TStackView>
  inline void ProductionProfile<TOutput>::doSecondaries(TStackView& vS) {
    auto projectile = vS.getProjectile();
    auto particle = vS.begin();
    while (particle != vS.end()) {
      Code pid = particle.getPID();
      if (is_muon(pid)) {
        this->write(particle.getPosition(), projectile.getPID(), particle.getWeight());
        CORSIKA_LOG_TRACE("muon produced from {}: Ekin={} GeV, weight = {}",
                          projectile.getPID(), particle.getKineticEnergy() / 1_GeV,
                          particle.getWeight());
        count_++;
      }
      ++particle; // next entry in SecondaryView
    }
  }

  template <typename TOutput>
  inline YAML::Node ProductionProfile<TOutput>::getConfig() const {
    YAML::Node node;
    node["type"] = "ProductionProfile";

    return node;
  }

} // namespace corsika
