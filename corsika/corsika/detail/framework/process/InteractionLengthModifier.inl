/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <utility>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  template <class TUnderlyingProcess>
  inline InteractionLengthModifier<TUnderlyingProcess>::InteractionLengthModifier(
      TUnderlyingProcess&& process,
      std::function<InteractionLengthModifier::functor_signature> modifier)
      : process_{std::move(process)}
      , modifier_{std::move(modifier)} {}

  template <class TUnderlyingProcess>
  template <typename TSecondaryView>
  inline void InteractionLengthModifier<TUnderlyingProcess>::doInteraction(
      TSecondaryView& view) {
    process_.doInteraction(view);
  }

  template <class TUnderlyingProcess>
  template <typename TParticle>
  inline GrammageType InteractionLengthModifier<TUnderlyingProcess>::getInteractionLength(
      TParticle const& particle) {
    GrammageType const original = process_.getInteractionLength(particle);
    Code const pid = particle.getPID();
    HEPEnergyType const energy = particle.getEnergy();

    return modifier_(original, pid, energy);
  }

  template <class TUnderlyingProcess>
  inline TUnderlyingProcess const&
  InteractionLengthModifier<TUnderlyingProcess>::getProcess() const {
    return process_;
  }

  template <class TUnderlyingProcess>
  inline TUnderlyingProcess& InteractionLengthModifier<TUnderlyingProcess>::getProcess() {
    return process_;
  }

} // namespace corsika
