/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/random/RNGManager.hpp>

#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  /**
   * A simple model for elastic hadronic interactions based on the formulas
   * in Gaisser, Engel, Resconi, Cosmic Rays and Particle Physics (Cambridge Univ. Press,
   * 2016), section 4.2 and Donnachie, Landshoff, Phys. Lett. B 296, 227 (1992)
   *
   * Currently only \f$p\f$ projectiles are supported and cross-sections are assumed to be
   * \f$pp\f$-like even for nuclei.
   *
   * \todo add unit test
   */
  class HadronicElasticInteraction
      : public InteractionProcess<HadronicElasticInteraction> {
  private:
    using SquaredHEPEnergyType = decltype(HEPEnergyType() * HEPEnergyType());

    using eV2 = decltype(square(electronvolt));
    using inveV2 = decltype(1 / square(electronvolt));

    inveV2 B(eV2 s) const;
    CrossSectionType getCrossSection(SquaredHEPEnergyType s) const;

  public:
    HadronicElasticInteraction( // x & y values taken from DL for pp collisions
        CrossSectionType x = 0.0217 * barn, CrossSectionType y = 0.05608 * barn);

    template <typename TParticle>
    GrammageType getInteractionLength(TParticle const& p);

    template <typename TParticle>
    ProcessReturn doInteraction(TParticle&);

  private:
    CrossSectionType parX_, parY_;

    static double constexpr gfEpsilon = 0.0808;
    static double constexpr gfEta = 0.4525;
    // Froissart-Martin is not violated up for sqrt s < 10^32 eV with these values [DL].

    default_prng_type& RNG_ =
        RNGManager<>::getInstance().getRandomStream("HadronicElasticModel");
  };

} // namespace corsika

#include <corsika/detail/modules/HadronicElasticModel.inl>
