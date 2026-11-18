/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/modules/sibyll/HadronInteractionModel.hpp>
#include <corsika/modules/sibyll/NuclearInteractionModel.hpp>

namespace corsika::sibyll {
  /**
   * This class combines sibyll::HadronInteractionModel, which can only handle hadron
   * projectiles, and sibyll::NuclearInteractionModel, which can handle only nucleus
   * projectile, into a single process. The getCrossSection() and doInteraction() methods
   * forward to the underlying models, depending on the projectile type.
   */
  class InteractionModel {
  public:
    using nuclear_model_type = NuclearInteractionModel<HadronInteractionModel>;

    InteractionModel(std::set<Code> const&, std::set<Code> const&);

    CrossSectionType getCrossSection(Code, Code, FourMomentum const&,
                                     FourMomentum const&) const;

    std::tuple<CrossSectionType, CrossSectionType> getCrossSectionInelEla(
        Code, Code, FourMomentum const&, FourMomentum const&) const;

    template <typename TSecondaries>
    void doInteraction(TSecondaries&, Code, Code, FourMomentum const&,
                       FourMomentum const&);

    HadronInteractionModel& getHadronInteractionModel();
    HadronInteractionModel const& getHadronInteractionModel() const;
    nuclear_model_type& getNuclearInteractionModel();
    nuclear_model_type const& getNuclearInteractionModel() const;

  private:
    HadronInteractionModel hadronSibyll_;
    nuclear_model_type nuclearSibyll_;
  };
} // namespace corsika::sibyll

#include <corsika/detail/modules/sibyll/InteractionModel.inl>
