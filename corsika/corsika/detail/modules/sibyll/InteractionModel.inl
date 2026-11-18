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

  inline InteractionModel::InteractionModel(std::set<Code> const& nuccomp,
                                            std::set<Code> const& list)
      : hadronSibyll_{list}
      , nuclearSibyll_{hadronSibyll_, nuccomp} {}

  inline HadronInteractionModel& InteractionModel::getHadronInteractionModel() {
    return hadronSibyll_;
  }

  inline HadronInteractionModel const& InteractionModel::getHadronInteractionModel()
      const {
    return hadronSibyll_;
  }

  inline typename InteractionModel::nuclear_model_type&
  InteractionModel::getNuclearInteractionModel() {
    return nuclearSibyll_;
  }

  inline typename InteractionModel::nuclear_model_type const&
  InteractionModel::getNuclearInteractionModel() const {
    return nuclearSibyll_;
  }

  inline CrossSectionType InteractionModel::getCrossSection(
      Code projCode, Code targetCode, FourMomentum const& proj4mom,
      FourMomentum const& target4mom) const {
    if (is_nucleus(projCode))
      return getNuclearInteractionModel().getCrossSection(projCode, targetCode, proj4mom,
                                                          target4mom);
    else
      return getHadronInteractionModel().getCrossSection(projCode, targetCode, proj4mom,
                                                         target4mom);
  }

  inline std::tuple<CrossSectionType, CrossSectionType>
  InteractionModel::getCrossSectionInelEla(Code projCode, Code targetCode,
                                           FourMomentum const& proj4mom,
                                           FourMomentum const& target4mom) const {
    if (is_nucleus(projCode))
      return {getNuclearInteractionModel().getCrossSection(projCode, targetCode, proj4mom,
                                                           target4mom),
              CrossSectionType::zero()};
    else
      return getHadronInteractionModel().getCrossSectionInelEla(projCode, targetCode,
                                                                proj4mom, target4mom);
  }

  template <typename TSecondaries>
  inline void InteractionModel::doInteraction(TSecondaries& view, Code projCode,
                                              Code targetCode,
                                              FourMomentum const& proj4mom,
                                              FourMomentum const& target4mom) {
    if (is_nucleus(projCode))
      return getNuclearInteractionModel().doInteraction(view, projCode, targetCode,
                                                        proj4mom, target4mom);
    else
      return getHadronInteractionModel().doInteraction(view, projCode, targetCode,
                                                       proj4mom, target4mom);
  }
} // namespace corsika::sibyll
