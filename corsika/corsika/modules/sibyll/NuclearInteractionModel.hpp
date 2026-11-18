/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <set>

namespace corsika::sibyll {

  /**
   * The sibyll::NuclearInteractionModel provides the SIBYLL semi superposition model.
   *
   * can transform a proton-nucleus interaction model into a nucleus-nucleus interaction
   * model.
   *
   * pass list of nuclei in the environment
   *
   * @tparam TNucleonModel
   */
  template <class TNucleonModel>
  class NuclearInteractionModel {

  public:
    NuclearInteractionModel(TNucleonModel&, std::set<Code> const&);

    ~NuclearInteractionModel();

    bool constexpr isValid(Code const projectileId, Code const targetId,
                           HEPEnergyType const sqrtSnn) const;

    void initializeNuclearCrossSections(std::set<Code> const&);

    void printCrossSectionTable(Code) const;
    CrossSectionType readCrossSectionTable(int const, Code const,
                                           HEPEnergyType const) const;
    HEPEnergyType getMinEnergyPerNucleonCoM() const { return gMinEnergyPerNucleonCoM_; }
    HEPEnergyType getMaxEnergyPerNucleonCoM() const { return gMaxEnergyPerNucleonCoM_; }
    unsigned int constexpr getMaxNucleusAProjectile() const {
      return gMaxNucleusAProjectile_;
    }
    unsigned int constexpr getMaxNFragments() const { return gMaxNFragments_; }
    unsigned int constexpr getNEnergyBins() const { return gNEnBins_; }

    CrossSectionType getCrossSection(Code const, Code const,
                                     FourMomentum const& projectileP4,
                                     FourMomentum const& targetP4) const;

    template <typename TSecondaryView>
    void doInteraction(TSecondaryView&, Code const, Code const,
                       FourMomentum const& projectileP4, FourMomentum const& targetP4);

  private:
    int count_ = 0;

    TNucleonModel& hadronicInteraction_;
    std::map<Code, int> targetComponentsIndex_;
    std::shared_ptr<spdlog::logger> logger_ =
        get_logger("corsika_sibyll_NuclearInteractionModel");
    static unsigned int constexpr gNSample_ =
        500; // number of samples in MC estimation of cross section
    static unsigned int constexpr gMaxNucleusAProjectile_ = 56;
    static unsigned int constexpr gNEnBins_ = 6;
    static unsigned int constexpr gMaxNFragments_ = 60;
    // energy limits defined by table used for cross section in signuc.f
    // 10**1 GeV to 10**6 GeV
    static HEPEnergyType constexpr gMinEnergyPerNucleonCoM_ = 10. * 1e9 * electronvolt;
    static HEPEnergyType constexpr gMaxEnergyPerNucleonCoM_ = 1.e6 * 1e9 * electronvolt;
  };

} // namespace corsika::sibyll

#include <corsika/detail/modules/sibyll/NuclearInteractionModel.inl>
