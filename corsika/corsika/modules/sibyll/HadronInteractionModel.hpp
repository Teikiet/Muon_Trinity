/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/geometry/FourVector.hpp>

#include <set>
#include <tuple>

namespace corsika::sibyll {

  /**
   * @brief Provides the SIBYLL hadron-nucleus interaction model.
   *
   * This is a TModel argument for InteractionProcess<TModel>.
   */

  class HadronInteractionModel {

  public:
    HadronInteractionModel();
    /**
     * @brief construct hadron interaction model with SIBYLL
     *
     * @param list: list of particles that should be stable inside SIBYLL. all other
     * particles should decay (if they can)
     */
    HadronInteractionModel(std::set<Code> list);
    ~HadronInteractionModel();

    /**
     * @brief Set the Verbose flag.
     *
     * If flag is true, SIBYLL will printout additional secondary particle information
     * lists, etc.
     *
     * @param flag to switch.
     */
    void setVerbose(bool const flag);

    /**
     * @brief evaluated validity of collision system.
     *
     * sibyll only accepts nuclei with 4<=A<=18 as targets, or protons aka Hydrogen or
     * neutrons (p,n == nucleon).
     */
    bool constexpr isValid(Code const projectileId, Code const targetId,
                           HEPEnergyType const sqrtSnn) const;

    /**
     * Returns inelastic AND elastic cross sections.
     *
     * These cross sections must correspond to the process described in doInteraction
     * AND elastic scattering (sigma_tot = sigma_inel + sigma_el). Allowed targets are:
     * nuclei or single nucleons (p,n,hydrogen). This "InelEla" method is used since
     * Sibyll must be useful inside the NuclearInteraction model, which requires that.
     *
     * @param projectile is the Code of the projectile
     * @param target is the Code of the target
     * @param projectileP4: four-momentum of projectile
     * @param targetP4: four-momentum of target
     *
     * @return a tuple of: inelastic cross section, elastic cross section
     */
    std::tuple<CrossSectionType, CrossSectionType> getCrossSectionInelEla(
        Code const projectile, Code const target, FourMomentum const& projectileP4,
        FourMomentum const& targetP4) const;

    /**
     * Returns inelastic (production) cross section.
     *
     * This cross section must correspond to the process described in doInteraction.
     * Allowed targets are: nuclei or single nucleons (p,n,hydrogen).
     *
     * @param projectile is the Code of the projectile
     * @param target is the Code of the target
     * @param projectileP4: four-momentum of projectile
     * @param targetP4: four-momentum of target
     *
     * @return inelastic cross section
     * elastic cross section
     */
    CrossSectionType getCrossSection(Code const projectile, Code const target,
                                     FourMomentum const& projectileP4,
                                     FourMomentum const& targetP4) const {
      return std::get<0>(
          getCrossSectionInelEla(projectile, target, projectileP4, targetP4));
    }

    /**
     * In this function SIBYLL is called to produce one event. The
     * event is copied (and boosted) into the shower lab frame.
     *
     * @param view is the stack object for the secondaries
     * @param projectile is the Code of the projectile
     * @param target is the Code of the target
     * @param projectileP4: four-momentum of projectile
     * @param targetP4: four-momentum of target
     */

    template <typename TSecondaries>
    void doInteraction(TSecondaries& view, Code const projectile, Code const target,
                       FourMomentum const& projectileP4, FourMomentum const& targetP4);

  private:
    void setParticleListStable(std::set<Code>);
    void setAllParticlesUnstable();
    HEPEnergyType constexpr getMinEnergyCoM() const { return minEnergyCoM_; }
    HEPEnergyType constexpr getMaxEnergyCoM() const { return maxEnergyCoM_; }

    // hard model limits
    static HEPEnergyType constexpr minEnergyCoM_ = 10. * 1e9 * electronvolt;
    static HEPEnergyType constexpr maxEnergyCoM_ = 1.e6 * 1e9 * electronvolt;
    static unsigned int constexpr maxTargetMassNumber_ = 18;
    static unsigned int constexpr minNuclearTargetA_ = 4;

    std::shared_ptr<spdlog::logger> logger_ =
        get_logger("corsika_sibyll_HadronInteractionModel");

    // data members
    int count_ = 0;
    int nucCount_ = 0;
    bool sibyll_listing_;
    bool const internal_decays_;
    std::set<Code> stable_particles_;
  };

} // namespace corsika::sibyll

#include <corsika/detail/modules/sibyll/HadronInteractionModel.inl>
