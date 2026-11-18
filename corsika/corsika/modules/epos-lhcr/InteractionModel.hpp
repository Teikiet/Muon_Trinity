/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <tuple>
#include <set>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/random/RNGManager.hpp>

namespace corsika::EPOS_LHCR {

  class InteractionModel {

  public:
    InteractionModel(std::set<Code> = {}, std::string const& dataPath = "",
                     bool const epos_printout_on = false);
    ~InteractionModel();

    /**
     * Returns production and elastic cross section for hadrons in epos.
     *
     * Allowed targets are: nuclei or single nucleons (p,n,hydrogen). This routine
     * calculates the cross sections from scratch.
     *
     * Note: **Very slow!**, use tabulation for any performance application.
     *
     * @param corsikaId - PID of beam particle,
     * @param targetId - PID of target particle
     * @param sqrtS - center-of-mass energy.
     */
    std::tuple<CrossSectionType, CrossSectionType> calcCrossSectionCoM(
        Code const corsikaId, int const, int const, Code const targetId, int const,
        int const, HEPEnergyType const sqrtS) const;

    /**
     * Returns production and elastic cross section for hadrons in epos by reading
     * pre-calculated tables from epos.
     */
    std::tuple<CrossSectionType, CrossSectionType> readCrossSectionTableLab(
        Code const, int const, int const, Code const, HEPEnergyType const) const;

    /**
     * Returns production and elastic cross section.
     *
     * Allowed configurations are
     * hadron-nucleon, hadron-nucleus and nucleus-nucleus. Inputs are particle id's mass
     * and charge numbers and total energy in the lab.
     */
    std::tuple<CrossSectionType, CrossSectionType> getCrossSectionInelEla(
        Code const projectileId, Code const targetId, FourMomentum const& projectileP4,
        FourMomentum const& targetP4) const;

    /**
     * Checks validity of projectile, target and energy combination.
     *
     * EPOSLHC only accepts nuclei with X<=A<=Y as targets, or protons aka Hydrogen or
     * neutrons (p,n == nucleon).
     */
    bool isValid(Code const projectileId, Code const targetId,
                 HEPEnergyType const sqrtS) const;

    /**
     * Get the inelatic/production cross section.
     *
     * @param projectileId
     * @param targetId
     * @param projectileP4
     * @param targetP4
     * @return CrossSectionType
     */
    CrossSectionType getCrossSection(Code const projectileId, Code const targetId,
                                     FourMomentum const& projectileP4,
                                     FourMomentum const& targetP4) const {
      return std::get<0>(
          getCrossSectionInelEla(projectileId, targetId, projectileP4, targetP4));
    }

    /**
     * Calculate one hadron-hadron interaction.
     */
    template <typename TSecondaries>
    void doInteraction(TSecondaries&, Code const projectileId, Code const targetId,
                       FourMomentum const& projectileP4, FourMomentum const& targetP4);

    void initializeEventCoM(Code const, int const, int const, Code const, int const,
                            int const, HEPEnergyType const) const;
    void initializeEventLab(Code const, int const, int const, Code const, int const,
                            int const, HEPEnergyType const) const;
    void configureParticles(Code const, int const, int const, Code const, int const,
                            int const) const;

  private:
    // initialize and setParticlesStable are private since they can only be called once at
    // the beginning and are already called in the constructor!
    void initialize() const;
    void setParticleListStable(std::set<Code>) const;
    inline static bool isInitialized_ = false;

    std::string data_path_;
    unsigned int count_ = 0;
    bool epos_listing_;

    default_prng_type& RNG_ = RNGManager<>::getInstance().getRandomStream("epos-lhcr");
    std::shared_ptr<spdlog::logger> logger_ = get_logger("corsika_epos_Interaction");
    HEPEnergyType const minEnergyCoM_ = 6 * 1e9 * electronvolt;
    HEPEnergyType const maxEnergyCoM_ = 2.e6 * 1e9 * electronvolt;
    static Code constexpr maxNucleus_ = Code::Lead;
  };

} // namespace corsika::EPOS_LHCR

#include <corsika/detail/modules/epos-lhcr/InteractionModel.inl>
