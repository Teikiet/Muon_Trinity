/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <tuple>
#include <boost/filesystem/path.hpp>

#include <corsika/framework/utility/CorsikaData.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/modules/pythia8/Pythia8.hpp>

namespace corsika::pythia8 {

  using HEPEnergyTypeSqr = decltype(1_GeV * 1_GeV);

  /**
   * @brief Defines the interface to PYTHIA8. Configured for
   * DIS neutrino - nucleon interactions.
   *
   */

  class NeutrinoInteraction : public InteractionProcess<NeutrinoInteraction> {

  public:
    /**
     * Constructs the interface for PYTHIA8 neutrino interactions
     *
     * @param listOfStableParticles - Switch on Pythia internal decay for all hadrons
     * except those in the list
     * @param handleNC - Switch on/off neutral current interactions
     * @param handleCC - Switch on/off charged current interactions
     */
    NeutrinoInteraction(std::set<Code> const& listOfStableParticles,
                        bool const& handleNC = true, bool const& handleCC = true);
    ~NeutrinoInteraction();
    /**
     * Returns inelastic (production) cross section.
     *
     * This cross section must correspond to the process described in doInteraction.
     * Allowed targets are: nuclei or single nucleons (p,n,hydrogen).
     *
     * @param projectile is the Code of the projectile
     * @param target is the Code of the target
     * @param projectileP4 is the four momentum of the projectile
     * @param targetP4 is the four momentum of the target
     *
     * @return inelastic cross section
     */
    /* NOTE: the cross section for neutrino interactions is fixed to an arbitrary small
     * value. This is needed just so the interaction process is not skipped in the process
     * loop when the interaction is forced. If needed the exact value could be made
     * accessible. */
    CrossSectionType getCrossSection(Code const projectile, Code const target,
                                     FourMomentum const& projectileP4,
                                     FourMomentum const& targetP4) const {
      if (isValid(projectile, target, projectileP4, targetP4))
        return cross_section_;
      else
        return CrossSectionType::zero();
    };

    bool isValid(Code const projectileId, [[maybe_unused]] Code const targetId,
                 FourMomentum const& projectileP4, FourMomentum const& targetP4) const {
      auto const S = (projectileP4 + targetP4).getNormSqr();
      return is_neutrino(projectileId) &&
             (is_nucleus(targetId) || targetId == Code::Proton ||
              targetId == Code::Neutron) &&
             S >= minQ2_;
    };

    /**
     * In this function PYTHIA is called to produce one event. The
     * event is copied (and boosted) into the shower lab frame.
     */
    template <typename TView>
    void doInteraction(TView& output, Code const projectileId, Code const targetId,
                       FourMomentum const& projectileP4, FourMomentum const& targetP4);

  private:
    std::shared_ptr<spdlog::logger> logger_ =
        get_logger("corsika_pythia8_NeutrinoInteraction");
    CrossSectionType const cross_section_ = 4_yb;
    int count_ = 0;
    std::set<Code> const stable_particles_ = {};
    bool const handle_nc_ = true;
    bool const handle_cc_ = true;
    HEPEnergyTypeSqr minQ2_ = 25_GeV * 1_GeV;
    Pythia8::Pythia pythiaMain_;
  };

  /**
   * @brief Wraps NeutrinoInteraction for a single Pythia call per shower.
   *
   * NeutrinoInteraction::doInteraction re-inits Pythia on every call, which fails
   * on reuse. Use this wrapper when tracked secondary neutrinos must not trigger
   * a second interaction.
   */
  class SingleUseNeutrinoInteraction
      : public InteractionProcess<SingleUseNeutrinoInteraction> {

  public:
    SingleUseNeutrinoInteraction(std::set<Code> const& stableParticles, bool const NC,
                                 bool const CC)
        : inner_(stableParticles, NC, CC) {}

    void reset() { used_ = false; }

    CrossSectionType getCrossSection(Code const projectileId, Code const targetId,
                                     FourMomentum const& projectileP4,
                                     FourMomentum const& targetP4) const {
      if (used_) return CrossSectionType::zero();
      return inner_.getCrossSection(projectileId, targetId, projectileP4, targetP4);
    }

    template <typename TSecondaryView>
    void doInteraction(TSecondaryView& view, Code const projectileId, Code const targetId,
                       FourMomentum const& projectileP4, FourMomentum const& targetP4) {
      inner_.doInteraction(view, projectileId, targetId, projectileP4, targetP4);
      used_ = true;
    }

  private:
    NeutrinoInteraction inner_;
    bool used_ = false;
  };
} // namespace corsika::pythia8

#include <corsika/detail/modules/pythia8/NeutrinoInteraction.inl>