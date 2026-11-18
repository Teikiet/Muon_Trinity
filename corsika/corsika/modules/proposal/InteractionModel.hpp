/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <PROPOSAL/PROPOSAL.h>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/process/ProcessReturn.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/random/UniformRealDistribution.hpp>
#include <corsika/modules/proposal/ProposalProcessBase.hpp>
#include <corsika/modules/proposal/HadronicPhotonModel.hpp>

namespace corsika::proposal {

  //!
  //! Electro-magnetic and photon stochastic losses produced by proposal. It makes
  //! use of interpolation tables which are runtime intensive calculation, but can be
  //! reused by setting the \param PROPOSAL::InterpolationDef::path_to_tables variable.
  //! Hadroninc interactions of photons with nuclei are included. The cross section is
  //! calculated by PROPOSAL. For the production of hadronic secondaries an external model
  //! is needed that implements the
  //! doInteraction(TSecondaries& view, Code const projectile, Code const
  //! target,FourMomentum const& projectileP4, FourMomentum const& targetP4) routine.
  //! @tparam THadronicModel
  //!

  template <class THadronicLEModel, class THadronicHEModel>
  class InteractionModel
      : public ProposalProcessBase,
        public HadronicPhotonModel<THadronicLEModel, THadronicHEModel> {

    enum { eSECONDARIES, eINTERACTION, eLPM_SUPPRESSION };
    struct LPM_calculator;
    using calculator_t = std::tuple<std::unique_ptr<PROPOSAL::SecondariesCalculator>,
                                    std::unique_ptr<PROPOSAL::Interaction>,
                                    std::unique_ptr<LPM_calculator>>;

    std::unordered_map<calc_key_t, calculator_t, hash>
        calc_; //!< Stores the secondaries and interaction calculators.

    //!
    //! Build the secondaries and interaction calculators and add it to calc.
    //!
    void buildCalculator(Code, size_t const&) final;

    inline static auto logger_{get_logger("corsika_proposal_InteractionModel")};

    // Calculators for the LPM effect
    struct LPM_calculator {
      std::unique_ptr<PROPOSAL::crosssection::PhotoPairLPM> photo_pair_lpm_ = nullptr;
      std::unique_ptr<PROPOSAL::crosssection::BremsLPM> brems_lpm_ = nullptr;
      std::unique_ptr<PROPOSAL::crosssection::EpairLPM> epair_lpm_ = nullptr;

      const double mass_density_baseline_; // base mass density, note PROPOSAL units
      const double particle_mass_;         // particle mass, note PROPOSAL units

      LPM_calculator(const PROPOSAL::Medium& medium, const Code code,
                     std::vector<PROPOSAL::InteractionType> inter_types)
          : mass_density_baseline_(medium.GetMassDensity())
          , particle_mass_(particle[code].mass) {
        // at the moment, we always calculate the LPM effect based on the default cross
        // sections. one could check for the parametrizations used in the other
        // calculators.
        if (std::find(inter_types.begin(), inter_types.end(),
                      PROPOSAL::InteractionType::Photopair) != inter_types.end())
          photo_pair_lpm_ = std::make_unique<PROPOSAL::crosssection::PhotoPairLPM>(
              particle[code], medium, PROPOSAL::crosssection::PhotoPairKochMotz());
        if (std::find(inter_types.begin(), inter_types.end(),
                      PROPOSAL::InteractionType::Brems) != inter_types.end())
          brems_lpm_ = std::make_unique<PROPOSAL::crosssection::BremsLPM>(
              particle[code], medium, PROPOSAL::crosssection::BremsElectronScreening());
        if (std::find(inter_types.begin(), inter_types.end(),
                      PROPOSAL::InteractionType::Epair) != inter_types.end())
          epair_lpm_ =
              std::make_unique<PROPOSAL::crosssection::EpairLPM>(particle[code], medium);
      }
    };

    //!
    //! Checks whether the interaction is suppressed by the LPM effect
    //!
    bool CheckForLPM(const LPM_calculator&, const HEPEnergyType,
                     const PROPOSAL::InteractionType,
                     const std::vector<PROPOSAL::ParticleState>&, const MassDensityType,
                     const PROPOSAL::Component&, const double v);

  public:
    //!
    //! Produces the stoachastic loss calculator for leptons based on nuclear
    //! compositions and stochastic description limited by the particle cut.
    //!
    template <typename TEnvironment>
    InteractionModel(TEnvironment const& env, THadronicLEModel&, THadronicHEModel&,
                     HEPEnergyType const&);

    //!
    //! Calculate the rates for the different targets and interactions. Sample a
    //! pair of interaction-type, component and rate, followed by sampling a loss and
    //! produce the corresponding secondaries and store them on the particle stack.
    //! interactions in PROPOSAL are:
    //!
    //! InteractionType::Particle
    //! InteractionType::Brems
    //! InteractionType::Ioniz
    //! InteractionType::Epair
    //! InteractionType::Photonuclear
    //! InteractionType::MuPair
    //! InteractionType::Hadrons
    //! InteractionType::ContinuousEnergyLoss
    //! InteractionType::WeakInt
    //! InteractionType::Compton
    //! InteractionType::Decay
    //! InteractionType::Annihilation
    //! InteractionType::Photopair
    //! InteractionType::Photoproduction
    //! InteractionType::Photoeffect
    //!
    //! more information can be found at:
    //! https://github.com/tudo-astroparticlephysics/PROPOSAL
    template <typename TSecondaryView>
    ProcessReturn doInteraction(TSecondaryView&, Code const projectileId,
                                FourMomentum const& projectileP4);

    //!
    //! Calculates and returns the cross section.
    //!
    template <typename TParticle>
    CrossSectionType getCrossSection(TParticle const& p, Code const projectileId,
                                     FourMomentum const& projectileP4);
  };

} // namespace corsika::proposal

#include <corsika/detail/modules/proposal/InteractionModel.inl>
