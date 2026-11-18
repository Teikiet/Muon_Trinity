/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/framework/utility/COMBoost.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>

#include <limits>
#include <memory>
#include <random>
#include <tuple>
#include <PROPOSAL/particle/Particle.h>
#include <corsika/modules/proposal/InteractionModel.hpp>

namespace corsika::proposal {

  template <typename THadronicLEModel, typename THadronicHEModel>
  template <typename TEnvironment>
  inline InteractionModel<THadronicLEModel, THadronicHEModel>::InteractionModel(
      TEnvironment const& _env, THadronicLEModel& _hadintLE, THadronicHEModel& _hadintHE,
      HEPEnergyType const& _enthreshold)
      : ProposalProcessBase(_env)
      , HadronicPhotonModel<THadronicLEModel, THadronicHEModel>(_hadintLE, _hadintHE,
                                                                _enthreshold) {
    //! Initialize PROPOSAL tables for all media and all particles
    for (auto const& medium : media) {
      for (auto const particle_code : tracked) {
        buildCalculator(particle_code, medium.first);
      }
    }
  }

  template <typename THadronicLEModel, typename THadronicHEModel>
  inline void InteractionModel<THadronicLEModel, THadronicHEModel>::buildCalculator(
      Code code, size_t const& component_hash) {
    // search crosssection builder for given particle
    auto p_cross = cross.find(code);
    if (p_cross == cross.end())
      throw std::runtime_error("PROPOSAL could not find corresponding builder");

    // interpolate the crosssection for given media and energy cut. These may
    // take some minutes if you have to build the tables and cannot read the tables
    // from disk
    auto c = p_cross->second(media.at(component_hash), proposal_energycutsettings[code]);

    // Look which interactions take place and build the corresponding
    // interaction and secondary builder. The interaction integral will
    // interpolated too and saved in the calc map by a key build out of a hash
    // of composed of the component and particle code.
    auto inter_types = PROPOSAL::CrossSectionVector::GetInteractionTypes(c);
    calc_[std::make_pair(component_hash, code)] = std::make_tuple(
        PROPOSAL::make_secondaries(inter_types, particle[code], media.at(component_hash)),
        PROPOSAL::make_interaction(c, true, true),
        std::make_unique<LPM_calculator>(media.at(component_hash), code, inter_types));
  }

  template <typename THadronicLEModel, typename THadronicHEModel>
  template <typename TStackView>
  inline ProcessReturn
  InteractionModel<THadronicLEModel, THadronicHEModel>::doInteraction(
      TStackView& view, Code const projectileId,
      [[maybe_unused]] FourMomentum const& projectileP4) {

    auto const& rootCS = get_root_CoordinateSystem();
    auto const projectile = view.getProjectile();

    if (canInteract(projectileId)) {

      // get or build corresponding calculators
      auto c = getCalculator(projectile, calc_);

      // get the rates of the interaction types for every component.
      std::uniform_real_distribution<double> distr(0., 1.);

      // sample a interaction-type, loss and component
      auto rates =
          std::get<eINTERACTION>(c->second)->Rates(projectile.getEnergy() / 1_MeV);
      auto [type, target_hash, v] = std::get<eINTERACTION>(c->second)->SampleLoss(
          projectile.getEnergy() / 1_MeV, rates, distr(RNG_));

      // TODO: This should become obsolete as soon #482 is fixed
      if (type == PROPOSAL::InteractionType::Undefined) {
        CORSIKA_LOG_WARN(
            "PROPOSAL: No particle interaction possible. "
            "Put initial particle back on stack.");
        view.addSecondary(std::make_tuple(projectileId,
                                          projectile.getEnergy() - get_mass(projectileId),
                                          projectile.getDirection()));
        return ProcessReturn::Ok;
      }

      // Read how much random numbers are required to calculate the secondaries.
      // Calculate the secondaries and deploy them on the corsika stack.
      auto rnd = std::vector<double>(
          std::get<eSECONDARIES>(c->second)->RequiredRandomNumbers(type));
      for (auto& it : rnd) it = distr(RNG_);
      Point const& place = projectile.getPosition();

      auto point = PROPOSAL::Cartesian3D(place.getX() / 1_cm, place.getY() / 1_cm,
                                         place.getZ() / 1_cm);
      auto projectile_dir = projectile.getDirection();
      auto d = projectile_dir.getComponents();
      auto direction = PROPOSAL::Cartesian3D(d.getX().magnitude(), d.getY().magnitude(),
                                             d.getZ().magnitude());

      auto loss = PROPOSAL::StochasticLoss(
          static_cast<int>(type), v * projectile.getEnergy() / 1_MeV, point, direction,
          projectile.getTime() / 1_s, 0., projectile.getEnergy() / 1_MeV);
      PROPOSAL::Component target;
      if (type != PROPOSAL::InteractionType::Ioniz)
        target = PROPOSAL::Component::GetComponentForHash(target_hash);

      auto sec =
          std::get<eSECONDARIES>(c->second)->CalculateSecondaries(loss, target, rnd);

      // Check for LPM effect suppression!

      if (std::get<eLPM_SUPPRESSION>(c->second)) {
        auto lpm_suppression = CheckForLPM(
            *std::get<eLPM_SUPPRESSION>(c->second), projectile.getEnergy(), type, sec,
            projectile.getNode()->getModelProperties().getMassDensity(place), target, v);
        if (lpm_suppression) {
          // Discard interaction - put initial particle back on stack
          CORSIKA_LOG_DEBUG("LPM suppression detected!");
          view.addSecondary(std::make_tuple(
              projectileId, projectile.getEnergy() - get_mass(projectileId),
              projectile.getDirection()));
          return ProcessReturn::Ok;
        }
      }

      for (auto& s : sec) {
        auto E = s.energy * 1_MeV;
        auto vecProposal = s.direction;
        auto dir =
            DirectionVector({vecProposal.GetX(), vecProposal.GetY(), vecProposal.GetZ()});

        if (static_cast<PROPOSAL::ParticleType>(s.type) ==
            PROPOSAL::ParticleType::Hadron) {
          FourMomentum const photonP4(E, E * dir);
          // for PROPOSAL media target.GetAtomicNum() returns the atomic number in units
          // of atomic mass, so Nitrogen is 14.0067. When using media in CORSIKA the same
          // field is filled with the number of nucleons (ie. 14 for Nitrogen) To be sure
          // we explicitly cast to int
          auto const A = int(target.GetAtomicNum());
          auto const Z = int(target.GetNucCharge());
          Code const targetId = get_nucleus_code(A, Z);
          CORSIKA_LOGGER_DEBUG(
              logger_,
              "photo-hadronic interaction of projectile={} with target={}! Energy={} GeV",
              projectileId, targetId, E / 1_GeV);
          this->doHadronicPhotonInteraction(view, rootCS, photonP4, targetId);
        } else {
          auto sec_code = convert_from_PDG(static_cast<PDGCode>(s.type));
          // use mass provided by PROPOSAL to ensure correct conversion to kinetic energy
          auto massProposal =
              PROPOSAL::ParticleDef::GetParticleDefForType(s.type).mass * 1_MeV;
          view.addSecondary(std::make_tuple(sec_code, E - massProposal, dir));
        }
      }
    }
    return ProcessReturn::Ok;
  }

  template <typename THadronicLEModel, typename THadronicHEModel>
  template <typename TParticle>
  inline CrossSectionType
  InteractionModel<THadronicLEModel, THadronicHEModel>::getCrossSection(
      TParticle const& projectile, Code const projectileId,
      FourMomentum const& projectileP4) {

    // ==============================================
    // this block better diappears. RU 26.10.2021
    //
    // determine the volume where the particle is (last) known to be
    auto const* currentLogicalNode = projectile.getNode();
    media::NuclearComposition const& composition =
        currentLogicalNode->getModelProperties().getNuclearComposition();
    auto const meanMass = composition.getAverageMassNumber() * constants::u;
    // ==============================================

    if (canInteract(projectileId)) {
      auto c = getCalculator(projectile, calc_);
      return meanMass / (std::get<eINTERACTION>(c->second)->MeanFreePath(
                             projectileP4.getTimeLikeComponent() / 1_MeV) *
                         1_g / (1_cm * 1_cm));
    }

    return CrossSectionType::zero();
  }

  template <typename THadronicLEModel, typename THadronicHEModel>
  bool InteractionModel<THadronicLEModel, THadronicHEModel>::CheckForLPM(
      const LPM_calculator& lpm_calculator, const HEPEnergyType projectile_energy,
      const PROPOSAL::InteractionType type,
      const std::vector<PROPOSAL::ParticleState>& sec, const MassDensityType mass_density,
      const PROPOSAL::Component& target, const double v) {

    double suppression_factor;
    double current_density = mass_density * ((1_cm * 1_cm * 1_cm) / 1_g);

    if (type == PROPOSAL::InteractionType::Photopair && lpm_calculator.photo_pair_lpm_) {
      if (sec.size() != 2)
        throw std::runtime_error(
            "Invalid number of secondaries for Photopair in CheckForLPM");
      auto x =
          sec[0].energy / (sec[0].energy + sec[1].energy); // particle energy asymmetry
      double density_correction = current_density / lpm_calculator.mass_density_baseline_;
      suppression_factor = lpm_calculator.photo_pair_lpm_->suppression_factor(
          projectile_energy / 1_MeV, x, target, density_correction);
    } else if (type == PROPOSAL::InteractionType::Brems && lpm_calculator.brems_lpm_) {
      double density_correction = current_density / lpm_calculator.mass_density_baseline_;
      suppression_factor = lpm_calculator.brems_lpm_->suppression_factor(
          projectile_energy / 1_MeV, v, target, density_correction);
    } else if (type == PROPOSAL::InteractionType::Epair && lpm_calculator.epair_lpm_) {
      if (sec.size() != 3)
        throw std::runtime_error(
            "Invalid number of secondaries for EPair in CheckForLPM");
      double rho = (sec[1].energy - sec[2].energy) /
                   (sec[1].energy + sec[2].energy); // todo: check
      // it would be better if these variables are calculated within PROPOSAL
      double beta = (v * v) / (2 * (1 - v));
      double xi = (lpm_calculator.particle_mass_ * lpm_calculator.particle_mass_) /
                  (PROPOSAL::ME * PROPOSAL::ME) * v * v / 4 * (1 - rho * rho) / (1 - v);
      double density_correction = current_density / lpm_calculator.mass_density_baseline_;
      suppression_factor = lpm_calculator.epair_lpm_->suppression_factor(
          projectile_energy / 1_MeV, v, rho * rho, beta, xi, density_correction);
    } else {
      return false;
    }

    // rejection sampling
    std::uniform_real_distribution<double> distr(0., 1.);
    auto rnd_num = distr(RNG_);

    CORSIKA_LOGGER_TRACE(logger_,
                         "Checking for LPM suppression! energy: {}, v: {}, type: {}, "
                         "target: {}, mass_density: {}, mass_density_baseline: {}, "
                         "suppression_factor: {}, rnd: {}, result: {}, pmass: {} ",
                         projectile_energy, v, type, target.GetName(),
                         mass_density / (1_g / (1_cm * 1_cm * 1_cm)),
                         lpm_calculator.mass_density_baseline_, suppression_factor,
                         rnd_num, (rnd_num > suppression_factor),
                         lpm_calculator.particle_mass_);

    if (rnd_num > suppression_factor)
      return true;
    else
      return false;
  }
} // namespace corsika::proposal
