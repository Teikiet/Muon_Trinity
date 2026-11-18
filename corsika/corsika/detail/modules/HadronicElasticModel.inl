/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/HadronicElasticModel.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/random/ExponentialDistribution.hpp>
#include <corsika/framework/utility/COMBoost.hpp>

#include <iomanip>
#include <iostream>

namespace corsika {

  inline HadronicElasticInteraction::HadronicElasticInteraction(CrossSectionType x,
                                                                CrossSectionType y)
      : parX_(x)
      , parY_(y) {}

  template <typename TParticle>
  inline GrammageType HadronicElasticInteraction::getInteractionLength(
      TParticle const& p) {
    if (p.getPID() == Code::Proton) {
      auto const* currentNode = p.getNode();
      auto const& mediumComposition =
          currentNode->getModelProperties().getNuclearComposition();

      auto const& components = mediumComposition.getComponents();
      auto const& fractions = mediumComposition.getFractions();

      auto const projectileMomentum = p.getMomentum();
      auto const projectileMomentumSquaredNorm = projectileMomentum.getSquaredNorm();
      auto const projectileEnergy = p.getEnergy();

      auto const avgCrossSection = [&]() {
        CrossSectionType avgCrossSection = 0_b;

        for (size_t i = 0; i < fractions.size(); ++i) {
          auto const targetMass = get_mass(components[i]);
          auto const s = static_pow<2>(projectileEnergy + targetMass) -
                         projectileMomentumSquaredNorm;
          avgCrossSection += getCrossSection(s) * fractions[i];
        }

        CORSIKA_LOG_DEBUG("avgCrossSection: {} mb", avgCrossSection / 1_mb);

        return avgCrossSection;
      }();

      auto const avgTargetMassNumber = mediumComposition.getAverageMassNumber();

      GrammageType const interactionLength =
          avgTargetMassNumber * constants::u / avgCrossSection;

      return interactionLength;
    } else {
      return std::numeric_limits<double>::infinity() * 1_g / (1_cm * 1_cm);
    }
  }

  template <typename TParticle>
  inline ProcessReturn HadronicElasticInteraction::doInteraction(TParticle& p) {
    if (p.getPID() != Code::Proton) { return ProcessReturn::Ok; }

    const auto* currentNode = p.getNode();
    const auto& composition = currentNode->getModelProperties().getNuclearComposition();
    const auto& components = composition.getComponents();

    std::vector<CrossSectionType> cross_section_of_components(
        composition.getComponents().size());

    auto const projectileMomentum = p.getMomentum();
    auto const projectileMomentumSquaredNorm = projectileMomentum.getSquaredNorm();
    auto const projectileEnergy = p.getEnergy();

    for (size_t i = 0; i < components.size(); ++i) {
      auto const targetMass = get_mass(components[i]);
      auto const s =
          static_pow<2>(projectileEnergy + targetMass) - projectileMomentumSquaredNorm;
      cross_section_of_components[i] = CrossSection(s);
    }

    const auto targetCode = composition.SampleTarget(cross_section_of_components, RNG_);

    auto const targetMass = get_mass(targetCode);

    std::uniform_real_distribution phiDist(0., 2 * M_PI);

    FourVector const projectileLab(projectileEnergy, projectileMomentum);
    FourVector const targetLab(targetMass, MomentumVector({0_eV, 0_eV, 0_eV}));
    COMBoost const boost(projectileLab, targetMass);

    auto const projectileCoM = boost.toCoM(projectileLab);
    auto const targetCoM = boost.toCoM(targetLab);

    auto const pProjectileCoMSqNorm =
        projectileCoM.getSpaceLikeComponents().getSquaredNorm();
    auto const pProjectileCoMNorm = sqrt(pProjectileCoMSqNorm);

    auto const eProjectileCoM = projectileCoM.getTimeLikeComponent();
    auto const eTargetCoM = targetCoM.getTimeLikeComponent();

    auto const sqrtS = eProjectileCoM + eTargetCoM;
    auto const s = static_pow<2>(sqrtS);

    auto const B = this->B(s);
    CORSIKA_LOG_DEBUG(B);

    ExponentialDistribution tDist(1 / B);
    auto const absT = [&]() {
      decltype(tDist(RNG_)) absT;
      auto const maxT = 4 * pProjectileCoMSqNorm;

      do {
        // |t| cannot become arbitrarily large, max. given by GER eq. (4.16), so we just
        // throw again until we have an acceptable value. Note that the formula holds in
        // any frame despite of what is stated in the book.
        absT = tDist(RNG_);
      } while (absT >= maxT);

      return absT;
    }();

    CORSIKA_LOG_DEBUG(
        "HadronicElasticInteraction: s = {}"
        " GeV²; absT = {} "
        " GeV² (max./GeV² = {})",
        s * constants::invGeVsq, absT * constants::invGeVsq,
        4 * constants::invGeVsq * projectileMomentumSquaredNorm);

    auto const theta = 2 * asin(sqrt(absT / (4 * pProjectileCoMSqNorm)));
    auto const phi = phiDist(RNG_);

    auto const projectileScatteredLab =
        boost.fromCoM(FourVector<HEPEnergyType, MomentumVector>(
            eProjectileCoM, boost.makeCoMVector(QuantityVector<hepmomentum_d>{
                                pProjectileCoMNorm * sin(theta) * cos(phi),
                                pProjectileCoMNorm * sin(theta) * sin(phi),
                                pProjectileCoMNorm * cos(theta)})));

    p.setMomentum(projectileScatteredLab.getSpaceLikeComponents());
    p.setEnergy(
        sqrt(projectileScatteredLab.getSpaceLikeComponents().getSquaredNorm() +
             static_pow<2>(get_mass(
                 p.getPID())))); // Don't use energy from boost. It can be smaller than
                                 // the momentum due to limited numerical accuracy.

    return ProcessReturn::Ok;
  }

  inline HadronicElasticInteraction::inveV2 HadronicElasticInteraction::B(eV2 s) const {
    auto constexpr b_p = 2.3;
    auto const result =
        (2 * b_p + 2 * b_p + 4 * pow(s * constants::invGeVsq, gfEpsilon) - 4.2) *
        constants::invGeVsq;
    CORSIKA_LOG_DEBUG("B({}) = {}  GeV¯²", s, result / constants::invGeVsq);

    return result;
  }

  inline CrossSectionType HadronicElasticInteraction::getCrossSection(
      SquaredHEPEnergyType s) const {
    // assuming every target behaves like a proton, parX_ and parY_ are universal
    CrossSectionType const sigmaTotal = parX_ * pow(s * constants::invGeVsq, gfEpsilon) +
                                        parY_ * pow(s * constants::invGeVsq, -gfEta);

    // according to Schuler & Sjöstrand, PRD 49, 2257 (1994)
    // (we ignore rho because rho^2 is just ~2 %)
    auto const sigmaElastic =
        static_pow<2>(sigmaTotal) /
        (16 * constants::pi * convert_HEP_to_SI<CrossSectionType::dimension_type>(B(s)));

    CORSIKA_LOG_DEBUG("HEM sigmaTot = {} mb", sigmaTotal / 1_mb);
    CORSIKA_LOG_DEBUG("HEM sigmaElastic = {} mb", sigmaElastic / 1_mb);
    return sigmaElastic;
  }

} // namespace corsika
