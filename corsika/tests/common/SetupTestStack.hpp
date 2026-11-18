/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>

#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <tests/common/SetupStack.hpp>
#include <SetupStack.hpp>

/**
 * \file SetupTestStack
 *
 * standard stack setup for unit tests.
 */

namespace corsika {

  using DummyEnvironmentInterface =
      media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
  using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

  namespace setup::testing {

    /**
     * \function setup_stack
     *
     * standard stack setup for unit tests.
     *
     * \return a tuple with element 0 being a Stack object filled with
     * one particle, and element 1 the StackView on it.
     */

    inline std::tuple<std::unique_ptr<test::Stack>, std::unique_ptr<test::StackView>>
    setup_stack(Code const vProjectileType, HEPEnergyType const vMomentum,
                DummyEnvironment::BaseNodeType* const vNodePtr,
                CoordinateSystem const& cs) {

      auto stack = std::make_unique<test::Stack>();

      Point const origin(cs, {0_m, 0_m, 0_m});
      MomentumVector const pLab(cs, {vMomentum, 0_GeV, 0_GeV});
      HEPMassType const mass = get_mass(vProjectileType);
      HEPEnergyType const Ekin = sqrt(vMomentum * vMomentum + mass * mass) - mass;

      auto particle = stack->addParticle(
          std::make_tuple(vProjectileType, Ekin, pLab.normalized(), origin, 0_ns));
      particle.setNode(vNodePtr);
      return std::make_tuple(std::move(stack),
                             std::make_unique<test::StackView>(particle));
    }

  } // namespace setup::testing

} // namespace corsika