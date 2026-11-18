/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/interfaces/IMediumPropertyModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>

#include <limits>

namespace corsika {

  using DummyEnvironmentInterface =
      media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
  using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

  namespace setup::testing {

    /**
     * \function setup_environment
     *
     * standard environment for unit testing.
     *
     */

    inline std::tuple<std::unique_ptr<DummyEnvironment>, CoordinateSystem const*,
                      DummyEnvironment::BaseNodeType*>
    setup_environment(Code const vTargetCode,
                      MagneticFluxType const& BfieldZ = MagneticFluxType::zero()) {

      auto env = std::make_unique<DummyEnvironment>();
      auto& universe = *(env->getUniverse());
      auto const& cs = env->getCoordinateSystem();

      /**
       * our world is a sphere at 0,0,0 with R=infty
       */
      auto world = DummyEnvironment::createNode<Sphere>(Point{cs, 0_m, 0_m, 0_m}, 100_km);

      /**
       * construct suited environment medium model:
       */
      using MyHomogeneousModel = media::MediumPropertyModel<media::UniformMagneticField<
          media::HomogeneousMedium<DummyEnvironmentInterface>>>;

      world->setModelProperties<MyHomogeneousModel>(
          media::Medium::AirDry1Atm, Vector(cs, 0_T, 0_T, BfieldZ),
          1_kg / (1_m * 1_m * 1_m),
          media::NuclearComposition(std::vector<Code>{vTargetCode},
                                    std::vector<double>{1.}));

      DummyEnvironment::BaseNodeType* nodePtr = world.get();
      universe.addChild(std::move(world));

      return std::make_tuple(std::move(env), &cs, nodePtr);
    }

  } // namespace setup::testing

} // namespace corsika
