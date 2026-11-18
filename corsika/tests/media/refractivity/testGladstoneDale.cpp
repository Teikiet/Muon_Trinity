/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/refractivity/UniformRefractiveIndex.hpp>
#include <corsika/media/refractivity/ExponentialRefractiveIndex.hpp>
#include <corsika/media/refractivity/GladstoneDaleRefractiveIndex.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>

#include <SetupTestTrajectory.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

template <typename TInterface2>
using MyExtraEnv2 =
    GladstoneDaleRefractiveIndex<MediumPropertyModel<UniformMagneticField<TInterface2>>>;

TEST_CASE("GladstoneDaleRefractiveIndex w/ Homogeneous medium") {

  logging::set_level(logging::level::info);

  // get a CS and a point
  auto const& gCS = get_root_CoordinateSystem();

  Point const gOrigin(gCS, {0_m, 0_m, 0_m});

  // setup interface types
  using IModelInterface = IRefractiveIndexModel<IMediumModel>;
  using AtmModel =
      GladstoneDaleRefractiveIndex<media::HomogeneousMedium<IModelInterface>>;

  // the constant density
  const auto density{19.2_g / cube(1_cm)};

  // the composition we use for the homogenous medium
  media::NuclearComposition const protonComposition({Code::Proton}, {1.});

  // the refractive index at sea level
  const double n0{1.000327};

  // a point at the surface of the earth
  Point const surface_{gCS, 0_m, 0_m, constants::EarthRadius::Mean};

  // a random point in the atmosphere
  Point const p1_{gCS, 1_km, 1_km, constants::EarthRadius::Mean + 10_km};

  // create the atmospheric model and check refractive index
  AtmModel medium(n0, surface_, density, protonComposition);

  CHECK(n0 - medium.getRefractiveIndex(surface_) == Approx(0));
  CHECK(n0 - medium.getRefractiveIndex(p1_) == Approx(0));
}

TEST_CASE("GladstoneDaleRefractiveIndex w/ 5-layered atmosphere") {

  logging::set_level(logging::level::info);

  // get a CS
  auto const& gCS = get_root_CoordinateSystem();

  // the center of the earth
  Point const center_{gCS, 0_m, 0_m, 0_m};

  // a point at the surface of the earth
  Point const surface_{gCS, 0_m, 0_m, constants::EarthRadius::Mean};

  // the refractive index at sea level
  const double n0{1.000327};

  // a reference point to calculate the refractive index there
  Point const ref_{gCS, 0_km, 0_km, constants::EarthRadius::Mean + 10_km};

  // setup a 5-layered environment
  using EnvironmentInterface =
      IRefractiveIndexModel<IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>>;
  using EnvType = Environment<EnvironmentInterface>;
  EnvType env;

  create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv2>(
      env, AtmosphereId::LinsleyUSStd, center_, n0, surface_, Medium::AirDry1Atm,
      MagneticFieldVector{gCS, 0_T, 50_uT, 0_T});

  // get the universe for this environment
  auto const* const universe{env.getUniverse().get()};
  auto const* node{universe->getContainingNode(ref_)};
  // get the refractive index
  auto const rIndex1{node->getModelProperties().getRefractiveIndex(ref_)};
  auto const rIndex2{node->getModelProperties().getRefractiveIndex(surface_)};

  CHECK(rIndex1 - n0 == Approx(-0.0002591034));
  CHECK(rIndex2 - n0 == Approx(0));
}
