/*
 * (c) Copyright 2026 CORSIKA Project, corsika8@kit.edu
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

template <typename TInterface>
using MyExtraEnv =
    ExponentialRefractiveIndex<MediumPropertyModel<UniformMagneticField<TInterface>>>;

TEST_CASE("ExponentialRefractiveIndex w/ Homogeneous medium") {

  logging::set_level(logging::level::info);

  // get a CS and a point
  auto const& gCS = get_root_CoordinateSystem();

  Point const gOrigin(gCS, {0_m, 0_m, 0_m});

  // setup interface types
  using IModelInterface = IRefractiveIndexModel<IMediumModel>;
  using AtmModel = ExponentialRefractiveIndex<media::HomogeneousMedium<IModelInterface>>;

  // the constant density
  const auto density{19.2_g / cube(1_cm)};

  // the composition we use for the homogenous medium
  media::NuclearComposition const protonComposition({Code::Proton}, {1.});

  // a new refractive index
  const double n0{1};
  const InverseLengthType lambda{6 / 1_m};

  // the center of the earth
  Point const center_{gCS, 0_m, 0_m, 0_m};
  // earth's radius
  LengthType const radius_{constants::EarthRadius::Mean};

  // create the atmospheric model and check refractive index
  AtmModel medium(n0, lambda, center_, constants::EarthRadius::Mean, density,
                  protonComposition);
  CHECK(n0 - medium.getRefractiveIndex(
                 Point(gCS, 0_m, 0_m, constants::EarthRadius::Mean)) ==
        Approx(0));

  // another refractive index
  const double n0_{1};
  const InverseLengthType lambda_{1 / 1_km};

  // distance from the center
  LengthType const dist_{4_km};

  // create the atmospheric model and check refractive index
  AtmModel medium_(n0_, lambda_, center_, dist_, density, protonComposition);
  CHECK(medium_.getRefractiveIndex(Point(gCS, 4_km, 3_km, 0_km)) == Approx(0.3678794412));

  // define axis vector
  Vector const axis(gCS, QuantityVector<dimensionless_d>(0, 0, 1));

  // check the density and nuclear composition
  REQUIRE(density == medium.getMassDensity(Point(gCS, 0_m, 0_m, 0_m)));
  medium.getNuclearComposition();
  REQUIRE(density == medium_.getMassDensity(Point(gCS, 0_m, 0_m, 0_m)));
  medium_.getNuclearComposition();

  SpeedType const velocity = 1_m / second;

  // the end time of our line
  TimeType const tEnd = 1_s;

  LengthType const length = tEnd * velocity;

  // create a line of length 1 m
  Line const line(gOrigin, Vector<SpeedType::dimension_type>(
                               gCS, {velocity, 0_m / second, 0_m / second}));

  // and the associated trajectory
  setup::Trajectory const track =
      setup::testing::make_track<setup::Trajectory>(line, tEnd);

  // and check the integrated grammage
  REQUIRE((medium.getIntegratedGrammage(track) / (density * length)) == Approx(1));
  REQUIRE((medium.getArclengthFromGrammage(track, density * 5_m) / 5_m) == Approx(1));
  REQUIRE((medium_.getIntegratedGrammage(track) / (density * length)) == Approx(1));
  REQUIRE((medium_.getArclengthFromGrammage(track, density * 5_m) / 5_m) == Approx(1));
}

TEST_CASE("ExponentialRefractiveIndex w/ 5-layered atmosphere") {

  logging::set_level(logging::level::info);

  // get a CS
  auto const& gCS = get_root_CoordinateSystem();

  // the center of the earth
  Point const center_{gCS, 0_m, 0_m, 0_m};

  // another refractive index
  const double n0{2};
  const InverseLengthType lambda{1 / 1_km};

  // a reference point to calculate the refractive index there
  Point const ref_{gCS, 0_m, 0_m, constants::EarthRadius::Mean};

  // setup a 5-layered environment
  using EnvironmentInterface =
      IRefractiveIndexModel<IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>>;
  using EnvType = Environment<EnvironmentInterface>;
  EnvType env;

  create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, AtmosphereId::LinsleyUSStd, center_, n0, lambda, center_,
      constants::EarthRadius::Mean, Medium::AirDry1Atm,
      MagneticFieldVector{gCS, 0_T, 50_uT, 0_T});

  // get the universe for this environment
  auto const* const universe{env.getUniverse().get()};
  auto const* node{universe->getContainingNode(ref_)};
  // get the refractive index
  auto const rIndex{node->getModelProperties().getRefractiveIndex(ref_)};

  CHECK(rIndex - n0 == Approx(0));
}