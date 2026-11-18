/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/CorsikaData.hpp>
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
#include <corsika/media/refractivity/LinearTabulatedRefractiveIndex.hpp>

#include <corsika/media/CORSIKA7Atmospheres.hpp>

#include <SetupTestTrajectory.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <catch2/catch_all.hpp>

#include <corsika/framework/utility/CorsikaFenv.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

TEST_CASE("Linear Tabulated Refractive Index Profile") {

  logging::set_level(logging::level::info);

  // get a CS and a point
  auto const& gCS = get_root_CoordinateSystem();

  Point const gOrigin(gCS, {0_m, 0_m, 0_m});

  // setup interface types
  using IModelInterface = IRefractiveIndexModel<IMediumModel>;
  using AtmModel =
      LinearTabulatedRefractiveIndex<media::HomogeneousMedium<IModelInterface>>;

  // the constant density
  const auto density{19.2_g / cube(1_cm)};

  // the composition we use for the homogenous medium
  media::NuclearComposition const protonComposition({Code::Proton}, {1.});

  // the center of the earth
  Point const center_{gCS, 0_m, 0_m, 0_m};
  // earth's radius

  // create the atmospheric model and check refractive index
  const auto transformer = PlanarTransformer(make_translation(
      gCS, {0_m, 0_m, constants::EarthRadius::Mean})); //, center_, radius_);
  AtmModel medium(corsika_data("CHERENKOV/atmosphere/atmprof1.dat").string(), transformer,
                  density, protonComposition);

  CHECK(medium.getRefractiveIndex(Point(gCS, 0_m, 0_m, constants::EarthRadius::Mean)) ==
        Approx(1.0 + 0.27072E-03));

  // another refractive index test with the same file

  // create the atmospheric model with a different transformer and check refractive index
  auto transformer2 =
      PlanarTransformer(make_translation(gCS, {0_m, 0_m, constants::EarthRadius::Mean}));
  AtmModel medium_(corsika_data("CHERENKOV/atmosphere/atmprof1.dat").string(),
                   transformer2, density, protonComposition);

  // Note: the refractive index depends on altitude from the file, not on the custom
  // parameters
  CHECK(medium_.getRefractiveIndex(
            Point(gCS, 0_m, 0_m, constants::EarthRadius::Mean + 120_km)) == Approx(1));

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
