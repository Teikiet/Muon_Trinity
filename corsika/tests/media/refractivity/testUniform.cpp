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

TEST_CASE("UniformRefractiveIndex w/ Homogeneous medium") {

  logging::set_level(logging::level::info);

  auto const& gCS = get_root_CoordinateSystem();

  Point const gOrigin(gCS, {0_m, 0_m, 0_m});

  // set up our interface types
  using IModelInterface = IRefractiveIndexModel<IMediumModel>;
  using AtmModel = UniformRefractiveIndex<media::HomogeneousMedium<IModelInterface>>;

  // the constant density
  const auto density{19.2_g / cube(1_cm)};

  // the composition we use for the homogenous medium
  media::NuclearComposition const protonComposition({Code::Proton}, {1.});

  // the refractive index that we use
  const double n{1.000327};

  // create the atmospheric model
  AtmModel medium(n, density, protonComposition);

  // and require that it is constant
  CHECK(n == medium.getRefractiveIndex(Point(gCS, -10_m, 4_m, 35_km)));
  CHECK(n == medium.getRefractiveIndex(Point(gCS, +210_m, 0_m, 7_km)));
  CHECK(n == medium.getRefractiveIndex(Point(gCS, 0_m, 0_m, 0_km)));
  CHECK(n == medium.getRefractiveIndex(Point(gCS, 100_km, 400_km, 350_km)));

  // a new refractive index
  const double n2{2.3472123};

  // update the refractive index of this atmospheric model
  medium.setRefractiveIndex(n2);

  // check that the returned refractive index is correct
  CHECK(n2 == medium.getRefractiveIndex(Point(gCS, -10_m, 4_m, 35_km)));
  CHECK(n2 == medium.getRefractiveIndex(Point(gCS, +210_m, 0_m, 7_km)));
  CHECK(n2 == medium.getRefractiveIndex(Point(gCS, 0_m, 0_m, 0_km)));
  CHECK(n2 == medium.getRefractiveIndex(Point(gCS, 100_km, 400_km, 350_km)));

  // define our axis vector
  Vector const axis(gCS, QuantityVector<dimensionless_d>(0, 0, 1));

  // check the density and nuclear composition
  CHECK(density == medium.getMassDensity(Point(gCS, 0_m, 0_m, 0_m)));
  medium.getNuclearComposition();

  SpeedType const speed = 1_m / second;

  // create a line of length 1 m
  Line const line(gOrigin, VelocityVector(gCS, {speed, 0_m / second, 0_m / second}));

  // the end time of our line
  auto const tEnd = 1_s;

  LengthType const length = tEnd * speed;

  // and the associated trajectory
  setup::Trajectory const track =
      setup::testing::make_track<setup::Trajectory>(line, tEnd);

  // and check the integrated grammage
  CHECK((medium.getIntegratedGrammage(track) / (density * length)) == Approx(1));
  CHECK((medium.getArclengthFromGrammage(track, density * 5_m) / 5_m) == Approx(1));
}
