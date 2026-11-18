/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/media/density_and_composition/FlatExponential.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/density_and_composition/InhomogeneousMedium.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/medium/MediumProperties.hpp>

#include <catch2/catch_all.hpp>

#include <SetupTestTrajectory.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

TEST_CASE("MediumProperties") {

  logging::set_level(logging::level::info);

  // test access of medium properties via enum and class types

  const Medium type = Medium::AirDry1Atm;
  const MediumData& air = mediumData(type);
  CHECK(air.getIeff() == 85.7);
  CHECK(air.getCbar() == 10.5961);
  CHECK(air.getX0() == 1.7418);
  CHECK(air.getX1() == 4.2759);
  CHECK(air.getAA() == 0.10914);
  CHECK(air.getSK() == 3.3994);
  CHECK(air.getDlt0() == 0.0);
}

TEST_CASE("MediumPropertyModel w/ Homogeneous") {

  logging::set_level(logging::level::info);

  Point const gOrigin(0_m, 0_m, 0_m);

  // setup our interface types
  using IModelInterface = IMediumPropertyModel<IMediumModel>;
  using AtmModel = media::MediumPropertyModel<media::HomogeneousMedium<IModelInterface>>;

  // the constant density
  const auto density{19.2_g / cube(1_cm)};

  // the composition we use for the homogenous medium
  media::NuclearComposition const protonComposition({Code::Proton}, {1.});

  // the refrative index that we use
  const Medium type = Medium::AirDry1Atm;

  // create the atmospheric model
  AtmModel medium(type, density, protonComposition);

  // and require that it is constant
  CHECK(type == medium.getMedium());

  // a new refractive index
  const Medium type2 = Medium::StandardRock;

  // update the refractive index of this atmospheric model
  medium.setMedium(type2);

  // check that the returned refractive index is correct
  CHECK(type2 == medium.getMedium());

  // define our axis vector
  Vector const axis(QuantityVector<dimensionless_d>(0, 0, 1));

  // check the density and nuclear composition
  CHECK(density == medium.getMassDensity(Point(0_m, 0_m, 0_m)));
  medium.getNuclearComposition();

  SpeedType const speed = 1_m / second;

  // create a line of length 1 m
  Line const line(gOrigin, VelocityVector(speed, 0_m / second, 0_m / second));

  // the end time of our line
  auto const tEnd = 1_s;

  LengthType const length = tEnd * speed;

  // and the associated trajectory
  setup::Trajectory const trajectory =
      corsika::setup::testing::make_track<setup::Trajectory>(line, tEnd);

  // and check the integrated grammage
  CHECK((medium.getIntegratedGrammage(trajectory) / (density * length)) == Approx(1));
  CHECK((medium.getArclengthFromGrammage(trajectory, density * 5_m) / 5_m) == Approx(1));
}
