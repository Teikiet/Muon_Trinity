/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/media/AtmosphereTable.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/FlatInterpolatedMedium.hpp>
#include <corsika/setup/SetupTrajectory.hpp>
#include <SetupTestTrajectory.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using Catch::Approx;

const std::string refDataDir =
    std::string(REFDATADIR) + "/atmosphere_table"; // from cmake

TEST_CASE("AtmosphereTable") {
  logging::set_level(logging::level::debug);

  // origin
  Point const gOrigin(0_m, 0_m, 0_m);

  // not-existing table
  CHECK_THROWS(media::AtmosphereTable(gOrigin, "fake_path.dat"));

  // test a table with duplicate entries
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_duplicates.dat"));

  // test a table with only one row
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_fewRows.dat"));

  // test a table with inconsistent rows
  CHECK_THROWS(
      media::AtmosphereTable(gOrigin, refDataDir + "/badTable_inconsistent.dat"));

  // test a table with missing mandatory parameters
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_missingPar.dat"));

  // test a table with a non-numerical value
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_NaN.dat"));

  // test a table with a negative density value
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_negDens.dat"));

  // test a table with a negative thickness value
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_negThick.dat"));

  // test a table with unsorted table rows
  CHECK_THROWS(media::AtmosphereTable(gOrigin, refDataDir + "/badTable_unsorted.dat"));

  // test a table with undefined values at zero altitude - tests that the first table row
  // was replicated at 0_m
  auto atm = media::AtmosphereTable(gOrigin, refDataDir + "/table_zeroAlt.dat");
  CHECK(atm.getDensityValues()[0] == 0.10651E-02_g / static_pow<3>(1_cm));
  CHECK(atm.getDensityValues()[1] == 0.10651E-02_g / static_pow<3>(1_cm));

  // query an undefined parameter
  atm = media::AtmosphereTable(gOrigin, refDataDir + "/table_undefPar.dat");
  CHECK_THROWS(atm.getRefractiveIndexValues());
}

TEST_CASE("FlatInterpolatedMedium") {
  logging::set_level(logging::level::info);

  // coordinate system and origin
  auto const& gCS = get_root_CoordinateSystem();
  Point const gOrigin(0_m, 0_m, 0_m);

  // environment
  media::Environment<corsika::media::IMediumModel> env;

  // create atmosphere from atmprof1.dat
  auto atm =
      media::AtmosphereTable(gOrigin, corsika_data("CHERENKOV/atmosphere/atmprof1.dat"));
  media::create_interpolated_atmosphere<corsika::media::IMediumModel>(env, atm);

  auto earthRadius = constants::EarthRadius::Mean;

  SECTION("density") {
    // tests density values between two table entries in atmprof1.dat
    auto const alt0 = 50_km + earthRadius;
    auto const alt1 = 55_km + earthRadius;
    auto const pt0 = Point(gCS, 0_m, 0_m, alt0);
    auto const pt1 = Point(gCS, 0_m, 0_m, alt1);
    auto const rho0 =
        env.getUniverse()->getContainingNode(pt0)->getModelProperties().getMassDensity(
            pt0);
    auto const rho1 =
        env.getUniverse()->getContainingNode(pt1)->getModelProperties().getMassDensity(
            pt1);

    // compare density with exact tabulated values
    CHECK(rho0 / 1_g * static_pow<3>(1_cm) / 1.1019E-06 == Approx(1));
    CHECK(rho1 / 1_g * static_pow<3>(1_cm) / 6.0359E-07 == Approx(1));

    // calculate density slope in logarithmic space
    auto const logAlt0 = log10(alt0 / 1_m);
    auto const logAlt1 = log10(alt1 / 1_m);
    auto const logRho0 = log10(rho0 / 1_g * static_pow<3>(1_cm));
    auto const logRho1 = log10(rho1 / 1_g * static_pow<3>(1_cm));
    auto const logRhoSlope = (logRho1 - logRho0) / (logAlt1 - logAlt0);

    // interpolate to midpoint between tabulated values
    auto alt = 52.5_km + earthRadius;
    auto pt = Point(gCS, 0_m, 0_m, alt);
    auto const logAlt = log10(alt / 1_m);
    auto const logRho = logRho0 + logRhoSlope * (logAlt - logAlt0);
    auto const rhoCalc = std::pow(10, logRho) * 1_g / static_pow<3>(1_cm);

    // get density from the medium
    auto const rho =
        env.getUniverse()->getContainingNode(pt)->getModelProperties().getMassDensity(pt);

    // compare density with interpolated value
    CHECK(rho / rhoCalc == Approx(1));
  }

  SECTION("grammage and arclength") {
    // test medium by creating tracks, obtaining the integrated grammage, reverting that
    // back to track length and comparing to the real length

    // test a horizontal track
    auto lengthExact = 500_m;
    auto alt = 30_km;
    auto pt = Point(gCS, 0_m, 0_m, alt + earthRadius);
    auto line = Line(pt, Vector<SpeedType::dimension_type>(
                             gCS, {lengthExact / second, 0_m / second, 0_m / second}));
    auto trajectory1 = setup::testing::make_track<setup::Trajectory>(line, 1_s);

    // integrated grammage obtained from the medium
    auto grammage = env.getUniverse()
                        ->getContainingNode(pt)
                        ->getModelProperties()
                        .getIntegratedGrammage(trajectory1);
    // invert grammage to track length
    auto length = env.getUniverse()
                      ->getContainingNode(pt)
                      ->getModelProperties()
                      .getArclengthFromGrammage(trajectory1, grammage);
    // compare to exact length
    CHECK(length / lengthExact == Approx(1));

    // test a long vertical track
    lengthExact = 20_km;
    alt = 10_km;
    pt = Point(gCS, 0_m, 0_m, alt + earthRadius);
    line = Line(pt, Vector<SpeedType::dimension_type>(
                        gCS, {0_m / second, 0_m / second, lengthExact / second}));
    auto trajectory2 = setup::testing::make_track<setup::Trajectory>(line, 1_s);

    // integrated grammage obtained from the medium
    grammage = env.getUniverse()
                   ->getContainingNode(pt)
                   ->getModelProperties()
                   .getIntegratedGrammage(trajectory2);
    // invert grammage to track length
    length = env.getUniverse()
                 ->getContainingNode(pt)
                 ->getModelProperties()
                 .getArclengthFromGrammage(trajectory2, grammage);
    // compare to exact length
    CHECK(length / lengthExact == Approx(1));

    // test a short inclined track
    auto incl = 60.0 * (M_PI / 180);
    auto zDiff = 10_m;
    auto xDiff = zDiff / cos(incl);
    lengthExact = sqrt(xDiff * xDiff + zDiff * zDiff);
    alt = 1_km;
    pt = Point(gCS, 0_m, 0_m, alt + earthRadius);
    line = Line(pt, Vector<SpeedType::dimension_type>(
                        gCS, {xDiff / second, 0_m / second, zDiff / second}));
    auto trajectory3 = setup::testing::make_track<setup::Trajectory>(line, 1_s);

    // integrated grammage obtained from the medium
    grammage = env.getUniverse()
                   ->getContainingNode(pt)
                   ->getModelProperties()
                   .getIntegratedGrammage(trajectory3);
    // invert grammage to track length
    length = env.getUniverse()
                 ->getContainingNode(pt)
                 ->getModelProperties()
                 .getArclengthFromGrammage(trajectory3, grammage);
    // compare to exact length
    CHECK(length / lengthExact == Approx(1));
  }
}
