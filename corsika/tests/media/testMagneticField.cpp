/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/utility/CorsikaData.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/magnetic/GeomagneticModel.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

const std::string refDataDir = std::string(REFDATADIR); // from cmake

TEST_CASE("UniformMagneticField w/ Homogeneous Medium") {

  logging::set_level(logging::level::trace);

  auto const& gCS = get_root_CoordinateSystem();

  Point const gOrigin(0_m, 0_m, 0_m);

  SECTION("UniformMagneticField interface") {

    // setup our interface types
    using IModelInterface = media::IMagneticFieldModel<media::IMediumModel>;
    using AtmModel = UniformMagneticField<media::HomogeneousMedium<IModelInterface>>;

    // the composition we use for the homogenous medium
    media::NuclearComposition const protonComposition({Code::Proton}, {1.});

    // create a magnetic field vector
    Vector B0(0_T, 0_T, 0_T);

    // the constant density
    const auto density{19.2_g / cube(1_cm)};

    // create our atmospheric model
    AtmModel medium(B0, density, protonComposition);

    // and test at several locations
    CHECK(B0.getComponents() ==
          medium.getMagneticField(Point(-10_m, 4_m, 35_km)).getComponents());
    CHECK(B0.getComponents() ==
          medium.getMagneticField(Point(1000_km, -1000_km, 1000_km)).getComponents());
    CHECK(B0.getComponents() ==
          medium.getMagneticField(Point(0_m, 0_m, 0_m)).getComponents());

    // create a new magnetic field vector
    Vector B1(23_T, 57_T, -4_T);

    // and update this atmospheric model
    medium.setMagneticField(B1);

    // and test at several locations
    CHECK(B1.getComponents() ==
          medium.getMagneticField(Point(-10_m, 4_m, 35_km)).getComponents());
    CHECK(B1.getComponents() ==
          medium.getMagneticField(Point(1000_km, -1000_km, 1000_km)).getComponents());
    CHECK(B1.getComponents() ==
          medium.getMagneticField(Point(0_m, 0_m, 0_m)).getComponents());
  }

  SECTION("GeomagneticModel") {

    // non existing data file
    CHECK_THROWS(GeomagneticModel(gCS, gOrigin, corsika_data("GeoMag/IDoNotExist.COF")));
    // and empty data file
    CHECK_THROWS(
        GeomagneticModel(gCS, gOrigin, refDataDir + "/EmptyGeomagneticData.COF"));

    {
      media::GeomagneticModel wmm(gCS, gOrigin, corsika_data("GeoMag/WMM.COF"));

      // create earth magnetic field vector
      MagneticFieldVector WMM_B_1 = wmm.getField(2022.5, 100_km, -80, -120);
      MagneticFieldVector WMM_B_2 = wmm.getField(2022.5, 0_km, 0, 120);
      MagneticFieldVector WMM_B_3 = wmm.getField(2020, 100_km, 80, 0);
      MagneticFieldVector WMM_B_4 = wmm.getField(2019.999, 100_km, 80, 0);
      CHECK(WMM_B_1.getX() / 1_nT == Approx(5815).margin(0.05));
      CHECK(WMM_B_1.getY() / 1_nT == Approx(-14803).margin(0.05));
      CHECK(WMM_B_1.getZ() / 1_nT == Approx(49755.3).margin(0.05));
      CHECK(WMM_B_2.getX() / 1_nT == Approx(39684.7).margin(0.05));
      CHECK(WMM_B_2.getY() / 1_nT == Approx(42.2).margin(0.05));
      CHECK(WMM_B_2.getZ() / 1_nT == Approx(10809.5).margin(0.05));
      CHECK(WMM_B_3.getX() / 1_nT == Approx(6261.8).margin(0.05));
      CHECK(WMM_B_3.getY() / 1_nT == Approx(185.5).margin(0.05));
      CHECK(WMM_B_3.getZ() / 1_nT == Approx(-52429.1).margin(0.05));
      CHECK(WMM_B_4.getX() / 1_nT == Approx(6262).margin(0.5));
      CHECK(WMM_B_4.getY() / 1_nT == Approx(186).margin(0.5));
      CHECK(WMM_B_4.getZ() / 1_nT == Approx(-52429).margin(0.5));
    }

    {
      media::GeomagneticModel igrf(gCS, gOrigin, corsika_data("GeoMag/IGRF13.COF"));

      // create earth magnetic field vector
      MagneticFieldVector IGRF_B_1 = igrf.getField(2022.5, 100_km, -80, -120);
      MagneticFieldVector IGRF_B_2 = igrf.getField(1937.5, 10_km, 80, 120);
      MagneticFieldVector IGRF_B_3 = igrf.getField(1990, 50_km, 80, 0);
      MagneticFieldVector IGRF_B_4 = igrf.getField(2012, 0_km, 0, -120);
      CHECK(IGRF_B_1.getX() / 1_nT == Approx(5821).margin(0.5));
      CHECK(IGRF_B_1.getY() / 1_nT == Approx(-14796).margin(0.5));
      CHECK(IGRF_B_1.getZ() / 1_nT == Approx(49776).margin(0.5));
      CHECK(IGRF_B_2.getX() / 1_nT == Approx(2957).margin(0.5));
      CHECK(IGRF_B_2.getY() / 1_nT == Approx(-1014).margin(0.5));
      CHECK(IGRF_B_2.getZ() / 1_nT == Approx(-57610).margin(0.5));
      CHECK(IGRF_B_3.getX() / 1_nT == Approx(6545).margin(0.5));
      CHECK(IGRF_B_3.getY() / 1_nT == Approx(1473).margin(0.5));
      CHECK(IGRF_B_3.getZ() / 1_nT == Approx(-52884).margin(0.5));
      CHECK(IGRF_B_4.getX() / 1_nT == Approx(30052).margin(0.5));
      CHECK(IGRF_B_4.getY() / 1_nT == Approx(-4738).margin(0.5));
      CHECK(IGRF_B_4.getZ() / 1_nT == Approx(-5190).margin(0.5));
    }
  }
}
