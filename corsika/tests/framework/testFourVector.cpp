/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <cmath>

using namespace corsika;
using Catch::Approx;

TEST_CASE("four vectors") {

  logging::set_level(logging::level::info);

  /*
    Test: P2 = E2 - p2 all in [GeV]
    This is the typical HEP application
   */
  SECTION("Energy momentum in hep-units") {

    HEPEnergyType E0 = 10_GeV;
    Vector<hepmomentum_d> P0(10_GeV, 10_GeV, 10_GeV);

    FourVector p0(E0, P0);

    CHECK(p0.getNormSqr() == -200_GeV * 1_GeV);
    CHECK(p0.getNorm() == sqrt(200_GeV * 1_GeV));
  }

  /*
    Check space/time-like
   */
  SECTION("Space/time likeness") {

    HEPEnergyType E0 = 20_GeV;
    Vector<hepmomentum_d> P0(10_GeV, 0_GeV, 0_GeV);
    Vector<hepmomentum_d> P1(10_GeV, 10_GeV, 20_GeV);
    Vector<hepmomentum_d> P2(0_GeV, 20_GeV, 0_GeV);

    FourVector p0(E0, P0);
    FourVector p1(E0, P1);
    FourVector p2(E0, P2);

    CHECK(p0.isSpacelike());
    CHECK(!p0.isTimelike());

    CHECK(!p1.isSpacelike());
    CHECK(p1.isTimelike());

    CHECK(!p2.isSpacelike());
    CHECK(!p2.isTimelike());
  }

  /*
    Test: P2 = E2/c2 - p2 with E in [GeV/c] and P in [GeV]
    This requires additional factors of c
   */
  SECTION("Energy momentum in SI-units") {

    auto E1 = 100_GeV / constants::c;
    Vector<hepmomentum_d> P1(10_GeV, 5_GeV, 15_GeV);

    FourVector p1(E1, P1);

    double const check = 100 * 100 - 10 * 10 - 5 * 5 - 15 * 15; // for dummies...

    CHECK(p1.getNormSqr() / 1_GeV / 1_GeV == Approx(check));
    CHECK(p1.getNorm() / 1_GeV == Approx(sqrt(check)));
  }

  /**
    Test: P2 = T2/c2 - r2 with T in [s] and r in [m]
    This requires additional factors of c
   */
  SECTION("Spacetime in SI-units") {

    TimeType T2 = 10_m / constants::c;
    Vector<length_d> P2(10_m, 5_m, 5_m);

    double const check = 10 * 10 - 10 * 10 - 5 * 5 - 5 * 5; // for dummies...

    FourVector p2(T2, P2);

    CHECK(p2.getNormSqr() == check * 1_m * 1_m);
    CHECK(p2.getNorm() == sqrt(abs(check)) * 1_m);
  }

  /**
     Testing the math operators
   */

  SECTION("Operators and comutions") {

    HEPEnergyType E1 = 100_GeV;
    Vector<hepmomentum_d> P1(0_GeV, 0_GeV, 0_GeV);

    HEPEnergyType E2 = 0_GeV;
    Vector<hepmomentum_d> P2(10_GeV, 0_GeV, 0_GeV);

    FourVector const p1(E1, P1);
    FourVector const p2(E2, P2);

    CHECK(p1.getNorm() / 1_GeV == Approx(100.));
    CHECK(p2.getNorm() / 1_GeV == Approx(10.));

    SECTION("product") {
      FourVector p3 = p1 + p2;
      CHECK(p3.getNorm() / 1_GeV == Approx(sqrt(100. * 100. - 100.)));
      p3 -= p2;
      CHECK(p3.getNorm() / 1_GeV == Approx(100.));
      CHECK(p1.getNorm() / 1_GeV == Approx(100.));
      CHECK(p2.getNorm() / 1_GeV == Approx(10.));
    }

    SECTION("difference") {
      FourVector p3 = p1 - p2;
      CHECK(p3.getNorm() / 1_GeV == Approx(sqrt(100. * 100. - 100.)));
      p3 += p2;
      CHECK(p3.getNorm() / 1_GeV == Approx(100.));
      CHECK(p1.getNorm() / 1_GeV == Approx(100.));
      CHECK(p2.getNorm() / 1_GeV == Approx(10.));
    }

    SECTION("scale") {
      double s = 10;
      FourVector p3 = p1 * s;
      CHECK(p3.getNorm() / 1_GeV == Approx(sqrt(100. * 100. * s * s)));
      p3 /= 10;
      CHECK(p3.getNorm() / 1_GeV == Approx(sqrt(100. * 100.)));
      CHECK(p1.getNorm() / 1_GeV == Approx(100.));
      CHECK(p2.getNorm() / 1_GeV == Approx(10.));
    }
  }

  /**
     The FourVector class can be used with reference template
     arguments. In this configuration it does not hold any data
     itself, but rather just refers to data located elsewhere. Thus,
     it merely provides the physical/mathematical wrapper around the
     data.
   */

  SECTION("Use as wrapper") {

    TimeType T = 10_m / constants::c;
    Vector<length_d> P(10_m, 5_m, 5_m);

    TimeType const T_c = 10_m / constants::c;
    Vector<length_d> const P_c(10_m, 5_m, 5_m);

    /*
      this does not compile, and it shoudn't!
      FourVector<TimeType&, Vector<length_d>&> p0(T_c, P_c);
    */
    FourVector<TimeType&, Vector<length_d>&> p1(T, P);
    FourVector<TimeType const&, Vector<length_d> const&> p2(T, P);
    FourVector<TimeType const&, Vector<length_d> const&> p3(T_c, P_c);

    p1 *= 10;
    // p2 *= 10; // this does not compile, and it shoudn't !
    // p3 *= 10; // this does not compile, and it shoudn't !!

    double const check = 10 * 10 - 10 * 10 - 5 * 5 - 5 * 5; // for dummies...
    CHECK(p1.getNormSqr() / (1_m * 1_m) == Approx(10. * 10. * check));
    CHECK(p2.getNorm() / 1_m == Approx(10 * sqrt(abs(check))));
    CHECK(p3.getNorm() / 1_m == Approx(sqrt(abs(check))));
  }
}
