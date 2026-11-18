/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/geometry/Helix.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>

using namespace corsika;
using Catch::Approx;

double constexpr absMargin = 1.0e-8;

TEST_CASE("Helix class") {

  logging::set_level(logging::level::info);

  Point r0(0_m, 0_m, 0_m);

  SECTION("Helix") {
    Vector<SpeedType::dimension_type> const vPar(0_m / second, 0_m / second,
                                                 4_m / second);

    Vector<SpeedType::dimension_type> const vPerp(3_m / second, 0_m / second,
                                                  0_m / second);

    auto const T = 1_s;
    auto const omegaC = 2 * M_PI / T;

    Helix const helix(r0, omegaC, vPar, vPerp);

    CHECK((helix.getPosition(1_s).getCoordinates() -
           QuantityVector<length_d>(0_m, 0_m, 4_m))
              .getNorm()
              .magnitude() == Approx(0).margin(absMargin));

    CHECK((helix.getPosition(0.25_s).getCoordinates() -
           QuantityVector<length_d>(-3_m / (2 * M_PI), -3_m / (2 * M_PI), 1_m))
              .getNorm()
              .magnitude() == Approx(0).margin(absMargin));

    CHECK((helix.getPosition(7_s) -
           helix.getPositionFromArclength(helix.getArcLength(0_s, 7_s)))
              .getNorm()
              .magnitude() == Approx(0).margin(absMargin));

    /*
    // we have to consider this, if we need it
    auto const t = 1234_s;
    Trajectory<Helix> const base(helix, t);
    CHECK(helix.getPosition(t).GetCoordinates() == base.GetPosition(1.).GetCoordinates());

    CHECK(base.ArcLength(0_s, 1_s) / 1_m == Approx(5));
    */
  }
}
