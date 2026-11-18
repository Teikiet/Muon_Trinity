/*
 * (c) Copyright 2026 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/utility/FindXmax.hpp>

using namespace corsika;
using Catch::Approx;

TEST_CASE("find_xmax") {

  SECTION("interpolate_profile success") {
    auto const result = find_xmax::interpolate_profile(1., 2., 3., 0., 1., 0.);
    REQUIRE(result.has_value());
    auto const [Xmax, Ymax] = *result;
    CHECK(Xmax == Approx(2.).margin(1e-9));
    CHECK(Ymax == Approx(1.).margin(1e-9));
  }

  SECTION("interpolate_profile solve failure") {
    auto const result = find_xmax::interpolate_profile(1., 1., 1., 0., 1., 0.);
    CHECK_FALSE(result.has_value());
  }

  SECTION("interpolate_profile no maximum") {
    // upward-opening parabola: minimum at x=2
    auto const result = find_xmax::interpolate_profile(1., 2., 3., 1., 0., 1.);
    CHECK_FALSE(result.has_value());
  }

  SECTION("estimate_parabola_parameters solve failure") {
    std::vector<double> const x = {5., 5., 5.};
    std::vector<double> const y = {1., 2., 3.};
    auto const result = find_xmax::estimate_parabola_parameters(x, y);
    CHECK_FALSE(result.has_value());
  }

  SECTION("estimate_parabola_parameters success") {
    std::vector<double> const x = {1., 2., 3., 4., 5.};
    std::vector<double> const y = {0., 1., 0., -3., -8.};
    auto const result = find_xmax::estimate_parabola_parameters(x, y);
    REQUIRE(result.has_value());
    auto const [Xmax, Ymax] = *result;
    CHECK(Xmax == Approx(2.).margin(1e-6));
    CHECK(Ymax == Approx(1.).margin(1e-6));
  }

  SECTION("estimate_parabola_parameters large x offset") {
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i < 21; ++i) {
      double const xi = 990. + i;
      x.push_back(xi);
      y.push_back(1. - (xi - 1000.) * (xi - 1000.));
    }
    auto const result = find_xmax::estimate_parabola_parameters(x, y);
    REQUIRE(result.has_value());
    auto const [Xmax, Ymax] = *result;
    CHECK(Xmax == Approx(1000.).margin(1e-3));
    CHECK(Ymax == Approx(1.).margin(1e-3));
  }
}
