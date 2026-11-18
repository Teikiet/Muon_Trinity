/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/random/ExponentialDistribution.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/random/UniformRealDistribution.hpp>
#include <iostream>
#include <limits>
#include <random>

using namespace corsika;
using Catch::Approx;

SCENARIO("random-number streams can be registered and retrieved") {
  GIVEN("a RNGManager") {
    RNGManager<>& rngManager = RNGManager<>::getInstance();

    WHEN("the sequence name is not registered") {
      CHECK(rngManager.isRegistered("stream_A") == false);

      THEN("a sequence is registered by name") {
        rngManager.registerRandomStream("stream_A");

        THEN("the sequence can be retrieved") {
          CHECK_NOTHROW(rngManager.getRandomStream("stream_A"));

          THEN("we can check that the sequence exists") {
            CHECK_NOTHROW(rngManager.getRandomStream("stream_A"));

            THEN("an unknown sequence cannot be retrieved") {
              CHECK(rngManager.isRegistered("stream_A") == true);

              THEN("an unknown sequence cannot be retrieved") {
                CHECK_THROWS(rngManager.getRandomStream("stream_UNKNOWN"));

                THEN("an unknown sequence is not registered") {
                  CHECK(rngManager.isRegistered("stream_UNKNOWN") == false);
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("UniformRealDistribution") {

  logging::set_level(logging::level::info);

  std::mt19937 rng;

  corsika::UniformRealDistribution<LengthType> dist(1_m, 2_m);

  SECTION("range") {
    corsika::UniformRealDistribution<LengthType> dist(1_m, 2_m);

    LengthType min =
        +1_m * std::numeric_limits<typename LengthType::value_type>::infinity();
    LengthType max =
        -1_m * std::numeric_limits<typename LengthType::value_type>::infinity();

    for (int i{0}; i < 1'000'000; ++i) {
      LengthType x = dist(rng);
      min = std::min(min, x);
      max = std::max(max, x);
    }

    CHECK(min / 1_m == Approx(1.));
    CHECK(max / 2_m == Approx(1.));
  }

  SECTION("range") {
    corsika::UniformRealDistribution<LengthType> dist(18_cm);

    LengthType min =
        +1_m * std::numeric_limits<typename LengthType::value_type>::infinity();
    LengthType max =
        -1_m * std::numeric_limits<typename LengthType::value_type>::infinity();

    for (int i{0}; i < 1'000'000; ++i) {
      LengthType x = dist(rng);
      min = std::min(min, x);
      max = std::max(max, x);
    }

    CHECK(min / 1_m == Approx(0.).margin(1e-3));
    CHECK(max / 18_cm == Approx(1.));
  }
}

TEST_CASE("ExponentialDistribution") {

  logging::set_level(logging::level::info);

  std::mt19937 rng;

  auto const beta = 15_m;

  corsika::ExponentialDistribution dist(beta);

  SECTION("mean") {
    std::remove_const<decltype(beta)>::type mean = beta * 0;

    int constexpr N = 1'000'000;

    for (int i{0}; i < N; ++i) {
      decltype(beta) x = dist(rng);
      mean += x / N;
    }

    CHECK(mean / beta == Approx(1).margin(1e-2));
  }
}
