/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/modules/thinning/EMThinning.hpp>

#include <SetupStack.hpp>
#include <SetupTestTrajectory.hpp>
#include <SetupTestEnvironment.hpp>

#include <catch2/catch_all.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

using namespace corsika;
using Catch::Approx;

using DummyEnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

TEST_CASE("EMThinning", "process,secondary") {
  RNGManager<>::getInstance().registerRandomStream("thinning");
  logging::set_level(logging::level::info);

  feenableexcept(FE_INVALID);
  using EnvType = DummyEnvironment;

  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();

  // setup empty particle stack
  test::Stack stack;
  stack.clear();

  SECTION("non-EM primary") {
    auto prim = stack.addParticle(std::make_tuple(Code::Proton, 10_GeV,
                                                  DirectionVector(rootCS, {1, 0, 0}),
                                                  Point(rootCS, 0_m, 0_m, 0_m), 0_ns));
    test::StackView view{prim};
    auto projectile = view.getProjectile();

    auto sec1 = projectile.addSecondary(
        std::make_tuple(Code::PiPlus, 2_GeV, DirectionVector(rootCS, {1, 0, 0})));
    auto sec2 = projectile.addSecondary(
        std::make_tuple(Code::KMinus, 8_GeV, DirectionVector(rootCS, {1, 0, 0})));

    EMThinning thinning{100_GeV, 1e20};
    thinning.doSecondaries(view);
    REQUIRE(view.getEntries() == 2);
    REQUIRE(sec1.getWeight() == 1.);
    REQUIRE(sec2.getWeight() == 1.);
  }

  SECTION("non-EM interaction") {
    auto prim = stack.addParticle(std::make_tuple(Code::Photon, 10_GeV,
                                                  DirectionVector(rootCS, {1, 0, 0}),
                                                  Point(rootCS, 0_m, 0_m, 0_m), 0_ns));
    test::StackView view{prim};
    auto projectile = view.getProjectile();

    // mimic photohadronic interaction, more than 2 secondaries
    auto sec1 = projectile.addSecondary(
        std::make_tuple(Code::PiPlus, 2_GeV, DirectionVector(rootCS, {1, 0, 0})));
    auto sec2 = projectile.addSecondary(
        std::make_tuple(Code::KMinus, 3_GeV, DirectionVector(rootCS, {1, 0, 0})));
    auto sec3 = projectile.addSecondary(
        std::make_tuple(Code::Proton, 5_GeV, DirectionVector(rootCS, {1, 0, 0})));

    EMThinning thinning{100_GeV, 1e20};
    thinning.doSecondaries(view);
    REQUIRE(view.getEntries() == 3);
    REQUIRE(sec1.getWeight() == 1.);
    REQUIRE(sec2.getWeight() == 1.);
    REQUIRE(sec3.getWeight() == 1.);
  }

  SECTION("above threshold") {
    auto prim = stack.addParticle(std::make_tuple(Code::Photon, 10_GeV,
                                                  DirectionVector(rootCS, {1, 0, 0}),
                                                  Point(rootCS, 0_m, 0_m, 0_m), 0_ns));

    test::StackView view{prim};
    auto projectile = view.getProjectile();

    auto sec1 = projectile.addSecondary(
        std::make_tuple(Code::Electron, 2_GeV, DirectionVector(rootCS, {1, 0, 0})));
    auto sec2 = projectile.addSecondary(
        std::make_tuple(Code::Positron, 8_GeV, DirectionVector(rootCS, {1, 0, 0})));

    EMThinning thinning{1_MeV, 1e20};
    thinning.doSecondaries(view);
    REQUIRE(view.getEntries() == 2);
    REQUIRE(sec1.getWeight() == 1.);
    REQUIRE(sec2.getWeight() == 1.);
  }

  SECTION("below threshold") {
    EMThinning thinning{100_GeV, 1e20};

    boost::accumulators::accumulator_set<
        double, boost::accumulators::stats<boost::accumulators::tag::mean>>
        acc;
    boost::accumulators::accumulator_set<
        HEPEnergyType, boost::accumulators::stats<boost::accumulators::tag::mean>>
        accE;

    // statistical test, 10000 iterations
    for (int i = 0; i < 10'000; ++i) {
      stack.clear();
      auto prim = stack.addParticle(std::make_tuple(Code::Photon, 10_GeV,
                                                    DirectionVector(rootCS, {1, 0, 0}),
                                                    Point(rootCS, 0_m, 0_m, 0_m), 0_ns));

      test::StackView view{prim};
      auto projectile = view.getProjectile();

      projectile.addSecondary(
          std::make_tuple(Code::Electron, 2_GeV, DirectionVector(rootCS, {1, 0, 0})));
      projectile.addSecondary(
          std::make_tuple(Code::Positron, 8_GeV, DirectionVector(rootCS, {1, 0, 0})));

      thinning.doSecondaries(view);
      REQUIRE(view.getEntries() == 1);

      double sumWeight{};
      HEPEnergyType sumE{};

      for (auto& p : view) {
        sumWeight += p.getWeight();
        sumE += p.getWeight() * p.getEnergy();
      }

      acc(sumWeight);
      accE(sumE);
    }

    REQUIRE(boost::accumulators::mean(acc) == Approx(2).epsilon(2e-2));
    REQUIRE(boost::accumulators::mean(accE) / 1_GeV == Approx(10).epsilon(2e-2));
  }

  SECTION("weight limitation") {
    double constexpr maxWeight = 3; // min acceptance prob. >= 0.3333
    EMThinning thinning{100_GeV, maxWeight};

    boost::accumulators::accumulator_set<
        double, boost::accumulators::stats<boost::accumulators::tag::mean>>
        acc;
    boost::accumulators::accumulator_set<
        HEPEnergyType, boost::accumulators::stats<boost::accumulators::tag::mean>>
        accE;

    // statistical test, 10000 iterations
    for (int i = 0; i < 10'000; ++i) {
      stack.clear();
      auto prim = stack.addParticle(std::make_tuple(Code::Photon, 10_GeV,
                                                    DirectionVector(rootCS, {1, 0, 0}),
                                                    Point(rootCS, 0_m, 0_m, 0_m), 0_ns));

      test::StackView view{prim};
      auto projectile = view.getProjectile();

      projectile.addSecondary(
          std::make_tuple(Code::Electron, 2_GeV, DirectionVector(rootCS, {1, 0, 0})));
      projectile.addSecondary(
          std::make_tuple(Code::Positron, 8_GeV, DirectionVector(rootCS, {1, 0, 0})));

      thinning.doSecondaries(view);

      double sumWeight{};
      HEPEnergyType sumE{};

      for (auto& p : view) {
        auto const w = p.getWeight();
        REQUIRE(w <= maxWeight);
        sumWeight += w;
        sumE += p.getEnergy() * w;
      }

      acc(sumWeight);
      accE(sumE);
    }

    REQUIRE(boost::accumulators::mean(acc) == Approx(2).epsilon(2e-2));
    REQUIRE(boost::accumulators::mean(accE) / 1_GeV == Approx(10).epsilon(2e-2));
  }

  SECTION("max. weight reached") {
    double constexpr maxWeight = 3; // min acceptance prob. >= 0.3333
    EMThinning thinning{100_GeV, maxWeight};

    for (int i = 0; i < 10'000; ++i) {
      stack.clear();
      auto prim = stack.addParticle(std::make_tuple(Code::Photon, 10_GeV,
                                                    DirectionVector(rootCS, {1, 0, 0}),
                                                    Point(rootCS, 0_m, 0_m, 0_m), 0_ns));
      prim.setWeight(3.5);

      test::StackView view{prim};
      auto projectile = view.getProjectile();

      auto sec1 = projectile.addSecondary(
          std::make_tuple(Code::Electron, 2_GeV, DirectionVector(rootCS, {1, 0, 0})));
      auto sec2 = projectile.addSecondary(
          std::make_tuple(Code::Positron, 8_GeV, DirectionVector(rootCS, {1, 0, 0})));

      // cannot perform thinning anymore, always keep all secondaries, weight unchanged
      thinning.doSecondaries(view);
      REQUIRE(view.getEntries() == 2);
      REQUIRE(sec1.getWeight() == 3.5);
      REQUIRE(sec2.getWeight() == 3.5);
    }
  }
}
