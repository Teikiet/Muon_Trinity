/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace corsika::media;
using Catch::Approx;

struct DummyRNG {
  double v_;
  DummyRNG(double v)
      : v_(v) {}
  static constexpr int max() { return 10; }
  static constexpr int min() { return 0; }
  double operator()() const { return v_; }
};

TEST_CASE("NuclearComposition") {

  logging::set_level(logging::level::info);

  // incompatible input: wrong vectors
  CHECK_THROWS(media::NuclearComposition({Code::Oxygen, Code::Carbon},
                                         {0.20, 0.05, 1 - 0.20 - 0.05}));
  // incompatible input: wrong fractions
  CHECK_THROWS(media::NuclearComposition({Code::Oxygen, Code::Carbon},
                                         {0.21, 0.05, 1 - 0.20 - 0.05}));
  // incompatible input: wrong fractions
  CHECK_THROWS(media::NuclearComposition({Code::Oxygen, Code::Carbon},
                                         {0.19, 0.05, 1 - 0.20 - 0.05}));

  media::NuclearComposition const testComposition(
      {Code::Oxygen, Code::Carbon, Code::Nitrogen}, {0.20, 0.05, 1 - 0.20 - 0.05});

  CHECK(testComposition.getSize() == 3);
  CHECK(testComposition.getFractions() == std::vector<double>{0.2, 0.05, 1 - 0.2 - 0.05});
  CHECK(testComposition.getComponents() ==
        std::vector<Code>{Code::Oxygen, Code::Carbon, Code::Nitrogen});

  CHECK(testComposition.getHash() ==
        18183071370474897160U); // we need a stable hasing algorithm
  CHECK(testComposition.getAverageMassNumber() == 14.3);

  CHECK(testComposition.getWeighted([](Code) -> double { return 1; }) ==
        std::vector<double>{0.2, 0.05, 1 - 0.2 - 0.05});

  std::vector<CrossSectionType> const testCX =
      testComposition.getWeighted([](Code) -> CrossSectionType { return 1_mb; });
  std::vector<CrossSectionType> const checkCX{0.2_mb, 0.05_mb, 1_mb - 0.2_mb - 0.05_mb};
  for (auto i1 = testCX.begin(), i2 = checkCX.begin(); i1 != testCX.end(); ++i1, ++i2) {
    CHECK(*i1 / 1_mb == Approx(*i2 / 1_mb));
  }

  CHECK(testComposition.getWeightedSum([](Code) -> double { return 1; }) == 1);

  CHECK(testComposition.getWeightedSum([](Code) -> CrossSectionType { return 1_mb; }) ==
        1_mb);

  CHECK(testComposition.sampleTarget(testCX, DummyRNG(0.1)) == Code::Oxygen);
}
