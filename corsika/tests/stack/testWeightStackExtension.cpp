/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/stack/CombinedStack.hpp>
#include <corsika/stack/DummyStack.hpp>
#include <corsika/stack/WeightStackExtension.hpp>

using namespace corsika;

#include <catch2/catch_all.hpp>

#include <iostream>
using namespace std;
using Catch::Approx;

// the Weight stack:
template <typename TStackIter>
using DummyWeightDataInterface =
    typename weights::MakeWeightDataInterface<TStackIter>::type;

// combine dummy stack with geometry information for tracking
template <typename TStackIter>
using StackWithGeometryInterface =
    CombinedParticleInterface<dummy_stack::DummyStack::pi_type, DummyWeightDataInterface,
                              TStackIter>;

using TestStack = CombinedStack<typename dummy_stack::DummyStack::stack_data_type,
                                weights::WeightData, StackWithGeometryInterface>;

TEST_CASE("WeightStackExtension", "[stack]") {

  logging::set_level(logging::level::info);

  dummy_stack::NoData noData;

  SECTION("write weights") {

    double const weight = 5.1;

    TestStack s;
    auto const p0 = s.addParticle(std::make_tuple(noData), std::tuple<double>{weight});

    CHECK(s.getEntries() == 1);
    CHECK(p0.getWeight() == weight);
  }

  SECTION("write/read weights") {
    double const weight = 15;

    TestStack s;
    auto p = s.addParticle(std::make_tuple(noData));
    p.setWeight(weight);
    CHECK(s.getEntries() == 1);

    const auto pout = s.getNextParticle();
    CHECK(pout.getWeight() == weight);
  }

  SECTION("stack fill and cleanup") {

    double const weight = 16;

    TestStack s;
    // add 99 particles, each 10th particle has weigth 100
    for (int i = 0; i < 99; ++i) {
      auto p = s.addParticle(std::tuple<dummy_stack::NoData>{noData});
      p.setWeight(((i + 1) % 10 == 0) ? 100 : weight);
    }

    CHECK(s.getEntries() == 99);
    double v = 0;
    for (int i = 0; i < 99; ++i) {
      auto p = s.getNextParticle();
      v += p.getWeight();
      p.erase();
    }
    CHECK(v == Approx(90 * weight + 9 * 100));
    CHECK(s.getEntries() == 0);
  }

  SECTION("stack operations") {

    TestStack s;

    double const weight = 16;
    // add 99 particles, each 10th particle is a positron
    // i=9, 19, 29, etc. are positron
    // add 99 particles, each 10th particle has weigth 100
    for (int i = 0; i < 99; ++i) {
      auto p = s.addParticle(std::tuple<dummy_stack::NoData>{noData});
      p.setWeight(((i + 1) % 10 == 0) ? 100 : weight);
    }

    // copy
    {
      s.copy(s.begin() + 9, s.begin() + 10); // 100 to 16
      const auto& p9 = s.cbegin() + 9;
      const auto& p10 = s.cbegin() + 10;

      CHECK(p9.getWeight() == 100);
      CHECK(p10.getWeight() == 100);
    }

    // copy
    {
      s.copy(s.begin() + 93, s.begin() + 9); // 16 to 100
      const auto& p93 = s.cbegin() + 93;
      const auto& p9 = s.cbegin() + 9;

      CHECK(p9.getWeight() == 16);
      CHECK(p93.getWeight() == 16);
    }

    // copy
    {
      s.copy(s.begin() + 89, s.begin() + 79); // 100 to 100
      const auto& p89 = s.cbegin() + 89;
      const auto& p79 = s.cbegin() + 79;

      CHECK(p89.getWeight() == 100);
      CHECK(p79.getWeight() == 100);
    }

    // invalid copy
    { CHECK_THROWS(s.copy(s.begin(), s.end() + 1000)); }

    // swap
    {
      s.swap(s.begin() + 11, s.begin() + 10);
      const auto& p11 = s.cbegin() + 11; // now: positron
      const auto& p10 = s.cbegin() + 10; // now: electron

      CHECK(p11.getWeight() == 100);
      CHECK(p10.getWeight() == 16);
    }

    // swap two positrons
    {
      s.swap(s.begin() + 29, s.begin() + 59);
      const auto& p29 = s.cbegin() + 29;
      const auto& p59 = s.cbegin() + 59;

      CHECK(p29.getWeight() == 100);
      CHECK(p59.getWeight() == 100);
    }

    for (int i = 0; i < 99; ++i) s.last().erase();
    CHECK(s.getEntries() == 0);
  }
}
