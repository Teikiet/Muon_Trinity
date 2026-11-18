/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/stack/DummyStack.hpp>

using namespace corsika;

#include <catch2/catch_all.hpp>

#include <tuple>

TEST_CASE("DummyStack", "[stack]") {

  logging::set_level(logging::level::info);

  using TestStack = dummy_stack::DummyStack;

  dummy_stack::NoData noData;

  SECTION("write node") {

    TestStack s;
    s.addParticle(std::tuple<dummy_stack::NoData>{noData});
    CHECK(s.getEntries() == 1);
  }

  SECTION("stack fill and cleanup") {

    TestStack s;
    // add 99 particles, each 10th particle is a nucleus with A=i and Z=A/2!

    for (int i = 0; i < 99; ++i) {
      s.addParticle(std::tuple<dummy_stack::NoData>{noData});
    }

    CHECK(s.getEntries() == 99);
    for (int i = 0; i < 99; ++i) s.getNextParticle().erase();
    CHECK(s.getEntries() == 0);
  }
}
