/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/stack/CombinedStack.hpp>
#include <corsika/stack/DummyStack.hpp>
#include <corsika/stack/GeometryNodeStackExtension.hpp>

using namespace corsika;

#include <catch2/catch_all.hpp>

#include <iostream>
using namespace std;

// this is our dummy environment, it only knows its trivial BaseNodeType
class DummyEnv {
public:
  typedef int BaseNodeType;
};

// the GeometryNode stack needs to know the type of geometry-nodes from the DummyEnv:
template <typename TStackIter>
using DummyGeometryDataInterface =
    typename node::MakeGeometryDataInterface<TStackIter, DummyEnv>::type;

// combine dummy stack with geometry information for tracking
template <typename TStackIter>
using StackWithGeometryInterface =
    CombinedParticleInterface<dummy_stack::DummyStack::pi_type,
                              DummyGeometryDataInterface, TStackIter>;

using TestStack = CombinedStack<typename dummy_stack::DummyStack::stack_data_type,
                                node::GeometryData<DummyEnv>, StackWithGeometryInterface>;

TEST_CASE("GeometryNodeStackExtension", "stack") {

  logging::set_level(logging::level::info);

  dummy_stack::NoData noData;

  SECTION("write node") {

    const int data = 5;

    TestStack s;
    s.addParticle(std::make_tuple(noData), std::tuple<const int*>{&data});

    CHECK(s.getEntries() == 1);
  }

  SECTION("write/read node") {
    const int data = 15;

    TestStack s;
    auto p = s.addParticle(std::make_tuple(noData));
    p.setNode(&data);
    CHECK(s.getEntries() == 1);

    const auto pout = s.getNextParticle();
    CHECK(*(pout.getNode()) == 15);
  }

  SECTION("stack fill and cleanup") {

    const int data = 16;

    TestStack s;
    // add 99 particles, each 10th particle is a nucleus with A=i and Z=A/2!
    for (int i = 0; i < 99; ++i) {
      auto p = s.addParticle(std::tuple<dummy_stack::NoData>{noData});
      p.setNode(&data);
    }

    CHECK(s.getEntries() == 99);
    double v = 0;
    for (int i = 0; i < 99; ++i) {
      auto p = s.getNextParticle();
      v += *(p.getNode());
      p.erase();
    }
    CHECK(v == 99 * data);
    CHECK(s.getEntries() == 0);
  }

  SECTION("copy and swap") {

    const int data1 = 16;
    const int data2 = 17;

    TestStack s;
    // add 99 particles, each 10th particle is a nucleus with A=i and Z=A/2!
    for (int i = 0; i < 4; ++i) {
      auto p = s.addParticle(std::tuple<dummy_stack::NoData>{noData});
      p.setNode(i % 2 == 0 ? &data1 : &data2);
    }

    CHECK(*((s.first() + 0)).getNode() == 16);
    CHECK(*((s.first() + 1)).getNode() == 17);
    CHECK(*((s.first() + 2)).getNode() == 16);
    CHECK(*((s.first() + 3)).getNode() == 17);

    s.copy(s.first() + 0, s.first() + 1);
    CHECK(*((s.first() + 0)).getNode() == 16);
    CHECK(*((s.first() + 1)).getNode() == 16);

    s.swap(s.first() + 3, s.first() + 1);
    CHECK(*((s.first() + 1)).getNode() == 17);
    CHECK(*((s.first() + 3)).getNode() == 16);
  }
}
