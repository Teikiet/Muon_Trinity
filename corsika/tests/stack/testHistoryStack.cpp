/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/stack/CombinedStack.hpp>
#include <corsika/stack/DummyStack.hpp>
#include <corsika/stack/history/HistoryStackExtension.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;

// this is our dummy environment, it only knows its trivial BaseNodeType
class DummyEvent {
private:
  size_t parent_;
  std::vector<int> secondaries_;

public:
  DummyEvent() {}
  DummyEvent(const size_t parent) { parent_ = parent; }

  size_t getParentIndex() { return parent_; }
  void addSecondary(const int particle) { secondaries_.push_back(particle); }
  int multiplicity() const { return secondaries_.size(); }
};

// the GeometryNode stack needs to know the type of geometry-nodes from the DummyEnv:
template <typename TStackIter>
using DummyHistoryDataInterface =
    typename history::MakeHistoryDataInterface<TStackIter, DummyEvent>::type;

// combine dummy stack with geometry information for tracking
template <typename TStackIter>
using StackWithHistoryInterface =
    CombinedParticleInterface<dummy_stack::DummyStack::pi_type, DummyHistoryDataInterface,
                              TStackIter>;

using TestStack =
    CombinedStack<typename dummy_stack::DummyStack::stack_data_type,
                  history::HistoryData<DummyEvent>, StackWithHistoryInterface>;

using EvtPtr = std::shared_ptr<DummyEvent>;

TEST_CASE("HistoryStackExtension", "stack") {

  logging::set_level(logging::level::info);

  const dummy_stack::NoData noData;
  TestStack s;

  auto p = s.addParticle(std::tuple<dummy_stack::NoData>{noData});

  SECTION("add lone particle") {
    CHECK(s.getEntries() == 1);

    EvtPtr evt = p.getEvent();
    CHECK(evt == nullptr);
  }

  SECTION("add and remove particles") {

    auto p = s.addParticle(std::tuple<dummy_stack::NoData>{noData});
    CHECK(s.getEntries() == 2);
    p.erase();
    CHECK(s.getEntries() == 1);
    s.purge();
    CHECK(s.getEntries() == 1);
  }
}
