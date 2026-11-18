/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/stack/history/Event.hpp>
#include <corsika/stack/history/HistorySecondaryProducer.hpp>
#include <corsika/stack/history/HistoryStackExtension.hpp>
#include <corsika/stack/history/HistorySecondaryProducer.hpp>
#include <corsika/stack/DummyStack.hpp>
#include <corsika/framework/stack/CombinedStack.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;

// the GeometryNode stack needs to know the type of geometry-nodes from the DummyEnv:
template <typename TStackIter>
using DummyHistoryDataInterface =
    typename history::MakeHistoryDataInterface<TStackIter, history::Event>::type;

// combine dummy stack with geometry information for tracking
template <typename TStackIter>
using StackWithHistoryInterface =
    CombinedParticleInterface<dummy_stack::DummyStack::pi_type, DummyHistoryDataInterface,
                              TStackIter>;

using TestStack =
    CombinedStack<typename dummy_stack::DummyStack::stack_data_type,
                  history::HistoryData<history::Event>, StackWithHistoryInterface,
                  history::HistorySecondaryProducer>;

// the correct secondary stack view
using TestStackView = typename TestStack::stack_view_type;
using EvtPtr = std::shared_ptr<history::Event>;

template <typename Event>
int count_generations(Event const* event) {
  int genCounter = 0;
  while (event) {
    event = event->parentEvent().get();
    genCounter++;
  }

  return genCounter;
}

TEST_CASE("HistoryStackExtensionView", "[stack]") {

  logging::set_level(logging::level::info);

  // in this test we only use one singel stack !
  const dummy_stack::NoData noData;
  TestStack stack;

  // add primary particle
  auto p0 = stack.addParticle(std::make_tuple(noData));

  CHECK(stack.getEntries() == 1);
  EvtPtr evt = p0.getEvent();
  CHECK(evt == nullptr);
  CHECK(count_generations(evt.get()) == 0);

  SECTION("interface test, view") {

    // add secondaries, 1st generation
    auto pnext = stack.getNextParticle();
    TestStackView hview0(pnext);

    auto const ev0 = p0.getEvent();
    CHECK(ev0 == nullptr);

    CORSIKA_LOG_DEBUG("loop VIEW");

    // add 5 secondaries
    for (int i = 0; i < 5; ++i) {
      auto sec = hview0.addSecondary(std::make_tuple(noData));

      CHECK(sec.getParentEventIndex() == i);
      CHECK(sec.getEvent().get() != nullptr);
      CHECK(sec.getEvent()->parentEvent() == nullptr);
      CHECK(count_generations(sec.getEvent().get()) == 1);
    }

    // read 1st genertion particle particle
    auto p1 = stack.getNextParticle();
    CHECK(count_generations(p1.getEvent().get()) == 1);

    TestStackView hview1(p1);

    auto const ev1 = p1.getEvent();

    // add second generation of secondaries
    // add 10 secondaries
    for (int i = 0; i < 10; ++i) {
      auto sec = hview1.addSecondary(std::make_tuple(noData));

      CHECK(sec.getParentEventIndex() == i);
      CHECK(sec.getEvent()->parentEvent() == ev1);
      CHECK(sec.getEvent()->parentEvent()->parentEvent() == ev0);

      CHECK(count_generations(sec.getEvent().get()) == 2);

      const auto org_projectile = stack.at(sec.getEvent()->projectileIndex());
      CHECK(org_projectile.getEvent() == sec.getEvent()->parentEvent());
    }

    // read 2nd genertion particle particle
    auto p2 = stack.getNextParticle();

    TestStackView hview2(p2);

    auto const ev2 = p2.getEvent();

    // add third generation of secondaries
    // add 15 secondaries
    for (int i = 0; i < 15; ++i) {
      CORSIKA_LOG_TRACE("loop, view: " + std::to_string(i));

      auto sec = hview2.addSecondary(std::make_tuple(noData));
      CORSIKA_LOG_TRACE("loop, ---- ");

      CHECK(sec.getParentEventIndex() == i);
      CHECK(sec.getEvent()->parentEvent() == ev2);
      CHECK(sec.getEvent()->parentEvent()->parentEvent() == ev1);
      CHECK(sec.getEvent()->parentEvent()->parentEvent()->parentEvent() == ev0);

      CHECK(count_generations(sec.getEvent().get()) == 3);
    }
  }

  SECTION("also test projectile access") {

    CORSIKA_LOG_TRACE("projectile test");

    // add secondaries, 1st generation
    TestStackView hview0(p0);
    auto proj0 = hview0.getProjectile();
    auto const ev0 = p0.getEvent();
    CHECK(ev0 == nullptr);

    CORSIKA_LOG_TRACE("loop");

    // add 5 secondaries
    for (int i = 0; i < 5; ++i) {
      CORSIKA_LOG_TRACE("loop " + std::to_string(i));
      auto sec = proj0.addSecondary(std::make_tuple(noData));

      CHECK(sec.getParentEventIndex() == i);
      CHECK(sec.getEvent() != nullptr);
      CHECK(sec.getEvent()->parentEvent() == nullptr);
    }
    CHECK(stack.getEntries() == 6);

    // read 1st genertion particle particle
    auto p1 = stack.getNextParticle();

    TestStackView hview1(p1);
    auto proj1 = hview1.getProjectile();
    auto const ev1 = p1.getEvent();

    // add second generation of secondaries
    // add 10 secondaries
    for (unsigned int i = 0; i < 10; ++i) {
      auto sec = proj1.addSecondary(std::make_tuple(noData));

      CHECK(sec.getParentEventIndex() == int(i));
      CHECK(sec.getEvent()->parentEvent() == ev1);
      CHECK(sec.getEvent()->secondaries().size() == i + 1);
      CHECK(sec.getEvent()->parentEvent()->parentEvent() == ev0);

      const auto org_projectile = stack.at(sec.getEvent()->projectileIndex());
      CHECK(org_projectile.getEvent() == sec.getEvent()->parentEvent());
    }
    CHECK(stack.getEntries() == 16);
  }
}
