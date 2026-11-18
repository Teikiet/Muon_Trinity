/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/process/DynamicInteractionProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using Catch::Approx;

struct DummyStack {
  using stack_view_type = int&;
};

struct DummyProcessA : public InteractionProcess<DummyProcessA> {
  int id_;

  DummyProcessA(int id)
      : id_{id} {}

  // mockup to "identify" the process via its ID
  CrossSectionType getCrossSection(Code, Code, FourMomentum const&,
                                   FourMomentum const&) const {
    return id_ * 10_mb;
  }

  // mockup: set argument to our ID
  template <typename TArgument>
  void doInteraction(TArgument& argument, Code, Code, FourMomentum const&,
                     FourMomentum const&) {
    argument = id_ * 10;
  }

  // to make sure we can test proper handling of dangling references and
  // stuff like that
  ~DummyProcessA() { id_ = 0; }

  // prevent copying
  DummyProcessA(DummyProcessA const&) = delete;
  DummyProcessA& operator=(DummyProcessA const&) = delete;

  // allow moving
  DummyProcessA(DummyProcessA&&) = default;
  DummyProcessA& operator=(DummyProcessA&&) = default;
};

// same as above to have another class
struct DummyProcessB : public InteractionProcess<DummyProcessB> {
  int id_;

  DummyProcessB(int id)
      : id_{id} {}

  CrossSectionType getCrossSection(Code, Code, FourMomentum const&,
                                   FourMomentum const&) const {
    return id_ * 1_mb;
  }

  template <typename TArgument>
  void doInteraction(TArgument& argument, Code, Code, FourMomentum const&,
                     FourMomentum const&) {
    argument = id_;
  }

  // to make sure we can test proper handling of dangling references and
  // stuff like that
  ~DummyProcessB() { id_ = 0; }

  // prevent copying
  DummyProcessB(DummyProcessB const&) = delete;
  DummyProcessB& operator=(DummyProcessB const&) = delete;

  // allow moving
  DummyProcessB(DummyProcessB&&) = default;
  DummyProcessB& operator=(DummyProcessB&&) = default;
};

TEST_CASE("DynamicInteractionProcess", "[process]") {
  FourMomentum const fourVec{1_GeV, {1_GeV, 1_GeV, 1_GeV}};

  std::shared_ptr<DummyProcessA> a = std::make_shared<DummyProcessA>(1);
  std::shared_ptr<DummyProcessB> b = std::make_shared<DummyProcessB>(2);

  SECTION("getCrossSection") {
    DynamicInteractionProcess<DummyStack> dynamic;

    dynamic = a;
    REQUIRE(dynamic.getCrossSection(Code::Electron, Code::Electron, fourVec, fourVec) /
                10_mb ==
            Approx(1).epsilon(1e-6));

    dynamic = b;
    REQUIRE(dynamic.getCrossSection(Code::Electron, Code::Electron, fourVec, fourVec) /
                2_mb ==
            Approx(1).epsilon(1e-6));

    dynamic = DynamicInteractionProcess<DummyStack>{std::make_shared<DummyProcessA>(3)};
    REQUIRE(dynamic.getCrossSection(Code::Electron, Code::Electron, fourVec, fourVec) /
                30_mb ==
            Approx(1).epsilon(1e-6));

    dynamic = std::make_shared<DummyProcessB>(4);
    REQUIRE(dynamic.getCrossSection(Code::Electron, Code::Electron, fourVec, fourVec) /
                4_mb ==
            Approx(1).epsilon(1e-6));

    // test process going out of scope
    { dynamic = std::make_shared<DummyProcessB>(5); }
    REQUIRE(dynamic.getCrossSection(Code::Electron, Code::Electron, fourVec, fourVec) /
                5_mb ==
            Approx(1).epsilon(1e-6));
  }

  SECTION("doInteraction") {
    int value{};

    DynamicInteractionProcess<DummyStack> dynamic;

    dynamic = a;
    dynamic.doInteraction(value, Code::Electron, Code::Electron, fourVec, fourVec);
    REQUIRE(value == 10);

    dynamic = b;
    dynamic.doInteraction(value, Code::Electron, Code::Electron, fourVec, fourVec);
    REQUIRE(value == 2);

    dynamic = DynamicInteractionProcess<DummyStack>{std::make_shared<DummyProcessA>(3)};
    dynamic.doInteraction(value, Code::Electron, Code::Electron, fourVec, fourVec);
    REQUIRE(value == 30);

    dynamic = std::make_shared<DummyProcessB>(4);
    dynamic.doInteraction(value, Code::Electron, Code::Electron, fourVec, fourVec);
    REQUIRE(value == 4);
  }
}
