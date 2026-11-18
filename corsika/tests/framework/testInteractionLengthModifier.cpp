/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/process/InteractionLengthModifier.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;

struct DummyProcess : public InteractionProcess<DummyProcess> {
  int id{0};

  template <typename TParticle>
  GrammageType getInteractionLength(TParticle const&) {
    return 100_g / 1_cm / 1_cm;
  }

  template <typename TArgument>
  void doInteraction(TArgument& arg) {
    arg = 3;
  }

  DummyProcess() = default;
  DummyProcess(DummyProcess&&) = default;

  // prevent copying
  DummyProcess(DummyProcess const&) = delete;
  DummyProcess& operator=(DummyProcess const&) = delete;
};

struct DummyParticle {
  HEPEnergyType getEnergy() const { return 47_TeV; }

  Code getPID() const { return Code::MuPlus; }
};

TEST_CASE("InteractionLengthModifier", "[process]") {
  DummyProcess u;
  u.id = 38;

  auto const modifier = [](GrammageType orig, Code, HEPEnergyType) -> GrammageType {
    return orig * 2;
  };

  InteractionLengthModifier mod{std::move(u), modifier};
  REQUIRE(std::is_same_v<decltype(mod), InteractionLengthModifier<DummyProcess>>);

  SECTION("getInteractionLength") {
    DummyParticle const p;
    REQUIRE(mod.getInteractionLength(p) == 200_g / 1_cm / 1_cm);
  }

  SECTION("doInteraction") {
    int k = 0;
    mod.doInteraction(k);
    REQUIRE(k == 3);
  }

  SECTION("getProcess") {
    DummyProcess& uRef = mod.getProcess();
    REQUIRE(uRef.id == 38);

    decltype(mod) const& modConstRef = mod;
    DummyProcess const& uConstRef = modConstRef.getProcess();
    REQUIRE(uConstRef.id == 38);
  }
}
