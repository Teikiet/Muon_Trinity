/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/Pythia8.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/random/RNGManager.hpp>

#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using Catch::Approx;

using DummyEnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

TEST_CASE("Pythia8", "modules") {

  logging::set_level(logging::level::info);

  SECTION("linking pythia") {
    using namespace Pythia8;

    // Generator. Process selection. LHC initialization. Histogram.
    Pythia pythia;

    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    pythia.readString("ProcessLevel:all = off");

    pythia.init();

    Event& event = pythia.event;
    event.reset();

    pythia.particleData.mayDecay(321, true);
    double pz = 100.;
    double m = 0.49368;
    event.append(321, 1, 0, 0, 0., 0., 100., sqrt(pz * pz + m * m), m);

    if (!pythia.next())
      CORSIKA_LOG_CRITICAL("decay failed!");
    else
      CORSIKA_LOG_DEBUG("particles after decay: {}", event.size());
    event.list();

    // loop over final state
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal()) {
        CORSIKA_LOG_DEBUG("particle: id= {}", pythia.event[i].id());
      }
  }

  SECTION("pythia interface") {

    std::set<Code> const particleList = {Code::PiPlus, Code::PiMinus, Code::KPlus,
                                         Code::KMinus, Code::K0Long,  Code::K0Short};
    RNGManager<>::getInstance().registerRandomStream("pythia");
    corsika::pythia8::Decay decay(particleList);
  }
}

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

using namespace corsika;

template <typename TStackView>
auto sumMomentum(TStackView const& view, CoordinateSystem const& vCS) {
  MomentumVector sum{vCS, 0_eV, 0_eV, 0_eV};
  for (auto const& p : view) { sum += p.getMomentum(); }
  return sum;
}

TEST_CASE("Pythia8Interface", "modules") {

  logging::set_level(logging::level::info);

  auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Proton);
  auto const& cs = *csPtr;
  {
    [[maybe_unused]] auto const& env_dummy = env;
    [[maybe_unused]] auto const& node_dummy = nodePtr;
  }

  SECTION("pythia decay") {
    HEPEnergyType const P0 = 10_GeV;
    auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
        Code::PiPlus, P0, (DummyEnvironment::BaseNodeType* const)nodePtr, cs);
    auto& stack = *stackPtr;
    auto& view = *secViewPtr;

    auto const& particle = stack.getNextParticle();
    auto const plab = MomentumVector(cs, {P0, 0_GeV, 0_GeV});

    std::set<Code> const particleList = {Code::PiPlus, Code::PiMinus, Code::KPlus,
                                         Code::KMinus, Code::K0Long,  Code::K0Short};

    RNGManager<>::getInstance().registerRandomStream("pythia");

    corsika::pythia8::Decay decay(particleList);

    CORSIKA_LOG_INFO("stack: {} {}", stack.asString(), particle.asString());
    [[maybe_unused]] const TimeType time = decay.getLifetime(particle);
    double const gamma = particle.getEnergy() / get_mass(Code::PiPlus);
    REQUIRE(time == get_lifetime(Code::PiPlus) * gamma);
    decay.doDecay(*secViewPtr);
    CORSIKA_LOG_INFO("piplus->{}", stack.asString());
    REQUIRE(stack.getEntries() == 3); // piplus, muplu, numu
    auto const pSum = sumMomentum(view, cs);
    REQUIRE((pSum - plab).getNorm() / 1_GeV == Approx(0).margin(1e-4));
    REQUIRE((pSum.getNorm() - plab.getNorm()) / 1_GeV == Approx(0).margin(1e-4));
  }

  SECTION("pythia decay config") {
    corsika::pythia8::Decay decay({Code::PiPlus, Code::PiMinus});
    REQUIRE(decay.isDecayHandled(Code::PiPlus));
    REQUIRE(decay.isDecayHandled(Code::PiMinus));
    REQUIRE_FALSE(decay.isDecayHandled(Code::KPlus));

    const std::vector<Code> particleTestList = {Code::PiPlus, Code::PiMinus, Code::KPlus,
                                                Code::LambdaBar, Code::D0Bar};

    // setup decays
    decay.setHandleDecay(particleTestList);
    for (auto& pCode : particleTestList) REQUIRE(decay.isDecayHandled(pCode));

    // individually
    decay.setHandleDecay(Code::KMinus);

    // impossible
    REQUIRE_THROWS(decay.setHandleDecay(Code::Photon));

    REQUIRE(decay.isDecayHandled(Code::PiPlus));
    REQUIRE_FALSE(decay.isDecayHandled(Code::Photon));

    // possible decays
    REQUIRE_FALSE(decay.canHandleDecay(Code::Photon));
    REQUIRE_FALSE(decay.canHandleDecay(Code::Proton));
    REQUIRE_FALSE(decay.canHandleDecay(Code::Electron));
    REQUIRE(decay.canHandleDecay(Code::PiPlus));
    REQUIRE(decay.canHandleDecay(Code::MuPlus));
  }

  SECTION("pythia neutrino interaction") {

    // this will be a nu-p collision at sqrts=3.5TeV -> no problem for pythia
    auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
        Code::NuE, 7_TeV, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
    auto& view = *secViewPtr;

    corsika::pythia8::NeutrinoInteraction collision({}, true, false);

    REQUIRE(collision.isValid(
        Code::NuE, Code::Proton,
        {calculate_total_energy(100_TeV, NuE::mass), {0_eV, 0_eV, 100_TeV}},
        {Nitrogen::mass, {0_eV, 0_eV, 0_eV}}));

    // non-zero cross section for neutrino projectile
    REQUIRE(collision.getCrossSection(
                Code::NuE, Code::Proton, {100_GeV, {0_eV, 0_eV, 100_GeV}},
                {Proton::mass, {0_eV, 0_eV, 0_eV}}) != CrossSectionType::zero());

    HEPEnergyType const P0 = 100_PeV;
    collision.doInteraction(view, Code::NuE, Code::Nitrogen,
                            {calculate_total_energy(P0, NuE::mass), {0_eV, 0_eV, P0}},
                            {Nitrogen::mass, {0_eV, 0_eV, 0_eV}});
    REQUIRE(view.getSize() >= 2);

    // no interaction for protons
    REQUIRE(collision.getCrossSection(
                Code::Proton, Code::Proton, {100_GeV, {0_eV, 0_eV, 100_GeV}},
                {Proton::mass, {0_eV, 0_eV, 0_eV}}) == CrossSectionType::zero());
  }

  SECTION("single-use neutrino interaction") {
    auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
        Code::NuE, 7_TeV, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
    auto& view = *secViewPtr;

    RNGManager<>::getInstance().registerRandomStream("pythia");

    corsika::pythia8::SingleUseNeutrinoInteraction collision({}, true, false);

    FourMomentum const projectileP4{calculate_total_energy(100_PeV, NuE::mass),
                                    {cs, {0_eV, 0_eV, 100_PeV}}};
    FourMomentum const targetP4{Nitrogen::mass, {cs, {0_eV, 0_eV, 0_eV}}};

    REQUIRE(collision.getCrossSection(Code::NuE, Code::Nitrogen, projectileP4,
                                      targetP4) != CrossSectionType::zero());

    collision.doInteraction(view, Code::NuE, Code::Nitrogen, projectileP4, targetP4);
    REQUIRE(view.getSize() >= 2);

    REQUIRE(collision.getCrossSection(Code::NuE, Code::Nitrogen, projectileP4,
                                      targetP4) == CrossSectionType::zero());

    collision.reset();
    REQUIRE(collision.getCrossSection(Code::NuE, Code::Nitrogen, projectileP4,
                                      targetP4) != CrossSectionType::zero());
  }

#include <tests/modules/testPythia8Interaction.inl>
}