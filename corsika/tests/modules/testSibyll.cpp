/*
 * (c) Copyright 2019 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/sibyll/ParticleConversion.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/COMBoost.hpp>

#include <SetupTestEnvironment.hpp>
#include <catch2/catch_all.hpp>
#include <tuple>

using namespace corsika;
using namespace corsika::sibyll;
using Catch::Approx;

using DummyEnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

TEST_CASE("Sibyll", "modules") {

  logging::set_level(logging::level::info);

  SECTION("Sibyll -> Corsika") {
    CHECK(Code::Electron ==
          corsika::sibyll::convertFromSibyll(corsika::sibyll::SibyllCode::Electron));
  }

  SECTION("Corsika -> Sibyll") {
    CHECK(corsika::sibyll::convertToSibyll(Electron::code) ==
          corsika::sibyll::SibyllCode::Electron);
    CHECK(corsika::sibyll::convertToSibyllRaw(Proton::code) == 13);
    CHECK(corsika::sibyll::convertToSibyll(XiC_2645_0::code) ==
          corsika::sibyll::SibyllCode::XiC_2645_0);
  }

  SECTION("canInteractInSibyll") {

    CHECK(corsika::sibyll::canInteract(Code::Proton));
    CHECK(corsika::sibyll::canInteract(Code::XiCPlus));

    CHECK_FALSE(corsika::sibyll::canInteract(Code::Electron));
    CHECK_FALSE(corsika::sibyll::canInteract(Code::SigmaC_2455_0));

    CHECK_FALSE(corsika::sibyll::canInteract(Code::Iron));
    CHECK_FALSE(corsika::sibyll::canInteract(Code::K0Bar));
    CHECK_FALSE(corsika::sibyll::canInteract(Code::Helium));
  }

  SECTION("cross-section type") {
    CHECK(corsika::sibyll::getSibyllXSCode(Code::Proton) == 1);
    CHECK(corsika::sibyll::getSibyllXSCode(Code::Electron) == 0);
    CHECK(corsika::sibyll::getSibyllXSCode(Code::K0Long) == 3);
    CHECK(corsika::sibyll::getSibyllXSCode(Code::K0Bar) == 0);
    CHECK(corsika::sibyll::getSibyllXSCode(Code::SigmaPlus) == 1);
    CHECK(corsika::sibyll::getSibyllXSCode(Code::PiMinus) == 2);
    CHECK(corsika::sibyll::getSibyllXSCode(Code::Helium) == 0);
  }

  SECTION("sibyll mass") {
    CHECK_FALSE(corsika::sibyll::getSibyllMass(Code::Electron) == 0_GeV);
    // Nucleus not a particle
    CHECK_THROWS(corsika::sibyll::getSibyllMass(Code::Iron));
    // Higgs not a particle in Sibyll
    CHECK_THROWS(corsika::sibyll::getSibyllMass(Code::H0));
  }
}

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>

#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>

template <typename TStackView>
auto sumMomentum(TStackView const& view, CoordinateSystem const& vCS) {
  Vector<hepenergy_d> sum{vCS, 0_eV, 0_eV, 0_eV};
  for (auto const& p : view) { sum += p.getMomentum(); }
  return sum;
}

TEST_CASE("SibyllInterface", "modules") {

  logging::set_level(logging::level::info);

  // the environment and stack should eventually disappear from here
  auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Oxygen);
  auto const& cs = *csPtr;
  { [[maybe_unused]] auto const& env_dummy = env; }

  auto [stack, viewPtr] = setup::testing::setup_stack(
      Code::Proton, 10_GeV, (DummyEnvironment::BaseNodeType* const)nodePtr, cs);
  test::StackView& view = *viewPtr;

  RNGManager<>::getInstance().registerRandomStream("sibyll");

  std::set<Code> const nuclearcomp = {Code::Hydrogen, Code::Nitrogen, Code::Oxygen,
                                      Code::Carbon};

  std::set<Code> stable_particles = {Code::PiPlus, Code::PiMinus};
  SECTION("InteractionInterface - valid targets") {

    corsika::sibyll::HadronInteractionModel model(stable_particles);
    // sibyll only accepts protons or nuclei with 4<=A<=18 as targets
    CHECK_FALSE(model.isValid(Code::Proton, Code::Electron, 100_GeV));
    CHECK(model.isValid(Code::Proton, Code::Hydrogen, 100_GeV));
    CHECK_FALSE(model.isValid(Code::Proton, Code::Deuterium, 100_GeV));
    CHECK(model.isValid(Code::Proton, Code::Helium, 100_GeV));
    CHECK_FALSE(model.isValid(Code::Proton, Code::Helium3, 100_GeV));
    CHECK_FALSE(model.isValid(Code::Proton, Code::Iron, 100_GeV));
    CHECK(model.isValid(Code::Proton, Code::Oxygen, 100_GeV));
    CHECK(model.isValid(Code::Rho0, Code::Oxygen, 100_GeV));
    // beam particles
    CHECK_FALSE(model.isValid(Code::Electron, Code::Oxygen, 100_GeV));
    CHECK_FALSE(model.isValid(Code::Iron, Code::Oxygen, 100_GeV));
    CHECK_FALSE(model.isValid(Code::K0Bar, Code::Oxygen, 100_GeV));
    // energy too low
    CHECK_FALSE(model.isValid(Code::Proton, Code::Proton, 9_GeV));
    CHECK(model.isValid(Code::Proton, Code::Proton, 11_GeV));
    // energy too high
    CHECK_FALSE(model.isValid(Code::Proton, Code::Proton, 1000001_GeV));
    CHECK(model.isValid(Code::Proton, Code::Proton, 999999_GeV));

    //  hydrogen target == proton target == neutron target
    FourMomentum const aP4(100_GeV, {cs, 99_GeV, 0_GeV, 0_GeV});
    FourMomentum const bP4(1_GeV, {cs, 0_GeV, 0_GeV, 0_GeV});
    auto const [xs_prod_pp, xs_ela_pp] =
        model.getCrossSectionInelEla(Code::Proton, Code::Proton, aP4, bP4);
    auto const [xs_prod_pn, xs_ela_pn] =
        model.getCrossSectionInelEla(Code::Proton, Code::Neutron, aP4, bP4);
    auto const [xs_prod_pHydrogen, xs_ela_pHydrogen] =
        model.getCrossSectionInelEla(Code::Proton, Code::Hydrogen, aP4, bP4);
    CHECK(xs_prod_pp == xs_prod_pHydrogen);
    CHECK(xs_prod_pp == xs_prod_pn);
    CHECK(xs_ela_pp == xs_ela_pHydrogen);
    CHECK(xs_ela_pn == xs_ela_pHydrogen);

    // invalids
    auto const xs_prod_0 = model.getCrossSection(Code::Electron, Code::Proton, aP4, bP4);
    CHECK(xs_prod_0 / 1_mb == Approx(0));
    CHECK_THROWS(model.doInteraction(view, Code::Electron, Code::Proton, aP4, bP4));

    CHECK_THROWS(convertFromSibyll(corsika::sibyll::SibyllCode::Unknown));
  }

  SECTION("InteractionInterface - low energy") {

    const HEPEnergyType P0 = 60_GeV;
    MomentumVector const plab = MomentumVector(cs, {P0, 0_eV, 0_eV});
    // also print particles after sibyll was called
    corsika::sibyll::HadronInteractionModel model(stable_particles);
    model.setVerbose(true);
    HEPEnergyType const Elab = sqrt(static_pow<2>(P0) + static_pow<2>(Proton::mass));
    FourMomentum const projectileP4(Elab, plab);
    FourMomentum const nucleusP4(Oxygen::mass, MomentumVector(cs, {0_eV, 0_eV, 0_eV}));
    view.clear();
    model.doInteraction(view, Code::Proton, Code::Oxygen, projectileP4, nucleusP4);
    auto const pSum = sumMomentum(view, cs);

    /*
      Interactions between hadrons (h) and nuclei (A) in Sibyll are treated in the
      hadron-nucleon center-of-mass frame (hnCoM). The incoming hadron (h) and
      nucleon (N) are assumed massless, such that the energy and momentum in the hnCoM are
      : E_i_cm = 0.5 * SQS and P_i_cm = +- 0.5 * SQS  where i is either the projectile
      hadron or the target nucleon and SQS is the hadron-nucleon center-of-mass energy.

      The true energies and momenta, accounting for the hadron masses, are: E_i = ( S +
      m_i**2 - m_j**2 ) / (2 * SQS) and Pcm = +-
      sqrt( (S-(m_j+m_i)**2) * (s-(m_j-m_i)**2) ) / (2*SQS) where m_i is the projectiles
      mass and m_j is the target particles mass. In terms of lab. frame variables Pcm =
      m_j * Plab_i / SQS, where Plab_i is the momentum of the projectile (i) in the lab.
      and m_j is the mass of the target, i.e. the particle at rest (usually a nucleon).

      Any hadron-nucleus event can contain several nucleon interactions. In case of Nw
      (number of wounded nucleons) nucleons interacting in the hadron-nucleus interaction,
      the total energy and momentum in the hadron(i)-nucleon(N) center-of-mass frame are:
      momentum: p_projectile + p_nucleon_1 + p_nucleon_2 + .... p_nucleon_Nw = -(Nw-1) *
      Pcm with center-of-mass momentum Pcm = p_projectile = - p_nucleon_i. For the energy:
      E_projectile + E_nucleon_1 + ... E_nucleon_Nw = E_projectile + Nw * E_nucleon.

      Using the above definitions of center-of-mass energies and momenta this leads to the
      total energy: E_tot = SQS/2 * (1+Nw) + (m_N**2-m_i**2)/(2*SQS) * (Nw-1) and P_tot
      = -m_N * Plab_i / SQS * (Nw-1).

      A Lorentz transformation of these quantities to the lab. frame recovers Plab_i for
      the total momentum, so momentum is exactly conserved, and Elab_i + Nw * m_N for the
      total energy. Not surprisingly, the total energy differs from the total energy
      before the collision by the mass of the additional nucleons (Nw-1)*m_N. In relative
      terms the additional energy is entirely negligible and as it is not kinetic energy
      there is zero influence on the shower development.

      Due to the ommission of the hadron masses in Sibyll, the total energy and momentum
      in the center-of-mass system after the collision are just: E_tot = SQS/2 * (1+Nw)
      and P_tot = SQS/2 * (1-Nw). After the Lorentztransformation the total momentum in
      the lab. thus differs from the initial value by (1-Nw)/2 * ( m_N + m_i**2 / (2 *
      Plab_i) ) and momentum is NOT conserved. Note however that the second term quickly
      vanishes as the lab. momentum of the projectile increases. The first term is fixed
      as it depends only on the number of additional nucleons, in relative terms it is
      always small at high energies.

      For this reason the numerical precision in these tests is limited to 5% to still
      pass at low energies and no absolute check is implemented, e.g.

          CHECK(pSum.getComponents(cs).getX() / P0 == Approx(1).margin(0.05));
          CHECK((pSum - plab).norm()/1_GeV == Approx(0).margin(plab.norm() * 0.05/1_GeV));

      /FR'2020

      See also:

      Issue 272 / MR 204
      https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/merge_requests/204

    */

    CHECK(pSum.getComponents(cs).getX() / P0 == Approx(1).margin(0.05));
    CHECK(pSum.getComponents(cs).getY() / 1_GeV == Approx(0).margin(1e-3));
    CHECK(pSum.getComponents(cs).getZ() / 1_GeV == Approx(0).margin(1e-3));

    CHECK((pSum - plab).getNorm() / 1_GeV ==
          Approx(0).margin(plab.getNorm() * 0.05 / 1_GeV));
    CHECK(pSum.getNorm() / P0 == Approx(1).margin(0.05));
    [[maybe_unused]] CrossSectionType const cx =
        model.getCrossSection(Code::Proton, Code::Oxygen, projectileP4, nucleusP4);
    CHECK(cx / 1_mb == Approx(300).margin(1));
    // CHECK(view.getEntries() == 9); //! \todo: this was 20 before refactory-2020: check
    //                                           "also sibyll not stable wrt. to compiler
    //                                           changes"
  }

  SECTION("NuclearInteractionInterface") {

    HEPMomentumType const P0 = 50_TeV;
    MomentumVector const plab = MomentumVector(cs, {P0, 0_eV, 0_eV});
    corsika::sibyll::HadronInteractionModel hmodel(stable_particles);

    NuclearInteractionModel nuclearModel(hmodel, nuclearcomp);

    CHECK(nuclearModel.isValid(Code::Helium, Code::Oxygen, 100_GeV));
    CHECK_FALSE(nuclearModel.isValid(Code::PiPlus, Code::Oxygen, 100_GeV));
    CHECK_FALSE(nuclearModel.isValid(Code::Electron, Code::Oxygen, 100_GeV));

    Code const pid = Code::Oxygen;
    HEPEnergyType const Elab = sqrt(static_pow<2>(P0) + static_pow<2>(get_mass(pid)));
    FourMomentum const P4(Elab, plab);
    FourMomentum const targetP4(get_mass(Code::Oxygen),
                                MomentumVector(cs, {0_eV, 0_eV, 0_eV}));
    nuclearModel.doInteraction(view, pid, Code::Oxygen, P4, targetP4);
    CrossSectionType const cx =
        nuclearModel.getCrossSection(pid, Code::Oxygen, P4, targetP4);
    CHECK(cx > 0_mb);           // this is not physics validation
    CHECK(view.getSize() != 0); // this is not physics validation

    // invalid to underlying model
    FourMomentum P4mu(
        100_GeV,
        {cs, {sqrt(static_pow<2>(100_GeV) - static_pow<2>(MuPlus::mass)), 0_eV, 0_eV}});
    CrossSectionType const cx0 =
        nuclearModel.getCrossSection(Code::MuPlus, Code::Oxygen, P4mu, targetP4);
    CHECK(cx0 / 1_mb == Approx(0));

    CHECK_THROWS(
        nuclearModel.doInteraction(view, Code::MuPlus, Code::Oxygen, P4mu, targetP4));
  }

  SECTION("CombinedInterface") {
    corsika::sibyll::InteractionModel combinedModel{nuclearcomp, stable_particles};
    corsika::sibyll::HadronInteractionModel const& hmodel =
        combinedModel.getHadronInteractionModel();
    auto const& nuclearModel = combinedModel.getNuclearInteractionModel();

    FourMomentum pP{
        1_TeV, {cs, {calculate_momentum(1_TeV, get_mass(Code::Proton)), 0_eV, 0_eV}}};
    FourMomentum pT{get_mass(Code::Oxygen), {cs, {0_eV, 0_eV, 0_eV}}};

    // pion projectiles go to hadron model
    CHECK(combinedModel.getCrossSection(Code::PiPlus, Code::Oxygen, pP, pT) ==
          hmodel.getCrossSection(Code::PiPlus, Code::Oxygen, pP, pT));
    CHECK_FALSE(combinedModel.getCrossSection(Code::PiPlus, Code::Oxygen, pP, pT) ==
                nuclearModel.getCrossSection(Code::PiPlus, Code::Oxygen, pP, pT));

    // nuclear projectiles go to nuclear model
    CHECK(combinedModel.getCrossSection(Code::Helium, Code::Oxygen, pP, pT) ==
          nuclearModel.getCrossSection(Code::Helium, Code::Oxygen, pP, pT));
    CHECK_FALSE(combinedModel.getCrossSection(Code::Helium, Code::Oxygen, pP, pT) ==
                hmodel.getCrossSection(Code::Helium, Code::Oxygen, pP, pT));
  }
}

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>

#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>

TEST_CASE("SibyllDecayInterface", "modules") {
  logging::set_level(logging::level::info);

  auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Oxygen);
  auto const& cs = *csPtr;
  { [[maybe_unused]] auto const& env_dummy = env; }

  RNGManager<>::getInstance().registerRandomStream("sibyll");

  SECTION("DecayInterface") {

    auto [stackPtr, viewPtr] = setup::testing::setup_stack(
        Code::Lambda, 10_GeV, (DummyEnvironment::BaseNodeType* const)nodePtr, cs);
    test::StackView& view = *viewPtr;
    auto& stack = *stackPtr;
    auto particle = stack.first();

    Decay model;
    model.printDecayConfig();
    [[maybe_unused]] TimeType const time = model.getLifetime(particle);
    auto const gamma = particle.getEnergy() / particle.getMass();
    CHECK(time == get_lifetime(Code::Lambda) * gamma);
    model.doDecay(view);
    // run checks
    // not physics validation, just check doDecay finished with something
    CHECK(stack.getEntries() > 1);
  }

  SECTION("DecayInterface - decay not handled") {
    // sibyll does not know the higgs for example
    auto [stackPtr, viewPtr] = setup::testing::setup_stack(
        Code::H0, 10_GeV, (DummyEnvironment::BaseNodeType* const)nodePtr, cs);
    test::StackView& view = *viewPtr;
    auto& stack = *stackPtr;
    auto particle = stack.first();

    Decay model;

    CHECK(model.getLifetime(particle) == std::numeric_limits<double>::infinity() * 1_s);
    CHECK_THROWS(model.doDecay(view));
  }

  SECTION("DecayConfiguration") {

    Decay model({Code::PiPlus, Code::PiMinus});
    model.printDecayConfig();
    CHECK(model.isDecayHandled(Code::PiPlus));
    CHECK(model.isDecayHandled(Code::PiMinus));
    CHECK_FALSE(model.isDecayHandled(Code::KPlus));

    std::vector<Code> const particleTestList = {Code::PiPlus, Code::PiMinus, Code::KPlus,
                                                Code::LambdaBar, Code::D0Bar};

    // setup decays
    model.setHandleDecay(particleTestList);
    for (auto& pCode : particleTestList) CHECK(model.isDecayHandled(pCode));

    // set decay individually
    model.setHandleDecay(Code::KMinus);
    // fail
    CHECK_THROWS(model.setHandleDecay(Code::H0));

    // possible decays
    CHECK_FALSE(model.canHandleDecay(Code::H0));
    CHECK_FALSE(model.canHandleDecay(Code::Proton));
    CHECK_FALSE(model.canHandleDecay(Code::Electron));
    CHECK(model.canHandleDecay(Code::PiPlus));
    CHECK(model.canHandleDecay(Code::MuPlus));
  }
}
