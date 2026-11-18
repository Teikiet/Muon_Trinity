/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/urqmd/UrQMD.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <SetupTestStack.hpp>
#include <SetupTestEnvironment.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

#include <tuple>
#include <utility>

#include <catch2/catch_all.hpp>

using namespace corsika;
using namespace corsika::urqmd;
using Catch::Approx;

using DummyEnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

template <typename TStackView>
auto sumCharge(TStackView const& view) {
  int totalCharge = 0;
  for (auto const& p : view) { totalCharge += get_charge_number(p.getPID()); }
  return totalCharge;
}

template <typename TStackView>
auto sumMomentum(TStackView const& view, CoordinateSystem const& vCS) {
  MomentumVector sum{vCS, 0_eV, 0_eV, 0_eV};
  for (auto const& p : view) { sum += p.getMomentum(); }
  return sum;
}

TEST_CASE("UrQMD") {

  logging::set_level(logging::level::info);

  SECTION("conversion") {
    CHECK_THROWS(corsika::urqmd::convertFromUrQMD(106, 0));
    CHECK(corsika::urqmd::convertFromUrQMD(101, 0) == Code::Pi0);
    CHECK(corsika::urqmd::convertToUrQMD(Code::PiPlus) ==
          std::make_pair<int, int>(101, 2));
  }

  feenableexcept(FE_INVALID);
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  UrQMD urqmd;

  auto const& rootCS = get_root_CoordinateSystem();

  SECTION("valid") {
    // this is how it is currently done
    CHECK_FALSE(urqmd.isValid(Code::K0, Code::Proton));
    CHECK_FALSE(urqmd.isValid(Code::DPlus, Code::Proton));
    CHECK_FALSE(urqmd.isValid(Code::Electron, Code::Proton));
    CHECK_FALSE(urqmd.isValid(Code::Proton, Code::Electron));
    CHECK_FALSE(urqmd.isValid(Code::Oxygen, Code::Oxygen));
    CHECK_FALSE(urqmd.isValid(Code::PiPlus, Code::Omega));
    CHECK_FALSE(
        urqmd.isValid(Code::PiPlus, Code::Proton)); // Proton is not a valid target....

    CHECK_NOTHROW(urqmd.isValid(Code::Proton, Code::Oxygen));
    CHECK_NOTHROW(urqmd.isValid(Code::PiPlus, Code::Argon));
  }

  SECTION("cross sections") {

    FourMomentum const targetP4{Nitrogen::mass, {{0_eV, 0_eV, 0_eV}}};

    HEPMomentumType const P0 = 100_GeV;
    Code const validProjectileCodes[] = {
        Code::PiPlus,  Code::PiMinus, Code::Proton, Code::AntiProton, Code::AntiNeutron,
        Code::Neutron, Code::KPlus,   Code::KMinus, Code::K0Long};
    // Code::K0, Code::K0Bar  are not valid projectiles (no mass eigenstates)
    CrossSectionType const checkCX[] = {219_mb, 222_mb, 303_mb, 324_mb, 324_mb,
                                        303_mb, 189_mb, 198_mb, 172_mb};

    int i = 0;
    for (auto const code : validProjectileCodes) {
      FourMomentum const projectileP4{
          sqrt(static_pow<2>(get_mass(code)) + static_pow<2>(P0)), {{0_GeV, 0_GeV, P0}}};
      auto const cx = urqmd.getCrossSection(code, Code::Nitrogen, projectileP4, targetP4);
      CORSIKA_LOG_INFO("UrQMD cross seciton for {} is {} mb", code, cx / 1_mb);
      CHECK(cx / 1_mb == Approx(checkCX[i++] / 1_mb).margin(1));
    }

    // invalid
    CHECK_THROWS(urqmd.getTabulatedCrossSection(Code::Proton, Code::Proton, 100_GeV));
  }

  SECTION("pion+ projectile") {
    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Oxygen);
    auto const& cs = *csPtr;
    [[maybe_unused]] auto const& env_dummy = env;      // against warnings
    [[maybe_unused]] auto const& node_dummy = nodePtr; // against warnings

    auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
        Code::PiPlus, 40_GeV, (DummyEnvironment::BaseNodeType* const)nodePtr, cs);

    // must be assigned to variable, cannot be used as rvalue?!
    auto projectile = secViewPtr->getProjectile();
    auto const projectileMomentum = projectile.getMomentum();

    FourMomentum const projectileP4{
        sqrt(static_pow<2>(PiPlus::mass) + static_pow<2>(40_GeV)),
        {40_GeV, 0_GeV, 0_GeV}};
    FourMomentum const targetP4{Oxygen::mass, {0_GeV, 0_GeV, 0_GeV}};
    urqmd.doInteraction(*secViewPtr, Code::PiPlus, Code::Oxygen, projectileP4, targetP4);

    CHECK(sumCharge(*secViewPtr) ==
          get_charge_number(Code::PiPlus) + get_charge_number(Code::Oxygen));

    auto const secMomSum = sumMomentum(*secViewPtr, rootCS);
    CHECK((secMomSum - projectileMomentum).getNorm() / projectileMomentum.getNorm() ==
          Approx(0).margin(1e-2));
  }

  SECTION("K0Long projectile") {
    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Oxygen);
    auto const& cs = *csPtr;
    [[maybe_unused]] auto const& env_dummy = env;      // against warnings
    [[maybe_unused]] auto const& node_dummy = nodePtr; // against warnings

    auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
        Code::K0Long, 40_GeV, (DummyEnvironment::BaseNodeType* const)nodePtr, cs);
    CHECK(stackPtr->getEntries() == 1);
    CHECK(secViewPtr->getEntries() == 0);

    // must be assigned to variable, cannot be used as rvalue?!
    auto projectile = secViewPtr->getProjectile();
    auto const projectileMomentum = projectile.getMomentum();

    FourMomentum const projectileP4{
        sqrt(static_pow<2>(K0Long::mass) + static_pow<2>(40_GeV)),
        {40_GeV, 0_GeV, 0_GeV}};
    FourMomentum const targetP4{Oxygen::mass, {0_GeV, 0_GeV, 0_GeV}};
    urqmd.doInteraction(*secViewPtr, Code::K0Long, Code::Oxygen, projectileP4, targetP4);

    CHECK(sumCharge(*secViewPtr) ==
          get_charge_number(Code::K0Long) + get_charge_number(Code::Oxygen));

    auto const secMomSum = sumMomentum(*secViewPtr, rootCS);
    CHECK((secMomSum - projectileMomentum).getNorm() / projectileMomentum.getNorm() ==
          Approx(0).margin(1e-2));
  }
}
