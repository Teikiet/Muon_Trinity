/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/OnShellCheck.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <corsika/setup/SetupStack.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;

TEST_CASE("OnShellCheck", "[processes]") {

  logging::set_level(logging::level::info);

  feenableexcept(FE_INVALID);
  using EnvType = setup::Environment;
  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();

  // setup empty particle stack
  setup::Stack stack;
  stack.clear();
  // two energies
  const HEPEnergyType E = 10_GeV;
  // list of arbitrary particles
  std::array const particleList{Code::PiPlus, Code::PiMinus,  Code::Helium,
                                Code::Photon, Code::Electron, Code::MuPlus};

  std::array const mass_shifts{1.1, 1.001, 1.0, 1.0, 1.01, 1.0};

  SECTION("check particle masses") {

    OnShellCheck check(1.e-2, 0.01, false);

    // add primary particle to stack
    auto particle = stack.addParticle(
        std::make_tuple(Code::Proton, E, MomentumVector(rootCS, {0_GeV, 0_GeV, 0_GeV}),
                        Point(rootCS, 0_m, 0_m, 0_m), 0_ns));
    // view on secondary particles
    setup::StackView view{particle};
    // ref. to primary particle through the secondary view.
    // only this way the secondary view is populated
    auto projectile = view.getProjectile();
    // add secondaries, all with energies above the threshold
    // only cut is by species
    int count = -1;
    for (auto const proType : particleList) {
      count++;
      const auto pz = sqrt((E - get_mass(proType) * mass_shifts[count]) *
                           (E + get_mass(proType) * mass_shifts[count]));
      auto const momentum = MomentumVector(rootCS, {0_GeV, 0_GeV, pz});
      projectile.addSecondary(
          std::make_tuple(proType, E, momentum, Point(rootCS, 0_m, 0_m, 0_m), 0_ns));
    }
    check.doSecondaries(view);
    int i = -1;
    for (auto const& p : view) {
      i++;
      auto const Plab = FourVector(p.getEnergy(), p.getMomentum());
      auto const m_kinetic = Plab.getNorm();
      if (i == 0)
        CHECK(m_kinetic / PiPlus::mass == Approx(1));
      else if (i == 1)
        CHECK_FALSE(m_kinetic / PiMinus::mass == Approx(1));
    }
  }
}
