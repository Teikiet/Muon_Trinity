/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <catch2/catch_all.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>

using namespace corsika;
using Catch::Approx;

TEST_CASE("EnergyMomentumOperations") {

  logging::set_level(logging::level::info);

  HEPMomentumType const p0 = 30_GeV;
  HEPMassType const m0 = 40_GeV;
  HEPEnergyType const e0 = 50_GeV;

  CHECK(calculate_mass_sqr(e0, p0) / (m0 * m0) == Approx(1));
  CHECK(calculate_mass(e0, p0) / m0 == Approx(1));

  CHECK(calculate_momentum_sqr(e0, m0) / (p0 * p0) == Approx(1));
  CHECK(calculate_momentum(e0, m0) / p0 == Approx(1));

  CHECK(calculate_total_energy_sqr(p0, m0) / (e0 * e0) == Approx(1));
  CHECK(calculate_total_energy(p0, m0) / e0 == Approx(1));

  CHECK(calculate_kinetic_energy(p0, m0) / (e0 - m0) == Approx(1));
  CHECK(calculate_lab_energy((100_GeV * 100_GeV), 1_GeV, 100_MeV) / 49999.99_GeV ==
        Approx(1).margin(0.1));
}
