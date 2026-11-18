/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using Catch::Approx;

TEST_CASE("ParticleProperties", "[Particles]") {

  logging::set_level(logging::level::info);

  SECTION("Types") {
    CHECK(Electron::code == Code::Electron);
    CHECK(Positron::code == Code::Positron);
    CHECK(Proton::code == Code::Proton);
    CHECK(Neutron::code == Code::Neutron);
    CHECK(Photon::code == Code::Photon);
    CHECK(PiPlus::code == Code::PiPlus);
  }

  SECTION("Masses") {
    CHECK(Electron::mass / (511_keV) == Approx(1));
    CHECK(Electron::mass / get_mass(Code::Electron) == 1.);
    CHECK(Photon::mass / (1_eV) == 0.);
    CHECK(Photon::mass == get_mass(Code::Photon));

    CHECK((Proton::mass + Neutron::mass) / constants::nucleonMass == Approx(2));
  }

  SECTION("Charges") {
    CHECK(Electron::charge / constants::e == Approx(-1));
    CHECK(Positron::charge / constants::e == Approx(+1));
    CHECK(get_charge(Positron::anti_code) / constants::e == Approx(-1));
    CHECK(Photon::charge / constants::e == 0.);
    CHECK_FALSE(is_charged(Code::Photon));
    CHECK(is_charged(Code::Iron));
    CHECK(is_charged(Code::PiPlus));
  }

  SECTION("Names") {
    CHECK(Electron::name == "Electron");
    CHECK(get_name(Code::Electron) == "Electron");
    CHECK(PiMinus::name == "PiMinus");
    CHECK(Iron::name == "nucleus");
    CHECK(Photon::name == "Photon");
  }

  SECTION("PDG") {
    CHECK(get_PDG(Code::PiPlus) == PDGCode::PiPlus);
    CHECK(get_PDG(Code::DPlus) == PDGCode::DPlus);
    CHECK(get_PDG(Code::NuMu) == PDGCode::NuMu);
    CHECK(get_PDG(Code::NuE) == PDGCode::NuE);
    CHECK(get_PDG(Code::MuMinus) == PDGCode::MuMinus);
    CHECK(get_PDG(Code::Photon) == PDGCode::Photon);

    CHECK(static_cast<int>(get_PDG(Code::PiPlus)) == 211);
    CHECK(static_cast<int>(get_PDG(Code::DPlus)) == 411);
    CHECK(static_cast<int>(get_PDG(Code::NuMu)) == 14);
    CHECK(static_cast<int>(get_PDG(Code::NuEBar)) == -12);
    CHECK(static_cast<int>(get_PDG(Code::MuMinus)) == 13);
    CHECK(static_cast<int>(get_PDG(Code::Photon)) == 22);
  }

  SECTION("Conversion PDG -> internal") {
    CHECK(convert_from_PDG(PDGCode::KStarMinus) == Code::KStarMinus);
    CHECK(convert_from_PDG(PDGCode::MuPlus) == Code::MuPlus);

    CHECK(static_cast<int>(Code::SigmaC_2455_MinusBar) == 116);
    CHECK(static_cast<int>(get_PDG(Code::SigmaC_2455_MinusBar)) == -4212);
    CHECK(static_cast<int>(convert_from_PDG(PDGCode(-4212))) == 116);

    CHECK(convert_from_PDG(PDGCode::SigmaC_2455_MinusBar) == Code::SigmaC_2455_MinusBar);
  }

  SECTION("Lifetimes") {
    CHECK(get_lifetime(Code::Electron) ==
          std::numeric_limits<double>::infinity() * si::second);
    CHECK(get_lifetime(Code::DPlus) < get_lifetime(Code::Photon));
    CHECK(get_lifetime(Code::RhoPlus) / si::second ==
          (Approx(4.414566727909413e-24).epsilon(1e-3)));
    CHECK(get_lifetime(Code::SigmaMinusBar) / si::second ==
          (Approx(8.018175e-11).epsilon(1e-5)));
    CHECK(get_lifetime(Code::MuPlus) / si::second ==
          (Approx(2.196981e-06).epsilon(1e-5)));
  }

  SECTION("Energy thresholds") {
    //! by default energy thresholds are set to zero
    CHECK(get_kinetic_energy_propagation_threshold(Electron::code) == 1_GeV);
    set_kinetic_energy_propagation_threshold(Electron::code, 10_GeV);
    CHECK_FALSE(get_kinetic_energy_propagation_threshold(Code::Electron) == 1_GeV);
    CHECK(get_kinetic_energy_propagation_threshold(Code::Electron) == 10_GeV);

    //! by default energy thresholds are set to zero
    CHECK(get_energy_production_threshold(Neutron::code) == 1_MeV);
    set_energy_production_threshold(Neutron::code, 1_GeV);
    CHECK_FALSE(get_energy_production_threshold(Code::Neutron) == 1_MeV);
    CHECK(get_energy_production_threshold(Code::Neutron) == 1_GeV);
  }

  SECTION("Particle groups: electromagnetic") {
    CHECK(is_em(Code::Photon));
    CHECK(is_em(Code::Electron));
    CHECK_FALSE(is_em(Code::MuPlus));
    CHECK_FALSE(is_em(Code::NuE));
    CHECK_FALSE(is_em(Code::Proton));
    CHECK_FALSE(is_em(Code::PiPlus));
    CHECK_FALSE(is_em(Code::Oxygen));
  }

  SECTION("Particle groups: hadrons") {
    CHECK_FALSE(is_hadron(Code::Photon));
    CHECK_FALSE(is_hadron(Code::Electron));
    CHECK_FALSE(is_hadron(Code::MuPlus));
    CHECK_FALSE(is_hadron(Code::NuE));
    CHECK(is_hadron(Code::Proton));
    CHECK(is_hadron(Code::PiPlus));
    CHECK(is_hadron(Code::Oxygen));
  }

  SECTION("Particle groups: muons") {
    CHECK_FALSE(is_muon(Code::Photon));
    CHECK_FALSE(is_muon(Code::Electron));
    CHECK(is_muon(Code::MuPlus));
    CHECK(is_muon(Code::MuMinus));
    CHECK_FALSE(is_muon(Code::NuE));
    CHECK_FALSE(is_muon(Code::Proton));
    CHECK_FALSE(is_muon(Code::PiPlus));
    CHECK_FALSE(is_muon(Code::Oxygen));
  }

  SECTION("Particle groups: neutrinos") {
    CHECK_FALSE(is_neutrino(Code::Photon));
    CHECK_FALSE(is_neutrino(Code::Electron));
    CHECK_FALSE(is_neutrino(Code::MuPlus));
    CHECK_FALSE(is_neutrino(Code::Proton));
    CHECK_FALSE(is_neutrino(Code::PiPlus));
    CHECK_FALSE(is_neutrino(Code::Oxygen));

    CHECK(is_neutrino(Code::NuE));
    CHECK(is_neutrino(Code::NuMu));
    CHECK(is_neutrino(Code::NuTau));
    CHECK(is_neutrino(Code::NuEBar));
    CHECK(is_neutrino(Code::NuMuBar));
    CHECK(is_neutrino(Code::NuTauBar));
  }

  SECTION("Nuclei") {
    CHECK_FALSE(is_nucleus(Code::Photon));
    CHECK(is_nucleus(Code::Argon));
    CHECK_FALSE(is_nucleus(Code::Proton));
    CHECK(is_nucleus(Code::Hydrogen));
    CHECK(Argon::is_nucleus);
    CHECK_FALSE(EtaC::is_nucleus);

    CHECK(get_nucleus_A(Code::Hydrogen) == 1);
    CHECK(get_nucleus_A(Code::Tritium) == 3);
    CHECK(Hydrogen::nucleus_Z == 1);
    CHECK(Tritium::nucleus_A == 3);

    CHECK(is_nucleus(get_nucleus_code(1, 1)));
    CHECK(is_nucleus(get_nucleus_code(100, 100)));
    CHECK(get_nucleus_code(208, 82) == Code::Lead);
    CHECK_FALSE(is_nucleus(Code::Electron));
    CHECK(is_nucleus(Code::Lead));
    CHECK(get_nucleus_Z(Code::Lead) == 82);
    CHECK(get_nucleus_A(Code::Lead) == 208);

    // impossible nucleus
    CHECK_THROWS(get_nucleus_code(20, 40));

    // getters
    auto const testId = get_nucleus_code(40, 20);
    CHECK(get_nucleus_A(testId) == 40);
    CHECK(get_nucleus_Z(testId) == 20);
    CHECK(is_nucleus(testId));
    CHECK(get_nucleus_mass(testId) == 20 * Proton::mass + 20 * Neutron::mass);
    CHECK(get_name(testId) == "nucleus");
    CHECK(get_charge(testId) == 20 * constants::e);
  }
}
