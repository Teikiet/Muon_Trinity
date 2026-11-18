/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <cmath>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>

namespace corsika {

  inline HEPEnergyType get_kinetic_energy_propagation_threshold(Code const code) {
    if (is_nucleus(code)) return particle::detail::threshold_nuclei;
    return particle::detail::propagation_thresholds[static_cast<CodeIntType>(code)];
  }

  inline void set_kinetic_energy_propagation_threshold(Code const code,
                                                       HEPEnergyType const val) {
    if (is_nucleus(code))
      particle::detail::threshold_nuclei = val;
    else
      particle::detail::propagation_thresholds[static_cast<CodeIntType>(code)] = val;
  }

  inline HEPMassType constexpr get_mass(Code const code) {
    if (is_nucleus(code)) { return get_nucleus_mass(code); }
    return particle::detail::masses[static_cast<CodeIntType>(code)];
  }

  inline HEPEnergyType get_energy_production_threshold(Code const p) {
    return particle::detail::production_thresholds[static_cast<CodeIntType>(p)];
  }

  inline void set_energy_production_threshold(Code const p, HEPEnergyType const val) {
    particle::detail::production_thresholds[static_cast<CodeIntType>(p)] = val;
  }

  inline bool constexpr is_charged(Code const c) { return get_charge_number(c) != 0; }

  inline bool constexpr is_nucleus(Code const code) { return code >= Code::Nucleus; }

  inline Code constexpr get_nucleus_code(size_t const A,
                                         size_t const Z) { // 10LZZZAAAI
    if (Z > A) { throw std::runtime_error("Z cannot be larger than A in nucleus."); }
    return static_cast<Code>(static_cast<CodeIntType>(Code::Nucleus) + Z * 10000 +
                             A * 10);
  }

  inline size_t constexpr get_nucleus_Z(Code const code) {
    return (static_cast<CodeIntType>(code) % static_cast<CodeIntType>(Code::Nucleus)) /
           10000;
  }

  inline size_t constexpr get_nucleus_A(Code const code) {
    return (static_cast<CodeIntType>(code) % 10000) / 10;
  }

  inline PDGCode constexpr get_PDG(Code const code) {
    if (code < Code::Nucleus) {
      return particle::detail::pdg_codes[static_cast<CodeIntType>(code)];
    }
    size_t const Z = get_nucleus_Z(code);
    size_t const A = get_nucleus_A(code);
    return static_cast<PDGCode>(static_cast<CodeIntType>(Code::Nucleus) + Z * 10000 +
                                A * 10); // 10LZZZAAAI
  }

  inline PDGCode constexpr get_PDG(unsigned int const A, unsigned int const Z) {
    return PDGCode(1000000000 + Z * 10000 + A * 10);
  }

  inline int16_t constexpr get_charge_number(Code const code) {
    if (is_nucleus(code)) return get_nucleus_Z(code);
    return particle::detail::electric_charges[static_cast<CodeIntType>(code)];
  }

  inline ElectricChargeType constexpr get_charge(Code const code) {
    if (code == Code::Nucleus)
      throw std::runtime_error("charge of particle::Nucleus undefined");
    return get_charge_number(code) * constants::e;
  }

  inline std::string_view constexpr get_name(Code const code) {
    if (is_nucleus(code)) { return "nucleus"; }
    return particle::detail::names[static_cast<CodeIntType>(code)];
  }

  inline std::string get_name(Code code, full_name) {
    if (is_nucleus(code)) {
      return fmt::format("nucleus ({},{})", get_nucleus_A(code), get_nucleus_Z(code));
    }

    return std::string{get_name(code)};
  }

  inline TimeType constexpr get_lifetime(Code const p) {
    return particle::detail::lifetime[static_cast<CodeIntType>(p)] * second;
  }

  inline bool constexpr is_hadron(Code const code) {
    if (is_nucleus(code)) return true;
    return particle::detail::isHadron[static_cast<CodeIntType>(code)];
  }

  inline bool constexpr is_kaon(Code const c) {
    return c == Code::KPlus || c == Code::KMinus || c == Code::K0Long ||
           c == Code::K0Short;
  }

  inline bool constexpr is_em(Code const c) {
    return c == Code::Electron || c == Code::Positron || c == Code::Photon;
  }

  inline bool constexpr is_muon(Code const c) {
    return c == Code::MuPlus || c == Code::MuMinus;
  }

  inline bool constexpr is_neutrino(Code const c) {
    return c == Code::NuE || c == Code::NuMu || c == Code::NuTau || c == Code::NuEBar ||
           c == Code::NuMuBar || c == Code::NuTauBar;
  }

  inline std::ostream& operator<<(std::ostream& stream, corsika::Code const code) {
    return stream << get_name(code);
  }

  inline Code convert_from_PDG(PDGCode const p) {
    static_assert(particle::detail::conversionArray.size() % 2 == 1);
    // this will fail, for the strange case where the maxPDG is negative...
    int constexpr maxPDG{(particle::detail::conversionArray.size() - 1) >> 1};
    auto const k = static_cast<PDGCodeIntType>(p);
    if (std::abs(k) <= maxPDG) {
      return particle::detail::conversionArray[k + maxPDG];
    } else {
      if (1000000000 <= k && k <= 1009999990) { // nucleus (no L or I)
        int const Z = (k - 1000000000) / 10000;
        int const A = (k - 1000000000 - 10000 * Z) / 10;
        return get_nucleus_code(A, Z);
      }
      return particle::detail::conversionMap.at(p);
    }
  }

  inline HEPMassType constexpr get_nucleus_mass(Code const code) {
    unsigned int const A = get_nucleus_A(code);
    unsigned int const Z = get_nucleus_Z(code);
    return get_nucleus_mass(A, Z);
  }

  inline HEPMassType constexpr get_nucleus_mass(unsigned int const A,
                                                unsigned int const Z) {
    return get_mass(Code::Proton) * Z + (A - Z) * get_mass(Code::Neutron);
  }

  inline std::string get_nucleus_name(Code const code) {
    size_t const A = get_nucleus_A(code);
    size_t const Z = get_nucleus_Z(code);
    return fmt::format("Nucleus_A{}_Z{}", A, Z);
  }

  inline std::initializer_list<Code> constexpr get_all_particles() {
    return particle::detail::all_particles;
  }

} // namespace corsika
