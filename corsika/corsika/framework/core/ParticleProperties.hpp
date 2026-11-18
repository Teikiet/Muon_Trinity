/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 *   @file ParticleProperties.hpp
 *
 * Interface to particle properties.
 */

#pragma once

#include <array>
#include <cstdint>
#include <cmath>
#include <iosfwd>
#include <string_view>
#include <type_traits>
#include <unordered_map>

#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  /**
   * @defgroup Particles Particle Properties
   *
   * The properties of all particles are saved in static and flat
   * arrays. There is a enum corsika::Code to identify each
   * particle, and each individual particle has its own static class,
   * which can be used to retrieve its physical properties.
   *
   * The properties of all elementary particles are accessible here. The data
   * are taken from the Pythia ParticleData.xml file.
   *
   * Particle data can be accessed via global function in namespace corsika, or via
   * static classes for each particle type. These classes all have the interface (example
   * for the class corsika::Electron):
   *
   * @code{.cpp}
   *   static constexpr Code code{Code::Electron};
   *   static constexpr Code anti_code{Code::Positron};
   *   static constexpr HEPMassType mass{corsika::get_mass(code)};
   *   static constexpr ElectricChargeType charge{corsika::get_charge(code)};
   *   static constexpr int charge_number{corsika::get_charge_number(code)};
   *   static constexpr std::string_view name{corsika::get_name(code)};
   *   static constexpr bool is_nucleus{corsika::is_nucleus(code)};
   * @endcode
   *
   * The names, relations and properties of all particles known to CORSIKA 8 are listed
   * below.
   *
   * **Note** on energy threshold on particle production as well as particle propagation.
   * The functions:
   * @code {.cpp}
   * HEPEnergyType constexpr get_energy_production_threshold(Code const);
   * void constexpr set_energy_production_threshold(Code const, HEPEnergyType const);
   * @endcode
   * can be used to tune the transition where explicit production of new particles, e.g.
   * in Bremsstrahlung, is simulated versus a continuous handling of low-energy particles
   * as generic energy losses. The default value for all particle types is 1 MeV.
   *
   * Furthermore, the functions:
   * @code {.cpp}
   * HEPEnergyType constexpr get_kinetic_energy_propagation_threshold(Code const);
   * void constexpr set_kinetic_energy_propagation_threshold(Code const, HEPEnergyType
   *                                                         const);
   * @endcode
   * are used to discard low energy particle during tracking. The default value for all
   * particle types is 1 GeV.
   *
   * @addtogroup Particles
   * @{
   */

  /**
   * @enum Code
   *
   * The Code enum is the actual place to define CORSIKA 8 particle codes.
   */
  enum class Code : std::int32_t;

  /**
   * @enum PDGCode
   *
   * Specifically for PDG ids.
   */
  enum class PDGCode : std::int32_t;
  /**
   * Internal integer type for enum Code.
   */
  typedef std::underlying_type<Code>::type CodeIntType;

  /**
   * Internal integer type for enum PDGCode.
   */
  typedef std::underlying_type<PDGCode>::type PDGCodeIntType;
} // namespace corsika

// data arrays, etc., as generated automatically
#include <corsika/framework/core/GeneratedParticleProperties.inc>

namespace corsika {

  // forward declarations to be used in GeneratedParticleProperties
  struct full_name {}; //!< tag class for get_name()

  int16_t constexpr get_charge_number(Code const);     //!< electric charge in units of e
  ElectricChargeType constexpr get_charge(Code const); //!< electric charge
  HEPMassType constexpr get_mass(Code const);          //!< mass

  /**
   * Get the kinetic energy propagation threshold.
   *
   * Particles are tracked only above the kinetic energy propagation threshold. Below
   * this, they are discarded and removed. Sensible default values must be configured for
   * a simulation.
   */
  HEPEnergyType get_kinetic_energy_propagation_threshold(Code const);

  /**
   * Set the kinetic energy propagation threshold object.
   */
  void set_kinetic_energy_propagation_threshold(Code const, HEPEnergyType const);

  /**
   * Get the particle production energy threshold.
   *
   * The (total) energy below which a particle is only handled stoachastically (no
   * production below this energy). This is for example important for stochastic discrete
   * Bremsstrahlung versus low-energy Bremsstrahlung as part of continuous energy losses.
   */
  HEPEnergyType get_energy_production_threshold(Code const); //!<

  /**
   * Set the particle production energy threshold in total energies.
   */
  void set_energy_production_threshold(Code const, HEPEnergyType const);

  //! Particle code according to PDG, "Monte Carlo Particle Numbering Scheme"
  PDGCode constexpr get_PDG(Code const);
  PDGCode constexpr get_PDG(unsigned int const A, unsigned int const Z);
  std::string_view constexpr get_name(Code const); //!< name of the particle as string
  std::string get_name(Code,
                       full_name); //!< get name of particle, including (A,Z) for nuclei
  TimeType constexpr get_lifetime(Code const); //!< lifetime

  bool constexpr is_hadron(Code const); //!< true if particle is hadron
  bool constexpr is_kaon(Code const);   //!< true if particle is strange meson
  bool constexpr is_em(Code const); //!< true if particle is electron, positron or photon
  bool constexpr is_muon(Code const);     //!< true if particle is mu+ or mu-
  bool constexpr is_neutrino(Code const); //!< true if particle is (anti-) neutrino
  bool constexpr is_charged(Code const);  //!< true if particle is charged

  /**
   * @brief Creates the Code for a nucleus of type 10LZZZAAAI.
   *
   * @return internal nucleus Code
   */
  Code constexpr get_nucleus_code(size_t const A, size_t const Z);

  /**
   * Checks if Code corresponds to a nucleus.
   *
   * @return true if nucleus.
   * @return false  if not nucleus.
   */
  bool constexpr is_nucleus(Code const);

  /**
   * Get the mass number A for nucleus.
   *
   * @return int size of nucleus.
   */
  size_t constexpr get_nucleus_A(
      Code const); //!< returns A for hard-coded nucleus, otherwise 0

  /**
   * Get the charge number Z for nucleus.
   *
   * @return int charge of nucleus.
   */
  size_t constexpr get_nucleus_Z(
      Code const); //!< returns Z for hard-coded nucleus, otherwise 0

  /**
   * @brief Calculates the mass of nucleus.
   *
   * @return HEPMassType the mass of (A,Z) nucleus, disregarding binding energy.
   */
  HEPMassType constexpr get_nucleus_mass(Code const code);

  /**
   * @brief Calculates the mass of nucleus.
   *
   * @return HEPMassType the mass of (A,Z) nucleus, disregarding binding energy.
   */
  HEPMassType constexpr get_nucleus_mass(unsigned int const A, unsigned int const Z);

  /**
   * @brief Get the nucleus name.
   *
   * @param code
   * @return std::string_view
   */
  inline std::string get_nucleus_name(Code const code);

  /**
   * @brief convert PDG code to CORSIKA 8 internal code.
   *
   * @return Code internal code.
   */
  Code convert_from_PDG(PDGCode const);

  /**
   * @brief Returns list of all non-nuclei particles.
   *
   * @return std::initializer_list<Code> constexpr
   */
  std::initializer_list<Code> constexpr get_all_particles();

  /**
   * @brief Code output operator.
   *
   * The output stream operator for human-readable particle codes.
   *
   * @return std::ostream&
   */
  std::ostream& operator<<(std::ostream&, corsika::Code);

  /** @}*/

} // namespace corsika

#include <corsika/detail/framework/core/ParticleProperties.inl>

// constants in namespaces-like static classes, generated automatically
#include <corsika/framework/core/GeneratedParticleClasses.inc>
