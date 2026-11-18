/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/PhysicalGeometry.hpp>
#include <corsika/framework/core/PhysicalConstants.hpp>

/**
 * \file EnergyMomentumOperations.hpp
 *
 * Relativistic energy momentum calculations.
 */

namespace corsika {

  namespace detail {
    using HEPEnergyTypeSqr = decltype(1_GeV * 1_GeV);
  }

  /**
   * \f[m^2=E_{tot}^2-p^2\f]
   *
   * @param E total energy.
   * @param p particle momentum.
   * @return HEPEnergyType-squared
   */
  auto constexpr calculate_mass_sqr(HEPEnergyType const E, HEPMomentumType const p) {
    return (E - p) * (E + p);
  }

  /**
   * \f[m=sqrt(E_{tot}^2-p^2) \f]
   *
   * @param E total energy.
   * @param p particle momentum.
   * @return HEPEnergyType
   */
  HEPEnergyType constexpr calculate_mass(HEPEnergyType const E, HEPMomentumType const p) {
    return sqrt(calculate_mass_sqr(E, p));
  }

  /**
   * \f[p^2=E_{tot}^2-m^2\f]
   *
   * @param E total energy.
   * @param m particle mass.
   * @return HEPEnergyType-squared
   */
  auto constexpr calculate_momentum_sqr(HEPEnergyType const E, HEPMassType const m) {
    return (E - m) * (E + m);
  }

  /**
   * \f[p=sqrt(E_{tot}^2-m^2) \f]
   *
   * @param E total energy.
   * @param m particle mass.
   * @return HEPEnergyType
   */
  HEPEnergyType constexpr calculate_momentum(HEPEnergyType const E, HEPMassType const m) {
    return sqrt(calculate_momentum_sqr(E, m));
  }

  /**
   * \f[E_{tot}^2=p^2+m^2\f]
   *
   * @param p momentum.
   * @param m particle mass.
   * @return HEPEnergyType-squared
   */
  auto constexpr calculate_total_energy_sqr(HEPMomentumType const p,
                                            HEPMassType const m) {
    return p * p + m * m;
  }

  /**
   * \f[E_{tot}=sqrt(p^2+m^2)\f]
   *
   * @param p momentum.
   * @param m particle mass.
   * @return HEPEnergyType
   */
  HEPEnergyType constexpr calculate_total_energy(HEPMomentumType const p,
                                                 HEPMassType const m) {
    return sqrt(calculate_total_energy_sqr(p, m));
  }

  /**
   * \f[E_{kin}=sqrt(p^2+m^2) - m \f]
   *
   * @param p momentum.
   * @param m particle mass.
   * @return HEPEnergyType
   */
  HEPEnergyType constexpr calculate_kinetic_energy(HEPMomentumType const p,
                                                   HEPMassType const m) {
    return calculate_total_energy(p, m) - m;
  }

  /**
   * \f[E_{lab}=(\sqrt{s}^2 - m_{proj}^2 - m_{targ}^2) / (2 m_{targ}) \f]
   *
   * @param p momentum.
   * @param m particle mass.
   * @return HEPEnergyType
   */
  HEPEnergyType constexpr calculate_lab_energy(detail::HEPEnergyTypeSqr sqrtS_sqr,
                                               HEPMassType const m_proj,
                                               HEPMassType const m_targ) {
    return (sqrtS_sqr - static_pow<2>(m_proj) - static_pow<2>(m_targ)) / (2 * m_targ);
  }

  /**
   * \f[E_{com}=sqrt{2 * m_{proj} * m_{targ} * E_{lab} + m_{proj}^2 + m_{targ}^2} \f]
   *
   * @param E lab. energy.
   * @param m particle mass.
   * @return HEPEnergyType
   */
  HEPEnergyType constexpr calculate_com_energy(HEPEnergyType Elab,
                                               HEPMassType const m_proj,
                                               HEPMassType const m_targ) {
    return sqrt(2 * Elab * m_targ + static_pow<2>(m_proj) + static_pow<2>(m_targ));
  }

} // namespace corsika
