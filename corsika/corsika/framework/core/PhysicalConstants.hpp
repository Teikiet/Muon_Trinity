/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>

/**
 * \file PhysicalConstants.hpp
 *
 * Constants are defined with static units, based on the package
 * (namespace) phys::units, imported in PhysicsUnits.hpp
 *
 * \namespace corsika::constants
 *
 * Physical and mathematical constants with units.
 */

namespace corsika::constants {

  using namespace phys::units;

  // acceleration of free-fall, standard
  constexpr quantity<acceleration_d> g_sub_n{Rep(9.80665L) * meter / square(second)};

  // Avogadro constant
  constexpr quantity<dimensions<0, 0, 0, 0, 0, -1>> N_sub_A{Rep(6.02214199e+23L) / mole};

  // elementary charge
  constexpr quantity<electric_charge_d> e{Rep(1.6021766208e-19L) * coulomb};

  // vacuum permittivity
  constexpr quantity<dimensions<-3, -1, 4, 2>> epsilonZero{Rep(8.8541878128e-12L) *
                                                           farad / meter};

  // electronvolt
  // constexpr quantity<hepenergy_d> eV{e / coulomb * joule};

  // Planck constant
  constexpr quantity<dimensions<2, 1, -1>> h{Rep(6.62606876e-34L) * joule * second};
  constexpr quantity<dimensions<2, 1, -1>> hBar{h / (2 * M_PI)};

  // speed of light in a vacuum
  constexpr quantity<speed_d> c{Rep(299792458L) * meter / second};
  constexpr auto cSquared = c * c;

  // hbar * c
  constexpr quantity<dimensions<1, 0, 0, 0, 0, 0, 0, 1>> hBarC{
      Rep(1.973'269'78e-7L) * electronvolt * meter}; // from RPP 2018

  auto constexpr invGeVsq = 1e-18 / (electronvolt * electronvolt);

  // unified atomic mass unit
  constexpr quantity<mass_d> u{Rep(1.6605402e-27L) * kilogram};

  auto constexpr nucleonMass = 0.5 * (0.93827 + 0.93957) * 1e9 * electronvolt;

  // molar gas constant
  auto constexpr R = Rep(8.314'459'8) * joule / (mole * kelvin);

  /**
   * A namespace containing various Earth radii.
   */
  namespace EarthRadius {
    static constexpr auto Mean{6'371'000 * meter};
    static constexpr auto Geomagnetic_reference{6'371'200 * meter};
    static constexpr auto Equatorial{6'378'137 * meter};
    static constexpr auto Polar{6'356'752 * meter};
    static constexpr auto PolarCurvature{6'399'593 * meter};
  } // namespace EarthRadius

  // etc.

} // namespace corsika::constants
