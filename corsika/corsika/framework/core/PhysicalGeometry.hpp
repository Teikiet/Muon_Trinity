/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Vector.hpp>

/**
 * \file PhysicalUnits.hpp
 *
 * Import and extend the phys::units package. The SI units are also imported into the
 * `\namespace corsika`, since they are used everywhere as integral part of the framework.
 */

namespace corsika {

  /**
   * A 3D vector defined in the root coordinate system with units HEPMomentumType
   **/
  using MomentumVector = Vector<hepmomentum_d>;

  /**
   * A 3D vector defined in the root coordinate system with no units. But, note, this is
   * not automatically normaliyed! It is not a "NormalVector".
   **/
  using DirectionVector = Vector<dimensionless_d>;

  /**
   * A 3D vector defined in the root coordinate system with units "velocity_t".
   *
   **/
  using VelocityVector = Vector<SpeedType::dimension_type>;

  /**
   * A 3D vector defined in the root coordinate system with units "length_t".
   *
   **/
  using LengthVector = Vector<length_d>;

  /**
   * A 3D vector defined in the root coordinate system with units ElectricFieldType
   **/
  typedef Vector<ElectricFieldType::dimension_type> ElectricFieldVector;

  /**
   * A 3D vector defined in the root coordinate system with units VectorPotentialType
   **/
  typedef Vector<VectorPotentialType::dimension_type> VectorPotential;

} // namespace corsika
