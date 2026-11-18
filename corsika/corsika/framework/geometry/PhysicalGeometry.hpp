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

  typedef Vector<hepmomentum_d> MomentumVector;
  typedef Vector<dimensionless_d> DirectionVector;
  typedef Vector<magnetic_flux_density_d> MagneticFieldVector;
  typedef Vector<SpeedType::dimension_type> VelocityVector;

} // namespace corsika
