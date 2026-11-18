/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <vector>

namespace corsika::history {

  /**
   * This class stores the non-common properties of secondaries in an event. All
   * other (common) properties are available via the event itself or its projectile.
   */
  struct SecondaryParticle {
    HEPEnergyType const energy_;
    Vector<hepmomentum_d> const momentum_;
    Code const pid_;

  public:
    SecondaryParticle(HEPEnergyType energy, Vector<hepmomentum_d> momentum, Code pid)
        : energy_{energy}
        , momentum_{momentum}
        , pid_{pid} {}
  };

} // namespace corsika::history
