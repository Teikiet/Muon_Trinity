/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>

#include <array>
#include <utility>
#include <string>

namespace corsika::urqmd {

  /**
   * Checks if particle with Code can interaction in UrQMD.
   */
  bool canInteract(Code const);

  /**
   * convert CORSIKA code to UrQMD code tuple.
   *
   * In the current implementation a detour via the PDG code is made.
   */
  std::pair<int, int> convertToUrQMD(Code const);

  /**
   * convert UrQMD code to CORSIKA code.
   */
  Code convertFromUrQMD(int vItyp, int vIso3);

} // namespace corsika::urqmd

#include <corsika/detail/modules/urqmd/ParticleConversion.inl>
