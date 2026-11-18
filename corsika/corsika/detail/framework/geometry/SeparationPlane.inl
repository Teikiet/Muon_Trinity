/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {

  inline SeparationPlane::SeparationPlane(Plane const& plane)
      : plane_(plane) {}

  inline bool SeparationPlane::contains(Point const& p) const {
    return !plane_.isAbove(p);
  }

  inline std::string SeparationPlane::asString() const { return plane_.asString(); }

} // namespace corsika