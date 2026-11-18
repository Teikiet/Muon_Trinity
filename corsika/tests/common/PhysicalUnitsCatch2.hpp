/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>

#include <catch2/catch_all.hpp>

using Catch::Approx;

namespace corsika::testing {

  /**
   * comparion of DimensionlessType quantities with the Catch2 Approx
   * class.

   * This allows code/checks of this type
   * `CHECK(v.normalize().norm() == Approx(1),margin(0)) `
   *
   **/
  inline bool operator==(DimensionlessType const a, Approx const& b) {
    return a.magnitude() == b;
  }

} // namespace corsika::testing
