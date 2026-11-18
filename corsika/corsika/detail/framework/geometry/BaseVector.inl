/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>

namespace corsika {

  template <typename TDimension>
  inline QuantityVector<TDimension> const& BaseVector<TDimension>::getQuantityVector()
      const {
    return quantityVector_;
  }

  template <typename TDimension>
  inline QuantityVector<TDimension>& BaseVector<TDimension>::getQuantityVector() {
    return quantityVector_;
  }

} // namespace corsika