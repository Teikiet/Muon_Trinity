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

  /**
   * Common base class for Vector and Point.
   *
   * Stores a QuantityVector with coordinates always in root coordinate system.
   * No coordinate system member is stored, reducing memory footprint by 25%.
   * Transformations are performed on-demand when accessing coordinates in
   * different coordinate systems.
   */
  template <typename TDimension>
  class BaseVector {

  public:
    /**
     * Construct from QuantityVector in root CS (direct storage, no transformation)
     */
    BaseVector(QuantityVector<TDimension> const& pQVector)
        : quantityVector_(pQVector) {}

    BaseVector() = delete;
    BaseVector(BaseVector const&) = default;
    BaseVector(BaseVector&& a) = default;
    BaseVector& operator=(BaseVector const&) = default;
    ~BaseVector() = default;

  protected:
    QuantityVector<TDimension> const& getQuantityVector() const;
    QuantityVector<TDimension>& getQuantityVector();
    void setQuantityVector(QuantityVector<TDimension> const& v) { quantityVector_ = v; }

  private:
    QuantityVector<TDimension> quantityVector_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/BaseVector.inl>