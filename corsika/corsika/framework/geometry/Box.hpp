/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/IVolume.hpp>

namespace corsika {

  /**
   * Describes a box in space
   *
   *  The center point and the orientation of the Box is set by
   *  a CoordinateSystem at construction and the sides extend
   *  by x, y, z in both directions.
   **/
  class Box : public IVolume {

  public:
    Box(CoordinateSystem const& cs, LengthType const x, LengthType const y,
        LengthType const z);

    Box(CoordinateSystem const& cs, LengthType const side);

    // Convenience constructors assuming root coordinate system
    Box(LengthType const x, LengthType const y, LengthType const z);
    Box(LengthType const side);

    //! returns true if the Point p is within the box
    bool contains(Point const& p) const override;

    Point const& getCenter() const { return center_; };

    CoordinateSystem const& getCoordinateSystem() const { return *cs_; }

    LengthType const getX() const { return x_; }
    LengthType const getY() const { return y_; }
    LengthType const getZ() const { return z_; }

    std::string asString() const;

    template <typename TDim>
    void rotate(QuantityVector<TDim> const& axis, double const angle);

  protected:
    Point center_;
    CoordinateSystem const* cs_; // local coordinate system with center_ in coordinate (0,
                                 // 0, 0) and user defined orientation
    LengthType x_;
    LengthType y_;
    LengthType z_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/Box.inl>