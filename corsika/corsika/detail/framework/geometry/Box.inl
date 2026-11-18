/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

namespace corsika {

  // Constructor with explicit CS for orientation
  inline Box::Box(const CoordinateSystem& cs, LengthType const x, LengthType const y,
                  LengthType const z)
      : center_(Point(
            cs, {0_m, 0_m, 0_m})) // Center at (0,0,0) in box's CS, transformed to root
      , cs_(&cs)
      , x_(x)
      , y_(y)
      , z_(z) {}

  inline Box::Box(const CoordinateSystem& cs, LengthType const side)
      : center_(Point(cs, {0_m, 0_m, 0_m}))
      , cs_(&cs)
      , x_(side / 2)
      , y_(side / 2)
      , z_(side / 2) {}

  // Convenience constructors using root CS
  inline Box::Box(LengthType const x, LengthType const y, LengthType const z)
      : Box(get_root_CoordinateSystem(), x, y, z) {}

  inline Box::Box(LengthType const side)
      : Box(get_root_CoordinateSystem(), side) {}

  inline bool Box::contains(Point const& p) const {
    // Transform point to box's local CS to check containment
    // Point p is in root CS, we get its coordinates in box's CS
    auto const px = p.getX(*cs_);
    auto const py = p.getY(*cs_);
    auto const pz = p.getZ(*cs_);

    return (abs(px) < x_) && (abs(py) < y_) && (abs(pz) < z_);
  }

  inline std::string Box::asString() const {
    std::ostringstream txt;
    txt << "center=" << center_ << " (root CS)"
        << ", x-axis=" << DirectionVector{*cs_, {1, 0, 0}}
        << ", y-axis: " << DirectionVector{*cs_, {0, 1, 0}}
        << ", z-axis: " << DirectionVector{*cs_, {0, 0, 1}};
    return txt.str();
  }

  template <typename TDim>
  inline void Box::rotate(QuantityVector<TDim> const& axis, double const angle) {
    // Create new rotated CS in registry and update pointer
    // make_rotation returns const reference from registry
    cs_ = &make_rotation(*cs_, axis, angle);

    // Note: center_ remains unchanged as it's stored in root CS
    // The box orientation has changed, but not its position
  }

} // namespace corsika