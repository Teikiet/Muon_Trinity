/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  // Constructor with CS explicit (transforms center to root)
  inline Sphere::Sphere(const CoordinateSystem& cs, Point const& pCenter,
                        LengthType const pRadius)
      : center_(Point(cs, pCenter.getCoordinates()))
      , radius_(pRadius) {}

  inline bool Sphere::contains(Point const& p) const {
    // Both center_ and p are in root CS internally
    // Vector subtraction works directly in root CS
    return radius_ * radius_ > (center_ - p).getSquaredNorm();
  }

  // Direct access in root CS
  inline Point const& Sphere::getCenter() const { return center_; }
  inline Point& Sphere::getCenter() { return center_; }

  inline void Sphere::setCenter(Point const& p) { center_ = p; }

  inline LengthType Sphere::getRadius() const { return radius_; }

  inline void Sphere::setRadius(LengthType const r) { radius_ = r; }

  inline std::string Sphere::asString() const {
    std::ostringstream txt;
    txt << "center=" << center_ << " (root CS), radius=" << radius_;
    return txt.str();
  }

} // namespace corsika