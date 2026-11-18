/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <string>

namespace corsika {

  class Plane {

  public:
    Plane(Point const& vCenter, DirectionVector const& vNormal)
        : center_(vCenter)
        , normal_(vNormal.normalized()) {}

    bool isAbove(Point const& vP) const;

    LengthType getDistanceTo(Point const& vP) const;

    Point const& getCenter() const;

    DirectionVector const& getNormal() const;

    std::string asString() const;

  public:
    Point const center_;
    DirectionVector const normal_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/Plane.inl>
