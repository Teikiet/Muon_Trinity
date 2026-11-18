/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/BaseVector.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Vector.hpp>

namespace corsika {

  // Constructors that transform to root CS (full transformation: rotation + translation)
  inline Point::Point(const CoordinateSystem& cs,
                      QuantityVector<length_d> const& pQVector)
      : BaseVector<length_d>(cs.isRoot()
                                 ? pQVector
                                 : QuantityVector<length_d>(cs.getTransformToRoot() *
                                                            pQVector.eigenVector_)) {}

  inline Point::Point(const CoordinateSystem& cs, LengthType x, LengthType y,
                      LengthType z)
      : Point(cs, QuantityVector<length_d>{x, y, z}) {}

  // getCoordinates: direct access in root CS
  inline QuantityVector<length_d> const& Point::getCoordinates() const {
    return BaseVector<length_d>::getQuantityVector();
  }

  inline QuantityVector<length_d>& Point::getCoordinates() {
    return BaseVector<length_d>::getQuantityVector();
  }

  // getCoordinates(cs): transform from root to target CS (full transformation)
  inline QuantityVector<length_d> Point::getCoordinates(
      const CoordinateSystem& pCS) const {
    if (pCS.isRoot()) { return BaseVector<length_d>::getQuantityVector(); }
    // Transform: root -> pCS (rotation + translation for points)
    return QuantityVector<length_d>(
        pCS.getTransformFromRoot() *
        BaseVector<length_d>::getQuantityVector().eigenVector_);
  }

  // getX/Y/Z: direct access to root CS coordinates
  inline LengthType Point::getX() const { return getCoordinates().getX(); }

  inline LengthType Point::getY() const { return getCoordinates().getY(); }

  inline LengthType Point::getZ() const { return getCoordinates().getZ(); }

  // getX/Y/Z(cs): transform individual coordinates to target CS
  inline LengthType Point::getX(const CoordinateSystem& pCS) const {
    return getCoordinates(pCS).getX();
  }

  inline LengthType Point::getY(const CoordinateSystem& pCS) const {
    return getCoordinates(pCS).getY();
  }

  inline LengthType Point::getZ(const CoordinateSystem& pCS) const {
    return getCoordinates(pCS).getZ();
  }

  // Point + Vector: both in root CS, direct addition
  inline Point Point::operator+(Vector<length_d> const& pVec) const {
    return Point(getCoordinates() + pVec.getComponents());
  }

  // Point - Vector: both in root CS, direct subtraction
  inline Point Point::operator-(Vector<length_d> const& pVec) const {
    return Point(getCoordinates() - pVec.getComponents());
  }

  // Point - Point: displacement vector between two points
  inline Vector<length_d> Point::operator-(Point const& pB) const {
    return Vector<length_d>(getCoordinates() - pB.getCoordinates());
  }

  inline std::ostream& operator<<(std::ostream& os, corsika::Point const& p) {
    auto const& qv = p.getCoordinates();
    os << qv << " (root CS)";
    return os;
  }

  inline LengthType distance(Point const& p1, Point const& p2) {
    return (p1 - p2).getNorm();
  }

} // namespace corsika