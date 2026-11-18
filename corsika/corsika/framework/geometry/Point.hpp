/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/BaseVector.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/geometry/Vector.hpp>

namespace corsika {

  /**
   * A Point represents a position in space.
   *
   * Points are stored internally with coordinates in the root coordinate system.
   * Unlike Vectors, Points are affected by both rotation and translation when
   * transforming between coordinate systems.
   */
  class Point : public BaseVector<length_d> {

  public:
    /**
     * Construct Point with coordinates already in root CS (direct storage)
     */
    Point(QuantityVector<length_d> const& pQVector)
        : BaseVector<length_d>(pQVector) {}

    /**
     * Construct Point with coordinates in root CS from individual values
     */
    Point(LengthType x, LengthType y, LengthType z)
        : BaseVector<length_d>({x, y, z}) {}

    /**
     * Construct Point from coordinates in arbitrary CS (transforms to root)
     */
    Point(const CoordinateSystem& cs, QuantityVector<length_d> const& pQVector);

    /**
     * Construct Point from individual coordinates in arbitrary CS (transforms to root)
     */
    Point(const CoordinateSystem& cs, LengthType x, LengthType y, LengthType z);

    /**
     * Get coordinates in root CS
     */
    inline QuantityVector<length_d> const& getCoordinates() const;

    /**
     * Get mutable coordinates in root CS
     */
    inline QuantityVector<length_d>& getCoordinates();

    /**
     * Get coordinates transformed to arbitrary CS
     */
    inline QuantityVector<length_d> getCoordinates(const CoordinateSystem& pCS) const;

    /**
     * Get X, Y, Z coordinate in root CS
     */
    inline LengthType getX() const;
    inline LengthType getY() const;
    inline LengthType getZ() const;

    /**
     * Get X, Y, Z coordinate in arbitrary CS
     *
     *
     * Note: If you need multiple coordinates (X, Y, Z),
     * call getCoordinates(cs) once and access components individually:
     *
     *   // Instead of:
     *   auto x = point.getX(cs);
     *   auto y = point.getY(cs);
     *   auto z = point.getZ(cs);
     *
     *   // Do:
     *   auto coords = point.getCoordinates(cs);
     *   auto x = coords.getX();
     *   auto y = coords.getY();
     *   auto z = coords.getZ();
     *
     */
    inline LengthType getX(const CoordinateSystem& pCS) const;
    inline LengthType getY(const CoordinateSystem& pCS) const;
    inline LengthType getZ(const CoordinateSystem& pCS) const;

    inline Point operator+(Vector<length_d> const& pVec) const;
    inline Point operator-(Vector<length_d> const& pVec) const;

    /**
     * Returns the displacement vector between two points
     */
    inline Vector<length_d> operator-(Point const& pB) const;
  };

  /*
   * calculates the distance between two points
   */
  inline LengthType distance(Point const& p1, Point const& p2);

} // namespace corsika

#include <corsika/detail/framework/geometry/Point.inl>