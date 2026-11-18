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
#include <corsika/framework/geometry/CoordinateSystem.hpp>

namespace corsika {
  /**
   * Describes a sphere in space
   *
   **/
  class Sphere : public IVolume {
  public:
    Sphere(Point const& pCenter, LengthType const pRadius)
        : center_(pCenter)
        , radius_(pRadius) {}
    // Construct Sphere with center in arbitrary CS (transforms to root)
    Sphere(const CoordinateSystem& cs, Point const& pCenter, LengthType const pRadius);
    //! returns true if the Point p is within the sphere
    bool contains(Point const& p) const override;
    Point const& getCenter() const;
    Point& getCenter();
    void setCenter(Point const&);
    LengthType getRadius() const;
    void setRadius(LengthType const);

    /**
     * Get coordinate system (always root CS for Sphere since it has no orientation).
     * This method provides API compatibility with Box for use in templated code.
     */
    const CoordinateSystem& getCoordinateSystem() const {
      return get_root_CoordinateSystem();
    }

    std::string asString() const;

  private:
    Point center_;
    LengthType radius_;
  };
} // namespace corsika
#include <corsika/detail/framework/geometry/Sphere.inl>