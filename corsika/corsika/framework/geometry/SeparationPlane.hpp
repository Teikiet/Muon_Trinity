/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/IVolume.hpp>

namespace corsika {

  /**
   *  Describes which exists on one side of a plane
   *
   *  The separation line is given by a plane where the unit
   *  vector of the plane points outside of the volume
   **/
  class SeparationPlane : public IVolume {

  public:
    SeparationPlane(Plane const& plane);

    ~SeparationPlane() {}

    Plane getPlane() const { return plane_; }

    //! returns true if the Point p is below the plane
    bool contains(Point const& p) const override;

    std::string asString() const;

  protected:
    Plane const plane_;
  };

} // namespace corsika

#include <corsika/detail/framework/geometry/SeparationPlane.inl>
