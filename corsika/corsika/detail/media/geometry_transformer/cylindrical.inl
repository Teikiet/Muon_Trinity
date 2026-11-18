/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    CylindricalTransformer::CylindricalTransformer(const CoordinateSystem& cs,
                                                   LengthType radius)
        : IGeometryTransformer(cs) {
      radius_ = radius;
    }

    LengthType CylindricalTransformer::getEffectiveHeight(Point const& point) const {

      double x = point.getX(cs_); // Surface Vector
      double y = point.getY(cs_); // Surface Vector
      double z = point.getZ(cs_); // Up down

      double ref = std::sqrt(x * x + y * y);
      Eigen::Vector3d normal = {ref, 0, z + radius_};

      return normal.norm() - radius_;
    }

  } // namespace media
} // namespace corsika
