/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    SphericalTransformer::SphericalTransformer(const CoordinateSystem& cs,
                                               LengthType radius)
        : IGeometryTransformer(cs)
        , radius_(radius) {}

    LengthType SphericalTransformer::getEffectiveHeight(Point const&) const {
      Point const center(cs_, 0_m, 0_m, 0_m);
      return distance(point, center_) - radius_; // If from surface a simple position.z()
                                                 // would work and be much faster;
    }

  } // namespace media
} // namespace corsika