/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    inline PlanarTransformer::PlanarTransformer(const CoordinateSystem& cs)
        : IGeometryTransformer(cs) {}

    inline LengthType PlanarTransformer::getEffectiveHeight(Point const& point) const {
      return point.getZ(cs_);
    }

  } // namespace media
} // namespace corsika