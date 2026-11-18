/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IGeometryTransformer.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/units/quantity.hpp>

namespace corsika {
  namespace media {
    class SphericalTransformer : public IGeometryTransformer {
    private:
      const LengthType radius_; // Radius of the Earth

    public:
      SphericalTransformer(const CoordinateSystem& cs, LengthType radius);
      ~SphericalTransformer() = default;

      LengthType getEffectiveHeight(Point const& position) const;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/geometry_transformer/spherical.inl>