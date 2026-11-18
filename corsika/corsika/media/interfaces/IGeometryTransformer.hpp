/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    /**
     * Interface for geometry transformers.
     *
     * This interface defines the methods required for geometry transformers
     * used in media models.
     */

    class IGeometryTransformer {
    protected:
      const CoordinateSystem& cs_;

    public:
      IGeometryTransformer(const CoordinateSystem& cs);

      virtual ~IGeometryTransformer() = default;

      virtual LengthType getEffectiveHeight(const Point& position) const = 0;
    };

    using IGeometryTransformerPtr = std::shared_ptr<IGeometryTransformer>;

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/interfaces/IGeometryTransformer.inl>
