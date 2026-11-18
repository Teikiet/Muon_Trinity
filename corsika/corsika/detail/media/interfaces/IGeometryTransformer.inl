/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    inline IGeometryTransformer::IGeometryTransformer(const CoordinateSystem& cs)
        : cs_(cs) {}

  } // namespace media
} // namespace corsika
