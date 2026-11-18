/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Sphere.hpp>
#include <limits>

namespace corsika {
  namespace media {

    struct Universe : public corsika::Sphere {
      Universe(corsika::CoordinateSystem const& pCS);
      bool contains(corsika::Point const&) const override;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/Universe.inl>
