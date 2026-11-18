/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/Environment.hpp>

namespace corsika {
  namespace media {

    inline Universe::Universe(CoordinateSystem const& pCS)
        : corsika::Sphere(pCS, Point{0 * meter, 0 * meter, 0 * meter},
                          meter * std::numeric_limits<double>::infinity()) {}

    inline bool Universe::contains(corsika::Point const&) const { return true; }

  } // namespace media
} // namespace corsika
