/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Point.hpp>

namespace corsika {

  class IVolume {

  public:
    //! returns true if the Point p is within the volume
    virtual bool contains(Point const& p) const = 0;

    virtual ~IVolume() = default;
  };

} // namespace corsika
