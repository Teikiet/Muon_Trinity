/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <string>

namespace corsika {
  namespace media {

    template <typename T>
    struct NameModel : public T {
      virtual std::string const& GetName() const = 0;
      virtual ~NameModel() = default;
    };

  } // namespace media
} // namespace corsika
