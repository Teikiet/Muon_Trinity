/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <functional>

namespace conex {
  void set_rng_function(std::function<void(double*, std::size_t)>);
}
