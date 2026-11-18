/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <urqmd.hpp>

#include <rng_impl.hpp>

IMPLEMENT_RNG(urqmd)

namespace urqmd {
  extern "C" double ranf_(int&) { return draw_std_rnd(); }
} // namespace urqmd
