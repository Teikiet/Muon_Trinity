/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <functional>

/** Declare the set_rng_function() function inside namespace NAME and corresponding C
 * function set_NAME_rng_function(). To be used inside the public header of a module
 * exposing its interface.
 */
#define DECLARE_RNG(NAME)                                             \
  namespace NAME {                                                    \
    void set_rng_function(std::function<void(double*, std::size_t)>); \
  }                                                                   \
  extern "C" void set_##NAME##_rng_function(void (*rng_function)(double*, std::size_t));
