/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <cfenv>

extern "C" {
#warning No enabling/disabling of floating point exceptions - platform needs better implementation

inline int feenableexcept(int /*excepts*/) noexcept { return -1; }

inline int fedisableexcept(int /*excepts*/) noexcept { return -1; }
}
