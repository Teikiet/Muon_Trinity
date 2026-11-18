/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/pythia8/Random.hpp>

namespace corsika::pythia8 {

  inline double Random::flat() { return Dist_(RNG_); }

} // namespace corsika::pythia8
