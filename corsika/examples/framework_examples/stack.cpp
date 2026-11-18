/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/stack/VectorStack.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>

#include <cassert>
#include <iomanip>
#include <iostream>

using namespace corsika;
using namespace std;

void fill(VectorStack& s) {
  CoordinateSystem const& rootCS = get_root_CoordinateSystem();
  for (int i = 0; i < 11; ++i) {
    s.addParticle(std::make_tuple(Code::Electron, 1.5_GeV * i - Electron::mass,
                                  DirectionVector(rootCS, {1, 0, 0}),
                                  Point(rootCS, 0_m, 0_m, 0_m), 0_ns));
  }
}

void read(VectorStack& s) {
  assert(s.getEntries() == 11); // stack has 11 particles

  HEPEnergyType total_energy;
  int i = 0;
  for (auto& p : s) {
    total_energy += p.getEnergy();
    // particles are electrons with 1.5 GeV energy times i
    assert(p.getPID() == Code::Electron);
    assert(p.getEnergy() == 1.5_GeV * (i++));
  }
}

int main() {

  logging::set_level(logging::level::warn);

  CORSIKA_LOG_INFO("stack_example");
  VectorStack s;
  fill(s);
  read(s);
  CORSIKA_LOG_INFO("done");
  return 0;
}
