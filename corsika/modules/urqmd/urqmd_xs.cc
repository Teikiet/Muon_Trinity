/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

// a little helper to dump UrQMD cross-sections

#include <corsika/process/urqmd/UrQMD.h>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/random/RNGManager.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace corsika;
using namespace corsika::units::si;

int main() {
  random::RNGManager::getInstance().registerRandomStream("UrQMD");
  corsika::UrQMD::UrQMD urqmd;

  std::vector<Code> const projectiles{{Code::Proton, Code::AntiProton, Code::Neutron,
                                       Code::AntiNeutron, Code::PiPlus, Code::PiMinus,
                                       Code::KPlus, Code::KMinus, Code::K0Short}};

  for (auto const& p : projectiles) {
    std::ofstream file(std::string("xs_") + GetName(p) + ".dat");
    for (auto Elab = GetMass(p) + 200_MeV; Elab <= 10_TeV; Elab *= 1.02) {
      file << Elab / 1_GeV << '\t'
           << urqmd.GetTabulatedCrossSection(p, Code::Nitrogen, Elab) / 1_mb << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
