/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <iomanip>
#include <string>

using namespace corsika;
using namespace std;

//
// This example prints out all of the particles that are
// known to CORSIKA 8 and the corresponding IDs in the
// hadronic model codes.
//

int main() {

  logging::set_level(logging::level::info);
  corsika_logger->set_pattern("[%n:%^%-8l%$] %v");

  logging::info(
      "\n"
      "------------------------------------------\n"
      "        particles in CORSIKA\n"
      "------------------------------------------\n");
  int const width = 25 + 10 + 10 + 10 + 10 + 16 + 16 + 17;
  logging::info(
      "Name                      | "
      "C8 code    | "
      "PDG-id     | "
      "SIBYLL-id  | "
      "QGSJETII-id| "
      "PDG-mass (GeV)   | "
      "SIBYLL-mass (GeV)|");
  logging::info("{:->{}}", ' ', width);
  for (auto const p : get_all_particles()) {
    if (!is_nucleus(p)) {
      corsika::sibyll::SibyllCode sib_id = corsika::sibyll::convertToSibyll(p);
      auto const sib_mass = (sib_id != corsika::sibyll::SibyllCode::Unknown
                                 ? to_string(corsika::sibyll::getSibyllMass(p) / 1_GeV)
                                 : "");
      auto const qgs_id = corsika::qgsjetII::convertToQgsjetII(p);
      logging::info("{:25} | {:10} | {:10} | {:10} | {:10} | {:>16.5} | {:>16.5} |",
                    get_name(p), p, static_cast<int>(get_PDG(p)),
                    (sib_id != corsika::sibyll::SibyllCode::Unknown
                         ? to_string(static_cast<int>(sib_id))
                         : ""),
                    (qgs_id != corsika::qgsjetII::QgsjetIICode::Unknown
                         ? to_string(static_cast<int>(qgs_id))
                         : ""),
                    get_mass(p) / 1_GeV, sib_mass);
    }
  }
  logging::info("{:->{}}", ' ', width);
}
