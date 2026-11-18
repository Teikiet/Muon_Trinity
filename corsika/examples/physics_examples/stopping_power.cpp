/*
 * (c) Copyright 2019 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>

#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/setup/SetupStack.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <fstream>
#include <iostream>
#include <limits>

using namespace corsika;
using namespace std;

//
// This example demonstrates the energy loss of muons as function of beta*gamma (=p/m)
//
int main() {

  logging::set_level(logging::level::info);

  CORSIKA_LOG_INFO("stopping_power");

  feenableexcept(FE_INVALID);

  // setup environment, geometry
  using EnvType = media::Environment<media::IMediumModel>;
  EnvType env;
  env.getUniverse()->setModelProperties<media::HomogeneousMedium<media::IMediumModel>>(
      1_g / cube(1_cm), media::NuclearComposition{{Code::Unknown}, {1.f}});

  CoordinateSystem const& rootCS = env.getCoordinateSystem();

  Point const injectionPos(
      rootCS, 0_m, 0_m,
      112.8_km); // this is the CORSIKA 7 start of atmosphere/universe

  BetheBlochPDG eLoss;

  setup::Stack<EnvType> stack;

  std::string fileName = "dEdX.dat";
  CORSIKA_LOG_INFO("Writing to file {}", fileName);
  std::ofstream file(fileName);
  file << "# beta*gamma, dE/dX / MeV/(g/cm²)" << std::endl;

  for (HEPEnergyType E0 = 200_MeV; E0 < 1_PeV; E0 *= 1.05) {
    stack.clear();
    const Code beamCode = Code::MuPlus;
    const HEPMassType mass = get_mass(beamCode);
    double theta = 0.;
    double phi = 0.;

    HEPMomentumType P0 = calculate_momentum(E0, mass);
    auto momentumComponents = [](double theta, double phi, HEPMomentumType ptot) {
      return std::make_tuple(ptot * sin(theta) * cos(phi), ptot * sin(theta) * sin(phi),
                             -ptot * cos(theta));
    };
    auto const [px, py, pz] =
        momentumComponents(theta / 180. * constants::pi, phi / 180. * constants::pi, P0);
    auto plab = MomentumVector(rootCS, {px, py, pz});

    stack.addParticle(std::make_tuple(
        beamCode, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
        plab.normalized(), injectionPos, 0_ns));

    auto const p = stack.getNextParticle();
    HEPEnergyType dE = eLoss.getTotalEnergyLoss(p, 1_g / square(1_cm));
    file << P0 / mass << "\t" << -dE / 1_MeV << std::endl;
  }
  CORSIKA_LOG_INFO("finished writing {}", fileName);
}
