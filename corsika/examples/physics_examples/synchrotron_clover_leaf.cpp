/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <corsika/output/OutputManager.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/refractivity/UniformRefractiveIndex.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <corsika/modules/radio/RadioProcess.hpp>
#include <corsika/modules/radio/CoREAS.hpp>
#include <corsika/modules/radio/ZHS.hpp>
#include <corsika/modules/radio/observers/Observer.hpp>
#include <corsika/modules/radio/observers/TimeDomainObserver.hpp>
#include <corsika/modules/radio/detectors/ObserverCollection.hpp>
#include <corsika/modules/radio/propagators/DummyTestPropagator.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>

#include <corsika/modules/TimeCut.hpp>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <typeinfo>

using namespace corsika;
using namespace std;

//
// A simple shower to get the electric field trace of an electron (& a positron)
// in order to reveal the "clover leaf" pattern of energy fluence to the ground.
//

int main() {

  logging::set_level(logging::level::info);
  corsika_logger->set_pattern("[%n:%^%-8l%$] custom pattern: %v");

  CORSIKA_LOG_INFO("Synchrotron radiation");

  feenableexcept(FE_INVALID);
  RNGManager<>::getInstance().registerRandomStream("cascade");
  std::random_device rd;
  auto seed = rd();
  RNGManager<>::getInstance().setSeed(seed);

  OutputManager output("clover_leaf_outputs");

  // create a suitable environment
  using IModelInterface = media::IRefractiveIndexModel<
      media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>>;
  using AtmModel = media::UniformRefractiveIndex<media::MediumPropertyModel<
      media::UniformMagneticField<media::HomogeneousMedium<IModelInterface>>>>;
  using EnvType = media::Environment<AtmModel>;
  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  // a refractive index for the vacuum
  const double ri_{1.};
  // the composition we use for the homogeneous medium
  media::NuclearComposition const Composition({Code::Nitrogen}, {1.});
  // create magnetic field vector
  auto const Bmag{0.00005_T};
  MagneticFieldVector B(rootCS, 0_T, Bmag, 0_T);
  // create a Sphere for the medium
  auto world = EnvType::createNode<Sphere>(center, 150_km);
  // set the environment properties
  world->setModelProperties<AtmModel>(ri_, media::Medium::AirDry1Atm, B,
                                      1_kg / (1_m * 1_m * 1_m), Composition);
  // bind things together
  env.getUniverse()->addChild(std::move(world));

  Point injectionPos(rootCS, 0_m, 0_m, 5_km);

  auto const centerX{center.getCoordinates().getX()};
  auto const centerY{center.getCoordinates().getY()};
  auto const centerZ{center.getCoordinates().getZ()};
  auto const injectionPosX_{injectionPos.getCoordinates().getX()};
  auto const injectionPosY_{injectionPos.getCoordinates().getY()};
  auto const injectionPosZ_{injectionPos.getCoordinates().getZ()};
  auto const triggerpoint_{Point(rootCS, injectionPosX_, injectionPosY_, injectionPosZ_)};

  // the observer characteristics
  const TimeType duration_{2e-6_s};            // 0.8e-4_s
  const InverseTimeType sampleRate_{1e+11_Hz}; // 1e+9_Hz

  // the detectors
  ObserverCollection<TimeDomainObserver> detectorCoREAS;

  std::string name_center = "CoREAS_R=0_m--Phi=0degrees";
  auto triggertime_center{((triggerpoint_ - center).getNorm() / constants::c) - 500_ns};
  TimeDomainObserver observer_center(name_center, center, rootCS, triggertime_center,
                                     duration_, sampleRate_, triggertime_center);
  detectorCoREAS.addObserver(observer_center);

  for (auto radius_1 = 25_m; radius_1 <= 1000_m; radius_1 += 25_m) {
    for (auto phi_1 = 0; phi_1 <= 315; phi_1 += 45) {
      auto phiRad_1 = phi_1 / 180. * M_PI;
      auto rr_1 = static_cast<int>(radius_1 / 1_m);
      auto const point_1{Point(rootCS, centerX + radius_1 * cos(phiRad_1),
                               centerY + radius_1 * sin(phiRad_1), centerZ)};
      CORSIKA_LOG_INFO("Observer point: {}", point_1);
      auto triggertime_1{((triggerpoint_ - point_1).getNorm() / constants::c) - 500_ns};
      std::string name_1 = "CoREAS_R=" + std::to_string(rr_1) +
                           "_m--Phi=" + std::to_string(phi_1) + "degrees";
      TimeDomainObserver observer_1(name_1, point_1, rootCS, triggertime_1, duration_,
                                    sampleRate_, triggertime_1);
      detectorCoREAS.addObserver(observer_1);
    }
  }

  // setup particle stack, and add primary particle
  setup::Stack<EnvType> stack;
  stack.clear();

  // electron
  const Code beamCode = Code::Electron;
  auto const charge = get_charge(beamCode);
  auto const mass = get_mass(beamCode);
  auto const gyroradius = 100_m;
  auto const pLabMag = convert_SI_to_HEP(charge * Bmag * gyroradius);
  auto const omega_inv =
      convert_HEP_to_SI<MassType::dimension_type>(mass) / (abs(charge) * Bmag);
  MomentumVector const plab{rootCS, 0_MeV, 0_MeV, -10_MeV};
  auto const Elab = sqrt(plab.getSquaredNorm() + static_pow<2>(mass));
  auto gamma = Elab / mass;
  TimeType const period = 2 * M_PI * omega_inv * gamma;

  // positron
  const Code beamCodeP = Code::Positron;
  auto const chargeP = get_charge(beamCodeP);
  auto const massP = get_mass(beamCodeP);
  // auto const gyroradius = 100_m;
  auto const pLabMagP = convert_SI_to_HEP(chargeP * Bmag * gyroradius);
  auto const omega_invP =
      convert_HEP_to_SI<MassType::dimension_type>(massP) / (abs(chargeP) * Bmag);
  MomentumVector const plabP{rootCS, 0_MeV, 0_MeV, -10_MeV};
  auto const ElabP = sqrt(plabP.getSquaredNorm() + static_pow<2>(massP));
  auto gammaP = ElabP / massP;
  TimeType const periodP = 2 * M_PI * omega_invP * gammaP;

  stack.addParticle(std::make_tuple(beamCode,
                                    calculate_kinetic_energy(plab.getNorm(), mass),
                                    plab.normalized(), injectionPos, 0_ns));
  stack.addParticle(std::make_tuple(beamCodeP,
                                    calculate_kinetic_energy(plabP.getNorm(), massP),
                                    plabP.normalized(), injectionPos, 0_ns));

  // setup relevant processes
  setup::Tracking tracking;

  // the radio signal propagator
  auto SP = make_dummy_test_radio_propagator(env);

  // put radio processes here
  RadioProcess<decltype(detectorCoREAS), CoREAS<decltype(detectorCoREAS), decltype(SP)>,
               decltype(SP)>
      coreas(detectorCoREAS, SP);
  output.add("CoREAS", coreas);

  TimeCut cut(period / 4);

  // assemble all processes into an ordered process list
  auto sequence = make_sequence(coreas, cut);

  output.startOfLibrary();
  // define air shower object, run simulation
  Cascade EAS(env, tracking, sequence, output, stack);
  EAS.run();

  CORSIKA_LOG_INFO("|p| electron = {} and E electron = {}", plab.getNorm(), Elab);
  CORSIKA_LOG_INFO("period: {}", period);
  CORSIKA_LOG_INFO("gamma: {}", gamma);

  CORSIKA_LOG_INFO("|p| positron = {} and E positron = {}", plabP.getNorm(), ElabP);
  CORSIKA_LOG_INFO("period: {}", periodP);
  CORSIKA_LOG_INFO("gamma: {}", gammaP);

  output.endOfLibrary();
}
