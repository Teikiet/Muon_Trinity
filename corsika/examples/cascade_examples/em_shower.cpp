/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/process/InteractionCounter.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>

#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/TrackWriter.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <typeinfo>

using namespace corsika;
using namespace std;

void registerRandomStreams(int seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {} ", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {} ", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
}

using EnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using EnvType = media::Environment<EnvironmentInterface>;
template <typename T>
using MyExtraEnv = media::MediumPropertyModel<media::UniformMagneticField<T>>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;

int main(int argc, char** argv) {

  logging::set_level(logging::level::warn);

  if (!(argc == 2 || argc == 3)) {
    std::cerr << "usage: em_shower <energy/GeV> [seed]" << std::endl
              << "seed = 0 for randomized seed" << std::endl;
    return 1;
  }
  feenableexcept(FE_INVALID);
  int seed = 0;

  if (argc >= 3) { seed = std::stoi(std::string(argv[2])); }
  // initialize random number sequence(s)
  registerRandomStreams(seed);

  // setup environment, geometry
  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};

  // build a Linsley US Standard atmosphere into `env`
  MagneticFieldVector bField{rootCS, 50_uT, 0_T, 0_T};
  media::create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, media::AtmosphereId::LinsleyUSStd, center, media::Medium::AirDry1Atm, bField);

  std::unordered_map<Code, HEPEnergyType> energy_resolution = {
      {Code::Electron, 5_MeV},
      {Code::Positron, 5_MeV},
      {Code::Photon, 5_MeV},
  };
  for (auto const& [pcode, energy] : energy_resolution)
    set_energy_production_threshold(pcode, energy);

  const Code beamCode = Code::Electron;
  auto const mass = get_mass(beamCode);
  const HEPEnergyType E0 = 1_GeV * std::stof(std::string(argv[1]));
  double theta = 0.;
  auto const thetaRad = theta / 180. * M_PI;

  HEPMomentumType P0 = calculate_momentum(E0, mass);
  auto momentumComponents = [](double thetaRad, HEPMomentumType ptot) {
    return std::make_tuple(ptot * sin(thetaRad), 0_eV, -ptot * cos(thetaRad));
  };

  auto const [px, py, pz] = momentumComponents(thetaRad, P0);
  auto plab = MomentumVector(rootCS, {px, py, pz});

  auto const observationHeight = 0.0_km + constants::EarthRadius::Mean;
  auto const injectionHeight = 112.75_km + constants::EarthRadius::Mean;
  auto const t = -observationHeight * cos(thetaRad) +
                 sqrt(-static_pow<2>(sin(thetaRad) * observationHeight) +
                      static_pow<2>(injectionHeight));
  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
  Point const injectionPos =
      showerCore + DirectionVector{rootCS, {-sin(thetaRad), 0, cos(thetaRad)}} * t;

  media::ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos) * 1.02,
                                     env, false, 1000};
  auto const dX = 10_g / square(1_cm); // Binning of the writers along the shower axis

  CORSIKA_LOG_INFO("Primary particle:   {}", beamCode);
  CORSIKA_LOG_INFO("Zenith angle:       {} (rad)", theta);
  CORSIKA_LOG_INFO("Momentum:           {} (GeV)", plab.getComponents() / 1_GeV);
  CORSIKA_LOG_INFO("Propagation dir:    {}", plab.getNorm());
  CORSIKA_LOG_INFO("Injection point:    {}", injectionPos.getCoordinates());
  CORSIKA_LOG_INFO("shower axis length: {} ",
                   (showerCore - injectionPos).getNorm() * 1.02);

  // setup processes, decays and interactions
  EnergyLossWriter energyloss{showerAxis, dX};
  ParticleCut<SubWriter<decltype(energyloss)>> cut(5_MeV, 5_MeV, 100_GeV, 100_GeV,
                                                   100_GeV, true, energyloss);

  corsika::sibyll::Interaction sibyll(corsika::media::get_all_elements_in_universe(env),
                                      corsika::setup::C7trackedParticles);
  corsika::sophia::InteractionModel sophia;
  HEPEnergyType heThresholdNN = 80_GeV;
  corsika::proposal::Interaction emCascade(
      env, sophia, sibyll.getHadronInteractionModel(), heThresholdNN);
  corsika::proposal::ContinuousProcess<SubWriter<decltype(energyloss)>> emContinuous(
      env, energyloss);

  //  NOT possible right now, due to interface differenc in PROPOSAL
  //  InteractionCounter emCascadeCounted(emCascade);

  OutputManager output("em_shower_outputs");

  output.add("energyloss", energyloss);

  TrackWriter tracks;
  output.add("tracks", tracks);

  LongitudinalWriter profile{showerAxis, dX};
  output.add("profile", profile);
  LongitudinalProfile<SubWriter<decltype(profile)>> longprof{profile};

  Plane const obsPlane(showerCore, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType, ParticleWriterParquet> observationLevel{
      obsPlane, DirectionVector(rootCS, {1., 0., 0.})};
  output.add("particles", observationLevel);

  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

  auto sequence = make_sequence(emCascade, emContinuous, longprof, observationLevel, cut);
  // define air shower object, run simulation
  TrackingType tracking;

  output.startOfLibrary();

  auto const primaryProperties = std::make_tuple(
      beamCode, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
      plab.normalized(), injectionPos, 0_ns);

  // setup particle stack, and add primary particle
  StackType stack;
  stack.clear();
  stack.addParticle(primaryProperties);
  primaryWriter.recordPrimary(primaryProperties);

  Cascade EAS(env, tracking, sequence, output, stack);

  // to fix the point of first interaction, uncomment the following two lines:
  //  EAS.forceInteraction();

  EAS.run();

  HEPEnergyType const Efinal =
      energyloss.getEnergyLost() + observationLevel.getEnergyGround();

  CORSIKA_LOG_INFO(
      "total energy budget (GeV): {}, "
      "relative difference (%): {}",
      Efinal / 1_GeV, (Efinal / E0 - 1) * 100);

  output.endOfLibrary();

  return EXIT_SUCCESS;
}
