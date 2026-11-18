/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/process/InteractionCounter.hpp>
#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>

#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/refractivity/UniformRefractiveIndex.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>

#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/TrackWriter.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/PROPOSAL.hpp>

#include <corsika/modules/radio/RadioProcess.hpp>
#include <corsika/modules/radio/CoREAS.hpp>
#include <corsika/modules/radio/ZHS.hpp>
#include <corsika/modules/radio/observers/Observer.hpp>
#include <corsika/modules/radio/observers/TimeDomainObserver.hpp>
#include <corsika/modules/radio/detectors/ObserverCollection.hpp>
#include <corsika/modules/radio/propagators/NumericalIntegratingPropagator.hpp>
#include <corsika/modules/radio/propagators/DummyTestPropagator.hpp>

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

//
// An example of running an casacde which generates radio input for a
// cascade that has hadronic interactions turned off. Currently, this
// will produce output for only one observer, but there are lines below,
// which are commented out, to enable a full star-shape pattern of observers
//

void registerRandomStreams(int seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {} ", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {} ", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
}

using EnvironmentInterface = media::IRefractiveIndexModel<
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>>;
using EnvType = media::Environment<EnvironmentInterface>;
template <typename TInterface>
using MyExtraEnv = media::UniformRefractiveIndex<
    media::MediumPropertyModel<media::UniformMagneticField<TInterface>>>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;

int main(int argc, char** argv) {

  logging::set_level(logging::level::warn);

  if (argc != 3) {
    std::cerr
        << "usage: radio_em_shower <energy/GeV> <seed> - put seed=0 to use random seed"
        << std::endl;
    return 1;
  }

  int seed{static_cast<int>(std::stof(std::string(argv[2])))};
  CORSIKA_LOG_INFO("Seed {}", seed);
  feenableexcept(FE_INVALID);
  // initialize random number sequence(s)
  registerRandomStreams(seed);

  // setup environment, geometry
  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};

  double const refractive_index = 1.000327;
  MagneticFieldVector const bField{rootCS, 50_uT, 0_T, 0_T};
  media::create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, media::AtmosphereId::LinsleyUSStd, center, refractive_index,
      media::Medium::AirDry1Atm, bField);

  std::unordered_map<Code, HEPEnergyType> energy_resolution = {
      {Code::Electron, 5_MeV},
      {Code::Positron, 5_MeV},
      {Code::Photon, 5_MeV},
  };
  for (auto const& [pcode, energy] : energy_resolution)
    set_energy_production_threshold(pcode, energy);

  Code const beamCode = Code::Electron;
  auto const mass = get_mass(beamCode);
  HEPEnergyType const E0 = 1_GeV * std::stof(std::string(argv[1]));
  double const theta = 0.;
  auto const thetaRad = theta / 180. * M_PI;

  HEPMomentumType const P0 = calculate_momentum(E0, mass);
  auto momentumComponents = [](double thetaRad, HEPMomentumType ptot) {
    return std::make_tuple(ptot * sin(thetaRad), 0_eV, -ptot * cos(thetaRad));
  };

  auto const [px, py, pz] = momentumComponents(thetaRad, P0);
  auto const plab = MomentumVector(rootCS, {px, py, pz});

  auto const observationHeight = 1.4_km + constants::EarthRadius::Mean;
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

  // setup the radio observers
  TimeType const groundHitTime{(showerCore - injectionPos).getNorm() / constants::c};

  // Radio observers and relevant information
  // the observer time variables
  TimeType const duration{1e-6_s};
  InverseTimeType const sampleRate{1e+9_Hz};

  // the detector (aka observer collection) for CoREAS and ZHS
  ObserverCollection<TimeDomainObserver> detectorCoREAS;
  ObserverCollection<TimeDomainObserver> detectorZHS;

  auto const showerCoreX{showerCore.getCoordinates().getX()};
  auto const showerCoreY{showerCore.getCoordinates().getY()};
  auto const injectionPosX{injectionPos.getCoordinates().getX()};
  auto const injectionPosY{injectionPos.getCoordinates().getY()};
  auto const injectionPosZ{injectionPos.getCoordinates().getZ()};
  auto const triggerpoint{Point(rootCS, injectionPosX, injectionPosY, injectionPosZ)};

  CORSIKA_LOG_INFO("Primary particle:   {}", beamCode);
  CORSIKA_LOG_INFO("Zenith angle:       {} (rad)", theta);
  CORSIKA_LOG_INFO("Momentum:           {} (GeV)", plab.getComponents() / 1_GeV);
  CORSIKA_LOG_INFO("Propagation dir:    {}", plab.getNorm());
  CORSIKA_LOG_INFO("Injection point:    {}", injectionPos.getCoordinates());
  CORSIKA_LOG_INFO("shower axis length: {} ",
                   (showerCore - injectionPos).getNorm() * 1.02);
  CORSIKA_LOG_INFO("Trigger Point is:   {}", triggerpoint);

  // // setup CoREAS observers - use the for loop for star shape pattern
  // for (auto radius_coreas = 25_m; radius_coreas <= 500_m; radius_coreas += 25_m) {
  //   for (auto phi_coreas = 0; phi_coreas <= 315; phi_coreas += 45) {
  auto radius_coreas = 200_m;
  auto phi_coreas = 45;
  auto phiRad_coreas = phi_coreas / 180. * M_PI;
  auto rr_coreas = static_cast<int>(radius_coreas / 1_m);
  auto const point_coreas{Point(rootCS, showerCoreX + radius_coreas * cos(phiRad_coreas),
                                showerCoreY + radius_coreas * sin(phiRad_coreas),
                                constants::EarthRadius::Mean)};
  std::cout << "Observer point: " << point_coreas << std::endl;
  CORSIKA_LOG_INFO("Observer point    {}", injectionPos.getCoordinates());
  auto triggertime_coreas{(triggerpoint - point_coreas).getNorm() / constants::c};
  std::string name_coreas = "CoREAS_R=" + std::to_string(rr_coreas) +
                            "_m--Phi=" + std::to_string(phi_coreas) + "degrees";
  TimeDomainObserver observer_coreas(name_coreas, point_coreas, rootCS,
                                     triggertime_coreas, duration, sampleRate,
                                     triggertime_coreas);
  detectorCoREAS.addObserver(observer_coreas);
  //   }
  // }

  // // setup ZHS observers - use the for loop for star shape pattern
  // for (auto radius_zhs = 25_m; radius_zhs <= 500_m; radius_zhs += 25_m) {
  //   for (auto phi_zhs = 0; phi_zhs <= 315; phi_zhs += 45) {
  auto radius_zhs = 200_m;
  auto phi_zhs = 45;
  auto phiRad_zhs = phi_zhs / 180. * M_PI;
  auto rr_zhs = static_cast<int>(radius_zhs / 1_m);
  auto const point_zhs{Point(rootCS, showerCoreX + radius_zhs * cos(phiRad_zhs),
                             showerCoreY + radius_zhs * sin(phiRad_zhs),
                             constants::EarthRadius::Mean)};
  auto triggertime_zhs{(triggerpoint - point_zhs).getNorm() / constants::c};
  std::string name_zhs = "ZHS_R=" + std::to_string(rr_zhs) +
                         "_m--Phi=" + std::to_string(phi_zhs) + "degrees";
  TimeDomainObserver observer_zhs(name_zhs, point_zhs, rootCS, triggertime_zhs, duration,
                                  sampleRate, triggertime_zhs);
  detectorZHS.addObserver(observer_zhs);
  //   }
  // }

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

  OutputManager output("radio_em_shower_outputs");

  output.add("energyloss", energyloss);

  TrackWriter tracks;
  output.add("tracks", tracks);

  LongitudinalWriter profile{showerAxis, dX};
  output.add("profile", profile);
  LongitudinalProfile<SubWriter<decltype(profile)>> longprof{profile};

  // the radio signal propagator
  auto SP = make_dummy_test_radio_propagator(env);

  // initiate CoREAS
  RadioProcess<decltype(detectorCoREAS), CoREAS<decltype(detectorCoREAS), decltype(SP)>,
               decltype(SP)>
      coreas(detectorCoREAS, SP);

  // register CoREAS with the output manager
  output.add("CoREAS", coreas);

  // initiate ZHS
  RadioProcess<decltype(detectorZHS), ZHS<decltype(detectorZHS), decltype(SP)>,
               decltype(SP)>
      zhs(detectorZHS, SP);

  // register ZHS with the output manager
  output.add("ZHS", zhs);

  Plane const obsPlane(showerCore, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType, ParticleWriterParquet> observationLevel{
      obsPlane, DirectionVector(rootCS, {1., 0., 0.})};
  output.add("particles", observationLevel);

  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

  auto sequence = make_sequence(emCascade, emContinuous, longprof, coreas, zhs, tracks,
                                observationLevel, cut);
  // define air shower object, run simulation
  TrackingType tracking;

  auto const primaryProperties = std::make_tuple(
      beamCode, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
      plab.normalized(), injectionPos, 0_ns);

  output.startOfLibrary();

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
