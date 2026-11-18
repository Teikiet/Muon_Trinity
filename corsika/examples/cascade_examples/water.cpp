/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/* clang-format off */
// InteractionCounter used boost/histogram, which
// fails if boost/type_traits have been included before. Thus, we have
// to include it first...
#include <corsika/framework/process/InteractionCounter.hpp>
/* clang-format on */
#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/medium/MediumProperties.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/FLUKA.hpp>
#include <corsika/modules/tracking/TrackingStraight.hpp>

#include <corsika/output/OutputManager.hpp>
#include <corsika/stack/GeometryNodeStackExtension.hpp>
#include <corsika/stack/WeightStackExtension.hpp>
#include <corsika/stack/VectorStack.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

using namespace corsika;

using IMediumType = media::IMediumPropertyModel<media::IMediumModel>;
using EnvType = media::Environment<IMediumType>;
using StackType = setup::Stack<EnvType>;
using TrackingType = tracking_line::Tracking;
using Particle = StackType::particle_type;

void registerRandomStreams(int seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("epos");
  RNGManager<>::getInstance().registerRandomStream("pythia");
  RNGManager<>::getInstance().registerRandomStream("fluka");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {} ", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {} ", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
}

int main(int argc, char** argv) {
  // * process input
  Code beamCode;
  HEPEnergyType E0, eCut;
  int A, Z, n_event;
  int randomSeed;
  std::string output_dir;
  CLI::App app{"Cascade in water"};
  // we start by definining a sub-group for the primary ID
  auto opt_Z = app.add_option("-Z", Z, "Atomic number for primary")
                   ->check(CLI::Range(0, 26))
                   ->group("Primary");
  auto opt_A = app.add_option("-A", A, "Atomic mass number for primary")
                   ->needs(opt_Z)
                   ->check(CLI::Range(1, 58))
                   ->group("Primary");
  app.add_option("-p,--pdg", "PDG code for primary.")
      ->excludes(opt_A)
      ->excludes(opt_Z)
      ->group("Primary");
  app.add_option("-E,--energy", "Primary energy in GeV")
      ->required()
      ->check(CLI::PositiveNumber)
      ->group("Primary");
  app.add_option("--eCut", "Cut energy in GeV")->default_val(1.);
  app.add_option("-N,--nevent", n_event, "The number of events/showers to run.")
      ->default_val(1)
      ->check(CLI::PositiveNumber);
  app.add_option("-f,--filename", output_dir, "Filename for output library")
      ->check(CLI::NonexistentPath)
      ->default_val("output");
  app.add_option("-v", "Verbosity level: warn, info, debug, trace.")
      ->default_val("info")
      ->check(CLI::IsMember({"warn", "info", "debug", "trace"}));
  app.add_option("-s", randomSeed, "Seed for random number")
      ->check(CLI::NonNegativeNumber)
      ->default_val(0);

  // parse the command line options into the variables
  CLI11_PARSE(app, argc, argv);

  std::string_view const loglevel = app["-v"]->as<std::string_view>();
  if (loglevel == "warn") {
    logging::set_level(logging::level::warn);
  } else if (loglevel == "info") {
    logging::set_level(logging::level::info);
  } else if (loglevel == "debug") {
    logging::set_level(logging::level::debug);
  } else if (loglevel == "trace") {
#ifndef _C8_DEBUG_
    CORSIKA_LOG_ERROR("trace log level requires a Debug build.");
    return 1;
#endif
    logging::set_level(logging::level::trace);
  }

  // check that we got either PDG or A/Z
  // this can be done with option_groups but the ordering
  // gets all messed up
  if (app.count("--pdg") == 0) {
    if ((app.count("-A") == 0) || (app.count("-Z") == 0)) {
      CORSIKA_LOG_ERROR("If --pdg is not provided, then both -A and -Z are required.");
      return 1;
    }
  }

  // initialize random number sequence(s)
  registerRandomStreams(randomSeed);

  // check if we want to use a PDG code instead
  if (app.count("--pdg") > 0) {
    beamCode = convert_from_PDG(PDGCode(app["--pdg"]->as<int>()));
  } else {
    // check manually for proton and neutrons
    if ((A == 1) && (Z == 1))
      beamCode = Code::Proton;
    else if ((A == 1) && (Z == 0))
      beamCode = Code::Neutron;
    else
      beamCode = get_nucleus_code(A, Z);
  }

  eCut = app["--eCut"]->as<double>() * 1_GeV;

  // * environment and universe
  EnvType env;
  auto& universe = env.getUniverse();
  auto const& rootCS = env.getCoordinateSystem();

  // * Water geometry
  {
    Point const center{rootCS, 0_m, 0_m, 0_m};
    auto sphere = std::make_unique<Sphere>(center, 100_m);
    auto node = std::make_unique<media::VolumeTreeNode<IMediumType>>(std::move(sphere));
    media::NuclearComposition const nuclearComposition{{Code::Hydrogen, Code::Oxygen},
                                                       {2.0 / 3.0, 1.0 / 3.0}};
    // density of sea water
    auto density = 1.02_g / (1_cm * 1_cm * 1_cm);
    auto water_medium = std::make_shared<
        media::MediumPropertyModel<media::HomogeneousMedium<IMediumType>>>(
        media::Medium::WaterLiquid, density, nuclearComposition);
    node->setModelProperties(water_medium);
    universe->addChild(std::move(node));
  }

  // * make downward-going shower axis and a observation plane in x-y-plane
  auto injectorLength = 50_m;
  Point const injectionPos = Point(rootCS, {0_m, 0_m, injectorLength});
  auto const& injectCS = make_translation(rootCS, injectionPos.getCoordinates());
  DirectionVector upVec(rootCS, {0., 0., 1.});
  DirectionVector leftVec(rootCS, {1., 0., 0.});
  DirectionVector downVec(rootCS, {0., 0., -1.});

  std::vector<ObservationPlane<TrackingType, ParticleWriterParquet>> obsPlanes;
  const int nPlane = 5;
  for (int i = 0; i < nPlane - 1; i++) {
    Point planeCenter{injectCS, {0_m, 0_m, -(i + 1) * 3_m}};
    obsPlanes.push_back({Plane(planeCenter, upVec), leftVec, false});
  }
  auto& obsPlaneFinal = obsPlanes.emplace_back(
      Plane{Point{injectCS, {0_m, 0_m, -50_m}}, upVec}, leftVec, true);

  // * longitutional profile
  media::ShowerAxis const showerAxis{injectionPos, 1.2 * injectorLength * downVec, env};
  auto const dX = 1_g / square(1_cm); // Binning of the writers along the shower axis
  LongitudinalWriter longiWriter{showerAxis, dX};
  LongitudinalProfile<SubWriter<decltype(longiWriter)>> longprof{longiWriter};

  // * energy loss profile
  EnergyLossWriter dEdX{showerAxis, dX};

  // * physical process list
  // particle production threshold
  HEPEnergyType const emCut = eCut;
  HEPEnergyType const hadCut = eCut;
  ParticleCut<SubWriter<decltype(dEdX)>> cut(emCut, emCut, hadCut, hadCut, hadCut, true,
                                             dEdX);

  // tell proposal that we are interested in all energy losses above the particle cut
  set_energy_production_threshold(Code::Electron, std::min({emCut, hadCut}));
  set_energy_production_threshold(Code::Positron, std::min({emCut, hadCut}));
  set_energy_production_threshold(Code::Photon, std::min({emCut, hadCut}));
  set_energy_production_threshold(Code::MuMinus, std::min({emCut, hadCut}));
  set_energy_production_threshold(Code::MuPlus, std::min({emCut, hadCut}));
  set_energy_production_threshold(Code::TauMinus, std::min({emCut, hadCut}));
  set_energy_production_threshold(Code::TauPlus, std::min({emCut, hadCut}));

  // hadronic interactions
  HEPEnergyType heHadronModelThreshold = std::pow(10, 1.9) * 1_GeV;
  corsika::sibyll::Interaction sibyll(corsika::media::get_all_elements_in_universe(env),
                                      corsika::setup::C7trackedParticles);

  auto const all_elements = corsika::media::get_all_elements_in_universe(env);
  corsika::fluka::Interaction leIntModel{all_elements};
  InteractionCounter leIntCounted{leIntModel};
  struct EnergySwitch {
    HEPEnergyType cutE_;
    EnergySwitch(HEPEnergyType cutE)
        : cutE_(cutE) {}
    bool operator()(const Particle& p) const { return (p.getKineticEnergy() < cutE_); }
  };
  auto hadronSequence =
      make_select(EnergySwitch(heHadronModelThreshold), leIntCounted, sibyll);

  // decay process
  corsika::pythia8::Decay decayPythia;

  corsika::sophia::InteractionModel sophia;

  // EM process
  corsika::proposal::Interaction emCascade(
      env, sophia, sibyll.getHadronInteractionModel(), heHadronModelThreshold);
  corsika::proposal::ContinuousProcess<SubWriter<decltype(dEdX)>> emContinuous(env, dEdX);

  // total physics list
  auto physics_sequence =
      make_sequence(emCascade, emContinuous, hadronSequence, decayPythia);

  // * output module
  OutputManager output(output_dir);
  for (int i = 0; i < nPlane; i++) {
    output.add(fmt::format("particles_{:}", i), obsPlanes[i]);
  }
  // hard coded
  auto obsPlaneSequence =
      make_sequence(obsPlanes[0], obsPlanes[1], obsPlanes[2], obsPlanes[3], obsPlanes[4]);

  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(obsPlanes.back());
  output.add("primary", primaryWriter);

  output.add("profile", longiWriter);
  output.add("energyloss", dEdX);

  // * the final process sequence
  auto sequence = make_sequence(physics_sequence, longprof, obsPlaneSequence, cut);

  // * tracking and stack
  TrackingType tracking;
  StackType stack;

  // * cascade manager
  Cascade EAS(env, tracking, sequence, output, stack);

  E0 = app["-E"]->as<double>() * 1_GeV;
  HEPEnergyType mass = get_mass(beamCode);
  // convert Elab to Plab
  HEPMomentumType P0 = calculate_momentum(E0, mass);
  auto plab = MomentumVector(rootCS, {0_GeV, 0_GeV, -P0});

  // print our primary parameters all in one place
  if (app["--pdg"]->count() > 0) {
    CORSIKA_LOG_INFO("Primary PDG ID:     {}", app["--pdg"]->as<int>());
  } else {
    CORSIKA_LOG_INFO("Primary Z/A:        {}/{}", Z, A);
  }
  CORSIKA_LOG_INFO("Primary Energy:     {}", E0);
  CORSIKA_LOG_INFO("Primary Momentum:   {}", P0);
  CORSIKA_LOG_INFO("Primary Direction:  {}", plab.getNorm());
  CORSIKA_LOG_INFO("Point of Injection: {}", injectionPos.getCoordinates());
  CORSIKA_LOG_INFO("Shower Axis Length: {}", injectorLength);

  // * main loop
  output.startOfLibrary();
  for (int i_shower = 0; i_shower < n_event; i_shower++) {
    stack.clear();
    CORSIKA_LOG_INFO("Event: {} / {}", i_shower, n_event);

    // * inject primary
    auto const primaryProperties = std::make_tuple(
        beamCode, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
        plab.normalized(), injectionPos, 0_ns);
    auto primary = stack.addParticle(primaryProperties);
    stack.addParticle(primaryProperties);

    EAS.run();

    // * report energy loss result
    HEPEnergyType const Efinal = dEdX.getEnergyLost() + obsPlaneFinal.getEnergyGround();
    CORSIKA_LOG_INFO(
        "total energy budget (TeV): {:.2f} (dEdX={:.2f} ground={:.2f}), "
        "relative difference (%): {:.3f}",
        E0 / 1_TeV, dEdX.getEnergyLost() / 1_TeV, obsPlaneFinal.getEnergyGround() / 1_TeV,
        (Efinal / E0 - 1.) * 100.);
  }
  output.endOfLibrary();

  return 0;
}
