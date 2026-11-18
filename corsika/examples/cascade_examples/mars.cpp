/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
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
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>

#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/framework/process/InteractionCounter.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>

#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/FlatExponential.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/density_and_composition/SlidingPlanarExponential.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/StackInspector.hpp>
#include <corsika/modules/TrackWriter.hpp>
#include <corsika/modules/FLUKA.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace corsika;
using namespace std;

using EnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using EnvType = media::Environment<EnvironmentInterface>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;

typedef decltype(1 * pascal) PressureType;
typedef decltype(1 * degree_celsius) TemperatureType;

class MarsAtmModel {
public:
  MarsAtmModel() = delete;
  MarsAtmModel(PressureType a, InverseLengthType b, TemperatureType c,
               decltype(1 * degree_celsius / 1_m) d)
      : a_(a)
      , b_(b)
      , c_(c)
      , d_(d) {}

  MassDensityType operator()(LengthType height) const {
    PressureType const pressure = a_ * exp(-b_ * height);
    TemperatureType const temperature = -c_ - d_ * height + 273.1_K; // in K
    constexpr decltype(square(1_m) / (square(1_s) * 1_K)) constant =
        1000 * 0.1921 * square(1_m) / (square(1_s) * 1_K);
    return pressure / (constant * temperature);
  }

private:
  PressureType a_;
  InverseLengthType b_;
  TemperatureType c_;
  decltype(1_K / 1_m) d_;
};

void registerRandomStreams(int seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
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

template <typename T>
using MyExtraEnv = media::MediumPropertyModel<media::UniformMagneticField<T>>;

int main(int argc, char** argv) {

  // the main command line description
  CLI::App app{"Simulate standard (downgoing) showers with CORSIKA 8."};

  // some options that we want to fill in
  int A, Z, nevent = 0;

  // the following section adds the options to the parser

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
  // the remainding options
  app.add_option("-E,--energy", "Primary energy in GeV")
      ->required()
      ->check(CLI::PositiveNumber)
      ->group("Primary");
  app.add_option("-z,--zenith", "Primary zenith angle (deg)")
      ->required()
      ->default_val(0.)
      ->check(CLI::Range(0, 90))
      ->group("Primary");
  app.add_option("-a,--azimuth", "Primary azimuth angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0, 360))
      ->group("Primary");
  app.add_option("-N,--nevent", nevent, "The number of events/showers to run.")
      ->required()
      ->check(CLI::PositiveNumber)
      ->group("Library/Output");
  app.add_option("-f,--filename", "Filename for output library.")
      ->required()
      ->default_val("corsika_library")
      ->check(CLI::NonexistentPath)
      ->group("Library/Output");
  app.add_option("-s,--seed", "The random number seed.")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");
  app.add_option("-v,--verbosity", "Verbosity level: warn, info, debug, trace.")
      ->default_val("info")
      ->check(CLI::IsMember({"warn", "info", "debug", "trace"}))
      ->group("Misc.");

  // parse the command line options into the variables
  CLI11_PARSE(app, argc, argv);

  if (app.count("--verbosity")) {
    auto const loglevel = app["--verbosity"]->as<std::string>();
    if (loglevel == "warn") {
      logging::set_level(logging::level::warn);
    } else if (loglevel == "info") {
      logging::set_level(logging::level::info);
    } else if (loglevel == "debug") {
      logging::set_level(logging::level::debug);
    } else if (loglevel == "trace") {

#ifndef C8_DEBUG
      CORSIKA_LOG_ERROR("trace log level requires a Debug build.");
      return 1;
#endif
      logging::set_level(logging::level::trace);
    }
  }

  // check that we got either PDG or A/Z
  if (app.count("--pdg") == 0) {
    if ((app.count("-A") == 0) || (app.count("-Z") == 0)) {
      CORSIKA_LOG_ERROR("If --pdg is not provided, then both -A and -Z are required.");
      return 1;
    }
  }

  // initialize random number sequence(s)
  registerRandomStreams(app["--seed"]->as<int>());

  /* === START: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */
  EnvType env;
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  LengthType const radiusMars = 3389.5_km;
  auto builder =
      media::make_layered_spherical_atmosphere_builder<EnvironmentInterface, MyExtraEnv>::
          create(center,
                 radiusMars,                                   // Mars
                 media::Medium::AirDry1Atm,                    // Mars, close enough
                 MagneticFieldVector{rootCS, 0_T, 0_uT, 0_T}); // Mars

  builder.setNuclearComposition(                   // Mars
      {{Code::Carbon, Code::Oxygen,                // 95.97 CO2
        Code::Nitrogen},                           // 1.89 N2 + 1.93 Argon + 0.146 O2
       {0.9597 / 3, 0.9597 * 2 / 3, 1 - 0.9597}}); // values taken from AIRES manual

  MarsAtmModel layer1(0.699e3 * pascal, 0.00009 / 1_m, 31.0 * degree_celsius,
                      0.000998 * 1 * degree_celsius / 1_m);
  MarsAtmModel layer2(0.699e3 * pascal, 0.00009 / 1_m, 23.4 * degree_celsius,
                      0.00222 * 1 * degree_celsius / 1_m);

  builder.addTabularLayer(layer1, 100, 100_m, 7_km);
  builder.addTabularLayer(layer2, 300, 500_m, 100_km);
  builder.addLinearLayer(1_g / square(1_cm), 1e9_cm, 112.8_km);
  builder.assemble(env);
  /* === END: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */

  ofstream atmout("mars.dat");
  for (LengthType h = 0_m; h < 110_km; h += 10_m) {
    Point const ptest{rootCS, 0_m, 0_m, builder.getPlanetRadius() + h};
    auto rho =
        env.getUniverse()->getContainingNode(ptest)->getModelProperties().getMassDensity(
            ptest);
    atmout << h / 1_m << " " << rho / 1_kg * cube(1_m) << "\n";
  }
  atmout.close();

  /* === START: CONSTRUCT PRIMARY PARTICLE === */

  // parse the primary ID as a PDG or A/Z code
  Code beamCode;

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
  HEPEnergyType const mass = get_mass(beamCode);

  // particle energy
  HEPEnergyType const E0 = 1_GeV * app["--energy"]->as<double>();

  // direction of the shower in (theta, phi) space
  auto const thetaRad = app["--zenith"]->as<double>() / 180. * M_PI;
  auto const phiRad = app["--azimuth"]->as<double>() / 180. * M_PI;

  // convert Elab to Plab
  HEPMomentumType P0 = sqrt((E0 - mass) * (E0 + mass));

  // convert the momentum to the zenith and azimuth angle of the primary
  auto const [px, py, pz] =
      std::make_tuple(P0 * sin(thetaRad) * cos(phiRad), P0 * sin(thetaRad) * sin(phiRad),
                      -P0 * cos(thetaRad));
  auto plab = MomentumVector(rootCS, {px, py, pz});
  /* === END: CONSTRUCT PRIMARY PARTICLE === */

  /* === START: CONSTRUCT GEOMETRY === */
  auto const observationHeight = 0_km + builder.getPlanetRadius();
  auto const injectionHeight = 111.75_km + builder.getPlanetRadius();
  auto const t = -observationHeight * cos(thetaRad) +
                 sqrt(-static_pow<2>(sin(thetaRad) * observationHeight) +
                      static_pow<2>(injectionHeight));
  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
  Point const injectionPos =
      showerCore + DirectionVector{rootCS,
                                   {-sin(thetaRad) * cos(phiRad),
                                    -sin(thetaRad) * sin(phiRad), cos(thetaRad)}} *
                       t;

  // we make the axis much longer than the inj-core distance since the
  // profile will go beyond the core, depending on zenith angle
  /* === END: CONSTRUCT GEOMETRY === */

  // create the output manager that we then register outputs with
  OutputManager output(app["--filename"]->as<std::string>());

  media::ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos) * 1.2,
                                     env};
  auto const dX = 10_g / square(1_cm); // Binning of the writers along the shower axis

  EnergyLossWriter dEdX{showerAxis, dX};
  output.add("energyloss", dEdX);

  HEPEnergyType const emcut = 1_GeV;
  HEPEnergyType const hadcut = 1_GeV;
  ParticleCut<SubWriter<decltype(dEdX)>> cut(emcut, emcut, hadcut, hadcut, hadcut, true,
                                             dEdX);

  // tell proposal that we are interested in all energy losses above the particle cut
  set_energy_production_threshold(Code::Electron, std::min({emcut, hadcut}));
  set_energy_production_threshold(Code::Positron, std::min({emcut, hadcut}));
  set_energy_production_threshold(Code::Photon, std::min({emcut, hadcut}));
  set_energy_production_threshold(Code::MuMinus, std::min({emcut, hadcut}));
  set_energy_production_threshold(Code::MuPlus, std::min({emcut, hadcut}));
  set_energy_production_threshold(Code::TauMinus, std::min({emcut, hadcut}));
  set_energy_production_threshold(Code::TauPlus, std::min({emcut, hadcut}));

  /* === START: SETUP PROCESS LIST === */
  auto const all_elements = corsika::media::get_all_elements_in_universe(env);
  corsika::sibyll::Interaction sibyll(all_elements, corsika::setup::C7trackedParticles);
  InteractionCounter sibyllCounted(sibyll);

  corsika::pythia8::Decay decayPythia;

  // energy threshold for high energy hadronic model. Affects LE/HE switch for hadron
  // interactions and the hadronic photon model in proposal
  HEPEnergyType heHadronModelThreshold = 63.1_GeV;

  corsika::sophia::InteractionModel sophia;
  corsika::proposal::Interaction emCascade(
      env, sophia, sibyll.getHadronInteractionModel(), heHadronModelThreshold);

  // use BetheBlochPDG for hadronic continuous losses, and proposal otherwise
  corsika::proposal::ContinuousProcess<SubWriter<decltype(dEdX)>> emContinuousProposal(
      env, dEdX);
  BetheBlochPDG<SubWriter<decltype(dEdX)>> emContinuousBethe{dEdX};
  struct EMHadronSwitch {
    EMHadronSwitch() = default;
    bool operator()(const Particle& p) const { return is_hadron(p.getPID()); }
  };
  auto emContinuous =
      make_select(EMHadronSwitch(), emContinuousBethe, emContinuousProposal);

  LongitudinalWriter longprof{showerAxis, dX};
  output.add("profile", longprof);
  LongitudinalProfile<SubWriter<decltype(longprof)>> profile{longprof};

  corsika::fluka::Interaction leIntModel{all_elements};
  InteractionCounter leIntCounted{leIntModel};
  StackInspector<StackType> stackInspect(5000, false, E0);

  // assemble all processes into an ordered process list
  struct EnergySwitch {
    HEPEnergyType cutE_;
    EnergySwitch(HEPEnergyType cutE)
        : cutE_(cutE) {}
    bool operator()(Particle const& p) const { return (p.getKineticEnergy() < cutE_); }
  };
  auto hadronSequence =
      make_select(EnergySwitch(heHadronModelThreshold), leIntCounted, sibyllCounted);

  // track writer
  TrackWriter trackWriter;
  output.add("tracks", trackWriter); // register TrackWriter

  // observation plane
  Plane const obsPlane(showerCore, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType> observationLevel(obsPlane,
                                                  DirectionVector(rootCS, {1., 0., 0.}));
  // register the observation plane with the output
  output.add("particles", observationLevel);

  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

  // assemble the final process sequence
  auto sequence =
      make_sequence(stackInspect, hadronSequence, decayPythia, emCascade, emContinuous,
                    trackWriter, profile, observationLevel, cut);
  /* === END: SETUP PROCESS LIST === */

  // create the cascade object using the default stack and tracking
  // implementation
  TrackingType tracking;
  StackType stack;
  Cascade EAS(env, tracking, sequence, output, stack);

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
  CORSIKA_LOG_INFO("Shower Axis Length: {}", (showerCore - injectionPos).getNorm() * 1.2);

  // trigger the output manager to open the library for writing
  output.startOfLibrary();

  // loop over each shower
  for (int i_shower = 1; i_shower < nevent + 1; i_shower++) {
    CORSIKA_LOG_INFO("Shower {} / {} ", i_shower, nevent);

    // trigger the start of the outputs for this shower
    output.startOfShower();

    // directory for outputs
    string const outdir(app["--filename"]->as<std::string>());
    string const labHist_file = outdir + "/inthist_lab_" + to_string(i_shower) + ".npz";
    string const cMSHist_file = outdir + "/inthist_cms_" + to_string(i_shower) + ".npz";

    // setup particle stack, and add primary particle
    stack.clear();

    // add the desired particle to the stack
    auto const primaryProperties = std::make_tuple(
        beamCode, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
        plab.normalized(), injectionPos, 0_ns);
    stack.addParticle(primaryProperties);

    primaryWriter.recordPrimary(primaryProperties);

    // run the shower
    EAS.run();

    HEPEnergyType const Efinal =
        dEdX.getEnergyLost() + observationLevel.getEnergyGround();

    CORSIKA_LOG_INFO(
        "total energy budget (GeV): {}, "
        "relative difference (%): {}",
        Efinal / 1_GeV, (Efinal / E0 - 1) * 100);

    // trigger the output manager to save this shower to disk
    output.endOfShower();
  }

  // and finalize the output on disk
  output.endOfLibrary();
}
