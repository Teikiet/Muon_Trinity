/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/process/InteractionCounter.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>
// #include <corsika/framework/utility/CorsikaFenv.hpp>

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

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/TrackWriter.hpp>
#include <corsika/modules/UrQMD.hpp>
#include <corsika/modules/CONEX.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace corsika;
using namespace std;

//
// An example of running an EAS where the hadronic cascade is
// handled by sibyll+URQMD and the EM cascade is treated with
// CONEX + Bethe Bloch (as opposed to PROPOSAL).
//

/**
 * Random number stream initialization
 *
 * @param seed
 */
void registerRandomStreams(uint64_t seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("conex");
  RNGManager<>::getInstance().registerRandomStream("epos");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("pythia");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {} ", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {} ", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
}

/**
 * New (for demonstration) ContinuousProcess which will check if a particles has traversed
 * below the observation level.
 */
class TrackCheck : public ContinuousProcess<TrackCheck> {

public:
  /**
   * Construct a new Track Check object.
   *
   * @param plane -- the actual observation level
   */
  TrackCheck(Plane const& plane)
      : plane_(plane) {}

  /**
   * The doContinous method to check a particular particle.
   *
   * @tparam TParticle
   * @tparam TTrack
   * @param particle
   * @return ProcessReturn
   */
  template <typename TParticle>
  ProcessReturn doContinuous(Step<TParticle> const& step, bool const) {
    auto const delta = step.getParticlePre().getPosition() - plane_.getCenter();
    auto const n = plane_.getNormal();
    auto const proj = n.dot(delta);
    if (proj < -1_m) {
      CORSIKA_LOG_INFO("particle {} failes: proj={}, delta={}, p={}",
                       step.getParticlePre().asString(), proj, delta,
                       step.getPositionPost());
      throw std::runtime_error("particle below obs level");
    }
    return ProcessReturn::Ok;
  }

  /**
   * No limit on tracking step length imposed here, of course.
   *
   * @tparam TParticle
   * @tparam TTrack
   * @return LengthType
   */
  template <typename TParticle, typename TTrack>
  LengthType getMaxStepLength(TParticle const&, TTrack const&) const {
    return std::numeric_limits<double>::infinity() * 1_m;
  }

private:
  Plane plane_;
};

/**
 * Selection of environment interface implementation:
 */
using EnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using EnvType = media::Environment<EnvironmentInterface>;
template <typename T>
using MyExtraEnv = media::MediumPropertyModel<media::UniformMagneticField<T>>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;

int main(int argc, char** argv) {

  logging::set_level(logging::level::info);

  CORSIKA_LOG_INFO("hybrid_MC");

  if (argc < 4) {
    CORSIKA_LOG_ERROR(
        "\n"
        "usage: hybrid_MC <A> <Z> <energy/GeV> [seed] \n"
        "       if no seed is given, a random seed is chosen");
    return 1;
  }

  uint64_t seed = 0;
  if (argc > 4) seed = std::stol(std::string(argv[4]));
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

  unsigned short const A = std::stoi(std::string(argv[1]));
  unsigned short const Z = std::stoi(std::string(argv[2]));
  Code const beamCode = get_nucleus_code(A, Z);
  auto const mass = get_mass(beamCode);
  HEPEnergyType const E0 = 1_GeV * std::stof(std::string(argv[3]));
  double theta = 0.;
  auto const thetaRad = theta / 180. * M_PI;

  HEPMomentumType P0 = calculate_momentum(E0, mass);
  auto momentumComponents = [](double thetaRad, HEPMomentumType ptot) {
    return std::make_tuple(ptot * sin(thetaRad), 0_eV, -ptot * cos(thetaRad));
  };

  auto const [px, py, pz] = momentumComponents(thetaRad, P0);
  auto plab = MomentumVector(rootCS, {px, py, pz});

  auto const observationHeight = 0_km + constants::EarthRadius::Mean;
  auto const injectionHeight = 112.75_km + constants::EarthRadius::Mean;
  auto const t = -observationHeight * cos(thetaRad) +
                 sqrt(-static_pow<2>(sin(thetaRad) * observationHeight) +
                      static_pow<2>(injectionHeight));
  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
  Point const injectionPos =
      showerCore +
      Vector<dimensionless_d>{rootCS, {-sin(thetaRad), 0, cos(thetaRad)}} * t;

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

  // SETUP WRITERS

  OutputManager output("hybrid_MC_outputs");

  // register energy losses as output
  EnergyLossWriter dEdX{showerAxis, dX};
  output.add("energyloss", dEdX);

  // create a track writer and register it with the output manager
  TrackWriter<TrackWriterParquet> tracks;
  output.add("tracks", tracks);

  ParticleCut<SubWriter<decltype(dEdX)>> cut(3_GeV, false, dEdX);
  BetheBlochPDG<SubWriter<decltype(dEdX)>> eLoss(dEdX);

  LongitudinalWriter profile{showerAxis, dX};
  output.add("profile", profile);
  LongitudinalProfile<SubWriter<decltype(profile)>> longprof{profile};

  Plane const obsPlane(showerCore, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType> observationLevel(obsPlane,
                                                  DirectionVector(rootCS, {1., 0., 0.}));
  output.add("particles", observationLevel);

  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

  // SETUP PROCESSES, DECAYS, INTERACTIONS

  corsika::sibyll::Interaction sibyll(corsika::media::get_all_elements_in_universe(env),
                                      corsika::setup::C7trackedParticles);
  InteractionCounter sibyllCounted{sibyll};

  corsika::pythia8::Decay decayPythia;

  CONEXhybrid // SubWriter<decltype(dEdX>, SubWriter<decltype(profile)>>
      conex_model(center, showerAxis, t, injectionHeight, E0, get_PDG(Code::Proton), dEdX,
                  profile);

  corsika::urqmd::UrQMD urqmd_model;
  InteractionCounter urqmdCounted{urqmd_model};

  TrackCheck trackCheck(obsPlane);

  // assemble all processes into an ordered process list
  struct EnergySwitch {
    HEPEnergyType cutE_;
    EnergySwitch(HEPEnergyType cutE)
        : cutE_(cutE) {}
    bool operator()(const StackType::particle_type& p) const {
      return (p.getEnergy() < cutE_);
    }
  };
  auto hadronSequence = make_select(EnergySwitch(55_GeV), urqmdCounted, sibyllCounted);
  auto sequence = make_sequence(hadronSequence, decayPythia, eLoss, cut, conex_model,
                                longprof, observationLevel, trackCheck);

  output.startOfLibrary();

  StackType stack;
  stack.clear();

  // define air shower object, run simulation
  TrackingType tracking;
  Cascade EAS(env, tracking, sequence, output, stack);

  auto const primaryProperties = std::make_tuple(
      Code::Proton, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
      plab.normalized(), injectionPos, 0_ns);

  stack.addParticle(primaryProperties);
  primaryWriter.recordPrimary(primaryProperties);

  // to fix the point of first interaction, uncomment the following two lines:
  //  EAS.SetNodes();
  //  EAS.forceInteraction();

  EAS.run();

  const HEPEnergyType Efinal = dEdX.getEnergyLost() + observationLevel.getEnergyGround();
  CORSIKA_LOG_INFO(
      "total cut energy (GeV): {}, "
      "relative difference (%): {}",
      Efinal / 1_GeV, (Efinal / E0 - 1) * 100);

  auto const hists = sibyllCounted.getHistogram() + urqmdCounted.getHistogram();

  save_hist(hists.labHist(), "inthist_lab_hybrid.npz", true);
  save_hist(hists.CMSHist(), "inthist_cms_hybrid.npz", true);

  output.endOfLibrary();

  CORSIKA_LOG_INFO("done");
}
