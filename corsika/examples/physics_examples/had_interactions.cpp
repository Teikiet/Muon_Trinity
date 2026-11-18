/*
 * (c) Copyright 2019 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/Epos.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Pythia8.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Vector.hpp>

#include <corsika/framework/random/RNGManager.hpp>

#include <corsika/framework/core/EnergyMomentumOperations.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <corsika/modules/Random.hpp>

#include <fstream>
#include <iostream>
#include <tuple>

using namespace corsika;

using EnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using EnvType = media::Environment<EnvironmentInterface>;

template <typename TModel>
void create_events(TModel& model, std::string const& model_name, Code const projectileId,
                   FourMomentum projectileP4, Code const targetId, FourMomentum targetP4,
                   int const nEvents) {

  logging::set_level(logging::level::info);

  // calculate Ecm per nucleon (sqrtSNN). Used only for output!
  auto const projectileP4NN =
      projectileP4 / (is_nucleus(projectileId) ? get_nucleus_A(projectileId) : 1);
  auto const targetP4NN = targetP4 / (is_nucleus(targetId) ? get_nucleus_A(targetId) : 1);
  auto const SNN = (projectileP4NN + targetP4NN).getNormSqr();
  HEPEnergyType const sqrtSNN = sqrt(SNN);
  auto const sqrtS = (projectileP4 + targetP4).getNorm();
  CORSIKA_LOG_INFO("{} + {} interactions at sqrtSNN= {} GeV, sqrtS = {} GeV",
                   projectileId, targetId, sqrtSNN / 1_GeV, sqrtS / 1_GeV);

  // target material and environment (not used just there to create the C8 stack)
  EnvType env;
  auto& universe = *(env.getUniverse());
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  auto world = EnvType::createNode<Sphere>(Point{rootCS, 0_m, 0_m, 0_m}, 150_km);
  using MyHomogeneousModel = media::MediumPropertyModel<
      media::UniformMagneticField<media::HomogeneousMedium<EnvironmentInterface>>>;
  auto const props = world->setModelProperties<MyHomogeneousModel>(
      media::Medium::AirDry1Atm, MagneticFieldVector(rootCS, 0_T, 0_T, 0_T),
      1_kg / (1_m * 1_m * 1_m), media::NuclearComposition({targetId}, {1.}));
  world->setModelProperties(props);
  universe.addChild(std::move(world));

  // projectile (dummy)
  auto const plab3 = projectileP4.getSpaceLikeComponents();
  auto p0 = Point(rootCS, 0_m, 0_m, 0_m);
  setup::Stack<EnvType> stack;
  stack.addParticle(std::make_tuple(
      projectileId, calculate_kinetic_energy(plab3.getNorm(), get_mass(projectileId)),
      plab3.normalized(), p0, 0_ns));
  auto particle = stack.first();
  particle.setNode(universe.getContainingNode(p0));

  std::stringstream pname;
  pname << "particles-" << model_name << "-" << projectileId << "-" << targetId << "-"
        << sqrtSNN / 1_GeV << "GeV"
        << ".txt";

  unsigned int count_central_charged = 0;
  HEPMomentumType total_pt = 0_GeV;

  std::ofstream pout(pname.str());
  // add header
  pout << "# particle production in " << model_name << " for " << projectileId << " + "
       << targetId << " interactions."
       << "\n";
  pout << "# lab. momentum of projectile: "
       << projectileP4.getSpaceLikeComponents().getNorm() / 1_GeV << " (GeV)"
       << "\n";
  pout << "# CM energy per nucleon: " << sqrtSNN / 1_GeV << " (GeV)"
       << "\n";
  pout << "# i, j, particle name, PDG id, p.rap, rap., energy (GeV), px, py, pz (GeV), "
          "pt (GeV), charge"
       << "\n";
  for (int i = 0; i < nEvents; ++i) {

    // empty stack from previous interaction
    stack.clear();
    // add particle to get StackView to write secondaries to C8 stack
    stack.addParticle(std::make_tuple(
        projectileId, calculate_kinetic_energy(plab3.getNorm(), get_mass(projectileId)),
        plab3.normalized(), p0, 0_ns)); // irrelevant
    auto particle = stack.first();
    particle.setNode(universe.getContainingNode(p0));
    setup::Stack<EnvType>::stack_view_type output(particle);

    // generate secondaries and fill output
    model->doInteraction(output, projectileId, targetId, projectileP4, targetP4);

    // loop through secondaries and calculate higher level variables (rapidity, pt, ...),
    // write to file
    int j = -1;
    HEPMomentumType tpt = 0_GeV;
    for (auto& secondary : output) {
      ++j;
      // pseudorapidity, rapidity
      auto const p = secondary.getMomentum();
      auto const pz = p.getZ(rootCS);
      auto const px = p.getX(rootCS);
      auto const py = p.getY(rootCS);
      auto const pabs = p.getNorm();
      auto const pt = sqrt(p.getSquaredNorm() - square(pz));
      auto const en = secondary.getEnergy();
      auto const mt = sqrt(square(pt) + square(get_mass(secondary.getPID())));
      double prap, rap;
      if (pz / 1_GeV <= 0) {
        rap = log(mt / (en - pz));
        prap = log((pt + 1.0_eV) / (pabs - pz));
      } else {
        rap = log((en + pz) / mt);
        prap = log((pabs + pz) / (pt + 1.0_eV));
      }
      pout << i << " " << j << " " << secondary.getPID() << " "
           << static_cast<int>(get_PDG(secondary.getPID())) << " " << prap << " " << rap
           << " " << en / 1_GeV << " " << px / 1_GeV << " " << py / 1_GeV << " "
           << pz / 1_GeV << " " << pt / 1_GeV << " "
           << get_charge_number(secondary.getPID()) << "\n";
      // calculate plateau (CMS experiment)
      if (abs(prap) < 2.5 && get_charge_number(secondary.getPID()) != 0) {
        count_central_charged++;
      }
      tpt += pt;
    }
    total_pt += tpt / float(j);
  }
  std::cout << "average plateau height: " << count_central_charged / float(nEvents)
            << std::endl;
  std::cout << "average pt (GeV): " << total_pt / float(nEvents) / 1_GeV << std::endl;
  std::cout << "wrote hadronic events for " << model_name << " in file: " << pname.str()
            << std::endl;
};

int main(int argc, char** argv) {

  std::string int_model, projectile, target;
  HEPMomentumType P0;
  int Nevents = 10;
  if (argc < 2 || argc > 7) {
    std::cout << "usage: check <interaction model> \n valid models are: sibyll, "
                 "epos, qgsjet, pythia8. Other options are (in that order): <projectile> "
                 "<target> <projectile momentum (lab in GeV)> <no. of events>. Defaults "
                 "are: proton, proton, 100, 10"
              << std::endl;
    return 1;
  } else if (argc == 2) {
    int_model = std::string(argv[1]);
    projectile = "proton";
    target = "proton";
    P0 = 100_GeV;
  } else if (argc == 3) {
    int_model = std::string(argv[1]);
    projectile = std::string(argv[2]);
    target = "proton";
    P0 = 100_GeV;
  } else if (argc == 4) {
    int_model = std::string(argv[1]);
    projectile = std::string(argv[2]);
    target = std::string(argv[3]);
    P0 = 100_GeV;
  } else if (argc == 5) {
    int_model = std::string(argv[1]);
    projectile = std::string(argv[2]);
    target = std::string(argv[3]);
    P0 = std::stoi(std::string(argv[4])) * 1_GeV;
  } else if (argc == 6) {
    int_model = std::string(argv[1]);
    projectile = std::string(argv[2]);
    target = std::string(argv[3]);
    P0 = std::stoi(std::string(argv[4])) * 1_GeV;
    Nevents = std::stoi(std::string(argv[5]));
  }

  Code projectileId;
  if (projectile == "proton")
    projectileId = Code::Proton;
  else if (projectile == "antiproton")
    projectileId = Code::AntiProton;
  else if (projectile == "pion")
    projectileId = Code::PiPlus;
  else if (projectile == "oxygen")
    projectileId = Code::Oxygen;
  else {
    std::cout << "unknown projectile. pick proton, antiproton, pion or oxygen"
              << std::endl;
    return 0;
  }
  Code targetId;
  if (target == "proton")
    targetId = Code::Proton;
  else if (target == "nitrogen")
    targetId = Code::Nitrogen;
  else {
    std::cout << "unknown target. pick proton or nitrogen" << std::endl;
    return 0;
  }

  auto const& rootCS = get_root_CoordinateSystem();
  FourMomentum const P4proj{
      calculate_total_energy(get_mass(projectileId), P0),
      {rootCS, {HEPMomentumType::zero(), HEPMomentumType::zero(), P0}}};

  FourMomentum const P4targ{
      calculate_total_energy(get_mass(targetId), HEPMomentumType::zero()),
      {rootCS,
       {HEPMomentumType::zero(), HEPMomentumType::zero(), HEPMomentumType::zero()}}};

  /*
   the interface to hadronic models in C8 has two layers. The class InteractionModel
   (HadronInteractionModel in the case of sibyll) implements the interface to the
   interaction model. It provides the basic interface functions getCrossSection and
   doInteraction. The class Interaction extends the interface in InteractionModel to an
   element C8 process sequence. Here we use only the interface class.
  */
  if (int_model == "sibyll") {
    RNGManager<>::getInstance().registerRandomStream("sibyll");
    // this is the hadron model in sibyll. can only do hadron-nucleus interactions
    auto model = std::make_shared<corsika::sibyll::HadronInteractionModel>(
        corsika::setup::C7trackedParticles);
    create_events(model, int_model, projectileId, P4proj, targetId, P4targ, Nevents);
  } else if (int_model == "epos") {
    RNGManager<>::getInstance().registerRandomStream("epos");
    auto model = std::make_shared<corsika::epos::InteractionModel>(
        corsika::setup::C7trackedParticles);
    create_events(model, int_model, projectileId, P4proj, targetId, P4targ, Nevents);
  } else if (int_model == "pythia8") {
    RNGManager<>::getInstance().registerRandomStream("pythia");
    auto model = std::make_shared<corsika::pythia8::InteractionModel>(
        corsika::setup::C7trackedParticles);

    create_events(model, int_model, projectileId, P4proj, targetId, P4targ, Nevents);
  } else {
    RNGManager<>::getInstance().registerRandomStream("qgsjet");
    auto model = std::make_shared<corsika::qgsjetII::InteractionModel>();
    create_events(model, int_model, projectileId, P4proj, targetId, P4targ, Nevents);
  }
}
