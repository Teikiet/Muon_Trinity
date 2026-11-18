/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <boost/filesystem/path.hpp>
#include "corsika/framework/core/EnergyMomentumOperations.hpp"

namespace corsika::pythia8 {

  inline NeutrinoInteraction::NeutrinoInteraction(std::set<Code> const& list,
                                                  bool const& handleNC,
                                                  bool const& handleCC)
      : stable_particles_(list)
      , handle_nc_(handleNC)
      , handle_cc_(handleCC)
      , pythiaMain_{CORSIKA_Pythia8_XML_DIR, false} {

    pythiaMain_.setRndmEnginePtr(std::make_shared<corsika::pythia8::Random>());
  }

  inline NeutrinoInteraction::~NeutrinoInteraction() {
    CORSIKA_LOGGER_DEBUG(logger_, "Pythia::NeutrinoInteraction n= {}", count_);
  }

  template <class TView>
  void NeutrinoInteraction::doInteraction(TView& view, Code const projectileId,
                                          Code const targetId,
                                          FourMomentum const& projectileP4,
                                          FourMomentum const& targetP4) {

    CORSIKA_LOGGER_DEBUG(logger_, "Primary {} - {} interaction with E_nu = {} GeV",
                         projectileId, targetId,
                         projectileP4.getTimeLikeComponent() / 1_GeV);
    CORSIKA_LOGGER_DEBUG(
        logger_, "configure Pythia for primary neutrino interactions. NC={}, CC={}",
        handle_nc_, handle_cc_);
    if (!handle_nc_ && !handle_cc_) {
      CORSIKA_LOGGER_ERROR(
          logger_,
          "no neutrino interaction channel configured! Select either NC, CC or both!");
      throw std::runtime_error("Configuration error!");
    }
    CORSIKA_LOGGER_DEBUG(logger_, "minimal Q2 in DIS: {} GeV2", minQ2_ / 1_GeV / 1_GeV);

    if (!isValid(projectileId, targetId, projectileP4, targetP4)) {
      CORSIKA_LOGGER_ERROR(logger_, "wrong projectile, target or energy configuration!");
      throw std::runtime_error("Configuration error!");
    }
    // sample nucleon from nucleus A,Z
    double const fProtons = get_nucleus_Z(targetId) / double(get_nucleus_A(targetId));
    double const fNeutrons = 1. - fProtons;
    std::discrete_distribution<int> nucleonChannelDist{fProtons, fNeutrons};
    corsika::default_prng_type& rng =
        corsika::RNGManager<>::getInstance().getRandomStream("pythia");
    Code const targetNucleonId = (nucleonChannelDist(rng) ? Code::Neutron : Code::Proton);
    int const idTarget_pythia = static_cast<int>(get_PDG(targetNucleonId));
    auto const targetP4NN =
        targetP4 / (is_nucleus(targetId) ? get_nucleus_A(targetId) : 1);
    CORSIKA_LOGGER_DEBUG(logger_, "selected {} target. P4NN={}", targetNucleonId,
                         targetP4NN.getSpaceLikeComponents());

    // set projectile
    double const Q2min = minQ2_ / 1_GeV / 1_GeV;
    int const idProjectile_pythia = static_cast<int>(get_PDG(projectileId));

    // calculate CoM energy of the neutrino nucleon interaction
    double const ecm_pythia =
        calculate_com_energy(projectileP4.getTimeLikeComponent(), get_mass(projectileId),
                             get_mass(targetNucleonId)) /
        1_GeV;

    CORSIKA_LOGGER_DEBUG(logger_, "center-of-mass energy: {} GeV", ecm_pythia);
    // Set up CoM interaction
    pythiaMain_.readString("Beams:frameType = 1");
    // center-of-mass energy
    pythiaMain_.settings.parm("Beams:eCM", ecm_pythia);
    // BeamA (+z in CoM)
    pythiaMain_.settings.mode("Beams:idA", idProjectile_pythia);
    // BeamB (-z direction in CoM)
    pythiaMain_.settings.mode("Beams:idB", idTarget_pythia);

    // Set up DIS process within some phase space.
    // Neutral current (with gamma/Z interference).
    if (handle_nc_) pythiaMain_.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
    // Charged current.
    if (handle_cc_) pythiaMain_.readString("WeakBosonExchange:ff2ff(t:W) = on");
    // Phase-space cut: minimal Q2 of process.
    pythiaMain_.settings.parm("PhaseSpace:Q2Min", Q2min);

    // Set dipole recoil on. Necessary for DIS + shower.
    pythiaMain_.readString("SpaceShower:dipoleRecoil = on");

    // Allow emissions up to the kinematical limit,
    // since rate known to match well to matrix elements everywhere.
    pythiaMain_.readString("SpaceShower:pTmaxMatch = 2");

    // QED radiation off lepton not handled yet by the new procedure.
    pythiaMain_.readString("PDF:lepton = off");
    pythiaMain_.readString("TimeShower:QEDshowerByL = off");

    // switch on decays for all hadrons except the ones defined as tracked by C8
    if (!stable_particles_.empty()) {
      pythiaMain_.readString("HadronLevel:Decay = on");
      for (auto pCode : stable_particles_)
        pythiaMain_.particleData.mayDecay(static_cast<int>(get_PDG(pCode)), false);
    } else {
      // all hadrons stable
      pythiaMain_.readString("HadronLevel:Decay = off");
    }

    pythiaMain_.readString("Stat:showProcessLevel = off");
    pythiaMain_.readString("Stat:showPartonLevel = off");

    pythiaMain_.readString("Print:quiet = on");
    pythiaMain_.readString("Check:epTolErr = 0.1");
    pythiaMain_.readString("Check:epTolWarn = 0.0001");
    pythiaMain_.readString("Check:mTolErr = 0.01");

    pythiaMain_.init();

    // References to the event record
    Pythia8::Event& eventMain = pythiaMain_.event;
    // define boost assuming target nucleon is at rest
    COMBoost const boost{projectileP4, targetP4NN};
    // the boost is along the total momentum axis and does not include a rotation
    // get rotated frame where momentum of projectile in the CoM is along +z

    if (!pythiaMain_.next()) {
      throw std::runtime_error("Pythia neutrino collision failed ");
    } else {
      CORSIKA_LOGGER_DEBUG(logger_, "pythia neutrino interaction done!");
    }

    MomentumVector Plab_final{};
    auto Elab_final = HEPEnergyType::zero();
    CORSIKA_LOGGER_DEBUG(logger_, "particles generated in neutrino interaction:");
    for (int i = 0; i < eventMain.size(); ++i) {
      auto const& p8p = eventMain[i];
      if (p8p.isFinal()) {
        try {
          // get particle ids from pythia stack and convert to corsika id
          auto const volatile id = static_cast<PDGCode>(p8p.id());
          auto const pyId = convert_from_PDG(id);
          // get pythia momentum in CoM (rotated cs!)
          MomentumVector const pyPcom(boost.makeCoMVector(QuantityVector<hepmomentum_d>{
              p8p.px() * 1_GeV, p8p.py() * 1_GeV, p8p.pz() * 1_GeV}));
          // boost pythia 4Momentum in CoM to lab. frame. This includes the rotation.
          // calculate 3-momentum and kinetic energy
          CORSIKA_LOGGER_DEBUG(logger_, "momentum in CoM: {} GeV",
                               pyPcom.getComponents() / 1_GeV);
          auto const pyP4lab = boost.fromCoM(FourVector{p8p.e() * 1_GeV, pyPcom});
          auto const pyPlab = pyP4lab.getSpaceLikeComponents();
          HEPEnergyType const pyEkinlab =
              calculate_kinetic_energy(pyPlab.getNorm(), get_mass(pyId));

          // add to corsika stack
          auto pnew =
              view.addSecondary(std::make_tuple(pyId, pyEkinlab, pyPlab.normalized()));

          CORSIKA_LOGGER_DEBUG(logger_, "id = {}, E = {} GeV, p = {} GeV", pyId,
                               pyEkinlab / 1_GeV, pyPlab.getComponents() / 1_GeV);
          Plab_final += pnew.getMomentum();
          Elab_final += pnew.getEnergy();
        }
        // irreproducible in tests, LCOV_EXCL_START
        catch (std::out_of_range const& ex) {
          CORSIKA_LOGGER_CRITICAL(logger_, "Pythia ID {} unknown in C8", p8p.id());
          throw ex;
        }
        // LCOV_EXCL_STOP
      }
    }

    CORSIKA_LOGGER_DEBUG(logger_,
                         "conservation (all GeV): "
                         "Elab_final= {}"
                         ", Plab_final= {}",
                         Elab_final / 1_GeV, (Plab_final / 1_GeV).getComponents());

    count_++;
  }

} // namespace corsika::pythia8