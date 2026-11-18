/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/TAUOLA.hpp>
#include <SetupTestEnvironment.hpp>
#include <SetupTestStack.hpp>
#include <catch2/catch_all.hpp>
#include <boost/histogram.hpp>

#include <boost/format.hpp>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <corsika/modules/Pythia8.hpp>
#include <corsika/framework/random/RNGManager.hpp>

using namespace corsika::tauola;
using namespace corsika;
using namespace boost::histogram;
using Catch::Approx;

using DummyEnvironmentInterface =
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

template <typename TStackView>
auto sumMomentum(TStackView const& view, CoordinateSystem const& vCS) {
  MomentumVector sum{vCS, 0_eV, 0_eV, 0_eV};
  for (auto const& p : view) { sum += p.getMomentum(); }
  return sum;
}

TEST_CASE("TAUOLA", "modules") {
  //  Tpp::Tauola::setSeed(0,0,0);
  logging::set_level(logging::level::info);
  RNGManager<>::getInstance().registerRandomStream("tauola");

  // test that we can create and initialize the TAUOLA library.
  SECTION("Create TAUOLA Process") {
    Decay();
    Decay(Helicity::Unpolarized, false);
    Decay(Helicity::RightHanded, false);
    Decay(Helicity::LeftHanded, false);
    Decay(Helicity::Unpolarized, true);
    Decay(Helicity::Unpolarized, true, false);
    Decay(Helicity::Unpolarized, true, true);
    Decay(Helicity::Unpolarized, false, false);
  }

  // test the link with TAUOLA
  SECTION("linking TAUOLA") {
    double mass = Tpp::parmas_.amtau, px = 0, py = 0, pz = 1000;
    double p = sqrt(px * px + py * py + pz * pz);
    double E = sqrt(mass * mass + p * p);
    auto particle = std::make_unique<TauolaInterfaceParticle>(15, 1, mass, px, py, pz, E);
    Tpp::Tauola::decayOne(particle.get(), false, 0, 0, 0);

    particle->print();

    if ((particle->getDaughters()).size() == 0)
      CORSIKA_LOG_CRITICAL("decay failed!");
    else
      CORSIKA_LOG_DEBUG("tau decay into {} particles:",
                        (particle->getDaughters()).size());
    CHECK(particle->getDaughters().size() > 0);
    for (auto* daughter : particle->getDaughters()) {
      CHECK(daughter->getMothers().size() > 0);
    }

    // just checking the interface here
    CHECK(particle->getBarcode() == -1);
    particle->setMass(particle->findMassTauola(Code::RhoPlus) / 1_GeV);
    particle->setPdgID(static_cast<PDGCodeIntType>(get_PDG(Code::RhoPlus)));
    particle->setStatus(0);
    CHECK(particle->getStatus() == 0);
    CHECK(particle->getPdgID() == static_cast<PDGCodeIntType>(get_PDG(Code::RhoPlus)));

    particle->setMass(particle->findMassTauola(Code::KStarPlus) / 1_GeV);

    CHECK(particle->findMassTauola(Code::PiPlus) == get_mass(Code::PiPlus));
  }

  SECTION("TAUOLA lifetime") {
    auto const projectileCode = GENERATE(Code::PiMinus, Code::TauMinus, Code::TauPlus);
    //  setup environment and coordinate system
    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Proton);

    // Momentum of injection
    HEPEnergyType const P0 = 10_GeV;

    corsika::tauola::Decay decay;

    // Setup injected particle
    auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
        projectileCode, P0, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
    auto& stack = *stackPtr;
    auto const& particle = stack.getNextParticle();

    // basic sanity check lifetime is non-zero
    CHECK(decay.getLifetime(particle) > 0_ns);
  }
  // test corsika particles decay with TAUOLA
  SECTION("TAUOLA decay") {
    //  setup environment and coordinate system
    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Proton);
    [[maybe_unused]] auto const& env_dummy = env;
    auto const& cs = *csPtr;

    // Momentum of injection
    HEPEnergyType const P0 = 10_GeV;

    // Check decay of tau+ and tau-
    std::vector<Code> pids;
    pids.push_back(Code::TauMinus);
    pids.push_back(Code::TauPlus);

    for (auto& pid : pids) {
      CORSIKA_LOG_INFO("\n\nDecay test: {}", pid);

      // Instantiate tauola decay module
      corsika::tauola::Decay decay;

      // Setup injected particle
      auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
          pid, P0, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
      auto& stack = *stackPtr;
      auto& view = *secViewPtr;
      auto const& particle = stack.getNextParticle();
      auto const plab = MomentumVector(cs, {P0, 0_GeV, 0_GeV});

      // Check boosted lifetime
      CORSIKA_LOG_INFO("stack: {} {}", stack.asString(), particle.asString());
      const TimeType lifetime = decay.getLifetime(particle);
      double const gamma = particle.getEnergy() / get_mass(pid);
      CHECK(lifetime == get_lifetime(pid) * gamma);

      // Decay
      decay.doDecay(view);
      CHECK(stack.getEntries() > 1);
      CORSIKA_LOG_INFO("{}->{}", pid, stack.asString());

      // Total energy and momentum of decay products
      MomentumVector pSum{cs, 0_eV, 0_eV, 0_eV};
      HEPEnergyType eSum{0_GeV};
      for (auto const& p : view) {
        // Check time is right
        CHECK(p.getTime() == 0_ns);
        pSum += p.getMomentum();
        eSum += p.getEnergy();
      }

      CORSIKA_LOG_INFO("Momentum before decay: {:.3f} GeV, after decay: {:.3f} GeV",
                       plab.getNorm() / 1_GeV, pSum.getNorm() / 1_GeV);
      CORSIKA_LOG_INFO("Energy before decay: {:.3f} GeV, after decay: {:.3f} GeV",
                       particle.getEnergy() / 1_GeV, eSum / 1_GeV);

      // Check momentum and energy conservation
      CHECK((pSum - plab).getNorm() / 1_GeV == Approx(0).margin(1e-3));
      CHECK((pSum.getNorm() - plab.getNorm()) / 1_GeV == Approx(0).margin(1e-3));
      CHECK((particle.getEnergy() - eSum) / 1_GeV == Approx(0).margin(1e-3));
    }
  }

  // Compare performance of TAUOLA and PYTHIA on decay of tau-
  SECTION("TAUOLA vs. PYTHIA") {
    //  setup environment and coordinate system
    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Proton);
    [[maybe_unused]] auto const& env_dummy = env;
    auto const& cs = *csPtr;

    // Num of events
    const int N{1000};

    // Bin number for histgrams
    const unsigned int Nbins{50};

    // ========== TAUOLA ========

    // Create histograms
    // pid distribution
    auto tauola_pid = make_histogram(axis::regular<>(8000, -4000, 4000));
    // Histograms of energy for different particle
    auto tauola_electron = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_muon = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_nue = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_numu = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_nutau = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_photon = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_pizero = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_piplus = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto tauola_piminus = make_histogram(axis::regular<>(Nbins, 0.0, 1000));

    // Record time consumption
    std::chrono::steady_clock::time_point tauola_begin = std::chrono::steady_clock::now();

    // Test tau- decays
    const auto pcode{Code::TauMinus};
    const auto mass{get_mass(pcode)};

    // Energy of injection
    HEPEnergyType const E0 = 1000_GeV;
    HEPEnergyType const P0 = sqrt(E0 * E0 - mass * mass);

    // Instantiate tauola decay module
    corsika::tauola::Decay decay(Helicity::Unpolarized);

    for (int iEvent = 0; iEvent < N; iEvent++) {

      // Setup injected particle
      auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
          pcode, P0, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
      auto& stack = *stackPtr;
      auto& view = *secViewPtr;
      auto const& particle = stack.getNextParticle();
      auto const plab = MomentumVector(cs, {P0, 0_GeV, 0_GeV});

      // Decay
      decay.doDecay(*secViewPtr);
      CHECK(stack.getEntries() > 1);

      // Total energy and momentum of decay products
      MomentumVector pSum{cs, 0_eV, 0_eV, 0_eV};
      HEPEnergyType eSum{0_GeV};
      for (auto const& child : view) {
        // Check time is right
        CHECK(child.getTime() == 0_ns);

        auto P = child.getMomentum();
        auto E = child.getEnergy();

        // Record information for child particle
        tauola_pid(static_cast<int>(get_PDG(child.getPID())));

        switch (child.getPID()) {
          case Code::Electron:
            tauola_electron(E / 1_GeV);
            break;

          case Code::MuMinus:
            tauola_muon(E / 1_GeV);
            break;

          case Code::NuE:
          case Code::NuEBar:
            tauola_nue(E / 1_GeV);
            break;

          case Code::NuMu:
          case Code::NuMuBar:
            tauola_numu(E / 1_GeV);
            break;

          case Code::NuTau:
          case Code::NuTauBar:
            tauola_nutau(E / 1_GeV);
            break;

          case Code::Photon:
            tauola_photon(E / 1_GeV);
            break;

          case Code::PiMinus:
            tauola_piminus(E / 1_GeV);
            break;

          case Code::PiPlus:
            tauola_piplus(E / 1_GeV);
            break;

          case Code::Pi0:
            tauola_pizero(E / 1_GeV);
            break;

          default:
            break;
        }
        eSum += E;
        pSum += P;
      }

      // Check momentum and energy conservation
      CHECK((pSum - plab).getNorm() / 1_GeV == Approx(0).margin(1e-3));
      CHECK((pSum.getNorm() - plab.getNorm()) / 1_GeV == Approx(0).margin(1e-3));
      CHECK((particle.getEnergy() - eSum) / 1_GeV == Approx(0).margin(1e-3));
    }

    // check that no particles overflowed our histograms
    CHECK(tauola_pid.at(tauola_pid.axis().index(1e6)) == 0);

    // Record time consumption
    std::chrono::steady_clock::time_point tauola_end = std::chrono::steady_clock::now();

    // write the decay output file
    std::ofstream oftauola("tauola_decays.dat");
    for (int i = 0; i < tauola_electron.axis().size(); ++i) {
      oftauola << boost::format("%.1f %.1f %i %i %i %i %i %i %i %i %i\n") %
                      tauola_electron.axis().bin(i).lower() %
                      tauola_electron.axis().bin(i).upper() % tauola_electron.at(i) %
                      tauola_muon.at(i) % tauola_nue.at(i) % tauola_numu.at(i) %
                      tauola_nutau.at(i) % tauola_photon.at(i) % tauola_pizero.at(i) %
                      tauola_piplus.at(i) % tauola_piminus.at(i);
    }

    // write the TAUOLA PID hist
    std::ofstream oftauola_pid("tauola_pid.dat");
    for (auto&& x : indexed(tauola_pid)) {
      oftauola_pid << boost::format("%.1f %.1f %i\n") % x.bin().lower() %
                          x.bin().upper() % *x;
    }

    CORSIKA_LOG_INFO(
        "Time used by TAUOLA with {} events: {}ms", N,
        std::chrono::duration_cast<std::chrono::microseconds>(tauola_end - tauola_begin)
                .count() /
            1e3);

    // ========== PYTHIA ========
    using namespace Pythia8;
    RNGManager<>::getInstance().registerRandomStream("pythia");

    // create a pythia decay instance
    corsika::pythia8::Decay pythia;

    // check that the tau minus will be decayed by Pythia
    CHECK(pythia.isDecayHandled(Code::TauMinus));

    // make sure that Pythia can handle the tau decay
    CHECK(pythia.canHandleDecay(Code::TauMinus));

    // Create histograms
    // pid distribution
    auto pythia_pid = make_histogram(axis::regular<>(8000, -4000, 4000));
    // Histograms of energy for different particle
    auto pythia_electron = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_muon = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_nue = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_numu = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_nutau = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_photon = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_pizero = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_piplus = make_histogram(axis::regular<>(Nbins, 0.0, 1000));
    auto pythia_piminus = make_histogram(axis::regular<>(Nbins, 0.0, 1000));

    // Record time consumption
    std::chrono::steady_clock::time_point pythia_begin = std::chrono::steady_clock::now();

    for (int iEvent = 0; iEvent < N; iEvent++) {
      // Setup injected particle
      auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
          pcode, P0, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
      auto& stack = *stackPtr;
      auto& view = *secViewPtr;
      auto const& particle = stack.getNextParticle();
      auto const plab = MomentumVector(cs, {P0, 0_GeV, 0_GeV});

      // Decay by pythia
      pythia.doDecay(*secViewPtr);
      CHECK(stack.getEntries() > 1);

      // Total energy and momentum of decay products
      MomentumVector pSum{cs, 0_eV, 0_eV, 0_eV};
      HEPEnergyType eSum{0_GeV};
      for (auto const& child : view) {
        // Check time is right
        CHECK(child.getTime() == 0_ns);

        auto P = child.getMomentum();
        auto E = child.getEnergy();

        // Record information for child particle
        pythia_pid(static_cast<int>(get_PDG(child.getPID())));

        switch (child.getPID()) {
          case Code::Electron:
            pythia_electron(E / 1_GeV);
            break;

          case Code::MuMinus:
            pythia_muon(E / 1_GeV);
            break;

          case Code::NuE:
          case Code::NuEBar:
            pythia_nue(E / 1_GeV);
            break;

          case Code::NuMu:
          case Code::NuMuBar:
            pythia_numu(E / 1_GeV);
            break;

          case Code::NuTau:
          case Code::NuTauBar:
            pythia_nutau(E / 1_GeV);
            break;

          case Code::Photon:
            pythia_photon(E / 1_GeV);
            break;

          case Code::PiMinus:
            pythia_piminus(E / 1_GeV);
            break;

          case Code::PiPlus:
            pythia_piplus(E / 1_GeV);
            break;

          case Code::Pi0:
            pythia_pizero(E / 1_GeV);
            break;

          default:
            break;
        }
        eSum += E;
        pSum += P;
      }

      // Check momentum and energy conservation
      CHECK((pSum - plab).getNorm() / 1_GeV == Approx(0).margin(1e-4));
      CHECK((pSum.getNorm() - plab.getNorm()) / 1_GeV == Approx(0).margin(1e-4));
      CHECK((particle.getEnergy() - eSum) / 1_GeV == Approx(0).margin(1e-4));
    }

    // Record time consumption
    std::chrono::steady_clock::time_point pythia_end = std::chrono::steady_clock::now();

    // check that no particles overflowed our histograms
    CHECK(pythia_pid.at(pythia_pid.axis().index(1e6)) == 0);

    // write the decay output file
    std::ofstream ofpythia("pythia_decays.dat");
    for (int i = 0; i < pythia_electron.axis().size(); ++i) {
      ofpythia << boost::format("%.1f %.1f %i %i %i %i %i %i %i %i %i\n") %
                      pythia_electron.axis().bin(i).lower() %
                      pythia_electron.axis().bin(i).upper() % pythia_electron.at(i) %
                      pythia_muon.at(i) % pythia_nue.at(i) % pythia_numu.at(i) %
                      pythia_nutau.at(i) % pythia_photon.at(i) % pythia_pizero.at(i) %
                      pythia_piplus.at(i) % pythia_piminus.at(i);
    }

    // and write the Pythia PID file
    std::ofstream ofpythia_pid("pythia_pid.dat");
    for (auto&& x : indexed(pythia_pid)) {
      ofpythia_pid << boost::format("%.1f %.1f %i\n") % x.bin().lower() %
                          x.bin().upper() % *x;
    }
    CORSIKA_LOG_INFO(
        "Time used by PYTHIA with {} events: {}ms", N,
        std::chrono::duration_cast<std::chrono::microseconds>(pythia_end - pythia_begin)
                .count() /
            1e3);

    // calculate the number of microseconds per trial
    const auto tauola_per_trial{
        std::chrono::duration_cast<std::chrono::microseconds>(tauola_end - tauola_begin)
            .count() /
        N};
    const auto pythia_per_trial{
        std::chrono::duration_cast<std::chrono::microseconds>(pythia_end - pythia_begin)
            .count() /
        N};

    // print out our benchmarking information
    std::cerr << "========== START BENCHMARK INFO ==========\n";
    std::cerr << "TAUOLA: " << tauola_per_trial << "[µs/trial]\n";
    std::cerr << "PYTHIA: " << pythia_per_trial << "[µs/trial]\n";
    std::cerr << "TAUOLA/PYTHIA: " << float(tauola_per_trial) / float(pythia_per_trial)
              << "\n";
    std::cerr << "=========== END BENCHMARK INFO ===========\n";
  }

  SECTION("TAUOLA polarization") {
    //  setup environment and coordinate system
    auto [env, csPtr, nodePtr] = setup::testing::setup_environment(Code::Proton);
    [[maybe_unused]] auto const& env_dummy = env;
    auto const& cs = *csPtr;

    // Test tau- decays
    Code pid{Code::TauMinus};
    const auto mass{get_mass(pid)};

    // Energy of injection
    double normE = 100;
    HEPEnergyType const E0 = normE * 1_GeV;
    HEPEnergyType const P0 = sqrt(E0 * E0 - mass * mass);

    // Check decay of tau+ and tau-

    // Number of events
    const int N{1000};

    // Bin number for histograms
    const unsigned int Nbins{100};

    // Different states
    const auto states = {
        std::make_tuple(Helicity::RightHanded, true, "tauola_decays_rh.dat"),
        std::make_tuple(Helicity::RightHanded, false, "tauola_decays_anti_rh.dat"),
        std::make_tuple(Helicity::Unpolarized, true, "tauola_decays_nh.dat"),
        std::make_tuple(Helicity::Unpolarized, false, "tauola_decays_anti_nh.dat"),
        std::make_tuple(Helicity::LeftHanded, true, "tauola_decays_lh.dat"),
        std::make_tuple(Helicity::LeftHanded, false, "tauola_decays_anti_lh.dat")};

    for (const auto& [helicity, not_anti_tau, name] : states) {
      // Instantiate tauola decay module
      corsika::tauola::Decay decay(helicity, false, true);

      // Create histograms
      // Histograms of energy for different particle
      auto tauola_electron = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_muon = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_nue = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_numu = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_nutau = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_photon = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_pizero = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_piplus = make_histogram(axis::regular<>(Nbins, 0.0, normE));
      auto tauola_piminus = make_histogram(axis::regular<>(Nbins, 0.0, normE));

      for (int iEvent = 0; iEvent < N; iEvent++) {
        // Setup injected particle
        auto pid_tau = not_anti_tau ? Code::TauMinus : Code::TauPlus;
        auto [stackPtr, secViewPtr] = setup::testing::setup_stack(
            pid_tau, P0, (DummyEnvironment::BaseNodeType* const)nodePtr, *csPtr);
        auto& stack = *stackPtr;
        auto& view = *secViewPtr;

        auto particle = stack.first();
        auto const plab = MomentumVector(cs, {P0, 0_GeV, 0_GeV});

        // Perform decay
        decay.doDecay(*secViewPtr);

        // Check decay
        CHECK(stack.getEntries() > 1);

        // Total energy and momentum of decay products
        MomentumVector pSum{cs, 0_eV, 0_eV, 0_eV};
        HEPEnergyType eSum{0_GeV};
        for (auto const& child : view) {
          // Check time is right
          CHECK(child.getTime() == 0_ns);

          auto P = child.getMomentum();
          auto E = child.getEnergy();

          // Record information for child particle
          switch (child.getPID()) {
            case Code::Electron:
              tauola_electron(E / 1_GeV);
              break;

            case Code::MuMinus:
              tauola_muon(E / 1_GeV);
              break;

            case Code::NuE:
            case Code::NuEBar:
              tauola_nue(E / 1_GeV);
              break;

            case Code::NuMu:
            case Code::NuMuBar:
              tauola_numu(E / 1_GeV);
              break;

            case Code::NuTau:
            case Code::NuTauBar:
              tauola_nutau(E / 1_GeV);
              break;

            case Code::Photon:
              tauola_photon(E / 1_GeV);
              break;

            case Code::PiMinus:
              tauola_piminus(E / 1_GeV);
              break;

            case Code::PiPlus:
              tauola_piplus(E / 1_GeV);
              break;

            case Code::Pi0:
              tauola_pizero(E / 1_GeV);
              break;

            default:
              break;
          }
          eSum += E;
          pSum += P;
        }

        // Check momentum and energy conservation
        CHECK((pSum - plab).getNorm() / 1_GeV == Approx(0).margin(1e-2));
        CHECK((pSum.getNorm() - plab.getNorm()) / 1_GeV == Approx(0).margin(1e-2));
        CHECK((particle.getEnergy() - eSum) / 1_GeV == Approx(0).margin(1e-2));
      }

      // write the decay output file
      std::ofstream oftauola(name);
      for (int i = 0; i < tauola_electron.axis().size(); ++i) {
        oftauola << boost::format("%.1f %.1f %i %i %i %i %i %i %i %i %i\n") %
                        tauola_electron.axis().bin(i).lower() %
                        tauola_electron.axis().bin(i).upper() % tauola_electron.at(i) %
                        tauola_muon.at(i) % tauola_nue.at(i) % tauola_numu.at(i) %
                        tauola_nutau.at(i) % tauola_photon.at(i) % tauola_pizero.at(i) %
                        tauola_piplus.at(i) % tauola_piminus.at(i);
      }
    }
  }

  SECTION("TAUOLA decay config") {
    corsika::tauola::Decay decay;
    CHECK(decay.isDecayHandled(Code::WPlus));
    CHECK(decay.isDecayHandled(Code::WMinus));
    CHECK(decay.isDecayHandled(Code::Z0));
    CHECK(decay.isDecayHandled(Code::TauPlus));
    CHECK(decay.isDecayHandled(Code::TauMinus));
    CHECK(decay.isDecayHandled(Code::Eta));
    CHECK(decay.isDecayHandled(Code::K0Short));
    CHECK(decay.isDecayHandled(Code::RhoPlus));
    CHECK(decay.isDecayHandled(Code::RhoMinus));
    CHECK(decay.isDecayHandled(Code::KStarPlus));
    CHECK(decay.isDecayHandled(Code::KStarMinus));
  }
}
