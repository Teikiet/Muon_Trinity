/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <Tauola/Tauola.h>

#include <corsika/framework/utility/COMBoost.hpp>

namespace Tpp = Tauolapp;

namespace corsika::tauola {

  // default constructor
  Decay::Decay(const Helicity helicity, const bool radiative, const bool spincorr)
      : helicity_(helicity) {

    // if TAUOLA has not already been initialized,
    if (!Tpp::Tauola::getIsTauolaIni()) {

      // set the default units explicitly.
      // we use GeV from momentum and cm for length
      Tpp::Tauola::setUnits(Tpp::Tauola::GEV, Tpp::Tauola::CM);

      // we set the lifetime as zero so that the tau decays
      // at the vertex point given by CORSIKA since
      // CORSIKA handles proapgating the tau until the decay point.
      Tpp::Tauola::setTauLifetime(0.0);

      // don't decay any intermediaries - add these to the stack
      // in particular, eta's, K0, Pi-0's.
      // See Pg. 46 in the TAUOLA manual.
      Tpp::Tauola::setEtaK0sPi(1, 1, 0);

      // use the new currents (Pg. 40 in the user manual).
      Tpp::Tauola::setNewCurrents(1);

      // whether to use QCD radiative corrections in calculating
      // leptonic tau decay. If this is enabled, TAUOLA generates
      // an order of magnitude more photons than Pythia (for example)
      Tpp::Tauola::setRadiation(radiative);

      // we have to set the spin correlation state after
      // we have initialized TAUOLA
      Tpp::Tauola::spin_correlation.setAll(spincorr);

      // initialize the TAUOLA library
      Tpp::Tauola::initialize();

      // connect the C8 RNG manager to the TAUOLA RNG manager
      corsika::connect_random_stream(
          RNGManager<>::getInstance().getRandomStream("tauola"),
          [](corsika::rng_function_type f) { Decay::rngFcn_ = std::move(f); });
      Tpp::Tauola::setRandomGenerator(Decay::WrappedRNG);
    }
  }

  Decay::~Decay() {
    CORSIKA_LOGGER_INFO(logger_, "Number of particles decayed using TAUOLA: {}", count_);
  }

  bool Decay::isDecayHandled(const Code code) {
    // if known to tauola, it can decay
    if (code == Code::WPlus || code == Code::WMinus || code == Code::Z0 ||
        code == Code::TauPlus || code == Code::TauMinus || code == Code::Eta ||
        code == Code::K0Short || code == Code::RhoPlus || code == Code::RhoMinus ||
        code == Code::KStarPlus || code == Code::KStarMinus)
      return true;
    else
      return false;
  }

  // perform the decay
  template <typename TSecondaryView>
  void Decay::doDecay(TSecondaryView& view) {
    CORSIKA_LOGGER_DEBUG(logger_, "Decaying primary particle using TAUOLA...");
    count_++;

    auto projectile = view.getProjectile();

    // get the TAUOLA particle from this projectile
    auto particle{TauolaInterfaceParticle::from_projectile(projectile)};

    CORSIKA_LOGGER_TRACE(
        logger_, "[decay primary] code: {}, momentum: ({}, {}, {}) GeV, energy: {} GeV",
        convert_from_PDG(static_cast<PDGCode>(particle->getPdgID())), particle->getPx(),
        particle->getPy(), particle->getPz(), particle->getE());

    // decay this particle using TAUOLA
    decay(particle.get());

    std::vector<TauolaInterfaceParticle*> daughters;

    // we iterate over the daughter particles and check for
    // short-lived secondaries that we have to decay again using TAUOLA.
    // We only ever have to go one layer down in the hierarchy.
    for (const auto& daughter : particle->getDaughters()) {

      // check if this is decay can be handled by TAUOLA
      if (abs(daughter->getPdgID()) == 20213 ||
          isDecayHandled(convert_from_PDG(static_cast<PDGCode>(daughter->getPdgID())))) {
        // decay this secondary
        decay(daughter);
        // and now add all the daughters of this particle to our list
        for (const auto& child : daughter->getDaughters()) {

          // check if this is further decayable (K0 and Eta's typically)
          if (abs(child->getPdgID()) == 20213 ||
              isDecayHandled(convert_from_PDG(static_cast<PDGCode>(child->getPdgID())))) {
            // decay the cild
            decay(child);

            // and add it's daughters to the top-level daughters vector
            for (const auto& grandchild : child->getDaughters()) {
              daughters.push_back(dynamic_cast<TauolaInterfaceParticle*>(grandchild));
            }
          } else {
            // if we are here, the original daughter is non-decayable so add it to our
            // struct.
            daughters.push_back(dynamic_cast<TauolaInterfaceParticle*>(child));
          }
        }

      } else {
        // if we get here, `daughter` is long-lived and doesn't require a decay.
        daughters.push_back(dynamic_cast<TauolaInterfaceParticle*>(daughter));
      }
    }

    // we now have a complete list of the daughter particles of this decay
    // we now iterate over them and add them to CORSIKA stack.
    for (const auto& daughter : daughters) {

      // get particle's PDG ID
      const auto pdg{daughter->getPdgID()};

      // get the CORSIKA PID from our PDG ID
      const auto pid{convert_from_PDG(static_cast<PDGCode>(pdg))};

      // construct the momentum in the ROOT coordinate system
      const MomentumVector P(daughter->getPx() * 1_GeV, daughter->getPy() * 1_GeV,
                             daughter->getPz() * 1_GeV);

      HEPEnergyType const mass = get_mass(pid);
      HEPEnergyType const Ekin = sqrt(P.getSquaredNorm() + mass * mass) - mass;

      CORSIKA_LOGGER_TRACE(logger_,
                           "[decay result] code: {}, momentum: {} GeV, energy: {} GeV",
                           pid, P.getComponents() / 1_GeV, daughter->getE());

      // and add the particle to the CORSIKA stack
      projectile.addSecondary(std::make_tuple(pid, Ekin, P.normalized()));
    }
  }

  template <typename InterfaceParticle>
  void Decay::decay(InterfaceParticle* particle) {

    // return if we get a nullptr here
    // this should never happen (but to be safe!)
    // LCOV_EXCL_START
    if (particle == nullptr) {
      CORSIKA_LOGGER_WARN(logger_, "Got nullptr in decay routine. Skipping decay.");
      return;
    }
    // LCOV_EXCL_STOP

    // get the particle's PDG
    const auto pdg{static_cast<PDGCode>(particle->getPdgID())};

    // if this is a tau- or tau+, use the helicity in the decay
    if (pdg == PDGCode::TauMinus || pdg == PDGCode::TauPlus) {

      // get helicity from enum class
      double helic;
      switch (helicity_) {
        case Helicity::LeftHanded:
          helic = -1.;
          break;
        case Helicity::RightHanded:
          helic = 1.;
          break;
        case Helicity::Unpolarized:
          helic = 0.;
          break;
          // NOT testing undefined
        default: // LCOV_EXCL_START
          CORSIKA_LOGGER_ERROR(logger_, "Unknown helicity setting: {}", helicity_);
          throw std::runtime_error("Unknown helicity setting!");
          // LCOV_EXCL_STOP
      }

      double Polx{0};
      double Poly{0};
      double Polz{-helic};

      // and call the decay routine with the polarization
      Tpp::Tauola::decayOne(particle, false, Polx, Poly, Polz);

    } else { // this is not a tau, so just do a standard decay
      Tpp::Tauola::decayOne(particle);
    }
  }

  template <typename TParticle>
  TimeType Decay::getLifetime(TParticle const& particle) {
    // Can be moved into Cascade.inl, see Issue #271
    auto const pid = particle.getPID();
    if ((pid == Code::TauMinus) || (pid == Code::TauPlus)) {
      HEPEnergyType E = particle.getEnergy();
      HEPMassType m = particle.getMass();
      double const gamma = E / m;
      TimeType const t0 = get_lifetime(pid);
      auto const lifetime = gamma * t0;
      CORSIKA_LOGGER_TRACE(logger_,
                           "[decay time] code: {}, t0: {} ns, momentum: {} GeV, energy: "
                           "{} GeV, gamma: {}, tau: {} ns",
                           particle.getPID(), t0 / 1_ns,
                           particle.getMomentum().getComponents() / 1_GeV, E / 1_GeV,
                           gamma, lifetime / 1_ns);
      return lifetime;
    } else
      return std::numeric_limits<double>::infinity() * 1_s;
  }

  void Decay::PrintSummary() const { Tpp::Tauola::summary(); }

  inline double Decay::WrappedRNG() {
    double val;
    rngFcn_(&val, 1); // Fill one value
    return val;
  }

} // namespace corsika::tauola
