/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/modules/pythia8/Pythia8.hpp>
#include <corsika/modules/pythia8/Decay.hpp>
#include <corsika/modules/pythia8/Random.hpp>

#include <corsika/framework/utility/COMBoost.hpp>

namespace corsika::pythia8 {

  inline Decay::Decay(bool const print_listing)
      : Pythia8::Pythia(CORSIKA_Pythia8_XML_DIR)
      , print_listing_(print_listing) {
    init();
  }

  inline Decay::Decay(std::set<Code> const& those)
      : handleAllDecays_(false)
      , handledDecays_(those) {
    init();
  }

  inline Decay::~Decay() { CORSIKA_LOG_INFO("Pythia::Decay n={}", count_); }

  inline void Decay::init() {

    // run this only once during construction

    // link random number generator in pythia to CORSIKA8
    Pythia8::Pythia::setRndmEnginePtr(std::make_shared<corsika::pythia8::Random>());

    Pythia8::Pythia::readString("Next:numberShowInfo = 0");
    Pythia8::Pythia::readString("Next:numberShowProcess = 0");
    Pythia8::Pythia::readString("Next:numberShowEvent = 0");

    Pythia8::Pythia::readString("Print:quiet = on");
    Pythia8::Pythia::readString("Check:particleData = off");

    /*
       switching off event check in pythia is needed to allow decays that are off-shell
       according to the mass definition in pythia.
       the consistency of particle masses between event generators is an unsolved issues
    */
    CORSIKA_LOG_INFO("Pythia::Init: switching off event checking in pythia..");
    Pythia8::Pythia::readString("Check:event = 1");

    Pythia8::Pythia::readString("ProcessLevel:all = off");
    Pythia8::Pythia::readString("ProcessLevel:resonanceDecays = on");

    // making sure
    setStable(Code::Pi0);

    //    Pythia8::Pythia::particleData.readString("59:m0 = 101.00");

    // LCOV_EXCL_START, we don't validate pythia8 internals
    if (!Pythia8::Pythia::init())
      throw std::runtime_error("Pythia::Decay: Initialization failed!");
    // LCOV_EXCL_STOP
  }

  inline bool Decay::canHandleDecay(Code const vParticleCode) {
    // if known to pythia and not proton, electron or neutrino it can decay
    if (vParticleCode == Code::Proton || vParticleCode == Code::AntiProton ||
        vParticleCode == Code::NuE || vParticleCode == Code::NuMu ||
        vParticleCode == Code::NuTau || vParticleCode == Code::NuEBar ||
        vParticleCode == Code::NuMuBar || vParticleCode == Code::NuTauBar ||
        vParticleCode == Code::Electron || vParticleCode == Code::Positron)
      return false;
    else if (canDecay(vParticleCode)) // check pythia8 internal
      return true;
    else
      return false;
  }

  inline void Decay::setHandleDecay(Code const vParticleCode) {
    handleAllDecays_ = false;
    CORSIKA_LOG_DEBUG("Pythia::Decay: set to handle decay of {} ", vParticleCode);
    if (Decay::canHandleDecay(vParticleCode))
      handledDecays_.insert(vParticleCode);
    else
      throw std::runtime_error("this decay can not be handled by pythia!");
  }

  inline void Decay::setHandleDecay(std::vector<Code> const& vParticleList) {
    handleAllDecays_ = false;
    for (auto const p : vParticleList) setHandleDecay(p);
  }

  inline bool Decay::isDecayHandled(Code const vParticleCode) {
    if (handleAllDecays_ && canHandleDecay(vParticleCode))
      return true;
    else
      return handledDecays_.find(vParticleCode) != Decay::handledDecays_.end();
  }

  inline void Decay::setStable(std::vector<Code> const& particleList) {
    for (auto const p : particleList) Decay::setStable(p);
  }

  inline void Decay::setUnstable(Code const pCode) {
    CORSIKA_LOG_DEBUG("Pythia::Decay: setting {} unstable..", pCode);
    Pythia8::Pythia::particleData.mayDecay(static_cast<int>(get_PDG(pCode)), true);
  }

  inline void Decay::setStable(Code const pCode) {
    CORSIKA_LOG_DEBUG("Pythia::Decay: setting {} stable..", pCode);
    Pythia8::Pythia::particleData.mayDecay(static_cast<int>(get_PDG(pCode)), false);
  }

  inline bool Decay::isStable(Code const vCode) {
    return Pythia8::Pythia::particleData.canDecay(static_cast<int>(get_PDG(vCode)));
  }

  inline bool Decay::canDecay(Code const pCode) {
    bool const ans =
        Pythia8::Pythia::particleData.canDecay(static_cast<int>(get_PDG(pCode)));
    CORSIKA_LOG_DEBUG("Pythia::Decay: checking if particle: {} can decay in PYTHIA? {} ",
                      pCode, ans);
    return ans;
  }

  inline void Decay::printDecayConfig(Code const vCode) {
    CORSIKA_LOG_INFO("Decay: Pythia decay configuration:");
    CORSIKA_LOG_INFO(" {} is {} ", vCode, (isStable(vCode) ? "stable" : "unstable"));
  }

  inline void Decay::printDecayConfig() {
    CORSIKA_LOG_INFO("Pythia::Decay: decay configuration:");
    if (handleAllDecays_)
      CORSIKA_LOG_INFO(" all particles known to Pythia are handled by Pythia::Decay!");
    else
      for (auto const pCode : handledDecays_)
        CORSIKA_LOG_INFO("Decay of {} is handled by Pythia!", pCode);
  }

  template <typename TParticle>
  inline TimeType Decay::getLifetime(TParticle const& particle) {

    auto const pid = particle.getPID();
    if (canDecay(pid)) {
      HEPEnergyType E = particle.getEnergy();
      HEPMassType m = particle.getMass();

      double const gamma = E / m;

      TimeType const t0 = get_lifetime(pid);
      auto const lifetime = gamma * t0;
      CORSIKA_LOG_TRACE("Pythia::Decay: code: {}", particle.getPID());
      CORSIKA_LOG_TRACE("Pythia::Decay: MinStep: t0: {}", t0);
      CORSIKA_LOG_TRACE("Pythia::Decay: MinStep: energy: {} GeV", E / 1_GeV);
      CORSIKA_LOG_TRACE("Pythia::Decay: momentum: {} GeV",
                        particle.getMomentum().getComponents() / 1_GeV);
      CORSIKA_LOG_TRACE("Pythia::Decay: MinStep: gamma: {}", gamma);
      CORSIKA_LOG_TRACE("Pythia::Decay: MinStep: tau: {} ", lifetime);

      return lifetime;
    } else
      return std::numeric_limits<double>::infinity() * 1_s;
  }

  template <typename TView>
  inline void Decay::doDecay(TView& view) {

    auto projectile = view.getProjectile();

    auto const& labMomentum = projectile.getMomentum();

    // define target kinematics in lab frame
    // define boost to and from CoM frame
    // CoM frame definition in Pythia projectile: +z
    COMBoost const boost(labMomentum, projectile.getMass());

    count_++;

    // pythia stack
    Pythia8::Event& event = Pythia8::Pythia::event;
    event.reset();

    auto const particleId = projectile.getPID();

    // set particle unstable
    Decay::setUnstable(particleId);

    // input particle PDG
    auto const pdgCode = static_cast<int>(get_PDG(particleId));

    double constexpr px = 0;
    double constexpr py = 0;
    double constexpr pz = 0;
    double const en = projectile.getMass() / 1_GeV;
    double const m = en;

    // add particle to pythia stack
    event.append(pdgCode, 1, 0, 0, // PID, status, col, acol
                 px, py, pz, en, m);

    // LCOV_EXCL_START, we don't validate pythia8 internals
    if (!Pythia8::Pythia::next())
      throw std::runtime_error("Pythia::Decay: decay failed!");
    // LCOV_EXCL_STOP

    CORSIKA_LOG_DEBUG("Pythia::Decay: particles after decay: {} ", event.size());

    // LCOV_EXCL_START, we don't validate pythia8 internals
    if (print_listing_) {
      // list final state
      event.list();
    }
    // LCOV_EXCL_STOP

    // loop over final state
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal()) {
        try {
          auto const id = static_cast<PDGCode>(event[i].id());
          auto const pyId = convert_from_PDG(id);
          HEPEnergyType const Erest = event[i].e() * 1_GeV;
          MomentumVector const pRest(boost.makeCoMVector(QuantityVector<hepmomentum_d>{
              event[i].px() * 1_GeV, event[i].py() * 1_GeV, event[i].pz() * 1_GeV}));
          FourVector const fourMomRest{Erest, pRest};
          auto const fourMomLab = boost.fromCoM(fourMomRest);
          auto const p3 = fourMomLab.getSpaceLikeComponents();

          HEPEnergyType const mass = get_mass(pyId);
          HEPEnergyType const Ekin = sqrt(p3.getSquaredNorm() + mass * mass) - mass;

          CORSIKA_LOG_TRACE("particle: id={} momentum={} energy={} ", pyId,
                            fourMomLab.getSpaceLikeComponents().getComponents() / 1_GeV,
                            fourMomLab.getTimeLikeComponent());

          view.addSecondary(std::make_tuple(pyId, Ekin, p3.normalized()));
        } catch (std::out_of_range const& ex) {
          CORSIKA_LOG_CRITICAL("Pythia ID {} unknown in C8", event[i].id());
          throw ex;
        }
      }
    }

    // set particle stable
    Decay::setStable(particleId);
  }
} // namespace corsika::pythia8
