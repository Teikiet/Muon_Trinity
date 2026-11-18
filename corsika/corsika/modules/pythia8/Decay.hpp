/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/process/DecayProcess.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/modules/pythia8/Pythia8.hpp>

namespace corsika::pythia8 {

  class Decay : public DecayProcess<Decay>, public Pythia8::Pythia {

  public:
    Decay(bool const print_listing = false);
    Decay(std::set<Code> const&);
    ~Decay();

    // is Pythia::Decay set to handle the decay of this particle?
    bool isDecayHandled(Code const);

    //! is decay possible in principle?
    bool canHandleDecay(Code const);

    //! set Pythia::Decay to handle the decay of this particle!
    void setHandleDecay(Code const);
    //! set Pythia::Decay to handle the decay of this list of particles!
    void setHandleDecay(std::vector<Code> const&);
    //! set Pythia::Decay to handle all particle decays
    void setHandleAllDecays();

    //! print internal configuration for this particle
    void printDecayConfig(Code const);
    //! print configuration of decays in corsika
    void printDecayConfig();

    bool canDecay(Code const);

    template <typename TParticle>
    TimeType getLifetime(TParticle const&);

    template <typename TView>
    void doDecay(TView&);

  private:
    void init();

    bool isStable(Code const vCode);
    void setStable(std::vector<Code> const&);
    void setUnstable(Code const);
    void setStable(Code const);

    // data members
    int count_ = 0;
    bool handleAllDecays_ = true;
    std::set<Code> handledDecays_;
    bool print_listing_ = false;
  };

} // namespace corsika::pythia8

#include <corsika/detail/modules/pythia8/Decay.inl>
