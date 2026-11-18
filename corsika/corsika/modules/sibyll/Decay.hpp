/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/process/DecayProcess.hpp>

#include <set>
#include <vector>

namespace corsika::sibyll {

  class Decay : public DecayProcess<Decay> {

  public:
    Decay(const bool sibyll_listing = false);
    Decay(std::set<Code> const&);
    ~Decay();

    void printDecayConfig(const Code);
    void printDecayConfig();
    void setHadronsUnstable();

    // is Sibyll::Decay set to handle the decay of this particle?
    bool isDecayHandled(const Code);

    // is decay possible in principle?
    bool canHandleDecay(const Code);

    // set Sibyll::Decay to handle the decay of this particle!
    void setHandleDecay(const Code);
    // set Sibyll::Decay to handle the decay of this list of particles!
    void setHandleDecay(std::vector<Code> const&);
    // set Sibyll::Decay to handle all particle decays
    void setHandleAllDecay();

    template <typename TParticle>
    TimeType getLifetime(TParticle const&);

    /**
     In this function SIBYLL is called to produce to decay the input particle.
   */

    template <typename TSecondaryView>
    void doDecay(TSecondaryView&);

  private:
    // internal routines to set particles stable and unstable in the COMMON blocks in
    // sibyll
    void setStable(std::vector<Code> const&);
    void setUnstable(std::vector<Code> const&);

    void setStable(Code const);
    void setUnstable(Code const);

    // internally set all particles to decay/not to decay
    void setAllUnstable();
    void setAllStable();

    // will this particle be stable in sibyll ?
    bool isStable(Code const);
    // will this particle decay in sibyll ?
    bool isUnstable(Code const);
    // set particle with input code to decay or not
    void setDecay(Code const, bool const);

    // data members
    int count_ = 0;
    bool handleAllDecays_ = true;
    bool sibyll_listing_ = false;
    std::set<Code> handledDecays_;

    std::shared_ptr<spdlog::logger> logger_ = get_logger("corsika_sibyll_Decay");
  };

} // namespace corsika::sibyll

#include <corsika/detail/modules/sibyll/Decay.inl>
