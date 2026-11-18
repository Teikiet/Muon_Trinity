/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/modules/radio/RadioProcess.hpp>

namespace corsika {
  /*!
   * Wrapper for a radio process which turns the production of signal on and off depending
   * on where the currently-emitting particle is in space. A function with signature `bool
   * Fcn(Point const&)` must be given which turns this on and off
   */

  template <class TUnderlyingProcess>
  class LimitedRadioProcess
      : public ContinuousProcess<LimitedRadioProcess<TUnderlyingProcess>> {

    //! define the default constructor for the function signature
    static auto constexpr default_functor = [](Point const&) { return true; };

  public:
    //! The signature of the if-statement that decides if this RadioProcess should be
    //! called
    using functor_signature = bool(Point const&);

    LimitedRadioProcess(TUnderlyingProcess&& process,
                        std::function<functor_signature> runThis);

    // This function wraps the one from RadioProcess and runs it depending on output of
    // `runThis_`
    template <typename Particle>
    ProcessReturn doContinuous(Step<Particle> const& step, bool const);

    // This process just returns the one from RadioProcess since it would be more work to
    // turn this on and off than to just feed-through
    template <typename Particle, typename Track>
    LengthType getMaxStepLength(Particle const& vParticle, Track const& vTrack) const;

  private:
    TUnderlyingProcess process_;
    std::function<functor_signature> const runThis_{default_functor};
  };
} // namespace corsika

#include <corsika/detail/modules/radio/LimitedRadioProcess.inl>
