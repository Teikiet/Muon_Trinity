/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

namespace corsika {

  template <class TUnderlyingProcess>
  inline LimitedRadioProcess<TUnderlyingProcess>::LimitedRadioProcess(
      TUnderlyingProcess&& process,
      std::function<LimitedRadioProcess::functor_signature> runThis)
      : process_{std::forward<TUnderlyingProcess>(process)}
      , runThis_{std::move(runThis)} {}

  template <class TUnderlyingProcess>
  template <typename Particle>
  inline ProcessReturn LimitedRadioProcess<TUnderlyingProcess>::doContinuous(
      const Step<Particle>& step, const bool arg) {

    if ((step.getParticlePre().getPID() != Code::Electron) &&
        (step.getParticlePre().getPID() != Code::Positron)) {
      CORSIKA_LOG_DEBUG("Not e+/-, will not run RadioProcess");
      return ProcessReturn::Ok;
    }

    // Check for running the wrapped RadioProcess function
    if (runThis_(step.getPositionPre())) {
      CORSIKA_LOG_TRACE("Will run RadioProcess for particle at {}",
                        step.getPositionPre());
      return process_.doContinuous(step, arg);
    }
    CORSIKA_LOG_DEBUG("Particle at {} is not selected by LimitedRadioProcess",
                      step.getPositionPre());
    return ProcessReturn::Ok;
  }

  template <class TUnderlyingProcess>
  template <typename Particle, typename Track>
  inline LengthType LimitedRadioProcess<TUnderlyingProcess>::getMaxStepLength(
      Particle const& vParticle, Track const& vTrack) const {
    return process_.getMaxStepLength(vParticle, vTrack);
  }

} // namespace corsika