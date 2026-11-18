/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/output/BaseOutput.hpp>
#include <corsika/output/ParquetStreamer.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

namespace corsika {

  /**
   * The base interface for radio emission processes.
   *
   * TRadioImpl is the concrete implementation of the radio algorithm.
   * TObserverCollection is the detector instance that stores observers
   * and is responsible for managing the output writing.
   */
  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  class RadioProcess : public ContinuousProcess<
                           RadioProcess<TObserverCollection, TRadioImpl, TPropagator>>,
                       public BaseOutput {

    /*
     * A collection of filter objects for deciding on valid particles and tracks.
     */
    // std::vector<std::function<bool(ParticleType&, TrackType const&)>> filters_;

    /**
     * Get a reference to the underlying radio implementation.
     */
    TRadioImpl& implementation();

    /**
     *  Get a const reference to the underlying implementation.
     */
    TRadioImpl const& implementation() const;

  protected:
    TObserverCollection& observers_; ///< The radio observers we store into.
    TPropagator propagator_;         ///< The propagator implementation.
    unsigned int showerId_{0};       ///< The current event ID.
    ParquetStreamer output_;         //!< The parquet streamer for this process.

  public:
    using axistype = std::vector<long double>;
    /**
     * Construct a new RadioProcess.
     */
    RadioProcess(TObserverCollection& observers, TPropagator& propagator);

    /**
     * Perform the continuous process (radio emission).
     *
     * This handles filtering individual particle tracks
     * before passing them to `Simulate`.`
     *
     * @param particle    The current particle.
     * @param track       The current track.
     */
    template <typename Particle>
    ProcessReturn doContinuous(Step<Particle> const& step, bool const);

    /**
     * Return the maximum step length for this particle and track.
     *
     * This must be provided by the TRadioImpl.
     *
     * @param particle    The current particle.
     * @param track       The current track.
     *
     * @returns The maximum length of this track.
     */
    template <typename Particle, typename Track>
    LengthType getMaxStepLength(Particle const& vParticle, Track const& vTrack) const;

    /**
     * Called at the start of each library.
     */
    void startOfLibrary(boost::filesystem::path const& directory) final override;

    /**
     * Called at the end of each shower.
     */
    virtual void endOfShower(unsigned int const) final override;

    /**
     * Called at the end of each library.
     *
     */
    void endOfLibrary() final override;

    /**
     * Get the configuration of this output.
     */
    YAML::Node getConfig() const final;

  }; // END: class RadioProcess

} // namespace corsika

#include <corsika/detail/modules/radio/RadioProcess.inl>
