/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/modules/radio/RadioProcess.hpp>

namespace corsika {

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline TRadioImpl&
  RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::implementation() {
    return static_cast<TRadioImpl&>(*this);
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline TRadioImpl const&
  RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::implementation() const {
    return static_cast<TRadioImpl const&>(*this);
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::RadioProcess(
      TObserverCollection& observers, TPropagator& propagator)
      : observers_(observers)
      , propagator_(propagator) {}

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  template <typename Particle>
  inline ProcessReturn RadioProcess<TObserverCollection, TRadioImpl,
                                    TPropagator>::doContinuous(const Step<Particle>& step,
                                                               const bool) {

    // return immediately if radio process does not have any observers
    if (observers_.size() == 0) return ProcessReturn::Ok;

    // we want the following particles:
    // Code::Electron & Code::Positron

    // we wrap Simulate() in doContinuous as the plan is to add particle level
    // filtering or thinning for calculation of the radio emission. This is
    // important for controlling the runtime of radio (by ignoring particles
    // that aren't going to contribute i.e. heavy hadrons)
    // if (valid(step)) {
    auto const particleID_{step.getParticlePre().getPID()};
    if ((particleID_ == Code::Electron) || (particleID_ == Code::Positron)) {
      CORSIKA_LOG_DEBUG("Particle for radio calculation: {} ", particleID_);
      return this->implementation().simulate(step);
    } else {
      CORSIKA_LOG_DEBUG("Particle {} is irrelevant for radio", particleID_);
      return ProcessReturn::Ok;
    }
    //}
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  template <typename Particle, typename Track>
  inline LengthType
  RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::getMaxStepLength(
      [[maybe_unused]] const Particle& vParticle,
      [[maybe_unused]] const Track& vTrack) const {
    return meter * std::numeric_limits<double>::infinity();
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline void RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::startOfLibrary(
      const boost::filesystem::path& directory) {

    // setup the streamer
    output_.initStreamer((directory / ("observers.parquet")).string());

    // enable compression with the default level
    output_.enableCompression();

    // LCOV_EXCL_START
    // build the schema
    output_.addField("Time", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);

    output_.addField("Ex", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);

    output_.addField("Ey", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);

    output_.addField("Ez", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);
    // LCOV_EXCL_STOP
    // and build the streamer
    output_.buildStreamer();
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline void RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::endOfShower(
      const unsigned int) {

    // loop over every observer and instruct them to
    // flush data to disk, and then reset the observer
    // before the next event
    for (auto& observer : observers_.getObservers()) {

      auto const sampleRate = observer.getSampleRate() * 1_s;
      auto const radioImplementation =
          static_cast<std::string>(this->implementation().algorithm);

      // get the axis labels for this observer and write the first row.
      axistype axis = observer.implementation().getAxis();

      // get the copy of the waveform data for this event
      std::vector<double> const& dataX = observer.implementation().getWaveformX();
      std::vector<double> const& dataY = observer.implementation().getWaveformY();
      std::vector<double> const& dataZ = observer.implementation().getWaveformZ();

      // check for the axis name
      std::string label = "Unknown";
      if (observer.getDomainLabel() == "Time") {
        label = "Time";
      }
      // LCOV_EXCL_START
      else if (observer.getDomainLabel() == "Frequency") {
        label = "Frequency";
      }
      // LCOV_EXCL_STOP
      if (radioImplementation == "ZHS" && label == "Time") {
        for (size_t i = 0; i < axis.size() - 1; i++) {
          auto time = (axis.at(i + 1) + axis.at(i)) / 2.;
          auto Ex = -(dataX.at(i + 1) - dataX.at(i)) * sampleRate;
          auto Ey = -(dataY.at(i + 1) - dataY.at(i)) * sampleRate;
          auto Ez = -(dataZ.at(i + 1) - dataZ.at(i)) * sampleRate;

          *(output_.getWriter())
              << showerId_ << static_cast<double>(time) << static_cast<double>(Ex)
              << static_cast<double>(Ey) << static_cast<double>(Ez) << parquet::EndRow;
        }
      } else if (radioImplementation == "CoREAS" && label == "Time") {
        for (size_t i = 0; i < axis.size() - 1; i++) {
          *(output_.getWriter())
              << showerId_ << static_cast<double>(axis[i])
              << static_cast<double>(dataX[i]) << static_cast<double>(dataY[i])
              << static_cast<double>(dataZ[i]) << parquet::EndRow;
        }
      }

      observer.reset();
    }

    // increment our event counter
    showerId_++;
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline void RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::endOfLibrary() {
    output_.closeStreamer();
  }

  template <typename TObserverCollection, typename TRadioImpl, typename TPropagator>
  inline YAML::Node
  RadioProcess<TObserverCollection, TRadioImpl, TPropagator>::getConfig() const {

    // top-level YAML node
    YAML::Node config;

    // fill in some basics
    config["type"] = "RadioProcess";
    config["algorithm"] = this->implementation().algorithm;
    config["units"]["time"] = "ns";
    config["units"]["frequency"] = "GHz";
    config["units"]["electric field"] = "V/m";
    config["units"]["distance"] = "m";

    for (auto& observer : observers_.getObservers()) {
      // get the name/location of this observer
      auto name = observer.getName();
      auto location = observer.getLocation().getCoordinates();

      // get the observers config
      config["observers"][name] = observer.getConfig();

      // write the location of this observer
      config["observers"][name]["location"].push_back(location.getX() / 1_m);
      config["observers"][name]["location"].push_back(location.getY() / 1_m);
      config["observers"][name]["location"].push_back(location.getZ() / 1_m);
    }

    return config;
  }

} // namespace corsika
