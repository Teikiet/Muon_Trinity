/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <limits>

namespace corsika {

  template <typename TOutput>
  inline TrackWriter<TOutput>::TrackWriter() {}

  template <typename TOutput>
  template <typename TParticle>
  inline ProcessReturn TrackWriter<TOutput>::doContinuous(Step<TParticle> const& step,
                                                          bool const) {

    auto const start = step.getPositionPre().getCoordinates();
    auto const end = step.getPositionPost().getCoordinates();

    // write the track to the file
    TOutput::write(step.getParticlePre().getPID(), step.getEkinPre(),
                   step.getParticlePre().getWeight(), start, step.getTimePre(), end,
                   step.getEkinPost(), step.getTimePost());

    return ProcessReturn::Ok;
  }

  template <typename TOutput>
  template <typename TParticle, typename TTrack>
  inline LengthType TrackWriter<TOutput>::getMaxStepLength(TParticle const&,
                                                           TTrack const&) {
    return meter * std::numeric_limits<double>::infinity();
  }

  template <typename TOutput>
  YAML::Node TrackWriter<TOutput>::getConfig() const {
    using namespace units::si;

    YAML::Node node;

    // add default units for values
    node["type"] = "TrackWriter";
    node["units"] = "GeV | m | ns";

    return node;
  }
} // namespace corsika
