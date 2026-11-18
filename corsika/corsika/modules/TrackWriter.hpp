/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/modules/writers/TrackWriterParquet.hpp>
#include <corsika/modules/writers/WriterOff.hpp>
#include <corsika/framework/core/Step.hpp>

namespace corsika {

  /**
   * @ingroup Modules
   * @{
   *
   * To write 3D track data to disk.
   *
   * Since the only sole purpose of this module is to generate track
   * output on disk, the default output mode is not "WriterOff" but
   * directly TrackWriterParquet. It can of course be changed.
   *
   * @tparam TOutput with default TrackWriterParquet
   */

  template <typename TOutput = TrackWriterParquet>
  class TrackWriter : public ContinuousProcess<TrackWriter<TOutput>>, public TOutput {

  public:
    TrackWriter();

    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle> const&, bool const limitFlag);

    template <typename TParticle, typename TTrack>
    LengthType getMaxStepLength(TParticle const&, TTrack const&);

    YAML::Node getConfig() const;

  private:
  };

  //! @}

} // namespace corsika

#include <corsika/detail/modules/TrackWriter.inl>
