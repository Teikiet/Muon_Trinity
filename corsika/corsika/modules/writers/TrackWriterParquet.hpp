/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/output/BaseOutput.hpp>
#include <corsika/output/ParquetStreamer.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>

namespace corsika {

  class TrackWriterParquet : public BaseOutput {

  public:
    /**
     * Construct a new writer.
     *
     * @param name    The name of this output.
     */
    TrackWriterParquet();

    /**
     * Called at the start of each library.
     */
    void startOfLibrary(boost::filesystem::path const& directory) final override;

    /**
     * Called at the start of each shower.
     */
    void startOfShower(unsigned int const showerId) final override;

    /**
     * Called at the end of each shower.
     */
    void endOfShower(unsigned int const showerId) final override;

    /**
     * Called at the end of each library.
     *
     * This must also increment the run number since we override
     * the default behaviour of BaseOutput.
     */
    void endOfLibrary() final override;

    /**
     * Write a track to the file.
     */
    void write(Code const pid, HEPEnergyType const KinenergyPre, double const weight,
               QuantityVector<length_d> const& start, TimeType const t_start,
               QuantityVector<length_d> const& end, HEPEnergyType const KinenergyPost,
               TimeType const t_end);

  private:
    ParquetStreamer output_; ///< The primary output file.
    unsigned int showerId_;  ///< event Id counter

  }; // class TrackWriterParquet

} // namespace corsika

#include <corsika/detail/modules/writers/TrackWriterParquet.inl>
