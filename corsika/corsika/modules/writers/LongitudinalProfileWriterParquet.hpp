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

#include <array>

namespace corsika {

  /**
   * The actual writer to save longitudinal profile data to disk.
   *
   * The purpose of this class is not to collect single-particle-level profile data,
   * but to write entire binned profiles to disk at the end of a shower.
   * To fill the shower data, you have to use ProfileWriter in combination with
   * SubWriter. The LongitudinalProfileWriterParquet is the default output mode of the
   * LongitudinalWriter.
   */

  template <size_t NColumns>
  class LongitudinalProfileWriterParquet : public BaseOutput {

  public:
    /**
     * Construct a new writer.
     */
    LongitudinalProfileWriterParquet(std::array<const char*, NColumns> const& colNames);

    /**
     * Called at the start of each library.
     */
    void startOfLibrary(boost::filesystem::path const& directory) override;

    /**
     * Called at the start of each shower.
     */
    void startOfShower(unsigned int const showerId) override;

    /**
     * Called at the end of each shower.
     */
    void endOfShower(unsigned int const showerId) override;

    /**
     * Called at the end of each library.
     *
     * This must also increment the run number since we override
     * the default behaviour of BaseOutput.
     */
    void endOfLibrary() override;

    /**
     * Add profile to disk.
     */
    void write(unsigned int const showerId, GrammageType const grammage,
               std::array<double, NColumns> const& data);

  private:
    std::array<const char*, NColumns> columns_; //!< column names
    ParquetStreamer output_; //!< The parquet streamer for this process.

  }; // namespace corsika

} // namespace corsika

#include <corsika/detail/modules/writers/LongitudinalProfileWriterParquet.inl>
