/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/InteractionHistogram.hpp>
#include <corsika/modules/writers/InteractionHistogramWriterParquet.hpp>
#include <corsika/output/BaseOutput.hpp>
#include <corsika/output/YAMLStreamer.hpp>

#include <boost/filesystem.hpp>
#include <functional>
#include <yaml-cpp/yaml.h>

namespace corsika {

  /**
   * Registers combined InteractionHistogram output with OutputManager.
   *
   * Parquet layout is implemented by InteractionHistogramWriterParquet: separate files
   * for laboratory-frame energy and for \f$\sqrt{s}\f$, each a wide table (one row per
   * energy bin, one FLOAT column per projectile PDG).
   *
   * For per-shower histograms, pass @a resetCountersForShower to clear the underlying
   * InteractionCounter histograms at each startOfShower (before the cascade runs).
   */
  class InteractionHistogramWriter : public BaseOutput {

  public:
    /**
     * @param getCombined            Called at each endOfShower; return value is written.
     * @param resetCountersForShower If set, called at each startOfShower before tracking
     *                               (use InteractionCounter::resetHistogram()).
     */
    explicit InteractionHistogramWriter(
        std::function<InteractionHistogram()> getCombined,
        std::function<void()> resetCountersForShower = {});

    void startOfLibrary(boost::filesystem::path const& directory) final override;

    void startOfShower(unsigned int const showerId) final override;

    void endOfShower(unsigned int const showerId) final override;

    void endOfLibrary() final override;

    YAML::Node getSummary() const final override;

    YAML::Node getConfig() const final override;

  private:
    std::function<InteractionHistogram()> getCombined_;
    std::function<void()> resetCountersForShower_;
    InteractionHistogramWriterParquet parquet_;
    YAML::Node summary_;
  };

} // namespace corsika

#include <corsika/detail/modules/writers/InteractionHistogramWriter.inl>
