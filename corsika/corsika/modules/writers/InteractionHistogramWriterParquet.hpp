/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/InteractionHistogram.hpp>

#include <boost/filesystem.hpp>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace corsika {

  /**
   * Parquet output for InteractionHistogram: **one pair of files per shower** —
   * `lab_energy_<shower_id>.parquet` and `cms_energy_<shower_id>.parquet` — each a wide
   * table (one row per energy bin, one FLOAT column per projectile PDG on axis 0).
   * Schema is local to that shower, so new projectile types do not require column
   * backfill across showers.
   */
  class InteractionHistogramWriterParquet {

  public:
    InteractionHistogramWriterParquet() = default;

    void startOfLibrary(boost::filesystem::path const& directory);

    void endOfShower(unsigned int showerId, InteractionHistogram const& h);

    void endOfLibrary();

    std::uint64_t rowsLastShower() const { return rowsLastShower_; }
    std::uint64_t totalRows() const { return totalRows_; }

    std::vector<int> const& lastShowerLabProjectilePdgs() const { return last_lab_pdgs_; }

    std::vector<int> const& lastShowerCmsProjectilePdgs() const { return last_cms_pdgs_; }

    /**
     * Sum of bin counts over all energy bins per projectile PDG (from the lab histogram;
     * same totals as summing the CMS histogram per PDG).
     */
    std::map<int, double> const& lastShowerIntegratedPerPdg() const {
      return last_integrated_per_pdg_;
    }

  private:
    boost::filesystem::path directory_;
    std::uint64_t totalRows_{0};
    std::uint64_t rowsLastShower_{0};
    std::vector<int> last_lab_pdgs_;
    std::vector<int> last_cms_pdgs_;
    std::map<int, double> last_integrated_per_pdg_;
  };

} // namespace corsika

#include <corsika/detail/modules/writers/InteractionHistogramWriterParquet.inl>
