/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <string>

namespace corsika {

  inline InteractionHistogramWriter::InteractionHistogramWriter(
      std::function<InteractionHistogram()> getCombined,
      std::function<void()> resetCountersForShower)
      : getCombined_(std::move(getCombined))
      , resetCountersForShower_(std::move(resetCountersForShower)) {}

  inline void InteractionHistogramWriter::startOfLibrary(
      boost::filesystem::path const& directory) {
    parquet_.startOfLibrary(directory);
    // Always write per-output config so Python/corsika8.io (and standalone tests)
    // see ``name`` and ``type``. Directory basename matches OutputManager's
    // registered name when using output.add("interaction_hist", ...).
    YAML::Node cfg = getConfig();
    cfg["name"] = directory.filename().string();
    YAMLStreamer yaml;
    yaml.writeYAML(cfg, directory / "config.yaml");
    setInit(true);
  }

  inline void InteractionHistogramWriter::startOfShower(unsigned int const /*showerId*/) {
    if (resetCountersForShower_) { resetCountersForShower_(); }
  }

  inline void InteractionHistogramWriter::endOfShower(unsigned int const showerId) {
    parquet_.endOfShower(showerId, getCombined_());
    auto const key = "shower_" + std::to_string(showerId);
    summary_[key]["rows"] = parquet_.rowsLastShower();
    YAML::Node lab_pdgs;
    for (int pdg : parquet_.lastShowerLabProjectilePdgs()) { lab_pdgs.push_back(pdg); }
    summary_[key]["projectile_pdgs_lab_energy"] = lab_pdgs;
    YAML::Node cms_pdgs;
    for (int pdg : parquet_.lastShowerCmsProjectilePdgs()) { cms_pdgs.push_back(pdg); }
    summary_[key]["projectile_pdgs_cms_energy"] = cms_pdgs;

    YAML::Node integrated;
    for (auto const& entry : parquet_.lastShowerIntegratedPerPdg()) {
      integrated[std::to_string(entry.first)] = entry.second;
    }
    summary_[key]["integrated_interactions_per_pdg"] = integrated;
  }

  inline void InteractionHistogramWriter::endOfLibrary() {
    parquet_.endOfLibrary();
    summary_["total_rows"] = parquet_.totalRows();
  }

  inline YAML::Node InteractionHistogramWriter::getSummary() const { return summary_; }

  inline YAML::Node InteractionHistogramWriter::getConfig() const {
    YAML::Node cfg;
    cfg["type"] = "InteractionHistogram";
    cfg["format"] = "parquet";
    cfg["layout"] = "wide_per_pdg_one_file_pair_per_shower";
    YAML::Node per_shower;
    per_shower["lab"] = "lab_energy_<shower_id>.parquet";
    per_shower["cms"] = "cms_energy_<shower_id>.parquet";
    cfg["per_shower_files"] = per_shower;

    YAML::Node units;
    units["energy_low_eV_energy_high_eV"] =
        "eV; laboratory-frame projectile kinetic energy (bin edges)";
    units["sqrt_s_low_eV_sqrt_s_high_eV"] = "eV; center-of-mass sqrt(s) (bin edges)";
    units["pdg_<PDG>"] = "interaction counts per energy bin (dimensionless)";
    cfg["column_units"] = units;

    cfg["per_shower_histogram"] = static_cast<bool>(resetCountersForShower_);
    return cfg;
  }

} // namespace corsika
