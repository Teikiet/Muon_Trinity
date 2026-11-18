/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/FindXmax.hpp>

#include <corsika/media/ShowerAxis.hpp>

#include <exception>
#include <algorithm>
#include <iostream>

namespace corsika {

  template <typename TOutput>
  inline ProductionWriter<TOutput>::ProductionWriter(media::ShowerAxis const& axis,
                                                     GrammageType dX)
      : ProductionWriter<TOutput>{
            axis, static_cast<unsigned int>(std::ceil(axis.getMaximumX() / dX)) + 1,
            dX} {}

  template <typename TOutput>
  inline ProductionWriter<TOutput>::ProductionWriter(media::ShowerAxis const& axis,
                                                     size_t nbins, GrammageType dX)
      : TOutput(production_profile::ProjectileIndexNames)
      , showerAxis_(axis)
      , dX_(dX)
      , nBins_(nbins)
      , profile_{nbins} {
    summary_ = YAML::Node();
  }

  template <typename TOutput>
  inline void ProductionWriter<TOutput>::startOfLibrary(
      boost::filesystem::path const& directory) {
    TOutput::startOfLibrary(directory);
  }

  template <typename TOutput>
  inline void ProductionWriter<TOutput>::startOfShower(unsigned int const showerId) {
    profile_.clear();
    for (size_t i = 0; i < nBins_; ++i) { profile_.emplace_back(); }
    TOutput::startOfShower(showerId);
  }

  template <typename TOutput>
  inline void ProductionWriter<TOutput>::endOfLibrary() {
    TOutput::endOfLibrary();
  }

  template <typename TOutput>
  inline void ProductionWriter<TOutput>::endOfShower(unsigned int const showerId) {

    size_t constexpr window_size = 20;

    if (profile_.size() < window_size) {
      CORSIKA_LOGGER_WARN(
          TOutput::getLogger(),
          "MPD profile too small (size={}) for window {}. Skipping Xmax estimation.",
          profile_.size(), window_size);
    } else {
      // determined XMumax and dNdXmax from quadratic interpolation
      double maximum = 0;
      size_t iMaximum = 0;

      int constexpr HadronIdx =
          static_cast<int>(corsika::production_profile::ProjectileIndex::Hadron);

      auto x_at = [&](size_t i) { return (i + 0.5) * (dX_ / 1_g * square(1_cm)); };
      auto y_at = [&](size_t i) { return profile_[i].at(HadronIdx); };

      for (size_t i = 0; i <= profile_.size() - window_size; ++i) {
        double value = 0;
        for (size_t j = 0; j < window_size; ++j) { value += y_at(i + j); }
        if (value > maximum) {
          maximum = value;
          iMaximum = i;
        }
      }

      std::vector<double> xs, ys;
      for (size_t j = 0; j < window_size; ++j) {
        xs.push_back(x_at(iMaximum + j));
        ys.push_back(y_at(iMaximum + j));
      }

      if (auto result = find_xmax::estimate_parabola_parameters(xs, ys)) {
        auto [Xmumax, dNdXmumax] = *result;
        summary_["shower_" + std::to_string(showerId)]["XmuMax"] = Xmumax;
        summary_["shower_" + std::to_string(showerId)]["dNdXmuMax"] = dNdXmumax;
      } else {
        CORSIKA_LOGGER_WARN(TOutput::getLogger(),
                            "XmuMax estimation failed. Skipping XmuMax.");
      }
    }

    // write profile to file
    int iRow{0};
    for (production_profile::ProfileData const& row : profile_) {
      // here: write to underlying writer (e.g. parquet)
      TOutput::write(showerId, iRow * dX_, row);
      iRow++;
    }
    TOutput::endOfShower(showerId);
  }

  template <typename TOutput>
  inline void ProductionWriter<TOutput>::write(Point const& p0, Code const projectile_pid,
                                               double const weight) {
    GrammageType const grammage = showerAxis_.getProjectedX(p0);

    // Note: particle may go also "upward", thus, grammageEnd<grammageStart
    size_t const bin = std::ceil(grammage / dX_);

    // Guard against a bin outside the profile range. The ShowerAxis is anchored
    // on the nominal axis, so a particle injected far off-axis (e.g. after the
    // magnetic-field back-tracking) can project beyond the axis' maximum
    // grammage and yield bin >= profile_.size(). Mirror LongitudinalWriter.
    if (bin >= profile_.size()) {
      CORSIKA_LOGGER_WARN(TOutput::getLogger(),
                          "Grammage bin {} outside of profile {}. skipping.", bin,
                          profile_.size());
      return;
    }

    CORSIKA_LOGGER_TRACE(TOutput::getLogger(), "grammage={} bin={}",
                         grammage / 1_g * square(1_cm), bin);
    if (projectile_pid == Code::Photon) {
      profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::Photon)] +=
          weight;
    } else if (projectile_pid == Code::Positron || projectile_pid == Code::Electron) {
      profile_.at(
          bin)[static_cast<int>(production_profile::ProjectileIndex::ElectronPositron)] +=
          weight;
    } else if (is_muon(projectile_pid)) {
      profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::Muon)] +=
          weight;
    } else if (is_hadron(projectile_pid)) {
      profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::Hadron)] +=
          weight;
      if (projectile_pid == Code::PiPlus || projectile_pid == Code::PiMinus) {
        profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::Pion)] +=
            weight;
      } else if (is_kaon(projectile_pid)) {
        profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::Kaon)] +=
            weight;
      } else {
        profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::Heavy)] +=
            weight;
      }
    }
    profile_.at(bin)[static_cast<int>(production_profile::ProjectileIndex::All)] +=
        weight;
  }

  template <typename TOutput>
  inline YAML::Node ProductionWriter<TOutput>::getConfig() const {

    YAML::Node node;

    node["type"] = "ProductionProfile";
    node["units"]["grammage"] = "g/cm^2";
    node["bin-size"] = dX_ / (1_g / square(1_cm));
    node["nbins"] = nBins_;

    return node;
  }

  template <typename TOutput>
  inline YAML::Node ProductionWriter<TOutput>::getSummary() const {
    return summary_;
  }

} // namespace corsika
