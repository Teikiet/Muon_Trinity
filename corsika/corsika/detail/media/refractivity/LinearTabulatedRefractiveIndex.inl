/*
 * (c) Copyright 2026 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <cctype>
#include <algorithm>
#include <cmath>
#include <limits>

namespace corsika {
  namespace media {

    template <typename T, typename TGeometry>
    template <typename... Args>
    LinearTabulatedRefractiveIndex<T, TGeometry>::LinearTabulatedRefractiveIndex(
        std::filesystem::path const& tabulated_amosphere_path,
        TGeometry const& transformer, Args&&... args)
        : T(std::forward<Args>(args)...)
        , tabulated_amosphere_(tabulated_amosphere_path)
        , transformer_(transformer) {

      // Interpolate between the tabulated points in higher precision
      auto file_data = this->parse_file(tabulated_amosphere_.string());
      this->interpolate(file_data);
    }

    template <typename T, typename TGeometry>
    std::vector<typename LinearTabulatedRefractiveIndex<T, TGeometry>::table_row_data>
    LinearTabulatedRefractiveIndex<T, TGeometry>::parse_file(
        std::string_view const& tabulated_amosphere_path) {
      // CSV read from file
      std::ifstream file(tabulated_amosphere_path.data());

      if (!file) {
        throw std::runtime_error("Failed to open file: " +
                                 std::string(tabulated_amosphere_path));
      }

      std::string line;
      std::vector<table_row_data> refractive_index_data;

      for (int lineNumber = 0; std::getline(file, line); lineNumber++) {

        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

        std::istringstream iss(line);
        table_row_data row;

        if (!(iss >> row.altitude >> row.density >> row.grammage >>
              row.refractive_index)) {
          throw std::runtime_error("Parse error at line " + std::to_string(lineNumber) +
                                   ": expected 4 numeric columns");
        }

        if (iss.fail()) {
          throw std::runtime_error("Could not parse value at line " +
                                   std::to_string(lineNumber));
        }

        if (!std::isfinite(row.altitude) || !std::isfinite(row.density) ||
            !std::isfinite(row.grammage) || !std::isfinite(row.refractive_index)) {
          throw std::runtime_error("Non-finite value at line " +
                                   std::to_string(lineNumber));
        }

        row.altitude = row.altitude * 1000.0; // convert km to m
        row.refractive_index =
            row.refractive_index; // refractive index is already in correct units

        refractive_index_data.push_back(row);
      }

      // Check if we have at least two data points
      if (refractive_index_data.size() < 2) {
        throw std::runtime_error("Less than two rows found in: " +
                                 tabulated_amosphere_.string());
      }

      std::sort(refractive_index_data.begin(), refractive_index_data.end(),
                [](table_row_data const& a, table_row_data const& b) {
                  return a.altitude < b.altitude;
                });

      // Check for duplicate altitudes
      for (size_t i = 1; i < refractive_index_data.size(); ++i) {
        if (refractive_index_data[i].altitude == refractive_index_data[i - 1].altitude) {
          throw std::runtime_error("Duplicate altitude at line " + std::to_string(i) +
                                   " in file: " + tabulated_amosphere_.string());
        }
      }

      // Check if first entry starts at 0 m
      if (refractive_index_data.front().altitude > 0) {
        CORSIKA_LOG_WARN(
            "First altitude in tabulated atmosphere is greater than 0 m. "
            "Filling lower altitudes with first entry's refractive index.");
      }

      // Initialize max_height_
      max_height_ = refractive_index_data.back().altitude;

      return refractive_index_data;
    }

    template <typename T, typename TGeometry>
    void LinearTabulatedRefractiveIndex<T, TGeometry>::interpolate(
        std::vector<
            typename LinearTabulatedRefractiveIndex<T, TGeometry>::table_row_data> const&
            refractive_index_data) {

      if (refractive_index_data[0].altitude > refractive_index_data[1].altitude) {
        logging::error(
            "Error in tabulated atmosphere file: heights are not sorted in ascending "
            "order");

        throw(std::runtime_error(
            "Error in tabulated atmosphere file: heights are not sorted "
            "in ascending order"));
      }

      auto data1 = refractive_index_data[0];
      auto data2 = refractive_index_data[1];
      size_t data_idx = 2;

      refractive_index_profile_.resize(static_cast<size_t>(this->max_height_ + 0.5f));

      for (unsigned int interpol_height = 0;
           interpol_height < static_cast<size_t>(this->max_height_ + 0.5f);
           interpol_height++) {

        if (interpol_height > data2.altitude) {
          if (data_idx <
              refractive_index_data.size()) { // Load new data points if possible
            data1 = data2;
            data2 = refractive_index_data[data_idx];
            data_idx++;
          } else // No new data points available, use last known refractive index
          {
            refractive_index_profile_[interpol_height] = data2.refractive_index;
            continue;
          }
        }

        if (interpol_height <
            data1.altitude) { // Before first data point, fill with first value
          refractive_index_profile_[interpol_height] = data1.refractive_index;
          continue;
        }

        // Default case: interpolate between data1 and data2

        double interpolated_value = this->log_interpolate_between(
            data1.altitude, data1.refractive_index, data2.altitude,
            data2.refractive_index, static_cast<double>(interpol_height));

        refractive_index_profile_[interpol_height] = 1.0 + interpolated_value;
      }
    }

    template <typename T, typename TGeometry>
    inline double LinearTabulatedRefractiveIndex<T, TGeometry>::log_interpolate_between(
        double heigh1, double value1, double height2, double value2,
        double query_height) const {

      double log_value1 = std::log(value1);
      double log_value2 = std::log(value2);

      double height_diff = height2 - heigh1;
      double log_value_diff = log_value2 - log_value1;

      double fraction = (query_height - heigh1) / height_diff;

      return std::exp(log_value1 + fraction * log_value_diff);
    }

    template <typename T, typename TGeometry>
    inline double LinearTabulatedRefractiveIndex<T, TGeometry>::getRefractiveIndex(
        Point const& point) const {
      double altitude = (transformer_.getEffectiveHeight(point) / 1_m);
      if (altitude < 0) {
        throw(std::runtime_error(
            "LinearTabulatedRefractiveIndex: point is inside the earth"));
      }
      if (altitude > this->max_height_) { return this->refractive_index_profile_.back(); }

      int lower = std::floor(altitude);
      int upper = std::ceil(altitude);
      int max_index = refractive_index_profile_.size() - 1;
      lower = std::max(0, std::min(max_index, lower));
      upper = std::max(0, std::min(max_index, upper));
      double lower_index = refractive_index_profile_[lower];
      double upper_index = refractive_index_profile_[upper];
      double fraction = altitude - lower;
      return lower_index + fraction * (upper_index - lower_index);
    }

  } // namespace media
} // namespace corsika
