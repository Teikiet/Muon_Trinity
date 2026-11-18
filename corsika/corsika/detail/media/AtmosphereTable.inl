/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <algorithm>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <boost/math/tr1.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <cmath>

namespace corsika {
  namespace media {
    inline AtmosphereTable::AtmosphereTable(Point const& center,
                                            boost::filesystem::path const& path)
        : center_(center)
        , tablePath_(path) {

      // open input table file
      boost::filesystem::ifstream file(tablePath_, std::ios::in);

      // exit if file opening failed
      if (!file.is_open()) {
        CORSIKA_LOG_ERROR("Failed opening atmosphere table '{}'", tablePath_);
        throw std::runtime_error("Cannot load AtmosphereTable data.");
      }

      std::string line;
      std::vector<std::string> lineSplit;
      std::vector<double> params;
      int lineId = 1;

      // read table
      while (getline(file >> std::ws, line)) {
        // skip empty lines and comments
        if (line.size() == 0 || line.at(0) == '#') {
          lineId++;
          continue;
        }

        // split line string by whitespaces
        boost::split(lineSplit, line, boost::is_any_of(" \t\n\r"));

        // extract parameters from the line
        for (const auto& subStr : lineSplit) {
          // skip empty strings
          if (subStr.empty()) continue;

          std::istringstream iss(subStr);
          double val;

          // check if line substring contains a proper value
          if (iss >> val && iss.eof() && std::isfinite(val)) {
            params.push_back(val);
          } else {
            CORSIKA_LOG_ERROR(
                "Atmosphere table '{}' contains NaN ('{}') on line {}, column {}",
                tablePath_, subStr, lineId, params.size() + 1);
            throw std::runtime_error("Atmosphere table contains NaN");
          }
        }

        // check number of parameters - at least two (altitude and density),
        // optionally two more (thickness, refractive index)
        if (params.size() < 2) {
          CORSIKA_LOG_ERROR(
              "Atmosphere table '{}' has too few parameters on line {} "
              "(expected at least 2, got {})",
              tablePath_, lineId, lineSplit.size());
          CORSIKA_LOG_ERROR(
              "Expected table format is: altitude [km], density [g/cm^3], "
              "OPTIONALLY thickness [g/cm^2], refractive index [as n-1]");
          throw std::runtime_error(
              "Atmosphere table has an unexpected number of parameters");
        }

        // store altitude and density values
        alt_.push_back(params[0] * 1_km);
        dens_.push_back(params[1] * 1_g / cube(1_cm));

        // if provided, also store thickness and refractive index values
        if (params.size() > 2) { thick_.push_back(params[2] * 1_g / square(1_cm)); }

        if (params.size() > 3) {
          // values in the table expected as n-1, so add 1
          refr_.push_back(1 + params[3]);
        }

        params.clear();
        lineId++;
      }

      file.close();

      // test whether the parsed table is valid
      check_table();
    }

    void AtmosphereTable::check_table() {
      // check number of parameter consistency across rows
      if (alt_.size() != dens_.size() ||
          (thick_.size() != 0 && alt_.size() != thick_.size()) ||
          (refr_.size() != 0 && alt_.size() != refr_.size())) {
        CORSIKA_LOG_ERROR(
            "Atmosphere table '{}' has inconsistent numbers of parameters in rows",
            tablePath_);
        throw std::runtime_error(
            "Atmosphere table has inconsistent numbers of parameters in rows");
      }

      // number of rows
      auto nRows = alt_.size();

      // check size of table (at least 2 rows)
      if (nRows < 2) {
        CORSIKA_LOG_ERROR("Atmosphere table '{}' has fewer than two rows", tablePath_);
        throw std::runtime_error("Atmosphere table has fewer than two rows");
      }

      // check atmosphere is defined at zero altitude
      if (alt_[0] > 0_m) {
        CORSIKA_LOG_WARN(
            "Atmosphere undefined at zero altitude in table '{}', using parameters of "
            "the first point above ground (altitude {})",
            tablePath_, alt_[0]);

        // add entry at zero altitude and use parameters of the first point above
        alt_.insert(alt_.begin(), 0_m);
        dens_.insert(dens_.begin(), dens_[0]);
        thick_.insert(thick_.begin(), thick_[0]);
        refr_.insert(refr_.begin(), refr_[0]);
      }

      // check table for negative density/thickness
      for (std::size_t n = 0; n < nRows; n++) {
        if (dens_[n] < MassDensityType::zero()) {
          CORSIKA_LOG_ERROR(
              "Negative atmosphere density defined in table '{}' at altitude {}",
              tablePath_, alt_[n]);
          throw std::runtime_error("Negative atmosphere density defined at altitude " +
                                   std::to_string(alt_[n] / 1_km) + " km");
        }

        if (thick_.size() != 0 && thick_[n] < GrammageType::zero()) {
          CORSIKA_LOG_ERROR(
              "Negative atmosphere thickness defined in table '{}' at altitude {}",
              tablePath_, alt_[n]);
          throw std::runtime_error("Negative atmosphere thickness defined at altitude " +
                                   std::to_string(alt_[n] / 1_km) + " km");
        }
      }

      // check whether table rows are sorted by ascending altitude
      if (!std::is_sorted(alt_.begin(), alt_.end())) {
        CORSIKA_LOG_ERROR("Table rows in table '{}' are not sorted by ascending altitude",
                          tablePath_);
        throw std::runtime_error("Table rows are not sorted by ascending altitude");
      }

      // duplicate entries
      for (unsigned int n = 1; n < nRows; n++) {
        if (alt_[n - 1] == alt_[n]) {
          CORSIKA_LOG_ERROR("Atmosphere table '{}' has duplicate entries for altitude {}",
                            tablePath_, alt_[n]);
          throw std::runtime_error("Atmosphere table has duplicate entries for " +
                                   std::to_string(alt_[n] / 1_km) + " km altitude");
        }
      }
    }

    const std::vector<LengthType>& AtmosphereTable::getAltitudeValues() const {
      return alt_;
    }

    const std::vector<MassDensityType>& AtmosphereTable::getDensityValues() const {
      return dens_;
    }

    const std::vector<GrammageType>& AtmosphereTable::getThicknessValues() const {
      if (thick_.size() != 0) {
        return thick_;
      } else {
        CORSIKA_LOG_ERROR("Atmosphere thickness queried, but not defined in table '{}'",
                          tablePath_);
        throw(std::runtime_error(
            "Atmosphere thickness queried, but not defined in table '" +
            tablePath_.string() + "'"));
      }
    }

    const std::vector<double>& AtmosphereTable::getRefractiveIndexValues() const {
      if (refr_.size() != 0) {
        return refr_;
      } else {
        CORSIKA_LOG_ERROR("Refractive index queried, but not defined in table '{}'",
                          tablePath_);
        throw(std::runtime_error("Refractive index queried, but not defined in table '" +
                                 tablePath_.string() + "'"));
      }
    }

    template <typename TEnvironmentInterface, template <typename> typename TExtraEnv,
              typename TEnvironment, typename... TArgs>
    void create_interpolated_atmosphere(TEnvironment& env, AtmosphereTable const& table,
                                        TArgs... args) {

      // construct the atmosphere builder
      auto builder = media::make_layered_spherical_atmosphere_builder<
          TEnvironmentInterface, TExtraEnv>::create(table.getCenter(),
                                                    constants::EarthRadius::Mean,
                                                    std::forward<TArgs>(args)...);

      // air nuclear composition
      builder.setNuclearComposition(corsika::media::standardAirComposition);

      // get vectors of altitudes and densities
      auto const altitudeVec = table.getAltitudeValues();
      auto const densityVec = table.getDensityValues();

      builder.addInterpolatedLayer(altitudeVec, densityVec);

      builder.assemble(env);
    }
  } // namespace media
} // namespace corsika
