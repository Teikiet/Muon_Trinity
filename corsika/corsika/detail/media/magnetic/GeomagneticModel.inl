/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/Logging.hpp>
#include <boost/math/tr1.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <stdexcept>
#include <string>
#include <cmath>

namespace corsika {
  namespace media {

    inline GeomagneticModel::GeomagneticModel(const CoordinateSystem& centerCS,
                                              Point const& center,
                                              boost::filesystem::path const path)
        : center_(center)
        , centerCS_(centerCS) {

      // Read in coefficients
      boost::filesystem::ifstream file(path, std::ios::in);

      // Exit if file opening failed
      if (!file.is_open()) {
        CORSIKA_LOG_ERROR("Failed opening data file {}", path);
        throw std::runtime_error("Cannot load media::GeomagneticModel data.");
      }

      // media::GeomagneticModel supports two types of input data: WMM.COF and IGRF.COF
      // They have only slightly different format and content and can be easily
      // differentiated here.

      std::string line;
      while (getline(file >> std::ws, line)) {

        double epoch;
        std::string model_name;
        std::string release_date; // just for WMM
        int nMax =
            12; // the spherical max n (l) shell (for IGRF), for WMM this is 12
                // Note that n=l=0 is the monopole and is not included in the model.
        int dummyInt;
        double dummyDouble;

        std::istringstream firstLine(line);

        // check comments and ignore:
        if (firstLine.peek() == '#' ||                     // normal comment
            line.size() == 0 ||                            // empty line
            line.find("9999999999999999999999999") == 0) { // crazy WMM comment
          continue;
        }

        // check IGRF format:
        if (firstLine >> model_name >> epoch >> nMax >> dummyInt >> dummyInt >>
            dummyDouble >> dummyDouble >> dummyDouble >> dummyDouble >> model_name >>
            dummyInt) {
          static bool info = false;
          if (!info) {
            CORSIKA_LOG_INFO("Reading IGRF input data format.");
            info = true;
          }
        } else {
          // check WMM format:
          firstLine.clear();
          firstLine.seekg(0, std::ios::beg);
          if (firstLine >> epoch >> model_name >> release_date) {
            CORSIKA_LOG_INFO("Reading WMM input data format.");
          } else {
            CORSIKA_LOG_ERROR("line: {}", line);
            throw std::runtime_error("Incompatible input data for GeomagneticModel");
          }
        }

        int nPar = 0;
        for (int i = 0; i < nMax; ++i) { nPar += i + 2; }
        int iEpoch = int(epoch);

        if (parameters_.count(iEpoch) != 0) {
          throw std::runtime_error(
              "media::GeomagneticModel input file has duplicate Epoch. Fix.");
        }
        parameters_[iEpoch] = std::vector<ParameterLine>(nPar);

        for (int i = 0; i < nPar; i++) {
          file >> parameters_[iEpoch][i].n >> parameters_[iEpoch][i].m >>
              parameters_[iEpoch][i].g >> parameters_[iEpoch][i].h >>
              parameters_[iEpoch][i].dg >> parameters_[iEpoch][i].dh;
          file.ignore(9999999, '\n');
        }
      }
      file.close();

      if (parameters_.size() == 0) {
        CORSIKA_LOG_ERROR("No input data read!");
        throw std::runtime_error("No input data read");
      }
    }

    inline MagneticFieldVector GeomagneticModel::getField(double const year,
                                                          LengthType const altitude,
                                                          double const latitude,
                                                          double const longitude) {

      int iYear = int(year);
      auto iEpoch = parameters_.rbegin();
      for (; iEpoch != parameters_.rend(); ++iEpoch) {
        if (iEpoch->first <= iYear) { break; }
      }
      CORSIKA_LOG_DEBUG("Found Epoch {} for year {}", iEpoch->first, year);
      if (iEpoch == parameters_.rend()) {
        CORSIKA_LOG_WARN("Year {} is before first EPOCH. Results unclear.", year);
        iEpoch--; // move one epoch back
      }
      if (altitude < -1_km || altitude > 600_km) {
        CORSIKA_LOG_WARN("Altitude should be between -1_km and 600_km.");
      }
      if (latitude <= -90 || latitude >= 90) {
        CORSIKA_LOG_ERROR("Latitude has to be between -90 and 90 degree.");
        throw std::runtime_error("Latitude has to be between -90 and 90 degree.");
      } else if (latitude < -89.992 || latitude > 89.992) {
        CORSIKA_LOG_WARN("Latitude is close to the poles.");
      }
      if (longitude < -180 || longitude > 180) {
        CORSIKA_LOG_WARN("Longitude should be between -180 and 180 degree.");
      }

      double const epoch = double(iEpoch->first);
      auto iNextEpoch = iEpoch; // next epoch for interpolation
      --iNextEpoch;
      bool const lastEpoch = (iEpoch == parameters_.rbegin());
      auto const delta_t = year - epoch;
      CORSIKA_LOG_DEBUG(
          "identified: t_epoch={}, delta_t={}, lastEpoch={} (false->interpolate)", epoch,
          delta_t, lastEpoch);

      double const lat_geo = latitude * constants::pi / 180;
      double const lon = longitude * constants::pi / 180;

      // Transform into spherical coordinates
      double constexpr f = 1 / 298.257223563;
      double constexpr e_squared = f * (2 - f);
      LengthType R_c =
          constants::EarthRadius::Equatorial / sqrt(1 - e_squared * pow(sin(lat_geo), 2));
      LengthType p = (R_c + altitude) * cos(lat_geo);
      LengthType z = sin(lat_geo) * (altitude + R_c * (1 - e_squared));
      LengthType r = sqrt(p * p + z * z);
      double lat_sph = asin(z / r);

      double legendre, next_legendre, derivate_legendre;
      double magneticfield[3] = {0, 0, 0};

      for (size_t j = 0; j < iEpoch->second.size(); j++) {

        ParameterLine p = iEpoch->second[j];

        // Time interpolation
        if (iEpoch == parameters_.rbegin()) {
          // this is the latest epoch in time, or time-dependence (dg/dh) was specified
          // we use the extrapolation factors dg/dh:
          p.g = p.g + delta_t * p.dg;
          p.h = p.h + delta_t * p.dh;
        } else {
          // we linearly interpolate between two epochs
          ParameterLine const next_p = iNextEpoch->second[j];
          double const length = iNextEpoch->first - epoch;
          double p_g = p.g + (next_p.g - p.g) * delta_t / length;
          double p_h = p.h + (next_p.h - p.h) * delta_t / length;
          CORSIKA_LOG_TRACE(
              "interpolation: delta-g={}, delta-h={}, delta-t={}, length={} g1={} g2={} "
              "g={} h={} ",
              next_p.g - p.g, next_p.h - p.h, year - epoch, length, next_p.g, p.g, p_g,
              p_h);
          p.g = p_g;
          p.h = p_h;
        }

        legendre = boost::math::tr1::assoc_legendre(p.n, p.m, sin(lat_sph));
        next_legendre = boost::math::tr1::assoc_legendre(p.n + 1, p.m, sin(lat_sph));

        // Schmidt semi-normalization
        if (p.m > 0) {
          // Note: n! = tgamma(n+1)
          legendre *= sqrt(2 * std::tgamma(p.n - p.m + 1) / std::tgamma(p.n + p.m + 1));
          next_legendre *=
              sqrt(2 * std::tgamma(p.n + 1 - p.m + 1) / std::tgamma(p.n + 1 + p.m + 1));
        }
        derivate_legendre =
            (p.n + 1) * tan(lat_sph) * legendre -
            sqrt(pow(p.n + 1, 2) - pow(p.m, 2)) / cos(lat_sph) * next_legendre;

        magneticfield[0] +=
            pow(constants::EarthRadius::Geomagnetic_reference / r, p.n + 2) *
            (p.g * cos(p.m * lon) + p.h * sin(p.m * lon)) * derivate_legendre;
        magneticfield[1] +=
            pow(constants::EarthRadius::Geomagnetic_reference / r, p.n + 2) * p.m *
            (p.g * sin(p.m * lon) - p.h * cos(p.m * lon)) * legendre;
        magneticfield[2] +=
            (p.n + 1) * pow(constants::EarthRadius::Geomagnetic_reference / r, p.n + 2) *
            (p.g * cos(p.m * lon) + p.h * sin(p.m * lon)) * legendre;
      }
      magneticfield[0] *= -1;
      magneticfield[1] /= cos(lat_sph);
      magneticfield[2] *= -1;

      // Transform back into geodetic coordinates
      double magneticfield_geo[3];
      magneticfield_geo[0] = magneticfield[0] * cos(lat_sph - lat_geo) -
                             magneticfield[2] * sin(lat_sph - lat_geo);
      magneticfield_geo[1] = magneticfield[1];
      magneticfield_geo[2] = magneticfield[0] * sin(lat_sph - lat_geo) +
                             magneticfield[2] * cos(lat_sph - lat_geo);

      return MagneticFieldVector{centerCS_, magneticfield_geo[0] * 1_nT,
                                 magneticfield_geo[1] * -1_nT,
                                 magneticfield_geo[2] * -1_nT};
    }

  } // namespace media
} // namespace corsika
