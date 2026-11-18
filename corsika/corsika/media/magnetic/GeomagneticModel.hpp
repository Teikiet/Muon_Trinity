/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/utility/CorsikaData.hpp>

#include <boost/filesystem.hpp>
#include <fstream>
#include <string>

namespace corsika {
  namespace media {

    /**
     * A magnetic field calculated with the WMM or IGRF model.
     * WMM: https://geomag.bgs.ac.uk/documents/WMM2020_Report.pdf
     * IGRF: https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
     */

    class GeomagneticModel {

      /**
       * Internal data structure for a single shell of the spherical harmonic
       * parameterization.
       */
      struct ParameterLine {
        int n;
        int m;
        double g;
        double h;
        double dg; //< time dependence of g in "per year"
        double dh; //< time dependence of h in "per year"
      };

    public:
      /**
       * Construct a new World Magnetic Model object.
       *
       * @param centerCS Coordinate system of the center.
       * @param center Center of Earth.
       * @param data Data table to read.
       */
      GeomagneticModel(
          CoordinateSystem const& centerCS, Point const& center,
          boost::filesystem::path const path = corsika::corsika_data("GeoMag/WMM.COF"));

      /**
       * Calculates the value of the magnetic field.
       *
       * @param  year        Year of the evaluation, between 2020 and 2025.
       * @param  altitude    Height of the location to evaluate the field at,
       *                    in km between -1 and 850.
       * @param  latitude    Latitude of the location to evaluate the field at,
       *                   in degrees between -90 and 90 (negative for southern
       * hemisphere).
       * @param  longitute   Longitude of the location to evaluate the field at,
       *                      in degrees between -180 and 180 (negative for western
       *   hemisphere).
       *
       * @returns    The magnetic field vector in nT.
       */
      MagneticFieldVector getField(double const year, LengthType const altitude,
                                   double const latitude, double const longitude);

    private:
      Point center_; //< Center of Earth.
      const CoordinateSystem& centerCS_;
      std::map<int, std::vector<ParameterLine>> parameters_; //< epoch to Parameters map
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/magnetic/GeomagneticModel.inl>