/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/utility/CorsikaData.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/detail/media/LayeredSphericalAtmosphereBuilder.hpp>

#include <boost/filesystem.hpp>
#include <string>

namespace corsika {
  namespace media {
    /**
     * A tabulated atmosphere model constructed from an external atmospheric profile table
     *
     * The expected format and units of column values are
     * Alt[km]  rho[g/cm^3]  thick[g/cm^2]  n-1
     *
     * Altitude and density values are required, thickness and refractive index
     * values are optional
     *
     * For examples, see tables in modules/data/CHERENKOV/atmosphere
     */

    class AtmosphereTable {

    public:
      /**
       * Construct a new AtmosphereTable object.
       *
       * @param center Atmosphere/planet center
       * @param path Data table to read
       */
      AtmosphereTable(Point const& center, boost::filesystem::path const& path);

      /**
       * Returns the center of the atmosphere model
       * @returns Center of the atmosphere model
       */
      inline Point const& getCenter() const { return center_; }

      /**
       * Returns a vector of altitudes of layer edges (measured from the surface)
       * @returns Vector of altitudes of layer edges (measured from the surface)
       */
      std::vector<LengthType> const& getAltitudeValues() const;

      /**
       * Returns a vector of densities at layer edges
       * @returns Vector of densities at layer edges
       */
      std::vector<MassDensityType> const& getDensityValues() const;

      /**
       * Returns a vector of atmosphere thicknesses (cumulated densities) at layer edges
       * @returns Vector of atmosphere thicknesses at layer edges
       */
      std::vector<GrammageType> const& getThicknessValues() const;

      /**
       * Returns a vector of refractive indices at layer edges
       * @returns Vector of refractive indices at layer edges
       */
      std::vector<double> const& getRefractiveIndexValues() const;

      /**
       * Checks the parsed table for consistency, sortedness and proper values
       */
      void check_table();

    private:
      // center of atmosphere/planet
      Point center_;

      // path to table
      boost::filesystem::path tablePath_;

      // vectors to hold the atmosphere parameters
      std::vector<LengthType> alt_;
      std::vector<MassDensityType> dens_;
      std::vector<GrammageType> thick_;
      std::vector<double> refr_;
    };

    /**
     * Function to create an atmosphere using tabulated values
     * @tparam TEnvironmentInterface
     * @tparam TExtraEnv
     * @tparam TEnvironment
     * @tparam TArgs
     * @param env Environment to create the atmosphere in
     * @param table Table describing the atmosphere
     * @param args
     */
    template <typename TEnvironmentInterface,
              template <typename> typename TExtraEnv = detail::NoExtraModel,
              typename TEnvironment, typename... TArgs>
    void create_interpolated_atmosphere(TEnvironment& env, AtmosphereTable const& table,
                                        TArgs... args);
  } // namespace media
} // namespace corsika

#include <corsika/detail/media/AtmosphereTable.inl>
