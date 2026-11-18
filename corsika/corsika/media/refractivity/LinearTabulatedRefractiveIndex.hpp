/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>
#include <corsika/media/interfaces/IGeometryTransformer.hpp>
#include <corsika/media/geometry_transformer/planar.hpp>
#include <filesystem>

namespace corsika {
  namespace media {
    /**
     * An exponential refractive index.
     *
     * This class returns the value of an exponential refractive index
     * for all evaluated locations.
     *
     */
    template <typename T, typename TGeometry = PlanarTransformer>
    class LinearTabulatedRefractiveIndex : public T {

      std::filesystem::path tabulated_amosphere_;

      struct table_row_data {
        double altitude;
        double density;
        double grammage;
        double refractive_index;
      };

      double max_height_;

      std::vector<double> refractive_index_profile_;

      const TGeometry transformer_;

    protected:
      std::vector<typename LinearTabulatedRefractiveIndex<T, TGeometry>::table_row_data>
      parse_file(std::string_view const& tabulated_amosphere_path);
      void interpolate(std::vector<typename LinearTabulatedRefractiveIndex<
                           T, TGeometry>::table_row_data> const& refractive_index_data);

      double log_interpolate_between(double heigh1, double value1, double height2,
                                     double value2, double query_height) const;

    public:
      /**
       * Construct an LinearTabulatedRefractiveIndex.
       *
       * This is initialized with the path to a tabulated atmosphere in csv format with
       * two columns (height, refractive index) and returns the value of the refractive
       * index at a given point location.
       *
       * @param field    The refractive index to return to a given point.
       */
      template <typename... Args>
      LinearTabulatedRefractiveIndex(
          std::filesystem::path const& tabulated_amosphere_path, TGeometry const&,
          Args&&... args);

      /**
       * Evaluate the refractive index at a given location using its z-coordinate.
       *
       * @param  point    The location to evaluate at.
       * @returns    The refractive index at this point.
       */
      double getRefractiveIndex(Point const& point) const final override;

    }; // END: class LinearTabulatedRefractiveIndex

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/refractivity/LinearTabulatedRefractiveIndex.inl>