/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>
#include <corsika/media/geometry_transformer/planar.hpp>

namespace corsika {
  namespace media {

    /**
     * An exponential refractive index.
     *
     * This class returns the value of an exponential refractive index
     * for all evaluated locations.
     *
     */
    template <typename T, class TGeometry = PlanarTransformer>
    class ExponentialRefractiveIndex : public T {

      double n0_;                ///< n0 constant.
      InverseLengthType lambda_; ///< lambda parameter.
      Point center_;             ///< center of the planet.
      LengthType radius_;        ///< the planet radius.

    public:
      /**
       * Construct an ExponentialRefractiveIndex.
       *
       * This is initialized with two parameters n_0 and lambda
       * and returns the value of the exponential refractive index
       * at a given point location.
       *
       * @param field    The refractive index to return to a given point.
       */
      template <typename... Args>
      ExponentialRefractiveIndex(double const n0, InverseLengthType const lambda,
                                 Point const center, LengthType const radius,
                                 Args&&... args);

      /**
       * Evaluate the refractive index at a given location using its z-coordinate.
       *
       * @param  point    The location to evaluate at.
       * @returns    The refractive index at this point.
       */
      double getRefractiveIndex(Point const& point) const final override;

    }; // END: class ExponentialRefractiveIndex

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/refractivity/ExponentialRefractiveIndex.inl>
