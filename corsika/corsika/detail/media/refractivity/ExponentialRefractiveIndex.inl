/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>

namespace corsika {
  namespace media {

    template <typename T, typename TGeometry>
    template <typename... TArgs>
    ExponentialRefractiveIndex<T, TGeometry>::ExponentialRefractiveIndex(
        double const n0, InverseLengthType const lambda, Point const center,
        LengthType const radius, TArgs&&... args)
        : T(std::forward<TArgs>(args)...)
        , n0_(n0)
        , lambda_(lambda)
        , center_(center)
        , radius_(radius) {}

    template <typename T, typename TGeometry>
    inline double ExponentialRefractiveIndex<T, TGeometry>::getRefractiveIndex(
        Point const& point) const {
      return n0_ * exp((-lambda_) * (distance(point, center_) - radius_));
    }

  } // namespace media
} // namespace corsika
