/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>

namespace corsika {
  namespace media {

    template <typename T>
    template <typename... Args>
    inline media::GladstoneDaleRefractiveIndex<T>::GladstoneDaleRefractiveIndex(
        double const referenceRefractiveIndex, Point const point, Args&&... args)
        : T(std::forward<Args>(args)...)
        , referenceRefractivity_(referenceRefractiveIndex - 1)
        , referenceInvDensity_(1 / this->getMassDensity(point)) {}

    template <typename T>
    inline double media::GladstoneDaleRefractiveIndex<T>::getRefractiveIndex(
        Point const& point) const {
      return referenceRefractivity_ *
                 (this->getMassDensity(point) * referenceInvDensity_) +
             1.;
    }

  } // namespace media
} // namespace corsika
