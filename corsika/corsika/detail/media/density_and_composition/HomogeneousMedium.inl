/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

namespace corsika {
  namespace media {

    template <typename T>
    inline HomogeneousMedium<T>::HomogeneousMedium(
        MassDensityType density, media::NuclearComposition const& nuclComp)
        : density_(density)
        , nuclComp_(nuclComp) {}

    template <typename T>
    inline MassDensityType HomogeneousMedium<T>::getMassDensity(Point const&) const {
      return density_;
    }

    template <typename T>
    inline media::NuclearComposition const& HomogeneousMedium<T>::getNuclearComposition()
        const {
      return nuclComp_;
    }

    template <typename T>
    inline GrammageType HomogeneousMedium<T>::getIntegratedGrammage(
        BaseTrajectory const& track) const {
      return track.getLength() * density_;
    }

    template <typename T>
    inline LengthType HomogeneousMedium<T>::getArclengthFromGrammage(
        BaseTrajectory const&, GrammageType grammage) const {
      return grammage / density_;
    }

  } // namespace media
} // namespace corsika
