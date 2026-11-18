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

    template <typename T, typename TDensityFunction>
    template <typename... TArgs>
    inline InhomogeneousMedium<T, TDensityFunction>::InhomogeneousMedium(
        media::NuclearComposition const& nuclComp, TArgs&&... rhoTArgs)
        : nuclComp_(nuclComp)
        , densityFunction_(rhoTArgs...) {}

    template <typename T, typename TDensityFunction>
    inline MassDensityType InhomogeneousMedium<T, TDensityFunction>::getMassDensity(
        Point const& point) const {
      return densityFunction_.evaluateAt(point);
    }

    template <typename T, typename TDensityFunction>
    inline media::NuclearComposition const&
    InhomogeneousMedium<T, TDensityFunction>::getNuclearComposition() const {
      return nuclComp_;
    }

    template <typename T, typename TDensityFunction>
    inline GrammageType InhomogeneousMedium<T, TDensityFunction>::getIntegratedGrammage(
        BaseTrajectory const& line) const {
      return densityFunction_.getIntegrateGrammage(line);
    }

    template <typename T, typename TDensityFunction>
    inline LengthType InhomogeneousMedium<T, TDensityFunction>::getArclengthFromGrammage(
        BaseTrajectory const& line, GrammageType grammage) const {
      return densityFunction_.getArclengthFromGrammage(line, grammage);
    }

  } // namespace media
} // namespace corsika
