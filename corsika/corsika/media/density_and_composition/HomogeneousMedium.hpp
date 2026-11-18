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
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {
  namespace media {

    /**
     * a homogeneous medium
     */

    template <typename T>
    class HomogeneousMedium : public T {

    public:
      HomogeneousMedium(MassDensityType density,
                        media::NuclearComposition const& nuclComp);

      MassDensityType getMassDensity(Point const&) const override;

      media::NuclearComposition const& getNuclearComposition() const override;

      GrammageType getIntegratedGrammage(BaseTrajectory const&) const override;

      LengthType getArclengthFromGrammage(BaseTrajectory const&,
                                          GrammageType grammage) const override;

    private:
      MassDensityType const density_;
      media::NuclearComposition const nuclComp_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density_and_composition/HomogeneousMedium.inl>
