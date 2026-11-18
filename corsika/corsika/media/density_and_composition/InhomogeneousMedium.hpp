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
     * A general inhomogeneous medium. The mass density distribution TDensityFunction must
     * be a \f$C^2\f$-function.
     */

    template <typename T, typename TDensityFunction>
    class InhomogeneousMedium : public T {

    public:
      template <typename... TArgs>
      InhomogeneousMedium(media::NuclearComposition const& nuclComp, TArgs&&... rhoTArgs);

      MassDensityType getMassDensity(Point const& point) const override;

      media::NuclearComposition const& getNuclearComposition() const override;

      GrammageType getIntegratedGrammage(BaseTrajectory const& line) const override;

      LengthType getArclengthFromGrammage(BaseTrajectory const& pLine,
                                          GrammageType grammage) const override;

    private:
      media::NuclearComposition const nuclComp_;
      TDensityFunction const densityFunction_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density_and_composition/InhomogeneousMedium.inl>
