/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

#include <limits>
#include <vector>

namespace corsika {
  namespace media {

    /**
     * This class provides the grammage/length conversion functionality for
     * (locally) flat tabulated atmospheres.
     */
    template <typename TDerived>
    class BaseTabular {

    public:
      BaseTabular(Point const& point, LengthType const referenceHeight,
                  std::function<MassDensityType(LengthType)> const& rho,
                  unsigned int const nBins, LengthType const deltaHeight);

      Point const& getAnchorPoint() const { return point_; }

    protected:
      auto const& getImplementation() const;

      MassDensityType getMassDensity(LengthType const height) const;

      GrammageType getIntegratedGrammage(BaseTrajectory const& line) const;

      LengthType getArclengthFromGrammage(BaseTrajectory const& line,
                                          GrammageType const grammage) const;

    private:
      unsigned int const nBins_;
      LengthType deltaHeight_;
      std::vector<MassDensityType> density_;
      Point const point_;
      LengthType const referenceHeight_;

    }; // class BaseTabular

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density/BaseTabular.inl>
