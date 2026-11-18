/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {
  namespace media {

    class IMediumModel {
    public:
      virtual ~IMediumModel() = default; // LCOV_EXCL_LINE

      virtual MassDensityType getMassDensity(Point const&) const = 0;

      /**
       * Integrate the matter density along trajectory.
       *
       * @return GrammageType as integrated matter density along the BaseTrajectory
       *
       * @todo think about the mixin inheritance of the trajectory vs the BaseTrajectory
       *       approach; for now, only lines are supported (?).
       */
      virtual GrammageType getIntegratedGrammage(BaseTrajectory const&) const = 0;

      /**
       * Calculates the length along the trajectory.
       *
       * The length along the trajectory is determined at which the integrated matter
       * density is reached. If the specified matter density cannot be reached (is too
       * large) the result becomes meaningless and could be "infinity" (discuss this).
       *
       * @return LengthType the length corresponding to grammage.
       */
      virtual LengthType getArclengthFromGrammage(BaseTrajectory const&,
                                                  GrammageType) const = 0;

      virtual media::NuclearComposition const& getNuclearComposition() const = 0;
    };

  } // namespace media
} // namespace corsika
