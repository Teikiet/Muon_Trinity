/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <limits>

#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {
  namespace media {

    /**
     * Helper class to integrate 1D density functions.
     *
     * Integrated density (grammage) can be calculated along a Trajectory. Also the
     * inverse (arclength) for a specified grammage can be calculated. In addition, also a
     * maximum possible length for a maximum allowed uncertainty can be evaluated.
     *
     * @tparam TDerived Must provide the **evaluateAt** and **getFirstDerivative**
     * methods.
     */

    template <typename TDerived>
    class LinearApproximationIntegrator {

      /**
       * Return the type TDerived of this object.
       *
       * @return TDerived const&
       */
      TDerived const& getImplementation() const;

    public:
      /**
       * Get the integrated Grammage along line.
       *
       * Calculated as \f[ X=(\varrho_0 + 0.5* d\varrho/dl \cdot l ) l \f].
       *
       * @param line
       * @return GrammageType
       */
      GrammageType getIntegrateGrammage(BaseTrajectory const& line) const;

      /**
       * Get the arclength from Grammage along line.
       *
       * Calculate as \f[ (1 - 0.5 * X * d\varrho/dl / (\varrho_0^2)) * X / \varrho_0 \f].
       *
       * @param line
       * @param grammage
       * @return LengthType
       */
      LengthType getArclengthFromGrammage(BaseTrajectory const& line,
                                          GrammageType grammage) const;

      /**
       * Get the maximum length l to keep uncertainty below \f[ relError \f].
       *
       * @param line
       * @param relError
       * @return LengthType
       */
      LengthType getMaximumLength(BaseTrajectory const& line, double relError) const;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density/integrator/LinearApproximationIntegrator.inl>
