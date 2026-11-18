/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/density/integrator/LinearApproximationIntegrator.hpp>

namespace corsika {

  namespace media {

    template <typename TDerived>
    inline TDerived const& LinearApproximationIntegrator<TDerived>::getImplementation()
        const {
      return *static_cast<TDerived const*>(this);
    }

    template <typename TDerived>
    inline GrammageType LinearApproximationIntegrator<TDerived>::getIntegrateGrammage(
        BaseTrajectory const& line) const {
      LengthType const length = line.getLength();
      auto const c0 = getImplementation().evaluateAt(line.getPosition(0));
      auto const c1 = getImplementation().rho_.getFirstDerivative(line.getPosition(0),
                                                                  line.getDirection(0));
      CORSIKA_LOG_INFO("length={} c0={} c1={} pos={} dir={} return={}", length, c0, c1,
                       line.getPosition(0), line.getDirection(0),
                       (c0 + 0.5 * c1 * length) * length);
      return (c0 + 0.5 * c1 * length) * length;
    }

    template <typename TDerived>
    inline LengthType LinearApproximationIntegrator<TDerived>::getArclengthFromGrammage(
        BaseTrajectory const& line, GrammageType grammage) const {
      auto const c0 = getImplementation().rho_(line.getPosition(0));
      auto const c1 = getImplementation().rho_.getFirstDerivative(line.getPosition(0),
                                                                  line.getDirection(0));

      return (1 - 0.5 * grammage * c1 / (c0 * c0)) * grammage / c0;
    }

    template <typename TDerived>
    inline LengthType LinearApproximationIntegrator<TDerived>::getMaximumLength(
        BaseTrajectory const& line, [[maybe_unused]] double relError) const {
      [[maybe_unused]] auto const c1 = getImplementation().rho_.getSecondDerivative(
          line.getPosition(0), line.getDirection(0));

      // todo: provide a real, working implementation
      return 1_m * std::numeric_limits<double>::infinity();
    }

  } // namespace media
} // namespace corsika
