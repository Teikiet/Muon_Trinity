/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/media/density/integrator/LinearApproximationIntegrator.hpp>

namespace corsika {
  namespace media {

    template <class TDerivableRho,
              template <typename> class TIntegrator = LinearApproximationIntegrator>
    class DensityFunction
        : public TIntegrator<DensityFunction<TDerivableRho, TIntegrator>> {
      friend class TIntegrator<DensityFunction<TDerivableRho, TIntegrator>>;

      TDerivableRho rho_; //!< functor for density

    public:
      DensityFunction(TDerivableRho rho)
          : rho_(rho) {}

      MassDensityType evaluateAt(Point const& p) const { return rho_(p); }
    };

  } // namespace media
} // namespace corsika
