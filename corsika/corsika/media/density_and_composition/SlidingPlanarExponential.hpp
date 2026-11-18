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
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {
  namespace media {

    // clang-format off
  /**
   * The SlidingPlanarExponential models mass density as:
   * \f[
   *   \varrho(r) = \varrho_0 \exp\left( \frac{|p_0 - r|}{\lambda} \right).
   * \f]
   * For grammage/length conversion, the density distribution is approximated as
   * locally flat at the starting point \f$ r_0 \f$ of the trajectory with the axis pointing
   * from \f$ p_0 \f$ to \f$ r_0 \f$.
   */
    // clang-format on

    template <typename T>
    class SlidingPlanarExponential : public BaseExponential<SlidingPlanarExponential<T>>,
                                     public T {

      using Base = BaseExponential<SlidingPlanarExponential<T>>;

    public:
      SlidingPlanarExponential(Point const& p0, MassDensityType const rho0,
                               LengthType const lambda,
                               media::NuclearComposition const& nuclComp,
                               LengthType const referenceHeight = LengthType::zero());

      MassDensityType getMassDensity(Point const& point) const override;

      media::NuclearComposition const& getNuclearComposition() const override;

      GrammageType getIntegratedGrammage(BaseTrajectory const& line) const override;

      LengthType getArclengthFromGrammage(BaseTrajectory const& line,
                                          GrammageType const grammage) const override;

    private:
      media::NuclearComposition const nuclComp_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density_and_composition/SlidingPlanarExponential.inl>
