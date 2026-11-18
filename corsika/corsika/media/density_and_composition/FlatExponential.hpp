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
#include <corsika/media/density/BaseExponential.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>

namespace corsika {
  namespace media {

    // clang-format off
  /**
   * flat exponential density distribution with
   * \f[
   *  \varrho(\vec{r}) = \varrho_0 \exp\left( \frac{1}{\lambda} (\vec{r} - \vec{p}) \cdot
   *    \vec{a} \right).
   * \f]
   * \f$ \vec{a} \f$ denotes the axis and should be normalized to avoid degeneracy
   * with the scale parameter \f$ \lambda \f$, \f$ \vec{r} \f$ is the location of 
   * the evaluation, \f$ \vec{p} \f$ is the anchor point at which \f$ \varrho_0 \f$ 
   * is given. Thus, the unit vector \f$ \vec{a} \f$ specifies the direction of
   * <em>decreasing</em> <b>height/altitude</b>.
   */
    // clang-format on
    template <typename T>
    class FlatExponential : public BaseExponential<FlatExponential<T>>, public T {
      using base_type = BaseExponential<FlatExponential<T>>;

    public:
      FlatExponential(Point const& point, Vector<dimensionless_d> const& axis,
                      MassDensityType const rho, LengthType const lambda,
                      media::NuclearComposition const& nuclComp);

      MassDensityType getMassDensity(Point const& point) const override;

      media::NuclearComposition const& getNuclearComposition() const override;

      GrammageType getIntegratedGrammage(BaseTrajectory const& line) const override;

      LengthType getArclengthFromGrammage(BaseTrajectory const& line,
                                          GrammageType const grammage) const override;

    private:
      DirectionVector const axis_;
      media::NuclearComposition const nuclComp_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density_and_composition/FlatExponential.inl>
