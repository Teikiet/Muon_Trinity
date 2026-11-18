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
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/density/BaseTabular.hpp>

namespace corsika {
  namespace media {

    // clang-format off
  /**
   * The SlidingPlanarTabular models mass density as
   * \f[
   *   \varrho(r) = \varrho_0 \rho\left( |p_0 - r| \right).
   * \f]
   * For grammage/length conversion, the density distribution is approximated as
   * locally flat at the starting point \f$ r_0 \f$ of the trajectory with the 
   * axis pointing rom \f$ p_0 \f$ to \f$ r_0 \f$ defining the local height.
   */
    // clang-format on

    template <typename TDerived>
    class SlidingPlanarTabular : public BaseTabular<SlidingPlanarTabular<TDerived>>,
                                 public TDerived {

      using Base = BaseTabular<SlidingPlanarTabular<TDerived>>;

    public:
      SlidingPlanarTabular(Point const& p0,
                           std::function<MassDensityType(LengthType)> const& rho,
                           unsigned int const nBins, LengthType const deltaHeight,
                           media::NuclearComposition const& nuclComp,
                           LengthType referenceHeight = LengthType::zero());

      MassDensityType getMassDensity(Point const& point) const override;

      media::NuclearComposition const& getNuclearComposition() const override;

      GrammageType getIntegratedGrammage(BaseTrajectory const& line) const override;

      LengthType getArclengthFromGrammage(BaseTrajectory const& line,
                                          GrammageType grammage) const override;

    private:
      media::NuclearComposition const nuclComp_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density_and_composition/SlidingPlanarTabular.inl>
