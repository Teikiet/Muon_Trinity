/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    template <typename T>
    inline SlidingPlanarTabular<T>::SlidingPlanarTabular(
        Point const& p0, std::function<MassDensityType(LengthType)> const& rho,
        unsigned int const nBins, LengthType const deltaHeight,
        media::NuclearComposition const& nuclComp, LengthType referenceHeight)
        : BaseTabular<SlidingPlanarTabular<T>>(p0, referenceHeight, rho, nBins,
                                               deltaHeight)
        , nuclComp_(nuclComp) {}

    template <typename T>
    inline MassDensityType SlidingPlanarTabular<T>::getMassDensity(
        Point const& point) const {
      auto const heightFull =
          (point - BaseTabular<SlidingPlanarTabular<T>>::getAnchorPoint()).getNorm();
      return BaseTabular<SlidingPlanarTabular<T>>::getMassDensity(heightFull);
    }

    template <typename T>
    inline media::NuclearComposition const&
    SlidingPlanarTabular<T>::getNuclearComposition() const {
      return nuclComp_;
    }

    template <typename T>
    inline GrammageType SlidingPlanarTabular<T>::getIntegratedGrammage(
        BaseTrajectory const& traj) const {
      return BaseTabular<SlidingPlanarTabular<T>>::getIntegratedGrammage(traj);
    }

    template <typename T>
    inline LengthType SlidingPlanarTabular<T>::getArclengthFromGrammage(
        BaseTrajectory const& traj, GrammageType const grammage) const {
      return BaseTabular<SlidingPlanarTabular<T>>::getArclengthFromGrammage(traj,
                                                                            grammage);
    }

  } // namespace media
} // namespace corsika
