/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>

#include <exception>

namespace corsika {
  namespace media {

    template <typename TDerived>
    inline BaseTabular<TDerived>::BaseTabular(
        Point const& point, LengthType const referenceHeight,
        std::function<MassDensityType(LengthType)> const& rho, unsigned int const nBins,
        LengthType const deltaHeight)
        : nBins_(nBins)
        , deltaHeight_(deltaHeight)
        , point_(point)
        , referenceHeight_(referenceHeight) {
      density_.resize(nBins_);
      for (unsigned int bin = 0; bin < nBins; ++bin) {
        density_[bin] = rho(deltaHeight_ * bin);
        CORSIKA_LOG_DEBUG("new tabulated atm bin={} h={} rho={}", bin, deltaHeight_ * bin,
                          density_[bin]);
      }
    }

    template <typename TDerived>
    inline auto const& BaseTabular<TDerived>::getImplementation() const {
      return *static_cast<TDerived const*>(this);
    }

    template <typename TDerived>
    inline MassDensityType BaseTabular<TDerived>::getMassDensity(
        LengthType const height) const {
      double const fbin = (height - referenceHeight_) / deltaHeight_;
      int const bin = int(fbin);
      if (bin < 0) return MassDensityType::zero();
      if (bin >= int(nBins_ - 1)) {
        CORSIKA_LOG_ERROR(
            "invalid height {} (corrected {}) in BaseTabular atmosphere. Min 0, max {}. "
            "If max is too low: increase!",
            height, height - referenceHeight_, nBins_ * deltaHeight_);
        throw std::runtime_error("invalid height");
      }
      return density_[bin] + (fbin - bin) * (density_[bin + 1] - density_[bin]);
    }

    template <typename TDerived>
    inline GrammageType BaseTabular<TDerived>::getIntegratedGrammage(
        BaseTrajectory const& traj) const {

      Point pCurr = traj.getPosition(0);
      DirectionVector dCurr = traj.getDirection(0);
      LengthType height1 = (traj.getPosition(0) - point_).getNorm() - referenceHeight_;
      LengthType height2 = (traj.getPosition(1) - point_).getNorm() - referenceHeight_;

      LengthType const fullLength = traj.getLength(1);

      int sign = +1; // normal direction
      if (height1 > height2) {
        std::swap(height1, height2);
        pCurr = traj.getPosition(1);
        dCurr = traj.getDirection(1);
        sign = -1; // inverted direction
      }

      double const fbin1 = height1 / deltaHeight_;
      unsigned int const bin1 = int(fbin1);

      double const fbin2 = height2 / deltaHeight_;
      unsigned int const bin2 = int(fbin2);

      if (fbin1 == fbin2) { return GrammageType::zero(); }

      if (bin1 >= nBins_ - 1 || bin2 >= nBins_ - 1) {
        CORSIKA_LOG_ERROR("invalid height {} {} in BaseTabular atmosphere integration",
                          height1, height2);
        throw std::runtime_error("invalid height");
      }

      // interpolated start/end densities
      MassDensityType const rho1 = getMassDensity(height1 + referenceHeight_);
      MassDensityType const rho2 = getMassDensity(height2 + referenceHeight_);

      // within first bin
      if (bin1 == bin2) { return fullLength * (rho2 + rho1) / 2; }

      // inclination of trajectory (local)
      DirectionVector axis((pCurr - point_).normalized()); // to gravity center
      double cosTheta = axis.dot(dCurr);

      // distance to next height bin
      unsigned int bin = bin1;
      LengthType dD = (deltaHeight_ * (bin + 1) - height1) / cosTheta * sign;
      LengthType distance = dD;

      GrammageType X = dD * (rho1 + density_[bin + 1]) / 2;
      double frac = (sign > 0 ? distance / fullLength : 1 - distance / fullLength);
      pCurr = traj.getPosition(frac);
      dCurr = traj.getDirection(frac);

      for (++bin; bin < bin2; ++bin) {
        // inclination of trajectory
        axis = (pCurr - point_).normalized();
        cosTheta = axis.dot(dCurr);
        // distance to next height bin
        dD = deltaHeight_ / cosTheta * sign;
        distance += dD;
        GrammageType const dX = dD * (density_[bin] + density_[bin + 1]) / 2;
        X += dX;
        frac = (sign > 0 ? distance / fullLength : 1 - distance / fullLength);
        pCurr = traj.getPosition(frac);
        dCurr = traj.getDirection(frac);
      }

      // inclination of trajectory
      axis = ((pCurr - point_).normalized());
      cosTheta = axis.dot(dCurr);
      // distance to next height bin
      dD = (height2 - deltaHeight_ * bin2) / cosTheta * sign;
      X += dD * (rho2 + density_[bin2]) / 2;
      distance += dD;
      return X;
    }

    template <typename TDerived>
    inline LengthType BaseTabular<TDerived>::getArclengthFromGrammage(
        BaseTrajectory const& traj, GrammageType const grammage) const {

      if (grammage < GrammageType::zero()) {
        CORSIKA_LOG_ERROR("cannot integrate negative grammage");
        throw std::runtime_error("negative grammage error");
      }
      LengthType const height =
          (traj.getPosition(0) - point_).getNorm() - referenceHeight_;

      double const fbin = height / deltaHeight_;
      int bin = int(fbin);

      if (bin >= int(nBins_ - 1)) {
        CORSIKA_LOG_ERROR("invalid height {} in BaseTabular atmosphere integration",
                          height);
        throw std::runtime_error("invalid height");
      }

      // interpolated start/end densities
      MassDensityType const rho = getMassDensity(height + referenceHeight_);

      // inclination of trajectory
      Point pCurr = traj.getPosition(0);
      DirectionVector dCurr = traj.getDirection(0);
      DirectionVector axis((pCurr - point_).normalized());
      double cosTheta = axis.dot(dCurr);
      int sign = +1; // height increasing along traj
      if (cosTheta < 0) {
        cosTheta = -cosTheta; // absolute value only
        sign = -1;            // height decreasing along traj
      }

      // height -> distance
      LengthType const deltaDistance = deltaHeight_ / cosTheta;

      // start with 0 g/cm2
      GrammageType X(GrammageType::zero());
      LengthType distance(LengthType::zero());

      // within first bin
      distance =
          (sign > 0 ? deltaDistance * (bin + 1 - fbin) : deltaDistance * (fbin - bin));
      GrammageType binGrammage = (sign > 0 ? distance * (rho + density_[bin + 1]) / 2
                                           : distance * (rho + density_[bin]) / 2);
      if (X + binGrammage > grammage) {
        double const binFraction = (grammage - X) / binGrammage;
        return distance * binFraction;
      }
      X += binGrammage;

      // the following bins (along trajectory)
      for (bin += sign; bin < int(nBins_ - 1) && bin >= 0; bin += sign) {

        binGrammage = deltaDistance * (density_[bin] + density_[bin + 1]) / 2;
        if (X + binGrammage > grammage) {
          double const binFraction = (grammage - X) / binGrammage;
          return distance + deltaDistance * binFraction;
        }
        X += binGrammage;
        distance += deltaDistance;
      }
      return std::numeric_limits<double>::infinity() * meter;
    }

  } // namespace media
} // namespace corsika
