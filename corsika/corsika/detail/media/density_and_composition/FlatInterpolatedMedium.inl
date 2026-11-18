/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

namespace corsika {
  namespace media {

    template <typename T>
    inline FlatInterpolatedMedium<T>::FlatInterpolatedMedium(
        Point const& point, LengthType altStart, std::vector<LengthType> const& altVec,
        std::vector<MassDensityType> const& densityVec,
        NuclearComposition const& nuclComp)
        : point_(point)
        , nuclComp_(nuclComp)
        , altStart_(altStart) {

      // get upper medium boundary
      LengthType const altEnd = altStart_ + altVec.at(altVec.size() - 1);

      // calculate length step size
      lengthStep_ = (altEnd - altStart_) / steps_;
      invLengthStep_ = 1.0 / lengthStep_;

      // initialize iterators for accessing altitudes and densities - start at the
      // atmosphere top
      auto iterAltHi = altVec.rbegin();
      auto iterAltLo = iterAltHi + 1;
      auto iterDensHi = densityVec.rbegin();
      auto iterDensLo = iterDensHi + 1;

      // log values of altitudes and densities
      double logAltHi = log10((*iterAltHi + altStart_) / 1_m);
      double logAltLo = log10((*iterAltLo + altStart_) / 1_m);
      double logDensHi = log10(*iterDensHi / 1_g * cube(1_cm));
      double logDensLo = log10(*iterDensLo / 1_g * cube(1_cm));
      double logDensSlope = (logDensLo - logDensHi) / (logAltHi - logAltLo);

      // generate lookup tables for density and grammage
      for (std::size_t idx = 0; idx < steps_; idx++) {
        // calculate altitude corresponding to this lookup index
        LengthType const alt = altEnd - (idx + 1) * lengthStep_;

        // ensure iterators point to closest table entries above and below this altitude
        while (alt < (*iterAltLo + altStart_)) {
          // current altitude not between altitude iterators, advance all iterators
          iterAltHi++;
          iterAltLo++;
          iterDensHi++;
          iterDensLo++;

          // recalculate log values of altitudes and densities
          logAltHi = log10((*iterAltHi + altStart_) / 1_m);
          logAltLo = log10((*iterAltLo + altStart_) / 1_m);
          logDensHi = log10(*iterDensHi / 1_g * cube(1_cm));
          logDensLo = log10(*iterDensLo / 1_g * cube(1_cm));
          logDensSlope = (logDensLo - logDensHi) / (logAltHi - logAltLo);
        }

        // logDistance from the top of the current layer
        double const logDistInLayer = (logAltHi - log10(alt / 1_m));

        // get logDensity by interpolating in log. space
        double const logDensity = logDensHi + logDensSlope * logDistInLayer;

        // convert from logDensity to density
        MassDensityType const density = std::pow(10, logDensity) * 1_g / cube(1_cm);

        // initialize thickness to zero
        auto thickness = GrammageType::zero();

        // calculate thickness for every step except the first one (first one = atmosphere
        // top, thickness is always zero)
        if (idx > 0)
          // thickness T for this step calculated using previous thickness T0:
          // T = T0 + avgDensity * lengthStep
          thickness = thicknessTable_[steps_ - idx] +
                      0.5 * (density + densityTable_[steps_ - idx]) * lengthStep_;

        // fill lookup tables
        densityTable_[steps_ - idx - 1] = density;
        thicknessTable_[steps_ - idx - 1] = thickness;
      }

      // total atmosphere logThickness
      double const logThickTotal = log10(thicknessTable_[0] / 1_g * square(1_cm));

      // atmosphere thickness after one step from the top (the top has zero thickness,
      // cannot use log)
      thickMin_ = thicknessTable_[steps_ - 2];
      logThickMin_ = log10(thickMin_ / 1_g * square(1_cm));

      // determine thickness step in log. space
      logThickStep_ = (logThickTotal - logThickMin_) / steps_;
      invLogThickStep_ = 1.0 / logThickStep_;

      // initialize thickness iterators
      auto iterThickHi = thicknessTable_.rbegin();
      auto iterThickLo = iterThickHi + 1;

      // distance from atmosphere top to current layer
      LengthType distToLayer = lengthStep_;

      // calculate arclength lookup table entries
      for (std::size_t idx = 0; idx < steps_; idx++) {
        // calculate logThickness corresponding to this lookup index
        double const logThickness = logThickMin_ + idx * logThickStep_;

        // convert logThickness to thickness
        GrammageType const thickness = pow(10, logThickness) * 1_g / square(1_cm);

        // ensure lookup index points to closest table entry to this thickness
        while (thickness > *iterThickLo && iterThickLo != (thicknessTable_.rend() - 1)) {
          // current altitude not between altitude iterators, advance all iterators
          iterThickHi++;
          iterThickLo++;
          // update distance from atmosphere top to current layer
          distToLayer += lengthStep_;
        }

        // get distance from the top of the current layer by interpolating
        LengthType const distInLayer =
            (thickness - *iterThickHi) / (*iterThickLo - *iterThickHi) * lengthStep_;

        // current altitude is atmosphere top minus distance to the current layer minus
        // distance inside the layer
        LengthType altitude = (altEnd - altStart_) - distToLayer - distInLayer;

        // special treatment for first or last element of the lookup table
        if (idx == 0) { altitude = altEnd - altStart_; }
        if (idx == steps_ - 1) { altitude = LengthType::zero(); }

        // fill lookup table
        altitudeTable_[idx] = altitude;
      }
    }

    template <typename T>
    inline MassDensityType FlatInterpolatedMedium<T>::getMassDensity(
        Point const& p0) const {
      return getMassDensity((p0 - point_).getNorm());
    }

    template <typename T>
    inline MassDensityType FlatInterpolatedMedium<T>::getMassDensity(
        LengthType const height) const {

      // lookup table index
      std::size_t idx = (height - altStart_) * invLengthStep_;

      // ensure index falls within the lookup table range
      idx = std::min(idx, steps_ - 1);

      return densityTable_[idx];
    }

    template <typename T>
    inline media::NuclearComposition const&
    FlatInterpolatedMedium<T>::getNuclearComposition() const {
      return nuclComp_;
    }

    template <typename T>
    inline GrammageType FlatInterpolatedMedium<T>::getIntegratedGrammage(
        BaseTrajectory const& track) const {

      // zero track length, zero grammage
      if (track.getLength() == LengthType::zero()) { return GrammageType::zero(); }

      // get altitudes of track endpoints
      LengthType const alt0 = (track.getPosition(0) - point_).getNorm();

      // correction factor for track inclination: |cos(angle)| between track and axis,
      // calculated using the dot product
      LengthType const uDotA = track.getDirection(0).dot(track.getPosition(0) - point_);
      double const corrFactor = std::abs(uDotA / alt0);

      // special treatment for horizontal or near-horizontal tracks
      // factor of 1.75e-4 corresponds to 0.01deg off horizontal
      if (corrFactor < 1.75e-4) { return getMassDensity(alt0) * track.getLength(); }

      // altitude of track endpoint in flat approximation - add or subtract track length
      // projection
      LengthType const alt1 = alt0 + (uDotA < 0_m ? -track.getLength() * corrFactor
                                                  : track.getLength() * corrFactor);

      // lookup indices for track endpoints
      double idx0 = (alt0 - altStart_) * invLengthStep_;
      double idx1 = (alt1 - altStart_) * invLengthStep_;

      // ensure indices fall within the lookup table range
      idx0 = std::clamp(idx0, 0.0, static_cast<double>(steps_ - 1));
      idx1 = std::clamp(idx1, 0.0, static_cast<double>(steps_ - 1));

      // integer indices for closest table entries
      std::size_t const idx0Lo = std::floor(idx0);
      std::size_t const idx0Hi = std::ceil(idx0);
      std::size_t const idx1Lo = std::floor(idx1);
      std::size_t const idx1Hi = std::ceil(idx1);

      // closest thickness table entries
      GrammageType const thick0Lo = thicknessTable_[idx0Lo];
      GrammageType const thick0Hi = thicknessTable_[idx0Hi];
      GrammageType const thick1Lo = thicknessTable_[idx1Lo];
      GrammageType const thick1Hi = thicknessTable_[idx1Hi];

      // interpolate start- and endpoint thicknesses
      GrammageType const thick0 = thick0Lo + (thick0Hi - thick0Lo) * (idx0 - idx0Lo);
      GrammageType const thick1 = thick1Lo + (thick1Hi - thick1Lo) * (idx1 - idx1Lo);

      // thickness difference between endpoints
      double const thickDiff = std::abs((thick0 - thick1) / 1_g * square(1_cm));

      // grammage along the track is the thickness difference scaled by the correction
      // factor
      return thickDiff / corrFactor * 1_g / square(1_cm);
    }

    template <typename T>
    inline LengthType FlatInterpolatedMedium<T>::getArclengthFromGrammage(
        BaseTrajectory const& track, GrammageType const grammage) const {

      // zero grammage, zero track length
      if (grammage == GrammageType::zero()) { return LengthType::zero(); }

      // get track startpoint altitude
      LengthType const alt0 = (track.getPosition(0) - point_).getNorm();

      // correction factor for track inclination: |cos(angle)| between track and axis,
      // calculated using the dot product
      LengthType const uDotA = track.getDirection(0).dot(track.getPosition(0) - point_);
      double const corrFactor = std::abs(uDotA / alt0);

      // special treatment for horizontal or near-horizontal tracks
      // factor of 1.75e-4 corresponds to 0.01deg off horizontal
      if (corrFactor < 1.75e-4) { return grammage / getMassDensity(alt0); }

      // lookup index for track startpoint
      double idx = (alt0 - altStart_) * invLengthStep_;

      // ensure index fall within the lookup table range
      idx = std::clamp(idx, 0.0, static_cast<double>(steps_ - 1));

      // integer indices for closest table entries
      std::size_t const idxHi = std::ceil(idx);
      std::size_t const idxLo = std::floor(idx);

      // closest thickness table entries
      GrammageType const thickHi = thicknessTable_[idxHi];
      GrammageType const thickLo = thicknessTable_[idxLo];

      // interpolate thickness at track startpoint
      GrammageType const thick0 = thickLo + (thickHi - thickLo) * (idx - idxLo);

      // determine atmosphere thickness at track endpoint by adding or subtracting
      // grammage along the track (downward track - grammage is added, upward
      // track - grammage is subtracted)
      GrammageType thick1 =
          thick0 +
          (uDotA < LengthType::zero() ? grammage * corrFactor : -grammage * corrFactor);

      // prevent negative atmosphere thickness
      thick1 = std::max(thick1, thickMin_);

      // convert to logThickness
      double const logThick0 = log10(thick0 / 1_g * square(1_cm));
      double const logThick1 = log10(thick1 / 1_g * square(1_cm));

      // lookup indices for track endpoints
      double idx0 = (logThick0 - logThickMin_) * invLogThickStep_;
      double idx1 = (logThick1 - logThickMin_) * invLogThickStep_;

      // ensure indices fall within the lookup table range
      idx0 = std::clamp(idx0, 0.0, static_cast<double>(steps_ - 1));
      idx1 = std::clamp(idx1, 0.0, static_cast<double>(steps_ - 1));

      // check index difference
      if (std::abs(idx0 - idx1) < 1) {
        // for close indices, density changes very little - take density at track
        // startpoint and simplify by trackLength = grammage / density
        return grammage / getMassDensity(alt0);
      }

      // integer indices for closest table entries
      std::size_t const idx1Lo = std::floor(idx1);
      std::size_t const idx1Hi = std::ceil(idx1);

      // closest altitude table entries
      LengthType const alt1Lo = altitudeTable_[idx1Lo];
      LengthType const alt1Hi = altitudeTable_[idx1Hi];

      // interpolate to get endpoint altitude
      LengthType const alt1 = alt1Lo + (alt1Hi - alt1Lo) * (idx1 - idx1Lo);

      // altitude difference between endpoints
      double const altDiff = std::abs((alt0 - alt1 - altStart_) / 1_m);

      // total track length is the difference between endpoint arclengths,
      // scaled by the track angle
      return altDiff / corrFactor * 1_m;
    }
  } // namespace media
} // namespace corsika
