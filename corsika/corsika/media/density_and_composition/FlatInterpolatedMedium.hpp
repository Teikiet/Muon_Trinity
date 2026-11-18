/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/BaseTrajectory.hpp>
#include <corsika/framework/geometry/Line.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

namespace corsika {
  namespace media {

    /**
     * A medium defined by mass density values at given altitudes, using interpolation to
     * generate lookup tables
     *
     * The medium uses flat atmosphere approximation - should not be used for very
     * inclined showers
     */
    template <typename T>
    class FlatInterpolatedMedium : public T {

    public:
      /**
       * Construct a new FlatInterpolatedMedium object.
       *
       * @param point Medium center point
       * @param altStart Starting altitude (measured from medium center), typically
       * corresponding to ground level
       * @param alt Vector of altitudes
       * @param density Vector of densities
       * @param nuclComp Nuclear composition object
       * @
       */
      FlatInterpolatedMedium(Point const& point, LengthType altStart,
                             std::vector<LengthType> const& alt,
                             std::vector<MassDensityType> const& density,
                             NuclearComposition const& nuclComp);

      /**
       * Returns the center of the medium
       * @returns Center of the medium
       */
      Point const& getAnchorPoint() const { return point_; }

      /**
       * Returns the mass density at a given point in space
       * @param point Point at which the density is evaluated
       * @returns Density at a given point in space
       */
      MassDensityType getMassDensity(Point const& point) const override;

      /**
       * Returns the mass density at a given height (measured from medium center)
       * @param height Height at which the mass density is evaluated
       * @returns Mass density at a given height
       */
      MassDensityType getMassDensity(LengthType const height) const;

      /**
       * Returns the nuclear composition object
       * @returns Nuclear composition object
       */
      NuclearComposition const& getNuclearComposition() const override;

      /**
       * Returns the integrated grammage along a given track
       * @param track Track along which the integrated grammage is calculated
       * @returns Integrated grammage along a given track
       */
      GrammageType getIntegratedGrammage(BaseTrajectory const& track) const override;

      /**
       * Returns the track length corresponding to a given integrated grammage
       * @param track Track whose length is evaluated
       * @param grammage Integrated grammage
       * @returns Length of the track corresponding to a given integrated grammage
       */
      LengthType getArclengthFromGrammage(BaseTrajectory const& track,
                                          GrammageType const grammage) const override;

    private:
      Point const point_;
      NuclearComposition const nuclComp_;

      // boundary altitude of the atmosphere
      LengthType const altStart_;

      // number of entries in lookup tables - roughly 1m spacing
      static constexpr std::size_t steps_ = 120000;

      // density lookup table - density at given altitude
      std::array<MassDensityType, steps_> densityTable_;
      // thickness lookup table - atmosphere thickness (cumulated density from the top of
      // the atmosphere) at given altitude
      std::array<GrammageType, steps_> thicknessTable_;
      // altitude lookup table - altitude at given atmosphere thickness (reverse of the
      // thickness table)
      std::array<LengthType, steps_> altitudeTable_;

      // lookup table step sizes in length and logThickness
      LengthType lengthStep_{};
      double logThickStep_;
      // inverse step sizes
      InverseLengthType invLengthStep_{};
      double invLogThickStep_;

      // thickness of the first non-zero table entry, used in lookups
      GrammageType thickMin_;
      double logThickMin_;
    };
  } // namespace media
} // namespace corsika

#include <corsika/detail/media/density_and_composition/FlatInterpolatedMedium.inl>
