/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

#include <corsika/media/Environment.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/modules/radio/propagators/RadioPropagator.hpp>

namespace corsika {

  /**
   * This class implements a tabulated propagator that approximates
   * the Earth's atmosphere as flat. Signal propagation is rectilinear
   * and this is intended to be used for vertical showers
   * (<60 degrees zenith angle) for fast simulations. The table needs
   * to have at least 10 elements.
   *
   */
  template <typename TEnvironment>
  class TabulatedFlatAtmospherePropagator final
      : public RadioPropagator<TabulatedFlatAtmospherePropagator<TEnvironment>,
                               TEnvironment> {

    using Base =
        RadioPropagator<TabulatedFlatAtmospherePropagator<TEnvironment>, TEnvironment>;
    using SignalPathCollection = typename Base::SignalPathCollection;

  public:
    /**
     * Construct a new FlatEarthPropagator with a given environment.
     *
     */
    TabulatedFlatAtmospherePropagator(TEnvironment const& env, Point const& upperLimit,
                                      Point const& lowerLimit, LengthType const step);

    /**
     * Return the collection of paths from `source` to `destination`.
     * Hence, the signal propagated from the
     * emission point to the oberserver location.
     *
     */
    template <typename Particle>
    SignalPathCollection propagate(Particle const& particle, Point const& source,
                                   Point const& destination);

  private:
    Point const upperLimit_; ///< the upper point of the table.
    Point const lowerLimit_; ///< the lowest point of the table (ideally earth's surface).
    LengthType const step_;  ///< the tabulation step.
    InverseLengthType const
        inverseStep_;            ///< inverse of the step used to speed up calculations.
    LengthType const maxHeight_; ///< z coordinate of upper limit (maximum height).
    LengthType const minHeight_; ///< z coordinate of lower limit (minimum height) - 1_km
                                 ///< for safety reasons.
    std::vector<double> refractivityTable_; ///< the table that stores refractivity.
    std::vector<LengthType>
        heightTable_; ///< the table that stores the height using the step above.
    std::vector<double>
        integratedRefractivityTable_; ///< the table that stores integrated refractivity.
    double slopeRefrLower_; ///< used to interpolate refractivity for particles below the
                            ///< lowest point.
    double slopeIntRefrLower_; ///< used to interpolate integrated refractivity for
                               ///< particles below the lowest point.
    double slopeRefrUpper_; ///< used to interpolate refractivity for particles above the
                            ///< highest point.
    double slopeIntRefrUpper_; ///< used to interpolate integrated refractivity for
                               ///< particles above the highest point.
    double lastElement_;       ///< last index of tables.
    std::deque<Point> points;  ///< the points that the signal has propagated through.
    std::vector<double>
        rindex; ///< the refractive index values along the signal propagation path.

  }; // End: FlatEarthPropagator

  template <typename TEnvironment>
  TabulatedFlatAtmospherePropagator<TEnvironment>
  make_tabulated_flat_atmosphere_radio_propagator(TEnvironment const& env,
                                                  Point const& upperLimit,
                                                  Point const& lowerLimit,
                                                  LengthType const step) {
    return TabulatedFlatAtmospherePropagator<TEnvironment>(env, upperLimit, lowerLimit,
                                                           step);
  }

} // namespace corsika

#include <corsika/detail/modules/radio/propagators/TabulatedFlatAtmospherePropagator.inl>
